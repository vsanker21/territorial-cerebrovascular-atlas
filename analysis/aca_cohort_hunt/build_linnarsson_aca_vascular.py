"""
Build external ACA cortex vascular cohort from Linnarsson / Siletti HBCA
CELLxGENE census dissections.

Primary ACA pool (TerritoryVoxels ACA-dominant HO ROIs):
  A24 anterior cingulate, A32 rostral MFC, A25 subcallosal, A23 caudal cingulate.

Exclude A14 (gyrus rectus / medial OFC) — MCA on TerritoryVoxels Frontal Orbital.

Outputs:
  data/processed/aca_cohort_hunt/linnarsson_aca_vascular.h5ad
  results/aca_cohort_hunt/linnarsson_aca_vascular_summary.json
  results/aca_cohort_hunt/linnarsson_aca_tv_priors.json
"""
from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import anndata as ad
import nibabel as nib
import numpy as np
import pandas as pd
import scanpy as sc
from nilearn import datasets, image

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/aca_cohort_hunt/linnarsson_acc"
OUT_H5AD = ROOT / "data/processed/aca_cohort_hunt/linnarsson_aca_vascular.h5ad"
OUT_DIR = ROOT / "results/aca_cohort_hunt"
OUT_SUMMARY = OUT_DIR / "linnarsson_aca_vascular_summary.json"
OUT_TV = OUT_DIR / "linnarsson_aca_tv_priors.json"

TV_PATH = ROOT / (
    "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/"
    "Atlas_182_MNI152/TerritoryVoxels_BMM.nii"
)

# Brodmann / dissection -> Harvard-Oxford cortical labels for TerritoryVoxels
DISSECTIONS = {
    "A24": {
        "file": "Cx_CgGr_A24.h5ad",
        "aca_pool": True,
        "ho_labels": ["Cingulate Gyrus, anterior division"],
        "role": "primary_ACA_cingulate",
    },
    "A32": {
        "file": "Cx_A32_RoG.h5ad",
        "aca_pool": True,
        "ho_labels": ["Paracingulate Gyrus", "Cingulate Gyrus, anterior division"],
        "role": "ACA_dorsal_MFC",
    },
    "A25": {
        "file": "Cx_A25_SCG.h5ad",
        "aca_pool": True,
        "ho_labels": ["Frontal Medial Cortex", "Cingulate Gyrus, anterior division"],
        "role": "ACA_subcallosal_borderline",
    },
    "A23": {
        "file": "Cx_A23_CgGC.h5ad",
        "aca_pool": True,
        "ho_labels": ["Cingulate Gyrus, posterior division"],
        "role": "ACA_caudal_cingulate",
    },
    "A14": {
        "file": "Cx_A14_ReG.h5ad",
        "aca_pool": False,
        "ho_labels": ["Frontal Orbital Cortex"],
        "role": "EXCLUDE_MCA_medial_OFC",
    },
}

VASC_TYPES = {
    "endothelial cell",
    "pericyte",
    "vascular associated smooth muscle cell",
    "fibroblast",
}
TV_MAP = {1: "ACA", 2: "MCA", 3: "PCA", 4: "VB"}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def tv_roi_vector(ho_labels: list[str]) -> dict:
    ho = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-1mm")
    labels = list(ho["labels"])
    ho_data = np.asanyarray(image.load_img(ho["maps"]).dataobj)
    mask = np.zeros(ho_data.shape, dtype=bool)
    missing = []
    for lab in ho_labels:
        if lab not in labels:
            missing.append(lab)
            continue
        mask |= ho_data == labels.index(lab)
    if not mask.any():
        raise KeyError(f"No HO voxels for {ho_labels}; missing={missing}")

    tv_img = nib.load(str(TV_PATH))
    if tv_img.shape[:3] != ho_data.shape:
        tv_img = image.resample_to_img(tv_img, ho["maps"], interpolation="nearest")
    tv = np.rint(np.asanyarray(tv_img.dataobj)).astype(np.int32)
    vals = tv[mask]
    assigned = vals[np.isin(vals, [1, 2, 3, 4])]
    frac_asg = float(len(assigned) / max(1, len(vals)))
    if len(assigned) == 0:
        return {
            "tv_dominant": "unassigned",
            "tv_p_ACA": 0.0,
            "tv_p_MCA": 0.0,
            "tv_p_PCA": 0.0,
            "tv_p_VB": 0.0,
            "tv_p_dominant": 0.0,
            "tv_assignment_fraction": frac_asg,
            "ho_missing": missing,
        }
    counts = {TV_MAP[k]: float((assigned == k).sum()) for k in TV_MAP}
    total = float(len(assigned))
    probs = {f"tv_p_{k}": counts[k] / total for k in counts}
    dom = max(counts, key=counts.get)
    return {
        "tv_dominant": dom,
        **probs,
        "tv_p_dominant": probs[f"tv_p_{dom}"],
        "tv_assignment_fraction": frac_asg,
        "ho_missing": missing,
        "tv_source": "TerritoryVoxels_HO",
    }


def extract_vascular(path: Path, dissection: str, meta: dict) -> ad.AnnData:
    log(f"Reading {path.name}")
    a = sc.read_h5ad(path)
    keep = a.obs["cell_type"].astype(str).isin(VASC_TYPES)
    b = a[keep].copy()
    b.obs["dissection_code"] = dissection
    b.obs["dissection_role"] = meta["role"]
    b.obs["aca_pool"] = meta["aca_pool"]
    b.obs["source_dataset"] = "linnarsson_siletti_hbca_cellxgene"
    b.obs["cohort_arm"] = "ACA_cortex" if meta["aca_pool"] else "MCA_medial_OFC"
    return b


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)

    tv_priors = {}
    for code, meta in DISSECTIONS.items():
        log(f"TerritoryVoxels for {code}")
        tv_priors[code] = {
            **tv_roi_vector(meta["ho_labels"]),
            "ho_labels": meta["ho_labels"],
            "aca_pool": meta["aca_pool"],
            "role": meta["role"],
        }

    OUT_TV.write_text(json.dumps(tv_priors, indent=2), encoding="utf-8")
    log(f"Wrote {OUT_TV}")

    pieces = []
    per_dissection = {}
    for code, meta in DISSECTIONS.items():
        path = RAW / meta["file"]
        if not path.exists():
            raise FileNotFoundError(path)
        b = extract_vascular(path, code, meta)
        tv = tv_priors[code]
        for k, v in tv.items():
            if k.startswith("tv_") or k in {"ho_labels", "role"}:
                b.obs[k] = json.dumps(v) if isinstance(v, list) else v
        per_dissection[code] = {
            "n_vascular": int(b.n_obs),
            "n_endo": int((b.obs["cell_type"].astype(str) == "endothelial cell").sum()),
            "n_peri": int((b.obs["cell_type"].astype(str) == "pericyte").sum()),
            "n_vSMC": int(
                (b.obs["cell_type"].astype(str) == "vascular associated smooth muscle cell").sum()
            ),
            "n_fib": int((b.obs["cell_type"].astype(str) == "fibroblast").sum()),
            "n_donors": int(b.obs["donor_id"].nunique()) if b.n_obs else 0,
            "aca_pool": meta["aca_pool"],
            "tv_dominant": tv["tv_dominant"],
            "tv_p_ACA": tv["tv_p_ACA"],
            "tv_assignment_fraction": tv["tv_assignment_fraction"],
        }
        pieces.append(b)

    all_v = ad.concat(pieces, join="outer", merge="same", index_unique="-")
    aca = all_v[all_v.obs["aca_pool"].astype(bool)].copy()
    aca.write_h5ad(OUT_H5AD)
    log(f"Wrote {OUT_H5AD} n={aca.n_obs}")

    summary = {
        "object": str(OUT_H5AD.relative_to(ROOT)),
        "n_vascular_aca_pool": int(aca.n_obs),
        "n_endo_aca_pool": int((aca.obs["cell_type"].astype(str) == "endothelial cell").sum()),
        "n_donors_aca_pool": int(aca.obs["donor_id"].nunique()),
        "cell_type_counts": aca.obs["cell_type"].astype(str).value_counts().to_dict(),
        "dissection_counts": aca.obs["dissection_code"].astype(str).value_counts().to_dict(),
        "per_dissection": per_dissection,
        "excluded": {
            "A14": "MCA-leaning Frontal Orbital on TerritoryVoxels; kept out of ACA pool"
        },
        "verdict": (
            "GO_as_external_ACA_label_anchor_but_n_sparse"
            if aca.n_obs >= 50
            else "NO_GO_too_sparse"
        ),
        "notes": [
            "Adult neurotypical only in these census slices.",
            "Not vessel-enriched; vascular yield is low per dissection.",
            "Use as TerritoryVoxels-mapped ACA cortex vascular seed; prefer LIBD dACC "
            "(GSE296731) once annotated vascular gate is available for power.",
        ],
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
