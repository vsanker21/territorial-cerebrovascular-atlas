"""
Attach arterial-territory probability vectors to Yang samples.

Each Yang sample inherits the ROI×territory distribution of its dissection
region (hippocampus or superior frontal cortex), computed in Liu atlas space.

Primary vectors use hard-atlas voxel fractions (same as feasibility probe).
If PROXY soft maps exist, also report mean soft probabilities and mean
borderzone fraction inside the ROI — clearly labeled PROXY_NOT_BMM.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn import datasets, image

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis" / "feasibility"))
from yang_territory_probe import (  # noqa: E402
    FSL_MNI152_1MM_AFFINE,
    LEVEL2_MAJOR,
    load_arterial_atlas,
    resample_mask_to_atlas,
    roi_mask_from_labels,
    shannon_entropy,
    territory_vector,
)

META = ROOT / "data" / "processed" / "yang" / "yang_sample_metadata.csv"
ATLAS_DIR = ROOT / "data" / "raw" / "arterial_atlas" / "Atlas_182"
PROXY_DIR = ROOT / "data" / "processed" / "arterial_proxy"
BMM_PATH = (
    ROOT
    / "data"
    / "raw"
    / "arterial_atlas"
    / "nitrc"
    / "extracted"
    / "Atlas_MNI152"
    / "Atlas_182_MNI152"
    / "ProbArterialAtlas_BMM.nii"
)
BORDER_PATH = BMM_PATH.parent / "BorderZone_ProbAve.nii"
TERR_PATH = BMM_PATH.parent / "TerritoryVoxels_BMM.nii"
OUT = ROOT / "data" / "processed" / "yang" / "yang_sample_territory_vectors.csv"
OUT_JSON = ROOT / "data" / "processed" / "yang" / "yang_region_territory_lookup.json"

REGION_TO_HO = {
    "hippocampus": {
        "atlas": "subcortical",
        "labels": ["Left Hippocampus", "Right Hippocampus"],
    },
    "superior_frontal_cortex": {
        "atlas": "cortical",
        "labels": ["Superior Frontal Gyrus"],
    },
}


def load_ho():
    cort = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-2mm")
    sub = datasets.fetch_atlas_harvard_oxford("sub-maxprob-thr25-2mm")
    return {
        "cortical": (cort.maps, list(cort.labels)),
        "subcortical": (sub.maps, list(sub.labels)),
    }


def region_hard_vector(region: str, atlases, art_l2) -> dict:
    spec = REGION_TO_HO[region]
    maps_img, labels = atlases[spec["atlas"]]
    mask = roi_mask_from_labels(maps_img, labels, spec["labels"])
    mask_rs = resample_mask_to_atlas(mask, maps_img.affine, art_l2)
    l2 = np.asanyarray(art_l2.dataobj)
    mask_rs = mask_rs & (l2 > 0) & (l2 < 9)
    return territory_vector(mask_rs, l2, LEVEL2_MAJOR)


def region_proxy_stats(region: str, atlases, art_l2) -> dict | None:
    prob_path = PROXY_DIR / "PROXY_soft_territory_probs.nii.gz"
    border_path = PROXY_DIR / "PROXY_borderzone_fraction.nii.gz"
    if not prob_path.exists():
        return None
    spec = REGION_TO_HO[region]
    maps_img, labels = atlases[spec["atlas"]]
    mask = roi_mask_from_labels(maps_img, labels, spec["labels"])
    mask_rs = resample_mask_to_atlas(mask, maps_img.affine, art_l2)
    probs = nib.load(str(prob_path)).get_fdata()
    border = nib.load(str(border_path)).get_fdata()
    if probs.ndim == 4 and probs.shape[-1] == 4:
        pass
    else:
        raise ValueError(probs.shape)
    # Restrict to voxels with any soft mass
    mass = probs.sum(axis=-1) > 0
    sel = mask_rs & mass
    if sel.sum() == 0:
        return None
    mean_p = probs[sel].mean(axis=0)
    mean_p = mean_p / mean_p.sum()
    return {
        "proxy_p_ACA": float(mean_p[0]),
        "proxy_p_MCA": float(mean_p[1]),
        "proxy_p_PCA": float(mean_p[2]),
        "proxy_p_VB": float(mean_p[3]),
        "proxy_borderzone_fraction_mean": float(border[sel].mean()),
        "proxy_entropy_bits": shannon_entropy(mean_p),
        "proxy_type": "PROXY_NOT_BMM",
    }


def main():
    if not META.exists():
        raise SystemExit(f"Missing {META}; run parse_geo_metadata.py first")
    samples = pd.read_csv(META)
    art_l2 = load_arterial_atlas(ATLAS_DIR / "ArterialAtlas_level2.nii")
    atlases = load_ho()

    lookup = {}
    for region in sorted(samples["region"].unique()):
        if region not in REGION_TO_HO:
            continue
        hard = region_hard_vector(region, atlases, art_l2)
        proxy = region_proxy_stats(region, atlases, art_l2)
        entry = {
            "region": region,
            "n_voxels": hard["n_voxels"],
            "p_ACA": hard["probs"].get("ACA", np.nan),
            "p_MCA": hard["probs"].get("MCA", np.nan),
            "p_PCA": hard["probs"].get("PCA", np.nan),
            "p_VB": hard["probs"].get("VB", np.nan),
            "dominant": hard["dominant"],
            "p_dominant": hard["p_dominant"],
            "entropy_bits": hard["entropy_bits"],
            "method": "hard_atlas_voxel_fraction",
            "high_specificity_primary": bool(
                hard["p_dominant"] >= 0.50 and hard["entropy_bits"] <= 1.0
            ),
        }
        if proxy:
            entry.update(proxy)
        lookup[region] = entry

    rows = []
    for _, s in samples.iterrows():
        region = s["region"]
        base = lookup.get(region, {})
        row = {
            "geo_accession": s.get("geo_accession"),
            "title": s.get("title"),
            "donor_id": s.get("donor_id"),
            "disease": s.get("disease"),
            "region": region,
            "territory_method": base.get("method"),
            "dominant_territory": base.get("dominant"),
            "p_dominant": base.get("p_dominant"),
            "entropy_bits": base.get("entropy_bits"),
            "p_ACA": base.get("p_ACA"),
            "p_MCA": base.get("p_MCA"),
            "p_PCA": base.get("p_PCA"),
            "p_VB": base.get("p_VB"),
            "high_specificity_primary": base.get("high_specificity_primary"),
            "analysis_stream": (
                "H1_primary"
                if base.get("high_specificity_primary")
                else ("Aim3_borderzone" if region in REGION_TO_HO else "unmapped")
            ),
        }
        for k in (
            "proxy_p_ACA",
            "proxy_p_MCA",
            "proxy_p_PCA",
            "proxy_p_VB",
            "proxy_borderzone_fraction_mean",
            "proxy_entropy_bits",
            "proxy_type",
        ):
            row[k] = base.get(k)
        rows.append(row)

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT, index=False)
    OUT_JSON.write_text(json.dumps(lookup, indent=2), encoding="utf-8")
    print(out_df.to_string(index=False))
    print("\nRegion lookup:")
    print(json.dumps(lookup, indent=2))
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
