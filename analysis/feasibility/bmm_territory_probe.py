"""
Re-run Yang / Class-B ROI territory assignability using real ProbArterialAtlas_BMM.

Prefer Atlas_182_MNI152 (FSL 182×218×182). Mean BMM probability vector inside each
Harvard-Oxford ROI is the confirmatory territory assignment (Aim 1 design).
"""
from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn import datasets, image

ROOT = Path(__file__).resolve().parents[2]
BMM_DIR = (
    ROOT
    / "data"
    / "raw"
    / "arterial_atlas"
    / "nitrc"
    / "extracted"
    / "Atlas_MNI152"
    / "Atlas_182_MNI152"
)
OUT = ROOT / "results" / "feasibility"
OUT.mkdir(parents=True, exist_ok=True)
ORDER = ["ACA", "MCA", "PCA", "VB"]

ROIS = {
    "hippocampus": {"atlas": "subcortical", "labels": ["Left Hippocampus", "Right Hippocampus"], "role": "Yang"},
    "superior_frontal_gyrus": {"atlas": "cortical", "labels": ["Superior Frontal Gyrus"], "role": "Yang"},
    "temporal_pole": {"atlas": "cortical", "labels": ["Temporal Pole"], "role": "Garcia_exvivo_MCA"},
    "middle_frontal_gyrus": {"atlas": "cortical", "labels": ["Middle Frontal Gyrus"], "role": "ClassB_MCA_candidate"},
    "frontal_medial_cortex": {"atlas": "cortical", "labels": ["Frontal Medial Cortex"], "role": "ClassB_ACA_candidate"},
    "angular_gyrus": {"atlas": "cortical", "labels": ["Angular Gyrus"], "role": "ClassB_MCA_candidate"},
    "occipital_pole": {"atlas": "cortical", "labels": ["Occipital Pole"], "role": "ClassA_PCA"},
}


def shannon(p):
    p = np.asarray(p, float)
    p = p[p > 0]
    p = p / p.sum()
    return float(-(p * np.log2(p)).sum())


def main():
    bmm_img = nib.load(str(BMM_DIR / "ProbArterialAtlas_BMM.nii"))
    border_img = nib.load(str(BMM_DIR / "BorderZone_ProbAve.nii"))
    terr_img = nib.load(str(BMM_DIR / "TerritoryVoxels_BMM.nii"))
    bmm = np.asanyarray(bmm_img.dataobj).astype(np.float64)
    border = np.asanyarray(border_img.dataobj).astype(np.float64)
    terr = np.asanyarray(terr_img.dataobj)
    if terr.ndim == 4:
        terr = terr[..., 0]

    # Inspect affine — MNI152 package should have proper affine
    affine_info = {
        "bmm_affine": bmm_img.affine.tolist(),
        "bmm_shape": list(bmm.shape),
        "border_shape": list(border.shape),
        "territory_voxels_unique": sorted(int(x) for x in np.unique(terr)),
    }

    cort = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-2mm")
    sub = datasets.fetch_atlas_harvard_oxford("sub-maxprob-thr25-2mm")
    atlases = {
        "cortical": (cort.maps, list(cort.labels)),
        "subcortical": (sub.maps, list(sub.labels)),
    }

    rows = []
    for name, spec in ROIS.items():
        maps_img, labels = atlases[spec["atlas"]]
        data = np.asanyarray(maps_img.dataobj)
        mask = np.zeros(data.shape, dtype=bool)
        for lab in spec["labels"]:
            mask |= data == labels.index(lab)
        mask_img = nib.Nifti1Image(mask.astype(np.uint8), maps_img.affine)
        mask_rs = image.resample_to_img(
            mask_img, bmm_img, interpolation="nearest", force_resample=True, copy_header=True
        )
        m = np.asanyarray(mask_rs.dataobj).astype(bool)
        # Voxels with any BMM mass
        mass = bmm.sum(axis=-1) > 0
        sel = m & mass
        n = int(sel.sum())
        if n == 0:
            rows.append({"roi": name, "n_voxels": 0, "role": spec["role"], "status": "empty"})
            continue
        mean_p = bmm[sel].mean(axis=0)
        # Renormalize in case of numerical drift
        mean_p = mean_p / mean_p.sum()
        # Borderzone: mean MCA/ACA and MCA/PCA ratios; also fraction of multi-claimed voxels (label 4)
        bz_ratio_mean = border[sel].mean(axis=0) if border.ndim == 4 else [np.nan, np.nan]
        multi = float((terr[sel] == 4).mean()) if n else np.nan
        none = float((terr[sel] == 5).mean()) if n else np.nan
        dom_i = int(np.argmax(mean_p))
        row = {
            "roi": name,
            "role": spec["role"],
            "n_voxels": n,
            "p_ACA": float(mean_p[0]),
            "p_MCA": float(mean_p[1]),
            "p_PCA": float(mean_p[2]),
            "p_VB": float(mean_p[3]),
            "dominant": ORDER[dom_i],
            "p_dominant": float(mean_p[dom_i]),
            "entropy_bits": shannon(mean_p),
            "mean_MCA_over_ACA_ratio": float(bz_ratio_mean[0]),
            "mean_MCA_over_PCA_ratio": float(bz_ratio_mean[1]),
            "frac_multi_claimed_voxels": multi,
            "frac_unclaimed_voxels": none,
            "high_specificity_pdom_ge_0.50": bool(mean_p[dom_i] >= 0.50),
            "high_specificity_H_le_1.0": bool(shannon(mean_p) <= 1.0),
            "method": "ProbArterialAtlas_BMM_mean_in_HO_ROI",
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "bmm_roi_territory_overlap.csv", index=False)

    yang = df[df["role"] == "Yang"]
    hip = yang[yang["roi"] == "hippocampus"].iloc[0]
    sfg = yang[yang["roi"] == "superior_frontal_gyrus"].iloc[0]
    verdict = {
        "method": "Liu ProbArterialAtlas_BMM (Atlas_182_MNI152)",
        "affine_info": affine_info,
        "hippocampus_pca_dominant": bool(hip["dominant"] == "PCA" and hip["p_dominant"] >= 0.50),
        "superior_frontal_aca_dominant": bool(sfg["dominant"] == "ACA" and sfg["p_dominant"] >= 0.50),
        "yang_within_donor_viable": bool(
            hip["dominant"] == "PCA"
            and sfg["dominant"] == "ACA"
            and hip["p_dominant"] >= 0.50
            and sfg["p_dominant"] >= 0.50
        ),
        "recommended_yang_contrast": f"{hip['dominant']} (hippocampus) vs {sfg['dominant']} (superior frontal)",
        "garcia_temporal_mca": None,
    }
    tp = df[df["roi"] == "temporal_pole"]
    if len(tp):
        t = tp.iloc[0]
        verdict["garcia_temporal_mca"] = bool(t["dominant"] == "MCA" and t["p_dominant"] >= 0.50)

    (OUT / "bmm_feasibility_verdict.json").write_text(json.dumps(verdict, indent=2), encoding="utf-8")
    md = [
        "# BMM territory feasibility — confirmatory",
        "",
        "Using real `ProbArterialAtlas_BMM` from NITRC Atlas_182_MNI152.",
        "",
        "## Verdict",
        "",
        "```json",
        json.dumps(verdict, indent=2),
        "```",
        "",
        "## ROI table",
        "",
        df.to_markdown(index=False),
        "",
    ]
    (OUT / "bmm_roi_territory_summary.md").write_text("\n".join(md), encoding="utf-8")
    print(json.dumps(verdict, indent=2))
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
