"""
Yang → arterial territory feasibility probe
==========================================

Scientific purpose
------------------
Section 9 of the Territorial Cerebrovascular Atlas plan requires a cheap
go/no-go test: map Yang et al. (Nature 2022) dissection regions onto the
Liu et al. (Sci Data 2023) arterial territory atlas and ask whether samples
are territorially assignable with usable specificity.

Critical accuracy corrections relative to the draft plan text
------------------------------------------------------------
1. Yang cortex is *superior frontal cortex*, not generic neocortex
   (GEO GSE163577; Nature 2022). Superior frontal gyrus straddles ACA
   (medial) and MCA frontal pars (lateral). Expect higher entropy /
   borderzone fraction than a pure MCA ROI.
2. Hippocampus is classically PCA / anterior-choroidal–associated, not a
   pure cytoarchitectonic match to "PCA cortex".
3. GitHub Chin-Fu-Liu/Arterial_Atlas currently ships hard parcellations
   (ArterialAtlas.nii, ArterialAtlas_level2.nii), not ProbArterialAtlas_BMM
   or BorderZone_ProbAve. This probe therefore estimates territory
   probability vectors as *voxel-fraction distributions within anatomical
   ROIs* after resampling ROIs into atlas space. That is a principled
   first-order proxy for the Bernoulli-mixture maps; BMM maps should
   replace this step once obtained from NITRC.

Go criteria (pre-specified for interpretation)
----------------------------------------------
- Hippocampus dominant territory is PCA (level-2, bilateral pooled) with
  P(dominant) >= 0.50 and Shannon entropy H <= 1.0 bits (4 major territories).
- Superior frontal: report whether it is ACA-dominant, MCA-dominant, or
  mixed; do NOT require MCA purity (that would be anatomically wrong).
- Same-donor Yang contrast viability: if hippocampus is PCA-dominant and
  superior frontal is ACA- or MCA-dominant (or ACA+MCA combined > 0.70 with
  PCA < 0.20), the within-donor territorial contrast remains usable, but
  the clean "PCA vs MCA" framing in the draft plan must be revised to
  "PCA (hippocampus) vs anterior-circulation superior frontal (ACA/MCA)".

Outputs
-------
results/feasibility/yang_roi_territory_overlap.csv
results/feasibility/yang_roi_territory_summary.md
results/feasibility/contrast_region_candidates.csv
"""
from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn import datasets, image

ROOT = Path(__file__).resolve().parents[2]
# Prefer 182×218×182 (classic FSL MNI152 1mm grid). GitHub NIfTIs ship a
# placeholder affine (origin at 0); we restore the standard FSL MNI152 1mm
# affine before any ROI resampling. Re-validate against NITRC MNI152 zip.
ATLAS_DIR = ROOT / "data" / "raw" / "arterial_atlas" / "Atlas_182"
FALLBACK_ATLAS_DIR = ROOT / "data" / "raw" / "arterial_atlas" / "Atlas"
OUT_DIR = ROOT / "results" / "feasibility"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Standard FSL MNI152 1mm (182×218×182), radiological (negative X)
FSL_MNI152_1MM_AFFINE = np.array(
    [
        [-1.0, 0.0, 0.0, 90.0],
        [0.0, 1.0, 0.0, -126.0],
        [0.0, 0.0, 1.0, -72.0],
        [0.0, 0.0, 0.0, 1.0],
    ],
    dtype=float,
)

# Level-2 intensities from ArterialAtlasLabels.txt
LEVEL2_MAJOR = {
    1: "ACA",  # ACAL
    2: "ACA",  # ACAR
    3: "MCA",  # MCAL
    4: "MCA",  # MCAR
    5: "PCA",  # PCAL
    6: "PCA",  # PCAR
    7: "VB",   # VBL
    8: "VB",   # VBR
    # 9,10 ventricles excluded
}

LEVEL1_TO_MAJOR = {
    1: "ACA", 2: "ACA", 3: "ACA", 4: "ACA",
    5: "MCA", 6: "MCA", 7: "MCA", 8: "MCA", 9: "MCA", 10: "MCA",
    11: "MCA", 12: "MCA", 13: "MCA", 14: "MCA", 15: "MCA", 16: "MCA",
    17: "PCA", 18: "PCA", 19: "PCA", 20: "PCA", 21: "PCA", 22: "PCA",
    23: "PCA", 24: "PCA",
    25: "VB", 26: "VB", 27: "VB", 28: "VB", 29: "VB", 30: "VB",
}

# Harvard-Oxford cortical / subcortical labels relevant to Yang + 2×2 design
# Names must match nilearn Harvard-Oxford label strings.
YANG_ROIS = {
    "hippocampus": {
        "atlas": "subcortical",
        "labels": ["Left Hippocampus", "Right Hippocampus"],
        "yang_role": "Yang primary: hippocampus samples",
        "a_priori_expectation": "PCA-dominant",
    },
    "superior_frontal_gyrus": {
        "atlas": "cortical",
        "labels": ["Superior Frontal Gyrus"],
        "yang_role": "Yang primary: superior frontal cortex samples",
        "a_priori_expectation": "ACA/MCA mixed (not pure MCA)",
    },
}

# Candidate ROIs for the Aim-2 2×2 confound-breaking design
CONTRAST_ROIS = {
    "occipital_pole": {
        "atlas": "cortical",
        "labels": ["Occipital Pole"],
        "contrast_role": "Territory-matched / region-differing vs hippocampus (both PCA-leaning)",
    },
    "lingual_gyrus": {
        "atlas": "cortical",
        "labels": ["Lingual Gyrus"],
        "contrast_role": "Additional PCA cortical ROI",
    },
    "frontal_medial_cortex": {
        "atlas": "cortical",
        "labels": ["Frontal Medial Cortex"],
        "contrast_role": "Region-matched cortex / ACA-leaning vs lateral frontal MCA",
    },
    "middle_frontal_gyrus": {
        "atlas": "cortical",
        "labels": ["Middle Frontal Gyrus"],
        "contrast_role": "MCA-leaning isocortex candidate (DLPFC-adjacent)",
    },
    "juxtapositional_lobule": {
        "atlas": "cortical",
        "labels": ["Juxtapositional Lobule Cortex (formerly Supplementary Motor Cortex)"],
        "contrast_role": "Medial frontal ACA-adjacent candidate",
    },
}


def shannon_entropy(probs: np.ndarray, eps: float = 1e-12) -> float:
    p = probs[probs > 0]
    if p.size == 0:
        return float("nan")
    p = p / p.sum()
    return float(max(0.0, -(p * np.log2(p + eps)).sum()))


def load_harvard_oxford():
    cort = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-2mm")
    sub = datasets.fetch_atlas_harvard_oxford("sub-maxprob-thr25-2mm")
    return {
        "cortical": (cort.maps, list(cort.labels)),
        "subcortical": (sub.maps, list(sub.labels)),
    }


def load_arterial_atlas(path: Path) -> nib.Nifti1Image:
    """Load atlas, squeeze trailing singleton, restore FSL MNI152 1mm affine."""
    img = nib.load(str(path))
    data = np.asanyarray(img.dataobj)
    if data.ndim == 4 and data.shape[-1] == 1:
        data = data[..., 0]
    if data.shape[:3] not in {(182, 218, 182), (181, 217, 181)}:
        raise ValueError(f"Unexpected atlas shape {data.shape} for {path}")
    affine = FSL_MNI152_1MM_AFFINE.copy()
    # 181-grid is a 1-voxel crop of the 182-grid; keep same voxel size/orientation.
    return nib.Nifti1Image(data.astype(np.int16), affine)


def roi_mask_from_labels(maps_img, label_names, wanted) -> np.ndarray:
    data = np.asanyarray(maps_img.dataobj)
    mask = np.zeros(data.shape, dtype=bool)
    for name in wanted:
        if name not in label_names:
            raise KeyError(f"Label '{name}' not in atlas. Available: {label_names}")
        idx = label_names.index(name)
        mask |= data == idx
    return mask


def resample_mask_to_atlas(mask_bool: np.ndarray, mask_affine, atlas_img) -> np.ndarray:
    mask_img = nib.Nifti1Image(mask_bool.astype(np.uint8), mask_affine)
    resampled = image.resample_to_img(
        mask_img,
        atlas_img,
        interpolation="nearest",
        force_resample=True,
        copy_header=True,
    )
    return np.asanyarray(resampled.dataobj).astype(bool)


def territory_vector(roi_mask: np.ndarray, territory_data: np.ndarray, mapping: dict) -> dict:
    vals = territory_data[roi_mask]
    vals = vals[vals > 0]
    # Drop ventricles (level2 9/10 or level1 31/32)
    vals = vals[~np.isin(vals, [9, 10, 31, 32])]
    n = vals.size
    if n == 0:
        return {
            "n_voxels": 0,
            "probs": {},
            "dominant": None,
            "p_dominant": float("nan"),
            "entropy_bits": float("nan"),
        }
    majors = np.array([mapping.get(int(v), "OTHER") for v in vals])
    majors = majors[majors != "OTHER"]
    n = majors.size
    counts = {k: int((majors == k).sum()) for k in ("ACA", "MCA", "PCA", "VB")}
    probs = {k: counts[k] / n for k in counts}
    dominant = max(probs, key=probs.get)
    p_vec = np.array([probs["ACA"], probs["MCA"], probs["PCA"], probs["VB"]], dtype=float)
    return {
        "n_voxels": int(n),
        "probs": probs,
        "dominant": dominant,
        "p_dominant": float(probs[dominant]),
        "entropy_bits": shannon_entropy(p_vec),
        "p_anterior": float(probs["ACA"] + probs["MCA"]),
    }


def score_roi(name: str, spec: dict, atlases, art_l1, art_l2) -> dict:
    maps_img, labels = atlases[spec["atlas"]]
    mask = roi_mask_from_labels(maps_img, labels, spec["labels"])
    mask_rs = resample_mask_to_atlas(mask, maps_img.affine, art_l2)
    # Restrict to voxels that have a territory label
    l2 = np.asanyarray(art_l2.dataobj)
    l1 = np.asanyarray(art_l1.dataobj)
    mask_rs = mask_rs & (l2 > 0) & (l2 < 9)
    v_l2 = territory_vector(mask_rs, l2, LEVEL2_MAJOR)
    v_l1 = territory_vector(mask_rs, l1, LEVEL1_TO_MAJOR)
    row = {
        "roi": name,
        "labels": ";".join(spec["labels"]),
        "n_voxels_in_atlas": v_l2["n_voxels"],
        "dominant_l2": v_l2["dominant"],
        "p_dominant_l2": v_l2["p_dominant"],
        "entropy_l2_bits": v_l2["entropy_bits"],
        "p_ACA": v_l2["probs"].get("ACA", 0.0) if v_l2["probs"] else np.nan,
        "p_MCA": v_l2["probs"].get("MCA", 0.0) if v_l2["probs"] else np.nan,
        "p_PCA": v_l2["probs"].get("PCA", 0.0) if v_l2["probs"] else np.nan,
        "p_VB": v_l2["probs"].get("VB", 0.0) if v_l2["probs"] else np.nan,
        "p_anterior_ACA_MCA": v_l2.get("p_anterior", np.nan),
        "dominant_l1_majorized": v_l1["dominant"],
        "p_dominant_l1_majorized": v_l1["p_dominant"],
        "entropy_l1_majorized_bits": v_l1["entropy_bits"],
        "high_specificity_pdom_ge_0.50": bool(v_l2["p_dominant"] >= 0.50) if v_l2["dominant"] else False,
        "high_specificity_H_le_1.0": bool(v_l2["entropy_bits"] <= 1.0) if np.isfinite(v_l2["entropy_bits"]) else False,
    }
    row.update({k: spec.get(k) for k in ("yang_role", "a_priori_expectation", "contrast_role") if k in spec})
    return row


def interpret_yang(rows: pd.DataFrame) -> dict:
    hip = rows.loc[rows["roi"] == "hippocampus"].iloc[0]
    sfg = rows.loc[rows["roi"] == "superior_frontal_gyrus"].iloc[0]
    hip_ok = hip["dominant_l2"] == "PCA" and hip["p_dominant_l2"] >= 0.50
    sfg_ant = sfg["p_anterior_ACA_MCA"] >= 0.70 and sfg["p_PCA"] < 0.20
    sfg_pure = sfg["dominant_l2"] in {"ACA", "MCA"} and sfg["p_dominant_l2"] >= 0.50
    verdict = {
        "hippocampus_pca_dominant": bool(hip_ok),
        "superior_frontal_anterior_circulation": bool(sfg_ant),
        "superior_frontal_single_territory_dominant": bool(sfg_pure),
        "within_donor_contrast_viable": bool(hip_ok and sfg_ant),
        "draft_plan_pca_vs_mca_framing_accurate": bool(
            hip_ok and sfg["dominant_l2"] == "MCA" and sfg["p_dominant_l2"] >= 0.50
        ),
        "recommended_yang_contrast_framing": (
            f"PCA (hippocampus) vs {sfg['dominant_l2']} (superior frontal)"
            if (hip_ok and sfg_ant)
            else "NOT YET VIABLE — inspect overlap table and consider alternate cortical ROIs / BMM maps"
        ),
    }
    return verdict


def write_summary(yang_df: pd.DataFrame, contrast_df: pd.DataFrame, verdict: dict) -> str:
    lines = [
        "# Yang → territory feasibility probe — results",
        "",
        "## Methods note",
        "",
        "- Territory vectors = voxel-fraction of Liu hard atlas labels inside Harvard-Oxford ROIs,",
        "  resampled (nearest) into ArterialAtlas space.",
        "- This is a proxy until NITRC `ProbArterialAtlas_BMM` / borderzone maps are ingested.",
        "- Yang cortex = **superior frontal cortex** (not generic cortex).",
        "",
        "## Pre-specified go criteria",
        "",
        "- Hippocampus: PCA-dominant with P≥0.50",
        "- Superior frontal: anterior-circulation (ACA+MCA) ≥0.70 and PCA <0.20",
        "- Within-donor Yang contrast viable if both hold",
        "",
        "## Verdict",
        "",
        "```json",
        json.dumps(verdict, indent=2),
        "```",
        "",
        "## Yang ROI territory overlap",
        "",
        yang_df.to_markdown(index=False),
        "",
        "## Aim-2 contrast ROI candidates",
        "",
        contrast_df.to_markdown(index=False),
        "",
        "## Interpretation guidance",
        "",
        "1. If hippocampus is PCA-dominant but superior frontal is mixed ACA/MCA with high entropy,",
        "   Yang still supports a **posterior vs anterior circulation** contrast, but **not** a clean",
        "   MCA territorial identity test. Use Winkler/Garcia regional metadata or additional atlases",
        "   for MCA-pure cortex.",
        "2. Prefer 2×2 contrasts from the candidate table where one ROI is PCA-dominant and another",
        "   is MCA- or ACA-dominant with P≥0.50 and H≤1.0.",
        "3. Replace hard-atlas fractions with BMM probabilities before locking the pre-registration",
        "   specificity threshold.",
        "",
    ]
    return "\n".join(lines)


def main():
    atlas_dir = ATLAS_DIR if (ATLAS_DIR / "ArterialAtlas.nii").exists() else FALLBACK_ATLAS_DIR
    art_l1 = load_arterial_atlas(atlas_dir / "ArterialAtlas.nii")
    art_l2 = load_arterial_atlas(atlas_dir / "ArterialAtlas_level2.nii")
    print(f"Using arterial atlas dir: {atlas_dir}")
    print(f"Restored affine:\n{art_l2.affine}")
    atlases = load_harvard_oxford()

    yang_rows = [score_roi(k, v, atlases, art_l1, art_l2) for k, v in YANG_ROIS.items()]
    contrast_rows = [score_roi(k, v, atlases, art_l1, art_l2) for k, v in CONTRAST_ROIS.items()]

    yang_df = pd.DataFrame(yang_rows)
    contrast_df = pd.DataFrame(contrast_rows)
    verdict = interpret_yang(yang_df)

    yang_path = OUT_DIR / "yang_roi_territory_overlap.csv"
    contrast_path = OUT_DIR / "contrast_region_candidates.csv"
    summary_path = OUT_DIR / "yang_roi_territory_summary.md"
    verdict_path = OUT_DIR / "yang_feasibility_verdict.json"

    yang_df.to_csv(yang_path, index=False)
    contrast_df.to_csv(contrast_path, index=False)
    summary_path.write_text(write_summary(yang_df, contrast_df, verdict), encoding="utf-8")
    verdict_path.write_text(json.dumps(verdict, indent=2), encoding="utf-8")

    print("Wrote:", yang_path)
    print("Wrote:", contrast_path)
    print("Wrote:", summary_path)
    print("Wrote:", verdict_path)
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
