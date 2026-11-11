"""
Attach TerritoryVoxels_BMM-based territorial annotations to Yang samples.

Scientific rules (from results/feasibility/bmm_confirmatory_interpretation.md):
  - Labels 1–4 = ACA, MCA, PCA, VB (empirical mapping for this NITRC release)
  - Labels 5–6 excluded from primary probability mass
  - Report assignment_fraction = (# voxels in 1–4) / (# ROI voxels)
  - Dominant territory among assigned; gate H1_primary if
      assignment_fraction >= 0.20 AND p_dominant >= 0.50
  - Superior frontal with mixed ACA/MCA → Aim3_borderzone stream
  - Hard-atlas vectors retained as secondary columns (not confirmatory alone)
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
    LEVEL2_MAJOR,
    load_arterial_atlas,
    resample_mask_to_atlas,
    roi_mask_from_labels,
    shannon_entropy,
    territory_vector,
)

META = ROOT / "data/processed/yang/yang_sample_metadata.csv"
NITRC_MNI = (
    ROOT / "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/Atlas_182_MNI152"
)
ATLAS_DIR = NITRC_MNI if (NITRC_MNI / "ArterialAtlas_level2.nii").exists() else (
    ROOT / "data/raw/arterial_atlas/Atlas_182"
)
TERR_PATH = NITRC_MNI / "TerritoryVoxels_BMM.nii"
BORDER_PATH = NITRC_MNI / "BorderZone_ProbAve.nii"
OUT = ROOT / "data/processed/yang/yang_sample_territory_voxels.csv"
OUT_JSON = ROOT / "data/processed/yang/yang_region_territory_voxels_lookup.json"
OUT_MD = ROOT / "results/feasibility/yang_territory_voxels_annotation.md"

TV_LABELS = {1: "ACA", 2: "MCA", 3: "PCA", 4: "VB"}
REGION_TO_HO = {
    "hippocampus": {"atlas": "subcortical", "labels": ["Left Hippocampus", "Right Hippocampus"]},
    "superior_frontal_cortex": {"atlas": "cortical", "labels": ["Superior Frontal Gyrus"]},
}


def load_ho():
    cort = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-2mm")
    sub = datasets.fetch_atlas_harvard_oxford("sub-maxprob-thr25-2mm")
    return {
        "cortical": (cort.maps, list(cort.labels)),
        "subcortical": (sub.maps, list(sub.labels)),
    }


def tv_vector(mask: np.ndarray, tv: np.ndarray) -> dict:
    roi_n = int(mask.sum())
    vals = tv[mask]
    assigned = vals[(vals >= 1) & (vals <= 4)]
    n_asg = int(assigned.size)
    probs = {name: 0.0 for name in ("ACA", "MCA", "PCA", "VB")}
    if n_asg:
        for lab, name in TV_LABELS.items():
            probs[name] = float((assigned == lab).sum() / n_asg)
    p = np.array([probs["ACA"], probs["MCA"], probs["PCA"], probs["VB"]], dtype=float)
    if p.sum() > 0:
        dom_i = int(p.argmax())
        dominant = ["ACA", "MCA", "PCA", "VB"][dom_i]
        p_dom = float(p[dom_i])
        ent = shannon_entropy(p)
    else:
        dominant, p_dom, ent = None, float("nan"), float("nan")
    # mixed anterior circulation (ACA+MCA) without PCA/VB dominance
    anterior = probs["ACA"] + probs["MCA"]
    is_borderish = bool(
        n_asg > 0
        and anterior >= 0.70
        and max(probs["ACA"], probs["MCA"]) < 0.70
        and probs["PCA"] < 0.25
    )
    return {
        "n_roi_voxels": roi_n,
        "n_assigned_1to4": n_asg,
        "assignment_fraction": float(n_asg / roi_n) if roi_n else float("nan"),
        "p_ACA": probs["ACA"],
        "p_MCA": probs["MCA"],
        "p_PCA": probs["PCA"],
        "p_VB": probs["VB"],
        "dominant": dominant,
        "p_dominant": p_dom,
        "entropy_bits": ent,
        "aca_mca_borderzone_flag": is_borderish,
        "n_label5": int((vals == 5).sum()),
        "n_label6": int((vals == 6).sum()),
    }


def border_stats(mask: np.ndarray, border_img) -> dict:
    if border_img is None:
        return {"borderzone_defined_fraction": float("nan"), "borderzone_mean_where_defined": float("nan")}
    b = np.asanyarray(border_img.dataobj, dtype=float)
    if b.ndim == 4:
        b = b[..., 0]
    sel = mask & np.isfinite(b)
    frac_def = float(sel.sum() / mask.sum()) if mask.sum() else float("nan")
    mean_def = float(np.nanmean(b[sel])) if sel.sum() else float("nan")
    return {
        "borderzone_defined_fraction": frac_def,
        "borderzone_mean_where_defined": mean_def,
    }


def main():
    if not META.exists():
        raise SystemExit(f"Missing {META}")
    if not TERR_PATH.exists():
        raise SystemExit(f"Missing TerritoryVoxels at {TERR_PATH}")

    samples = pd.read_csv(META)
    art_l2 = load_arterial_atlas(ATLAS_DIR / "ArterialAtlas_level2.nii")
    tv_img = nib.load(str(TERR_PATH))
    # Ensure TV on same grid as art_l2 (NITRC MNI152 should match Atlas_182)
    if tv_img.shape[:3] != art_l2.shape[:3]:
        tv_img = image.resample_to_img(tv_img, art_l2, interpolation="nearest")
    # NITRC int16 maps can carry scl_slope → float labels like 1.00004578.
    # Round to nearest integer before label logic.
    tv = np.rint(np.asanyarray(tv_img.dataobj, dtype=float)).astype(np.int32)
    if tv.ndim == 4:
        tv = tv[..., 0]

    border_img = None
    if BORDER_PATH.exists():
        border_img = nib.load(str(BORDER_PATH))
        if border_img.shape[:3] != art_l2.shape[:3]:
            border_img = image.resample_to_img(border_img, art_l2, interpolation="continuous")

    atlases = load_ho()
    lookup = {}
    for region, spec in REGION_TO_HO.items():
        maps_img, labels = atlases[spec["atlas"]]
        mask = roi_mask_from_labels(maps_img, labels, spec["labels"])
        mask_rs = resample_mask_to_atlas(mask, maps_img.affine, art_l2)
        hard = territory_vector(mask_rs & (np.asanyarray(art_l2.dataobj) > 0) & (np.asanyarray(art_l2.dataobj) < 9), np.asanyarray(art_l2.dataobj), LEVEL2_MAJOR)
        tv_stats = tv_vector(mask_rs, tv)
        bstats = border_stats(mask_rs, border_img)
        h1 = bool(
            tv_stats["assignment_fraction"] >= 0.20
            and tv_stats["p_dominant"] >= 0.50
            and not tv_stats["aca_mca_borderzone_flag"]
        )
        if tv_stats["aca_mca_borderzone_flag"]:
            stream = "Aim3_borderzone"
        elif h1:
            stream = "H1_primary"
        else:
            stream = "exploratory_low_assignment"
        lookup[region] = {
            **{f"tv_{k}": v for k, v in tv_stats.items()},
            **bstats,
            "hard_p_ACA": hard["probs"].get("ACA"),
            "hard_p_MCA": hard["probs"].get("MCA"),
            "hard_p_PCA": hard["probs"].get("PCA"),
            "hard_p_VB": hard["probs"].get("VB"),
            "hard_dominant": hard["dominant"],
            "analysis_stream": stream,
            "method": "TerritoryVoxels_BMM_labels_1to4",
            "label_map": TV_LABELS,
        }

    rows = []
    for _, s in samples.iterrows():
        region = s["region"]
        base = lookup.get(region, {})
        rows.append(
            {
                "geo_accession": s.get("geo_accession"),
                "title": s.get("title"),
                "donor_id": s.get("donor_id"),
                "disease": s.get("disease"),
                "region": region,
                "territory_method": base.get("method"),
                "analysis_stream": base.get("analysis_stream"),
                "tv_dominant": base.get("tv_dominant"),
                "tv_p_dominant": base.get("tv_p_dominant"),
                "tv_assignment_fraction": base.get("tv_assignment_fraction"),
                "tv_entropy_bits": base.get("tv_entropy_bits"),
                "tv_p_ACA": base.get("tv_p_ACA"),
                "tv_p_MCA": base.get("tv_p_MCA"),
                "tv_p_PCA": base.get("tv_p_PCA"),
                "tv_p_VB": base.get("tv_p_VB"),
                "tv_aca_mca_borderzone_flag": base.get("tv_aca_mca_borderzone_flag"),
                "borderzone_defined_fraction": base.get("borderzone_defined_fraction"),
                "borderzone_mean_where_defined": base.get("borderzone_mean_where_defined"),
                "hard_dominant": base.get("hard_dominant"),
                "hard_p_ACA": base.get("hard_p_ACA"),
                "hard_p_MCA": base.get("hard_p_MCA"),
                "hard_p_PCA": base.get("hard_p_PCA"),
                "hard_p_VB": base.get("hard_p_VB"),
            }
        )

    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    # JSON-safe
    def _safe(o):
        if isinstance(o, dict):
            return {k: _safe(v) for k, v in o.items()}
        if isinstance(o, float) and (np.isnan(o) or np.isinf(o)):
            return None
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        return o

    OUT_JSON.write_text(json.dumps(_safe(lookup), indent=2), encoding="utf-8")

    lines = [
        "# Yang TerritoryVoxels annotation",
        "",
        f"Method: `{lookup[next(iter(lookup))]['method']}`",
        "",
        "| Region | Stream | Dom | P_dom | Assign frac | ACA | MCA | PCA | VB | Border flag |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for region, e in lookup.items():
        lines.append(
            f"| {region} | {e['analysis_stream']} | {e['tv_dominant']} | "
            f"{e['tv_p_dominant']:.3f} | {e['tv_assignment_fraction']:.3f} | "
            f"{e['tv_p_ACA']:.3f} | {e['tv_p_MCA']:.3f} | {e['tv_p_PCA']:.3f} | {e['tv_p_VB']:.3f} | "
            f"{e['tv_aca_mca_borderzone_flag']} |"
        )
    lines += [
        "",
        f"Wrote `{OUT.relative_to(ROOT).as_posix()}`",
        "",
        "Contrast lock: hippocampus (PCA-dominant) vs superior_frontal_cortex (ACA/MCA borderzone).",
    ]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(df.to_string(index=False))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
