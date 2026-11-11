"""
ATLAS Aim 4 stroke arm: lesion masks → Liu arterial territories (bilateral).

Primary atlas: ArterialAtlas_level2.nii (ACAL/ACAR/MCAL/MCAR/PCAL/PCAR/VBL/VBR)
collapsed to major territories ACA/MCA/PCA/VB.

Rationale: TerritoryVoxels_BMM in this NITRC release is not usable for bilateral
stroke geography (territory labels 1-4 fall almost entirely in one hemisphere;
contralesional lesions land in label 6). ArterialAtlas_level2 is the hierarchical
major-territory map with explicit left/right parcels (ArterialAtlasLabels.txt).

ATLAS R3.0 masks share a common MNI152NLin2009aSym 1 mm grid; level2 is
nearest-neighbor resampled once onto that grid.

Usage:
  python analysis/aim4/atlas_lesion_tv_pipeline.py
  python analysis/aim4/atlas_lesion_tv_pipeline.py --n 50
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn import image

ROOT = Path(__file__).resolve().parents[2]
ATLAS_ROOT = ROOT / "data/raw/atlas_stroke/extracted/ATLAS3_Training_Preprocessed"
LEVEL2_PATH = ROOT / (
    "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/"
    "Atlas_182_MNI152/ArterialAtlas_level2.nii"
)
OUT = ROOT / "results/aim4/atlas"

# ArterialAtlasLabels.txt level2 → major territory
LEVEL2_TO_MAJOR = {
    1: "ACA",  # ACAL
    2: "ACA",  # ACAR
    3: "MCA",  # MCAL
    4: "MCA",  # MCAR
    5: "PCA",  # PCAL
    6: "PCA",  # PCAR
    7: "VB",  # VBL
    8: "VB",  # VBR
}
MAJOR = ["ACA", "MCA", "PCA", "VB"]


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def find_masks() -> list[Path]:
    return sorted(ATLAS_ROOT.rglob("*label-lesion*mask*.nii.gz"))


def parse_ids(mask: Path) -> dict:
    sub = None
    site = None
    for i, p in enumerate(mask.parts):
        if p.startswith("sub-") and i > 0:
            site = mask.parts[i - 1]
            sub = p
            break
    meta = {}
    if sub:
        meta_p = mask.parent / f"{sub}_ses-1_metadata.csv"
        if meta_p.exists():
            try:
                m = pd.read_csv(meta_p).iloc[0].to_dict()
                meta = {
                    "days_post_stroke": m.get("DAYS_POST_STROKE"),
                    "chronicity": m.get("CHRONICITY"),
                    "atlas2_dataset": m.get("ATLAS2_DATASET"),
                    "session_id": m.get("SESSION_ID"),
                    "site_meta": m.get("SITE"),
                }
            except Exception:
                pass
    return {"subject_id": sub, "site": site, "mask_path": str(mask), **meta}


def collapse_level2(level2: np.ndarray) -> np.ndarray:
    """Map level2 labels to 1=ACA,2=MCA,3=PCA,4=VB; 0 elsewhere."""
    out = np.zeros(level2.shape, dtype=np.int16)
    code = {"ACA": 1, "MCA": 2, "PCA": 3, "VB": 4}
    for lab, maj in LEVEL2_TO_MAJOR.items():
        out[level2 == lab] = code[maj]
    return out


def tv_fractions(lesion: np.ndarray, tv: np.ndarray) -> dict:
    lesion = lesion.astype(bool)
    row = {"n_lesion_voxels": int(lesion.sum())}
    assigned = 0
    code = {"ACA": 1, "MCA": 2, "PCA": 3, "VB": 4}
    for name in MAJOR:
        m = tv == code[name]
        v = int((lesion & m).sum())
        assigned += v
        row[f"lesion_vox_{name}"] = v
        row[f"tv_vox_{name}"] = int(m.sum())
        row[f"lesion_frac_of_lesion_{name}"] = (
            v / row["n_lesion_voxels"] if row["n_lesion_voxels"] else 0.0
        )
        row[f"lesion_burden_in_{name}"] = (
            v / row[f"tv_vox_{name}"] if row[f"tv_vox_{name}"] else 0.0
        )
    row["lesion_vox_unassigned"] = row["n_lesion_voxels"] - assigned
    row["lesion_frac_assigned_tv"] = (
        assigned / row["n_lesion_voxels"] if row["n_lesion_voxels"] else 0.0
    )
    shares = {n: row[f"lesion_frac_of_lesion_{n}"] for n in MAJOR}
    if row["n_lesion_voxels"] > 0 and any(shares.values()):
        dom = max(shares, key=shares.get)
        row["dominant_territory"] = dom
        row["dominant_share"] = shares[dom]
    else:
        row["dominant_territory"] = None
        row["dominant_share"] = 0.0
    return row


def apply_qc(row: dict) -> dict:
    reasons = []
    flags = []
    n = int(row.get("n_lesion_voxels") or 0)
    frac = float(row.get("lesion_frac_assigned_tv") or 0.0)
    if n <= 0:
        reasons.append("FAIL_ZERO_LESION")
    if n > 0 and frac < 0.10:
        flags.append("FLAG_LOW_TV")
    # empty file / corrupt caught as error status upstream
    row["qc_reasons"] = ";".join(reasons)
    row["qc_flags"] = ";".join(flags)
    row["qc_primary"] = len(reasons) == 0
    row["qc_sensitivity_strict"] = len(reasons) == 0 and "FLAG_LOW_TV" not in flags
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=0, help="max subjects (0=all)")
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    out_csv = OUT / "atlas_lesion_tv_fractions.csv"
    masks = find_masks()
    if not masks:
        raise SystemExit(f"No lesion masks under {ATLAS_ROOT}")
    if args.n and args.n > 0:
        masks = masks[: args.n]
    log(f"ATLAS lesion masks: {len(masks)}")

    done = {}
    if args.resume and out_csv.exists():
        prev = pd.read_csv(out_csv)
        if "atlas_map" in prev.columns and (prev["atlas_map"] == "ArterialAtlas_level2").all():
            for _, r in prev.iterrows():
                if str(r.get("status", "")) == "ok":
                    done[str(r["subject_id"])] = r.to_dict()
            log(f"Resume: keeping {len(done)} completed")
        else:
            log("Resume skipped: prior CSV used different atlas_map")

    ref = nib.load(str(masks[0]))
    log(f"Resampling ArterialAtlas_level2 -> ATLAS grid {ref.shape} ...")
    l2_img = nib.load(str(LEVEL2_PATH))
    l2_r = image.resample_to_img(
        l2_img, ref, interpolation="nearest", force_resample=True, copy_header=True
    )
    l2 = np.rint(np.asanyarray(l2_r.dataobj, dtype=float)).astype(np.int32)
    if l2.ndim == 4:
        l2 = l2[..., 0]
    tv = collapse_level2(l2)
    tv_out = OUT / "ArterialAtlas_level2_major_on_ATLAS_MNI152NLin2009aSym.nii.gz"
    nib.save(nib.Nifti1Image(tv.astype(np.int16), ref.affine), str(tv_out))
    log(
        "Wrote "
        + tv_out.name
        + "; counts "
        + ", ".join(f"{n}={(tv == i).sum()}" for i, n in enumerate(MAJOR, 1))
    )

    rows = list(done.values())
    todo = [m for m in masks if parse_ids(m)["subject_id"] not in done]
    log(f"To process: {len(todo)}")

    for i, mask in enumerate(todo, 1):
        ids = parse_ids(mask)
        sid = ids["subject_id"]
        if i == 1 or i % 100 == 0 or i == len(todo):
            log(f"[{i}/{len(todo)}] {sid}")
        try:
            if Path(mask).stat().st_size < 100:
                raise RuntimeError("Empty or tiny file")
            img = nib.load(str(mask))
            if img.shape != ref.shape:
                img = image.resample_to_img(
                    img, ref, interpolation="nearest", force_resample=True, copy_header=True
                )
            lesion = np.asanyarray(img.dataobj) > 0
            if lesion.shape != tv.shape:
                raise RuntimeError(f"shape mismatch lesion{lesion.shape} tv{tv.shape}")
            fr = tv_fractions(lesion, tv)
            row = {
                **ids,
                **fr,
                "status": "ok",
                "space": "MNI152NLin2009aSym",
                "atlas_map": "ArterialAtlas_level2",
            }
            row = apply_qc(row)
        except Exception as e:
            row = {
                **ids,
                "status": "error",
                "error": str(e),
                "qc_primary": False,
                "atlas_map": "ArterialAtlas_level2",
            }
            log(f"  ERROR {sid}: {e}")
        rows.append(row)
        if i % 50 == 0 or i == len(todo):
            pd.DataFrame(rows).to_csv(out_csv, index=False)

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)

    ok = df[df["status"] == "ok"].copy()
    primary = ok[ok["qc_primary"] == True]  # noqa: E712
    strict = ok[ok["qc_sensitivity_strict"] == True]  # noqa: E712

    def share_means(d):
        return {
            t: float(d[f"lesion_frac_of_lesion_{t}"].mean()) if len(d) else None
            for t in MAJOR
        }

    dom_counts = (
        primary["dominant_territory"].value_counts(dropna=False).to_dict()
        if len(primary)
        else {}
    )
    summary = {
        "atlas_map": "ArterialAtlas_level2 (L/R collapsed to ACA/MCA/PCA/VB)",
        "atlas_map_rationale": (
            "TerritoryVoxels_BMM labels 1-4 are effectively unilateral in this release; "
            "ArterialAtlas_level2 provides bilateral major territories per ArterialAtlasLabels.txt."
        ),
        "n_masks_found": int(len(masks)),
        "n_rows": int(len(df)),
        "n_ok": int(len(ok)),
        "n_error": int((df["status"] != "ok").sum()),
        "n_primary_qc": int(len(primary)),
        "n_strict_qc": int(len(strict)),
        "n_fail_zero": int((ok["n_lesion_voxels"] <= 0).sum()) if len(ok) else 0,
        "n_flag_low_tv_in_primary": int(
            primary["qc_flags"].fillna("").astype(str).str.contains("FLAG_LOW_TV").sum()
        )
        if len(primary)
        else 0,
        "median_lesion_voxels_primary": float(primary["n_lesion_voxels"].median())
        if len(primary)
        else None,
        "median_frac_assigned_tv_primary": float(primary["lesion_frac_assigned_tv"].median())
        if len(primary)
        else None,
        "mean_lesion_share_primary": share_means(primary),
        "mean_lesion_share_strict": share_means(strict),
        "dominant_territory_counts_primary": {str(k): int(v) for k, v in dom_counts.items()},
        "tv_resampled_to": str(tv_out.relative_to(ROOT)),
        "interpretation_locks": [
            "ATLAS lesions are independent of Liu atlas construction set (prereg anti-circularity).",
            "Stroke geography uses ArterialAtlas_level2 (bilateral), not TerritoryVoxels_BMM.",
            "Descriptive stroke geography + MRA geometry null; molecular H3 still needs control-trained features.",
            "Primary QC: drop zero-lesion; FLAG_LOW_TV retained in primary, dropped in strict.",
        ],
    }
    (OUT / "atlas_lesion_tv_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    md = [
        "# ATLAS lesion x arterial territories (ArterialAtlas_level2)",
        "",
        f"**Atlas map:** {summary['atlas_map']}",
        f"**Rationale:** {summary['atlas_map_rationale']}",
        "",
        f"**n_ok / primary / strict:** {summary['n_ok']} / {summary['n_primary_qc']} / {summary['n_strict_qc']}",
        f"- Median lesion voxels (primary): {summary['median_lesion_voxels_primary']}",
        f"- Median assigned TV frac (primary): {summary['median_frac_assigned_tv_primary']}",
        f"- Mean lesion share ACA/MCA/PCA/VB (primary): {summary['mean_lesion_share_primary']}",
        f"- Dominant territory counts (primary): {summary['dominant_territory_counts_primary']}",
        "",
        "## Locked interpretation",
        "- Independent stroke outcome arm for Aim 4 (not Liu construction lesions).",
        "- Use bilateral level2 major territories; do not use TerritoryVoxels_BMM for ATLAS.",
        "- Geography is descriptive until nested with geometry (+/- molecular features).",
        "- Do not claim molecular H3 from lesion shares alone.",
        "",
        f"Table: `{out_csv.relative_to(ROOT)}`",
        "",
    ]
    (OUT / "atlas_lesion_tv_report.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
