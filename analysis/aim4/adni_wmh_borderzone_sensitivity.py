"""ADNI WMH × BorderZone_ProbAve sensitivity (primary QC set).

Resamples BorderZone_ProbAve into each subject's WMH MNI grid (2 mm).
Reports fraction of local WMH voxels in borderzone (prob > thr).
"""
from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nibabel.processing import resample_from_to
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
NEST = ROOT / "results/aim4/adni/aim4_nest_wmh_mra_subject_table.csv"
MASK_ROOT = ROOT / "data/processed/adni/t1_mni"
BZ = ROOT / (
    "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/"
    "Atlas_MNI152/BorderZone_ProbAve.nii"
)
OUT = ROOT / "results/aim4/adni"
THRESHOLDS = (0.25, 0.5, 0.75)


def load_bz_3d() -> nib.Nifti1Image:
    img = nib.load(str(BZ))
    data = np.asanyarray(img.dataobj, dtype=np.float32)
    if data.ndim == 4:
        data = data[..., 0]
    return nib.Nifti1Image(data, img.affine)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    nest = pd.read_csv(NEST)
    bz_img = load_bz_3d()
    ref = None
    bz_r = None
    rows = []
    for ptid in nest["ptid"].astype(str):
        mask_p = MASK_ROOT / ptid / "wmh_mni.nii.gz"
        if not mask_p.exists():
            continue
        wmh_img = nib.load(str(mask_p))
        if ref is None or wmh_img.shape != ref.shape or not np.allclose(wmh_img.affine, ref.affine):
            ref = wmh_img
            bz_r = np.asanyarray(
                resample_from_to(bz_img, (wmh_img.shape, wmh_img.affine), order=1).dataobj
            )
        wmh = np.asanyarray(wmh_img.dataobj) > 0
        n = int(wmh.sum())
        if n == 0:
            continue
        vals = bz_r[wmh]
        row = {
            "ptid": ptid,
            "n_wmh": n,
            "mean_bz_prob_in_wmh": float(vals.mean()),
            "wmh_frac_assigned_tv": float(
                nest.loc[nest["ptid"].astype(str) == ptid, "wmh_frac_assigned_tv"].iloc[0]
            ),
        }
        for thr in THRESHOLDS:
            row[f"frac_wmh_bz_gt_{thr}"] = float((vals > thr).mean())
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "aim4_wmh_borderzone_subject.csv", index=False)
    summary = {
        "n_subjects": int(len(df)),
        "mean_frac_wmh_bz_gt_0.5": float(df["frac_wmh_bz_gt_0.5"].mean()) if len(df) else None,
        "mean_frac_wmh_bz_gt_0.25": float(df["frac_wmh_bz_gt_0.25"].mean()) if len(df) else None,
        "mean_bz_prob_in_wmh": float(df["mean_bz_prob_in_wmh"].mean()) if len(df) else None,
        "interpretation_locks": [
            "Sensitivity only; primary lock remains TerritoryVoxels major shares.",
            "High borderzone overlap can help explain unlabeled / FLAG_LOW_TV WMH.",
        ],
    }
    if len(df):
        unassigned = 1.0 - df["wmh_frac_assigned_tv"]
        rho, p = stats.spearmanr(df["frac_wmh_bz_gt_0.5"], unassigned)
        summary["spearman_bz_gt0.5_vs_unassigned"] = {"rho": float(rho), "p": float(p)}
    (OUT / "aim4_wmh_borderzone_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = [
        "# ADNI WMH × BorderZone_ProbAve sensitivity",
        "",
        f"n={summary['n_subjects']}; mean frac WMH with BZ prob>0.5 = {summary.get('mean_frac_wmh_bz_gt_0.5')}",
        f"mean frac WMH with BZ prob>0.25 = {summary.get('mean_frac_wmh_bz_gt_0.25')}",
        f"mean BZ prob in WMH voxels = {summary.get('mean_bz_prob_in_wmh')}",
        f"Spearman BZ>0.5 vs unassigned: {summary.get('spearman_bz_gt0.5_vs_unassigned')}",
        "",
        "Locked: sensitivity only; does not change MCA-dominant major-territory geography lock.",
        "",
    ]
    (OUT / "aim4_wmh_borderzone_report.md").write_text("\n".join(md), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
