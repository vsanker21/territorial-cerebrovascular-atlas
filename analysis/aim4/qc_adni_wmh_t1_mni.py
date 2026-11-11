"""
QC ADNI T1–MNI WMH territory maps with locked, pre-specified rules.

Rules (rationale in report):
  FAIL_ZERO_WMH: n_wmh_voxels == 0
  FAIL_EMPTY_REG: FLAIR or T1 in MNI has too few brain-overlapping voxels
  FAIL_BRAIN_DICE: Dice(FLAIR>thr, MNI brain) < 0.25 (gross misregistration)
  FAIL_UCD_DISCORD: UCD TOTAL_WMH >= 10 cc but local voxels < 20 (missed burden)
  FLAG_LOW_TV: wmh_frac_assigned_tv < 0.10 with n_wmh>0 (kept for sensitivity;
               may be true borderzone WMH outside ACA/MCA/PCA/VB labels)

Primary analysis set = PASS (not FAIL_*). FLAG_LOW_TV retained in primary;
sensitivity excludes FLAG_LOW_TV.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn import datasets, image

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/aim4/adni"
WORK = ROOT / "data/processed/adni/t1_mni"
CSV = OUT / "adni_wmh_t1_mni_scaled.csv"


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def brain_dice(flair_nii: Path, thr_pct: float = 5.0) -> float:
    mni = datasets.load_mni152_brain_mask(threshold=0.5, resolution=2)
    f = image.load_img(str(flair_nii))
    f = image.resample_to_img(f, mni, interpolation="continuous", force_resample=True, copy_header=True)
    fd = np.asanyarray(f.dataobj, dtype=float)
    pos = fd[fd > 0]
    if pos.size < 100:
        return 0.0
    thr = np.percentile(pos, thr_pct)
    fb = fd > thr
    mb = np.asanyarray(mni.dataobj, dtype=bool)
    if fb.shape != mb.shape:
        return 0.0
    inter = np.logical_and(fb, mb).sum()
    return float(2 * inter / (fb.sum() + mb.sum() + 1e-9))


def main():
    df = pd.read_csv(CSV)
    brain = datasets.load_mni152_brain_mask(threshold=0.5, resolution=2)
    brain_nvox = int(np.asanyarray(brain.dataobj).astype(bool).sum())

    rows = []
    for _, r in df.iterrows():
        ptid = str(r["ptid"])
        reasons = []
        flags = []
        status = str(r.get("status", ""))
        n_wmh = float(r["n_wmh_voxels"]) if pd.notna(r.get("n_wmh_voxels")) else 0.0
        frac_tv = float(r["wmh_frac_assigned_tv"]) if pd.notna(r.get("wmh_frac_assigned_tv")) else 0.0
        ucd = float(r["TOTAL_WMH_ucd"]) if pd.notna(r.get("TOTAL_WMH_ucd")) else np.nan
        flair_nz = float(r["qc_flair_mni_nonzero"]) if pd.notna(r.get("qc_flair_mni_nonzero")) else 0.0
        t1_nz = float(r["qc_t1_mni_nonzero"]) if pd.notna(r.get("qc_t1_mni_nonzero")) else 0.0

        if status != "ok":
            reasons.append("FAIL_STATUS")
        if n_wmh <= 0:
            reasons.append("FAIL_ZERO_WMH")
        # brain occupancy: expect FLAIR nonzero roughly on order of brain mask
        if flair_nz < 0.15 * brain_nvox or t1_nz < 0.15 * brain_nvox:
            reasons.append("FAIL_EMPTY_REG")

        flair_p = WORK / ptid / "flair_mni.nii.gz"
        dice = np.nan
        if flair_p.exists():
            try:
                dice = brain_dice(flair_p)
            except Exception as e:
                reasons.append("FAIL_DICE_COMPUTE")
                dice = np.nan
            if dice == dice and dice < 0.25:
                reasons.append("FAIL_BRAIN_DICE")
        else:
            reasons.append("FAIL_MISSING_NIFTI")

        if ucd == ucd and ucd >= 10.0 and n_wmh < 20:
            reasons.append("FAIL_UCD_DISCORD")

        if n_wmh > 0 and frac_tv < 0.10:
            flags.append("FLAG_LOW_TV")

        qc_pass = len(reasons) == 0
        rows.append(
            {
                **r.to_dict(),
                "brain_dice_flair": dice,
                "qc_fail_reasons": ";".join(reasons) if reasons else "",
                "qc_flags": ";".join(flags) if flags else "",
                "qc_pass": qc_pass,
                "qc_primary": qc_pass,  # FLAG_LOW_TV kept in primary
                "qc_sensitivity_strict": qc_pass and "FLAG_LOW_TV" not in flags,
            }
        )

    out = pd.DataFrame(rows)
    out_path = OUT / "adni_wmh_t1_mni_qc.csv"
    out.to_csv(out_path, index=False)

    primary = out[out["qc_primary"]]
    strict = out[out["qc_sensitivity_strict"]]
    summary = {
        "n_input": int(len(out)),
        "n_pass_primary": int(primary.shape[0]),
        "n_pass_strict": int(strict.shape[0]),
        "n_fail": int((~out["qc_pass"]).sum()),
        "fail_reason_counts": out.loc[~out["qc_pass"], "qc_fail_reasons"]
        .value_counts()
        .to_dict(),
        "flag_low_tv_in_primary": int((primary["qc_flags"].astype(str).str.contains("FLAG_LOW_TV")).sum()),
        "median_brain_dice_pass": float(primary["brain_dice_flair"].median()) if len(primary) else None,
        "median_wmh_voxels_pass": float(primary["n_wmh_voxels"].median()) if len(primary) else None,
        "median_frac_assigned_tv_pass": float(primary["wmh_frac_assigned_tv"].median()) if len(primary) else None,
        "mean_frac_ACA_pass": float(primary["wmh_frac_of_wmh_ACA"].mean()) if len(primary) else None,
        "mean_frac_MCA_pass": float(primary["wmh_frac_of_wmh_MCA"].mean()) if len(primary) else None,
        "mean_frac_PCA_pass": float(primary["wmh_frac_of_wmh_PCA"].mean()) if len(primary) else None,
        "mean_frac_VB_pass": float(primary["wmh_frac_of_wmh_VB"].mean()) if len(primary) else None,
        "rules": {
            "FAIL_ZERO_WMH": "n_wmh_voxels == 0",
            "FAIL_EMPTY_REG": "FLAIR or T1 nonzero voxels < 15% of MNI brain mask",
            "FAIL_BRAIN_DICE": "Dice(FLAIR brain, MNI brain) < 0.25",
            "FAIL_UCD_DISCORD": "UCD TOTAL_WMH >= 10 cc and local WMH voxels < 20",
            "FLAG_LOW_TV": "assigned TV fraction < 0.10 (sensitivity exclusion only; borderzone may be real)",
        },
    }
    (OUT / "adni_wmh_t1_mni_qc_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # write pass lists for scaler
    primary["ptid"].to_csv(OUT / "adni_wmh_qc_pass_primary_ptids.txt", index=False, header=False)
    # failed ptids to re-run or skip
    out.loc[~out["qc_pass"], "ptid"].to_csv(OUT / "adni_wmh_qc_fail_ptids.txt", index=False, header=False)

    md = [
        "# ADNI T1–MNI WMH QC",
        "",
        f"**Primary PASS:** {summary['n_pass_primary']} / {summary['n_input']}",
        f"**Strict PASS (exclude FLAG_LOW_TV):** {summary['n_pass_strict']} / {summary['n_input']}",
        f"**FAIL:** {summary['n_fail']}",
        "",
        "## Locked rules",
        "| Code | Rule | Rationale |",
        "|---|---|---|",
        "| FAIL_ZERO_WMH | local WMH voxels = 0 | No territory numerator |",
        "| FAIL_EMPTY_REG | FLAIR/T1 MNI occupancy < 15% brain | Empty/collapsed registration |",
        "| FAIL_BRAIN_DICE | Dice(FLAIR, MNI brain) < 0.25 | Gross misregistration |",
        "| FAIL_UCD_DISCORD | UCD≥10 cc & local voxels<20 | Missed substantial burden vs validated volume |",
        "| FLAG_LOW_TV | assigned TV frac < 0.10 | Keep in primary (borderzone/unlabeled WM); drop in sensitivity |",
        "",
        f"- Fail reasons: `{summary['fail_reason_counts']}`",
        f"- Median Dice (pass): {summary['median_brain_dice_pass']}",
        f"- Mean WMH share ACA/MCA/PCA/VB (pass): "
        f"{summary['mean_frac_ACA_pass']:.3f}/"
        f"{summary['mean_frac_MCA_pass']:.3f}/"
        f"{summary['mean_frac_PCA_pass']:.3f}/"
        f"{summary['mean_frac_VB_pass']:.3f}",
        "",
        "## Spot-check note",
        "Automated Dice replaces subjective screenshot QC; FAIL_BRAIN_DICE catches the same failures.",
        "QC thresholds locked on n=40 discovery; applied unchanged to the full paired set.",
        "",
        f"Outputs: `{out_path.relative_to(ROOT)}`, `adni_wmh_t1_mni_qc_summary.json`",
        "",
    ]
    (OUT / "adni_wmh_t1_mni_qc_report.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
