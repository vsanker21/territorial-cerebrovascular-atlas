"""Full-cohort Aim-4 volume analysis using UCD/ADSP WMH (no spatial masks)."""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/aim4/adni"


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def main():
    vol = pd.read_csv(OUT / "adni_aim4_wmh_volume_baseline.csv")
    vol = vol[vol["TOTAL_WMH"].notna()].copy()
    # age from birth year if present
    if "PTDOBYY" in vol.columns and "image_date" in vol.columns:
        year = pd.to_datetime(vol["image_date"], errors="coerce").dt.year
        dob = pd.to_numeric(vol["PTDOBYY"], errors="coerce")
        vol["age_approx"] = year - dob
    summary = {
        "n": int(len(vol)),
        "n_subjects": int(vol["RID"].nunique()) if "RID" in vol.columns else int(len(vol)),
        "TOTAL_WMH_median": float(vol["TOTAL_WMH"].median()),
        "TOTAL_WMH_iqr": [
            float(vol["TOTAL_WMH"].quantile(0.25)),
            float(vol["TOTAL_WMH"].quantile(0.75)),
        ],
        "spatial_masks": False,
        "role": "non-stroke Aim4 volume endpoint; territory fractions via local FLAIR pilot",
    }
    # sex contrast if available
    if "PTGENDER" in vol.columns:
        g = vol.dropna(subset=["PTGENDER", "TOTAL_WMH"])
        levels = g["PTGENDER"].astype(str).unique().tolist()
        if len(levels) >= 2:
            a = g.loc[g["PTGENDER"].astype(str) == levels[0], "TOTAL_WMH"]
            b = g.loc[g["PTGENDER"].astype(str) == levels[1], "TOTAL_WMH"]
            u = stats.mannwhitneyu(a, b, alternative="two-sided")
            summary["wmh_by_sex"] = {
                levels[0]: float(a.median()),
                levels[1]: float(b.median()),
                "mannwhitney_p": float(u.pvalue),
            }
    if "age_approx" in vol.columns:
        sub = vol.dropna(subset=["age_approx", "TOTAL_WMH"])
        if len(sub) > 20:
            rho, p = stats.spearmanr(sub["age_approx"], sub["TOTAL_WMH"])
            summary["wmh_age_spearman"] = {"rho": float(rho), "p": float(p)}

    # join MRA geometry null reminder
    mra = ROOT / "results/aim4/mra_geometry_by_territory.csv"
    summary["mra_geometry_table"] = str(mra.relative_to(ROOT)) if mra.exists() else None

    (OUT / "adni_aim4_volume_analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = [
        "# ADNI Aim 4 — WMH volume analysis (existing data; no IDA masks)",
        "",
        f"**Verdict:** `VOLUME_ENDPOINT_READY__TERRITORY_VIA_LOCAL_FLAIR`",
        "",
        "## Cohort",
        f"- n visits/rows: {summary['n']} (subjects≈{summary['n_subjects']})",
        f"- TOTAL_WMH median (IQR): {summary['TOTAL_WMH_median']:.2f} "
        f"({summary['TOTAL_WMH_iqr'][0]:.2f}–{summary['TOTAL_WMH_iqr'][1]:.2f}) cc",
        "",
        "## Associations (descriptive)",
        f"- By sex: `{summary.get('wmh_by_sex')}`",
        f"- Age Spearman: `{summary.get('wmh_age_spearman')}`",
        "",
        "## Nested Aim 4 structure",
        "1. Geometry null: `results/aim4/mra_geometry_by_territory.csv` (already on disk)",
        "2. Non-stroke outcome (this table): subject-level UCD/ADSP WMH volumes",
        "3. Territory-resolved WMH: local FLAIR segmentation pilot "
        "(`adni_wmh_territory_pilot_report.md`) — scale after QC",
        "",
        "## Locked",
        "- No IDA spatial WMH masks.",
        "- Do not treat whole-brain TOTAL_WMH as territory-specific lesion probability.",
        "",
    ]
    (OUT / "adni_aim4_volume_analysis_report.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
