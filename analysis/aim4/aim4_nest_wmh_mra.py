"""
Aim 4 nest: ADNI territory WMH fractions + MRA geometry null.

Design (anti-circularity):
  - Outcome: territorial WMH burden (local FLAIR seg in MNI ∩ TerritoryVoxels)
  - Geometry null: Mouches/Forkert MRA vessel probability / radius by territory
    (atlas-level constants — same for all subjects; tests composition vs geometry baseline)
  - Volume endpoint: UCD TOTAL_WMH (validated) as parallel non-spatial outcome
  - Analysis set: QC primary PASS from adni_wmh_t1_mni_qc.csv

Nested models (subject-level):
  M0: demographics (age proxy, sex) -> TOTAL_WMH (volume arm)
  M1: geometry features only (atlas TV geometry replicated per subject) — null for territory
  M2: M1 + territorial WMH burden vector (ACA/MCA/PCA/VB)
  Report: Spearman of local vs UCD; territory share summaries; OLS R^2 nest for TOTAL_WMH
          and for MCA burden (most powered).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/aim4/adni"
MRA = ROOT / "results/aim4/mra_geometry_by_territory.csv"
QC = OUT / "adni_wmh_t1_mni_qc.csv"
INV = OUT / "adni_aim4_cohort_inventory.csv"


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def r2_ols(X, y):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(X).all(axis=1) & np.isfinite(y)
    if mask.sum() < 8:
        return np.nan, int(mask.sum())
    Xs = StandardScaler().fit_transform(X[mask])
    reg = LinearRegression().fit(Xs, y[mask])
    return float(reg.score(Xs, y[mask])), int(mask.sum())


def _run_nest(analysis: pd.DataFrame, label: str, n_primary: int, n_strict: int, geo: dict):
    """Run nest models on `analysis` dataframe; write label-specific outputs."""
    for k, v in geo.items():
        analysis[k] = v

    burden_cols = [f"wmh_burden_in_{t}" for t in ["ACA", "MCA", "PCA", "VB"]]
    share_cols = [f"wmh_frac_of_wmh_{t}" for t in ["ACA", "MCA", "PCA", "VB"]]
    geo_feat = [
        "mra_vessel_probability_ACA",
        "mra_vessel_probability_MCA",
        "mra_vessel_probability_PCA",
        "mra_vessel_probability_VB",
        "mra_vessel_radius_ACA",
        "mra_vessel_radius_MCA",
        "mra_vessel_radius_PCA",
        "mra_vessel_radius_VB",
    ]
    geo_feat = [c for c in geo_feat if c in analysis.columns]

    y_vol = analysis["TOTAL_WMH_ucd"].astype(float).to_numpy()
    X0 = analysis[["age_approx", "sex_f"]].to_numpy()
    r2_0, n0 = r2_ols(X0, y_vol)
    X1 = analysis[geo_feat].to_numpy()
    r2_1, n1 = r2_ols(X1, y_vol)
    X2 = analysis[geo_feat + burden_cols].to_numpy()
    r2_2, n2 = r2_ols(X2, y_vol)
    X3 = analysis[["age_approx", "sex_f"] + burden_cols].to_numpy()
    r2_3, n3 = r2_ols(X3, y_vol)

    mca = analysis["wmh_burden_in_MCA"].astype(float)
    aca = analysis["wmh_burden_in_ACA"].astype(float)
    wil = stats.wilcoxon(mca, aca, zero_method="wilcox", alternative="two-sided") if len(analysis) >= 10 else None
    rho_ucd, p_ucd = stats.spearmanr(analysis["TOTAL_WMH_ucd"], analysis["n_wmh_voxels"])

    suffix = "" if label == "primary" else f"_{label}"
    keep = [
        "ptid",
        "RID",
        "TOTAL_WMH_ucd",
        "n_wmh_voxels",
        "wmh_frac_assigned_tv",
        "brain_dice_flair",
        "age_approx",
        "sex_f",
        "qc_flags",
    ] + burden_cols + share_cols + geo_feat
    keep = [c for c in keep if c in analysis.columns]
    nest_path = OUT / f"aim4_nest_wmh_mra_subject_table{suffix}.csv"
    analysis[keep].to_csv(nest_path, index=False)

    terr_summary = []
    for t in ["ACA", "MCA", "PCA", "VB"]:
        terr_summary.append(
            {
                "territory": t,
                "mean_wmh_burden": float(analysis[f"wmh_burden_in_{t}"].mean()),
                "median_wmh_burden": float(analysis[f"wmh_burden_in_{t}"].median()),
                "mean_wmh_share": float(analysis[f"wmh_frac_of_wmh_{t}"].mean()),
                "mra_vessel_probability_mean": geo.get(f"mra_vessel_probability_{t}"),
                "mra_vessel_radius_mean": geo.get(f"mra_vessel_radius_{t}"),
            }
        )
    pd.DataFrame(terr_summary).to_csv(OUT / f"aim4_nest_wmh_mra_territory_summary{suffix}.csv", index=False)

    summary = {
        "analysis_set": label,
        "n_analysis": int(len(analysis)),
        "n_primary_qc": n_primary,
        "n_strict_qc": n_strict,
        "spearman_local_vs_ucd": {"rho": float(rho_ucd), "p": float(p_ucd)},
        "ols_r2_TOTAL_WMH": {
            "M0_demographics": {"r2": r2_0, "n": n0},
            "M1_mra_geometry_null": {"r2": r2_1, "n": n1, "note": "atlas constants identical across subjects"},
            "M2_geometry_plus_tv_burden": {"r2": r2_2, "n": n2},
            "M3_demographics_plus_tv_burden": {"r2": r2_3, "n": n3},
            "delta_M3_minus_M0": (r2_3 - r2_0) if (r2_3 == r2_3 and r2_0 == r2_0) else None,
        },
        "wilcoxon_MCA_vs_ACA_burden": {
            "stat": float(wil.statistic) if wil else None,
            "p": float(wil.pvalue) if wil else None,
        },
        "territory_means": terr_summary,
        "interpretation_locks": [
            "MRA geometry is atlas-level (not subject-specific); M1 R^2 near 0 is expected.",
            "M3 vs M0 (territorial WMH burden predicting UCD TOTAL_WMH) is a partition/concordance check, not H3 molecular proof.",
            "Local WMH is exploratory vs UCD clinical volumes; report Spearman concordance.",
            "Descriptive geography: MCA-dominant WMH share in this QC-pass ADNI set.",
            "ATLAS stroke lesions remain a separate arm (passphrase-gated).",
            "Primary lock is qc_primary at z=1.5; strict excludes FLAG_LOW_TV (sensitivity only).",
        ],
    }
    (OUT / f"aim4_nest_wmh_mra_summary{suffix}.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    md = [
        f"# Aim 4 nest — ADNI territorial WMH + MRA geometry ({label})",
        "",
        f"**Analysis set:** QC {label} n={len(analysis)} "
        f"(registry: primary n={n_primary}, strict n={n_strict})",
        "",
        "## Outcomes",
        f"- Local vs UCD TOTAL_WMH Spearman ρ={rho_ucd:.3f}, p={p_ucd:.2e}",
        f"- Wilcoxon MCA vs ACA WMH burden: p="
        f"{(wil.pvalue if wil else float('nan')):.3g}",
        "",
        "## Nested OLS R² predicting UCD TOTAL_WMH",
        f"| Model | R² | n |",
        f"|---|---:|---:|",
        f"| M0 demographics (age, sex) | {r2_0:.3f} | {n0} |",
        f"| M1 MRA geometry null (atlas constants) | {r2_1:.3f} | {n1} |",
        f"| M2 geometry + TV WMH burden | {r2_2:.3f} | {n2} |",
        f"| M3 demographics + TV WMH burden | {r2_3:.3f} | {n3} |",
        f"| **Δ (M3 − M0)** | **{(r2_3 - r2_0) if r2_3==r2_3 and r2_0==r2_0 else float('nan'):.3f}** | |",
        "",
        f"## Territory means ({label})",
        pd.DataFrame(terr_summary).to_string(index=False),
        "",
        "## Locked interpretation",
        "- Atlas MRA geometry is a *shared* null covariate (not per-subject vascular anatomy); M1 R²≈0 is expected.",
        "- M3 predicting UCD TOTAL_WMH from territorial burdens of the *same* FLAIR WMH map is a "
        "**partition / concordance check** (parts vs clinical total), **not** H3 molecular proof.",
        "- Aim 4 molecular claim still requires control-trained vascular transcriptomic features "
        "predicting independent territorial WMH or ATLAS lesions.",
        "- Descriptive result: ADNI WMH mass is predominantly MCA-territory share in this QC set.",
        "- Do not claim stroke-lesion vulnerability until ATLAS decrypt + independent lesion maps.",
        "",
        f"Tables: `{nest_path.relative_to(ROOT)}`, `aim4_nest_wmh_mra_territory_summary{suffix}.csv`",
        "",
    ]
    report = OUT / ("aim4_nest_wmh_mra_report.md" if label == "primary" else f"aim4_nest_wmh_mra_report_{label}.md")
    report.write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(summary, indent=2, default=str))
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--analysis-set",
        choices=["primary", "strict", "both"],
        default="both",
        help="primary = qc_primary; strict = exclude FLAG_LOW_TV; both writes primary+strict",
    )
    args = ap.parse_args()

    if not QC.exists():
        raise SystemExit(f"Missing {QC} — run qc_adni_wmh_t1_mni.py first")
    qc = pd.read_csv(QC)
    primary = qc[qc["qc_primary"] == True].copy()  # noqa: E712
    strict = qc[qc["qc_sensitivity_strict"] == True].copy()  # noqa: E712
    if len(primary) < 5:
        raise SystemExit(f"Too few QC-pass subjects: {len(primary)}")

    mra = pd.read_csv(MRA)
    geo = {}
    for _, r in mra.iterrows():
        geo[f"mra_{r['metric']}_{r['territory']}"] = float(r["mean"])

    inv = pd.read_csv(INV)
    inv["ptid"] = inv["subject_id"].astype(str)
    dem = inv[["ptid", "PTGENDER", "PTDOBYY", "image_date", "TOTAL_WMH"]].drop_duplicates("ptid")

    def prep(df: pd.DataFrame) -> pd.DataFrame:
        out = df.merge(dem, on="ptid", how="left", suffixes=("", "_inv"))
        if "TOTAL_WMH_ucd" not in out.columns and "TOTAL_WMH" in out.columns:
            out["TOTAL_WMH_ucd"] = out["TOTAL_WMH"]
        scan_year = pd.to_datetime(out["image_date"], errors="coerce").dt.year
        dob = pd.to_datetime(out["PTDOBYY"], errors="coerce")
        birth_year = dob.dt.year.fillna(pd.to_numeric(out["PTDOBYY"], errors="coerce"))
        out["age_approx"] = scan_year - birth_year
        out["sex_f"] = (pd.to_numeric(out["PTGENDER"], errors="coerce") == 2).astype(float)
        return out

    primary = prep(primary)
    strict = prep(strict)
    n_primary, n_strict = len(primary), len(strict)

    if args.analysis_set in ("primary", "both"):
        _run_nest(primary, "primary", n_primary, n_strict, geo)
    if args.analysis_set in ("strict", "both"):
        if len(strict) < 5:
            raise SystemExit(f"Too few strict QC subjects: {len(strict)}")
        _run_nest(strict, "strict", n_primary, n_strict, geo)


if __name__ == "__main__":
    main()
