"""
Aim 4 nest — ATLAS territorial lesion fractions + MRA geometry null.

Design:
  - Outcome: total lesion voxels (1 mm³); descriptive geography of TV shares
  - Geometry null: Mouches/Forkert MRA atlas constants (identical across subjects)
  - Features: lesion_burden_in_{ACA,MCA,PCA,VB}
  - Analysis sets: QC primary / strict (FLAG_LOW_TV excluded)

Interpretation lock: M3−M0 with territorial burdens vs total lesion volume is a
partition/concordance-style check of the same lesion map — NOT molecular H3.
Molecular H3 still requires control-trained transcriptomic features.
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
OUT = ROOT / "results/aim4/atlas"
FRAC = OUT / "atlas_lesion_tv_fractions.csv"
MRA = ROOT / "results/aim4/mra_geometry_by_territory.csv"


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def r2_ols(X, y):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(X).all(axis=1) & np.isfinite(y)
    if mask.sum() < 8:
        return np.nan, int(mask.sum())
    Xs = StandardScaler().fit_transform(X[mask])
    return float(LinearRegression().fit(Xs, y[mask]).score(Xs, y[mask])), int(mask.sum())


def _run(analysis: pd.DataFrame, label: str, geo: dict, n_primary: int, n_strict: int):
    for k, v in geo.items():
        analysis[k] = v

    burden = [f"lesion_burden_in_{t}" for t in ["ACA", "MCA", "PCA", "VB"]]
    share = [f"lesion_frac_of_lesion_{t}" for t in ["ACA", "MCA", "PCA", "VB"]]
    geo_feat = [
        f"mra_vessel_probability_{t}" for t in ["ACA", "MCA", "PCA", "VB"]
    ] + [f"mra_vessel_radius_{t}" for t in ["ACA", "MCA", "PCA", "VB"]]
    geo_feat = [c for c in geo_feat if c in analysis.columns]

    # limited demog: chronicity (0/1/NaN)
    analysis["chronicity_num"] = pd.to_numeric(analysis.get("chronicity"), errors="coerce")
    y = analysis["n_lesion_voxels"].astype(float).to_numpy()

    X0 = analysis[["chronicity_num"]].to_numpy()
    # replace NaN chronicity with median for OLS completeness, track n
    if np.isnan(X0).all():
        r2_0, n0 = np.nan, 0
    else:
        med = np.nanmedian(X0)
        X0 = np.where(np.isnan(X0), med, X0)
        r2_0, n0 = r2_ols(X0, y)

    r2_1, n1 = r2_ols(analysis[geo_feat].to_numpy(), y)
    r2_2, n2 = r2_ols(analysis[geo_feat + burden].to_numpy(), y)
    if n0 > 0:
        X3 = np.column_stack([X0, analysis[burden].to_numpy()])
        r2_3, n3 = r2_ols(X3, y)
    else:
        r2_3, n3 = r2_ols(analysis[burden].to_numpy(), y)

    mca = analysis["lesion_burden_in_MCA"].astype(float)
    aca = analysis["lesion_burden_in_ACA"].astype(float)
    wil = (
        stats.wilcoxon(mca, aca, zero_method="wilcox", alternative="two-sided")
        if len(analysis) >= 10
        else None
    )

    terr = []
    for t in ["ACA", "MCA", "PCA", "VB"]:
        terr.append(
            {
                "territory": t,
                "mean_lesion_burden": float(analysis[f"lesion_burden_in_{t}"].mean()),
                "median_lesion_burden": float(analysis[f"lesion_burden_in_{t}"].median()),
                "mean_lesion_share": float(analysis[f"lesion_frac_of_lesion_{t}"].mean()),
                "mra_vessel_probability_mean": geo.get(f"mra_vessel_probability_{t}"),
                "mra_vessel_radius_mean": geo.get(f"mra_vessel_radius_{t}"),
            }
        )

    suffix = "" if label == "primary" else f"_{label}"
    keep = [
        "subject_id",
        "site",
        "n_lesion_voxels",
        "lesion_frac_assigned_tv",
        "dominant_territory",
        "dominant_share",
        "days_post_stroke",
        "chronicity",
        "qc_flags",
    ] + burden + share + geo_feat
    keep = [c for c in keep if c in analysis.columns]
    subj_path = OUT / f"aim4_nest_atlas_mra_subject_table{suffix}.csv"
    analysis[keep].to_csv(subj_path, index=False)
    pd.DataFrame(terr).to_csv(OUT / f"aim4_nest_atlas_mra_territory_summary{suffix}.csv", index=False)

    dom = analysis["dominant_territory"].value_counts(dropna=False).to_dict()
    summary = {
        "analysis_set": label,
        "n_analysis": int(len(analysis)),
        "n_primary_qc": n_primary,
        "n_strict_qc": n_strict,
        "median_lesion_voxels": float(analysis["n_lesion_voxels"].median()),
        "dominant_territory_counts": {str(k): int(v) for k, v in dom.items()},
        "ols_r2_n_lesion_voxels": {
            "M0_chronicity": {"r2": r2_0, "n": n0},
            "M1_mra_geometry_null": {
                "r2": r2_1,
                "n": n1,
                "note": "atlas constants identical across subjects",
            },
            "M2_geometry_plus_tv_burden": {"r2": r2_2, "n": n2},
            "M3_chronicity_or_burden": {"r2": r2_3, "n": n3},
            "delta_M3_minus_M0": (r2_3 - r2_0)
            if (r2_3 == r2_3 and r2_0 == r2_0)
            else (r2_3 if r2_3 == r2_3 else None),
        },
        "wilcoxon_MCA_vs_ACA_burden": {
            "stat": float(wil.statistic) if wil else None,
            "p": float(wil.pvalue) if wil else None,
        },
        "territory_means": terr,
        "interpretation_locks": [
            "MRA geometry is atlas-level (not subject-specific); M1 R^2 near 0 is expected.",
            "M3 with territorial lesion burdens predicting total lesion volume is a partition-style check, not molecular H3.",
            "ATLAS is the independent stroke outcome arm (anti-circularity vs Liu construction).",
            "Molecular H3 still requires control-trained vascular transcriptomic features.",
            "Descriptive result: report dominant-territory counts and mean lesion shares.",
        ],
    }
    (OUT / f"aim4_nest_atlas_mra_summary{suffix}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    md = [
        f"# Aim 4 nest — ATLAS territorial lesion + MRA geometry ({label})",
        "",
        f"**Analysis set:** QC {label} n={len(analysis)} "
        f"(registry primary n={n_primary}, strict n={n_strict})",
        "",
        f"- Median lesion voxels: {summary['median_lesion_voxels']:.0f}",
        f"- Dominant territory counts: {summary['dominant_territory_counts']}",
        f"- Wilcoxon MCA vs ACA lesion burden: p="
        f"{(wil.pvalue if wil else float('nan')):.3g}",
        "",
        "## Nested OLS R² predicting total lesion voxels",
        f"| Model | R² | n |",
        f"|---|---:|---:|",
        f"| M0 chronicity | {r2_0 if r2_0==r2_0 else float('nan'):.3f} | {n0} |",
        f"| M1 MRA geometry null | {r2_1 if r2_1==r2_1 else float('nan'):.3f} | {n1} |",
        f"| M2 geometry + TV lesion burden | {r2_2 if r2_2==r2_2 else float('nan'):.3f} | {n2} |",
        f"| M3 chronicity/burden | {r2_3 if r2_3==r2_3 else float('nan'):.3f} | {n3} |",
        f"| **Δ (M3 − M0)** | **{(r2_3 - r2_0) if r2_3==r2_3 and r2_0==r2_0 else float('nan'):.3f}** | |",
        "",
        f"## Territory means ({label})",
        pd.DataFrame(terr).to_string(index=False),
        "",
        "## Locked interpretation",
        "- Atlas MRA geometry is a shared null; M1 R²≈0 is expected.",
        "- M3 predicting total lesion volume from territorial burdens of the *same* lesion map "
        "is a **partition-style check**, **not** molecular H3.",
        "- Molecular Aim 4 claim still requires control-trained transcriptomic features × "
        "independent territorial lesions.",
        "- Report descriptive stroke geography (shares / dominant territory) with ATLAS citation.",
        "",
        f"Tables: `{subj_path.relative_to(ROOT)}`",
        "",
    ]
    report = OUT / (
        "aim4_nest_atlas_mra_report.md" if label == "primary" else f"aim4_nest_atlas_mra_report_{label}.md"
    )
    report.write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(summary, indent=2, default=str))
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis-set", choices=["primary", "strict", "both"], default="both")
    args = ap.parse_args()

    if not FRAC.exists():
        raise SystemExit(f"Missing {FRAC} — run atlas_lesion_tv_pipeline.py first")
    df = pd.read_csv(FRAC)
    ok = df[df["status"] == "ok"].copy()
    primary = ok[ok["qc_primary"] == True].copy()  # noqa: E712
    strict = ok[ok["qc_sensitivity_strict"] == True].copy()  # noqa: E712
    if len(primary) < 5:
        raise SystemExit(f"Too few primary QC subjects: {len(primary)}")

    mra = pd.read_csv(MRA)
    geo = {}
    for _, r in mra.iterrows():
        geo[f"mra_{r['metric']}_{r['territory']}"] = float(r["mean"])

    n_primary, n_strict = len(primary), len(strict)
    if args.analysis_set in ("primary", "both"):
        _run(primary, "primary", geo, n_primary, n_strict)
    if args.analysis_set in ("strict", "both"):
        if len(strict) < 5:
            raise SystemExit(f"Too few strict QC: {len(strict)}")
        _run(strict, "strict", geo, n_primary, n_strict)


if __name__ == "__main__":
    main()
