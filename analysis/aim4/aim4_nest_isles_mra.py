"""Aim 4 nest — ISLES'22 territorial lesion fractions + MRA geometry null."""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/aim4/isles"
FRAC = OUT / "isles_lesion_tv_fractions.csv"
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
    geo_feat = [f"mra_vessel_probability_{t}" for t in ["ACA", "MCA", "PCA", "VB"]] + [
        f"mra_vessel_radius_{t}" for t in ["ACA", "MCA", "PCA", "VB"]
    ]

    y = analysis["n_lesion_voxels"].astype(float).to_numpy()
    # demog: age + sex when available
    dem_cols = []
    if "age_years" in analysis.columns and analysis["age_years"].notna().sum() >= 8:
        analysis["_age"] = analysis["age_years"].astype(float)
        dem_cols.append("_age")
    if "sex_f" in analysis.columns and analysis["sex_f"].notna().sum() >= 8:
        analysis["_sex"] = analysis["sex_f"].astype(float)
        dem_cols.append("_sex")
    if dem_cols:
        X0 = analysis[dem_cols].to_numpy()
        # impute nan with median
        for j in range(X0.shape[1]):
            col = X0[:, j]
            med = np.nanmedian(col)
            X0[:, j] = np.where(np.isnan(col), med, col)
        r2_0, n0 = r2_ols(X0, y)
    else:
        r2_0, n0 = np.nan, 0

    r2_1, n1 = r2_ols(analysis[geo_feat].to_numpy(), y)
    r2_2, n2 = r2_ols(analysis[geo_feat + burden].to_numpy(), y)
    if n0 > 0:
        r2_3, n3 = r2_ols(np.column_stack([X0, analysis[burden].to_numpy()]), y)
    else:
        r2_3, n3 = r2_ols(analysis[burden].to_numpy(), y)

    mca = analysis["lesion_burden_in_MCA"].astype(float)
    aca = analysis["lesion_burden_in_ACA"].astype(float)
    w = stats.wilcoxon(mca, aca, alternative="two-sided") if len(analysis) >= 8 else None

    terr_rows = []
    for t in ["ACA", "MCA", "PCA", "VB"]:
        terr_rows.append(
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
    pd.DataFrame(terr_rows).to_csv(OUT / f"aim4_nest_isles_mra_territory_summary{suffix}.csv", index=False)
    analysis.to_csv(OUT / f"aim4_nest_isles_mra_subject_table{suffix}.csv", index=False)

    summary = {
        "analysis_set": label,
        "n_analysis": int(len(analysis)),
        "n_primary_qc": n_primary,
        "n_strict_qc": n_strict,
        "median_lesion_voxels": float(analysis["n_lesion_voxels"].median()),
        "dominant_territory_counts": {
            str(k): int(v) for k, v in analysis["dominant_territory"].value_counts(dropna=False).items()
        },
        "ols_r2_n_lesion_voxels": {
            "M0_demographics": {"r2": r2_0, "n": n0},
            "M1_mra_geometry_null": {
                "r2": r2_1,
                "n": n1,
                "note": "atlas constants identical across subjects",
            },
            "M2_geometry_plus_tv_burden": {"r2": r2_2, "n": n2},
            "M3_demog_or_burden": {"r2": r2_3, "n": n3},
            "delta_M3_minus_M0": (r2_3 - r2_0) if np.isfinite(r2_3) and np.isfinite(r2_0) else r2_3,
        },
        "wilcoxon_MCA_vs_ACA_burden": {"stat": float(w.statistic), "p": float(w.pvalue)} if w else None,
        "mean_lesion_share": {t: float(analysis[f"lesion_frac_of_lesion_{t}"].mean()) for t in ["ACA", "MCA", "PCA", "VB"]},
        "interpretation_locks": [
            "ISLES'22 independent stroke replication vs ATLAS.",
            "M3 predicting total lesion from territorial burdens is partition — NOT molecular H3.",
            "MRA geometry null uses atlas constants (R2~0 expected).",
            "ArterialAtlas_level2 bilateral map only.",
        ],
    }
    (OUT / f"aim4_nest_isles_mra_summary{suffix}.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8"
    )
    md = [
        f"# Aim 4 nest — ISLES'22 x ArterialAtlas_level2 ({label})",
        "",
        f"n={summary['n_analysis']}; dominant={summary['dominant_territory_counts']}",
        f"Mean lesion share ACA/MCA/PCA/VB={summary['mean_lesion_share']}",
        f"M0 R2={r2_0}; M1 MRA null R2={r2_1}; M3 R2={r2_3}; delta={summary['ols_r2_n_lesion_voxels']['delta_M3_minus_M0']}",
        f"Wilcoxon MCA vs ACA burden: {summary['wilcoxon_MCA_vs_ACA_burden']}",
        "",
        "Locked: geography + geometry null; M3 is partition not molecular H3.",
        "",
    ]
    (OUT / f"aim4_nest_isles_mra_report{suffix}.md").write_text("\n".join(md), encoding="utf-8")
    return summary


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(FRAC)
    ok = df[df["status"] == "ok"].copy()
    primary = ok[ok["qc_primary"] == True].copy()  # noqa: E712
    strict = ok[ok["qc_sensitivity_strict"] == True].copy()  # noqa: E712
    mra = pd.read_csv(MRA)
    geo = {}
    for _, r in mra.iterrows():
        t = r["territory"]
        if r["metric"] == "vessel_probability":
            geo[f"mra_vessel_probability_{t}"] = float(r["mean"])
        elif r["metric"] == "vessel_radius":
            geo[f"mra_vessel_radius_{t}"] = float(r["mean"])

    log("Primary nest…")
    s1 = _run(primary, "primary", geo, len(primary), len(strict))
    log("Strict nest…")
    s2 = _run(strict, "strict", geo, len(primary), len(strict))
    log(json.dumps({"primary": s1, "strict_delta_note": s2["ols_r2_n_lesion_voxels"]}, indent=2, default=str))


if __name__ == "__main__":
    main()
