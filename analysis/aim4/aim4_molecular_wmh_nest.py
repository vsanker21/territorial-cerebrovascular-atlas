"""
Exploratory Aim 4 molecular x ADNI WMH nest (ecological, control-trained).

Parallel to aim4_molecular_atlas_nest.py with ADNI territorial WMH shares.
LOCKED: EXPLORATORY_ECOLOGICAL_n3__NOT_CONFIRMATORY_H3.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
ALLEN = ROOT / "data/processed/allen/allen_endo_multiregion.h5ad"
ADNI_TERR = ROOT / "results/aim4/adni/aim4_nest_wmh_mra_territory_summary.csv"
MRA = ROOT / "results/aim4/mra_geometry_by_territory.csv"
OUT = ROOT / "results/aim4/adni"

MODULES = {
    "BBB_tight_junction": ["CLDN5", "OCLN", "TJP1", "TJP2", "F11R", "ESAM"],
    "angiogenesis": ["VEGFA", "FLT1", "KDR", "PECAM1", "ENG", "VWF", "CDH5"],
    "ECM_basement": ["COL4A1", "COL4A2", "LAMA4", "LAMB1", "HSPG2", "NID1"],
    "hypoxia_stress": ["HIF1A", "VEGFA", "BNIP3", "SLC2A1", "DDIT3", "ATF4"],
    "inflammation": ["ICAM1", "VCAM1", "SELE", "CCL2", "IL6", "NFKBIA"],
}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def dense(X):
    return X.toarray() if sparse.issparse(X) else np.asarray(X)


def mean_module(a: ad.AnnData, genes: list[str]) -> np.ndarray:
    present = [g for g in genes if g in a.var_names]
    if not present:
        return np.zeros(a.n_obs)
    return dense(a[:, present].X).mean(axis=1)


def r2_ols(X, y):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(y) < 3:
        return np.nan
    Xs = StandardScaler().fit_transform(X)
    return float(LinearRegression().fit(Xs, y).score(Xs, y))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(ALLEN)
    keep = a.obs["tv_arm"].astype(str).isin(["ACA", "MCA", "PCA"])
    a = a[keep].copy()
    a.obs["territory"] = a.obs["tv_arm"].astype(str)
    for name, genes in MODULES.items():
        a.obs[f"mod_{name}"] = mean_module(a, genes)
    mod_cols = [f"mod_{k}" for k in MODULES]
    donor_terr = (
        a.obs.groupby(["territory", "donor_id"], observed=True)[mod_cols].mean().reset_index()
    )
    terr_mol = donor_terr.groupby("territory", observed=True)[mod_cols].mean()
    terr_n_donors = donor_terr.groupby("territory", observed=True)["donor_id"].nunique()
    terr_n_cells = a.obs.groupby("territory", observed=True).size()

    adni = pd.read_csv(ADNI_TERR).set_index("territory")
    # tolerate column names
    share_col = next(c for c in adni.columns if "share" in c.lower() or "frac" in c.lower() and "mean" in c.lower() or c.startswith("mean_"))
    # prefer mean_wmh_share / mean_frac
    for cand in ("mean_wmh_share", "mean_lesion_share", "mean_frac_of_wmh", "mean_wmh_frac_of_wmh"):
        if cand in adni.columns:
            share_col = cand
            break
    else:
        share_col = [c for c in adni.columns if "mean" in c.lower() and ("share" in c.lower() or "frac" in c.lower())][0]

    mra = pd.read_csv(MRA)
    geo = {(r["territory"], r["metric"]): float(r["mean"]) for _, r in mra.iterrows()}

    rows = []
    for t in ["ACA", "MCA", "PCA"]:
        if t not in terr_mol.index or t not in adni.index:
            continue
        row = {
            "territory": t,
            "n_donors_allen": int(terr_n_donors.get(t, 0)),
            "n_cells_allen": int(terr_n_cells.get(t, 0)),
            "adni_mean_wmh_share": float(adni.loc[t, share_col]),
            "mra_vessel_probability": geo.get((t, "vessel_probability")),
            "mra_vessel_radius": geo.get((t, "vessel_radius")),
        }
        for c in mod_cols:
            row[c] = float(terr_mol.loc[t, c])
        rows.append(row)
    eco = pd.DataFrame(rows)
    eco.to_csv(OUT / "aim4_molecular_wmh_ecological_table.csv", index=False)

    y = eco["adni_mean_wmh_share"].to_numpy()
    spearman = []
    for c in mod_cols + ["mra_vessel_probability", "mra_vessel_radius"]:
        x = eco[c].to_numpy()
        rho, p = stats.spearmanr(x, y) if len(eco) >= 3 and np.nanstd(x) > 0 else (np.nan, np.nan)
        spearman.append({"feature": c, "spearman_rho_vs_wmh_share": float(rho) if np.isfinite(rho) else None, "p": float(p) if np.isfinite(p) else None, "n": len(eco)})
    spear_df = pd.DataFrame(spearman)
    spear_df.to_csv(OUT / "aim4_molecular_wmh_spearman.csv", index=False)

    Xg = eco[["mra_vessel_probability", "mra_vessel_radius"]].to_numpy()
    Xm = eco[mod_cols].to_numpy()
    nest = {
        "n_territories": int(len(eco)),
        "territories": eco["territory"].tolist(),
        "r2_geometry_only": r2_ols(Xg, y),
        "r2_modules_only": r2_ols(Xm, y),
        "r2_geometry_plus_modules": r2_ols(np.hstack([Xg, Xm]), y),
        "share_column_used": share_col,
        "note": "n=3 ecological points; descriptive only, not confirmatory H3",
    }
    summary = {
        "design": "ecological_allen_modules_x_ADNI_wmh_share",
        "verdict": "EXPLORATORY_ECOLOGICAL_n3__NOT_CONFIRMATORY_H3",
        "nested_r2": nest,
        "spearman": spearman,
        "interpretation_locks": [
            "Ecological nest only (n_territories=3); VB molecular arm absent.",
            "NOT confirmatory H3.",
            "ADNI geography nest remains the frozen imaging claim.",
        ],
    }
    (OUT / "aim4_molecular_wmh_nest_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    md = [
        "# Aim 4 exploratory molecular × ADNI WMH nest",
        "",
        f"**Verdict:** `{summary['verdict']}`",
        "",
        "## Territory table",
        eco.to_string(index=False),
        "",
        "## Spearman",
        spear_df.to_string(index=False),
        "",
        "## Nested R² (descriptive; n=3)",
        f"- Geometry: {nest['r2_geometry_only']}",
        f"- Modules: {nest['r2_modules_only']}",
        f"- Geometry+modules: {nest['r2_geometry_plus_modules']}",
        "",
        "Locked: exploratory only; not confirmatory molecular H3.",
        "",
    ]
    (OUT / "aim4_molecular_wmh_nest_report.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
