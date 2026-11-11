"""
Exploratory Aim 4 molecular x ATLAS nest (ecological, control-trained).

Design (anti-circularity):
  - Molecular features from Allen multi-region Endo ONLY (within-atlas, unsorted).
  - Arms: ACA=tv_arm ACA (ACC), MCA=tv_arm MCA (not MCA_leaning), PCA=tv_arm PCA (V1).
  - VB: no molecular arm on disk -> excluded from molecular models.
  - Modules: same non-neuronal panel as Class B hardening (harden_class_b.MODULES).
  - Aggregate: cell -> donor mean -> territory mean.
  - Outcome: ATLAS primary mean lesion_share / mean lesion_burden by territory.
  - Null: MRA geometry (vessel_probability, vessel_radius) by territory.

LOCKED: exploratory ecological nest (n_territories=3). NOT confirmatory H3.
NOT enrichment-matched. Soft DE genes are NOT used as features.
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
ATLAS_TERR = ROOT / "results/aim4/atlas/aim4_nest_atlas_mra_territory_summary.csv"
MRA = ROOT / "results/aim4/mra_geometry_by_territory.csv"
CLASS_B_MOD = ROOT / "results/class_b/class_b_module_scores_by_donor.csv"
OUT = ROOT / "results/aim4/atlas"

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
    if not ALLEN.exists():
        raise SystemExit(f"Missing {ALLEN}")
    if not ATLAS_TERR.exists():
        raise SystemExit(f"Missing {ATLAS_TERR}")

    log("Loading Allen multi-region Endo …")
    a = ad.read_h5ad(ALLEN)
    # restrict to clean TV arms
    arm_map = {"ACA": "ACA", "MCA": "MCA", "PCA": "PCA"}
    keep = a.obs["tv_arm"].astype(str).isin(arm_map)
    a = a[keep].copy()
    a.obs["territory"] = a.obs["tv_arm"].astype(str).map(arm_map)

    for name, genes in MODULES.items():
        a.obs[f"mod_{name}"] = mean_module(a, genes)

    mod_cols = [f"mod_{k}" for k in MODULES]
    # donor x territory means
    a.obs["donor_id"] = a.obs["donor_id"].astype(str)
    donor_terr = (
        a.obs.groupby(["territory", "donor_id"], observed=True)[mod_cols]
        .mean()
        .reset_index()
    )
    terr_mol = donor_terr.groupby("territory", observed=True)[mod_cols].mean()
    terr_n_donors = donor_terr.groupby("territory", observed=True)["donor_id"].nunique()
    terr_n_cells = a.obs.groupby("territory", observed=True).size()

    atlas = pd.read_csv(ATLAS_TERR).set_index("territory")
    mra = pd.read_csv(MRA)
    geo = {}
    for _, r in mra.iterrows():
        geo[(r["territory"], r["metric"])] = float(r["mean"])

    rows = []
    for t in ["ACA", "MCA", "PCA"]:
        if t not in terr_mol.index or t not in atlas.index:
            continue
        row = {
            "territory": t,
            "n_donors_allen": int(terr_n_donors.get(t, 0)),
            "n_cells_allen": int(terr_n_cells.get(t, 0)),
            "atlas_mean_lesion_share": float(atlas.loc[t, "mean_lesion_share"]),
            "atlas_mean_lesion_burden": float(atlas.loc[t, "mean_lesion_burden"]),
            "mra_vessel_probability": geo.get((t, "vessel_probability")),
            "mra_vessel_radius": geo.get((t, "vessel_radius")),
        }
        for c in mod_cols:
            row[c] = float(terr_mol.loc[t, c])
        rows.append(row)
    eco = pd.DataFrame(rows)
    eco_path = OUT / "aim4_molecular_atlas_ecological_table.csv"
    eco.to_csv(eco_path, index=False)
    donor_terr.to_csv(OUT / "aim4_allen_module_scores_by_donor_territory.csv", index=False)

    # Spearman: each module vs ATLAS lesion share
    spearman = []
    y = eco["atlas_mean_lesion_share"].to_numpy()
    for c in mod_cols + ["mra_vessel_probability", "mra_vessel_radius"]:
        x = eco[c].to_numpy()
        if len(eco) >= 3 and np.nanstd(x) > 0:
            rho, p = stats.spearmanr(x, y)
        else:
            rho, p = np.nan, np.nan
        spearman.append({"feature": c, "spearman_rho_vs_lesion_share": rho, "p": p, "n": len(eco)})
    spear_df = pd.DataFrame(spearman)
    spear_df.to_csv(OUT / "aim4_molecular_atlas_spearman.csv", index=False)

    # Nested R2 on n=3 (descriptive only)
    Xg = eco[["mra_vessel_probability", "mra_vessel_radius"]].to_numpy()
    Xm = eco[mod_cols].to_numpy()
    Xgm = np.hstack([Xg, Xm])
    nest = {
        "n_territories": int(len(eco)),
        "territories": eco["territory"].tolist(),
        "r2_geometry_only": r2_ols(Xg, y),
        "r2_modules_only": r2_ols(Xm, y),
        "r2_geometry_plus_modules": r2_ols(Xgm, y),
        "note": "n=3 ecological points; R2 is descriptive only, not confirmatory H3",
    }

    # Class B ACA vs MCA module means (2-point sensitivity; no PCA)
    classb_note = None
    if CLASS_B_MOD.exists():
        cb = pd.read_csv(CLASS_B_MOD)
        # Allen ACA only + Garcia MCA (prefer ROSMAP as closer to control cortex)
        aca = cb[cb["arm"] == "ACA"][mod_cols].mean()
        mca = cb[(cb["arm"] == "MCA") & (cb["batch"] == "Garcia_ROSMAP")][mod_cols].mean()
        if mca.isna().all():
            mca = cb[cb["arm"] == "MCA"][mod_cols].mean()
        classb = pd.DataFrame(
            [
                {"territory": "ACA", **aca.to_dict(), "atlas_mean_lesion_share": float(atlas.loc["ACA", "mean_lesion_share"])},
                {"territory": "MCA", **mca.to_dict(), "atlas_mean_lesion_share": float(atlas.loc["MCA", "mean_lesion_share"])},
            ]
        )
        classb.to_csv(OUT / "aim4_classb_modules_vs_atlas_share_2point.csv", index=False)
        # direction check: does module delta ACA-MCA align with share delta?
        deltas = []
        for c in mod_cols:
            deltas.append(
                {
                    "feature": c,
                    "delta_module_MCA_minus_ACA": float(classb.loc[1, c] - classb.loc[0, c]),
                    "delta_share_MCA_minus_ACA": float(
                        classb.loc[1, "atlas_mean_lesion_share"] - classb.loc[0, "atlas_mean_lesion_share"]
                    ),
                }
            )
        pd.DataFrame(deltas).to_csv(OUT / "aim4_classb_module_share_deltas.csv", index=False)
        classb_note = "2-point Class B module means (ACA pool vs Garcia MCA) saved; cannot estimate Spearman with n=2."

    summary = {
        "design": "ecological_allen_modules_x_ATLAS_lesion_share",
        "verdict": "EXPLORATORY_ECOLOGICAL_n3__NOT_CONFIRMATORY_H3",
        "allen_cells_by_territory": {t: int(terr_n_cells.get(t, 0)) for t in ["ACA", "MCA", "PCA"]},
        "allen_donors_by_territory": {t: int(terr_n_donors.get(t, 0)) for t in ["ACA", "MCA", "PCA"]},
        "vb_molecular_arm": None,
        "nested_r2": nest,
        "spearman": spearman,
        "classb_sensitivity": classb_note,
        "interpretation_locks": [
            "Ecological nest only (n_territories=3: ACA/MCA/PCA); VB has no molecular arm.",
            "Allen Endo is unsorted cortex — not enrichment-matched vessel-enriched ACC.",
            "NOT confirmatory H3 / Territorial Vulnerability Index.",
            "Do not use soft unmatched DE genes as molecular predictors.",
            "ATLAS lesion geography and MRA geometry null remain the frozen imaging layer.",
        ],
    }
    (OUT / "aim4_molecular_atlas_nest_summary.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8"
    )

    md = [
        "# Aim 4 exploratory molecular × ATLAS nest",
        "",
        f"**Verdict:** `{summary['verdict']}`",
        "",
        "## Design",
        "- Molecular: Allen multi-region Endo module means (ACA/MCA/PCA; VB absent).",
        "- Modules: BBB / angiogenesis / ECM / hypoxia / inflammation (Class B panel).",
        "- Outcome: ATLAS primary mean lesion share by territory.",
        "- Null: MRA vessel probability + radius.",
        "",
        "## Territory table",
        eco.to_string(index=False),
        "",
        "## Spearman (module/geometry vs ATLAS lesion share)",
        spear_df.to_string(index=False),
        "",
        "## Nested R² (descriptive only; n=3)",
        f"- Geometry only: {nest['r2_geometry_only']}",
        f"- Modules only: {nest['r2_modules_only']}",
        f"- Geometry + modules: {nest['r2_geometry_plus_modules']}",
        "",
        "## Locked interpretation",
        "- This is an **exploratory ecological** screen, **not** confirmatory molecular H3.",
        "- n_territories=3 precludes stable incremental inference; VB molecular arm missing.",
        "- Enrichment-matched vessel-enriched ACC still required for confirmatory Class B / H3.",
        "- Imaging layer (ATLAS geography + MRA null) remains frozen independently.",
        "",
        f"Tables: `{eco_path.relative_to(ROOT)}`, `aim4_molecular_atlas_spearman.csv`",
        "",
    ]
    (OUT / "aim4_molecular_atlas_nest_report.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
