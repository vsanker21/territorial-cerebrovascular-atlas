#!/usr/bin/env python3
"""Nature molecular rigor add-on: power tables + ambient-module bootstrap CIs.

Extends Midtemporal (MCA) vs Entorhinal (PCA) donor-aware contrasts:
  - Achieved power / MDE tables for paired and unpaired designs
  - Ambient score CIs and ambient-residualized module CIs (report all)
  - No soft DE marker lists; TRCI/BCI excluded
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(r"F:\Cerebrovascular SCS")
H5AD = ROOT / "data/processed/garcia/garcia_rosmap_vascular.h5ad"
OUT = ROOT / "results/tracks/nature_rigor"
OUT.mkdir(parents=True, exist_ok=True)

MCA_REG = "Midtemporal_cortex"
PCA_REG = "Entorhinal_cortex"
PAIR = (MCA_REG, PCA_REG)
MDE_D = 0.5

MODULES = {
    "BBB_tight_junction": ["CLDN5", "OCLN", "TJP1", "TJP2", "CDH5", "ESAM", "PECAM1"],
    "angiogenesis": ["VEGFA", "FLT1", "KDR", "ENG", "VWF", "ANGPT2", "TEK"],
    "ECM_basement": ["COL4A1", "COL4A2", "LAMA4", "LAMB1", "NID1", "HSPG2"],
    "hypoxia_stress": ["HIF1A", "VEGFA", "SLC2A1", "BNIP3", "DDIT4"],
    "inflammation": ["ICAM1", "VCAM1", "SELE", "CCL2", "IL6", "NFKB1"],
    "flow_shear": ["KLF2", "KLF4", "NOS3", "THBD", "PLAT", "MEF2C"],
    "mural_contractile": ["ACTA2", "TAGLN", "MYH11", "CNN1", "MYLK", "PDGFRB", "RGS5"],
}
AMBIENT = ["RBFOX1", "SNAP25", "SYT1", "SLC17A7", "GAD1", "MBP", "PLP1", "MOBP", "AQP4", "GFAP"]


def log(msg: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def log1p_norm(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    if X.max() > 50:
        lib = X.sum(axis=1, keepdims=True)
        lib[lib == 0] = 1.0
        X = np.log1p(X / lib * 1e4)
    return X


def module_mat(X: np.ndarray, var_names: list[str], genes: list[str]) -> np.ndarray:
    idx = [var_names.index(g) for g in genes if g in var_names]
    if not idx:
        return np.full(X.shape[0], np.nan)
    return X[:, idx].mean(axis=1)


def bootstrap_ci(x: np.ndarray, B: int = 2000, seed: int = 0) -> list[float]:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return [np.nan, np.nan]
    rng = np.random.default_rng(seed)
    boots = [float(rng.choice(x, size=len(x), replace=True).mean()) for _ in range(B)]
    return [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]


def paired_power_approx(n: int, d: float, alpha: float = 0.05) -> float:
    """Two-sided paired t approximate power for effect size d (mean/sd of paired diffs)."""
    if n < 3 or not np.isfinite(d):
        return np.nan
    from scipy.stats import norm, t as tdist

    df = n - 1
    ncp = d * np.sqrt(n)
    tcrit = tdist.ppf(1 - alpha / 2, df)
    # power ≈ 1 - F(tcrit; ncp) + F(-tcrit; ncp) using normal approx on noncentral t
    # Use scipy nct if available
    try:
        from scipy.stats import nct

        return float(1 - nct.cdf(tcrit, df, ncp) + nct.cdf(-tcrit, df, ncp))
    except Exception:
        # Normal approximation
        return float(1 - norm.cdf(tcrit - ncp) + norm.cdf(-tcrit - ncp))


def mde_at_power(n: int, power: float = 0.8, alpha: float = 0.05) -> float:
    """Approximate minimum |d| detectable at given power for paired design."""
    if n < 3:
        return np.nan
    # Binary search on d
    lo, hi = 0.01, 3.0
    for _ in range(40):
        mid = (lo + hi) / 2
        p = paired_power_approx(n, mid, alpha)
        if not np.isfinite(p):
            return np.nan
        if p < power:
            lo = mid
        else:
            hi = mid
    return float(hi)


def donor_means(meta: pd.DataFrame, score_col: str, celltype: str) -> pd.DataFrame:
    sub = meta[meta["celltype"] == celltype].copy()
    return (
        sub.groupby(["donor", "region"], observed=True)[score_col]
        .mean()
        .reset_index()
        .pivot(index="donor", columns="region", values=score_col)
    )


def main() -> None:
    log("Loading Garcia for power + ambient CI rigor")
    a = ad.read_h5ad(H5AD)
    obs = a.obs.copy()
    obs["donor"] = obs["subject"].astype(str)
    obs["region"] = obs["brain_region"].astype(str)
    obs["celltype"] = obs["celltype"].astype(str)
    mask = obs["region"].isin(PAIR).values
    idx = np.where(mask)[0]
    X = a.X[idx]
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = log1p_norm(np.asarray(X, dtype=float))
    var_names = list(a.var_names)
    meta = obs.iloc[idx].copy()

    amb = module_mat(X, var_names, AMBIENT)
    meta["ambient"] = amb

    for name, genes in MODULES.items():
        s = module_mat(X, var_names, genes)
        resid = s.copy()
        m = np.isfinite(s) & np.isfinite(amb)
        if m.sum() > 10 and np.nanstd(amb[m]) > 0:
            b1 = np.polyfit(amb[m], s[m], 1)
            resid = s - (b1[0] * amb + b1[1])
        meta[name] = s
        meta[f"{name}_resid"] = resid

    power_rows = []
    ambient_rows = []
    module_ci_rows = []

    for ctype in ["Endo", "Per", "SMC"]:
        # Ambient donor means by region
        piv_amb = donor_means(meta, "ambient", ctype)
        if MCA_REG in piv_amb.columns and PCA_REG in piv_amb.columns:
            both = piv_amb[[MCA_REG, PCA_REG]].dropna()
            a_vals = piv_amb[MCA_REG].dropna().values
            b_vals = piv_amb[PCA_REG].dropna().values
            if len(both) >= 3:
                d = both[PCA_REG].values - both[MCA_REG].values
                boots = []
                rng = np.random.default_rng(11)
                for _ in range(2000):
                    boots.append(float(rng.choice(d, size=len(d), replace=True).mean()))
                ambient_rows.append(
                    {
                        "celltype": ctype,
                        "score": "ambient",
                        "design": "paired",
                        "n_paired_donors": int(len(both)),
                        "mean_MCA": float(both[MCA_REG].mean()),
                        "mean_PCA": float(both[PCA_REG].mean()),
                        "delta_PCA_minus_MCA": float(np.mean(d)),
                        "ci95_lo": float(np.percentile(boots, 2.5)),
                        "ci95_hi": float(np.percentile(boots, 97.5)),
                        "wilcoxon_p": float(stats.wilcoxon(both[PCA_REG], both[MCA_REG]).pvalue),
                    }
                )
            ambient_rows.append(
                {
                    "celltype": ctype,
                    "score": "ambient",
                    "design": "region_marginal",
                    "n_MCA": int(len(a_vals)),
                    "n_PCA": int(len(b_vals)),
                    "mean_MCA": float(np.nanmean(a_vals)) if len(a_vals) else np.nan,
                    "mean_PCA": float(np.nanmean(b_vals)) if len(b_vals) else np.nan,
                    "ci95_MCA": bootstrap_ci(a_vals, seed=21),
                    "ci95_PCA": bootstrap_ci(b_vals, seed=22),
                }
            )

        for mod in MODULES:
            col = f"{mod}_resid"
            if ctype == "Endo" and mod == "mural_contractile":
                continue
            piv = donor_means(meta, col, ctype)
            if MCA_REG not in piv.columns or PCA_REG not in piv.columns:
                continue
            both = piv[[MCA_REG, PCA_REG]].dropna()
            n_paired = int(len(both))
            n_mca = int(piv[MCA_REG].notna().sum())
            n_pca = int(piv[PCA_REG].notna().sum())
            if n_paired >= 3:
                dvec = both[PCA_REG].values - both[MCA_REG].values
                sd = float(np.std(dvec, ddof=1)) if len(dvec) > 1 else np.nan
                d_hat = float(np.mean(dvec) / sd) if sd and sd > 0 else np.nan
                boots = []
                rng = np.random.default_rng(31)
                for _ in range(2000):
                    boots.append(float(rng.choice(dvec, size=len(dvec), replace=True).mean()))
                ci = [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]
                power_at_mde = paired_power_approx(n_paired, MDE_D)
                power_at_obs = paired_power_approx(n_paired, abs(d_hat)) if np.isfinite(d_hat) else np.nan
                mde80 = mde_at_power(n_paired, 0.8)
                module_ci_rows.append(
                    {
                        "celltype": ctype,
                        "module": mod,
                        "design": "paired_ambient_resid",
                        "n_paired_donors": n_paired,
                        "delta_PCA_minus_MCA": float(np.mean(dvec)),
                        "ci95_lo": ci[0],
                        "ci95_hi": ci[1],
                        "cohens_d_paired": d_hat,
                        "wilcoxon_p": float(stats.wilcoxon(both[PCA_REG], both[MCA_REG]).pvalue),
                        "ci_excludes_0": bool(ci[0] > 0 or ci[1] < 0),
                        "meets_MDE_d05": bool(abs(d_hat) >= MDE_D) if np.isfinite(d_hat) else False,
                    }
                )
                power_rows.append(
                    {
                        "celltype": ctype,
                        "module": mod,
                        "n_paired_donors": n_paired,
                        "n_MCA_donors": n_mca,
                        "n_PCA_donors": n_pca,
                        "observed_d": d_hat,
                        "prespec_MDE_d": MDE_D,
                        "power_at_prespec_MDE": power_at_mde,
                        "power_at_observed_d": power_at_obs,
                        "MDE_for_power_0.8": mde80,
                        "underpowered_vs_prespec": bool(
                            np.isfinite(power_at_mde) and power_at_mde < 0.8
                        ),
                        "reporting": (
                            "EXPLORATORY_UNDERPOWERED"
                            if n_paired < 10 or (np.isfinite(power_at_mde) and power_at_mde < 0.8)
                            else "POWER_OK_FOR_MDE_EXPLORATORY_UNTIL_REPLICATION"
                        ),
                    }
                )
            else:
                power_rows.append(
                    {
                        "celltype": ctype,
                        "module": mod,
                        "n_paired_donors": n_paired,
                        "n_MCA_donors": n_mca,
                        "n_PCA_donors": n_pca,
                        "observed_d": np.nan,
                        "prespec_MDE_d": MDE_D,
                        "power_at_prespec_MDE": paired_power_approx(max(n_paired, 3), MDE_D)
                        if n_paired >= 3
                        else np.nan,
                        "power_at_observed_d": np.nan,
                        "MDE_for_power_0.8": mde_at_power(n_paired, 0.8) if n_paired >= 3 else np.nan,
                        "underpowered_vs_prespec": True,
                        "reporting": "EXPLORATORY_UNDERPOWERED",
                    }
                )

    power_df = pd.DataFrame(power_rows)
    amb_df = pd.DataFrame(ambient_rows)
    mod_df = pd.DataFrame(module_ci_rows)
    power_df.to_csv(OUT / "nature_power_table.csv", index=False)
    amb_df.to_csv(OUT / "nature_ambient_cis.csv", index=False)
    mod_df.to_csv(OUT / "nature_module_ambient_resid_cis.csv", index=False)

    summary = {
        "date": time.strftime("%Y-%m-%d"),
        "pair": {"MCA": MCA_REG, "PCA": PCA_REG},
        "prespec_MDE_d": MDE_D,
        "n_power_rows": int(len(power_df)),
        "n_underpowered": int(power_df["underpowered_vs_prespec"].sum())
        if len(power_df)
        else 0,
        "modules_ci_excludes_0": int(mod_df["ci_excludes_0"].sum()) if len(mod_df) else 0,
        "verdict": "POWER_AND_AMBIENT_CI_TABLES_WRITTEN",
        "locks": [
            "No soft DE markers",
            "Donor-level primary",
            "Ambient residualization for modules",
            "TRCI/BCI excluded",
        ],
    }
    (OUT / "nature_power_ambient_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    lines = [
        "# Nature rigor — power tables + ambient module CIs",
        "",
        f"**Date:** {summary['date']}",
        f"**Verdict:** `{summary['verdict']}`",
        f"- Prespec MDE |d|≥{MDE_D}",
        f"- Power rows: {summary['n_power_rows']}; underpowered vs prespec: {summary['n_underpowered']}",
        f"- Ambient-resid modules with CI excluding 0: {summary['modules_ci_excludes_0']}",
        "",
        "## Files",
        "- `nature_power_table.csv`",
        "- `nature_ambient_cis.csv`",
        "- `nature_module_ambient_resid_cis.csv`",
        "",
        "## Interpretation lock",
        "Underpowered contrasts remain EXPLORATORY. Do not elevate to confirmatory DE/markers.",
        "Enrichment-matched ACC still required for Nature H1c.",
        "",
    ]
    # Highlight top module CI results
    if len(mod_df):
        lines.append("## Module Δ (PCA−MCA) with CIs (paired ambient-resid)")
        show = mod_df.sort_values("wilcoxon_p")
        for _, r in show.iterrows():
            lines.append(
                f"- **{r['celltype']} / {r['module']}**: Δ={r['delta_PCA_minus_MCA']:.4f} "
                f"CI[{r['ci95_lo']:.4f},{r['ci95_hi']:.4f}] d={r['cohens_d_paired']:.3f} "
                f"p={r['wilcoxon_p']:.4g} n={int(r['n_paired_donors'])}"
            )
        lines.append("")
    (OUT / "nature_power_ambient_report.md").write_text("\n".join(lines), encoding="utf-8")
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
