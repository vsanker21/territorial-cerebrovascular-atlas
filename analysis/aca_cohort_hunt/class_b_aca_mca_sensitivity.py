"""Class B sensitivity analyses for ACA vs MCA Harmony object.

Hardens SUPPORTIVE claim with:
  1. Leave-one-ACA-source-out donor OOF AUROC
  2. PLP1 ambient filter (oligo contamination check)
  3. Donor pseudobulk Wilcoxon DE (FDR)
  4. Balanced cell subsample (equal n per arm) cell CV

Uses existing integrated object when present; otherwise rebuilds via class_b loader.
"""
from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[2]
H5AD = ROOT / "data/processed/class_b/aca_mca_harmony_endo.h5ad"
OUT_DIR = ROOT / "results/class_b"
OUT_JSON = OUT_DIR / "aca_mca_sensitivity_summary.json"
OUT_MD = OUT_DIR / "aca_mca_sensitivity_report.md"
OUT_PB = OUT_DIR / "aca_mca_donor_pseudobulk_de.csv"
OUT_FOLDS = OUT_DIR / "aca_mca_sensitivity_donor_oof.csv"

OLIGO = ["PLP1", "MBP", "MOBP", "MAG", "OPALIN"]


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def dense(X):
    return X.toarray() if sparse.issparse(X) else np.asarray(X)


def donor_oof_auroc(X, y, donor, min_cells: int = 3):
    donors = pd.Series(donor)
    usable = [d for d, n in donors.value_counts().items() if n >= min_cells]
    rows = []
    oof = np.full(len(y), np.nan)
    for d in usable:
        te = donor == d
        tr = ~te
        if len(np.unique(y[tr])) < 2:
            continue
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], y[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        oof[te] = proba
        rows.append(
            {
                "held_out": d,
                "mean_mca_score": float(proba.mean()),
                "true_label": int(y[te][0]),
                "n_test": int(te.sum()),
            }
        )
    ddf = pd.DataFrame(rows)
    auc = np.nan
    if len(ddf) >= 4 and ddf["true_label"].nunique() > 1:
        auc = float(roc_auc_score(ddf["true_label"], ddf["mean_mca_score"]))
    return auc, ddf


def cell_cv_auroc(X, y, n_splits: int = 5):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=0)
    aucs = []
    for tr, te in skf.split(X, y):
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], clf.predict_proba(X[te])[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs))


def reharmony(a: ad.AnnData) -> ad.AnnData:
    b = a.copy()
    sc.pp.highly_variable_genes(b, n_top_genes=2000, flavor="seurat", subset=False)
    hvg = b[:, b.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=30, svd_solver="arpack")
    sc.external.pp.harmony_integrate(hvg, key="batch", basis="X_pca", adjusted_basis="X_pca_harmony")
    b.obsm["X_pca_harmony"] = hvg.obsm["X_pca_harmony"]
    return b


def oligo_ok(a: ad.AnnData, cap: float = 1.5) -> np.ndarray:
    present = [g for g in OLIGO if g in a.var_names]
    if not present:
        return np.ones(a.n_obs, dtype=bool)
    X = dense(a[:, present].X)
    return (X <= cap).all(axis=1)


def donor_pseudobulk_de(a: ad.AnnData, min_cells: int = 3, max_genes: int = 3000) -> pd.DataFrame:
    """Mean expression per donor_key within arm; Wilcoxon across donors."""
    ycol = "class_b_arm"
    # Restrict to most variable genes for speed / interpretability
    Xfull = dense(a.X)
    v = Xfull.var(axis=0)
    top = np.argsort(v)[::-1][:max_genes]
    genes = a.var_names.to_numpy()[top]
    X = Xfull[:, top]
    means = []
    meta = []
    for (arm, donor), idx in a.obs.groupby([ycol, "donor_key"]).groups.items():
        ix = a.obs_names.get_indexer(idx)
        if len(ix) < min_cells:
            continue
        means.append(X[ix].mean(axis=0))
        meta.append({"arm": str(arm), "donor_key": donor, "n": len(ix)})
    if len(means) < 6:
        return pd.DataFrame()
    M = np.vstack(means)
    meta = pd.DataFrame(meta)
    aca = (meta["arm"] == "ACA").to_numpy()
    mca = (meta["arm"] == "MCA").to_numpy()
    if aca.sum() < 3 or mca.sum() < 3:
        return pd.DataFrame()
    rows = []
    for j, g in enumerate(genes):
        xa = M[aca, j]
        xm = M[mca, j]
        try:
            stat, p = stats.mannwhitneyu(xm, xa, alternative="two-sided")
        except ValueError:
            continue
        rows.append(
            {
                "gene": str(g),
                "mean_ACA": float(xa.mean()),
                "mean_MCA": float(xm.mean()),
                "delta": float(xm.mean() - xa.mean()),
                "U": float(stat),
                "pval": float(p),
                "n_aca_donors": int(aca.sum()),
                "n_mca_donors": int(mca.sum()),
            }
        )
    de = pd.DataFrame(rows)
    if de.empty:
        return de
    de["padj"] = multipletests(de["pval"], method="fdr_bh")[1]
    return de.sort_values("pval")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log(f"Reading {H5AD}")
    a = sc.read_h5ad(H5AD)
    assert "X_pca_harmony" in a.obsm
    assert "class_b_arm" in a.obs

    X = a.obsm["X_pca_harmony"]
    y = (a.obs["class_b_arm"].astype(str) == "MCA").astype(int).to_numpy()
    donor = a.obs["donor_key"].astype(str).to_numpy()
    batch = a.obs["batch"].astype(str).to_numpy()

    results = {}
    fold_frames = []

    # Baseline
    base_auc, base_df = donor_oof_auroc(X, y, donor)
    base_df["analysis"] = "baseline"
    fold_frames.append(base_df)
    cell_m, cell_s = cell_cv_auroc(X, y)
    results["baseline"] = {
        "donor_oof_auroc": base_auc,
        "cell_cv_auroc_mean": cell_m,
        "cell_cv_auroc_std": cell_s,
        "n_aca": int((y == 0).sum()),
        "n_mca": int((y == 1).sum()),
        "n_aca_donors": int(base_df.loc[base_df.true_label == 0].shape[0]),
        "n_mca_donors": int(base_df.loc[base_df.true_label == 1].shape[0]),
    }
    log(f"baseline donor OOF={base_auc:.3f} cell={cell_m:.3f}")

    # Leave-one-ACA-source-out
    aca_sources = sorted(a.obs.loc[a.obs["class_b_arm"].astype(str) == "ACA", "batch"].astype(str).unique())
    for src in aca_sources:
        keep = ~((a.obs["class_b_arm"].astype(str) == "ACA") & (a.obs["batch"].astype(str) == src))
        # also drop if this leaves no ACA
        sub = a[keep].copy()
        if (sub.obs["class_b_arm"].astype(str) == "ACA").sum() < 20:
            log(f"skip drop_{src}: too few ACA left")
            continue
        sub = reharmony(sub)
        Xs = sub.obsm["X_pca_harmony"]
        ys = (sub.obs["class_b_arm"].astype(str) == "MCA").astype(int).to_numpy()
        ds = sub.obs["donor_key"].astype(str).to_numpy()
        auc, ddf = donor_oof_auroc(Xs, ys, ds)
        ddf["analysis"] = f"drop_ACA_source_{src}"
        fold_frames.append(ddf)
        cm, cs = cell_cv_auroc(Xs, ys)
        results[f"drop_{src}"] = {
            "donor_oof_auroc": auc,
            "cell_cv_auroc_mean": cm,
            "cell_cv_auroc_std": cs,
            "n_aca": int((ys == 0).sum()),
            "n_mca": int((ys == 1).sum()),
            "n_aca_donors": int((ddf.true_label == 0).sum()),
            "n_mca_donors": int((ddf.true_label == 1).sum()),
        }
        log(f"drop {src}: donor OOF={auc:.3f} cell={cm:.3f} n_aca={(ys==0).sum()}")

    # Single ACA source only vs MCA
    for src in aca_sources:
        keep = (a.obs["class_b_arm"].astype(str) == "MCA") | (
            (a.obs["class_b_arm"].astype(str) == "ACA") & (a.obs["batch"].astype(str) == src)
        )
        sub = a[keep].copy()
        if (sub.obs["class_b_arm"].astype(str) == "ACA").sum() < 15:
            continue
        sub = reharmony(sub)
        Xs = sub.obsm["X_pca_harmony"]
        ys = (sub.obs["class_b_arm"].astype(str) == "MCA").astype(int).to_numpy()
        ds = sub.obs["donor_key"].astype(str).to_numpy()
        auc, ddf = donor_oof_auroc(Xs, ys, ds)
        ddf["analysis"] = f"only_ACA_source_{src}"
        fold_frames.append(ddf)
        cm, cs = cell_cv_auroc(Xs, ys)
        results[f"only_{src}"] = {
            "donor_oof_auroc": auc,
            "cell_cv_auroc_mean": cm,
            "cell_cv_auroc_std": cs,
            "n_aca": int((ys == 0).sum()),
            "n_mca": int((ys == 1).sum()),
            "n_aca_donors": int((ddf.true_label == 0).sum()),
            "n_mca_donors": int((ddf.true_label == 1).sum()),
        }
        log(f"only {src}: donor OOF={auc:.3f} cell={cm:.3f} n_aca={(ys==0).sum()}")

    # PLP1 / oligo ambient filter
    ok = oligo_ok(a, cap=1.5)
    log(f"oligo filter retain {ok.sum()}/{a.n_obs}")
    sub = a[ok].copy()
    if (sub.obs["class_b_arm"].astype(str) == "ACA").sum() >= 20 and (
        sub.obs["class_b_arm"].astype(str) == "MCA"
    ).sum() >= 20:
        sub = reharmony(sub)
        Xs = sub.obsm["X_pca_harmony"]
        ys = (sub.obs["class_b_arm"].astype(str) == "MCA").astype(int).to_numpy()
        ds = sub.obs["donor_key"].astype(str).to_numpy()
        auc, ddf = donor_oof_auroc(Xs, ys, ds)
        ddf["analysis"] = "oligo_ambient_filter"
        fold_frames.append(ddf)
        cm, cs = cell_cv_auroc(Xs, ys)
        results["oligo_ambient_filter"] = {
            "donor_oof_auroc": auc,
            "cell_cv_auroc_mean": cm,
            "cell_cv_auroc_std": cs,
            "n_aca": int((ys == 0).sum()),
            "n_mca": int((ys == 1).sum()),
            "n_retained": int(ok.sum()),
            "n_dropped": int((~ok).sum()),
        }
        log(f"oligo filter: donor OOF={auc:.3f} cell={cm:.3f}")

    # Balanced subsample
    rng = np.random.default_rng(0)
    aca_ix = np.where(y == 0)[0]
    mca_ix = np.where(y == 1)[0]
    n = min(len(aca_ix), len(mca_ix))
    ix = np.concatenate([rng.choice(aca_ix, n, replace=False), rng.choice(mca_ix, n, replace=False)])
    cm, cs = cell_cv_auroc(X[ix], y[ix])
    auc, ddf = donor_oof_auroc(X[ix], y[ix], donor[ix], min_cells=2)
    ddf["analysis"] = "balanced_subsample"
    fold_frames.append(ddf)
    results["balanced_subsample"] = {
        "donor_oof_auroc": auc,
        "cell_cv_auroc_mean": cm,
        "cell_cv_auroc_std": cs,
        "n_per_arm": int(n),
    }
    log(f"balanced n={n}: donor OOF={auc:.3f} cell={cm:.3f}")

    # Donor pseudobulk DE on expression (not Harmony space)
    log("Donor pseudobulk DE")
    # Use log-normalized X; may already be log1p
    pb = donor_pseudobulk_de(a, min_cells=3)
    if len(pb):
        pb.to_csv(OUT_PB, index=False)
        results["pseudobulk_de"] = {
            "n_genes_tested": int(len(pb)),
            "n_fdr_05": int((pb["padj"] < 0.05).sum()),
            "n_fdr_10": int((pb["padj"] < 0.1).sum()),
            "top_mca": pb.sort_values("delta", ascending=False).head(10)[["gene", "delta", "padj"]].to_dict("records"),
            "top_aca": pb.sort_values("delta", ascending=True).head(10)[["gene", "delta", "padj"]].to_dict("records"),
        }
        log(f"pseudobulk FDR<0.05: {(pb['padj']<0.05).sum()}")
    else:
        results["pseudobulk_de"] = {"n_genes_tested": 0}

    pd.concat(fold_frames, ignore_index=True).to_csv(OUT_FOLDS, index=False)

    # Verdict aggregation
    donor_aucs = [v.get("donor_oof_auroc") for k, v in results.items() if isinstance(v, dict) and "donor_oof_auroc" in v]
    donor_aucs = [x for x in donor_aucs if x == x]
    supportive = sum(1 for x in donor_aucs if x >= 0.7)
    results["meta"] = {
        "n_sensitivity_donor_aucs": len(donor_aucs),
        "n_supportive_ge_0.7": supportive,
        "min_donor_oof": float(np.min(donor_aucs)) if donor_aucs else None,
        "median_donor_oof": float(np.median(donor_aucs)) if donor_aucs else None,
        "verdict": (
            "ROBUST_SUPPORTIVE"
            if supportive >= max(2, len(donor_aucs) - 1) and np.min(donor_aucs) >= 0.65
            else "SUPPORTIVE_WITH_FRAGILITY"
            if results["baseline"]["donor_oof_auroc"] >= 0.7
            else "MARGINAL"
        ),
    }

    OUT_JSON.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")

    # Markdown
    lines = [
        "# Class B sensitivity — ACA vs MCA",
        "",
        f"**Meta verdict:** `{results['meta']['verdict']}`",
        "",
        f"Baseline donor OOF AUROC **{results['baseline']['donor_oof_auroc']:.3f}**; "
        f"cell CV **{results['baseline']['cell_cv_auroc_mean']:.3f}**.",
        "",
        "## Donor OOF across sensitivities",
        "",
        "| Analysis | Donor OOF | Cell CV | n ACA | n MCA |",
        "|---|---:|---:|---:|---:|",
    ]
    for k, v in results.items():
        if k in ("meta", "pseudobulk_de") or not isinstance(v, dict) or "donor_oof_auroc" not in v:
            continue
        lines.append(
            f"| {k} | {v['donor_oof_auroc']:.3f} | {v.get('cell_cv_auroc_mean', float('nan')):.3f} | "
            f"{v.get('n_aca', v.get('n_per_arm', ''))} | {v.get('n_mca', v.get('n_per_arm', ''))} |"
        )
    pb = results.get("pseudobulk_de", {})
    lines += [
        "",
        "## Donor pseudobulk DE",
        "",
        f"- Genes tested: {pb.get('n_genes_tested', 0)}",
        f"- FDR < 0.05: **{pb.get('n_fdr_05', 0)}**",
        f"- FDR < 0.10: {pb.get('n_fdr_10', 0)}",
        "",
        "## Interpretation",
        "",
        "- Leave-one-ACA-source-out tests whether signal is driven by a single ACA cohort.",
        "- Oligo ambient filter addresses PLP1 among cell-level DE hits.",
        "- Pseudobulk DE is the stricter cross-donor molecular claim.",
        "",
        f"Objects: `{OUT_JSON.name}`, `{OUT_PB.name}`, `{OUT_FOLDS.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_MD}")
    log(json.dumps(results["meta"], indent=2))


if __name__ == "__main__":
    main()
