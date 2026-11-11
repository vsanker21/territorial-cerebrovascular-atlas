"""Garcia ROSMAP zonation-stratified same-territory null (PFC vs MCA isocortex).

Replication of weak H1 null with explicit capEndo / aEndo / vEndo strata.
"""
from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[2]
ROSMAP = ROOT / "data/processed/garcia/garcia_rosmap_vascular.h5ad"
OUT_DIR = ROOT / "results/h1"
OUT_JSON = OUT_DIR / "rosmap_zonation_null_summary.json"
OUT_MD = OUT_DIR / "rosmap_zonation_null_report.md"


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def dense(X):
    return X.toarray() if sparse.issparse(X) else np.asarray(X)


def cell_cv(X, y):
    if len(np.unique(y)) < 2 or len(y) < 40:
        return np.nan, np.nan
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    aucs = []
    for tr, te in skf.split(X, y):
        if len(np.unique(y[tr])) < 2:
            continue
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], clf.predict_proba(X[te])[:, 1]))
    return (float(np.mean(aucs)), float(np.std(aucs))) if aucs else (np.nan, np.nan)


def donor_oof(X, y, donor, min_cells=5):
    rows = []
    for d, n in pd.Series(donor).value_counts().items():
        if n < min_cells:
            continue
        te = donor == d
        tr = ~te
        if len(np.unique(y[tr])) < 2:
            continue
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], y[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        rows.append({"held_out": d, "score": float(proba.mean()), "label": int(y[te][0]), "n": int(te.sum())})
    df = pd.DataFrame(rows)
    auc = np.nan
    if len(df) >= 4 and df["label"].nunique() > 1:
        auc = float(roc_auc_score(df["label"], df["score"]))
    return auc, df


def arm_map(region: str) -> str:
    r = region.lower()
    if "prefrontal" in r or "pfc" in r or "dlpfc" in r:
        return "PFC_MCA_leaning"
    if "midtemporal" in r or "middle temporal" in r or "angular" in r:
        return "MCA_isocortex"
    if "hippocamp" in r or "entorhinal" in r:
        return "PCA_medial_temporal"
    return "other"


def run_stratum(endo: sc.AnnData, stratum: str):
    # PCA for features
    sc.pp.highly_variable_genes(endo, n_top_genes=1200, flavor="seurat", subset=False)
    h = endo[:, endo.var["highly_variable"]].copy()
    sc.pp.scale(h, max_value=10)
    sc.tl.pca(h, n_comps=20)
    X = h.obsm["X_pca"]
    # null: MCA_isocortex vs PFC
    sub = endo[endo.obs["arm"].isin(["MCA_isocortex", "PFC_MCA_leaning"])].copy()
    if sub.n_obs < 40:
        return {"stratum": stratum, "status": "too_few"}
    Xs = X[endo.obs_names.get_indexer(sub.obs_names)]
    y = (sub.obs["arm"].astype(str) == "MCA_isocortex").astype(int).to_numpy()
    donor = sub.obs["subject"].astype(str).to_numpy() if "subject" in sub.obs else sub.obs_names.astype(str)
    c_auc, c_sd = cell_cv(Xs, y)
    d_auc, _ = donor_oof(Xs, y, donor)
    return {
        "stratum": stratum,
        "contrast": "MCA_isocortex_vs_PFC_same_territory_null",
        "n": int(sub.n_obs),
        "n_mca_iso": int((y == 1).sum()),
        "n_pfc": int((y == 0).sum()),
        "cell_cv_auroc": c_auc,
        "cell_cv_std": c_sd,
        "donor_oof_auroc": d_auc,
        "interpretation": "Same-territory region null — high AUROC = region/batch, not territory",
    }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    a = sc.read_h5ad(ROSMAP)
    endo = a[a.obs["subtype"].astype(str).isin(["capEndo", "aEndo", "vEndo"])].copy()
    region_col = "brain_region" if "brain_region" in endo.obs else "region"
    endo.obs["arm"] = endo.obs[region_col].astype(str).map(arm_map)
    results = {
        "n_endo": int(endo.n_obs),
        "by_subtype": endo.obs["subtype"].astype(str).value_counts().to_dict(),
        "by_arm": endo.obs["arm"].astype(str).value_counts().to_dict(),
        "strata": [],
    }
    for st in ["all", "capEndo", "aEndo", "vEndo"]:
        if st == "all":
            sub = endo
        else:
            sub = endo[endo.obs["subtype"].astype(str) == st].copy()
        if sub.n_obs < 40:
            results["strata"].append({"stratum": st, "status": "too_few", "n": int(sub.n_obs)})
            continue
        # Need to compute PCA on this subset directly
        sc.pp.highly_variable_genes(sub, n_top_genes=1200, flavor="seurat", subset=False)
        h = sub[:, sub.var["highly_variable"]].copy()
        sc.pp.scale(h, max_value=10)
        sc.tl.pca(h, n_comps=20)
        keep = sub.obs["arm"].isin(["MCA_isocortex", "PFC_MCA_leaning"])
        s2 = sub[keep].copy()
        if s2.n_obs < 40:
            results["strata"].append({"stratum": st, "status": "too_few_after_arm_filter", "n": int(s2.n_obs)})
            continue
        X = h.obsm["X_pca"][keep.to_numpy()]
        y = (s2.obs["arm"].astype(str) == "MCA_isocortex").astype(int).to_numpy()
        donor = s2.obs["subject"].astype(str).to_numpy() if "subject" in s2.obs else s2.obs_names.astype(str)
        c_auc, c_sd = cell_cv(X, y)
        d_auc, _ = donor_oof(X, y, donor)
        rec = {
            "stratum": st,
            "contrast": "MCA_isocortex_vs_PFC_same_territory_null",
            "n": int(s2.n_obs),
            "n_mca_iso": int((y == 1).sum()),
            "n_pfc": int((y == 0).sum()),
            "cell_cv_auroc": c_auc,
            "cell_cv_std": c_sd,
            "donor_oof_auroc": d_auc,
        }
        results["strata"].append(rec)
        log(f"{st}: cell={c_auc} donor={d_auc} n={s2.n_obs}")

    results["verdict"] = "NULL_REPLICATION_ZONATION_STRATIFIED"
    results["claim"] = (
        "This is a same-territory region null (both MCA-leaning). "
        "It does not test ACA vs MCA. Use only to show residual region effects after zonation strata."
    )
    OUT_JSON.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
    lines = [
        "# ROSMAP zonation-stratified same-territory null",
        "",
        f"**Verdict:** `{results['verdict']}`",
        "",
        results["claim"],
        "",
        "| Stratum | n | Cell CV | Donor OOF |",
        "|---|---:|---:|---:|",
    ]
    for r in results["strata"]:
        if "cell_cv_auroc" not in r:
            lines.append(f"| {r['stratum']} | {r.get('n','')} | — | — |")
            continue
        lines.append(
            f"| {r['stratum']} | {r['n']} | {r['cell_cv_auroc']:.3f} | "
            f"{r['donor_oof_auroc'] if r['donor_oof_auroc']==r['donor_oof_auroc'] else float('nan'):.3f} |"
        )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
