"""Moonshot #2 — Counterfactual enrichment transport (SENSITIVITY, not H1c).

Learn BVE treatment from Matthews paired parenchyma vs vascular fractions;
transport unsorted ACC endothelium into pseudo-BVE latent space; compare to
Garcia MCA BVE. Never report gene markers as confirmatory territorial DE.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import anndata as ad
import anndata as ad_mod
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

from _paths import ALLEN_ACC, CLASS_B, GARCIA_TEMPORAL, GARCIA_TEMPORAL_FALLBACK, HONESTY, LINNARSSON, MATTHEWS_H5, RES

OUT = RES / "02_enrichment_transport"
OUT.mkdir(parents=True, exist_ok=True)

ENDO_MARKERS = ["CLDN5", "PECAM1", "FLT1", "KDR", "CDH5", "VWF", "EGFL7", "ECSCR"]


def _read_matthews() -> ad.AnnData:
    files = sorted(MATTHEWS_H5.glob("*_filtered_feature_bc_matrix.h5"))
    ads = []
    for fp in files:
        m = re.search(r"(GSM\d+)_(P|V)_(\d+)", fp.name)
        if not m:
            continue
        gsm, frac, num = m.group(1), m.group(2), m.group(3)
        a = sc.read_10x_h5(fp)
        a.var_names_make_unique()
        a.obs["gsm"] = gsm
        a.obs["fraction"] = "parenchyma" if frac == "P" else "vascular"
        a.obs["pair_id"] = num
        a.obs["sample"] = f"{frac}_{num}"
        ads.append(a)
    combined = ad_mod.concat(ads, join="outer", index_unique="-")
    return combined


def _score_endo(a: ad.AnnData) -> np.ndarray:
    genes = list(a.var_names.astype(str))
    idx = [genes.index(g) for g in ENDO_MARKERS if g in a.var_names]
    if not idx:
        return np.zeros(a.n_obs)
    X = a.X[:, idx]
    if sparse.issparse(X):
        return np.asarray(X.mean(axis=1)).ravel()
    return X.mean(axis=1)


def _to_dense_log(a: ad.AnnData, max_cells: int = 8000, rng=None) -> tuple[np.ndarray, ad.AnnData]:
    rng = rng or np.random.default_rng(0)
    if a.n_obs > max_cells:
        idx = rng.choice(a.n_obs, max_cells, replace=False)
        a = a[idx].copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat", subset=False)
    if "highly_variable" not in a.var.columns or int(a.var["highly_variable"].sum()) < 100:
        # variance fallback
        Xtmp = a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)
        v = Xtmp.var(axis=0)
        top = np.argsort(v)[-2000:]
        a = a[:, top].copy()
    else:
        a = a[:, a.var["highly_variable"]].copy()
    X = a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)
    return X, a


def _donor_oof(X, y, groups) -> float:
    mask = np.isfinite(y)
    X, y, groups = X[mask], y[mask], groups[mask]
    if len(np.unique(groups)) < 3:
        return float("nan")
    gkf = GroupKFold(n_splits=min(5, len(np.unique(groups))))
    preds = np.full(len(y), np.nan)
    for tr, te in gkf.split(X, y, groups):
        if len(np.unique(y[tr])) < 2:
            continue
        sc_ = StandardScaler()
        clf = LogisticRegression(max_iter=2000, class_weight="balanced")
        clf.fit(sc_.fit_transform(X[tr]), y[tr])
        preds[te] = clf.predict_proba(sc_.transform(X[te]))[:, 1]
    ok = ~np.isnan(preds)
    if ok.sum() < 10 or len(np.unique(y[ok])) < 2:
        return float("nan")
    return float(roc_auc_score(y[ok], preds[ok]))


def run() -> dict:
    rng = np.random.default_rng(42)
    print("Loading Matthews…")
    mat = _read_matthews()
    mat.obs["endo_score"] = _score_endo(mat)
    # Keep putative endothelium: top quintile within each sample or absolute
    thr = np.quantile(mat.obs["endo_score"], 0.80)
    endo = mat[mat.obs["endo_score"] >= thr].copy()
    print("Matthews endo-like cells", endo.n_obs, "thr", thr)

    # Shared PCA on Matthews endo
    X_m, endo_h = _to_dense_log(endo, max_cells=6000, rng=rng)
    pca = PCA(n_components=30, random_state=0)
    Z = pca.fit_transform(X_m)
    endo_h.obsm["X_pca"] = Z
    frac = endo_h.obs["fraction"].astype(str).values
    pair = endo_h.obs["pair_id"].astype(str).values

    # Learn mean transport parenchyma -> vascular in PCA, within pairs
    deltas = []
    for pid in np.unique(pair):
        pmask = (pair == pid) & (frac == "parenchyma")
        vmask = (pair == pid) & (frac == "vascular")
        if pmask.sum() < 5 or vmask.sum() < 5:
            continue
        deltas.append(Z[vmask].mean(axis=0) - Z[pmask].mean(axis=0))
    if not deltas:
        raise RuntimeError("No paired Matthews deltas")
    delta = np.mean(deltas, axis=0)
    # Also OT-lite: map each parenchyma cell to kNN vascular centroid shift
    P = Z[frac == "parenchyma"]
    V = Z[frac == "vascular"]
    nn = NearestNeighbors(n_neighbors=10).fit(V)
    _, idx = nn.kneighbors(P)
    ot_shift = (V[idx].mean(axis=1) - P).mean(axis=0)

    transport_stats = {
        "n_pairs_used": int(len(deltas)),
        "mean_delta_norm": float(np.linalg.norm(delta)),
        "ot_shift_norm": float(np.linalg.norm(ot_shift)),
        "delta_vs_ot_cosine": float(
            np.dot(delta, ot_shift) / (np.linalg.norm(delta) * np.linalg.norm(ot_shift) + 1e-12)
        ),
    }

    # Validation: can we classify V vs P better after? (should already be separable)
    y_vp = (frac == "vascular").astype(int)
    auroc_raw = _donor_oof(Z, y_vp.astype(float), pair)
    transport_stats["matthews_V_vs_P_donor_oof_auroc"] = auroc_raw

    # Load ACC unsorted endo + Garcia MCA
    acc_paths = [ALLEN_ACC, LINNARSSON]
    acc_ads = []
    for p in acc_paths:
        if p.exists():
            a = ad.read_h5ad(p)
            a.obs["acc_source"] = p.stem
            acc_ads.append(a)
    acc = ad_mod.concat(acc_ads, join="outer", index_unique="-") if acc_ads else None

    garcia_path = GARCIA_TEMPORAL if GARCIA_TEMPORAL.exists() else GARCIA_TEMPORAL_FALLBACK
    # Prefer class_b MCA arm + garcia tight if temporal too big
    results = {"honesty_lock": HONESTY["2"], "transport": transport_stats}

    if acc is None:
        results["error"] = "No ACC object"
        (OUT / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
        return results

    # Use Class B object: ACA vs MCA already harmony-aligned — apply transport in that PCA
    cb = ad.read_h5ad(CLASS_B)
    # X_pca_harmony exists
    if "X_pca_harmony" in cb.obsm:
        Zcb = np.asarray(cb.obsm["X_pca_harmony"])
    else:
        Xcb = cb.X.toarray() if sparse.issparse(cb.X) else np.asarray(cb.X)
        Zcb = PCA(n_components=min(30, Xcb.shape[1]), random_state=0).fit_transform(Xcb)

    # Map Matthews delta into Class B PCA dim by padding/truncating
    d = Zcb.shape[1]
    delta_cb = np.zeros(d)
    ot_cb = np.zeros(d)
    m = min(d, len(delta))
    delta_cb[:m] = delta[:m]
    ot_cb[:m] = ot_shift[:m]

    is_aca = cb.obs["tv_dominant"].astype(str).isin(["ACA"]).values
    is_mca = cb.obs["tv_dominant"].astype(str).isin(["MCA"]).values
    # Heuristic: class_b_arm may encode ACA sources
    if is_aca.sum() == 0 and "class_b_arm" in cb.obs.columns:
        is_aca = cb.obs["class_b_arm"].astype(str).str.contains("ACA|Allen|LIBD|Linn|ACC", case=False, na=False).values
        is_mca = ~is_aca

    Z_aca = Zcb[is_aca]
    Z_mca = Zcb[is_mca]
    Z_aca_t = Z_aca + delta_cb  # counterfactual transport
    Z_aca_ot = Z_aca + ot_cb

    def _sep(Za, Zm, label):
        X = np.vstack([Za, Zm])
        y = np.array([1] * len(Za) + [0] * len(Zm), dtype=float)
        # fake groups by chunking cells into pseudo-donors if needed
        if "donor" in cb.obs.columns:
            ga = cb.obs.loc[is_aca, "donor"].astype(str).values
            gm = cb.obs.loc[is_mca, "donor"].astype(str).values
            groups = np.concatenate([ga, gm])
        else:
            groups = np.array([f"g{i//20}" for i in range(len(y))])
        return _donor_oof(X, y, groups)

    sep_raw = _sep(Z_aca, Z_mca, "raw")
    sep_t = _sep(Z_aca_t, Z_mca, "transported")
    sep_ot = _sep(Z_aca_ot, Z_mca, "ot")

    # Distance to MCA centroid
    mca_c = Z_mca.mean(axis=0)
    dist_raw = float(np.linalg.norm(Z_aca.mean(0) - mca_c))
    dist_t = float(np.linalg.norm(Z_aca_t.mean(0) - mca_c))
    dist_ot = float(np.linalg.norm(Z_aca_ot.mean(0) - mca_c))

    results.update(
        {
            "n_aca_cells": int(is_aca.sum()),
            "n_mca_cells": int(is_mca.sum()),
            "donor_oof_auroc_raw_ACA_vs_MCA": sep_raw,
            "donor_oof_auroc_transported_mean_shift": sep_t,
            "donor_oof_auroc_transported_ot": sep_ot,
            "distance_ACA_to_MCA_centroid_raw": dist_raw,
            "distance_ACA_to_MCA_centroid_transported": dist_t,
            "distance_ACA_to_MCA_centroid_ot": dist_ot,
            "interpretation": (
                "SENSITIVITY ONLY. If transport reduces ACA–MCA distance / AUROC, "
                "enrichment batch may inflate separability; if AUROC remains high after "
                "transport toward BVE, residual territory-linked state is compatible with "
                "hypothesis — still NOT confirmatory H1c markers."
            ),
        }
    )

    pd.DataFrame(
        [
            {"method": "raw", "auroc": sep_raw, "dist_to_mca": dist_raw},
            {"method": "mean_shift_transport", "auroc": sep_t, "dist_to_mca": dist_t},
            {"method": "knn_ot_transport", "auroc": sep_ot, "dist_to_mca": dist_ot},
        ]
    ).to_csv(OUT / "transport_separability.csv", index=False)

    (OUT / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    (OUT / "report.md").write_text(
        "\n".join(
            [
                "# Moonshot 02 — Counterfactual enrichment transport",
                "",
                f"**Honesty lock:** {HONESTY['2']}",
                "",
                f"- Matthews pairs used: {transport_stats['n_pairs_used']}",
                f"- V vs P AUROC: {auroc_raw:.3f}",
                f"- Class B ACA vs MCA AUROC raw / transported / OT: "
                f"{sep_raw:.3f} / {sep_t:.3f} / {sep_ot:.3f}",
                f"- Distance to MCA centroid raw → transported: {dist_raw:.3f} → {dist_t:.3f}",
                "",
                results["interpretation"],
            ]
        ),
        encoding="utf-8",
    )
    print("DONE 02", results)
    return results


if __name__ == "__main__":
    run()
