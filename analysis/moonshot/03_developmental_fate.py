"""Moonshot #3 — Developmental fate memory without HOX RNA.

Align Wälchli fetal CNS EC with adult territorial endothelia in shared PCA/Harmony
space. Falsify Track B if adult territory is orthogonal to developmental state
after anatomy covariates.
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from _paths import CLASS_B, HONESTY, RES

OUT = RES / "03_developmental_fate"
OUT.mkdir(parents=True, exist_ok=True)

FETAL_MTX = Path(r"F:\Cerebrovascular SCS\data\processed\walchli\fetal_mtx")
ADULT_MTX = Path(r"F:\Cerebrovascular SCS\data\processed\walchli\adult_temporal_mtx")


def _read_mtx(folder: Path, label: str) -> ad.AnnData:
    a = sc.read_mtx(folder / "matrix.mtx").T
    genes = pd.read_csv(folder / "features.tsv", header=None)[0].astype(str).values
    bars = pd.read_csv(folder / "barcodes.tsv", header=None)[0].astype(str).values
    a.var_names = genes
    a.obs_names = bars
    a.var_names_make_unique()
    meta = pd.read_csv(folder / "meta.csv", index_col=0)
    # align meta
    meta.index = meta.index.astype(str)
    common = a.obs_names.intersection(meta.index)
    if len(common) == 0:
        # try barcode without prefix
        meta.index = bars[: len(meta)]
        common = a.obs_names.intersection(meta.index)
    a = a[common].copy()
    for c in meta.columns:
        a.obs[c] = meta.loc[a.obs_names, c].values
    a.obs["dev_stage"] = label
    return a


def _prep(a: ad.AnnData, n_cells: int, rng) -> ad.AnnData:
    if a.n_obs > n_cells:
        a = a[rng.choice(a.n_obs, n_cells, replace=False)].copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def run() -> dict:
    rng = np.random.default_rng(0)
    results = {"honesty_lock": HONESTY["3"]}

    if not FETAL_MTX.exists():
        results["status"] = "GATED_FETAL_MTX_MISSING"
        (OUT / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
        return results

    print("Loading fetal…")
    fetal = _read_mtx(FETAL_MTX, "fetal")
    fetal = _prep(fetal, 5000, rng)

    adult_parts = []
    if ADULT_MTX.exists():
        print("Loading adult Wälchli temporal…")
        adult_w = _read_mtx(ADULT_MTX, "adult_walchli_temporal_MCA")
        adult_w = _prep(adult_w, 5000, rng)
        adult_w.obs["territory"] = "MCA"
        adult_parts.append(adult_w)

    cb = ad.read_h5ad(CLASS_B)
    cb.obs["dev_stage"] = "adult_class_b"
    # territory
    if "tv_dominant" in cb.obs:
        cb.obs["territory"] = cb.obs["tv_dominant"].astype(str)
    elif "class_b_arm" in cb.obs:
        cb.obs["territory"] = cb.obs["class_b_arm"].astype(str)
    else:
        cb.obs["territory"] = "unknown"
    # Class B already log-normalized-ish; take expression matrix
    adult_parts.append(cb)

    # Ensure territory column exists on all adult parts
    for a in adult_parts:
        if "territory" not in a.obs.columns:
            a.obs["territory"] = "MCA" if "walchli" in str(a.obs.get("dev_stage", [""])[0]).lower() else "unknown"
    fetal.obs["territory"] = "developmental"

    # Intersect genes
    genes = set(fetal.var_names.astype(str))
    for a in adult_parts:
        genes &= set(a.var_names.astype(str))
    genes = sorted(genes)
    print("shared genes", len(genes))
    if len(genes) < 500:
        results["status"] = "TOO_FEW_SHARED_GENES"
        results["n_shared"] = len(genes)
        (OUT / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
        return results

    fetal_s = fetal[:, genes].copy()
    adults = [a[:, genes].copy() for a in adult_parts]
    combined = ad.concat([fetal_s] + adults, join="inner", index_unique="-")

    sc.pp.highly_variable_genes(combined, n_top_genes=2000, flavor="seurat", subset=True)
    X = combined.X.toarray() if sparse.issparse(combined.X) else np.asarray(combined.X)
    Z = PCA(n_components=30, random_state=0).fit_transform(StandardScaler().fit_transform(X))
    combined.obsm["X_pca"] = Z

    # Developmental axis: PC that best separates fetal vs adult
    y_fet = (combined.obs["dev_stage"].astype(str) == "fetal").astype(int).values
    # correlation of each PC with fetal label
    corrs = [abs(stats.pearsonr(Z[:, i], y_fet)[0]) for i in range(Z.shape[1])]
    dev_pc = int(np.argmax(corrs))
    combined.obs["dev_score"] = Z[:, dev_pc]
    results["developmental_pc"] = dev_pc
    results["developmental_pc_abs_corr"] = float(corrs[dev_pc])

    # Among adult cells: does territory associate with dev_score?
    adult_mask = combined.obs["dev_stage"].astype(str) != "fetal"
    adult = combined[adult_mask].copy()
    # Restrict to ACA/MCA if present
    terr = adult.obs["territory"].astype(str)
    keep = terr.isin(["ACA", "MCA"]) | terr.str.contains("ACA|MCA|Allen|Garcia|ACC", case=False, na=False)
    # Map
    def map_terr(t):
        t = str(t)
        if t == "ACA" or "ACA" in t or "Allen" in t or "LIBD" in t or "Linn" in t or "ACC" in t:
            return "ACA"
        if t == "MCA" or "MCA" in t or "Garcia" in t or "temporal" in t.lower():
            return "MCA"
        return "other"

    adult.obs["terr2"] = [map_terr(t) for t in adult.obs["territory"]]
    sub = adult[adult.obs["terr2"].isin(["ACA", "MCA"])].copy()
    if sub.n_obs < 50:
        results["status"] = "TOO_FEW_ADULT_TERRITORY_CELLS"
        (OUT / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
        return results

    # Unadjusted association
    y = (sub.obs["terr2"] == "ACA").astype(int).values
    r, p = stats.pointbiserialr(y, sub.obs["dev_score"].values)
    results["pointbiserial_terr_devscore_r"] = float(r)
    results["pointbiserial_terr_devscore_p"] = float(p)

    # Residualize dev_score on dataset batch
    batch = pd.get_dummies(sub.obs["dev_stage"].astype(str), drop_first=True)
    if batch.shape[1] > 0:
        lr = LinearRegression().fit(batch.values, sub.obs["dev_score"].values)
        resid = sub.obs["dev_score"].values - lr.predict(batch.values)
    else:
        resid = sub.obs["dev_score"].values
    r2, p2 = stats.pointbiserialr(y, resid)
    results["pointbiserial_terr_devscore_resid_batch_r"] = float(r2)
    results["pointbiserial_terr_devscore_resid_batch_p"] = float(p2)

    # Donor-OOF AUROC predicting territory from first 10 PCs (incl dev) vs without PC0/dev
    donor = sub.obs["donor"].astype(str).values if "donor" in sub.obs else np.array([f"c{i//15}" for i in range(sub.n_obs)])
    Zs = sub.obsm["X_pca"]

    def oof(X):
        gkf = GroupKFold(n_splits=min(5, len(np.unique(donor))))
        preds = np.full(len(y), np.nan)
        for tr, te in gkf.split(X, y, donor):
            if len(np.unique(y[tr])) < 2:
                continue
            sc_ = StandardScaler()
            clf = LogisticRegression(max_iter=2000, class_weight="balanced")
            clf.fit(sc_.fit_transform(X[tr]), y[tr])
            preds[te] = clf.predict_proba(sc_.transform(X[te]))[:, 1]
        ok = ~np.isnan(preds)
        return float(roc_auc_score(y[ok], preds[ok])) if ok.sum() > 10 and len(np.unique(y[ok])) > 1 else float("nan")

    auroc_all = oof(Zs[:, :15])
    # Remove developmental PC
    cols = [i for i in range(15) if i != dev_pc]
    auroc_nodev = oof(Zs[:, cols])
    results["donor_oof_auroc_with_devPC"] = auroc_all
    results["donor_oof_auroc_without_devPC"] = auroc_nodev
    results["delta_auroc_devPC"] = (
        float(auroc_all - auroc_nodev) if np.isfinite(auroc_all) and np.isfinite(auroc_nodev) else float("nan")
    )

    # Track B decision
    # Kill if residual association |r|<0.1 and delta_auroc < 0.02
    kill = (abs(r2) < 0.10 and (not np.isfinite(results["delta_auroc_devPC"]) or abs(results["delta_auroc_devPC"]) < 0.02))
    results["track_B_language"] = "KILL" if kill else "KEEP_GATED_EXPLORATORY"
    results["n_fetal"] = int(fetal.n_obs)
    results["n_adult_territory"] = int(sub.n_obs)

    pd.DataFrame(
        {
            "dev_score": sub.obs["dev_score"].values,
            "dev_score_resid_batch": resid,
            "territory": sub.obs["terr2"].values,
            "dev_stage": sub.obs["dev_stage"].values,
        }
    ).to_csv(OUT / "adult_devscore_by_territory.csv", index=False)

    (OUT / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    (OUT / "report.md").write_text(
        "\n".join(
            [
                "# Moonshot 03 — Developmental fate memory",
                "",
                f"**Honesty lock:** {HONESTY['3']}",
                "",
                f"- Developmental PC: {dev_pc} (|corr|={results['developmental_pc_abs_corr']:.3f})",
                f"- Territory↔dev_score r (raw / batch-resid): {r:.3f} / {r2:.3f}",
                f"- Donor OOF AUROC with/without dev PC: {auroc_all:.3f} / {auroc_nodev:.3f}",
                f"- **Track B language: {results['track_B_language']}**",
            ]
        ),
        encoding="utf-8",
    )
    print("DONE 03", results)
    return results


if __name__ == "__main__":
    run()
