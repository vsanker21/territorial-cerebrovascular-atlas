"""
True Class B path: external ACA cortex endo vs Garcia MCA endo (tight).

Arms
  ACA: Linnarsson A24/A32/A25/A23 annotated endo + LIBD dACC endothelial_seed
  MCA: Garcia temporal + ROSMAP MCA tight endo (gate v2b)

Analyses
  1. Harmony integrate (batch = source dataset)
  2. Cell-level 5-fold CV AUROC on Harmony PCs
  3. Donor-level OOF AUROC (leave-one-donor mean score)
  4. Wilcoxon DE (caveat: cross-dataset)

Outputs under results/class_b/ and data/processed/class_b/
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
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[2]
LIN = ROOT / "data/processed/aca_cohort_hunt/linnarsson_aca_vascular.h5ad"
LIBD_GLOBUS = ROOT / "data/processed/aca_cohort_hunt/libd_dacc_globus_vascular.h5ad"
LIBD_AZ = ROOT / "data/processed/aca_cohort_hunt/libd_dacc_azimuth_vascular.h5ad"
LIBD_SEED = ROOT / "data/processed/aca_cohort_hunt/libd_dacc_vascular_seed.h5ad"
ALLEN = ROOT / "data/processed/aca_cohort_hunt/allen_acc_vascular.h5ad"
GARCIA_TEMP = ROOT / "data/processed/garcia/garcia_temporal_mca_qc_tight.h5ad"
GARCIA_ROSMAP = ROOT / "data/processed/garcia/garcia_rosmap_vascular_tight.h5ad"

OUT_DIR = ROOT / "results/class_b"
OUT_H5AD = ROOT / "data/processed/class_b/aca_mca_harmony_endo.h5ad"
OUT_ACA_POOL = ROOT / "data/processed/aca_cohort_hunt/aca_endo_pool.h5ad"
OUT_DE = OUT_DIR / "aca_mca_harmony_class_b_de.csv"
OUT_DONOR = OUT_DIR / "aca_mca_harmony_donor_aware.csv"
OUT_MD = OUT_DIR / "aca_mca_harmony_class_b_report.md"
OUT_SUMMARY = ROOT / "results/analysis_progress_summary.md"
OUT_JSON = OUT_DIR / "aca_mca_harmony_summary.json"

NEURONAL_CAPS = {"LINGO1": 0.5, "RBFOX1": 0.5, "KCNIP4": 0.5, "SNAP25": 0.5}
ENDO_MARKERS = ["CLDN5", "PECAM1", "FLT1", "VWF", "CDH5", "EGFL7"]


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def dense(X):
    return X.toarray() if sparse.issparse(X) else np.asarray(X)


def ensure_lognorm(a: ad.AnnData) -> ad.AnnData:
    """If X looks like raw counts, normalize+log1p; keep counts layer."""
    xmax = float(a.X.max()) if a.n_obs else 0.0
    if xmax > 20:
        if "counts" not in a.layers:
            a.layers["counts"] = a.X.copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        log(f"  log-normalized (was raw, xmax={xmax:.1f})")
    return a


def to_gene_symbols(a: ad.AnnData) -> ad.AnnData:
    """Map Ensembl var_names -> gene symbols when available (Linnarsson CELLxGENE)."""
    a = a.copy()
    symbol_col = None
    for c in ("feature_name", "Gene", "gene_symbols", "gene_symbol"):
        if c in a.var.columns:
            symbol_col = c
            break
    if symbol_col is None:
        # already symbols if CLDN5 present
        return a
    syms = a.var[symbol_col].astype(str).to_numpy()
    # drop missing / ensembl leftovers
    ok = np.array([s not in {"nan", "NA", ""} and not s.startswith("ENSG") for s in syms])
    a = a[:, ok].copy()
    syms = a.var[symbol_col].astype(str).to_numpy()
    # make unique
    a.var_names = pd.Index(syms).astype(str)
    a.var_names_make_unique()
    log(f"  mapped to gene symbols via {symbol_col}; n_vars={a.n_vars}")
    return a


def neuronal_cap_ok(a: ad.AnnData) -> np.ndarray:
    present = [g for g in NEURONAL_CAPS if g in a.var_names]
    if not present:
        return np.ones(a.n_obs, dtype=bool)
    X = dense(a[:, present].X)
    ok = np.ones(a.n_obs, dtype=bool)
    for i, g in enumerate(present):
        ok &= X[:, i] <= NEURONAL_CAPS[g]
    return ok


def endo_marker_ok(a: ad.AnnData) -> np.ndarray:
    present = [g for g in ENDO_MARKERS if g in a.var_names]
    if not present:
        return np.ones(a.n_obs, dtype=bool)
    X = dense(a[:, present].X)
    return (X > 0).any(axis=1)


def prep_arm(
    a: ad.AnnData,
    arm: str,
    batch: str,
    donor_col: str,
    *,
    neuronal_cap: float | None = 0.5,
) -> ad.AnnData:
    a = a.copy()
    a = ensure_lognorm(a)
    keep = endo_marker_ok(a)
    if neuronal_cap is not None:
        present = [g for g in NEURONAL_CAPS if g in a.var_names]
        if present:
            X = dense(a[:, present].X)
            neur_ok = np.ones(a.n_obs, dtype=bool)
            for i, _g in enumerate(present):
                neur_ok &= X[:, i] <= neuronal_cap
            keep &= neur_ok
    a = a[keep].copy()
    a.obs["class_b_arm"] = arm
    a.obs["batch"] = batch
    a.obs["donor"] = a.obs[donor_col].astype(str)
    a.obs["donor_key"] = a.obs["batch"].astype(str) + "::" + a.obs["donor"].astype(str)
    return a


def _prep_libd_seed(libd: ad.AnnData) -> ad.AnnData:
    """Marker-seed fallback — do not treat as published annotations."""
    libd = libd[libd.obs["cell_type_seed"].astype(str) == "endothelial_seed"].copy()
    libd = ensure_lognorm(libd)
    keep = endo_marker_ok(libd)
    present = [g for g in NEURONAL_CAPS if g in libd.var_names]
    if present:
        X = dense(libd[:, present].X)
        neur_ok = np.ones(libd.n_obs, dtype=bool)
        for i, g in enumerate(present):
            neur_ok &= X[:, i] <= 1.0
        keep &= neur_ok
    libd = libd[keep].copy()
    libd.obs["class_b_arm"] = "ACA"
    libd.obs["batch"] = "LIBD_dACC_seed"
    libd.obs["donor"] = libd.obs["donor_id"].astype(str)
    libd.obs["donor_key"] = libd.obs["batch"].astype(str) + "::" + libd.obs["donor"].astype(str)
    return libd


def load_aca_endo() -> ad.AnnData:
    parts: list[ad.AnnData] = []

    lin = sc.read_h5ad(LIN)
    lin = to_gene_symbols(lin)
    lin = lin[lin.obs["cell_type"].astype(str) == "endothelial cell"].copy()
    # Published taxonomy; mild ambient cap (unsorted cortex)
    lin = prep_arm(lin, "ACA", "Linnarsson", "donor_id", neuronal_cap=1.0)
    log(f"Linnarsson ACA endo after caps: {lin.n_obs}")
    parts.append(lin)

    if ALLEN.exists():
        al = sc.read_h5ad(ALLEN)
        al = to_gene_symbols(al)
        ct = al.obs["cell_type"].astype(str).str.lower()
        al = al[ct.str.contains("endo")].copy()
        donor_col = "donor_id" if "donor_id" in al.obs else al.obs.columns[0]
        # Published Allen Endo: ambient neuronal RNA is common in unsorted snRNA;
        # do not apply neuronal caps that discard author-labeled Endo.
        al = prep_arm(al, "ACA", "Allen_ACC", donor_col, neuronal_cap=None)
        log(f"Allen ACC endo after caps: {al.n_obs}")
        parts.append(al)

    # Prefer paper-identical GLOBUS sce_Final > method-matched Azimuth > marker seed
    if LIBD_GLOBUS.exists():
        libd = sc.read_h5ad(LIBD_GLOBUS)
        ct = libd.obs["cell_type"].astype(str).str.lower()
        libd = libd[ct.str.contains("endo")].copy()
        libd = prep_arm(libd, "ACA", "LIBD_dACC_globus", "donor_id", neuronal_cap=None)
        log(f"LIBD dACC GLOBUS published Endo after caps: {libd.n_obs}")
        parts.append(libd)
    elif LIBD_AZ.exists():
        libd = sc.read_h5ad(LIBD_AZ)
        ct = libd.obs["cell_type"].astype(str).str.lower()
        libd = libd[ct.str.contains("endo")].copy()
        libd = prep_arm(libd, "ACA", "LIBD_dACC_azimuth", "donor_id", neuronal_cap=None)
        log(f"LIBD dACC Azimuth Endo after caps: {libd.n_obs}")
        parts.append(libd)
    elif LIBD_SEED.exists():
        libd = sc.read_h5ad(LIBD_SEED)
        libd = _prep_libd_seed(libd)
        log(f"LIBD dACC endo seed (FALLBACK) after soft caps: {libd.n_obs}")
        parts.append(libd)
    else:
        log("WARN: no LIBD ACA object found")

    shared = parts[0].var_names
    for p in parts[1:]:
        shared = shared.intersection(p.var_names)
    if len(shared) < 1000:
        raise RuntimeError(f"ACA pool gene intersection too small: {len(shared)}")
    pool = ad.concat([p[:, shared].copy() for p in parts], join="inner", merge="same", index_unique="-")
    pool.obs["cohort_arm"] = "ACA_cortex"
    pool.obs["tv_dominant"] = "ACA"
    log(f"ACA endo pool total: {pool.n_obs} donors={pool.obs['donor_key'].nunique()}")
    return pool


def load_mca_endo() -> ad.AnnData:
    gt = sc.read_h5ad(GARCIA_TEMP)
    # PatientID is donor
    donor_col = "PatientID" if "PatientID" in gt.obs else "SampleID"
    gt = prep_arm(gt, "MCA", "Garcia_exvivo", donor_col)
    log(f"Garcia temporal MCA endo after caps: {gt.n_obs}")

    gr = sc.read_h5ad(GARCIA_ROSMAP)
    donor_col = "subject" if "subject" in gr.obs else "library_id"
    gr = prep_arm(gr, "MCA", "Garcia_ROSMAP", donor_col)
    log(f"Garcia ROSMAP MCA endo after caps: {gr.n_obs}")

    shared = gt.var_names.intersection(gr.var_names)
    return ad.concat(
        [gt[:, shared].copy(), gr[:, shared].copy()],
        join="inner",
        merge="same",
        index_unique="-",
    )


def harmony_integrate(comb: ad.AnnData) -> ad.AnnData:
    sc.pp.highly_variable_genes(comb, n_top_genes=2000, flavor="seurat", subset=False)
    hvg = comb[:, comb.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=30, svd_solver="arpack")
    log("Harmony integrate")
    sc.external.pp.harmony_integrate(
        hvg, key="batch", basis="X_pca", adjusted_basis="X_pca_harmony"
    )
    sc.pp.neighbors(hvg, use_rep="X_pca_harmony", n_neighbors=15)
    sc.tl.umap(hvg)
    comb.obsm["X_pca_harmony"] = hvg.obsm["X_pca_harmony"]
    comb.obsm["X_umap"] = hvg.obsm["X_umap"]
    return comb


def donor_aware(X: np.ndarray, y: np.ndarray, donor: np.ndarray, batch: np.ndarray):
    rows = []
    # Cell CV
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    cell_aucs = []
    for i, (tr, te) in enumerate(skf.split(X, y)):
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], y[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        auc = roc_auc_score(y[te], proba)
        cell_aucs.append(auc)
        rows.append({"fold_type": "cell_cv", "held_out": f"fold{i}", "auroc": auc, "n_test": len(te)})

    donors = pd.Series(donor)
    vc = donors.value_counts()
    usable = [d for d in vc.index if vc[d] >= 3]
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
                "fold_type": "leave_one_donor",
                "held_out": d,
                "mean_mca_score": float(proba.mean()),
                "true_label": int(y[te][0]),
                "n_test": int(te.sum()),
                "batch": batch[te][0],
                "auroc": float(roc_auc_score(y[te], proba)) if len(np.unique(y[te])) > 1 else np.nan,
            }
        )

    drows = []
    for d in usable:
        m = donor == d
        if np.isnan(oof[m]).all():
            continue
        drows.append(
            {
                "fold_type": "donor_oof_score",
                "held_out": d,
                "mean_mca_score": float(np.nanmean(oof[m])),
                "true_label": int(y[m][0]),
                "n_test": int(m.sum()),
                "batch": batch[m][0],
            }
        )
    rows.extend(drows)
    ddf = pd.DataFrame(drows)
    donor_auc = np.nan
    if len(ddf) >= 4 and ddf["true_label"].nunique() > 1:
        donor_auc = float(roc_auc_score(ddf["true_label"], ddf["mean_mca_score"]))

    # Leave-one-batch-out (dataset)
    batch_aucs = []
    for b in np.unique(batch):
        te = batch == b
        tr = ~te
        if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
            # single-arm batch: score only
            if len(np.unique(y[tr])) < 2:
                continue
            clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
            clf.fit(X[tr], y[tr])
            proba = clf.predict_proba(X[te])[:, 1]
            rows.append(
                {
                    "fold_type": "leave_one_batch",
                    "held_out": b,
                    "mean_mca_score": float(proba.mean()),
                    "true_label": int(y[te][0]),
                    "n_test": int(te.sum()),
                    "auroc": np.nan,
                    "note": "single-class batch",
                }
            )
            continue
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], y[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        auc = roc_auc_score(y[te], proba)
        batch_aucs.append(auc)
        rows.append(
            {
                "fold_type": "leave_one_batch",
                "held_out": b,
                "auroc": auc,
                "n_test": int(te.sum()),
            }
        )

    summary = {
        "cell_cv_auroc_mean": float(np.mean(cell_aucs)) if cell_aucs else np.nan,
        "cell_cv_auroc_std": float(np.std(cell_aucs)) if cell_aucs else np.nan,
        "donor_level_oof_auroc": donor_auc,
        "n_donors_scored": int(len(ddf)),
        "n_aca_donors": int((ddf["true_label"] == 0).sum()) if len(ddf) else 0,
        "n_mca_donors": int((ddf["true_label"] == 1).sum()) if len(ddf) else 0,
        "batch_lobo_auroc_mean": float(np.mean(batch_aucs)) if batch_aucs else np.nan,
        "n_batch_lobo": len(batch_aucs),
    }
    return pd.DataFrame(rows), summary


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    OUT_ACA_POOL.parent.mkdir(parents=True, exist_ok=True)

    log("Load ACA endo pool")
    aca = load_aca_endo()
    aca.write_h5ad(OUT_ACA_POOL, compression="gzip")
    log(f"Wrote ACA pool n={aca.n_obs} donors={aca.obs['donor_key'].nunique()}")

    log("Load MCA endo (Garcia tight)")
    mca = load_mca_endo()
    log(f"MCA endo n={mca.n_obs} donors={mca.obs['donor_key'].nunique()}")

    shared = aca.var_names.intersection(mca.var_names)
    log(f"Shared genes: {len(shared)}")
    if len(shared) < 1000:
        raise RuntimeError(f"ACA∩MCA gene intersection too small: {len(shared)}")
    aca = aca[:, shared].copy()
    mca = mca[:, shared].copy()
    comb = ad.concat([aca, mca], join="inner", merge="same", index_unique="-")
    comb.obs["class_b_arm"] = comb.obs["class_b_arm"].astype("category")

    comb = harmony_integrate(comb)

    # DE
    sc.tl.rank_genes_groups(
        comb,
        groupby="class_b_arm",
        groups=["MCA"],
        reference="ACA",
        method="wilcoxon",
        n_genes=200,
    )
    de = sc.get.rank_genes_groups_df(comb, group="MCA")
    de["contrast"] = "MCA_vs_ACA_harmony"
    de["caveat"] = "cross_dataset_Harmony; neuronal caps applied; not vessel-matched"
    de.to_csv(OUT_DE, index=False)

    X = comb.obsm["X_pca_harmony"]
    y = (comb.obs["class_b_arm"].astype(str) == "MCA").astype(int).to_numpy()
    donor = comb.obs["donor_key"].astype(str).to_numpy()
    batch = comb.obs["batch"].astype(str).to_numpy()
    donor_df, summ = donor_aware(X, y, donor, batch)
    donor_df.to_csv(OUT_DONOR, index=False)

    comb.write_h5ad(OUT_H5AD, compression="gzip")

    # Neuronal DE contamination check
    neur = ["LINGO1", "RBFOX1", "KCNIP4", "SNAP25", "SYT1", "RBFOX3"]
    neur_hits = []
    for g in neur:
        hit = de[de["names"] == g]
        if len(hit):
            r = hit.iloc[0]
            neur_hits.append(
                f"{g}: rank={int(hit.index[0])+1}, logFC={r['logfoldchanges']:.2f}, padj={r['pvals_adj']:.2e}"
            )
    neur_txt = "\n".join(f"- {h}" for h in neur_hits) if neur_hits else "- none in top 200"

    top = de.head(15)
    bot = de.tail(10)
    counts = {
        "n_aca": int((comb.obs["class_b_arm"] == "ACA").sum()),
        "n_mca": int((comb.obs["class_b_arm"] == "MCA").sum()),
        "n_aca_donors": int(comb.obs.loc[comb.obs["class_b_arm"] == "ACA", "donor_key"].nunique()),
        "n_mca_donors": int(comb.obs.loc[comb.obs["class_b_arm"] == "MCA", "donor_key"].nunique()),
        "by_batch": comb.obs["batch"].astype(str).value_counts().to_dict(),
    }
    summ.update(counts)

    # Verdict
    cell_auc = summ["cell_cv_auroc_mean"]
    donor_auc = summ["donor_level_oof_auroc"]
    if np.isnan(donor_auc):
        verdict = "INCONCLUSIVE_donor_AUROC_undefined"
    elif donor_auc >= 0.7 and cell_auc >= 0.7:
        verdict = "SUPPORTIVE_transferable_ClassB_signal"
    elif cell_auc >= 0.75 and donor_auc < 0.6:
        verdict = "CELL_LEVEL_ONLY_likely_dataset_or_donor_leakage"
    elif donor_auc < 0.55:
        verdict = "NO_TRANSFERABLE_ClassB_signal_on_current_n"
    else:
        verdict = "WEAK_MIXED_signal"
    summ["verdict"] = verdict

    OUT_JSON.write_text(json.dumps(summ, indent=2), encoding="utf-8")

    md = f"""# Class B — true ACA cortex vs Garcia MCA (Harmony)

First TerritoryVoxels-mapped **ACA isocortex vascular** arm vs **MCA** endo.

## Arms

| Arm | Sources | N endo (after caps) | N donors |
|---|---|---:|---:|
| **ACA** | Linnarsson A24/A32/A25/A23 + LIBD dACC | {counts['n_aca']} | {counts['n_aca_donors']} |
| **MCA** | Garcia temporal + ROSMAP MCA tight v2b | {counts['n_mca']} | {counts['n_mca_donors']} |

Batch keys: `{counts['by_batch']}`

Neuronal caps: LINGO1/RBFOX1/KCNIP4/SNAP25 ≤ 0.5; require ≥1 endo marker.

## Separability

| Metric | Value |
|---|---:|
| Cell 5-fold CV AUROC | {summ['cell_cv_auroc_mean']:.3f} ± {summ['cell_cv_auroc_std']:.3f} |
| **Donor-level OOF AUROC** | **{summ['donor_level_oof_auroc'] if not np.isnan(summ['donor_level_oof_auroc']) else 'NA'}** |
| Donors scored (ACA / MCA) | {summ['n_aca_donors']} / {summ['n_mca_donors']} |
| Batch LOBO AUROC (multi-class only) | {summ['batch_lobo_auroc_mean'] if not np.isnan(summ['batch_lobo_auroc_mean']) else 'NA'} |

**Verdict:** `{verdict}`

## Caveats

- Still **cross-dataset** (tech + processing confound); Harmony reduces but does not erase this.
- ACA arm is **not vessel-enriched**; endo n remains modest.
- Yang SFG is **not** used as ACA (borderzone).

## Neuronal contaminant check in DE top-200

{neur_txt}

## Top genes higher in MCA

{top[['names','scores','logfoldchanges','pvals_adj']].to_string(index=False)}

## Higher in ACA (bottom ranks)

{bot[['names','scores','logfoldchanges','pvals_adj']].to_string(index=False)}

## Objects

- `{OUT_ACA_POOL.relative_to(ROOT).as_posix()}`
- `{OUT_H5AD.relative_to(ROOT).as_posix()}`
- DE: `{OUT_DE.relative_to(ROOT).as_posix()}`
- Donor folds: `{OUT_DONOR.relative_to(ROOT).as_posix()}`
"""
    OUT_MD.write_text(md, encoding="utf-8")
    log(f"Wrote {OUT_MD}")
    log(json.dumps(summ, indent=2))

    # Patch progress summary
    prev = OUT_SUMMARY.read_text(encoding="utf-8") if OUT_SUMMARY.exists() else ""
    block = (
        f"\n## Class B true ACA vs MCA (this pass)\n"
        f"- Report: `{OUT_MD.relative_to(ROOT).as_posix()}`\n"
        f"- Cell CV AUROC {summ['cell_cv_auroc_mean']:.3f}; "
        f"donor OOF AUROC {summ['donor_level_oof_auroc']}; "
        f"verdict `{verdict}`\n"
        f"- N ACA/MCA endo: {counts['n_aca']}/{counts['n_mca']}\n"
    )
    if "## Class B true ACA vs MCA" not in prev:
        OUT_SUMMARY.write_text(prev.rstrip() + "\n" + block, encoding="utf-8")


if __name__ == "__main__":
    main()
