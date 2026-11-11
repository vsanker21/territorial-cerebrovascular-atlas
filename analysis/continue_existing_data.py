"""
Continue analyses on existing processed objects (no ATLAS/RADC/ADNI).

1. Yang Aim3 DE stratified by vascular lineage (endo, peri)
2. H4 first-pass: territory classifier trained on control cells → scored on AD
3. Class B exploratory: Garcia MCA endo vs Yang SFG borderzone endo (controls)

Outputs under results/yang/, results/h4/, results/class_b/
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
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[1]
YANG = ROOT / "data/processed/yang/yang_vascular_qc.h5ad"
GARCIA = ROOT / "data/processed/garcia/garcia_temporal_mca_qc.h5ad"
OUT_AIM3 = ROOT / "results/yang/aim3_lineage_stratified_de.csv"
OUT_H4 = ROOT / "results/h4/yang_territory_classifier_h4.csv"
OUT_H4_MD = ROOT / "results/h4/yang_territory_classifier_h4_report.md"
OUT_CB = ROOT / "results/class_b/garcia_mca_vs_yang_sfg_endo_de.csv"
OUT_CB_MD = ROOT / "results/class_b/garcia_mca_vs_yang_sfg_report.md"
OUT_SUMMARY = ROOT / "results/analysis_progress_summary.md"


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def subsample(adata: ad.AnnData, groupby: str, groups: list[str], n_max: int, seed: int = 0) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    parts = []
    for g in groups:
        idx = adata.obs_names[adata.obs[groupby].astype(str) == g]
        if len(idx) == 0:
            continue
        if len(idx) > n_max:
            idx = rng.choice(idx.to_numpy(), size=n_max, replace=False)
        parts.append(adata[idx])
    return ad.concat(parts, join="inner", merge="same")


def run_aim3_lineage(yang: ad.AnnData) -> pd.DataFrame:
    frames = []
    for lineage in ["endo", "peri"]:
        for disease in ["control", "AD"]:
            sub = yang[
                (yang.obs["vascular_lineage"].astype(str) == lineage)
                & (yang.obs["disease"].astype(str) == disease)
            ].copy()
            if sub.obs["territory_contrast"].nunique() < 2:
                continue
            sub = subsample(
                sub,
                "territory_contrast",
                ["PCA_hippocampus", "ACA_MCA_borderzone_SFG"],
                n_max=6000,
                seed=hash((lineage, disease)) % 10_000,
            )
            log(f"Aim3 DE {lineage}/{disease}: {sub.n_obs} cells")
            sc.tl.rank_genes_groups(
                sub,
                groupby="territory_contrast",
                groups=["ACA_MCA_borderzone_SFG"],
                reference="PCA_hippocampus",
                method="wilcoxon",
                n_genes=150,
            )
            df = sc.get.rank_genes_groups_df(sub, group="ACA_MCA_borderzone_SFG")
            df["lineage"] = lineage
            df["disease"] = disease
            df["contrast"] = "SFG_borderzone_vs_PCA_hippocampus"
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def highly_variable_genes(adata: ad.AnnData, n_top: int = 2000) -> list[str]:
    # adata already log-normalized
    a = adata.copy()
    sc.pp.highly_variable_genes(a, n_top_genes=n_top, flavor="seurat", subset=False)
    return a.var_names[a.var["highly_variable"]].tolist()


def run_h4_classifier(yang: ad.AnnData) -> tuple[pd.DataFrame, dict]:
    """
    Train territory classifier on control endo; evaluate within control (CV)
    and transfer to AD (same features). H4 expects lower AD performance /
    shifted decision scores if territorial identity degrades.
    """
    endo = yang[yang.obs["vascular_lineage"].astype(str) == "endo"].copy()
    ctrl = endo[endo.obs["disease"].astype(str) == "control"].copy()
    ad_ = endo[endo.obs["disease"].astype(str) == "AD"].copy()
    ctrl = subsample(
        ctrl,
        "territory_contrast",
        ["PCA_hippocampus", "ACA_MCA_borderzone_SFG"],
        n_max=8000,
        seed=1,
    )
    ad_ = subsample(
        ad_,
        "territory_contrast",
        ["PCA_hippocampus", "ACA_MCA_borderzone_SFG"],
        n_max=8000,
        seed=2,
    )

    genes = highly_variable_genes(ctrl, n_top=1500)
    # intersect with AD
    genes = [g for g in genes if g in ad_.var_names]
    log(f"H4 features: {len(genes)} HVGs")

    def xy(a: ad.AnnData):
        X = a[:, genes].X
        if hasattr(X, "toarray"):
            X = X.toarray()
        y = (a.obs["territory_contrast"].astype(str) == "ACA_MCA_borderzone_SFG").astype(int).to_numpy()
        return X, y

    Xc, yc = xy(ctrl)
    Xa, ya = xy(ad_)

    # Control CV
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    cv_rows = []
    oof = np.zeros(len(yc))
    for fold, (tr, te) in enumerate(skf.split(Xc, yc)):
        scaler = StandardScaler(with_mean=True)
        Xtr = scaler.fit_transform(Xc[tr])
        Xte = scaler.transform(Xc[te])
        clf = LogisticRegression(max_iter=500, C=0.1, solver="lbfgs")
        clf.fit(Xtr, yc[tr])
        proba = clf.predict_proba(Xte)[:, 1]
        oof[te] = proba
        cv_rows.append(
            {
                "split": f"control_cv_fold{fold}",
                "auroc": float(roc_auc_score(yc[te], proba)),
                "accuracy": float(accuracy_score(yc[te], (proba >= 0.5).astype(int))),
                "n": int(len(te)),
            }
        )

    # Fit on all controls → score AD
    scaler = StandardScaler()
    Xc_s = scaler.fit_transform(Xc)
    clf = LogisticRegression(max_iter=500, C=0.1, solver="lbfgs")
    clf.fit(Xc_s, yc)
    pa = clf.predict_proba(scaler.transform(Xa))[:, 1]
    ad_auroc = float(roc_auc_score(ya, pa))
    ad_acc = float(accuracy_score(ya, (pa >= 0.5).astype(int)))
    ctrl_auroc = float(roc_auc_score(yc, oof))
    ctrl_acc = float(accuracy_score(yc, (oof >= 0.5).astype(int)))

    # Per-arm mean score in AD vs control (decision confidence for true class)
    def mean_correct_score(y, p):
        # score toward true label
        return float(np.where(y == 1, p, 1 - p).mean())

    summary = {
        "control_cv_auroc": ctrl_auroc,
        "control_cv_accuracy": ctrl_acc,
        "ad_transfer_auroc": ad_auroc,
        "ad_transfer_accuracy": ad_acc,
        "delta_auroc_ad_minus_control": ad_auroc - ctrl_auroc,
        "control_mean_correct_score": mean_correct_score(yc, oof),
        "ad_mean_correct_score": mean_correct_score(ya, pa),
        "n_control": int(len(yc)),
        "n_ad": int(len(ya)),
        "n_features": len(genes),
        "lineage": "endo",
        "classes": "ACA_MCA_borderzone_SFG vs PCA_hippocampus",
        "interpretation": (
            "H4 predicts AD transfer AUROC < control CV AUROC "
            "(territorial identity degradation). First-pass exploratory."
        ),
    }
    cv_df = pd.DataFrame(cv_rows)
    cv_df.to_csv(ROOT / "results/h4/yang_h4_cv_folds.csv", index=False)
    return pd.DataFrame([summary]), summary


def run_class_b(yang: ad.AnnData, garcia: ad.AnnData) -> pd.DataFrame:
    """
    Exploratory cross-dataset DE: Garcia MCA temporal endo (controls/epilepsy)
    vs Yang SFG borderzone endo controls.
    Label dataset as group; report as exploratory pending batch harmonization.
    """
    y = yang[
        (yang.obs["vascular_lineage"].astype(str) == "endo")
        & (yang.obs["disease"].astype(str) == "control")
        & (yang.obs["territory_contrast"].astype(str) == "ACA_MCA_borderzone_SFG")
    ].copy()
    g = garcia[garcia.obs["vascular_lineage"].astype(str) == "endo"].copy()
    y.obs["class_b_arm"] = "Yang_SFG_borderzone"
    g.obs["class_b_arm"] = "Garcia_MCA_temporal"
    # Harmonize to shared genes
    shared = y.var_names.intersection(g.var_names)
    y = y[:, shared].copy()
    g = g[:, shared].copy()
    y = subsample(y, "class_b_arm", ["Yang_SFG_borderzone"], n_max=5000, seed=3)
    g = subsample(g, "class_b_arm", ["Garcia_MCA_temporal"], n_max=5000, seed=4)
    # Re-normalize lightly within concat after subset — already log-normed per dataset
    both = ad.concat([y, g], join="inner", merge="same", label=None)
    both.obs["class_b_arm"] = both.obs["class_b_arm"].astype("category")
    log(f"Class B DE: {both.n_obs} cells, {both.n_vars} shared genes")
    sc.tl.rank_genes_groups(
        both,
        groupby="class_b_arm",
        groups=["Garcia_MCA_temporal"],
        reference="Yang_SFG_borderzone",
        method="wilcoxon",
        n_genes=150,
    )
    df = sc.get.rank_genes_groups_df(both, group="Garcia_MCA_temporal")
    df["contrast"] = "Garcia_MCA_temporal_vs_Yang_SFG_borderzone"
    df["note"] = "exploratory_cross_dataset_no_batch_correction"
    return df


def main():
    for d in [
        ROOT / "results/yang",
        ROOT / "results/h4",
        ROOT / "results/class_b",
    ]:
        d.mkdir(parents=True, exist_ok=True)

    log("Loading Yang QC object")
    yang = sc.read_h5ad(YANG)
    log("Loading Garcia QC object")
    garcia = sc.read_h5ad(GARCIA)

    # 1) Aim3 lineage DE
    aim3 = run_aim3_lineage(yang)
    aim3.to_csv(OUT_AIM3, index=False)
    log(f"Wrote {OUT_AIM3}")

    # 2) H4 classifier
    h4_df, h4_sum = run_h4_classifier(yang)
    h4_df.to_csv(OUT_H4, index=False)
    h4_md = [
        "# H4 first-pass — Yang endothelial territory classifier",
        "",
        "Trained on **control** endo: PCA hippocampus vs ACA/MCA borderzone SFG.",
        "Transferred to **AD** endo without retraining.",
        "",
        "| Metric | Value |",
        "|---|---:|",
    ]
    for k, v in h4_sum.items():
        if isinstance(v, float):
            h4_md.append(f"| {k} | {v:.4f} |")
        else:
            h4_md.append(f"| {k} | {v} |")
    h4_md += [
        "",
        f"**Delta AUROC (AD - control CV) = {h4_sum['delta_auroc_ad_minus_control']:.4f}**",
        "",
        "Negative delta is consistent with H4 (degradation); treat as exploratory until "
        "donor-aware CV and lineage-cleaned markers are locked.",
    ]
    OUT_H4_MD.write_text("\n".join(h4_md) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_H4}")

    # 3) Class B exploratory
    cb = run_class_b(yang, garcia)
    cb.to_csv(OUT_CB, index=False)
    top = cb.head(20)
    bot = cb.tail(20)
    cb_md = [
        "# Class B exploratory — Garcia MCA temporal endo vs Yang SFG borderzone endo",
        "",
        "**Caveat:** cross-dataset; no batch correction / mutual nearest neighbors yet. "
        "Use for hypothesis generation only.",
        "",
        "## Top genes higher in Garcia MCA",
        "",
        top[["names", "scores", "logfoldchanges", "pvals_adj"]].to_string(index=False),
        "",
        "## Top genes higher in Yang SFG borderzone (bottom of ranked list)",
        "",
        bot[["names", "scores", "logfoldchanges", "pvals_adj"]].to_string(index=False),
    ]
    OUT_CB_MD.write_text("\n".join(cb_md) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_CB}")

    # Summary
    lines = [
        "# Analysis progress — existing data (ATLAS/RADC/ADNI pending)",
        "",
        "## Yang Aim3 (lineage-stratified)",
        f"- DE table: `{OUT_AIM3.relative_to(ROOT).as_posix()}`",
        f"- Rows: {len(aim3):,} (endo/peri × control/AD)",
        "",
        "## H4 first-pass",
        f"- Control CV AUROC: **{h4_sum['control_cv_auroc']:.3f}**",
        f"- AD transfer AUROC: **{h4_sum['ad_transfer_auroc']:.3f}**",
        f"- Delta (AD-control): **{h4_sum['delta_auroc_ad_minus_control']:.3f}**",
        f"- Report: `{OUT_H4_MD.relative_to(ROOT).as_posix()}`",
        "",
        "## Class B exploratory",
        f"- `{OUT_CB.relative_to(ROOT).as_posix()}`",
        "",
        "## Still waiting",
        "- ATLAS encryption code email",
        "- RADC clinical-pathologic approval",
        "- ADNI approval",
    ]
    OUT_SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")
    try:
        print("\n".join(lines))
    except UnicodeEncodeError:
        print(OUT_SUMMARY.read_text(encoding="utf-8").encode("ascii", "replace").decode("ascii"))


if __name__ == "__main__":
    main()
