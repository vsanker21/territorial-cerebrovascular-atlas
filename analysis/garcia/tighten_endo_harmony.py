"""
Tighten Garcia endo gates (Yang gate-v2 pattern) on:
  - garcia_temporal_mca_qc.h5ad (ex vivo MCA)
  - garcia_rosmap_vascular.h5ad (ROSMAP MCA isocortex)

Then re-run Harmony Class B DE: Yang tight SFG endo vs Garcia tight MCA endo.

Also expands contaminant caps with RBFOX1 / KCNIP4 (dominant prior Class B DE hits).
"""
from __future__ import annotations

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
YANG = ROOT / "data/processed/yang/yang_vascular_qc_tight.h5ad"
TEMP_IN = ROOT / "data/processed/garcia/garcia_temporal_mca_qc.h5ad"
ROSMAP_IN = ROOT / "data/processed/garcia/garcia_rosmap_vascular.h5ad"
TEMP_OUT = ROOT / "data/processed/garcia/garcia_temporal_mca_qc_tight.h5ad"
ROSMAP_OUT = ROOT / "data/processed/garcia/garcia_rosmap_vascular_tight.h5ad"
TEMP_OUT_STRICT = ROOT / "data/processed/garcia/garcia_temporal_mca_qc_tight_yangstrict.h5ad"

OUT_DIR = ROOT / "results/class_b"
OUT_GATE = OUT_DIR / "garcia_gate_v2_report.md"
OUT_DE = OUT_DIR / "yang_garcia_harmony_class_b_de_tight.csv"
OUT_MD = OUT_DIR / "yang_garcia_harmony_class_b_tight_report.md"
OUT_H5AD = ROOT / "data/processed/class_b/yang_garcia_harmony_endo_tight.h5ad"
OUT_COMPARE = OUT_DIR / "harmony_soft_vs_tight.md"
OUT_SUMMARY = ROOT / "results/analysis_progress_summary.md"
SOFT_DE = OUT_DIR / "yang_garcia_harmony_class_b_de.csv"
SOFT_MD = OUT_DIR / "yang_garcia_harmony_class_b_report.md"

ENDO_MARKERS = ["CLDN5", "PECAM1", "FLT1", "VWF", "CDH5"]
PERI_MARKERS = ["PDGFRB", "RGS5", "KCNJ8", "ABCC9"]
SMC_MARKERS = ["ACTA2", "TAGLN", "MYH11", "CNN1"]
FIB_MARKERS = ["COL1A1", "DCN", "LUM", "PDGFRA"]
# Classic neuron panel for module scores (LINGO1/RBFOX1/KCNIP4 handled as hard caps)
NEURON = ["SNAP25", "SLC17A7", "GAD1", "RBFOX3"]
ASTRO = ["AQP4", "GFAP", "ALDH1L1"]
OLIGO = ["MBP", "MOBP", "PLP1"]
MICRO = ["PTPRC", "CX3CR1", "P2RY12"]

# Hard caps = genes that dominated soft Harmony DE (neuronal ambient).
# Keep set minimal so clean endo remain usable (LINGO1+RBFOX1+KCNIP4 alone ~444 cells).
CONTAMINANT_MAX = {
    "LINGO1": 0.5,
    "RBFOX1": 0.5,
    "KCNIP4": 0.5,
    "SNAP25": 0.5,
}

# Garcia-calibrated gate v2b
# Neuronal hard caps are the Class B artifact fix; oligo ambient is common in Garcia
# snRNA and should not zero the cohort (median oligo score among neuron-clean endo ~0.48).
ENDO_MIN_SCORE = 0.25
LINEAGE_MARGIN = 0.20
PARENCHYMA_MAX = 0.50
NEURON_MAX = 0.20
ASTRO_MAX = 0.50
OLIGO_MAX = 2.00
MICRO_MAX = 0.50
GATE_LABEL = "v2b_garcia_calibrated"


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def present(adata: ad.AnnData, genes: list[str]) -> list[str]:
    return [g for g in genes if g in adata.var_names]


def expr_df(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    g = present(adata, genes)
    if not g:
        return pd.DataFrame(index=adata.obs_names)
    X = adata[:, g].X
    if sparse.issparse(X):
        X = X.toarray()
    return pd.DataFrame(X, index=adata.obs_names, columns=g)


def ensure_module_scores(adata: ad.AnnData) -> ad.AnnData:
    panels = {
        "endo": ENDO_MARKERS,
        "peri": PERI_MARKERS,
        "smc": SMC_MARKERS,
        "fib": FIB_MARKERS,
        "neuron": NEURON,
        "astro": ASTRO,
        "oligo": OLIGO,
        "micro": MICRO,
    }
    for name, genes in panels.items():
        col = f"score_{name}"
        g = present(adata, genes)
        if len(g) < 2:
            adata.obs[col] = 0.0
            continue
        # Recompute always for consistency (neuron panel expanded)
        sc.tl.score_genes(adata, g, score_name=col, use_raw=False)
    adata.obs["score_vascular"] = (
        adata.obs["score_endo"]
        + adata.obs["score_peri"]
        + adata.obs["score_smc"]
        + adata.obs["score_fib"]
    ) / 4.0
    adata.obs["score_parenchyma"] = (
        adata.obs["score_neuron"] + adata.obs["score_astro"] + adata.obs["score_oligo"]
    ) / 3.0
    return adata


def apply_tight_endo(adata: ad.AnnData) -> tuple[ad.AnnData, dict]:
    """Return endo-only tight object + funnel stats."""
    a = ensure_module_scores(adata)
    n0 = a.n_obs
    soft_endo = a[a.obs["vascular_lineage"].astype(str) == "endo"].copy()
    n_soft = soft_endo.n_obs

    expr = expr_df(a, list(CONTAMINANT_MAX) + ENDO_MARKERS)
    contam = pd.Series(False, index=a.obs_names)
    for g, thr in CONTAMINANT_MAX.items():
        if g in expr.columns:
            contam |= expr[g] > thr

    endo_pos = pd.Series(False, index=a.obs_names)
    for g in ENDO_MARKERS:
        if g in expr.columns:
            endo_pos |= expr[g] > 0

    endo_margin = a.obs["score_endo"] - a.obs[["score_peri", "score_smc", "score_fib"]].max(
        axis=1
    )
    clean = (
        (~contam)
        & endo_pos
        & (a.obs["score_endo"] >= ENDO_MIN_SCORE)
        & (endo_margin >= LINEAGE_MARGIN)
        & (a.obs["score_parenchyma"] <= PARENCHYMA_MAX)
        & (a.obs["score_neuron"] <= NEURON_MAX)
        & (a.obs["score_astro"] <= ASTRO_MAX)
        & (a.obs["score_oligo"] <= OLIGO_MAX)
        & (a.obs["score_micro"] <= MICRO_MAX)
    )
    a.obs["pass_endo_gate_v2"] = clean.astype(bool)
    a.obs["contam_gate_v2"] = contam.astype(bool)
    a.obs["vascular_lineage_soft"] = a.obs["vascular_lineage"].astype(str)
    a.obs["vascular_lineage"] = np.where(clean, "endo", "non_endo_or_ambiguous")
    a.obs["gate_version"] = GATE_LABEL

    tight = a[clean].copy()
    tight.obs["vascular_lineage"] = "endo"
    stats = {
        "n_input": int(n0),
        "n_soft_endo": int(n_soft),
        "n_tight_endo": int(tight.n_obs),
        "retention_of_soft_endo": float(tight.n_obs / n_soft) if n_soft else 0.0,
        "pct_soft_endo_LINGO1_gt_0.5": float(
            (expr_df(soft_endo, ["LINGO1"])["LINGO1"] > 0.5).mean()
        )
        if "LINGO1" in soft_endo.var_names
        else np.nan,
        "pct_tight_endo_LINGO1_gt_0.5": float(
            (expr_df(tight, ["LINGO1"])["LINGO1"] > 0.5).mean()
        )
        if tight.n_obs and "LINGO1" in tight.var_names
        else 0.0,
    }
    return tight, stats


def subsample(adata: ad.AnnData, n_max: int, seed: int) -> ad.AnnData:
    if adata.n_obs <= n_max:
        return adata.copy()
    rng = np.random.default_rng(seed)
    idx = rng.choice(adata.obs_names.to_numpy(), size=n_max, replace=False)
    return adata[idx].copy()


def gene_rank(de: pd.DataFrame, gene: str) -> str:
    if de is None or de.empty or "names" not in de.columns:
        return "n/a"
    hit = de.reset_index(drop=True)
    idxs = hit.index[hit["names"] == gene].tolist()
    if not idxs:
        return "not in top list"
    i = idxs[0]
    row = hit.iloc[i]
    return f"rank {i+1}, logFC={row['logfoldchanges']:.2f}, padj={row['pvals_adj']:.2e}"


def run_harmony(
    yang: ad.AnnData,
    garcia_parts: list[tuple[ad.AnnData, str, str]],
    yang_n_max: int = 4000,
):
    y = yang[
        (yang.obs["vascular_lineage"].astype(str) == "endo")
        & (yang.obs["disease"].astype(str) == "control")
        & (yang.obs["territory_contrast"].astype(str) == "ACA_MCA_borderzone_SFG")
    ].copy()
    y.obs["class_b_group"] = "Yang_SFG_borderzone"
    y.obs["batch"] = "Yang"
    y = subsample(y, yang_n_max, 10)

    parts = [y]
    for g, batch, seed_name in garcia_parts:
        if g.n_obs == 0:
            continue
        g = g.copy()
        g.obs["class_b_group"] = "Garcia_MCA"
        g.obs["batch"] = batch
        nmax = 4000 if batch == "Garcia_exvivo" else 2000
        parts.append(subsample(g, nmax, abs(hash(seed_name)) % 10_000))

    shared = parts[0].var_names
    for p in parts[1:]:
        shared = shared.intersection(p.var_names)
    log(f"Shared genes: {len(shared)}; parts={[p.n_obs for p in parts]}")
    parts = [p[:, shared].copy() for p in parts]
    comb = ad.concat(parts, join="inner", merge="same", index_unique="-")

    sc.pp.highly_variable_genes(comb, n_top_genes=2000, flavor="seurat", subset=False)
    hvg = comb[:, comb.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=30, svd_solver="arpack")
    log("Harmony")
    sc.external.pp.harmony_integrate(
        hvg, key="batch", basis="X_pca", adjusted_basis="X_pca_harmony"
    )
    sc.pp.neighbors(hvg, use_rep="X_pca_harmony", n_neighbors=15)
    sc.tl.umap(hvg)

    comb.obs["class_b_group"] = comb.obs["class_b_group"].astype("category")
    comb.obsm["X_pca_harmony"] = hvg.obsm["X_pca_harmony"]
    comb.obsm["X_umap"] = hvg.obsm["X_umap"]
    sc.tl.rank_genes_groups(
        comb,
        groupby="class_b_group",
        groups=["Garcia_MCA"],
        reference="Yang_SFG_borderzone",
        method="wilcoxon",
        n_genes=200,
    )
    de = sc.get.rank_genes_groups_df(comb, group="Garcia_MCA")
    de["contrast"] = "Garcia_MCA_vs_Yang_SFG_borderzone_harmony_tight"
    de["caveat"] = "Harmony batch correction; DE on expression; both arms endo gate v2"

    X = comb.obsm["X_pca_harmony"]
    ybin = (comb.obs["class_b_group"].astype(str) == "Garcia_MCA").astype(int).to_numpy()
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    aucs = []
    for tr, te in skf.split(X, ybin):
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], ybin[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(ybin[te], proba))
    summary = {
        "harmony_cv_auroc_mean": float(np.mean(aucs)),
        "harmony_cv_auroc_std": float(np.std(aucs)),
        "n_yang": int((comb.obs["class_b_group"] == "Yang_SFG_borderzone").sum()),
        "n_garcia_mca": int((comb.obs["class_b_group"] == "Garcia_MCA").sum()),
        "n_batches": int(comb.obs["batch"].nunique()),
        "n_cells": int(comb.n_obs),
    }
    return comb, de, summary


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)

    log("Tighten Garcia temporal")
    temp = sc.read_h5ad(TEMP_IN)
    temp_tight, temp_stats = apply_tight_endo(temp)
    temp_tight.write_h5ad(TEMP_OUT, compression="gzip")
    log(f"Temporal soft endo {temp_stats['n_soft_endo']} -> tight {temp_stats['n_tight_endo']}")

    log("Tighten Garcia ROSMAP")
    rosmap = sc.read_h5ad(ROSMAP_IN)
    # Restrict to MCA isocortex before gating
    rosmap_mca = rosmap[rosmap.obs["class_b_arm"].astype(str) == "MCA_isocortex"].copy()
    # Ensure lineage column exists
    if "vascular_lineage" not in rosmap_mca.obs:
        rosmap_mca.obs["vascular_lineage"] = (
            rosmap_mca.obs["celltype"]
            .astype(str)
            .map({"Endo": "endo", "Per": "peri", "SMC": "smc", "Fib": "fib"})
            .fillna("other")
        )
    # Score + gate on full ROSMAP then subset MCA — better score reference
    rosmap_all = rosmap.copy()
    if rosmap_all.obs["vascular_lineage"].dtype.name != "object":
        rosmap_all.obs["vascular_lineage"] = rosmap_all.obs["vascular_lineage"].astype(str)
    rosmap_tight_all, rosmap_stats_all = apply_tight_endo(rosmap_all)
    # Keep MCA isocortex tight endo
    mca_mask = rosmap_tight_all.obs_names.isin(
        rosmap.obs_names[rosmap.obs["class_b_arm"].astype(str) == "MCA_isocortex"]
    )
    rosmap_tight = rosmap_tight_all[mca_mask].copy()
    # Recompute stats for MCA subset
    soft_mca_endo = rosmap_mca[rosmap_mca.obs["vascular_lineage"].astype(str) == "endo"].n_obs
    rosmap_stats = {
        "n_input": int(rosmap_mca.n_obs),
        "n_soft_endo": int(soft_mca_endo),
        "n_tight_endo": int(rosmap_tight.n_obs),
        "retention_of_soft_endo": float(rosmap_tight.n_obs / soft_mca_endo)
        if soft_mca_endo
        else 0.0,
    }
    rosmap_tight.write_h5ad(ROSMAP_OUT, compression="gzip")
    log(
        f"ROSMAP MCA soft endo {rosmap_stats['n_soft_endo']} -> tight {rosmap_stats['n_tight_endo']}"
    )

    gate_md = f"""# Garcia endo gate v2b (Garcia-calibrated)

## Why not Yang-strict v2?
Yang gate v2 (score_endo>=0.5, margin>=0.4, parenchyma<=0, LINGO1/RBFOX1/KCNIP4 caps)
retained only **~0.5%** of Garcia temporal soft endo (~60 cells) — unusable for DE.
Garcia soft "endo" is heavily neuron-contaminated (~72% LINGO1+).

## Gate v2b rules (primary)
- Canonical endo marker > 0: {ENDO_MARKERS}
- score_endo >= {ENDO_MIN_SCORE}; margin vs other vascular >= {LINEAGE_MARGIN}
- Soft parenchyma caps: parenchyma<={PARENCHYMA_MAX}, neuron<={NEURON_MAX}
- Hard contaminant expression caps: {CONTAMINANT_MAX}
- Neuron module score uses classic panel only ({NEURON}); LINGO1/RBFOX1/KCNIP4 are hard caps

## Temporal ex vivo MCA
- Soft endo: **{temp_stats['n_soft_endo']:,}** -> tight v2b: **{temp_stats['n_tight_endo']:,}** ({100*temp_stats['retention_of_soft_endo']:.1f}%)
- Soft LINGO1>0.5: {100*temp_stats['pct_soft_endo_LINGO1_gt_0.5']:.1f}% -> tight: {100*temp_stats['pct_tight_endo_LINGO1_gt_0.5']:.1f}%
- Object: `{TEMP_OUT.relative_to(ROOT).as_posix()}`

## ROSMAP MCA isocortex (Midtemporal+Angular)
- Soft endo: **{rosmap_stats['n_soft_endo']:,}** -> tight v2b: **{rosmap_stats['n_tight_endo']:,}** ({100*rosmap_stats['retention_of_soft_endo']:.1f}%)
- Object: `{ROSMAP_OUT.relative_to(ROOT).as_posix()}`
"""
    OUT_GATE.write_text(gate_md, encoding="utf-8")
    log(f"Wrote {OUT_GATE}")

    if temp_tight.n_obs + rosmap_tight.n_obs < 100:
        raise SystemExit(
            f"Too few tight Garcia endo cells "
            f"(temp={temp_tight.n_obs}, rosmap={rosmap_tight.n_obs})"
        )

    log("Harmony Class B on tight arms")
    yang = sc.read_h5ad(YANG)
    # Balance Yang closer to Garcia n for fairer DE
    n_garcia = temp_tight.n_obs + rosmap_tight.n_obs
    yang_n = int(min(4000, max(1000, 3 * n_garcia)))
    comb, de, summ = run_harmony(
        yang,
        [
            (temp_tight, "Garcia_exvivo", "exvivo"),
            (rosmap_tight, "Garcia_ROSMAP", "rosmap"),
        ],
        yang_n_max=yang_n,
    )
    comb.write_h5ad(OUT_H5AD, compression="gzip")
    de.to_csv(OUT_DE, index=False)

    top = de.head(20)
    bot = de.tail(15)
    md = f"""# Harmonized Class B — tight endo gate v2b (Yang + Garcia)

Batches: Yang tight SFG | Garcia ex vivo MCA v2b | Garcia ROSMAP MCA v2b.
Hard neuronal contaminant caps (LINGO1/RBFOX1/KCNIP4); Garcia-calibrated score thresholds.

## Separability on Harmony embedding

| Metric | Value |
|---|---:|
| 5-fold CV AUROC | {summ['harmony_cv_auroc_mean']:.3f} +/- {summ['harmony_cv_auroc_std']:.3f} |
| N Yang SFG | {summ['n_yang']} |
| N Garcia MCA | {summ['n_garcia_mca']} |
| N cells / batches | {summ['n_cells']} / {summ['n_batches']} |

## Top genes higher in Garcia MCA

{top[['names','scores','logfoldchanges','pvals_adj']].to_string(index=False)}

## Higher in Yang SFG (bottom ranks)

{bot[['names','scores','logfoldchanges','pvals_adj']].to_string(index=False)}

Object: `{OUT_H5AD.relative_to(ROOT).as_posix()}`
"""
    OUT_MD.write_text(md, encoding="utf-8")
    log(f"Wrote {OUT_MD}")

    # Soft vs tight comparison
    soft_de = pd.read_csv(SOFT_DE) if SOFT_DE.exists() else pd.DataFrame()
    soft_auroc = "n/a"
    if SOFT_MD.exists():
        for line in SOFT_MD.read_text(encoding="utf-8").splitlines():
            if "5-fold CV AUROC" in line:
                soft_auroc = line.split("|")[-2].strip() if "|" in line else line
                break
    neuron_genes = ["LINGO1", "KCNIP4", "RBFOX1", "CSMD1", "SNAP25", "CLDN5", "PECAM1", "FLT1", "VWF"]
    cmp_rows = ["| Gene | Soft DE | Tight DE |", "|---|---|---|"]
    for g in neuron_genes:
        cmp_rows.append(f"| {g} | {gene_rank(soft_de, g)} | {gene_rank(de, g)} |")
    cmp = f"""# Harmony Class B — soft vs tight Garcia endo (gate v2b)

## Retention
- Temporal soft endo {temp_stats['n_soft_endo']:,} -> tight v2b {temp_stats['n_tight_endo']:,} ({100*temp_stats['retention_of_soft_endo']:.1f}%)
- ROSMAP MCA soft endo {rosmap_stats['n_soft_endo']:,} -> tight v2b {rosmap_stats['n_tight_endo']:,} ({100*rosmap_stats['retention_of_soft_endo']:.1f}%)
- Note: Yang-strict v2 retained ~60 temporal cells (0.5%) — too few; v2b is primary

## Harmony AUROC
| Gate | AUROC |
|---|---|
| Soft | {soft_auroc} |
| Tight v2b | {summ['harmony_cv_auroc_mean']:.3f} +/- {summ['harmony_cv_auroc_std']:.3f} |

## Key gene ranks in Garcia-higher DE list
{chr(10).join(cmp_rows)}

Lower/absent ranks for LINGO1/KCNIP4/RBFOX1 = cleaner Class B.
"""
    OUT_COMPARE.write_text(cmp, encoding="utf-8")
    log(f"Wrote {OUT_COMPARE}")

    summary = f"""# Analysis progress — existing data

## Garcia endo gate v2b + Harmony Class B
- Temporal: {temp_stats['n_soft_endo']:,} soft endo -> **{temp_stats['n_tight_endo']:,}** tight ({100*temp_stats['retention_of_soft_endo']:.1f}%)
- ROSMAP MCA: {rosmap_stats['n_soft_endo']:,} soft endo -> **{rosmap_stats['n_tight_endo']:,}** tight ({100*rosmap_stats['retention_of_soft_endo']:.1f}%)
- Yang-strict v2 was too harsh on Garcia (~0.5% retention); v2b = hard neuronal caps + calibrated scores
- Harmony CV AUROC: **{summ['harmony_cv_auroc_mean']:.3f}** +/- {summ['harmony_cv_auroc_std']:.3f}
- Cells: Yang {summ['n_yang']} vs Garcia MCA {summ['n_garcia_mca']}
- Gate report: `results/class_b/garcia_gate_v2_report.md`
- Soft vs tight: `results/class_b/harmony_soft_vs_tight.md`
- DE: `results/class_b/yang_garcia_harmony_class_b_de_tight.csv`
- Object: `data/processed/class_b/yang_garcia_harmony_endo_tight.h5ad`

## ATLAS
- Still blocked on passphrase — `results/atlas/atlas_passphrase_gate.md`

## Still waiting
- ATLAS encryption code email
- RADC clinical-pathologic approval
- ADNI approval
"""
    OUT_SUMMARY.write_text(summary, encoding="utf-8")
    try:
        print(cmp)
        print(summary)
    except UnicodeEncodeError:
        print(cmp.encode("ascii", "replace").decode("ascii"))


if __name__ == "__main__":
    main()
