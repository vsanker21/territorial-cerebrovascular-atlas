"""
Yang vascular QC + Aim3 (PCA hippocampus vs ACA/MCA borderzone SFG) marker DE.

Uses TerritoryVoxels streams already on yang_vascular.h5ad.
Writes a QC-filtered object and a first-pass within-Yang Aim3 DE table
(controls and AD reported separately; discovery emphasis on controls).
"""
from __future__ import annotations

import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]
IN_H5AD = ROOT / "data/processed/yang/yang_vascular.h5ad"
OUT_H5AD = ROOT / "data/processed/yang/yang_vascular_qc.h5ad"
OUT_DE = ROOT / "results/yang/aim3_borderzone_vs_pca_de.csv"
OUT_QC = ROOT / "results/yang/yang_vascular_qc_summary.csv"
OUT_MD = ROOT / "results/yang/yang_vascular_qc_aim3_report.md"

# Canonical vascular / contaminant markers (Yang Nature 2022 / Garcia-style panels)
MARKERS = {
    "endo": ["PECAM1", "CLDN5", "FLT1", "VWF", "CDH5", "KDR"],
    "peri": ["PDGFRB", "RGS5", "KCNJ8", "ABCC9", "CSPG4"],
    "smc": ["ACTA2", "TAGLN", "MYH11", "CNN1"],
    "fib": ["COL1A1", "DCN", "LUM", "PDGFRA"],
    "neuron": ["SNAP25", "SLC17A7", "GAD1", "RBFOX3"],
    "astro": ["AQP4", "GFAP", "ALDH1L1"],
    "oligo": ["MBP", "MOBP", "PLP1"],
    "micro": ["PTPRC", "CX3CR1", "P2RY12"],
}

MIN_GENES = 200
MAX_MT = 10.0
MIN_COUNTS = 300
# Soft vascular gate: mean log1p of endo+peri+smc+fib > neuron/astro/oligo
VAS_SCORE_MARGIN = 0.15


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def present(adata: ad.AnnData, genes: list[str]) -> list[str]:
    return [g for g in genes if g in adata.var_names]


def score_modules(adata: ad.AnnData) -> ad.AnnData:
    for name, genes in MARKERS.items():
        g = present(adata, genes)
        if len(g) < 2:
            adata.obs[f"score_{name}"] = 0.0
            continue
        sc.tl.score_genes(adata, g, score_name=f"score_{name}", use_raw=False)
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


def main():
    log(f"Loading {IN_H5AD}")
    adata = sc.read_h5ad(IN_H5AD)
    n0 = adata.n_obs

    # Basic QC
    keep = (
        (adata.obs["n_genes_by_counts"] >= MIN_GENES)
        & (adata.obs["total_counts"] >= MIN_COUNTS)
        & (adata.obs["pct_counts_mt"] <= MAX_MT)
    )
    adata = adata[keep].copy()
    log(f"After basic QC: {adata.n_obs}/{n0}")

    # Normalize for scoring / DE (counts preserved in layers)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    adata = score_modules(adata)
    vasc = adata.obs["score_vascular"] >= (adata.obs["score_parenchyma"] + VAS_SCORE_MARGIN)
    # Also require at least one vascular lineage score above 0
    vasc &= (
        (adata.obs["score_endo"] > 0)
        | (adata.obs["score_peri"] > 0)
        | (adata.obs["score_smc"] > 0)
        | (adata.obs["score_fib"] > 0)
    )
    adata.obs["pass_vascular_gate"] = vasc.astype(bool)
    adata_v = adata[vasc].copy()
    log(f"After vascular gate: {adata_v.n_obs}/{adata.n_obs}")

    # Assign coarse lineage by max score among vascular modules
    lin_scores = adata_v.obs[["score_endo", "score_peri", "score_smc", "score_fib"]]
    adata_v.obs["vascular_lineage"] = (
        lin_scores.idxmax(axis=1).str.replace("score_", "", regex=False).astype("category")
    )

    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    adata_v.write_h5ad(OUT_H5AD, compression="gzip")
    log(f"Wrote {OUT_H5AD}")

    # QC summary tables
    rows = []
    for keys, name in [
        (["analysis_stream", "disease"], "stream_disease"),
        (["territory_contrast", "disease"], "contrast_disease"),
        (["vascular_lineage", "territory_contrast"], "lineage_contrast"),
    ]:
        g = (
            adata_v.obs.groupby(keys, observed=True)
            .size()
            .reset_index(name="n_cells")
        )
        g["table"] = name
        rows.append(g)
    qc = pd.concat(rows, ignore_index=True)
    OUT_QC.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(OUT_QC, index=False)

    # Aim3 DE: PCA_hippocampus vs ACA_MCA_borderzone_SFG
    # Run separately in control and AD; focus report on control discovery
    de_frames = []
    for disease in ["control", "AD"]:
        sub = adata_v[adata_v.obs["disease"].astype(str) == disease].copy()
        # Need both contrasts
        if sub.obs["territory_contrast"].nunique() < 2:
            log(f"Skip DE for {disease}: missing contrast arm")
            continue
        # Cap cells per arm for speed/memory
        parts = []
        for arm in ["PCA_hippocampus", "ACA_MCA_borderzone_SFG"]:
            idx = sub.obs_names[sub.obs["territory_contrast"].astype(str) == arm]
            if len(idx) > 8000:
                rng = np.random.default_rng(42)
                idx = rng.choice(idx, size=8000, replace=False)
            parts.append(sub[idx])
        sub = ad.concat(parts, join="outer", merge="same")
        sub.obs["territory_contrast"] = sub.obs["territory_contrast"].astype("category")
        log(f"DE {disease}: {sub.n_obs} cells")
        sc.tl.rank_genes_groups(
            sub,
            groupby="territory_contrast",
            groups=["ACA_MCA_borderzone_SFG"],
            reference="PCA_hippocampus",
            method="wilcoxon",
            n_genes=200,
            use_raw=False,
        )
        df = sc.get.rank_genes_groups_df(sub, group="ACA_MCA_borderzone_SFG")
        df["disease"] = disease
        df["contrast"] = "ACA_MCA_borderzone_SFG_vs_PCA_hippocampus"
        de_frames.append(df)

    if de_frames:
        de = pd.concat(de_frames, ignore_index=True)
        de.to_csv(OUT_DE, index=False)
        log(f"Wrote {OUT_DE}")
    else:
        de = pd.DataFrame()

    # Report
    stream_counts = (
        adata_v.obs.groupby(["analysis_stream", "disease"], observed=True)
        .size()
        .reset_index(name="n")
    )
    lines = [
        "# Yang vascular QC + Aim3 first-pass DE",
        "",
        f"- Input cells: {n0:,}",
        f"- After basic QC (genes≥{MIN_GENES}, counts≥{MIN_COUNTS}, MT≤{MAX_MT}%): {int(keep.sum()):,}",
        f"- After vascular gate: **{adata_v.n_obs:,}**",
        f"- Object: `{OUT_H5AD.relative_to(ROOT).as_posix()}`",
        "",
        "## Cells by stream × disease",
        "",
        stream_counts.to_string(index=False),
        "",
        "## Lineage × territory_contrast",
        "",
        adata_v.obs.groupby(["vascular_lineage", "territory_contrast"], observed=True)
        .size()
        .reset_index(name="n")
        .to_string(index=False),
        "",
        "## Aim3 DE (borderzone SFG vs PCA hippocampus)",
        "",
        "Positive scores = higher in ACA/MCA borderzone SFG relative to PCA hippocampus.",
        "",
    ]
    if len(de):
        for disease in de["disease"].unique():
            top = de[de["disease"] == disease].head(15)
            lines += [f"### {disease} (top 15)", "", top[["names", "scores", "logfoldchanges", "pvals_adj"]].to_string(index=False), ""]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_MD}")
    print(stream_counts.to_string(index=False))


if __name__ == "__main__":
    main()
