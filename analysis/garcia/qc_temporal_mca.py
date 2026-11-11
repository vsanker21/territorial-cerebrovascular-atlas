"""Garcia temporal MCA QC summary (object already TerritoryVoxels-annotated)."""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]
IN_H5AD = ROOT / "data/processed/garcia/garcia_temporal_mca.h5ad"
OUT_H5AD = ROOT / "data/processed/garcia/garcia_temporal_mca_qc.h5ad"
OUT_CSV = ROOT / "results/garcia/garcia_temporal_mca_qc_summary.csv"
OUT_MD = ROOT / "results/garcia/garcia_temporal_mca_qc_report.md"

MIN_GENES = 200
MAX_MT = 15.0
MIN_COUNTS = 300


def main():
    adata = sc.read_h5ad(IN_H5AD)
    n0 = adata.n_obs
    keep = (
        (adata.obs["n_genes_by_counts"] >= MIN_GENES)
        & (adata.obs["total_counts"] >= MIN_COUNTS)
        & (adata.obs["pct_counts_mt"] <= MAX_MT)
    )
    adata = adata[keep].copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    for name, genes in {
        "endo": ["PECAM1", "CLDN5", "FLT1", "VWF"],
        "peri": ["PDGFRB", "RGS5", "KCNJ8"],
        "smc": ["ACTA2", "TAGLN", "MYH11"],
        "fib": ["COL1A1", "DCN", "LUM"],
    }.items():
        g = [x for x in genes if x in adata.var_names]
        if len(g) >= 2:
            sc.tl.score_genes(adata, g, score_name=f"score_{name}")
        else:
            adata.obs[f"score_{name}"] = 0.0
    scores = adata.obs[["score_endo", "score_peri", "score_smc", "score_fib"]]
    adata.obs["vascular_lineage"] = (
        scores.idxmax(axis=1).str.replace("score_", "", regex=False).astype("category")
    )
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    summ = (
        adata.obs.groupby(["SampleID", "tissue_location", "vascular_lineage"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    summ.to_csv(OUT_CSV, index=False)
    lines = [
        "# Garcia temporal MCA QC",
        "",
        f"- Input: {n0:,} → QC: **{adata.n_obs:,}**",
        f"- Territory: MCA (P=1.0), stream ClassB_MCA",
        f"- Object: `{OUT_H5AD.relative_to(ROOT).as_posix()}`",
        "",
        adata.obs["vascular_lineage"].value_counts().to_string(),
        "",
        summ.groupby("vascular_lineage", observed=True)["n_cells"].sum().to_string(),
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(adata.obs["vascular_lineage"].value_counts())


if __name__ == "__main__":
    main()
