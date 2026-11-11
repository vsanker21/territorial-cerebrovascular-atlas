"""
Load Yang GSE163577 10x filtered matrices and attach TerritoryVoxels streams.

Outputs
-------
data/processed/yang/yang_vascular.h5ad
data/processed/yang/yang_load_qc_summary.csv
results/yang/yang_load_report.md

QC (library-level, pre-integration; not a full cell-type call):
  - n_genes_by_counts, total_counts, pct_counts_mt
  - Soft filter report only; full object keeps all barcodes with metrics attached
    so downstream vascular subsetting can use paper annotations or marker gates.
"""
from __future__ import annotations

import json
import shutil
import tarfile
import tempfile
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data/raw/yang/GSE163577_RAW"
META = ROOT / "data/processed/yang/yang_sample_metadata.csv"
TV = ROOT / "data/processed/yang/yang_sample_territory_voxels.csv"
OUT_H5AD = ROOT / "data/processed/yang/yang_vascular.h5ad"
OUT_QC = ROOT / "data/processed/yang/yang_load_qc_summary.csv"
OUT_MD = ROOT / "results/yang/yang_load_report.md"
EXTRACT_DIR = ROOT / "data/raw/yang/GSE163577_RAW_extracted"


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def ensure_extracted(tar_path: Path) -> Path:
    """Extract a per-sample 10x tar.gz once; return directory containing matrix.mtx.gz."""
    EXTRACT_DIR.mkdir(parents=True, exist_ok=True)
    # Prefer finding an already-extracted matrix folder matching GSM
    gsm = tar_path.name.split("_")[0]
    existing = list(EXTRACT_DIR.glob(f"**/{gsm}*/matrix.mtx.gz")) + list(
        EXTRACT_DIR.glob(f"**/matrix.mtx.gz")
    )
    # Narrow: extract into EXTRACT_DIR / tar stem
    dest = EXTRACT_DIR / tar_path.name.replace(".tar.gz", "")
    mtx = list(dest.rglob("matrix.mtx.gz"))
    if mtx:
        return mtx[0].parent
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as tf:
        tf.extractall(dest)
    mtx = list(dest.rglob("matrix.mtx.gz"))
    if not mtx:
        raise FileNotFoundError(f"No matrix.mtx.gz in {tar_path}")
    return mtx[0].parent


def load_one(tar_path: Path, sample_row: pd.Series, tv_row: pd.Series) -> ad.AnnData:
    mtx_dir = ensure_extracted(tar_path)
    adata = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    gsm = sample_row["geo_accession"]
    adata.obs_names = [f"{gsm}_{bc}" for bc in adata.obs_names]
    # Sample metadata
    for col in ("geo_accession", "title", "donor_id", "disease", "region", "brain_region"):
        if col in sample_row.index:
            adata.obs[col] = sample_row[col]
    adata.obs["dataset"] = "yang_gse163577"
    adata.obs["library_id"] = gsm
    # TerritoryVoxels stream
    tv_cols = [
        "analysis_stream",
        "territory_method",
        "tv_dominant",
        "tv_p_dominant",
        "tv_assignment_fraction",
        "tv_entropy_bits",
        "tv_p_ACA",
        "tv_p_MCA",
        "tv_p_PCA",
        "tv_p_VB",
        "tv_aca_mca_borderzone_flag",
        "borderzone_defined_fraction",
        "borderzone_mean_where_defined",
        "hard_dominant",
        "hard_p_ACA",
        "hard_p_MCA",
        "hard_p_PCA",
        "hard_p_VB",
    ]
    for col in tv_cols:
        if col in tv_row.index:
            adata.obs[col] = tv_row[col]
    # Contrast labels for H1 / Aim3 design
    stream = str(tv_row.get("analysis_stream", ""))
    if stream == "H1_primary":
        adata.obs["territory_contrast"] = "PCA_hippocampus"
    elif stream == "Aim3_borderzone":
        adata.obs["territory_contrast"] = "ACA_MCA_borderzone_SFG"
    else:
        adata.obs["territory_contrast"] = "other"
    return adata


def qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    return adata


def main():
    if not META.exists() or not TV.exists():
        raise SystemExit("Missing metadata / TerritoryVoxels CSVs; run parse + attach scripts first")
    meta = pd.read_csv(META)
    tv = pd.read_csv(TV)
    merged = meta.merge(tv, on=["geo_accession", "donor_id", "disease", "region", "title"], how="left")
    if merged["analysis_stream"].isna().any():
        # fallback merge on geo only
        tv2 = tv.drop(columns=[c for c in tv.columns if c in meta.columns and c != "geo_accession"])
        merged = meta.merge(tv2, on="geo_accession", how="left")

    tars = {p.name.split("_")[0]: p for p in sorted(RAW_DIR.glob("GSM*_filtered_feature_bc_matrix.tar.gz"))}
    missing = [g for g in merged["geo_accession"] if g not in tars]
    if missing:
        raise SystemExit(f"Missing 10x tars for: {missing}")

    adatas = []
    lib_rows = []
    for _, row in merged.iterrows():
        gsm = row["geo_accession"]
        log(f"Loading {gsm} ({row['region']}, {row.get('analysis_stream')})")
        a = load_one(tars[gsm], row, row)
        a = qc_metrics(a)
        lib_rows.append(
            {
                "geo_accession": gsm,
                "donor_id": row["donor_id"],
                "disease": row["disease"],
                "region": row["region"],
                "analysis_stream": row.get("analysis_stream"),
                "territory_contrast": a.obs["territory_contrast"].iloc[0],
                "tv_dominant": row.get("tv_dominant"),
                "tv_p_ACA": row.get("tv_p_ACA"),
                "tv_p_MCA": row.get("tv_p_MCA"),
                "tv_p_PCA": row.get("tv_p_PCA"),
                "n_cells": a.n_obs,
                "n_genes": a.n_vars,
                "median_genes": float(np.median(a.obs["n_genes_by_counts"])),
                "median_counts": float(np.median(a.obs["total_counts"])),
                "median_pct_mt": float(np.median(a.obs["pct_counts_mt"])),
            }
        )
        adatas.append(a)

    log(f"Concatenating {len(adatas)} libraries…")
    adata = ad.concat(adatas, join="outer", merge="same", label="batch_idx", index_unique=None)
    # Ensure categorical streams
    for col in ("region", "disease", "analysis_stream", "territory_contrast", "tv_dominant", "donor_id"):
        if col in adata.obs:
            adata.obs[col] = adata.obs[col].astype("category")

    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    adata.uns["territory_annotation"] = {
        "method": "TerritoryVoxels_BMM_labels_1to4",
        "source": "Liu Sci Data 2023 / NITRC Atlas_182_MNI152",
        "yang_contrast": "PCA_hippocampus vs ACA_MCA_borderzone_SFG",
        "note": "Superior frontal is Aim3 borderzone, not clean MCA",
    }
    adata.uns["load_timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    log(f"Writing {OUT_H5AD} ({adata.n_obs} cells × {adata.n_vars} genes)")
    adata.write_h5ad(OUT_H5AD, compression="gzip")

    qc = pd.DataFrame(lib_rows)
    qc.to_csv(OUT_QC, index=False)

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    by_stream = (
        qc.groupby("analysis_stream", dropna=False)
        .agg(n_libs=("geo_accession", "count"), n_cells=("n_cells", "sum"))
        .reset_index()
    )
    lines = [
        "# Yang 10x load + TerritoryVoxels attach",
        "",
        f"- Cells: **{adata.n_obs:,}**",
        f"- Genes: **{adata.n_vars:,}**",
        f"- Libraries: **{len(qc)}**",
        f"- Object: `{OUT_H5AD.relative_to(ROOT).as_posix()}`",
        "",
        "## Streams",
        "",
        by_stream.to_string(index=False),
        "",
        "## Per-library QC (medians)",
        "",
        qc[
            [
                "geo_accession",
                "region",
                "disease",
                "analysis_stream",
                "n_cells",
                "median_genes",
                "median_counts",
                "median_pct_mt",
            ]
        ].to_string(index=False),
        "",
        "Contrast lock: `territory_contrast` = PCA_hippocampus vs ACA_MCA_borderzone_SFG.",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    log("Done")
    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()
