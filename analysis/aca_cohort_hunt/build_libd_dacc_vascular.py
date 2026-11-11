"""
Build LIBD spatialdACC (GSE296731) dACC vascular seed from CellRanger matrices.

True dorsal anterior cingulate (BA24/32/33), n=10 neurotypical adults.
Not vessel-enriched — gate by endothelial/mural markers + neuronal/oligo caps.

Outputs:
  data/processed/aca_cohort_hunt/libd_dacc_vascular_seed.h5ad
  results/aca_cohort_hunt/libd_dacc_vascular_summary.json
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

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[2]
EXTRACTED = ROOT / "data/raw/aca_cohort_hunt/GSE296731/extracted"
OUT_H5AD = ROOT / "data/processed/aca_cohort_hunt/libd_dacc_vascular_seed.h5ad"
OUT_DIR = ROOT / "results/aca_cohort_hunt"
OUT_SUMMARY = OUT_DIR / "libd_dacc_vascular_summary.json"

# TerritoryVoxels prior for anterior cingulate (from Linnarsson A24 probe)
TV_PRIOR = {
    "tv_dominant": "ACA",
    "tv_p_ACA": 1.0,
    "tv_p_MCA": 0.0,
    "tv_p_PCA": 0.0,
    "tv_p_VB": 0.0,
    "tv_p_dominant": 1.0,
    "tv_assignment_fraction": 0.724,
    "tv_source": "TerritoryVoxels_HO_cingulate_anterior",
    "anatomy": "dACC_BA24_32_33",
    "cohort_arm": "ACA_cortex",
}

ENDO = ["CLDN5", "FLT1", "PECAM1", "VWF", "EGFL7", "CDH5", "KDR", "ESAM"]
PERI = ["PDGFRB", "RGS5", "KCNJ8", "ABCC9", "CSPG4"]
VSMC = ["ACTA2", "TAGLN", "MYH11", "CNN1"]
FIB = ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"]
NEUR = ["SNAP25", "SYT1", "RBFOX3", "SLC17A7", "GAD1", "GAD2", "LINGO1", "KCNIP4", "RBFOX1"]
OLIGO = ["MBP", "PLP1", "MOBP", "MOG"]
ASTRO = ["AQP4", "GFAP", "ALDH1L1"]
MICRO = ["CX3CR1", "P2RY12", "CSF1R"]


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def score_genes(adata: ad.AnnData, genes: list[str]) -> np.ndarray:
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return np.zeros(adata.n_obs, dtype=np.float32)
    X = adata[:, present].X
    if sparse.issparse(X):
        X = X.toarray()
    # mean log1p already applied
    return np.asarray(X.mean(axis=1)).ravel().astype(np.float32)


def load_sample(prefix: Path, sample_id: str, donor: str) -> ad.AnnData:
    # GEO names: GSM..._BrXXXX-Nc_dACC_{matrix,barcodes,features}
    mtx = next(prefix.glob(f"*{sample_id}*matrix.mtx.gz"))
    stem = mtx.name.replace("_matrix.mtx.gz", "")
    barcodes = prefix / f"{stem}_barcodes.tsv.gz"
    features = prefix / f"{stem}_features.tsv.gz"
    log(f"Loading {stem}")
    a = sc.read_10x_mtx(prefix, prefix=f"{stem}_", var_names="gene_symbols", cache=False)
    # read_10x_mtx expects directory with files named matrix.mtx etc OR use prefix carefully.
    # Fallback manual:
    raise NotImplementedError


def load_mtx_triplet(mtx: Path) -> ad.AnnData:
    stem = mtx.name.replace("_matrix.mtx.gz", "")
    parent = mtx.parent
    barcodes = parent / f"{stem}_barcodes.tsv.gz"
    features = parent / f"{stem}_features.tsv.gz"
    log(f"Reading {stem}")
    from scipy.io import mmread

    mat = mmread(mtx).T.tocsr().astype(np.float32)  # cells x genes
    bc = pd.read_csv(barcodes, header=None, sep="\t")[0].astype(str).values
    feat = pd.read_csv(features, header=None, sep="\t")
    # 10x features: gene_id, gene_symbol, feature_type
    genes = feat[1].astype(str).values if feat.shape[1] >= 2 else feat[0].astype(str).values
    # make unique
    genes = pd.Index(genes).astype(str)
    genes = genes.where(~genes.duplicated(), genes + "-dup")
    a = ad.AnnData(mat)
    a.obs_names = pd.Index([f"{stem}-{b}" for b in bc])
    a.var_names = genes
    # parse donor from stem e.g. GSM8975421_Br2720-1c_dACC
    parts = stem.split("_")
    gsm = parts[0]
    donor = parts[1].split("-")[0] if len(parts) > 1 else "unknown"
    a.obs["geo_sample"] = gsm
    a.obs["donor_id"] = donor
    a.obs["library_id"] = stem
    a.obs["region"] = "dACC"
    return a


def gate_vascular(a: ad.AnnData) -> ad.AnnData:
    # GEO dumps are unfiltered CellRanger (includes empty droplets).
    a.obs["n_counts"] = np.asarray(a.X.sum(axis=1)).ravel()
    a.obs["n_genes"] = np.asarray((a.X > 0).sum(axis=1)).ravel()
    a = a[(a.obs["n_counts"] >= 800) & (a.obs["n_genes"] >= 400)].copy()
    sc.pp.filter_genes(a, min_cells=3)
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    a.obs["score_endo"] = score_genes(a, ENDO)
    a.obs["score_peri"] = score_genes(a, PERI)
    a.obs["score_vsmc"] = score_genes(a, VSMC)
    a.obs["score_fib"] = score_genes(a, FIB)
    a.obs["score_neur"] = score_genes(a, NEUR)
    a.obs["score_oligo"] = score_genes(a, OLIGO)
    a.obs["score_astro"] = score_genes(a, ASTRO)
    a.obs["score_micro"] = score_genes(a, MICRO)
    a.obs["score_vasc"] = np.maximum.reduce(
        [
            a.obs["score_endo"].to_numpy(),
            a.obs["score_peri"].to_numpy(),
            a.obs["score_vsmc"].to_numpy(),
            a.obs["score_fib"].to_numpy(),
        ]
    )

    # Among QC nuclei: top vascular quantile + hard neuronal/oligo caps
    v = a.obs["score_vasc"].to_numpy()
    thr = float(np.quantile(v, 0.97)) if len(v) else 0.0
    keep = (
        (a.obs["score_vasc"] >= thr)
        & (a.obs["score_neur"] < 0.30)
        & (a.obs["score_oligo"] < 0.70)
        & (a.obs["score_astro"] < 0.70)
        & (a.obs["score_micro"] < 0.70)
    )
    # assign provisional type by max vascular subscore among keepers
    b = a[keep].copy()
    if b.n_obs == 0:
        return b
    subs = b.obs[["score_endo", "score_peri", "score_vsmc", "score_fib"]]
    b.obs["cell_type_seed"] = subs.idxmax(axis=1).map(
        {
            "score_endo": "endothelial_seed",
            "score_peri": "pericyte_seed",
            "score_vsmc": "vSMC_seed",
            "score_fib": "fibroblast_seed",
        }
    )
    return b


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    mtcs = sorted(EXTRACTED.glob("*_matrix.mtx.gz"))
    if not mtcs:
        raise FileNotFoundError(f"No matrices in {EXTRACTED}")

    pieces = []
    per_lib = {}
    for mtx in mtcs:
        raw = load_mtx_triplet(mtx)
        n0 = raw.n_obs
        vasc = gate_vascular(raw)
        per_lib[raw.obs["library_id"].iloc[0]] = {
            "n_raw": int(n0),
            "n_vasc_seed": int(vasc.n_obs),
            "donor_id": str(raw.obs["donor_id"].iloc[0]),
            "type_counts": (
                vasc.obs["cell_type_seed"].value_counts().to_dict() if vasc.n_obs else {}
            ),
        }
        log(f"  -> vasc seed {vasc.n_obs}/{n0}")
        if vasc.n_obs:
            pieces.append(vasc)

    if not pieces:
        raise RuntimeError("No vascular seed cells retained")

    # intersect genes
    genes = set(pieces[0].var_names)
    for p in pieces[1:]:
        genes &= set(p.var_names)
    genes = sorted(genes)
    pieces = [p[:, genes].copy() for p in pieces]
    # restore counts layer preference for downstream: use log-normalized X
    pooled = ad.concat(pieces, join="inner", merge="same", index_unique=None)
    for k, v in TV_PRIOR.items():
        pooled.obs[k] = v
    pooled.obs["source_dataset"] = "LIBD_spatialdACC_GSE296731"
    pooled.write_h5ad(OUT_H5AD)

    summary = {
        "object": str(OUT_H5AD.relative_to(ROOT)),
        "n_vascular_seed": int(pooled.n_obs),
        "n_donors": int(pooled.obs["donor_id"].nunique()),
        "n_libraries": int(pooled.obs["library_id"].nunique()),
        "type_counts": pooled.obs["cell_type_seed"].value_counts().to_dict(),
        "per_library": per_lib,
        "tv_prior": TV_PRIOR,
        "gate": {
            "vasc_quantile": 0.98,
            "neur_cap": 0.35,
            "oligo_cap": 0.8,
            "note": "Unsorted cortex seed; tighten before Class B DE",
        },
        "verdict": (
            "GO_external_ACA_dACC_vascular_seed"
            if pooled.n_obs >= 100
            else "MARGINAL_n_low"
        ),
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
