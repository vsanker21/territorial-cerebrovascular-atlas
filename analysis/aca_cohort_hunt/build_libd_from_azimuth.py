"""Convert Azimuth-annotated LIBD dACC vascular MTX → h5ad (method-matched to Shah et al.).

Prefer Zenodo/GLOBUS sce_azimuth.Rdata when available (paper-identical cell set).
This path re-annotates GEO CellRanger matrices with humancortexref (same reference).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[2]
MTXDIR = ROOT / "data/processed/aca_cohort_hunt/libd_azimuth_vasc_mtx"
OUT = ROOT / "data/processed/aca_cohort_hunt/libd_dacc_azimuth_vascular.h5ad"
OUT_SUM = ROOT / "results/aca_cohort_hunt/libd_dacc_azimuth_summary.json"

# Published Shah et al. humancortexref counts (from their azimuth.txt log)
PUBLISHED_EXPECT = {"Endo": 89, "VLMC": 159}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def main():
    meta_path = MTXDIR / "metadata.csv"
    mtx_path = MTXDIR / "matrix.mtx"
    feat_path = MTXDIR / "features.tsv"
    bc_path = MTXDIR / "barcodes.tsv"
    for p in (meta_path, mtx_path, feat_path, bc_path):
        if not p.exists():
            raise FileNotFoundError(f"Missing Azimuth export: {p}")

    meta = pd.read_csv(meta_path)
    genes = pd.read_csv(feat_path, header=None)[0].astype(str).tolist()
    barcodes = pd.read_csv(bc_path, header=None)[0].astype(str).tolist()
    mat = mmread(mtx_path)
    if sparse.issparse(mat):
        mat = mat.tocsr()
    else:
        mat = sparse.csr_matrix(mat)
    if mat.shape != (len(genes), len(barcodes)):
        if mat.shape == (len(barcodes), len(genes)):
            mat = mat.T.tocsr()
        else:
            raise ValueError(f"shape mismatch {mat.shape} vs genes={len(genes)} bc={len(barcodes)}")

    obs = meta.copy()
    if "cell_barcode" in obs.columns:
        obs.index = obs["cell_barcode"].astype(str)
    else:
        obs.index = barcodes
    # align
    obs = obs.reindex(barcodes)
    obs["donor_id"] = obs.get("donor_id", pd.Series(index=obs.index, dtype=object))
    ct_col = None
    for c in ("predicted.subclass", "predicted.celltype.l2", "predicted.celltype.l1"):
        if c in obs.columns:
            ct_col = c
            break
    if ct_col is None:
        pred = [c for c in obs.columns if str(c).startswith("predicted")]
        ct_col = pred[0] if pred else None
    if ct_col is None:
        raise RuntimeError("No predicted.* annotation column in metadata")

    a = ad.AnnData(X=mat.T.tocsr(), obs=obs, var=pd.DataFrame(index=genes))
    a.obs["cell_type"] = a.obs[ct_col].astype(str)
    a.obs["annotation_source"] = "Azimuth_humancortexref_method_matched"
    a.obs["source_dataset"] = "LIBD_spatialdACC_GSE296731"
    a.obs["cohort_arm"] = "ACA_cortex"
    a.obs["region"] = "dACC"
    a.obs["tv_dominant"] = "ACA"
    a.obs["tv_p_ACA"] = 1.0
    a.obs["tv_assignment_fraction"] = 0.72
    a.obs["tv_source"] = "TerritoryVoxels_HO_cingulate_anterior"
    a.obs["vessel_enriched"] = False

    # keep Endo + VLMC primarily; optional peri/SMC if present
    labs = a.obs["cell_type"].str.lower()
    keep = labs.str.contains("endo|vlmc")
    a = a[keep].copy()

    counts = a.obs["cell_type"].value_counts().to_dict()
    summary = {
        "object": str(OUT.relative_to(ROOT)),
        "n_vascular": int(a.n_obs),
        "n_donors": int(a.obs["donor_id"].nunique()) if "donor_id" in a.obs else None,
        "cell_type_counts": counts,
        "annotation_column": ct_col,
        "annotation_source": "Azimuth_humancortexref_method_matched",
        "published_expect_from_paper_log": PUBLISHED_EXPECT,
        "note": "Paper-identical labels require Zenodo/GLOBUS sce_azimuth.Rdata; this is same reference re-run on GEO matrices.",
        "vessel_enriched": False,
        "verdict": "GO_method_matched_Azimuth_unsorted_dACC",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT_SUM.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(OUT, compression="gzip")
    OUT_SUM.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
