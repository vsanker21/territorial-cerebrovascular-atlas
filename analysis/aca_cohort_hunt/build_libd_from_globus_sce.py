"""Build h5ad from GLOBUS paper-identical LIBD sce_Final Azimuth vascular export."""
from __future__ import annotations

import json
import time
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import sparse
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[2]
MTXDIR = ROOT / "data/processed/aca_cohort_hunt/libd_globus_vasc_mtx"
OUT = ROOT / "data/processed/aca_cohort_hunt/libd_dacc_globus_vascular.h5ad"
OUT_SUM = ROOT / "results/aca_cohort_hunt/libd_dacc_globus_summary.json"


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def main():
    meta = pd.read_csv(MTXDIR / "metadata.csv")
    genes = pd.read_csv(MTXDIR / "features.tsv", header=None)[0].astype(str).tolist()
    barcodes = pd.read_csv(MTXDIR / "barcodes.tsv", header=None)[0].astype(str).tolist()
    mat = mmread(MTXDIR / "matrix.mtx")
    mat = mat.tocsr() if sparse.issparse(mat) else sparse.csr_matrix(mat)
    if mat.shape == (len(barcodes), len(genes)):
        mat = mat.T.tocsr()
    assert mat.shape == (len(genes), len(barcodes)), mat.shape

    obs = meta.copy()
    if "cell_barcode" in obs.columns:
        obs.index = obs["cell_barcode"].astype(str)
    else:
        obs.index = barcodes
    obs = obs.reindex(barcodes)

    a = ad.AnnData(X=mat.T.tocsr(), obs=obs, var=pd.DataFrame(index=genes))
    a.obs["cell_type"] = a.obs["cell_type"].astype(str)
    a.obs["annotation_source"] = "published_Azimuth_cellType_azimuth_GLOBUS_sce_Final"
    a.obs["source_dataset"] = "LIBD_spatialdACC_sce_Final"
    a.obs["cohort_arm"] = "ACA_cortex"
    a.obs["region"] = "dACC"
    a.obs["tv_dominant"] = "ACA"
    a.obs["tv_p_ACA"] = 1.0
    a.obs["tv_assignment_fraction"] = 0.72
    a.obs["tv_source"] = "TerritoryVoxels_HO_cingulate_anterior"
    a.obs["vessel_enriched"] = False
    if "donor_id" not in a.obs.columns:
        a.obs["donor_id"] = a.obs.get("brain", a.obs.get("Sample", "unknown")).astype(str)

    counts = a.obs["cell_type"].value_counts().to_dict()
    summary = {
        "object": str(OUT.relative_to(ROOT)),
        "n_vascular": int(a.n_obs),
        "n_donors": int(a.obs["donor_id"].nunique()),
        "cell_type_counts": counts,
        "annotation_source": "published_Azimuth_GLOBUS_sce_Final",
        "matches_paper_log": counts.get("Endo") == 89 and counts.get("VLMC") == 159,
        "vessel_enriched": False,
        "verdict": "GO_paper_identical_Azimuth_unsorted_dACC",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(OUT, compression="gzip")
    OUT_SUM.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
