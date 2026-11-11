"""Extract published-annotation ACC endothelial/VLMC cells from Allen Cortical Atlas."""
from __future__ import annotations

import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]
NN = ROOT / "data/raw/aca_cohort_hunt/allen_cortical_acc/NonNeuronal.h5ad"
ACC = ROOT / "data/raw/aca_cohort_hunt/allen_cortical_acc/ACC_dissection.h5ad"
OUT = ROOT / "data/processed/aca_cohort_hunt/allen_acc_vascular.h5ad"
OUT_SUM = ROOT / "results/aca_cohort_hunt/allen_acc_vascular_summary.json"

VASC = {
    "endothelial cell",
    "vascular leptomeningeal cell",
    "pericyte",
    "smooth muscle cell",
    "brain vascular cell",
}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT_SUM.parent.mkdir(parents=True, exist_ok=True)

    # Prefer NonNeuronal (already vascular-rich) filtered to ACC tissue
    log(f"Reading {NN.name}")
    a = sc.read_h5ad(NN)
    log(f"n={a.n_obs} cols={list(a.obs.columns)[:25]}")
    tissue_col = None
    for c in ["tissue", "tissue_ontology_term_id", "region", "dissection"]:
        if c in a.obs:
            tissue_col = c
            break
    ct_col = "cell_type" if "cell_type" in a.obs else None
    print("tissue_col", tissue_col, "ct_col", ct_col)
    if tissue_col:
        print(a.obs[tissue_col].astype(str).value_counts().head(20))
    if ct_col:
        print(a.obs[ct_col].astype(str).value_counts().head(30))

    # ACC filter
    t = a.obs[tissue_col].astype(str).str.lower()
    acc_mask = t.str.contains("anterior cingulate|cingulate")
    log(f"ACC tissue cells: {int(acc_mask.sum())}")
    b = a[acc_mask].copy()
    ct = b.obs[ct_col].astype(str).str.lower()
    # match flexible vascular names
    vasc_mask = (
        ct.str.contains("endo")
        | ct.str.contains("pericyt")
        | ct.str.contains("vlmc")
        | ct.str.contains("leptomeningeal")
        | ct.str.contains("smooth muscle")
        | ct.isin({x.lower() for x in VASC})
    )
    v = b[vasc_mask].copy()
    log(f"ACC vascular: {v.n_obs}")
    print(v.obs[ct_col].astype(str).value_counts())

    # If sparse, fall back to full ACC dissection Endo
    if v.n_obs < 50:
        log("Too few from NonNeuronal ACC; using ACC dissection")
        a2 = sc.read_h5ad(ACC)
        ct2 = a2.obs[ct_col].astype(str).str.lower()
        vasc2 = (
            ct2.str.contains("endo")
            | ct2.str.contains("pericyt")
            | ct2.str.contains("vlmc")
            | ct2.str.contains("leptomeningeal")
            | ct2.str.contains("smooth muscle")
        )
        v = a2[vasc2].copy()
        log(f"ACC dissection vascular: {v.n_obs}")
        print(v.obs[ct_col].astype(str).value_counts())

    v.obs["source_dataset"] = "Allen_Cortical_Cell_Atlas_Bakken"
    v.obs["cohort_arm"] = "ACA_cortex"
    v.obs["tv_dominant"] = "ACA"
    v.obs["tv_p_ACA"] = 1.0
    v.obs["tv_assignment_fraction"] = 0.72
    v.obs["tv_source"] = "TerritoryVoxels_HO_cingulate_anterior"
    v.obs["vessel_enriched"] = False
    v.obs["annotation_source"] = "published_Allen_taxonomy"

    donor_col = "donor_id" if "donor_id" in v.obs else None
    summary = {
        "object": str(OUT.relative_to(ROOT)),
        "n_vascular": int(v.n_obs),
        "n_donors": int(v.obs[donor_col].nunique()) if donor_col else None,
        "cell_type_counts": v.obs[ct_col].astype(str).value_counts().to_dict(),
        "vessel_enriched": False,
        "anatomy": "anterior_cingulate_cortex",
        "tv_dominant": "ACA",
        "verdict": "GO_annotated_ACA_unsorted_not_vessel_enriched",
    }
    v.write_h5ad(OUT, compression="gzip")
    OUT_SUM.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
