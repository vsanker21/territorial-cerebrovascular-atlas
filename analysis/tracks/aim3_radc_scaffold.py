"""
RADC Aim 3 linkage scaffolding — DUA-gated.

Builds the analysis skeleton against Garcia ROSMAP projids already local.
Runs ONLY when RADC Standard extract lands under Project #7508.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import anndata as ad
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "tracks" / "aim3_radc"
OUT.mkdir(parents=True, exist_ok=True)
H5AD = ROOT / "data/processed/garcia/garcia_rosmap_vascular.h5ad"
# Expected delivery path once DUA executes
RADC_CANDIDATES = [
    ROOT / "data/raw/radc/Project7508_Standard_Clinical_Pathologic.csv",
    ROOT / "data/raw/radc/standard_clinical_pathologic_variables.csv",
    ROOT / "data/raw/radc/rosmap_standard.csv",
]


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def main():
    a = ad.read_h5ad(H5AD, backed="r")
    # Recover projid from garcia sample table if present
    sample_csv = ROOT / "results/metadata_audit/garcia_rosmap_samples.csv"
    if sample_csv.exists():
        samp = pd.read_csv(sample_csv)
        donors = samp[["projid", "study"]].drop_duplicates() if "projid" in samp.columns else pd.DataFrame()
    else:
        donors = pd.DataFrame({"subject": a.obs["subject"].astype(str).unique()})

    donor_list = {
        "n_garcia_rosmap_subjects_in_h5ad": int(a.obs["subject"].nunique()),
        "n_projid_in_sample_table": int(donors["projid"].nunique()) if "projid" in donors.columns else None,
        "sample_table": str(sample_csv) if sample_csv.exists() else None,
    }
    if "projid" in donors.columns:
        donors.to_csv(OUT / "aim3_garcia_projid_key.csv", index=False)

    radc_path = next((p for p in RADC_CANDIDATES if p.exists()), None)
    status = {
        "date": time.strftime("%Y-%m-%d"),
        "dua_project": "RADC #7508",
        "rra": "RRA-539487",
        "data_received": radc_path is not None,
        "radc_path": str(radc_path) if radc_path else None,
        "garcia_key": donor_list,
        "planned_models": [
            "territorial_molecular_score ~ clinical diagnosis (NCI/MCI/AD) + age + sex + educ + region",
            "territorial_molecular_score ~ Braak/CERAD/NIA-Reagan + covariates",
            "territorial_molecular_score ~ arteriolosclerosis/atherosclerosis/CAA/infarcts + covariates",
            "APOE as modifier",
        ],
        "verdict": (
            "READY_TO_LINK"
            if radc_path
            else "GATED_AWAITING_RADC_DUA_DELIVERY__SCAFFOLD_ONLY"
        ),
        "not_identity_proof": (
            "Aim 3 is phenotype linkage / association — not developmental territorial identity (Track B)."
        ),
    }

    if radc_path:
        radc = pd.read_csv(radc_path, low_memory=False)
        status["radc_n_rows"] = int(len(radc))
        status["radc_columns_head"] = list(radc.columns)[:40]
        # attempt link
        idcol = next((c for c in ["projid", "individualID", "proj_id"] if c in radc.columns), None)
        if idcol and "projid" in donors.columns:
            merged = donors.merge(radc, left_on="projid", right_on=idcol, how="inner")
            merged.to_csv(OUT / "aim3_linked_table.csv", index=False)
            status["n_linked"] = int(len(merged))
            status["verdict"] = "LINKED__RUN_PHENOTYPE_MODELS"

    with open(OUT / "aim3_radc_status.json", "w", encoding="utf-8") as f:
        json.dump(status, f, indent=2)

    md = [
        "# Aim 3 — RADC phenotype linkage (scaffold)",
        "",
        f"**Verdict:** `{status['verdict']}`",
        f"- RRA-539487 / Project #7508",
        f"- Garcia subjects in h5ad: {donor_list['n_garcia_rosmap_subjects_in_h5ad']}",
        f"- projids in sample table: {donor_list['n_projid_in_sample_table']}",
        "",
        status["not_identity_proof"],
        "",
        "## When data lands",
        "1. Place Standard extract under `data/raw/radc/`.",
        "2. Re-run this script → `aim3_linked_table.csv`.",
        "3. Fit pre-specified ordinal/mixed models (diagnosis, AD path, CVD path).",
        "4. Do **not** claim Track B developmental identity from Aim 3 alone.",
        "",
    ]
    (OUT / "aim3_radc_report.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
