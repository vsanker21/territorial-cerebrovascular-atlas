"""Shared paths and honesty locks for the public-data moonshot stack."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(r"F:\Cerebrovascular SCS")
DATA = ROOT / "data"
RAW = DATA / "raw"
PROC = DATA / "processed"
RES = ROOT / "results" / "moonshot"
RES.mkdir(parents=True, exist_ok=True)

# Inputs
MULTIVINE = RAW / "multivine" / "multivine_snrna.h5ad"
GWAS_DIR = RAW / "gwas"
MATTHEWS_H5 = RAW / "matthews_bve" / "raw"
GARCIA_TEMPORAL = PROC / "garcia" / "garcia_temporal_mca_qc_tight.h5ad"
GARCIA_TEMPORAL_FALLBACK = PROC / "garcia" / "garcia_temporal_mca.h5ad"
CLASS_B = PROC / "class_b" / "aca_mca_harmony_endo.h5ad"
ALLEN_ACC = PROC / "aca_cohort_hunt" / "allen_acc_vascular.h5ad"
ALLEN_MULTI = PROC / "allen" / "allen_endo_multiregion.h5ad"
LINNARSSON = PROC / "aca_cohort_hunt" / "linnarsson_aca_vascular.h5ad"
YANG = PROC / "yang" / "yang_vascular_qc_tight.h5ad"
WALCHLI_FETAL = (
    RAW
    / "walchli"
    / "geo_GSE256493"
    / "suppl"
    / "GSE256493_Fetal_CNS_sorted_endothelial_cells_seurat_object.rds.gz"
)
WALCHLI_ADULT = (
    RAW
    / "walchli"
    / "geo_GSE256493"
    / "suppl"
    / "_decomp"
    / "adult_temporal_sorted.rds"
)
ISLES24_TABLE = (
    ROOT / "results" / "nature_gates" / "isles24_trackC" / "isles24_trackC_subject_table.csv"
)
ATLAS_FRAC = ROOT / "results" / "aim4" / "atlas" / "atlas_lesion_tv_fractions.csv"
ISLES_FRAC = ROOT / "results" / "aim4" / "isles" / "isles_lesion_tv_fractions.csv"
ADNI_WMH = ROOT / "results" / "aim4" / "adni" / "adni_wmh_t1_mni_scaled.csv"
MODULE_BY_DONOR = ROOT / "results" / "class_b" / "class_b_module_scores_by_donor.csv"
GARCIA_REGIONS = PROC / "garcia" / "garcia_rosmap_vascular.h5ad"  # Sun vascular regions in this build

HONESTY = {
    "1": "Falsify if territory adds nothing beyond cell type + A–V zonation",
    "2": "SENSITIVITY_NOT_CONFIRMATORY_H1c — never report transported ACC as markers",
    "3": "Kill Track B language if adult territory ⊥ developmental state after covariates",
    "4": "Only nested OOS / conformal metrics count; in-sample R² is descriptive",
    "5": "Theory figure — pair with genetics/imaging; not standalone claim",
    "6": "ECOLOGICAL_NOT_CONFIRMATORY_MOLECULAR_H3",
    "7": "Methods lever — with #1 as supplement; alone → Nat Methods",
}
