# Moonshot 01 — GWAS × MultiVINE territory filter

**Honesty lock:** Falsify if territory adds nothing beyond cell type + A–V zonation

- MultiVINE vascular nuclei: **39582** (dlPFC only)
- Genesets: {'AD_core': 20, 'cerebrovascular_combined': 132, 'cSVD_core': 7, 'ICH': 5, 'ischemic_stroke_HP': 53, 'stroke_core': 16, 'WMH': 38, 'WMH_core': 13}
- Falsifier: `{"class_b_gwas_score_donor_oof_auroc_ACA_vs_MCA": 0.9160373325733828, "null_random_modules_mean_auroc": 0.9157826303790255, "null_random_modules_sd_auroc": 0.009878438313378937, "gwas_beats_null": false, "garcia_MT_vs_EC_gwas_donor_oof_auroc": 0.48773694390715666}`

## Top enrichments
| cell_type            | geneset                  |   overlap_n |   hypergeom_p |        fdr |
|:---------------------|:-------------------------|------------:|--------------:|-----------:|
| Smooth_Muscle_Cell   | WMH_core                 |           4 |   0.00022115  | 0.0017692  |
| Smooth_Muscle_Cell   | cSVD_core                |           3 |   0.000493419 | 0.00394735 |
| Smooth_Muscle_Cell   | cerebrovascular_combined |           6 |   0.111035    | 0.888277   |
| Pericyte             | cSVD_core                |           1 |   0.162521    | 0.650085   |
| Pericyte             | WMH_core                 |           1 |   0.280812    | 1          |
| Smooth_Muscle_Cell   | stroke_core              |           1 |   0.333594    | 1          |
| Pericyte             | stroke_core              |           1 |   0.333594    | 1          |
| Arterial             | ischemic_stroke_HP       |           2 |   0.384164    | 1          |
| Pericyte             | ischemic_stroke_HP       |           2 |   0.384164    | 1          |
| Capillary            | AD_core                  |           1 |   0.398048    | 1          |
| Meningeal_Fibroblast | AD_core                  |           1 |   0.398048    | 1          |
| Arterial_vs_Venous   | AD_core                  |           1 |   0.398048    | 1          |

## Verdict notes
- FALSIFIED at Class B: GWAS modules do not beat random modules for ACA vs MCA — territory adds nothing beyond noise at this resolution.

Outputs: `gwas_enrichment_by_celltype.csv`, `gwas_projection_by_territory.csv`, `summary.json`