# B0 gates + Track C (flow programs)

**Date:** 2026-08-17
**Overall routing:** `B0i_PASS_COUNTS__B0ii_HOX_RNA_NOGO__TRACK_B_NEEDS_CHROMATIN__TRACK_C_PRIMARY`

## Data
- Source: `F:\Cerebrovascular SCS\data\processed\garcia\garcia_rosmap_vascular.h5ad`
- Note: Local Garcia ROSMAP multiregion vascular object used for B0/Track C. Sun 2023 MIT counts_full.rds download was incomplete/unreadable (exact 32MiB); full Sun release requires Dropbox/Synapse (syn51015750). Regions overlap Sun design; cell n is smaller than published 22,514.
- n_cells=11689, n_donors=214

## B0i — mural counts (Per+SMC) per region per donor
- Verdict: `PASS_MURAL_COUNTS_ADEQUATE_AMONG_REGION_PRESENT_DONORS`
- Medians: `{"Angular_gyrus": 9.0, "Anterior_thalamus": 12.0, "Entorhinal_cortex": 8.0, "Hippocampus": 8.0, "Mammillary_body": 1.5, "Midtemporal_cortex": 8.0, "Prefrontal_cortex": 8.0}`
- Matched pair: `{"Midtemporal_cortex": {"n_donors_with_region": 25, "n_donors_with_mural": 25, "median_mural": 8.0, "mean_mural": 9.0, "pct_donors_ge_5": 0.8, "pct_donors_ge_10": 0.4}, "Entorhinal_cortex": {"n_donors_with_region": 24, "n_donors_with_mural": 24, "median_mural": 8.0, "mean_mural": 8.166666666666666, "pct_donors_ge_5": 0.7083333333333334, "pct_donors_ge_10": 0.3333333333333333}, "Prefrontal_cortex": {"n_donors_with_region": 213, "n_donors_with_mural": 200, "median_mural": 8.0, "mean_mural": 12.23943661971831, "pct_donors_ge_5": 0.6854460093896714, "pct_donors_ge_10": 0.4272300469483568}, "Hippocampus": {"n_donors_with_region": 41, "n_donors_with_mural": 38, "median_mural": 8.0, "mean_mural": 7.804878048780488, "pct_donors_ge_5": 0.6585365853658537, "pct_donors_ge_10": 0.2926829268292683}, "Anterior_thalamus": {"n_donors_with_region": 25, "n_donors_with_mural": 24, "median_mural": 12.0, "mean_mural": 13.16, "pct_donors_ge_5": 0.8, "pct_donors_ge_10": 0.64}}`

## B0ii — HOX anterior vs posterior (lineage proxy)
- Verdict: `HOX_RNA_UNDETECTABLE_IN_ADULT_MURAL__LINEAGE_NEEDS_CHROMATIN__B0ii_RNA_NOGO`
- Donor Wilcoxon: `{"n_paired_donors": 39, "mean_anterior": 0.0, "mean_posterior": 0.0, "wilcoxon_stat": 0.0, "p": NaN, "direction_expected": "posterior_HOX_gt_anterior", "direction_observed": "anterior_ge_posterior"}`
- Matched pair: `{"pair": "Midtemporal_cortex(MCA) vs Entorhinal_cortex(PCA)", "n_paired_donors": 24, "mean_midtemporal": 0.0, "mean_entorhinal": 0.0, "wilcoxon_p": NaN, "direction_observed": "midtemporal_ge_entorhinal"}`

## Track C — endothelial flow/shear programs
- Verdict: `FLOW_PROGRAM_NULL_BETWEEN_POLES__TRACK_C_INCONCLUSIVE_ON_THIS_OBJECT`
- Pole: `{"mean_anterior": 0.48309700540389444, "mean_posterior": 0.4752120140656605, "p": 0.8246719895734236, "n_ant": 245, "n_post": 99}`
- MCA vs PCA regions: `{"mean_MCA_regions": 0.4963113734019209, "mean_PCA_regions": 0.467091489228464, "p": 0.3800882186181823}`

## Track B gate
- `{"vessel_enriched_ACC": "CLOSED_Heiman_Matthews__PENDING_Walchli", "sun_full_release": "BLOCKED_incomplete_MIT_RDS__need_Dropbox_or_Synapse", "chromatin": "NOT_AVAILABLE", "verdict": "TRACK_B_CLOSED__DO_NOT_SPEND"}`

## Locks
- TRCI/BCI not mixed into these tracks.
- Track A imaging partition reported separately.
- RADC Aim 3 still DUA-gated.
