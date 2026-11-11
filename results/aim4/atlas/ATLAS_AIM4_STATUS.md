# ATLAS Aim 4 stroke arm — status

**Date:** 2026-07-27

## Pipeline
1. Decrypt + extract ATLAS R3.0 preprocessed (done)
2. Lesion → **ArterialAtlas_level2** (L/R collapsed to ACA/MCA/PCA/VB)
3. Nest: MRA geometry null + territorial lesion burdens

**Why not TerritoryVoxels_BMM?** In this NITRC release, TV_BMM labels 1–4 are effectively unilateral; contralesional lesions fall in label 6. `ArterialAtlas_level2` is the bilateral hierarchical major-territory map (`ArterialAtlasLabels.txt`).

## Results (primary QC)

| Metric | Value |
|---|---|
| Masks | 1453 (1 empty/corrupt) |
| OK | 1452 |
| **Primary PASS** | **1445** |
| Strict PASS | 1444 |
| Zero-lesion fail | 7 |
| Median lesion voxels | 7449 |
| Median assigned frac | 0.997 |
| Mean lesion share MCA / ACA / PCA / VB | **0.695 / 0.096 / 0.092 / 0.097** |
| Dominant territory | MCA 1114, VB 139, PCA 113, ACA 78 |
| M1 MRA geometry null R² | 0.000 |
| Wilcoxon MCA vs ACA burden | p ≈ 1.6e-121 |

**Locked read:** M3 R²≈0.998 predicting total lesion volume from territorial burdens is a **partition check** (same lesion map), **not** molecular H3. Descriptive geography: MCA-dominant stroke lesions.

## Laterality sensitivity
L/R lesion shares ≈ balanced (mean L≈0.49 / R≈0.49; dominant L 728 / R 717). MCA remains modal within both L- and R-dominant subjects → collapsing L/R→major is not a unilateral-atlas artifact. Report: `atlas_laterality_report.md`.

## Exploratory molecular × ATLAS
Allen Endo modules (ACA/MCA/PCA) vs ATLAS lesion shares: `EXPLORATORY_ECOLOGICAL_n3__NOT_CONFIRMATORY_H3` (Spearman |ρ|≈0.5, p≈0.67). Report: `aim4_molecular_atlas_nest_report.md`. Confirmatory H3 still needs enrichment-matched ACC vascular features.

## Outputs
- `results/aim4/atlas/atlas_lesion_tv_fractions.csv`
- `results/aim4/atlas/aim4_nest_atlas_mra_report.md`
- `results/aim4/atlas/aim4_molecular_atlas_nest_report.md`
- `results/aim4/atlas/atlas_laterality_report.md`
- Scripts: `atlas_lesion_tv_pipeline.py`, `aim4_nest_atlas_mra.py`, `aim4_molecular_atlas_nest.py`, `atlas_laterality_sensitivity.py`

## Cite ATLAS
Liew et al., Scientific Data 2018; Liew et al. ATLAS R2/R3 releases as required by Terms of Use.
