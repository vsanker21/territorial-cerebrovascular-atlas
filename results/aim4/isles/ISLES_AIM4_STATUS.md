# ISLES'22 Aim 4 stroke replication — status

**Date:** 2026-07-27  
**DOI:** [10.5281/zenodo.7960856](https://doi.org/10.5281/zenodo.7960856)

## Pipeline
1. Download + extract ISLES-2022 (done)
2. DWI affine → MNI152 2 mm; nearest-neighbor lesion mask warp
3. Overlay **ArterialAtlas_level2** (L/R → ACA/MCA/PCA/VB)
4. Nest: MRA geometry null + territorial burdens

**Why not TerritoryVoxels_BMM?** Same rule as ATLAS — bilateral stroke geography uses level2.

## Results (primary QC)

| Metric | Value |
|---|---|
| Cases | 250 |
| OK registration | 250 |
| **Primary PASS** | **247** (3 `FAIL_ZERO_LESION` after warp) |
| Strict PASS | 247 |
| Median lesion voxels (2 mm) | 1472 |
| Median assigned frac | 0.959 |
| Mean lesion share MCA / ACA / PCA / VB | **0.486 / 0.065 / 0.126 / 0.238** |
| Dominant territory | MCA 137, VB 66, PCA 33, ACA 11 |
| M1 MRA geometry null R² | 0.000 |
| Wilcoxon MCA vs ACA burden | p ≈ 2.7e-22 |

**vs ATLAS (primary):** ATLAS mean MCA share ≈0.70 (dominant MCA 1114/1445); ISLES MCA share ≈0.49 with more VB-dominant cases (66/247). MCA remains the modal dominant territory and mean MCA > mean ACA in both cohorts. Higher ISLES VB share is reported descriptively (dataset mix / registration), not as molecular H3.

**Locked read:** M3 R²≈0.99 predicting total lesion from territorial burdens is a **partition check**, **not** molecular H3.

## Outputs
- `results/aim4/isles/isles_lesion_tv_fractions.csv`
- `results/aim4/isles/aim4_nest_isles_mra_report.md`
- Figure: `figures/open_agenda/08_isles_aim4_lesion_nest`
- Scripts: `isles_lesion_tv_pipeline.py`, `aim4_nest_isles_mra.py`, `download_isles22.py`

## Cite
Hernandez Petzsche et al., Scientific Data 2022 (ISLES'22); Zenodo 7960856.
