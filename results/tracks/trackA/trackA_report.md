# Track A — volume-normalized burden + geometric null + variance partition

**Date:** 2026-08-17
**Verdict:** `SHARED_GEOMETRIC_TERRITORIAL_DOMINATES__TRACK_C_ALIGNED`

## Design locks
- Outcome: volume-normalized territorial burden (lesion or WMH voxels / territory voxels).
- Geometry null (coarse): log territory volume, mean distance to ventricles, GM/WM means,
  MRA vessel probability & radius (atlas-level).
- Partition: unique geometry / unique territory / shared (commonality analysis).
- **Not** molecular H3. **Not** TRCI/BCI.

## Selection structure (A4)
- **ATLAS:** Chronic / rehabilitation-skewed sampling favors large MCA-territory lesions; agreement with acute cohorts is not fully independent confirmation of MCA biology.
- **ISLES22:** Acute DWI lesions; LVO/thrombectomy selection also MCA-skewed. Higher VB-dominant fraction than ATLAS is descriptive cohort mix.
- **ADNI_WMH:** Non-stroke SVD/WMH axis; TerritoryVoxels geography; not stroke selection. Mean MCA>ACA share claim; not majority-subject claim.

## Variance partition summary

| cohort   |   n_subjects |   R2_geometry |   R2_territory |   R2_joint |   unique_geometry |   unique_territory |    shared |
|:---------|-------------:|--------------:|---------------:|-----------:|------------------:|-------------------:|----------:|
| ATLAS    |         1445 |     0.236454  |      0.236454  |  0.236454  |       0           |                  0 | 0.236454  |
| ISLES22  |          247 |     0.0623765 |      0.0623765 |  0.0623765 |       0           |                  0 | 0.0623765 |
| ADNI_WMH |          155 |     0.122923  |      0.122923  |  0.122923  |      -1.11022e-16 |                  0 | 0.122923  |

## ATLAS
- n subjects: 1445
- geometry features: log_territory_voxels, mean_vent_dist_mm, gm_mean, wm_mean, gm_wm_ratio, mra_vessel_probability, mra_vessel_radius
- partition: `{"R2_geometry": 0.23645359457738124, "R2_territory": 0.23645359457738124, "R2_joint": 0.23645359457738124, "unique_geometry": 0.0, "unique_territory": 0.0, "shared": 0.23645359457738124}`
- MCA vs ACA: `{"stat": 1496288.0, "p": 1.0905871490624952e-91, "test": "mannwhitneyu_MCA_vs_ACA"}`

## ISLES22
- n subjects: 247
- geometry features: log_territory_voxels, mean_vent_dist_mm, gm_mean, wm_mean, gm_wm_ratio, mra_vessel_probability, mra_vessel_radius
- partition: `{"R2_geometry": 0.062376546833225244, "R2_territory": 0.062376546833225244, "R2_joint": 0.062376546833225244, "unique_geometry": 0.0, "unique_territory": 0.0, "shared": 0.062376546833225244}`
- MCA vs ACA: `{"stat": 43505.0, "p": 5.2977377629737584e-17, "test": "mannwhitneyu_MCA_vs_ACA"}`

## ADNI_WMH
- n subjects: 155
- geometry features: log_territory_voxels, mean_vent_dist_mm, gm_mean, wm_mean, gm_wm_ratio, mra_vessel_probability, mra_vessel_radius
- partition: `{"R2_geometry": 0.12292285848812134, "R2_territory": 0.12292285848812146, "R2_joint": 0.12292285848812134, "unique_geometry": -1.1102230246251565e-16, "unique_territory": 0.0, "shared": 0.12292285848812146}`
- MCA vs ACA: `{"stat": 12895.0, "p": 0.2280535448551504, "test": "mannwhitneyu_MCA_vs_ACA"}`

## Interpretation
- **Design lock:** Atlas-constant geometry is linearly equivalent to territory dummies (4 territories → rank ≤3). Identical R²(G)=R²(T)=R²(GT) and zero unique components are expected, not a coding error.
- **What Track A can claim:** volume-normalized MCA>ACA burden in stroke cohorts; shared geometric–territorial structure; selection bias documented.
- **What Track A cannot claim:** territory explains burden *beyond* geometry using atlas-only covariates — that requires subject-level architecture or molecular predictors.
- Large **shared** component aligns with Track C (architecture / developmental branching collinearity).
- **Not** molecular H3. **Not** TRCI/BCI.
