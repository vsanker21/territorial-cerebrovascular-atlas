# Aim 4 nest — ADNI territorial WMH + MRA geometry (strict)

**Analysis set:** QC strict n=100 (registry: primary n=155, strict n=100)

## Outcomes
- Local vs UCD TOTAL_WMH Spearman ρ=0.573, p=4.77e-10
- Wilcoxon MCA vs ACA WMH burden: p=0.249

## Nested OLS R² predicting UCD TOTAL_WMH
| Model | R² | n |
|---|---:|---:|
| M0 demographics (age, sex) | 0.029 | 100 |
| M1 MRA geometry null (atlas constants) | 0.000 | 100 |
| M2 geometry + TV WMH burden | 0.449 | 100 |
| M3 demographics + TV WMH burden | 0.477 | 100 |
| **Δ (M3 − M0)** | **0.448** | |

## Territory means (strict)
territory  mean_wmh_burden  median_wmh_burden  mean_wmh_share  mra_vessel_probability_mean  mra_vessel_radius_mean
      ACA         0.002328           0.000647        0.107356                     1.313527                0.239551
      MCA         0.001437           0.000474        0.318235                     0.499102                0.192529
      PCA         0.000602           0.000000        0.060871                     0.290947                0.177866
       VB         0.000209           0.000000        0.011313                     0.334738                0.072408

## Locked interpretation
- Atlas MRA geometry is a *shared* null covariate (not per-subject vascular anatomy); M1 R²≈0 is expected.
- M3 predicting UCD TOTAL_WMH from territorial burdens of the *same* FLAIR WMH map is a **partition / concordance check** (parts vs clinical total), **not** H3 molecular proof.
- Aim 4 molecular claim still requires control-trained vascular transcriptomic features predicting independent territorial WMH or ATLAS lesions.
- Descriptive result: ADNI WMH mass is predominantly MCA-territory share in this QC set.
- Do not claim stroke-lesion vulnerability until ATLAS decrypt + independent lesion maps.

Tables: `results\aim4\adni\aim4_nest_wmh_mra_subject_table_strict.csv`, `aim4_nest_wmh_mra_territory_summary_strict.csv`
