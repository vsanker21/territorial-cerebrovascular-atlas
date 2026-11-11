# Aim 4 nest — ADNI territorial WMH + MRA geometry (primary)

**Analysis set:** QC primary n=155 (registry: primary n=155, strict n=100)

## Outcomes
- Local vs UCD TOTAL_WMH Spearman ρ=0.521, p=3.57e-12
- Wilcoxon MCA vs ACA WMH burden: p=0.217

## Nested OLS R² predicting UCD TOTAL_WMH
| Model | R² | n |
|---|---:|---:|
| M0 demographics (age, sex) | 0.035 | 155 |
| M1 MRA geometry null (atlas constants) | 0.000 | 155 |
| M2 geometry + TV WMH burden | 0.434 | 155 |
| M3 demographics + TV WMH burden | 0.459 | 155 |
| **Δ (M3 − M0)** | **0.424** | |

## Territory means (primary)
territory  mean_wmh_burden  median_wmh_burden  mean_wmh_share  mra_vessel_probability_mean  mra_vessel_radius_mean
      ACA         0.001528           0.000000        0.070538                     1.313527                0.239551
      MCA         0.000932           0.000086        0.207277                     0.499102                0.192529
      PCA         0.000401           0.000000        0.039610                     0.290947                0.177866
       VB         0.000146           0.000000        0.007781                     0.334738                0.072408

## Locked interpretation
- Atlas MRA geometry is a *shared* null covariate (not per-subject vascular anatomy); M1 R²≈0 is expected.
- M3 predicting UCD TOTAL_WMH from territorial burdens of the *same* FLAIR WMH map is a **partition / concordance check** (parts vs clinical total), **not** H3 molecular proof.
- Aim 4 molecular claim still requires control-trained vascular transcriptomic features predicting independent territorial WMH or ATLAS lesions.
- Descriptive result: ADNI WMH mass is predominantly MCA-territory share in this QC set.
- Do not claim stroke-lesion vulnerability until ATLAS decrypt + independent lesion maps.

Tables: `results\aim4\adni\aim4_nest_wmh_mra_subject_table.csv`, `aim4_nest_wmh_mra_territory_summary.csv`
