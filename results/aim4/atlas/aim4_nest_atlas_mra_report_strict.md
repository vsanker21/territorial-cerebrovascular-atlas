# Aim 4 nest — ATLAS territorial lesion + MRA geometry (strict)

**Analysis set:** QC strict n=1444 (registry primary n=1445, strict n=1444)

- Median lesion voxels: 7455
- Dominant territory counts: {'MCA': 1114, 'VB': 139, 'PCA': 113, 'ACA': 78}
- Wilcoxon MCA vs ACA lesion burden: p=1.59e-121

## Nested OLS R² predicting total lesion voxels
| Model | R² | n |
|---|---:|---:|
| M0 chronicity | 0.000 | 1444 |
| M1 MRA geometry null | 0.000 | 1444 |
| M2 geometry + TV lesion burden | 0.998 | 1444 |
| M3 chronicity/burden | 0.998 | 1444 |
| **Δ (M3 − M0)** | **0.998** | |

## Territory means (strict)
territory  mean_lesion_burden  median_lesion_burden  mean_lesion_share  mra_vessel_probability_mean  mra_vessel_radius_mean
      ACA            0.011036              0.000187           0.095924                     1.313527                0.239551
      MCA            0.032699              0.005453           0.695637                     0.499102                0.192529
      PCA            0.006742              0.000068           0.091730                     0.290947                0.177866
       VB            0.002547              0.000000           0.097386                     0.334738                0.072408

## Locked interpretation
- Atlas MRA geometry is a shared null; M1 R²≈0 is expected.
- M3 predicting total lesion volume from territorial burdens of the *same* lesion map is a **partition-style check**, **not** molecular H3.
- Molecular Aim 4 claim still requires control-trained transcriptomic features × independent territorial lesions.
- Report descriptive stroke geography (shares / dominant territory) with ATLAS citation.

Tables: `results\aim4\atlas\aim4_nest_atlas_mra_subject_table_strict.csv`
