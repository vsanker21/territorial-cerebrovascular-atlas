# Aim 4 nest — ATLAS territorial lesion + MRA geometry (primary)

**Analysis set:** QC primary n=1445 (registry primary n=1445, strict n=1444)

- Median lesion voxels: 7449
- Dominant territory counts: {'MCA': 1114, 'VB': 139, 'PCA': 113, 'ACA': 78, 'nan': 1}
- Wilcoxon MCA vs ACA lesion burden: p=1.59e-121

## Nested OLS R² predicting total lesion voxels
| Model | R² | n |
|---|---:|---:|
| M0 chronicity | 0.000 | 1445 |
| M1 MRA geometry null | 0.000 | 1445 |
| M2 geometry + TV lesion burden | 0.998 | 1445 |
| M3 chronicity/burden | 0.998 | 1445 |
| **Δ (M3 − M0)** | **0.998** | |

## Territory means (primary)
territory  mean_lesion_burden  median_lesion_burden  mean_lesion_share  mra_vessel_probability_mean  mra_vessel_radius_mean
      ACA            0.011028              0.000187           0.095858                     1.313527                0.239551
      MCA            0.032676              0.005441           0.695156                     0.499102                0.192529
      PCA            0.006737              0.000067           0.091666                     0.290947                0.177866
       VB            0.002545              0.000000           0.097319                     0.334738                0.072408

## Locked interpretation
- Atlas MRA geometry is a shared null; M1 R²≈0 is expected.
- M3 predicting total lesion volume from territorial burdens of the *same* lesion map is a **partition-style check**, **not** molecular H3.
- Molecular Aim 4 claim still requires control-trained transcriptomic features × independent territorial lesions.
- Report descriptive stroke geography (shares / dominant territory) with ATLAS citation.

Tables: `results\aim4\atlas\aim4_nest_atlas_mra_subject_table.csv`
