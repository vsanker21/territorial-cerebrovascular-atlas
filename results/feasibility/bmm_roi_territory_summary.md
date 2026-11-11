# BMM territory feasibility — confirmatory

Using real `ProbArterialAtlas_BMM` from NITRC Atlas_182_MNI152.

## Verdict

```json
{
  "method": "Liu ProbArterialAtlas_BMM (Atlas_182_MNI152)",
  "affine_info": {
    "bmm_affine": [
      [
        -1.0,
        0.0,
        0.0,
        90.0
      ],
      [
        0.0,
        1.0,
        0.0,
        -126.0
      ],
      [
        0.0,
        0.0,
        1.0,
        -72.0
      ],
      [
        0.0,
        0.0,
        0.0,
        1.0
      ]
    ],
    "bmm_shape": [
      182,
      218,
      182,
      4
    ],
    "border_shape": [
      182,
      218,
      182,
      2
    ],
    "territory_voxels_unique": [
      0,
      1,
      1,
      3,
      3,
      5,
      5
    ]
  },
  "hippocampus_pca_dominant": true,
  "superior_frontal_aca_dominant": false,
  "yang_within_donor_viable": false,
  "recommended_yang_contrast": "PCA (hippocampus) vs MCA (superior frontal)",
  "garcia_temporal_mca": true
}
```

## ROI table

| roi                    | role                 |   n_voxels |      p_ACA |    p_MCA |       p_PCA |        p_VB | dominant   |   p_dominant |   entropy_bits |   mean_MCA_over_ACA_ratio |   mean_MCA_over_PCA_ratio |   frac_multi_claimed_voxels |   frac_unclaimed_voxels | high_specificity_pdom_ge_0.50   | high_specificity_H_le_1.0   | method                               |
|:-----------------------|:---------------------|-----------:|-----------:|---------:|------------:|------------:|:-----------|-------------:|---------------:|--------------------------:|--------------------------:|----------------------------:|------------------------:|:--------------------------------|:----------------------------|:-------------------------------------|
| hippocampus            | Yang                 |       5527 | 0          | 0.31278  | 0.686681    | 0.000538966 | PCA        |     0.686681 |     0.902694   |                       nan |                       nan |                           0 |                       0 | True                            | True                        | ProbArterialAtlas_BMM_mean_in_HO_ROI |
| superior_frontal_gyrus | Yang                 |      21887 | 0.456797   | 0.543203 | 0           | 0           | MCA        |     0.543203 |     0.994608   |                       nan |                       nan |                           0 |                       0 | True                            | True                        | ProbArterialAtlas_BMM_mean_in_HO_ROI |
| temporal_pole          | Garcia_exvivo_MCA    |      18571 | 0          | 1        | 0           | 0           | MCA        |     1        |    -0          |                       nan |                       nan |                           0 |                       0 | True                            | True                        | ProbArterialAtlas_BMM_mean_in_HO_ROI |
| middle_frontal_gyrus   | ClassB_MCA_candidate |      21956 | 0.00484456 | 0.995155 | 0           | 0           | MCA        |     0.995155 |     0.0442241  |                       nan |                       nan |                           0 |                       0 | True                            | True                        | ProbArterialAtlas_BMM_mean_in_HO_ROI |
| frontal_medial_cortex  | ClassB_ACA_candidate |       4802 | 0.0588822  | 0.941118 | 0           | 0           | MCA        |     0.941118 |     0.322992   |                       nan |                       nan |                           0 |                       0 | True                            | True                        | ProbArterialAtlas_BMM_mean_in_HO_ROI |
| angular_gyrus          | ClassB_MCA_candidate |       7584 | 0          | 0.999915 | 8.54678e-05 | 0           | MCA        |     0.999915 |     0.00127833 |                       nan |                       nan |                           0 |                       0 | True                            | True                        | ProbArterialAtlas_BMM_mean_in_HO_ROI |
| occipital_pole         | ClassA_PCA           |      20793 | 0          | 0.123185 | 0.801367    | 0.0754482   | PCA        |     0.801367 |     0.909461   |                       nan |                       nan |                           0 |                       0 | True                            | True                        | ProbArterialAtlas_BMM_mean_in_HO_ROI |
