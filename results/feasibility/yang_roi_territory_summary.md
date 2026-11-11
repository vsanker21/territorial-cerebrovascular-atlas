# Yang → territory feasibility probe — results

## Methods note

- Territory vectors = voxel-fraction of Liu hard atlas labels inside Harvard-Oxford ROIs,
  resampled (nearest) into ArterialAtlas space.
- This is a proxy until NITRC `ProbArterialAtlas_BMM` / borderzone maps are ingested.
- Yang cortex = **superior frontal cortex** (not generic cortex).

## Pre-specified go criteria

- Hippocampus: PCA-dominant with P≥0.50
- Superior frontal: anterior-circulation (ACA+MCA) ≥0.70 and PCA <0.20
- Within-donor Yang contrast viable if both hold

## Verdict

```json
{
  "hippocampus_pca_dominant": true,
  "superior_frontal_anterior_circulation": true,
  "superior_frontal_single_territory_dominant": true,
  "within_donor_contrast_viable": true,
  "draft_plan_pca_vs_mca_framing_accurate": false,
  "recommended_yang_contrast_framing": "PCA (hippocampus) vs ACA (superior frontal)"
}
```

## Yang ROI territory overlap

| roi                    | labels                             |   n_voxels_in_atlas | dominant_l2   |   p_dominant_l2 |   entropy_l2_bits |    p_ACA |     p_MCA |    p_PCA |       p_VB |   p_anterior_ACA_MCA | dominant_l1_majorized   |   p_dominant_l1_majorized |   entropy_l1_majorized_bits | high_specificity_pdom_ge_0.50   | high_specificity_H_le_1.0   | yang_role                                     | a_priori_expectation         |
|:-----------------------|:-----------------------------------|--------------------:|:--------------|----------------:|------------------:|---------:|----------:|---------:|-----------:|---------------------:|:------------------------|--------------------------:|----------------------------:|:--------------------------------|:----------------------------|:----------------------------------------------|:-----------------------------|
| hippocampus            | Left Hippocampus;Right Hippocampus |               11007 | PCA           |        0.757427 |          0.841281 | 0        | 0.236304  | 0.757427 | 0.00626874 |             0.236304 | PCA                     |                  0.757427 |                    0.841281 | True                            | True                        | Yang primary: hippocampus samples             | PCA-dominant                 |
| superior_frontal_gyrus | Superior Frontal Gyrus             |               35205 | ACA           |        0.926204 |          0.379934 | 0.926204 | 0.0737963 | 0        | 0          |             1        | ACA                     |                  0.926204 |                    0.379934 | True                            | True                        | Yang primary: superior frontal cortex samples | ACA/MCA mixed (not pure MCA) |

## Aim-2 contrast ROI candidates

| roi                    | labels                                                              |   n_voxels_in_atlas | dominant_l2   |   p_dominant_l2 |   entropy_l2_bits |     p_ACA |       p_MCA |    p_PCA |     p_VB |   p_anterior_ACA_MCA | dominant_l1_majorized   |   p_dominant_l1_majorized |   entropy_l1_majorized_bits | high_specificity_pdom_ge_0.50   | high_specificity_H_le_1.0   | contrast_role                                                          |
|:-----------------------|:--------------------------------------------------------------------|--------------------:|:--------------|----------------:|------------------:|----------:|------------:|---------:|---------:|---------------------:|:------------------------|--------------------------:|----------------------------:|:--------------------------------|:----------------------------|:-----------------------------------------------------------------------|
| occipital_pole         | Occipital Pole                                                      |               39110 | PCA           |        0.916211 |        0.415392   | 0         | 0.0837893   | 0.916211 | 0        |            0.0837893 | PCA                     |                  0.916211 |                  0.415392   | True                            | True                        | Territory-matched / region-differing vs hippocampus (both PCA-leaning) |
| lingual_gyrus          | Lingual Gyrus                                                       |               26873 | PCA           |        0.866409 |        0.567206   | 0         | 0           | 0.866409 | 0.133591 |            0         | PCA                     |                  0.866409 |                  0.567206   | True                            | True                        | Additional PCA cortical ROI                                            |
| frontal_medial_cortex  | Frontal Medial Cortex                                               |                7766 | ACA           |        0.999742 |        0.00344205 | 0.999742  | 0.000257533 | 0        | 0        |            1         | ACA                     |                  0.999742 |                  0.00344205 | True                            | True                        | Region-matched cortex / ACA-leaning vs lateral frontal MCA             |
| middle_frontal_gyrus   | Middle Frontal Gyrus                                                |               41847 | MCA           |        0.984419 |        0.115849   | 0.0155806 | 0.984419    | 0        | 0        |            1         | MCA                     |                  0.984419 |                  0.115849   | True                            | True                        | MCA-leaning isocortex candidate (DLPFC-adjacent)                       |
| juxtapositional_lobule | Juxtapositional Lobule Cortex (formerly Supplementary Motor Cortex) |               11676 | ACA           |        1        |        0          | 1         | 0           | 0        | 0        |            1         | ACA                     |                  1        |                  0          | True                            | True                        | Medial frontal ACA-adjacent candidate                                  |

## Interpretation guidance

1. If hippocampus is PCA-dominant but superior frontal is mixed ACA/MCA with high entropy,
   Yang still supports a **posterior vs anterior circulation** contrast, but **not** a clean
   MCA territorial identity test. Use Winkler/Garcia regional metadata or additional atlases
   for MCA-pure cortex.
2. Prefer 2×2 contrasts from the candidate table where one ROI is PCA-dominant and another
   is MCA- or ACA-dominant with P≥0.50 and H≤1.0.
3. Replace hard-atlas fractions with BMM probabilities before locking the pre-registration
   specificity threshold.
