# Garcia sample audit (MOESM2) — Class B readiness

## Verdict

Garcia supports Class B MCA arm via (i) ex vivo temporal epilepsy cortex and (ii) ROSMAP mid-temporal / angular gyrus if present in region counts. Pair with Yang superior frontal (ACA) for territory-differing isocortex contrast. Winkler normal cortex remains pending lobe labels but is likely temporal-enriched by clinical context.

## Ex vivo tissue locations

```json
{
  "Left temporal tip": 2,
  "Left temporal lobe": 2,
  "Left lateral temporal lobe": 2,
  "Right temporal tip": 1,
  "Left anterior temporal lobe, left medial temporal lobe": 1,
  "Right lateral temporal lobe,  medial temporal lobe": 1,
  "Posterior temporal lobe": 1,
  "Left anterior temporal lobe": 1,
  "Left Anterior Temporal Tip": 1,
  "Left Temporal Fissure": 1,
  "Anterior Lateral Temporal Lobe": 1,
  "Left Temporal Lobe": 1,
  "Left Temporal Tip": 1,
  "Repeat on left medial temporal lobe with hippocampus": 1
}
```

## ROSMAP brain_region counts

| brain_region      |   n_libraries |
|:------------------|--------------:|
| Prefrontal_cortex |           213 |
| Hippocampus       |             6 |

## Territory scores

| roi               | labels                             |   n_voxels_in_atlas | dominant_l2   |   p_dominant_l2 |   entropy_l2_bits |     p_ACA |    p_MCA |    p_PCA |       p_VB |   p_anterior_ACA_MCA | dominant_l1_majorized   |   p_dominant_l1_majorized |   entropy_l1_majorized_bits | high_specificity_pdom_ge_0.50   | high_specificity_H_le_1.0   | source      | label                                    |   n_libraries | class_b     | role                                                                         | tissue_locations                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
|:------------------|:-----------------------------------|--------------------:|:--------------|----------------:|------------------:|----------:|---------:|---------:|-----------:|---------------------:|:------------------------|--------------------------:|----------------------------:|:--------------------------------|:----------------------------|:------------|:-----------------------------------------|--------------:|:------------|:-----------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Prefrontal_cortex | Middle Frontal Gyrus               |               41847 | MCA           |        0.984419 |          0.115849 | 0.0155806 | 0.984419 | 0        | 0          |             1        | MCA                     |                  0.984419 |                    0.115849 | True                            | True                        | ROSMAP      | Prefrontal_cortex                        |           213 | conditional | mixed_or_MCA_leaning_proxy — PFC not medial-split; MFG used as lateral proxy | nan                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Hippocampus       | Left Hippocampus;Right Hippocampus |               11007 | PCA           |        0.757427 |          0.841281 | 0         | 0.236304 | 0.757427 | 0.00626874 |             0.236304 | PCA                     |                  0.757427 |                    0.841281 | True                            | True                        | ROSMAP      | Hippocampus                              |             6 | no          | PCA Class A / Yang replication                                               | nan                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| exvivo_temporal   | Temporal Pole                      |               34955 | MCA           |        1        |          0        | 0         | 1        | 0        | 0          |             1        | MCA                     |                  1        |                    0        | True                            | True                        | ex_vivo_BVE | temporal tip/lobe (from Tissue location) |            17 | yes_mca     | Garcia ex vivo epilepsy cortex — MCA prior                                   | {'Left temporal tip': 2, 'Left temporal lobe': 2, 'Left lateral temporal lobe': 2, 'Right temporal tip': 1, 'Left anterior temporal lobe, left medial temporal lobe': 1, 'Right lateral temporal lobe,  medial temporal lobe': 1, 'Posterior temporal lobe': 1, 'Left anterior temporal lobe': 1, 'Left Anterior Temporal Tip': 1, 'Left Temporal Fissure': 1, 'Anterior Lateral Temporal Lobe': 1, 'Left Temporal Lobe': 1, 'Left Temporal Tip': 1, 'Repeat on left medial temporal lobe with hippocampus': 1} |
