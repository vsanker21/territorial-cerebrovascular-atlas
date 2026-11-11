# Moonshot 04 — CoW topology IV + conformal

**Honesty lock:** Only nested OOS / conformal metrics count; in-sample R² is descriptive

- n QC-pass: 78
- Track C topology support: **FAIL_OOS**

| outcome             |   oos_r2_cow_features |   oos_r2_demo_only |   delta_oos_r2 |   oos_spearman_cow |   oos_spearman_p |   oos_spearman_demo |   conformal_q90 |   conformal_coverage_90 |   n |
|:--------------------|----------------------:|-------------------:|---------------:|-------------------:|-----------------:|--------------------:|----------------:|------------------------:|----:|
| lesion_share_MCA    |            -0.0509053 |         -0.0372271 |   -0.0136783   |          -0.280138 |        0.0129872 |          -0.171723  |       0.227851  |                0.764706 |  78 |
| lesion_share_ACA    |            -0.0439538 |         -0.0437787 |   -0.000175143 |          -0.126516 |        0.269697  |          -0.0644497 |       0.0877671 |                0.823529 |  78 |
| mca_minus_aca_share |            -0.0511043 |         -0.0385973 |   -0.0125069   |          -0.222474 |        0.0502626 |          -0.154752  |       0.299658  |                0.764706 |  78 |