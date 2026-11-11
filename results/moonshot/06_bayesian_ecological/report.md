# Moonshot 06 — Hierarchical Bayesian ecological vulnerability

**Honesty lock:** ECOLOGICAL_NOT_CONFIRMATORY_MOLECULAR_H3

## Mean territory shares
|          |       ACA |      MCA |       PCA |         VB |
|:---------|----------:|---------:|----------:|-----------:|
| ATLAS    | 0.0953958 | 0.691805 | 0.0912242 | 0.0968498  |
| ISLES22  | 0.0641878 | 0.479688 | 0.124389  | 0.234937   |
| ADNI_WMH | 0.0409487 | 0.131565 | 0.0229944 | 0.00451709 |

## Hierarchical models
| outcome        |   n_territories |   beta_mean |   beta_sd | beta_ci95                               |   P_beta_gt_0 |   spearman_raw |
|:---------------|----------------:|------------:|----------:|:----------------------------------------|--------------:|---------------:|
| share_ATLAS    |               4 |  -0.421834  |   2.83312 | [-5.668271999288554, 6.04444287705069]  |       0.40875 |      -0.632456 |
| share_ISLES22  |               4 |  -0.0708205 |   2.81491 | [-5.323060852423531, 6.212614336312661] |       0.47125 |      -0.948683 |
| share_ADNI_WMH |               4 |   0.587723  |   2.81624 | [-4.705546178632169, 6.837742727187886] |       0.58875 |      -0.316228 |

Ecological hierarchical posteriors estimated; do NOT interpret as confirmatory molecular H3. Use PPC and beta CI for honest uncertainty at n_territories≈3–4.