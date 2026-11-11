# H4 within-donor CV — Yang tight endo

For each donor with both arms, 5-fold stratified CV predicts `territory_contrast` from HVGs.
This tests whether territorial/regional arms are readable **within** a donor.

| Metric | Value |
|---|---:|
| Control within-donor AUROC mean | 0.953 |
| AD within-donor AUROC mean | 0.981 |
| Delta (AD - control) | **+0.028** |
| N control / AD donors | 4 / 4 |

## Per donor

| disease | donor_id | n_cells | auroc_mean | auroc_std | n_SFG | n_PCA |
|---|---|---:|---:|---:|---:|---:|
| control | Control_1 | 3000 | 0.828 | 0.009 | 4087 | 4234 |
| control | Control_2 | 2437 | 0.989 | 0.005 | 2257 | 937 |
| control | Control_3 | 2734 | 0.996 | 0.001 | 1234 | 1869 |
| control | Control_4 | 3000 | 0.998 | 0.001 | 2231 | 2769 |
| AD | AD_1 | 2338 | 0.989 | 0.006 | 838 | 3582 |
| AD | AD_2 | 1173 | 0.986 | 0.004 | 873 | 300 |
| AD | AD_3 | 766 | 0.956 | 0.018 | 422 | 344 |
| AD | AD_4 | 2720 | 0.993 | 0.002 | 1220 | 3474 |

## Interpretation

1. Within-donor separability is **strong in both control and AD** (AUROC ~0.95–0.98).
2. Delta is **not negative** — this pass does **not** support H4 degradation on Yang alone.
3. Caveat: in Yang, territory_contrast is aliased with anatomical region (hippocampus vs SFG),
   so high AUROC may reflect regional identity rather than arterial-territory identity.
4. Cross-donor LODO (~chance in controls) shows the program does **not** transfer across donors;
   see `yang_h4_lodo_report.md`.

H4 remains open pending Class B / multi-region designs and larger donor n.
