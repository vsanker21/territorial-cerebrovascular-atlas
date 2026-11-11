# H1 first-pass variance partition — Yang tight endo controls

**Critical caveat:** In Yang, `territory_contrast` is aliased with anatomical region
(PCA_hippocampus vs ACA/MCA borderzone SFG). This estimates territorial/regional
contrast after donor, **not** territory independent of region. Proper H1 requires
Class B (region-matched, territory-differing) and/or multi-region cohorts.

| Metric | Value |
|---|---:|
| Cells | 12000 |
| Donors | 8 |
| HVGs / PCs | 1000 / 30 |
| Mean PC R2 donor | 0.0567 |
| Mean PC R2 donor+territory | 0.0613 |
| Delta R2 territory after donor | 0.0046 |
| HVGs territory Bonferroni | 120 (12.0%) |

                      model  mean_pc_r2
                      donor    0.056725
   donor+territory_contrast    0.061332
               donor+region    0.061332
delta_territory_after_donor    0.004607
   delta_region_after_donor    0.004607
