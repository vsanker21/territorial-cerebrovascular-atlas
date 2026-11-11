# Wälchli adult temporal + pathology — formal exploratory A–V tables

**Verdict:** `EXPLORATORY_A-V_AND_PATHOLOGY__NOT_TERRITORY`

## Methods locks
- GEO GSE256493 Seurat objects are **double-gzipped**; decompressed to `data/raw/walchli/geo_GSE256493/suppl/_decomp/`.
- Adult control = **temporal lobe** FACS CD31+/CD45− EC → **MCA anatomy**, not ACA.
- Use for A–V zonation / pathology EC-state precedent only; **not** Class B territory proof.
- Enrichment yes relative to unsorted ACC, but wrong territory for ACA marker claims.

## Cohort table

```
                                     cohort                      file n_cells
adult_temporal_sorted adult_temporal_sorted adult_temporal_sorted.rds   76125
avm_sorted                       avm_sorted            avm_sorted.rds   20305
lgg_sorted                       lgg_sorted            lgg_sorted.rds   17373
                      n_genes path_values cluster_col n_ec_clusters
adult_temporal_sorted   26666          TL  ECclusters            12
avm_sorted               3000         AVM  ECclusters            12
lgg_sorted              26666         LGG  ECclusters            13
                      n_donors_proxy donor_col
adult_temporal_sorted              9   Patient
avm_sorted                         5   Patient
lgg_sorted                         6   Patient
```

## EC cluster counts

```
                                          cluster     n                cohort
adult_temporal_sorted.1                 Capillary 36268 adult_temporal_sorted
adult_temporal_sorted.2                 Arteriole 18487 adult_temporal_sorted
adult_temporal_sorted.3      Angiogenic capillary  4863 adult_temporal_sorted
adult_temporal_sorted.4                    Venule  4780 adult_temporal_sorted
adult_temporal_sorted.5                    Artery  3114 adult_temporal_sorted
adult_temporal_sorted.6              Large artery  2157 adult_temporal_sorted
adult_temporal_sorted.7                    EndoMT  2088 adult_temporal_sorted
adult_temporal_sorted.8                      Vein  1792 adult_temporal_sorted
adult_temporal_sorted.9                Large vein  1406 adult_temporal_sorted
adult_temporal_sorted.10            Mitochondrial   820 adult_temporal_sorted
adult_temporal_sorted.11       Proliferating cell   305 adult_temporal_sorted
adult_temporal_sorted.12               Stem-to-EC    45 adult_temporal_sorted
avm_sorted.1                            Capillary  5436            avm_sorted
avm_sorted.2                                 Vein  3440            avm_sorted
avm_sorted.3                 Angiogenic capillary  2275            avm_sorted
avm_sorted.4                               Venule  2160            avm_sorted
avm_sorted.5                           Large vein  1997            avm_sorted
avm_sorted.6                               Artery  1865            avm_sorted
avm_sorted.7                         Large artery  1129            avm_sorted
avm_sorted.8                            Arteriole  1071            avm_sorted
avm_sorted.9                        Mitochondrial   587            avm_sorted
avm_sorted.10                              EndoMT   231            avm_sorted
avm_sorted.11                  Proliferating cell    82            avm_sorted
avm_sorted.12                          Stem-to-EC    32            avm_sorted
lgg_sorted.1                            Capillary  6944            lgg_sorted
lgg_sorted.2                 Angiogenic capillary  4720            lgg_sorted
lgg_sorted.3                            Arteriole  2000            lgg_sorted
lgg_sorted.4                               Venule  1027            lgg_sorted
lgg_sorted.5                               Artery   854            lgg_sorted
lgg_sorted.6                           Large vein   513            lgg_sorted
lgg_sorted.7                                 Vein   418            lgg_sorted
lgg_sorted.8                        Mitochondrial   357            lgg_sorted
lgg_sorted.9                               EndoMT   306            lgg_sorted
lgg_sorted.10                        Large artery   157            lgg_sorted
lgg_sorted.11                  Proliferating cell    43            lgg_sorted
lgg_sorted.12                          Stem-to-EC    33            lgg_sorted
lgg_sorted.13            Proliferating Stem-to-EC     1            lgg_sorted
                                 frac
adult_temporal_sorted.1  4.764269e-01
adult_temporal_sorted.2  2.428506e-01
adult_temporal_sorted.3  6.388177e-02
adult_temporal_sorted.4  6.279146e-02
adult_temporal_sorted.5  4.090640e-02
adult_temporal_sorted.6  2.833498e-02
adult_temporal_sorted.7  2.742857e-02
adult_temporal_sorted.8  2.354023e-02
adult_temporal_sorted.9  1.846962e-02
adult_temporal_sorted.10 1.077176e-02
adult_temporal_sorted.11 4.006568e-03
adult_temporal_sorted.12 5.911330e-04
avm_sorted.1             2.677173e-01
avm_sorted.2             1.694164e-01
avm_sorted.3             1.120414e-01
avm_sorted.4             1.063777e-01
avm_sorted.5             9.835016e-02
avm_sorted.6             9.184930e-02
avm_sorted.7             5.560207e-02
avm_sorted.8             5.274563e-02
avm_sorted.9             2.890914e-02
avm_sorted.10            1.137651e-02
avm_sorted.11            4.038414e-03
avm_sorted.12            1.575967e-03
lgg_sorted.1             3.997007e-01
lgg_sorted.2             2.716859e-01
lgg_sorted.3             1.151212e-01
lgg_sorted.4             5.911472e-02
lgg_sorted.5             4.915674e-02
lgg_sorted.6             2.952858e-02
lgg_sorted.7             2.406032e-02
lgg_sorted.8             2.054913e-02
lgg_sorted.9             1.761354e-02
lgg_sorted.10            9.037011e-03
lgg_sorted.11            2.475105e-03
lgg_sorted.12            1.899499e-03
lgg_sorted.13            5.756058e-05
```

## Coarse A–V bins

```
                                      av_bin     n                cohort
adult_temporal_sorted.1            capillary 36268 adult_temporal_sorted
adult_temporal_sorted.2             arterial 23758 adult_temporal_sorted
adult_temporal_sorted.3               venous  7978 adult_temporal_sorted
adult_temporal_sorted.4 angiogenic_capillary  4863 adult_temporal_sorted
adult_temporal_sorted.5      activated_other  3258 adult_temporal_sorted
avm_sorted.1                          venous  7597            avm_sorted
avm_sorted.2                       capillary  5436            avm_sorted
avm_sorted.3                        arterial  4065            avm_sorted
avm_sorted.4            angiogenic_capillary  2275            avm_sorted
avm_sorted.5                 activated_other   932            avm_sorted
lgg_sorted.1                       capillary  6944            lgg_sorted
lgg_sorted.2            angiogenic_capillary  4720            lgg_sorted
lgg_sorted.3                        arterial  3011            lgg_sorted
lgg_sorted.4                          venous  1958            lgg_sorted
lgg_sorted.5                 activated_other   740            lgg_sorted
                              frac
adult_temporal_sorted.1 0.47642693
adult_temporal_sorted.2 0.31209195
adult_temporal_sorted.3 0.10480131
adult_temporal_sorted.4 0.06388177
adult_temporal_sorted.5 0.04279803
avm_sorted.1            0.37414430
avm_sorted.2            0.26771731
avm_sorted.3            0.20019700
avm_sorted.4            0.11204137
avm_sorted.5            0.04590002
lgg_sorted.1            0.39970068
lgg_sorted.2            0.27168595
lgg_sorted.3            0.17331491
lgg_sorted.4            0.11270362
lgg_sorted.5            0.04259483
```

## What this enables / does not
- Enables: adult temporal control A–V composition vs AVM / LGG (exploratory H4-adjacent / Aim2 zonation).
- Does **not** unlock enrichment-matched ACA↔MCA DE or territorial H1.
