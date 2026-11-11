# Harmonized Class B — tight endo gate v2b (Yang + Garcia)

Batches: Yang tight SFG | Garcia ex vivo MCA v2b | Garcia ROSMAP MCA v2b.
Hard neuronal contaminant caps (LINGO1/RBFOX1/KCNIP4); Garcia-calibrated score thresholds.

## Separability on Harmony embedding

| Metric | Value |
|---|---:|
| 5-fold CV AUROC | 0.682 +/- 0.037 |
| N Yang SFG | 1000 |
| N Garcia MCA | 319 |
| N cells / batches | 1319 / 3 |

## Top genes higher in Garcia MCA

    names    scores  logfoldchanges    pvals_adj
   MALAT1 19.143332        1.336469 3.686681e-77
      B2M 11.805666        2.535378 1.747000e-28
   ARGLU1 11.635586        2.542639 1.138300e-27
     MT2A 11.037309        2.649289 9.401869e-25
  SPARCL1 10.326012        1.935953 1.501130e-21
     DDX5  9.784878        2.054181 3.370802e-19
    IFI27  8.863323        2.685596 1.308338e-15
     PLP1  8.862732       30.195471 1.308338e-15
     RGCC  8.488050        3.047913 3.353540e-14
     EIF1  8.422380        2.769655 5.619833e-14
     FTH1  8.354433        3.719441 9.181105e-14
   TMSB4X  8.109653        2.931253 6.543544e-13
  RUNDC3B  8.006170        1.933060 1.468857e-12
    H3F3B  7.921172        2.360993 2.628410e-12
   TMSB10  7.851706        1.743808 4.437014e-12
HNRNPA2B1  7.734295        1.657950 1.058776e-11
    ITM2A  7.733958        2.449988 1.058776e-11
     SYT1  7.669724        3.602993 1.606754e-11
 RASGEF1B  7.669555        3.395756 1.606754e-11
     XAF1  7.662803        1.609378 1.647792e-11

## Higher in Yang SFG (bottom ranks)

   names   scores  logfoldchanges  pvals_adj
    RHOB 4.427737        1.628937   0.001238
 HNRNPH1 4.427567        0.956142   0.001238
     CA4 4.416257        4.049008   0.001299
    ICA1 4.394565        1.280808   0.001431
    CCNH 4.373800        0.926769   0.001568
  RNF213 4.372281        1.003191   0.001573
   SRSF4 4.362321        0.871089   0.001640
    EMP1 4.354218        1.525430   0.001695
   SERP1 4.348985        2.886106   0.001730
   TRA2B 4.346959        1.321778   0.001739
    CD46 4.336324        0.904432   0.001819
  PRKAG2 4.333454        1.400483   0.001836
   USP54 4.332188        0.860576   0.001839
MPHOSPH8 4.326364        1.433617   0.001882
    MT1M 4.322228        2.056873   0.001910

Object: `data/processed/class_b/yang_garcia_harmony_endo_tight.h5ad`
