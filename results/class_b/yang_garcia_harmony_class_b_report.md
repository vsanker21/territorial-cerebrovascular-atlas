# Harmonized Class B — Yang SFG borderzone vs Garcia MCA (Harmony)

Batches: Yang | Garcia ex vivo temporal | Garcia ROSMAP Midtemporal+Angular.
Endothelial only; Yang controls; Harmony on PCA.

## Separability on Harmony embedding

| Metric | Value |
|---|---:|
| 5-fold CV AUROC | 0.848 +/- 0.007 |
| N Yang SFG | 4000 |
| N Garcia MCA | 4559 |
| N cells / batches | 8559 / 3 |

High AUROC after Harmony = residual territorial/regional signal beyond batch.
Yang arm is ACA/MCA **borderzone** (not pure ACA) — Aim3/Class B hybrid.

## Top genes higher in Garcia MCA

   names    scores  logfoldchanges  pvals_adj
  KCNIP4 56.856842        6.160317        0.0
   CSMD1 51.521233        5.391779        0.0
  LINGO1 51.466785       31.163971        0.0
   PTPRD 50.132370        4.720232        0.0
  MALAT1 49.372131        0.837321        0.0
   PCDH9 49.365868        3.763452        0.0
  RBFOX1 49.332211        5.052557        0.0
   CADM2 48.700191        4.295936        0.0
RASGEF1B 46.913506        4.456766        0.0
   NRXN3 44.194412        4.073610        0.0
   LRP1B 44.150024        3.825002        0.0
  ADGRB3 41.895729        3.314057        0.0
 CNTNAP2 41.320900        4.213007        0.0
   NRXN1 41.015484        3.947368        0.0
  LRRTM4 41.008938        5.369564        0.0
  FRMD4A 40.216846        3.707581        0.0
IL1RAPL1 39.754230        3.766371        0.0
    DLG2 39.584908        2.609183        0.0
  NKAIN2 39.475853        4.695795        0.0
   NLGN1 38.972710        4.995085        0.0

## Higher in Yang SFG (bottom ranks)

  names    scores  logfoldchanges     pvals_adj
  KCNQ5 22.296158        4.482242 6.247074e-108
 LRRTM3 22.267439        6.301658 1.180691e-107
  CNTN5 22.260015        4.253316 1.386904e-107
  H3F3B 22.252840        1.959051 1.620184e-107
 SPTLC2 22.183939        1.876537 7.477903e-107
HNRNPDL 22.174807        1.576034 9.077932e-107
  CRLF3 22.069199        2.056731 9.391230e-106
   SYN3 22.060913        4.716933 1.122897e-105
   THRB 22.036123        2.879604 1.933204e-105
  PSEN1 22.014589        1.612172 3.095669e-105
  ENOX1 21.972733        2.922172 7.718996e-105
  GRIA4 21.862459        4.488698 8.659306e-104
   EIF1 21.852341        2.009295 1.076060e-103
  TCF25 21.816999        1.434826 2.311647e-103
  CDH12 21.747377        4.382061 1.052052e-102

Object: `data/processed/class_b/yang_garcia_harmony_endo.h5ad`
