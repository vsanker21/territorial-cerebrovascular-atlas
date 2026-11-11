# Yang vascular QC + Aim3 first-pass DE

- Input cells: 319,449
- After basic QC (genes≥200, counts≥300, MT≤10.0%): 309,402
- After vascular gate: **143,148**
- Object: `data/processed/yang/yang_vascular_qc.h5ad`

## Cells by stream × disease

analysis_stream disease     n
Aim3_borderzone      AD 16533
Aim3_borderzone control 39392
     H1_primary      AD 29668
     H1_primary control 57555

## Lineage × territory_contrast

vascular_lineage     territory_contrast     n
            endo ACA_MCA_borderzone_SFG 34435
            endo        PCA_hippocampus 54365
             fib ACA_MCA_borderzone_SFG  2662
             fib        PCA_hippocampus  4726
            peri ACA_MCA_borderzone_SFG 16492
            peri        PCA_hippocampus 19861
             smc ACA_MCA_borderzone_SFG  2336
             smc        PCA_hippocampus  8271

## Aim3 DE (borderzone SFG vs PCA hippocampus)

Positive scores = higher in ACA/MCA borderzone SFG relative to PCA hippocampus.

### control (top 15)

     names    scores  logfoldchanges     pvals_adj
     NEAT1 30.203161        0.574680 7.219526e-196
    LINGO1 28.570511        2.806369 2.620243e-175
FP671120.1 22.338930        0.763963 1.297126e-106
     PTH1R 19.907656        1.299763  2.342867e-84
     EGFL7 18.800737        0.900470  4.162848e-75
    ATP1A2 16.865749        1.069913  1.925492e-60
      DLC1 15.311142        0.755000  1.440417e-49
     AGAP1 13.856349        0.602742  2.055509e-40
      SHC2 13.682967        1.042630  2.049751e-39
    ZBTB46 13.543548        1.036101  1.261444e-38
    RIPOR3 13.481836        1.488409  2.795742e-38
     ITIH5 13.421058        0.731997  5.871319e-38
    COL5A3 12.917068        1.121519  4.319492e-35
     CTBP2 12.670524        0.663191  9.629604e-34
      HIP1 12.310199        0.713588  8.112811e-32

### AD (top 15)

   names    scores  logfoldchanges     pvals_adj
   PCDH9 36.219078        1.909088 1.022680e-282
 PITPNC1 26.941959        1.504595 1.188321e-155
   CADM2 20.724670        1.335150  1.392110e-91
   TACC1 19.418131        0.806890  3.031132e-80
  LINGO1 18.519400        2.226872  6.900864e-73
 PLEKHG1 18.317207        1.012679  2.529189e-71
IL1RAPL1 18.049246        1.080974  2.980224e-69
  AKAP12 18.025396        1.607412  4.129379e-69
   PDE4B 17.761202        1.270085  4.304360e-67
  ZSWIM6 17.218723        0.895712  5.371198e-63
  APOLD1 17.089939        1.109361  4.549698e-62
  ATP1B3 16.971935        2.035810  3.173969e-61
  ADGRB3 16.772501        1.187773  8.469802e-60
   EVA1C 16.770065        1.021927  8.469802e-60
  ARL17B 16.269300        1.773187  3.037205e-56

