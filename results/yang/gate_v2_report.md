# Yang gate v2 (tight) report

## Gate rules

- Endo: any of ['CLDN5', 'PECAM1', 'FLT1', 'VWF', 'CDH5'] > 0; score_endo ≥ 0.5; margin vs other vascular ≥ 0.4
- Peri: any of ['PDGFRB', 'RGS5', 'KCNJ8', 'ABCC9'] > 0; score_peri ≥ 0.25; margin ≥ 0.4
- Shared: parenchyma ≤ 0.0; neuron ≤ 0.05; astro ≤ 0.25; oligo ≤ 0.0; micro ≤ 0.05
- Contaminant expression caps: {'LINGO1': 0.5, 'SNAP25': 0.5, 'SLC17A7': 0.5, 'RBFOX3': 0.5, 'GAD1': 0.5, 'MBP': 1.0, 'MOBP': 1.0, 'PLP1': 1.0, 'AQP4': 1.0, 'GFAP': 1.5, 'PTPRC': 0.5}

## Cell counts

- Soft vascular gate: **143,148** → `{'endo': 88800, 'peri': 36353, 'smc': 10607, 'fib': 7388}`
- Tight gate v2: **65,382** → `{'endo': 44059, 'peri': 15121, 'smc': 4173, 'fib': 2029}`
- Retention: **45.7%** of soft-gated cells
- Object: `data/processed/yang/yang_vascular_qc_tight.h5ad`

## Lineage × contrast × disease (tight)

disease                                    AD  control
vascular_lineage territory_contrast                   
endo             ACA_MCA_borderzone_SFG  3353     9809
                 PCA_hippocampus         9562    21335
fib              ACA_MCA_borderzone_SFG   191      368
                 PCA_hippocampus          685      785
peri             ACA_MCA_borderzone_SFG  1561     4483
                 PCA_hippocampus         2909     6168
smc              ACA_MCA_borderzone_SFG   170      394
                 PCA_hippocampus         1439     2170

## Aim3 top hits (tight endo/peri)

### endo / control

     names    scores  logfoldchanges     pvals_adj
     NEAT1 26.956219        0.565570 5.391673e-156
FP671120.1 21.261410        0.792779  1.083827e-96
     EGFL7 20.302998        0.827497  4.508054e-88
FP236383.1 15.109479        0.512688  2.351759e-48
     PTH1R 12.083827        1.656171  1.002960e-30
     BTNL9 11.502553        0.702354  8.765716e-28
      SHC2 11.397617        0.874359  2.882374e-27
   RAPGEF1 11.283164        0.586165  9.690997e-27
      INSR 11.058380        0.559259  1.195817e-25
      DLC1 10.839881        0.635076  1.288095e-24

### endo / AD

  names    scores  logfoldchanges     pvals_adj
  PCDH9 30.058163        2.472634 5.727962e-194
 SORBS1 24.516682        1.413226 1.644686e-128
PLEKHG1 23.337561        1.336676 2.060870e-116
 PLXNA2 22.459352        1.671627 8.693615e-108
 INPP5D 22.323122        1.822727 1.478066e-106
PITPNC1 19.050711        1.495800  3.105439e-77
 AKAP12 18.051979        1.995281  3.190939e-69
  TACC1 17.969759        0.938557  1.252706e-68
SEC14L1 17.941866        1.124900  1.863272e-68
 APOLD1 17.703226        1.372944  1.207263e-66

### peri / control

     names    scores  logfoldchanges    pvals_adj
    RIPOR3 21.014698        1.436772 8.070693e-94
FP671120.1 20.851372        1.044388 1.655829e-92
     PTH1R 19.279827        0.897948 4.434833e-79
     ITIH5 17.562994        0.923927 1.763070e-65
    ATP1A2 17.440399        0.916868 1.365747e-64
     EGFL7 16.238947        1.997737 5.980344e-56
     ZFHX3 14.415057        0.580633 6.343207e-44
      HIP1 13.507624        0.883196 1.970348e-38
   DENND2A 13.459188        0.884915 3.647409e-38
   SLC20A2 12.799370        0.829744 1.979782e-34

### peri / AD

  names    scores  logfoldchanges    pvals_adj
 ATP1B3 16.474993        3.179218 1.861207e-56
  PCDH9 15.976871        1.752571 3.105854e-53
PITPNC1 15.516176        1.855858 3.014481e-50
  PCBP3 12.487050        1.431803 4.910560e-32
ARHGEF7 11.779386        1.311825 2.388622e-28
   ELL2 10.780917        1.782493 1.291638e-23
 PHLPP1 10.660015        1.130152 4.375913e-23
  CADM2  9.690709        1.268318 7.910876e-19
  PDE4D  9.664343        1.172113 9.555077e-19
   TJP1  9.572894        0.799996 2.179005e-18

## H4 (tight)

- Control CV AUROC: **0.878**
- AD transfer AUROC: **0.695**
- Delta: **-0.182**

