# Within-ROSMAP H1 path — published Endo subtypes + neuronal caps

Garcia ROSMAP is the only **within-dataset** multi-region vascular object on disk.
No ACA cortex arm; contrasts are:
1. MCA isocortex (Midtemporal+Angular) vs PFC (both MCA-leaning) — **region null**
2. MCA isocortex vs PCA medial temporal — territory differs but tissue class also differs

Neuronal caps: LINGO1/RBFOX1/KCNIP4/SNAP25 ≤ 0.5; require canonical endo marker.
capEndo stratum conditions on A–V zonation (H1 confound control).

## Counts

- Endo subtypes raw: **4,345**
- After neuronal caps: **928**
- By arm: `{'PFC_MCA_leaning': 497, 'PCA_medial_temporal': 210, 'MCA_isocortex': 123, 'other': 98}`
- By zonation: `{'capEndo': 648, 'vEndo': 206, 'aEndo': 74}`

## Variance partition (arm after subject)

- **all_endo**: delta R2 arm after subject = **0.0051** (subject R2=0.3304; n=830)
- **capEndo**: delta R2 arm after subject = **0.0047** (subject R2=0.4070; n=587)

## Top DE genes

### MCA_isocortex_vs_PFC_same_territory_region_null | all_endo_subtypes

     names  logfoldchanges  pvals_adj
FP671120.1        1.496138   0.000241
     PLCG2        1.053607   1.000000
   ANKRD17        0.837506   1.000000
      NKTR        0.704821   1.000000
     PDE4D        0.868062   1.000000
    PICALM        0.604284   1.000000
     TEAD1        0.751547   1.000000
    CCSER2        0.707962   1.000000
    MYCBP2        0.761049   1.000000
    PRRC2C        0.704510   1.000000

### MCA_isocortex_vs_PFC_same_territory_region_null | capEndo_only

     names  logfoldchanges  pvals_adj
FP671120.1        1.403457   0.240318
    ELOVL7        0.741743   0.406418
    RIMKLB        1.045727   1.000000
   ANKRD17        1.008557   1.000000
    PPP3CC        1.019054   1.000000
    COBLL1        0.720931   1.000000
     LAMA5        1.225387   1.000000
    CCSER2        0.845226   1.000000
    PICALM        0.639233   1.000000
    TBC1D5        0.867777   1.000000

### MCA_isocortex_vs_PCA_medial_temporal | all_endo_subtypes

     names  logfoldchanges  pvals_adj
    PICALM        0.444247        1.0
FP671120.1        0.739090        1.0
      OSMR        0.734256        1.0
     PDE4D        0.630868        1.0
   ST3GAL6        0.774285        1.0
      SRGN        0.554607        1.0
     ANXA3        0.679908        1.0
    MT-ND1        0.669132        1.0
     PDZD2        0.535764        1.0
      MT2A        0.687836        1.0

### MCA_isocortex_vs_PCA_medial_temporal | capEndo_only

  names  logfoldchanges  pvals_adj
   SRGN        0.805026        1.0
 PICALM        0.515211        1.0
   OSMR        0.929736        1.0
 ELOVL7        0.428700        1.0
  ANXA3        0.826680        1.0
  AGFG1        0.650727        1.0
 SORBS2        0.562590        1.0
ST3GAL6        0.935154        1.0
TSC22D1        0.603484        1.0
  PDE4D        0.802685        1.0

