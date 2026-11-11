# Within-Garcia Class B path — ROSMAP vascular

## Critical finding

Garcia has **no clean ACA isocortex arm** on disk.
- Midtemporal + Angular -> MCA (TerritoryVoxels)
- Prefrontal (MFG proxy) -> MCA-leaning
- Hippocampus / entorhinal -> PCA-leaning (tissue class differs)

True ACA↔MCA Class B therefore requires **Yang SFG × Garcia MCA** 
(harmonized), not within-Garcia alone.

## TerritoryVoxels by region

```json
{
  "Midtemporal_cortex": {
    "tv_dominant": "MCA",
    "tv_p_ACA": 0.0,
    "tv_p_MCA": 1.0,
    "tv_p_PCA": 0.0,
    "tv_p_VB": 0.0,
    "tv_p_dominant": 1.0,
    "tv_assignment_fraction": 0.5096273845605277,
    "tv_source": "TerritoryVoxels_HO"
  },
  "Angular_gyrus": {
    "tv_dominant": "MCA",
    "tv_p_ACA": 0.0,
    "tv_p_MCA": 1.0,
    "tv_p_PCA": 0.0,
    "tv_p_VB": 0.0,
    "tv_p_dominant": 1.0,
    "tv_assignment_fraction": 0.5989151570975335,
    "tv_source": "TerritoryVoxels_HO"
  },
  "Prefrontal_cortex": {
    "tv_dominant": "MCA",
    "tv_p_ACA": 0.0,
    "tv_p_MCA": 1.0,
    "tv_p_PCA": 0.0,
    "tv_p_VB": 0.0,
    "tv_p_dominant": 1.0,
    "tv_assignment_fraction": 0.4400975845622981,
    "tv_source": "TerritoryVoxels_HO"
  },
  "Hippocampus": {
    "tv_dominant": "PCA",
    "tv_p_ACA": 0.0,
    "tv_p_MCA": 0.25,
    "tv_p_PCA": 0.75,
    "tv_p_VB": 0.0,
    "tv_p_dominant": 0.75,
    "tv_assignment_fraction": 0.31,
    "tv_source": "prior_Yang_hippocampus_TerritoryVoxels"
  },
  "Entorhinal_cortex": {
    "tv_dominant": "PCA",
    "tv_p_ACA": 0.0,
    "tv_p_MCA": 0.4681107099879663,
    "tv_p_PCA": 0.5318892900120337,
    "tv_p_VB": 0.0,
    "tv_p_dominant": 0.5318892900120337,
    "tv_assignment_fraction": 0.2654273668072058,
    "tv_source": "TerritoryVoxels_HO"
  },
  "Mammillary_body": {
    "tv_dominant": "unknown",
    "tv_p_ACA": NaN,
    "tv_p_MCA": NaN,
    "tv_p_PCA": NaN,
    "tv_p_VB": NaN,
    "tv_p_dominant": NaN,
    "tv_assignment_fraction": NaN,
    "tv_source": "unmapped"
  },
  "Anterior_thalamus": {
    "tv_dominant": "unknown",
    "tv_p_ACA": NaN,
    "tv_p_MCA": NaN,
    "tv_p_PCA": NaN,
    "tv_p_VB": NaN,
    "tv_p_dominant": NaN,
    "tv_assignment_fraction": NaN,
    "tv_source": "unmapped"
  }
}
```

## Arm × lineage counts

        class_b_arm vascular_lineage    n
      MCA_isocortex             endo  559
      MCA_isocortex              fib  225
      MCA_isocortex             peri  400
      MCA_isocortex              smc   79
PCA_medial_temporal             endo  751
PCA_medial_temporal              fib  927
PCA_medial_temporal             peri  394
PCA_medial_temporal              smc  122
    PFC_MCA_leaning             endo 2524
    PFC_MCA_leaning              fib 1632
    PFC_MCA_leaning             peri 2302
    PFC_MCA_leaning              smc  305
              other             endo  511
              other              fib  580
              other             peri  262
              other              smc  116

## Contrasts run

MCA isocortex vs PFC: 2559 endo cells — both MCA-leaning (region contrast, not ACA↔MCA).
MCA isocortex vs PCA medial temporal: 1310 endo cells — territory differs but tissue class also differs.

### MCA_isocortex_vs_PFC_MCA_leaning

     names    scores  logfoldchanges    pvals_adj                                 class_b_status
FP671120.1 10.912179        1.244784 3.380969e-23 NOT_true_ClassB_same_territory_region_contrast
FP236383.1 10.052409        1.410401 1.503315e-19 NOT_true_ClassB_same_territory_region_contrast
     PLCG2  8.302241        1.071136 8.566064e-13 NOT_true_ClassB_same_territory_region_contrast
    RIMKLB  4.946819        0.568361 3.152534e-03 NOT_true_ClassB_same_territory_region_contrast
    SHANK3  4.938887        0.614364 3.152534e-03 NOT_true_ClassB_same_territory_region_contrast
   RAPGEF2  4.666024        0.379667 9.362687e-03 NOT_true_ClassB_same_territory_region_contrast
AC105402.3  4.617720        0.968316 1.037204e-02 NOT_true_ClassB_same_territory_region_contrast
     CPNE8  4.597129        0.455825 1.037204e-02 NOT_true_ClassB_same_territory_region_contrast
      AKT3  4.594895        0.357310 1.037204e-02 NOT_true_ClassB_same_territory_region_contrast
    PICALM  4.554360        0.383733 1.174840e-02 NOT_true_ClassB_same_territory_region_contrast
     HERC2  4.486501        0.387674 1.517645e-02 NOT_true_ClassB_same_territory_region_contrast
      SMG1  4.384323        0.507834 2.295321e-02 NOT_true_ClassB_same_territory_region_contrast

### MCA_isocortex_vs_PCA_medial_temporal

   names   scores  logfoldchanges  pvals_adj                           class_b_status
 ST3GAL6 3.869657        0.592589   0.722290 exploratory_ClassA_like_tissue_differing
  PICALM 3.753891        0.282589   0.737263 exploratory_ClassA_like_tissue_differing
    SRGN 3.694015        0.428573   0.822579 exploratory_ClassA_like_tissue_differing
   THSD4 3.574262        0.343790   1.000000 exploratory_ClassA_like_tissue_differing
   TACC1 3.562523        0.314115   1.000000 exploratory_ClassA_like_tissue_differing
SLC39A10 3.301460        0.383442   1.000000 exploratory_ClassA_like_tissue_differing
   MECOM 3.189239        0.234260   1.000000 exploratory_ClassA_like_tissue_differing
  KCNIP4 3.181708        0.703646   1.000000 exploratory_ClassA_like_tissue_differing
 SEC14L1 3.137853        0.281248   1.000000 exploratory_ClassA_like_tissue_differing
   SYNE1 3.102267        0.281261   1.000000 exploratory_ClassA_like_tissue_differing
    SYT1 3.068157        0.725763   1.000000 exploratory_ClassA_like_tissue_differing
  MT-CO1 3.049921        0.352145   1.000000 exploratory_ClassA_like_tissue_differing

