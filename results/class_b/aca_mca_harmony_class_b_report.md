# Class B — true ACA cortex vs Garcia MCA (Harmony)

First TerritoryVoxels-mapped **ACA isocortex vascular** arm vs **MCA** endo.

## Arms

| Arm | Sources | N endo (after caps) | N donors |
|---|---|---:|---:|
| **ACA** | Allen ACC + LIBD dACC GLOBUS (paper Azimuth) + Linnarsson A24/A32/A25/A23 | 264 | 18 |
| **MCA** | Garcia temporal + ROSMAP MCA tight v2b | 319 | 36 |

Batch keys: `{'Garcia_exvivo': 207, 'Allen_ACC': 165, 'Garcia_ROSMAP': 112, 'LIBD_dACC_globus': 87, 'Linnarsson': 12}`

Neuronal caps: published-annotation ACA arms keep author Endo labels (no harsh neuronal discard); MCA uses endo-marker + neuronal gate. Harmony batch = source dataset.

## Separability

| Metric | Value |
|---|---:|
| Cell 5-fold CV AUROC | 0.826 ± 0.039 |
| **Donor-level OOF AUROC** | **0.787** |
| Donors scored (ACA / MCA) | 18 / 36 |
| Batch LOBO AUROC (multi-class only) | NA (each batch is single-arm) |

**Verdict:** `SUPPORTIVE_transferable_ClassB_signal`

## Caveats

- Still **cross-dataset** (tech + processing confound); Harmony reduces but does not erase this.
- ACA arm is **not vessel-enriched**; MCA is tight/vessel-enriched mix — enrichment asymmetry remains.
- Batch LOBO undefined because no batch contains both arms.
- Top MCA genes include **PLP1** (possible oligo ambient) — interpret DE cautiously; donor OOF is the primary claim.
- Yang SFG is **not** used as ACA (borderzone).

## Neuronal contaminant check in DE top-200

- none in top 200

## FDR hits (MCA vs ACA Wilcoxon, padj < 0.05)

13 genes; top: MACF1, RASGEF1B, NPIPB5, RPGR, PLP1, CP, RBMS3-AS3, FKBP5, IL4R, EMP1, IGFBP3, SLCO4A1, LINC00623.

    names   scores  logfoldchanges    pvals_adj
    MACF1 7.794848        1.930411 4.360073e-13
 RASGEF1B 5.297661        1.784101 2.577006e-06
   NPIPB5 4.202564        1.271909 3.666834e-04
     RPGR 3.758993        0.365586 1.895735e-03
     PLP1 3.726392        1.397508 2.124030e-03
       CP 3.310483        1.217317 8.291064e-03
RBMS3-AS3 3.280104        1.945530 9.127281e-03
    FKBP5 3.241329        0.246236 1.029533e-02
     IL4R 2.817269       -0.008454 3.369667e-02
     EMP1 2.808625        0.357056 3.440288e-02
   IGFBP3 2.733050        1.015024 4.185277e-02
  SLCO4A1 2.693287       -0.178094 4.610006e-02
LINC00623 2.668342        1.644934 4.896969e-02
     FOSB 2.512252        2.323673 7.098722e-02
     CCNH 2.462610        0.037926 7.951156e-02

## Higher in ACA (bottom ranks)

      names   scores  logfoldchanges  pvals_adj
  LINC00320 0.547055        1.292502        1.0
    TMEM217 0.542115        0.988755        1.0
     ADAM28 0.539398        0.124431        1.0
MIR4435-2HG 0.535941       -0.456919        1.0
    TMEM156 0.535447       -0.249851        1.0
    FAM157C 0.532483        0.386578        1.0
       GNLY 0.521616       25.300734        1.0
        IL6 0.521616       26.178297        1.0
      SRSF4 0.520381       -0.758858        1.0
   CNTNAP3B 0.518405       -1.114234        1.0

## Objects

- `data/processed/aca_cohort_hunt/aca_endo_pool.h5ad`
- `data/processed/class_b/aca_mca_harmony_endo.h5ad`
- DE: `results/class_b/aca_mca_harmony_class_b_de.csv`
- Donor folds: `results/class_b/aca_mca_harmony_donor_aware.csv`
