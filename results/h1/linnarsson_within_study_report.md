# Linnarsson / Siletti within-study vascular (F/T acquired)

**Verdict:** `UNDERPOWERED_n_donors_3__slices_acquired`

## Design
- Same study: Siletti et al. HBCA CELLxGENE collection `283d65eb-…`
- **ACA:** A24, A32, A25, A23 (A14 excluded)
- **Lateral frontal (MCA):** MFG A46 + IFG A44–45
- **Temporal (MCA):** MTG + STG
- Metric: leave-one-**donor**-out **pooled cell-level** AUROC on Endo PCA
- Note: the same 3 adult donors contribute cells to ACA, lateral frontal, and temporal — cross-cohort donor-mean OOF is invalid here.

## Counts
- Vascular object: `data\processed\aca_cohort_hunt\linnarsson_ft_vascular.h5ad` — n=657
- Endo n=234; donors=3
- Endo by arm: `{'temporal': 157, 'ACA': 43, 'lateral_frontal': 34}`
- Endo by dissection: `{'MTG': 131, 'IFG': 28, 'STG': 26, 'A24': 17, 'A25': 11, 'A23': 9, 'A32': 6, 'MFG_A46': 6}`

## TerritoryVoxels priors (dominant)
- A24: dominant=ACA pACA=1.00 pMCA=0.00 assign=0.72
- A32: dominant=ACA pACA=1.00 pMCA=0.00 assign=0.56
- A25: dominant=ACA pACA=1.00 pMCA=0.00 assign=0.53
- A23: dominant=ACA pACA=0.87 pMCA=0.00 assign=0.44
- MFG_A46: dominant=MCA pACA=0.00 pMCA=1.00 assign=0.44
- IFG: dominant=MCA pACA=0.00 pMCA=1.00 assign=0.46
- MTG: dominant=MCA pACA=0.00 pMCA=1.00 assign=0.51
- STG: dominant=MCA pACA=0.00 pMCA=1.00 assign=0.49

## Contrasts (LODO pooled cell AUROC)
| contrast | n_pos | n_neg | donors | pooled AUROC | mean fold AUROC | ambient pos/neg |
|---|---:|---:|---:|---:|---:|---|
| ACA_vs_lateral_frontal | 43 | 34 | 3+3 | 0.496 | 0.557 | 0.555/0.589 |
| ACA_vs_temporal | 43 | 157 | 3+3 | 0.423 | 0.451 | 0.555/0.280 |
| lateral_frontal_vs_temporal_MCA_null | 34 | 157 | 3+3 | 0.310 | 0.293 | 0.589/0.280 |

## Interpretation locks
- Do **not** treat lateral frontal as ACA (TerritoryVoxels MCA-dominant).
- F vs T = MCA vs MCA null (regional, not territorial).
- Only 3 adult donors in the Siletti census → underpowered as sole H1 proof even if AUROC is high.
- Class B still requires external ACA arms; enrichment-matched DE needs vessel-enriched ACC (no public download).
