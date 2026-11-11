# Class B sensitivity — ACA vs MCA

**Meta verdict:** `ROBUST_SUPPORTIVE_separability__DE_contaminated`

Baseline donor OOF AUROC **0.787**; cell CV **0.826**.

## Donor OOF across sensitivities

| Analysis | Donor OOF | Cell CV | n ACA | n MCA |
|---|---:|---:|---:|---:|
| baseline | 0.787 | 0.826 | 264 | 319 |
| drop_Allen_ACC | 0.746 | 0.680 | 99 | 319 |
| drop_LIBD_dACC_globus | 0.667 | 0.911 | 177 | 319 |
| drop_Linnarsson | 0.893 | 0.838 | 252 | 319 |
| only_Allen_ACC | 1.000 | 0.944 | 165 | 319 |
| only_LIBD_dACC_globus | 0.723 | 0.672 | 87 | 319 |
| oligo_ambient_filter | 0.846 | 0.796 | 209 | 191 |
| balanced_subsample | 0.800 | 0.811 | 264 | 264 |

## Donor pseudobulk DE

- Genes tested: 3000
- FDR < 0.05: **1372** (inflated — cross-dataset)
- FDR < 0.10: 1831
- **Neuronal ambient hits (higher in ACA):** RBFOX1, KCNIP4, SNAP25, LINGO1 (all FDR≪0.05)

## Interpretation

- **Separability (donor OOF) is robust:** 7/8 sensitivities ≥0.67; median OOF **0.79**; oligo filter and balanced subsample remain supportive.
- Weakest leave-one-out: drop LIBD → OOF **0.67** (Allen+Linnarsson only; few ACA donors).
- `only_Allen_ACC` OOF=1.0 is **fragile** (4 ACA donors) — do not over-interpret.
- **Pseudobulk DE is NOT a clean territorial claim:** massive FDR count + neuronal genes higher in unsorted ACC vs vessel-enriched Garcia MCA → ambient/tech/enrichment confound.
- Primary Class B claim remains **donor-transferable arm separability**, not gene-level territorial markers from this contrast.

## Locked next scientific steps

1. Prefer within-dataset / enrichment-matched contrasts for confirmatory molecular H1.
2. Vessel-enriched ACC (outreach) still required for enrichment-matched Class B DE.
3. **Do not** re-run soft unmatched DE as territorial markers.
4. Waiting: ATLAS passphrase; RADC clinical delivery. ADNI non-stroke nest frozen separately (not Class B DE).

**Revised meta:** `ROBUST_SUPPORTIVE_separability__DE_contaminated`
