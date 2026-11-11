# Class B hardening (ambient + modules + zonation status)

**Verdict:** `SEPARABILITY_ROBUST_TO_AMBIENT`

| Analysis | Donor OOF | Cell CV |
|---|---:|---:|
| Baseline Harmony PCs | 0.787 | 0.826 |
| Ambient residualized PCs | 0.724 | 0.666 |
| Non-neuronal modules only | 0.599 | 0.863 |
| Harmony + ambient as features | 0.852 | — |

## Ambient imbalance

- Neuronal mean ACA=1.919 vs MCA=0.147 (MW p=2.62e-73)
- Oligo mean ACA=0.243 vs MCA=0.367

## Zonation

Class B object lacks zonation. Within Garcia ROSMAP MCA, capEndo vs a/vEndo cell CV AUROC=0.980 (labels real). Cannot zonation-condition cross-dataset Class B until ACA arm has A-V labels.

Claims paragraph: `docs\10_class_b_claims_paragraph.md`
