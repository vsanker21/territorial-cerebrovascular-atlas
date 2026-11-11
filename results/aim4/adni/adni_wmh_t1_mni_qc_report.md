# ADNI T1–MNI WMH QC

**Primary PASS:** 155 / 267
**Strict PASS (exclude FLAG_LOW_TV):** 100 / 267
**FAIL:** 112

## Locked rules
| Code | Rule | Rationale |
|---|---|---|
| FAIL_ZERO_WMH | local WMH voxels = 0 | No territory numerator |
| FAIL_EMPTY_REG | FLAIR/T1 MNI occupancy < 15% brain | Empty/collapsed registration |
| FAIL_BRAIN_DICE | Dice(FLAIR, MNI brain) < 0.25 | Gross misregistration |
| FAIL_UCD_DISCORD | UCD≥10 cc & local voxels<20 | Missed substantial burden vs validated volume |
| FLAG_LOW_TV | assigned TV frac < 0.10 | Keep in primary (borderzone/unlabeled WM); drop in sensitivity |

- Fail reasons: `{'FAIL_ZERO_WMH': 99, 'FAIL_ZERO_WMH;FAIL_UCD_DISCORD': 9, 'FAIL_UCD_DISCORD': 4}`
- Median Dice (pass): 0.49233326807124683
- Mean WMH share ACA/MCA/PCA/VB (pass): 0.071/0.207/0.040/0.008

## Spot-check note
Automated Dice replaces subjective screenshot QC; FAIL_BRAIN_DICE catches the same failures.
QC thresholds locked on n=40 discovery; applied unchanged to the full paired set.

Outputs: `results\aim4\adni\adni_wmh_t1_mni_qc.csv`, `adni_wmh_t1_mni_qc_summary.json`
