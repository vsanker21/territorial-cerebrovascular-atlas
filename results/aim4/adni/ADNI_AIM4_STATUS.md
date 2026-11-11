# ADNI Aim 4 — QC + scale + nest + robustness (frozen)

**Date:** 2026-07-27 (sensitivities added; primary nest unchanged)

## 1. Scale — done
267/267 paired FLAIR+T1+UCD≥3 cc OK. Local vs UCD Spearman (all) ρ≈0.508. Zero local WMH at z=1.5: 108.

## 2. QC — locked
| Set | n |
|---|---:|
| Input | 267 |
| **Primary PASS** | **155** |
| Strict PASS | 100 |
| FAIL | 112 (mostly ZERO_WMH) |

## 3. Nest — frozen
| | Primary n=155 | Strict n=100 |
|---|---:|---:|
| Spearman local↔UCD | 0.521 | 0.573 |
| M1 MRA null R² | 0.000 | 0.000 |
| Δ(M3−M0) | 0.424 | 0.448 |
| Mean share MCA / ACA | 0.207 / 0.071 | 0.318 / 0.107 |
| Subject-level MCA>ACA | 0.452 | — |

**Lock:** M3 = partition/concordance, not molecular H3. Geography = **mean** MCA share > mean ACA (not majority-subject claim).

## 4. z-threshold robustness — done
Primary z=1.5 justified: z=1.0 kills concordance (ρ≈0.28, Δ≈0.02); z=1.5 best all-subject Spearman + strongest Δ among mid thresholds; MCA>ACA share stable.  
Report: `adni_wmh_zthr_robustness_report.md`

## 5. Sensitivities (2026-07-27) — descriptive only
| Analysis | Result |
|---|---|
| Sex / APOE | MCA share unchanged (sex p≈0.97; e4 p≈0.94); APOE via ADNIMERGE2 APOERES.rda |
| DIAGNOSIS | CN/MCI mean MCA>ACA; CN vs MCI p≈0.15; Dementia n=9 exploratory |
| BorderZone_ProbAve | Mean frac WMH BZ>0.5 ≈0.036 — unlabeled WMH not mainly borderzone |
| Allen modules × WMH shares | `EXPLORATORY_ECOLOGICAL_n3__NOT_CONFIRMATORY_H3` |

## 6. Paper freeze
Integrated into claims, synthesis, prereg, captions (panel `06_adni_aim4_wmh_nest`).
