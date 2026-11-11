# ISLES'24 Track C — honest scale-up synthesis (2026-08-18)

## Scale dependence (scientific lock)

| Cohort | QC n | ΔR² (full angio − age) for MCA share | 5-fold CV R² (full) | Interpretation |
|---|---|---|---|---|
| First 40 | 40 | ≈0.17; bootstrap CI excludes 0 | **negative** | Unstable / overfit |
| Expanded 80 | 78 | ≈0.025; bootstrap CI still excludes 0 | **still negative** | Tiny in-sample angio signal; **not predictive** |

## Accurate Nature wording

Subject-level CoW/LVO features show a **small positive in-sample association** with MCA lesion share beyond age, but **cross-validated R² remains ≤0**. This does **not** yet support a confirmatory Track C claim that subject angio explains territorial lesion geography.

## What still helps the program

1. Imaging geography (ATLAS/ISLES'22 MCA-dominant shares) remains solid.  
2. Atlas-constant MRA null remains locked (collinear).  
3. ISLES'24 proves the **pipeline** (TopCoW + LVO + ArterialAtlas_level2) is executable.  
4. Mechanism for Nature still needs either stronger subject angio predictors (CTA calibre/collateral grades), larger n with CV+, or molecular gates (ACC / snATAC / RADC).

## Do not claim

- That n=40 ΔR²≈0.17 is confirmatory.  
- That CV-negative models triage clinical angio.  
- TRCI/BCI (excluded).

Outputs: `isles24_trackC_inference_report.md`, subject table n=80.
