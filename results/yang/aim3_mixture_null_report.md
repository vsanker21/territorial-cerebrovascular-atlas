# H2 / Aim 3 — convex mixture null

**Primary verdict:** `MOSTLY_MIXTURE_WITH_STRUCTURED_RESIDUAL_exploratory`

## Primary — within-Allen (flank-separating genes)
Borderzone: MCA_leaning Endo. Flanks: ACA (ACC) + MCA (MTG+ANG).
- status=ok
- w_ACA=0.21, w_MCA=0.79
- frac_var_explained=0.9405038247790525, cosine=0.9699324801240664
- p(FVE ≥ obs | gene permutation)=0.0
- donors border/ACA/MCA=4/4/6
- n_genes_selected=500
Residual modules:
  - BBB_tight_junction: resid=0.017560156171934405, p=0.17666666666666667
  - angiogenesis: resid=-0.3142307682894171, p=0.0
  - ECM_basement: resid=-0.4846575389005658, p=0.0
  - hypoxia_stress: resid=-0.2778393166874432, p=0.0
  - inflammation: resid=-0.015854101450079008, p=0.49666666666666665
  - autoregulation_proxy: resid=-0.39017549658194195, p=0.0

## Exploratory — Yang SFG (Aim3) vs Allen ACA + Garcia MCA
**Caveats:** cross-dataset; enrichment asymmetry; not confirmatory H1/H2.
- status=ok
- w_ACA=0.14, w_MCA=0.86
- frac_var_explained=0.6488945171100987, cosine=0.8056361606320045
- p(FVE ≥ obs)=0.0
- donors border/ACA/MCA=4/4/9
Residual modules:
  - BBB_tight_junction: resid=-0.1604948870690524, p=0.12666666666666668
  - angiogenesis: resid=0.007349250680041955, p=0.08
  - ECM_basement: resid=-0.17193729697541701, p=0.056666666666666664
  - hypoxia_stress: resid=-0.2855834120605065, p=0.01
  - inflammation: resid=-0.15859591934558417, p=0.10666666666666667
  - autoregulation_proxy: resid=-0.057350497847174366, p=0.6566666666666666

## Locked interpretation
- Primary within-Allen: flanks explain most borderzone variance (high FVE), with residual pathway modules that can remain non-null after mixture.
- Verdict `MOSTLY_MIXTURE_WITH_STRUCTURED_RESIDUAL_exploratory` = mixture-dominant + structured leftover (exploratory H2-compatible residual).
- Verdict `STRUCTURED_RESIDUAL_H2_COMPATIBLE_exploratory` = low FVE + coherent residual (stronger H2-compatible).
- Do not call Yang SFG clean ACA; Yang arm is exploratory (cross-dataset / enrichment asymmetry).
