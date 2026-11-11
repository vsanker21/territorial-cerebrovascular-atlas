# ISLES'24 Track C — formal inference (exploratory)

**Date:** 2026-08-18
**n (QC-pass):** 78
**Verdict:** `EXPLORATORY_ANGIO_DELTA_CI_EXCLUDES_0__NOT_CONFIRMATORY_ALONE`

## Primary: MCA lesion share — ΔR² (full angio − age)
`{"outcome": "lesion_share_MCA", "contrast": "age_full_minus_age", "n": 78, "R2_base": 0.008145226953355289, "R2_full": 0.03284376383876475, "delta": 0.02469853688540946, "ci95": [0.009610298421491106, 0.13749107365184077], "p_boot_delta_le_0": 0.0, "ci_excludes_0": true}`

## MCA share model R² (in-sample / 5-fold CV)
- **age**: R²_in=0.008; R²_cv5=-0.060
- **age_LVOvol**: R²_in=0.009; R²_cv5=-0.075
- **age_CoWcomplete**: R²_in=0.014; R²_cv5=-0.067
- **age_CoWMCA**: R²_in=0.032; R²_cv5=-0.034
- **age_angio**: R²_in=0.032; R²_cv5=-0.071
- **age_full**: R²_in=0.033; R²_cv5=-0.107

## Files
- `isles24_trackC_inference_models.csv`
- `isles24_trackC_inference_deltaR2.csv`
- `isles24_trackC_inference_spearman.csv`

## Interpretation
Positive angio ΔR² with CI excluding 0 supports Track C architecture mediation as an **exploratory** mechanism signal. Confirmatory package still needs larger n, registration audit, and (for Nature) linkage to molecular gates.
