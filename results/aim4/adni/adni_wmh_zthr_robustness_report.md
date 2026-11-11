# ADNI WMH z-threshold robustness

**Primary lock:** z_thr = **1.5** (unchanged).  
**Grid:** 1.0, 1.25, 1.5, 1.75, 2.0 on existing `flair_mni.nii.gz` (no re-registration).  
**Cross-check at z=1.5:** primary PASS 155 / strict 100 / zero-WMH 108 — matches locked QC.

## Summary

| z | zero-WMH frac | primary n | strict n | Spearman (all) | Spearman (primary) | Δ(M3−M0) primary | mean MCA / ACA share |
|---:|---:|---:|---:|---:|---:|---:|---|
| 1.00 | 0.00 | 267 | 263 | 0.283 | 0.283 | 0.021 | 0.242 / 0.040 |
| 1.25 | 0.12 | 233 | 185 | 0.492 | 0.461 | 0.162 | 0.220 / 0.043 |
| **1.50** | **0.40** | **155** | **100** | **0.508** | **0.521** | **0.424** | **0.207 / 0.071** |
| 1.75 | 0.61 | 98 | 61 | 0.378 | 0.599 | 0.415 | 0.217 / 0.102 |
| 2.00 | 0.65 | 86 | 52 | 0.341 | 0.580 | 0.182 | 0.221 / 0.101 |

## Rationale (why z=1.5 stays locked)

1. **z=1.0 eliminates zeros but destroys validity** — median WMH voxels explode (1178); Spearman vs UCD collapses to ~0.28; Δ(M3−M0)≈0.02. Consistent with false-positive WM hyperintensity dominating the map.
2. **z=1.25** improves sensitivity (12% zero) with moderate concordance (ρ≈0.46) but weaker partition ΔR² than z=1.5.
3. **z=1.5** maximizes all-subject Spearman among the grid (~0.51) with the strongest primary Δ(M3−M0) (~0.42) while retaining MCA-dominant share. Zero-WMH (~40%) is the cost of specificity vs UCD≥3 cc — QC drops those subjects by design.
4. **z≥1.75** further cuts n and worsens all-subject concordance; primary-pass ρ can look high by selection of remaining lesions.

MCA > ACA share is **qualitatively stable** across the grid (sensitivity, not a new claim).

## Locked interpretation

- Primary Aim 4 non-stroke table remains QC-pass at **z=1.5**.
- M3−M0 ΔR² is partition/concordance, **not** molecular H3.
- Gene-level Class B territorial markers remain blocked (unsorted ACC vs BVE MCA); do not re-run soft unmatched DE.
- Does not unlock ATLAS stroke Aim 4.

Detail: `results/aim4/adni/adni_wmh_zthr_robustness_detail.csv`  
JSON: `adni_wmh_zthr_robustness_summary.json`
