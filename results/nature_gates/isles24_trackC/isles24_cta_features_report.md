# ISLES'24 CTA features (registered to NCCT)

**Verdict:** `CTA_INSAMPLE_DELTA_BUT_CV_NONPOSITIVE__NOT_CONFIRMATORY`
- n=40
- method: SimpleITK_affine_CTA_to_NCCT_then_TopCoW_intensity

## Models
- **age**: R²_in=0.0017761706137338074; R²_cv5=-0.06471028975345527; n=40
- **age_CoW_LVO**: R²_in=0.12809473460082255; R²_cv5=-0.05247797674970034; n=40
- **age_CTA**: R²_in=0.12083154903065585; R²_cv5=-0.24818323805339348; n=31
- **age_CoW_LVO_CTA**: R²_in=0.27638251216338994; R²_cv5=-0.3960381064277494; n=31

## ΔR²
```
{
  "CTA_minus_age": {
    "n": 31,
    "R2_base": 0.007072845354598067,
    "R2_full": 0.12083154903065585,
    "delta": 0.11375870367605778,
    "ci95": [
      0.03636903068255265,
      0.49052936387114315
    ],
    "ci_excludes_0": true,
    "p_boot_delta_le_0": 0.0
  },
  "CTA_minus_CoWLVO": {
    "n": 31,
    "R2_base": 0.1307595118633489,
    "R2_full": 0.27638251216338994,
    "delta": 0.14562300030004105,
    "ci95": [
      0.03613276012413408,
      0.5145961620535829
    ],
    "ci_excludes_0": true,
    "p_boot_delta_le_0": 0.0
  },
  "full_minus_age": {
    "n": 31,
    "R2_base": 0.007072845354598067,
    "R2_full": 0.27638251216338994,
    "delta": 0.26930966680879187,
    "ci95": [
      0.12446682395301673,
      0.6852635211617579
    ],
    "ci_excludes_0": true,
    "p_boot_delta_le_0": 0.0
  }
}
```

## Spearman
```
[
  {
    "feature": "cta_mean_MCA",
    "rho": -0.2735459662288931,
    "p": 0.0876467696526608,
    "n": 40
  },
  {
    "feature": "cta_asym_MCA",
    "rho": -0.05524193548387097,
    "p": 0.7678716935426563,
    "n": 31
  },
  {
    "feature": "cta_LVO_over_CoW",
    "rho": -0.13189493433395877,
    "p": 0.41720583909785525,
    "n": 40
  },
  {
    "feature": "cta_cow_range_p90_p10",
    "rho": -0.15478424015009384,
    "p": 0.34025408521927125,
    "n": 40
  },
  {
    "feature": "cta_mean_LVO",
    "rho": -0.1853658536585366,
    "p": 0.25215150317852175,
    "n": 40
  }
]
```
