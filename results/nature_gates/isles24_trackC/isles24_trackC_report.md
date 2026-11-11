# ISLES'24 Track C expanded pilot

**Verdict:** `ISLES24_EXPANDED_WITH_COW_FEATURES`
- Cases: 80; reg OK: 80; QC pass: 78
- Nested: `{
  "outcome": "lesion_share_MCA",
  "qc_note": "QC-pass only",
  "n": 78,
  "models": {
    "age": {
      "R2": 0.008145226953355178,
      "n": 78
    },
    "age_LVOvol": {
      "R2": 0.008553516584506204,
      "n": 78
    },
    "age_CoWcomplete": {
      "R2": 0.01383044863186933,
      "n": 78
    },
    "age_CoWMCA": {
      "R2": 0.031747987829572666,
      "n": 78
    },
    "age_LVOvol_CoWcomplete_CoWMCA": {
      "R2": 0.03221610754668358,
      "n": 78
    },
    "age_full_angio": {
      "R2": 0.03284376383876497,
      "n": 78
    }
  },
  "delta_R2_vs_age": {
    "age_LVOvol": 0.00040828963115102646,
    "age_CoWcomplete": 0.005685221678514152,
    "age_CoWMCA": 0.023602760876217488,
    "age_LVOvol_CoWcomplete_CoWMCA": 0.024070880593328403,
    "age_full_angio": 0.024698536885409794
  },
  "mean_lesion_shares": {
    "ACA": 0.0788582889651224,
    "MCA": 0.772138744651266,
    "PCA": 0.04973352381529854,
    "VB": 0.016747212234711893
  },
  "verdict": "EXPLORATORY_EXPANDED__NOT_CONFIRMATORY",
  "note": "Subject-level CoW/LVO features; confirmatory needs larger n + registration QC audit."
}`

## Features
- TopCoW multi-label presence/volumes; MCA laterality; completeness
- LVO volume + centroid x (binary LVO often invariant)
- Lesion shares on ArterialAtlas_level2 majors

## Locks
- Not TV_BMM; exploratory; TRCI/BCI excluded
