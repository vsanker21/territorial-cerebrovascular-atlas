# Track C confirmatory path — ISLES'24 subject-level CTA + lesions

**Date:** 2026-08-18  
**Status:** Spec locked; full train archive ~**99 GB** (`train.7z` on Zenodo 17652035) — selective acquisition required  
**ADNI subject TOF:** **NO-GO** (protocol + local MRI key)

---

## Why ISLES'24 (highest scientific accuracy)

| Need for Nature Track C | ISLES'24 | ADNI | OpenNeuro ds005096 |
|---|---|---|---|
| Subject-level angio | **CTA** (+ CTP) | No TOF/MRA | TOF-MRA yes |
| Territorial lesion maps | Follow-up DWI **lesion masks** | WMH yes; no angio | Aneurysm only |
| Vessel architecture labels | **LVO masks** + **Circle of Willis** masks | No | Parent-artery aneurysm segs |
| Same subjects | Yes | N/A | N/A |

**Pre-specified confirmatory test (when ≥N subjects with CTA features + lesion territory shares):**

1. Map lesions to ArterialAtlas_level2 / TerritoryVoxels (not TV_BMM).  
2. Extract subject angio features from CTA (LVO presence/side, CoW completeness, proximal calibre proxies).  
3. Nested models: `territory_burden ~ angio_features + territory` vs reduced models.  
4. Confirmatory only if unique territory ΔR² ≥ MDE after subject angio; else architecture-mediated (Track C win) or null.

---

## Acquisition plan (do not download full 99 GB blindly)

1. Prefer challenge / Hugging Face selective case pulls (`hugging-science/isles24-stroke`, CC BY-NC-SA).  
2. Minimum pilot: **n=10–20** train cases with `cta`, `lvo-msk`, `cow-msk`, `lesion-msk`.  
3. Drop path: `data/raw/isles24/`  
4. License: non-commercial — compatible with academic Nature path; no commercial redistribution.

## Local stubs

- Acquisition status placeholder: extend `data/raw/isles/` notes when first cases land.  
- Analysis entrypoint (to implement on first cases): `analysis/tracks/isles24_trackC_angio_lesion.py` (pending data).

## Explicit NO-GOs

- Claiming ADNI atlas-MRA null as subject-level architecture.  
- Using OpenNeuro aneurysm TOF as confirmatory lesion–architecture mediation.  
- Folding TRCI/BCI into this mechanism arm.
