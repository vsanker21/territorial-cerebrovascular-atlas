# External ACA cortex vascular cohort — acquisition verdict

**Date:** 2026-07-24 (updated)  
**Goal:** Obtain an external adult ACA isocortex vascular (endo ± mural) arm for Class B, while waiting on RADC/ADNI/ATLAS.

## TerritoryVoxels prior (HO ROI ranking)

| ROI | Dominant | p(ACA) | Assignment frac | Use |
|---|---|---:|---:|---|
| Cingulate anterior | ACA | 1.00 | **0.72** | Best Class B ACA cortex |
| Paracingulate | ACA | ~1.00 | 0.56 | Include with A32 |
| Frontal medial / A25 proxy | ACA | 1.00 | 0.53 | Include (borderline assign) |
| Cingulate posterior (A23) | ACA | 0.87 | 0.44 | Include ACA-leaning |
| SFG (Yang) | ACA/MCA mix | ~0.53 | 0.24 | Borderzone — not clean ACA |
| Frontal orbital (A14) | **MCA** | 0.00 | 0.40 | **Exclude** |

---

## Vessel-enriched medial/ACC

**Public NO-GO.** No open cohort is both microvessel/BVE-enriched **and** medial ACC/dACC.  
See `results/aca_cohort_hunt/vessel_enriched_ACC_and_LIBD_annotations.md` and outreach draft `docs/outreach/request_vessel_enriched_ACC.md`.

Do **not** call GSE222007 / GSE252921 “ACA.”

---

## Annotated unsorted ACA cortex (GO — acquired)

### Linnarsson / Siletti (CELLxGENE)

- **Object:** `data/processed/aca_cohort_hunt/linnarsson_aca_vascular.h5ad`
- Pool: A24 + A32 + A25 + A23 (exclude A14)
- Vascular 134 / endo 43 / donors 3
- **Verdict:** GO annotated ACA cortex (not vessel-enriched)

### Allen Cortical Cell Atlas — ACC dissection (Bakken)

- **Object:** `data/processed/aca_cohort_hunt/allen_acc_vascular.h5ad`
- ACC vascular **278** (Endo **165** / VLMC **113**) / donors **5**
- Labels: published Allen taxonomy
- **Verdict:** **GO annotated ACA unsorted** (largest public Endo arm)

### LIBD spatialdACC GSE296731 (true dACC)

- **GLOBUS paper-identical labels ACQUIRED** (`sce_Final.Rdata` from `jhpce#dACC_spatial`)
  - Object: `libd_dacc_globus_vascular.h5ad` — Endo **89** / VLMC **159** / donors **10** (matches paper log exactly)
  - Method-matched Azimuth retained as secondary (`libd_dacc_azimuth_vascular.h5ad`)
  - Marker-seed deprecated for Class B

---

## Combined ACA endo pool (current)

`data/processed/aca_cohort_hunt/aca_endo_pool.h5ad`

| Batch | Endo n | Notes |
|---|---:|---|
| Allen_ACC | 165 | published taxonomy |
| LIBD_dACC_globus | 87 | paper-identical Azimuth (2 dropped by endo-marker gate) |
| Linnarsson | 12 | after mild caps |
| **Total** | **264** | **18 donor_keys** |

---

## Other sources (NO-GO for ACA)

| Source | Why |
|---|---|
| GSE222007 / GSE252921 / GSE226606 | Microvessel yes; prefrontal/frontal → MCA-leaning |
| Garcia BVE / Wälchli adult control | Temporal → MCA |
| Garcia ROSMAP in silico | No cingulate region |

---

## Recommended use now

1. Class B ACA arm = **Allen + Linnarsson + LIBD (Azimuth when ready)** vs Garcia MCA tight.
2. Prefer donor OOF over cell AUROC; do not claim vessel enrichment.
3. Continue Zenodo → extract `sce_azimuth.Rdata`; GLOBUS if auth available.
4. Vessel-enriched ACC only via outreach / RADC controlled access.
