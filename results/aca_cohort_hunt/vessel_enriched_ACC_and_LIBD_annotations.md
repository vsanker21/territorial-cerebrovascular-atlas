# Acquisition verdict — vessel-enriched medial/ACC + published LIBD annotations

**Date:** 2026-07-24 (updated)  
**Locked criteria:**
1. **Vessel-enriched** medial / anterior cingulate (ACC/dACC) adult human vascular sn/scRNA — **not** prefrontal-labeled GEO.
2. **Published LIBD spatialdACC cell-type annotations** to replace marker-seed Endo calls.

TerritoryVoxels prior: anterior cingulate = ACA (assign frac ≈0.72).

---

## 1. Vessel-enriched medial/ACC cohort

### Exhaustive public search (GO / NO-GO)

| Candidate | Enrichment | Anatomy | TerritoryVoxels | Verdict |
|---|---|---|---|---|
| GSE222007 / GSE252921 | Microvessel yes | Frontal / prefrontal only | MCA-leaning | **NO-GO as ACA** (on disk for MCA comparator) |
| GSE226606 | Microvessel pilot | Frontal n=2 | MCA-leaning | NO-GO |
| GSE173731 Garcia BVE | BVE yes | Temporal surgical | MCA | NO-GO for ACA |
| GSE173731 ROSMAP in silico | Not physical enrichment | No cingulate among 7 regions | — | NO-GO |
| Wälchli GSE256493 | Sorted/unsorted vascular | Adult control = temporal | MCA | NO-GO |
| LIBD GSE296731 dACC | **Not** vessel-enriched | True dACC | ACA | GO anatomy / NO-GO enrichment |
| Linnarsson A24/A32/A25/A23 | Not enriched | True ACC/MFC | ACA | GO anatomy / NO-GO enrichment |
| Allen Cortical ACC (Bakken) | Not enriched | Dedicated ACC dissection | ACA | **GO annotated unsorted** (acquired) |

### Scientific conclusion (locked)

**No publicly available adult human cohort meets both:** physical **microvessel/BVE enrichment** **and** **medial ACC/dACC** (TerritoryVoxels-ACA-dominant) anatomy.

Claiming prefrontal microvessel GEO or temporal BVE as “ACA Class B” would be anatomically false.

### External next actions only

1. Author outreach (draft: `docs/outreach/request_vessel_enriched_ACC.md`)
2. RADC/ROSMAP controlled access — ask for cingulate vascular-enriched libraries
3. Until then: Class B ACA arm = **anatomically correct unsorted ACC** with honest power limits

---

## 2. Published LIBD spatialdACC cell-type annotations

### What the paper actually published

Shah et al. annotated snRNA-seq with **Azimuth `humancortexref`** (Bakken Human Motor Cortex reference via Azimuth), column `predicted.subclass` → `cellType_azimuth`.

Published counts from their `azimuth.txt` log (n=35,161 QC nuclei):

| Cell type | n |
|---|---|
| **Endo** | **89** |
| **VLMC** | **159** |
| Oligo | 13256 |
| … | … |

Source code (paper-identical method):  
`https://github.com/LieberInstitute/spatialdACC/blob/main/code/snRNA-seq/05_azimuth/01_azimuth_annotation.R`

### Acquisition channels (all pursued)

| Channel | Status | Notes |
|---|---|---|
| **GLOBUS** `jhpce#dACC_spatial` | **DONE** | `sce_Final.Rdata` on disk → Endo 89 / VLMC 159 (paper-identical) |
| **Zenodo** `10.5281/zenodo.15830481` (12.58 GB) | Downloaded | Archive complete; SCE not in zip (symlink only) |
| **Method-matched Azimuth re-run** | DONE (secondary) | GEO matrices; superseded by GLOBUS for Class B |

Published marker tables acquired locally:

`data/raw/aca_cohort_hunt/libd_annotations/azimuth_markers/`  
(`top20_markers_Endo.csv`, `top20_markers_VLMC.csv`, `azimuth_DE_sig_genes_top30.csv`)

### Scientific ranking of label sources

1. **Zenodo/GLOBUS `sce_azimuth.Rdata`** — paper-identical barcodes + labels (**gold**)
2. **Re-run Azimuth `humancortexref` on same matrices** — method-matched, not identical cell set (**silver**)
3. **Marker-seed Endo** — interim only; **do not** call “published annotations” (**bronze / deprecate**)

---

## 3. Annotated ACA cortex vascular (unsorted) — acquired

| Object | n vascular | Endo | Donors | Labels | Path |
|---|---|---|---|---|---|
| Linnarsson ACC | 134 (prior) | 43→caps | 3 | published taxonomy | `linnarsson_aca_vascular.h5ad` |
| Allen Cortical ACC | **278** | **165** | **5** | published Allen taxonomy | `allen_acc_vascular.h5ad` |
| LIBD dACC | pending Azimuth / Zenodo | expect ~89 | 10 | Azimuth | GEO + annotation pipeline |

Allen build verdict: `GO_annotated_ACA_unsorted_not_vessel_enriched`

---

## Locked GO / NO-GO summary

| Need | Status |
|---|---|
| Vessel-enriched **medial/ACC** | **NO-GO publicly** — does not exist in surveyed open sources |
| Published LIBD annotations | **GO paper-identical** — GLOBUS `sce_Final` Endo 89 / VLMC 159 |
| Annotated ACC vascular (unsorted) | **GO** — Linnarsson + Allen ACC + LIBD Azimuth |

### Do not

- Call GSE222007/GSE252921 “ACA”
- Treat marker-seed LIBD as published annotations
- Inflate Class B confidence from cell AUROC when donor n is small
