# Exhaustive public hunt — vessel-enriched ACC for confirmatory DE (H1c)

**Date:** 2026-08-18  
**Criterion (locked):** Enrichment-matched **vessel-enriched / FACS EC / BVE / VINE** libraries from **ACC / medial / clearly ACA-dominant** tissue, usable vs Garcia MCA BVE for confirmatory DE.  
**NO-GO substitutes:** unenriched ACC endo DE; Yang SFG (borderzone); lateral/prefrontal “frontal” microvessel as ACC.

---

## Verdict

**No public vessel-enriched ACC / clean ACA vascular object was found across GEO/SRA, CELLxGENE, local Sun scADbbb, GitHub, BioStudies (noise), or Hugging Face.**

H1c remains gated on **collaborator share** (Wälchli pending) or **new generation**. MultiVINE (already emailed) is **dlPFC**, not ACC — useful for Track B chromatin, **not** confirmatory territorial DE vs MCA.

---

## Sources searched

| Source | Method | Result for vessel-enriched ACC |
|---|---|---|
| NCBI GEO / GDS | Broad + strict eutils queries (ACC/cingulate × BVE/microvessel/VINE/endothelial × snRNA) | Enrichment hits = frontal/prefrontal/hippocampus/temporal — **no ACC** |
| CELLxGENE Curation API | 388 collections; vascular filter | MultiVINE collection present — tissue = **dorsolateral prefrontal cortex** only |
| Sun 2023 scADbbb (local + MIT) | `ROSMAP.VascularCells.meta.rds` | Regions: PFC, Hip, Thal, Ent, AG, Midtemporal, Mammillary — **no ACC/cingulate** |
| Yang VINE GSE163577 | GEO brief | Hippocampus + **superior frontal** (SFG) — SFG = **NO-GO** as clean ACA |
| Matthews GSE222007 / GSE252921 / GSE226606 | GEO | **Frontal / prefrontal** microvessel — PI closed for medial/cingulate |
| Garcia / Heiman BVE | Prior + GEO Garcia∩ACC = 0 | Temporal BVE = **MCA** only |
| Wälchli GSE256493 | GEO / prior | Adult control **temporal** EC = MCA |
| Zenodo API | Query | HTTP 403 (blocked); no alternate ACC hit via web |
| GitHub | Repo search | Cingulate hits = imaging/behavior — **not** vascular snRNA |
| BioStudies | Query | Dominated by adenoid cystic carcinoma “ACC” false positives |
| Hugging Face datasets | Search | 0 vascular brain snRNA ACC hits |
| AWS Open Data | Not separately indexed for this niche | No known ACC-BVE bucket |

---

## Candidate scorecard (H1c)

| Dataset | Enrichment | Region → TerritoryVoxels | H1c confirmatory DE? |
|---|---|---|---|
| Garcia BVE (in hand) | Vessel-enriched | Temporal → **MCA** | MCA arm only |
| Allen / LIBD / Linnarsson ACC endo (in hand) | **Unenriched** | ACC → ACA-ish | Separability yes; **DE NO-GO** |
| Yang VINE SFG | Enriched | SFG → **borderzone** | **NO-GO** as ACA |
| Matthews frontal microvessel | Enriched | Lateral/prefrontal ≠ ACC | Closed by PI |
| MultiVINE CELLxGENE | Enriched multiome | **dlPFC** | **Not ACC**; Track B only |
| Sun scADbbb vascular | In silico / sparse vascular | No ACC region | **NO-GO** for H1c |
| Wälchli adult temporal | FACS EC | Temporal → MCA | MCA / A–V only |
| sgACC cortical snRNA (J Neurosci) | Unenriched | sgACC anatomy | Sparse endo; **NO-GO** |

---

## What *can* still be downloaded (not H1c)

1. **MultiVINE snRNA** (CELLxGENE H5AD, ~628 MB) — dlPFC; pair with snATAC ask to Yang for Track B.  
   Asset: `https://datasets.cellxgene.cziscience.com/25da7ccb-04ca-49d2-b00b-a794ef8c766d.h5ad`  
   Collection: `433700dc-e8a5-48b0-b5cd-beb22f3f88fe`

2. Sun public counts already local under `data/raw/sun2023_scADbbb/` (truncated; full Dropbox wall / Synapse DUC).

---

## Best remaining paths to true vessel-enriched ACC

1. **Wälchli reply** (final nudge sent) — share adult medial/ACC FACS EC.  
2. **Yang reply** — confirm whether any MultiVINE/VINE bank has cingulate (unlikely from public metadata; ask anyway).  
3. **RADC / Greg** — inventory ROSMAP for cingulate BVE beyond Garcia temporal.  
4. **Stanford generation** — ACC autopsy + BVE under IRB (longest).

---

## Scientific lock

Do **not** relabel prefrontal/SFG microvessel or unenriched ACC endo as confirmatory territorial markers. That would violate enrichment-match and TerritoryVoxels locks.
