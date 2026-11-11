# Master execution report — Tracks A / B0 / C / Aim3 / B gate

**Date:** 2026-08-18  
**Venue program:** Nature Neuroscience path locked in `docs/13_Nature_NatNeurosci_readiness.md`  
**Prereg amendment:** `docs/04b_preregistration_Nature_amendment_v2.md`  
**Concept note:** `territorial-identity-concept-note.docx`  
**TRCI/BCI:** application layer only — **out of** Nature identity package  

---

## Nature path decision (2026-08-18)

**Lead claim:** Track **C–primary** (architecture/perfusion-compatible territorial vascular state) + endothelial separability + territory-normalized geography.  
**Not lead:** Track B developmental mural identity (HOX RNA NO-GO; needs chromatin + ACC).  
**Not Nat Neurosci alone:** Track A atlas-null (Stroke-tier methods).  

**Readiness:** ~40% for Nat Neurosci submission; program is GO if external gates open.

### Nature rigor update — matched Midtemporal (MCA) vs Entorhinal (PCA)
- Endo donor-region OOF AUROC ≈ **0.56** (CI ≈0.45–0.67) — below confirmatory bar  
- Endo **hypoxia** module: PCA>MCA, d≈0.51, p≈0.025 — exploratory until replication  
- Endo **flow_shear**: MCA>PCA trend (CI excludes 0) but d≈0.47 **below MDE 0.5**  
- Pericyte **mural_contractile**: PCA>MCA, d≈0.67, p≈0.004 — exploratory  
- Most other modules null  

---

## Executive routing (locked)

| Track | Status | Verdict |
|---|---|---|
| **A — Imaging** | **Executed** | `SHARED_GEOMETRIC_TERRITORIAL_DOMINATES__TRACK_C_ALIGNED` |
| **B0i — mural counts** | **Executed** (local Garcia ROSMAP multiregion) | `PASS` among region-present donors (matched pair median mural = 8) |
| **B0ii — HOX RNA** | **Executed** | `HOX_RNA_UNDETECTABLE` → lineage needs chromatin; RNA gate NO-GO |
| **C — flow programs** | **Executed** | Flow module null anterior vs posterior / MCA vs PCA on this object |
| **B — developmental** | **Closed** | No vessel-enriched ACC; no chromatin; HOX RNA NO-GO |
| **Aim 3 — RADC** | **Scaffold only** | RRA-539487; awaiting DUA delivery |
| **TRCI/BCI** | **Frozen separate** | Not used in A/B/C claims |

**Recommended program path:** Finish **Track A** paper (volume-norm + honest null + selection note) → keep **Aim 3** warm for RADC → do **not** fund Track B until ACC + chromatin open → treat Track C as interpretive frame (architecture/shared variance), not as positive flow-module proof on current endo object.

---

## Track A (detail)

**Scripts:** `analysis/tracks/trackA_volume_norm_geometric_null.py`  
**Outputs:** `results/tracks/trackA/`

### What was done
1. **A1 Volume-normalized burden** — used existing `lesion_burden_in_*` / `wmh_burden_in_*` (lesion voxels ÷ territory voxels) for ATLAS (n=1445), ISLES’22 (n=247), ADNI WMH (n=155).
2. **A2 Geometric null** — atlas-derived per territory: log volume, mean distance to ventricles (HO), GM/WM means, MRA probability/radius.
3. **A3 Variance partition** — commonality analysis on subject×territory long tables.
4. **A4 Selection structure** — documented ATLAS rehab skew, ISLES LVO/thrombectomy skew, ADNI non-stroke.

### Key quantitative results
| Cohort | R² (between-territory) | Unique geometry | Unique territory | Shared |
|---|---:|---:|---:|---:|
| ATLAS | 0.236 | 0 | 0 | 0.236 |
| ISLES22 | 0.062 | 0 | 0 | 0.062 |
| ADNI_WMH | 0.123 | ~0 | 0 | 0.123 |

MCA vs ACA volume-normalized burden: ATLAS p≈1.1e−91; ISLES p≈5.3e−17; ADNI WMH p≈0.23 (ns).

### Scientific lock (critical)
Atlas-constant geometry features are **linearly equivalent** to territory indicator codes (same column space for 4 territories). Therefore unique geometry = unique territory = 0 and **all explained variance is shared**. This is not a bug — it is the concept note’s prediction that branching anatomy and territory labels are collinear. **Claiming “territory beyond geometry” requires subject-level architecture (individual MRA/angiography) or molecular features.** Track A today supports: (i) volume-normalized MCA>ACA in stroke cohorts; (ii) shared geometric–territorial structure; (iii) methods paper / Track C framing — **not** independent territorial biology after geometry.

---

## B0 gates (detail)

**Scripts:** `analysis/tracks/b0_trackC_mural_hox_flow.py`  
**Primary object:** `data/processed/garcia/garcia_rosmap_vascular.h5ad` (11,689 cells; 214 donors; Endo/Per/SMC/Fib; 7 regions)

### Data acquisition attempts (programmatic)
| Source | Result |
|---|---|
| `http://compbio.mit.edu/scADbbb/ROSMAP.VascularCells.counts_full.rds` | Downloaded **exact 32 MiB** — truncated/unreadable |
| Dropbox share (HTML) | Login wall; not machine-downloadable without auth |
| `http://compbio.mit.edu/scBBB/ROSMAP.VascularCells.counts.rds` (49 MB) + `meta.rds` | **Downloaded successfully** (Garcia-era / overlapping vascular release) |
| Synapse syn51015750 / syn52293417 | Requires DUC — not anonymous |
| GitHub `nasunmit/scadbbb` | Code only; points to MIT/Dropbox |

### B0i
Among donors **with that region profiled**, median mural (Per+SMC) ≈ 8–12 for Midtemporal, Entorhinal, PFC, Hippocampus, Thalamus. Matched pair Midtemporal vs Entorhinal: **PASS** (≥50% donors with ≥5 mural cells).

### B0ii
HOX cluster RNA in adult mural cells is **essentially undetectable** (module mean ≈ 0). This matches the concept note’s statement that lineage may persist as **chromatin**, not adult RNA — so the RNA HOX gate is a **NO-GO**, not a positive lineage test.

---

## Track C (detail)

Endothelial shear/flow module (KLF2/KLF4/NOS3/…) — anterior vs posterior p≈0.82; MCA-region vs PCA-region p≈0.38 → **inconclusive** on this object.  
**Interpretive support for Track C still comes from Track A shared variance + known mesoderm endothelium (ACA↔MCA endo separability ≠ crest lineage).** Do not claim flow-module territorialization from today’s null.

---

## Aim 3 RADC

**Script:** `analysis/tracks/aim3_radc_scaffold.py`  
**Status:** `GATED_AWAITING_RADC_DUA_DELIVERY__SCAFFOLD_ONLY`  
Garcia sample table has **219 projids** ready to link. Place Standard extract under `data/raw/radc/` when Rush delivers.

---

## Track B

**Do not spend.** Requirements unmet: vessel-enriched ACC (Heiman/Matthews closed; Wälchli pending), chromatin/snATAC, HOX RNA gate failed, full Sun 22.5k annotated release not anonymously downloadable.

---

## TRCI / BCI

Unchanged frozen application layer under `results/bci_trci/` and `docs/outreach/Perfusion_Territory_Neurointerfaces_concept_note.pdf`. **Excluded** from Track A/B/C identity claims per concept-note housekeeping.

---

## Immediate next actions

1. **Track A write-up** for Stroke-style methods/geography paper (volume-norm MCA>ACA + collinearity/shared-variance lock + selection bias).  
2. **Push RADC DUA** (RRA-539487) — only path to Aim 3.  
3. **Optional:** authenticate Dropbox/Synapse for full Sun scADbbb counts to re-run B0 on published 22,514-cell object.  
4. **Do not** open Track B experiments until ACC + chromatin exist.  
5. Keep TRCI/BCI collaborator ask separate.
