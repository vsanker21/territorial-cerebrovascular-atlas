# Moonshot stack — final scientific synthesis (public data)

**Date:** 2026-08-18  
**Code:** `analysis/moonshot/`  
**Results:** `results/moonshot/`  
**Claim reframe:** coordinate + genetics + geography + falsifiers — **not** confirmatory ACC markers.

---

## Venue gate (locked)

**Honest venue after this stack: *Brain* / *Nature Communications*.**  
**Not ready for *Nature Neuroscience* on public data alone.**

Reasons:
1. Territory-as-GWAS-filter **fails** Class B falsifier (GWAS modules ≈ random modules for ACA vs MCA).
2. CoW topology **fails** nested OOS (negative R² on all territorial share outcomes).
3. Track B developmental memory language is **KILL** (ΔAUROC from developmental PC = 0).
4. Vessel-enriched ACC for confirmatory H1c still **absent**.

Genetics **cell-type** biology is real (SMC × WMH/cSVD genes) and belongs in the paper — but it does **not** currently prove a territorial filter beyond cell type.

---

## Data acquired programmatically

| Object | Source | Path |
|---|---|---|
| MultiVINE snRNA (~628 MB, 65,479 nuclei) | CELLxGENE | `data/raw/multivine/multivine_snrna.h5ad` |
| GWAS Catalog + curated cores | EBI API + literature cores | `data/raw/gwas/` |
| Matthews BVE P/V (40×10x h5) | GEO GSE222007 RAW | `data/raw/matthews_bve/raw/` |
| Wälchli fetal CNS EC + adult temporal | GEO GSE256493 → mtx | `data/processed/walchli/` |
| Class B / Garcia / Allen / Yang / ISLES24 / ATLAS / ADNI | already local | (prior ingest) |

---

## Results by analysis (honesty locks enforced)

### 1 — GWAS × MultiVINE territory filter
- Vascular nuclei: **39,582** (dlPFC only — not ACC).
- **Best enrichment:** Smooth_Muscle_Cell × WMH_core — COL4A1, COL4A2, NBEAL1, NOTCH3 (hypergeom p≈2.2×10⁻⁴, FDR≈0.0018).
- **Falsifier (locked):** Class B ACA vs MCA donor-OOF AUROC from GWAS modules **0.916** vs random-module null **0.916±0.010** → **`gwas_beats_null = False`**.
- Matched Garcia Midtemporal vs Entorhinal GWAS AUROC ≈ **0.49** (null).
- **Verdict:** Cell-type genetics supported; **territory-as-extra-filter falsified** at current resolution.

### 2 — Counterfactual enrichment transport (NOT H1c)
- Matthews paired V vs P endo AUROC ≈ **0.93** (n_pairs=10); mean-shift and OT shifts align (cosine≈0.91).
- Class B ACA vs MCA AUROC: raw **0.82** → after transport **0.95–0.96**; distance to MCA centroid **increased**.
- **Interpretation (sensitivity only):** Applying a frontal BVE operator does **not** collapse ACA into MCA. Residual separability persists — compatible with a territory-linked state **or** with operator misspecification. **Forbidden:** calling this confirmatory territorial markers.

### 3 — Developmental fate memory
- Fetal CNS EC (n=5,000 subsample) aligned to adult Wälchli temporal + Class B (14,904 shared genes).
- Developmental PC explains fetal vs adult (|r|≈0.30), but territory↔dev_score |r|≈0.04 after batch residualization.
- Donor OOF AUROC with vs without dev PC: **0.980 = 0.980** (Δ=0).
- **Track B language: KILL.**

### 4 — CoW topology IV + conformal
- n=78 QC-pass ISLES’24.
- Nested OOS R² for MCA share / ACA share / MCA−ACA: **all negative**; CoW features do **not** beat demographics.
- 90% conformal coverage ≈0.76–0.82 (undercover slightly).
- **Track C topology: FAIL_OOS** (matches prior intensity-model failure; honest null).

### 5 — Wasserstein territorial manifold
- W₂(ACA,MCA)≈30 vs within-MCA≈2.2 (strong geometry in Harmony space).
- Conformal territory labels: accuracy≈0.82, set coverage≈0.92, mean set size≈1.55.
- Yang SFG OT mixture FVE≈**0.54** (prior module FVE≈0.94); SFG much closer to MCA (W₂≈10) than ACA (W₂≈56) — reinforces **SFG ≠ clean ACA**.

### 6 — Hierarchical Bayesian ecological H3
- Mean shares: ATLAS/ISLES/ADNI MCA-dominant (as before).
- Territory-level hierarchical β for module→share: **CIs include 0** at n_territories=4.
- **Label: ECOLOGICAL_NOT_CONFIRMATORY_MOLECULAR_H3.** Kills n=3 Spearman theater by replacing it with posterior uncertainty.

### 7 — Foundation-style probes
- Harmony PCA donor OOF AUROC ≈ **0.82** (donor-label shuffle null ≈ **0.43**).
- GWAS-module probe AUROC ≈ **0.92** (shuffle null ≈ **0.35**) — beats null, but not independent of #1 genesets.
- GenePT embeddings: download 404 → documented; Harmony+GWAS probes are primary.
- Alignment Spearman(coef vs MultiVINE enrichment) ≈0.14 (n.s.).
- **Methods supplement only**; does not rescue Nat Neurosci alone.

---

## What this means scientifically

The public-data moonshot **succeeded as a rigor program**: it produced publishable positives **and** decisive nulls.

| Claim | Status |
|---|---|
| Vascular cell-type GWAS enrichment (SMC–WMH/cSVD) | **Support** |
| Territory adds GWAS signal beyond cell type / noise | **Falsified** |
| Enrichment-matched confirmatory ACA markers | **Still gated** |
| Developmental HOX/fate memory in adult territory | **Kill Track B language** |
| Subject CoW topology predicts lesion share OOS | **Falsified** |
| Territory as OT coordinate + conformal labels | **Support (methods/theory)** |
| Ecological molecular→lesion vulnerability | **Inconclusive (wide posteriors)** |

**Manuscript spine that remains honest:**  
Tier-1 donor-transferable endothelial separability + MCA-dominant lesion/WMH geography + OT manifold + SMC GWAS programs + explicit falsifiers (HOX, CoW OOS, territory-GWAS filter, unmatched DE ambient).

**Still required for Nat Neurosci:** vessel-enriched ACC (H1c) and/or RADC phenotype linkage and/or subject-level angio that actually predicts OOS.

---

## Reproduce

```bash
cd "F:/Cerebrovascular SCS/analysis/moonshot"
python run_all_moonshot.py
```

Individual: `python 01_gwas_multivine.py` … `python 07_foundation_probes.py`
