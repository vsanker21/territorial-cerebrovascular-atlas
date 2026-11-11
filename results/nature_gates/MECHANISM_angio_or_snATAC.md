# Mechanism gate — subject-level angio/MRA **or** snATAC (2026-08-18)

**Nature lock:** Track C architecture/perfusion is the leading mechanism for endothelium. Track B developmental mural identity is **secondary and gated** — open **only if** chromatin (snATAC / methylation) supports lineage marks after HOX RNA NO-GO.

**Scientific rationale**
- Atlas-constant MRA geometry is **collinear** with territory labels → cannot claim “territory beyond geometry.”
- Endothelium is mesoderm-derived in all territories → endo separability ≠ crest lineage.
- HOX RNA undetectable in current mural object → lineage needs **chromatin accessibility**, not more RNA DE.

---

## Arm 1 — Subject-level angio / MRA (Track C decisive)

| Source | What it is | Fit for Nature mechanism | Verdict |
|---|---|---|---|
| ADNI TOF-MRA (controlled) | Subject-level angio + WMH | Best if linked to same donors with lesion/WMH territory shares | **Pursue under existing ADNI access** — subject calibre/branching/collateral |
| OpenNeuro ds005096 (RBWH TOF-MRA) | n=63 aneurysm TOF-MRA + segmentations | Subject-level angio geometry **without** territorial stroke/WMH lesions | **Methods pilot only** — calibrate vessel metrics; **not** lesion–architecture mediation |
| SOOP / ATLAS / ISLES (in hand) | Lesions ± outcomes | Strong geography; **no** subject angio | Keep for geography; mechanism incomplete |
| ICPSR 38464 clinical stroke MRI | Large DWI lesions | Geography; typically no TOF | Secondary geography only |

### Pre-specified analysis (when subject-level angio exists)

1. Extract subject-level features: proximal calibre, branching density, collateral score, territory supply volume (from angio, **not** atlas dummies).  
2. Nested models: `burden ~ angio_features + territory_dummies` vs `burden ~ angio_features` vs `burden ~ territory_dummies`.  
3. **Confirmatory** only if unique territory ΔR² ≥ MDE after subject angio; else interpret as architecture-mediated (Track C win) or null.

### Immediate executable step

- **ADNI inventory (2026-08-18):** **NO-GO** for subject TOF/MRA — 0 hits in 92,246 MRI key rows; ADNI core protocol has no TOF-MRA. Aim 4 nest remains atlas-constant only. Report: `results/nature_gates/adni_subject_tof_mra_inventory.md`.  
- **OpenNeuro ds005096 pilot:** 3 subjects parent-artery volumes extracted (methods only). Report: `results/nature_gates/openneuro_tof_mra_pilot.md`.  
- **Confirmatory Track C:** **ISLES'24** (CTA + CoW/LVO + lesion masks). Spec: `results/nature_gates/ISLES24_TrackC_path.md`. Full Zenodo train ≈99 GB — acquire selectively.

---

## Arm 2 — snATAC / MultiVINE (Track B only if chromatin works)

| Source | Modality | Access | Nature role |
|---|---|---|---|
| MultiVINE-Seq (Reid / Yang lab; Neuron report / Alzforum 2025) | snRNA + snATAC on vessel-enriched ROSMAP cortex | Likely Synapse / controlled; **not** found as public GEO Series under “MultiVINE” as of 2026-08-18 | **Primary Track B path** — HOX / lineage peak enrichment in mural cells |
| Standard ROSMAP snATAC (PFC) | Chromatin, not vessel-enriched | Synapse DUC | Underpowered for rare mural peaks unless sorted |
| Wälchli / other FACS EC | scRNA mainly | GEO GSE256493 | RNA only — not Track B chromatin |

### Track B GO criteria (all required)

1. Vessel-enriched or FACS mural/endo nuclei with snATAC.  
2. Region annotation mappable to TerritoryVoxels (prefer enrichment-matched ACA vs MCA).  
3. Pre-specified HOX / neural-crest mural peak set passes FDR; RNA-only HOX remains exploratory.  
4. If chromatin null → **permanently retire Track B lead**; keep Track C.

### Outreach (ready)

See `docs/outreach/email_Yang_MultiVINE_snATAC_ask_2026-08-18.md`.

---

## Decision tree

```
Subject-level angio available?
  YES → run nested architecture models (Track C)
  NO  → keep atlas geography as descriptive only

snATAC mural lineage peaks support HOX/crest?
  YES → elevate Track B as secondary Nature claim
  NO  → Track B stays gated / closed; do not claim developmental endo identity
```

**TRCI/BCI:** not a mechanism arm for the Nature identity package.
