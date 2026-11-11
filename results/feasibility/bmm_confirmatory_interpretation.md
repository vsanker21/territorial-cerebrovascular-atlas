# NITRC BMM confirmatory interpretation (2026-07-23)

## What was ingested
From authenticated NITRC download (user `vsanker`):

- `ArterialAtlas.zip` (9.95 MB)
- `Atlas_MNI152.zip` (9.80 MB)
- `Vascular_Probabilistic_Maps.nii` (227 MB)
- Labels file

Primary analysis space: `Atlas_182_MNI152` (FSL 182×218×182, proper MNI affine).

## Critical methodological findings

### 1. `ProbArterialAtlas_BMM` is not a simplex soft assignment
Across voxels with any mass, mean ∑p ≈ **0.35**; almost no voxels have ∑p ≈ 1.  
Therefore **ROI-mean of BMM channels followed by renormalization is not a valid territorial probability** — it overweight sparse channels and can invert anatomy (e.g. Frontal Medial Cortex).

### 2. `TerritoryVoxels_BMM` label scheme (this file) ≠ Sci Data README text
Empirically vs hard atlas:

| Label | Predominant hard major | Interpretation |
|---:|---|---|
| 1 | ACA (98.8%) | ACA |
| 2 | MCA (89.8%) | MCA |
| 3 | PCA (96.3%) | PCA |
| 4 | VB (97.2%) | **VB** (README said “multi” — incorrect for this release) |
| 5 | mixed | ambiguous / multi-claim candidate |
| 6 | mixed / large | unassigned / low-confidence (majority of brain voxels) |

### 3. BorderZone_ProbAve
~91% NaN — use only where defined; do not average blindly.

## Confirmatory ROI results (TerritoryVoxels labels 1–4 only)

| ROI | Role | Among assigned | Fraction of ROI assigned | Verdict |
|---|---|---|---:|---|
| Hippocampus | Yang | **PCA 75% / MCA 25%** | 0.31 | PCA-dominant — GO |
| Superior frontal gyrus | Yang | **ACA 54% / MCA 46%** | 0.25 | **ACA–MCA borderzone** |
| Temporal pole | Garcia ex vivo | **MCA 100%** | 0.43 | MCA — Class B GO |
| Middle frontal gyrus | Class B | **MCA 100%** | 0.49 | MCA — Class B GO |
| Occipital pole | Class A | **PCA 93%** | 0.34 | PCA — GO |
| Frontal medial cortex | Class B ACA | ACA 100% | **0.009** | Too sparsely assigned in TerritoryVoxels |

Hard-atlas voxel fractions (prior probe) remain useful where TerritoryVoxels assignment fraction is low: SFG hard ≈ ACA 93%, but that **overstates purity** relative to BMM-derived TerritoryVoxels.

## Revised scientific framing

1. **Yang within-donor contrast is not clean PCA↔ACA.**  
   Hippocampus is PCA-dominant; superior frontal is an **ACA/MCA watershed / mixed** territory.  
   This makes Yang exceptionally informative for **Aim 3 (molecular borderzone)** and still useful for posterior vs anterior-circulation contrasts, but **not** a pure Class B territory-differing test.

2. **Class B (ACA↔MCA) should use:**  
   - MCA: Garcia temporal ex vivo and/or MFG  
   - ACA: hard-atlas–anchored medial frontal / confirmed ACA TerritoryVoxels voxels; Frontal Medial is anatomically ACA but poorly covered by TerritoryVoxels assignment

3. **Primary H1 samples:** prefer ROIs with high TerritoryVoxels assignment fraction and P_dom ≥ 0.50 among labels 1–4 (temporal, MFG, occipital, hippocampus).

4. **Do not use raw BMM ROI-means for confirmatory claims** without mass thresholds and explicit simplex normalization checks.

## Feasibility status after NITRC
- **Still GO**, with sharper design: Yang = PCA vs borderzone; Class B MCA locked by Garcia temporal + MFG.
- Pre-registration should be updated to this TerritoryVoxels-based rule before OSF lock.
