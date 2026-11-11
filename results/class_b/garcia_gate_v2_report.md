# Garcia endo gate v2b (Garcia-calibrated)

## Why not Yang-strict v2?
Yang gate v2 (score_endo>=0.5, margin>=0.4, parenchyma<=0, LINGO1/RBFOX1/KCNIP4 caps)
retained only **~0.5%** of Garcia temporal soft endo (~60 cells) — unusable for DE.
Garcia soft "endo" is heavily neuron-contaminated (~72% LINGO1+).

## Gate v2b rules (primary)
- Canonical endo marker > 0: ['CLDN5', 'PECAM1', 'FLT1', 'VWF', 'CDH5']
- score_endo >= 0.25; margin vs other vascular >= 0.2
- Soft parenchyma caps: parenchyma<=0.5, neuron<=0.2
- Hard contaminant expression caps: {'LINGO1': 0.5, 'RBFOX1': 0.5, 'KCNIP4': 0.5, 'SNAP25': 0.5}
- Neuron module score uses classic panel only (['SNAP25', 'SLC17A7', 'GAD1', 'RBFOX3']); LINGO1/RBFOX1/KCNIP4 are hard caps

## Temporal ex vivo MCA
- Soft endo: **11,798** -> tight v2b: **207** (1.8%)
- Soft LINGO1>0.5: 72.2% -> tight: 0.0%
- Object: `data/processed/garcia/garcia_temporal_mca_qc_tight.h5ad`

## ROSMAP MCA isocortex (Midtemporal+Angular)
- Soft endo: **559** -> tight v2b: **112** (20.0%)
- Object: `data/processed/garcia/garcia_rosmap_vascular_tight.h5ad`
