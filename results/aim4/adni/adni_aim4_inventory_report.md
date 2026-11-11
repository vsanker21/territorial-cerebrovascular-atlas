# ADNI Aim 4 inventory — existing data (no IDA WMH masks)

**Date:** 2026-07-24

## Locked path
- Spatial WMH masks **not** available on LONI for this account.
- Proceed with: UCD/ADSP **volume** outcomes + **local FLAIR→WMH segmentation** for TerritoryVoxels fractions.

## Counts
- 3T FLAIR series: 8486 (2472 subjects)
- Primary cohort (1 FLAIR/subject): 2472
- With UCD TOTAL_WMH: 2420
- With ADSP PHC harmonized WMH: 2349
- ADSL join: False

## On-disk imaging packs
- T1: `Clinical_T1w_Imaging_Cohort_MRI*.zip`
- FLAIR: `2D T2 FLAIR.zip`, `2D T2 FLAIR_dataset.zip`

## Outputs
- `results\aim4\adni\adni_aim4_cohort_inventory.csv`
- `results/aim4/adni/adni_aim4_wmh_volume_baseline.csv`
