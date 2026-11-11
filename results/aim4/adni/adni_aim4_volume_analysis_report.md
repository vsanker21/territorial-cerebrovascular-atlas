# ADNI Aim 4 — WMH volume analysis (existing data; no IDA masks)

**Verdict:** `VOLUME_ENDPOINT_READY__TERRITORY_VIA_LOCAL_FLAIR`

## Cohort
- n visits/rows: 2420 (subjects≈2420)
- TOTAL_WMH median (IQR): 2.33 (0.76–5.83) cc

## Associations (descriptive)
- By sex: `{'1.0': 2.7408, '2.0': 2.031, 'mannwhitney_p': 8.899241079571421e-07}`
- Age Spearman: `None`

## Nested Aim 4 structure
1. Geometry null: `results/aim4/mra_geometry_by_territory.csv` (already on disk)
2. Non-stroke outcome (this table): subject-level UCD/ADSP WMH volumes
3. Territory-resolved WMH: local FLAIR segmentation pilot (`adni_wmh_territory_pilot_report.md`) — scale after QC

## Locked
- No IDA spatial WMH masks.
- Do not treat whole-brain TOTAL_WMH as territory-specific lesion probability.
