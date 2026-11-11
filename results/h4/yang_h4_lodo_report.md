# H4 donor-aware LODO — Yang tight endo

Leave-one-donor-out AUROC on **control** donors with both territorial arms;
classifier trained on all control both-arm cells then transferred per **AD** donor.

## Summary

| Metric | Value |
|---|---:|
| Control LODO AUROC mean | 0.472 |
| Control LODO AUROC std | 0.141 |
| AD transfer AUROC mean | 0.648 |
| AD transfer AUROC std | 0.083 |
| Delta (AD - control) | +0.177 |
| Control donors (both arms) | 4 |
| AD donors (both arms) | 4 |
| Features | 1500 |

## Per-fold

| fold_type | held_out_donor | n_train | n_test | auroc | accuracy |
|---|---|---:|---:|---:|---:|
| control_lodo | Control_1 | 4607 | 3393 | 0.627 | 0.563 |
| control_lodo | Control_2 | 6667 | 1333 | 0.358 | 0.429 |
| control_lodo | Control_3 | 6736 | 1264 | 0.596 | 0.566 |
| control_lodo | Control_4 | 5990 | 2010 | 0.307 | 0.352 |
| ad_transfer_by_donor | AD_1 | 8000 | 2708 | 0.735 | 0.679 |
| ad_transfer_by_donor | AD_2 | 8000 | 1036 | 0.525 | 0.474 |
| ad_transfer_by_donor | AD_3 | 8000 | 591 | 0.622 | 0.619 |
| ad_transfer_by_donor | AD_4 | 8000 | 3018 | 0.712 | 0.652 |

## Interpretation

Cross-donor LODO in **controls is near chance** (mean AUROC 0.47). The earlier
cell-level stratified CV (~0.87) was inflated by within-donor structure leaking
across random folds.

This does **not** support a strong shared cross-donor territorial classifier on Yang
alone. Use within-donor CV (`yang_h4_within_donor_report.md`) for the H4 degradation
test under Yang's paired design. Proper cross-donor H1/H4 needs Class B / multi-region
cohorts where territory is not aliased with region+donor sampling.
