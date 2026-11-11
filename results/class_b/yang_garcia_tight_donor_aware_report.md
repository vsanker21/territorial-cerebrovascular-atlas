# Donor-aware Class B — tight Harmony object (out-of-fold)

Donors recovered from barcodes (Yang GSM; Garcia D20-*; ROSMAP region_lib).
Scores are **leave-one-donor-out** (not in-sample).

| Metric | Value |
|---|---:|
| Cell-level CV AUROC | 0.682 +/- 0.037 |
| Donor-level OOF AUROC | 0.067 |
| N Yang / Garcia donors scored | 4 / 15 |
| Multiclass within-donor LODO folds | 0 |

Donor-level OOF AUROC tests whether a held-out donor's mean score ranks Yang vs Garcia correctly.
Near-chance = little transferable territorial signal beyond batch/donor structure.

## Per-fold / donor table

      fold_type             held_out    auroc  n_test  mean_garcia_score         batch  true_label                                note
        cell_cv                fold0 0.700703     264                NaN           NaN         NaN                                 NaN
        cell_cv                fold1 0.731875     264                NaN           NaN         NaN                                 NaN
        cell_cv                fold2 0.689453     264                NaN           NaN         NaN                                 NaN
        cell_cv                fold3 0.666719     264                NaN           NaN         NaN                                 NaN
        cell_cv                fold4 0.619921     263                NaN           NaN         NaN                                 NaN
leave_one_donor           GSM4982100      NaN     429           0.612880          Yang         0.0 single-class donor (OOF score only)
leave_one_donor           GSM4982106      NaN     221           0.589273          Yang         0.0 single-class donor (OOF score only)
leave_one_donor           GSM4982101      NaN     218           0.590925          Yang         0.0 single-class donor (OOF score only)
leave_one_donor           GSM4982102      NaN     132           0.629918          Yang         0.0 single-class donor (OOF score only)
leave_one_donor           D20-169016      NaN      74           0.273200 Garcia_exvivo         1.0 single-class donor (OOF score only)
leave_one_donor           D20-125067      NaN      51           0.470375 Garcia_exvivo         1.0 single-class donor (OOF score only)
leave_one_donor           D20-125071      NaN      19           0.503380 Garcia_exvivo         1.0 single-class donor (OOF score only)
leave_one_donor           D20-125065      NaN      12           0.374911 Garcia_exvivo         1.0 single-class donor (OOF score only)
leave_one_donor           D20-125070      NaN      11           0.463012 Garcia_exvivo         1.0 single-class donor (OOF score only)
leave_one_donor MidtemporalCortex_17      NaN      11           0.525125 Garcia_ROSMAP         1.0 single-class donor (OOF score only)
leave_one_donor           D20-125073      NaN      10           0.543600 Garcia_exvivo         1.0 single-class donor (OOF score only)
leave_one_donor           D20-125063      NaN      10           0.400133 Garcia_exvivo         1.0 single-class donor (OOF score only)
leave_one_donor      AngularGyrus_32      NaN       9           0.507582 Garcia_ROSMAP         1.0 single-class donor (OOF score only)
leave_one_donor           D20-125069      NaN       9           0.394438 Garcia_exvivo         1.0 single-class donor (OOF score only)
leave_one_donor      AngularGyrus_24      NaN       9           0.532858 Garcia_ROSMAP         1.0 single-class donor (OOF score only)
leave_one_donor      AngularGyrus_33      NaN       8           0.653470 Garcia_ROSMAP         1.0 single-class donor (OOF score only)
leave_one_donor      AngularGyrus_30      NaN       8           0.450992 Garcia_ROSMAP         1.0 single-class donor (OOF score only)
leave_one_donor      AngularGyrus_48      NaN       7           0.506785 Garcia_ROSMAP         1.0 single-class donor (OOF score only)
leave_one_donor MidtemporalCortex_38      NaN       5           0.316471 Garcia_ROSMAP         1.0 single-class donor (OOF score only)
donor_oof_score           GSM4982100      NaN     429           0.612880          Yang         0.0                                 NaN
donor_oof_score           GSM4982106      NaN     221           0.589273          Yang         0.0                                 NaN
donor_oof_score           GSM4982101      NaN     218           0.590925          Yang         0.0                                 NaN
donor_oof_score           GSM4982102      NaN     132           0.629918          Yang         0.0                                 NaN
donor_oof_score           D20-169016      NaN      74           0.273200 Garcia_exvivo         1.0                                 NaN
donor_oof_score           D20-125067      NaN      51           0.470375 Garcia_exvivo         1.0                                 NaN
donor_oof_score           D20-125071      NaN      19           0.503380 Garcia_exvivo         1.0                                 NaN
donor_oof_score           D20-125065      NaN      12           0.374911 Garcia_exvivo         1.0                                 NaN
donor_oof_score           D20-125070      NaN      11           0.463012 Garcia_exvivo         1.0                                 NaN
donor_oof_score MidtemporalCortex_17      NaN      11           0.525125 Garcia_ROSMAP         1.0                                 NaN
donor_oof_score           D20-125073      NaN      10           0.543600 Garcia_exvivo         1.0                                 NaN
donor_oof_score           D20-125063      NaN      10           0.400133 Garcia_exvivo         1.0                                 NaN
donor_oof_score      AngularGyrus_32      NaN       9           0.507582 Garcia_ROSMAP         1.0                                 NaN
donor_oof_score           D20-125069      NaN       9           0.394438 Garcia_exvivo         1.0                                 NaN
donor_oof_score      AngularGyrus_24      NaN       9           0.532858 Garcia_ROSMAP         1.0                                 NaN
donor_oof_score      AngularGyrus_33      NaN       8           0.653470 Garcia_ROSMAP         1.0                                 NaN
donor_oof_score      AngularGyrus_30      NaN       8           0.450992 Garcia_ROSMAP         1.0                                 NaN
donor_oof_score      AngularGyrus_48      NaN       7           0.506785 Garcia_ROSMAP         1.0                                 NaN
donor_oof_score MidtemporalCortex_38      NaN       5           0.316471 Garcia_ROSMAP         1.0                                 NaN
