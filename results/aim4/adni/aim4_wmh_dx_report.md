# ADNI Aim4 WMH geography — DIAGNOSIS strata

Primary QC n=155; DIAGNOSIS match rate=1.00

    stratum   n  mean_share_MCA  mean_share_ACA  frac_subjects_MCA_gt_ACA  median_TOTAL_WMH_ucd  median_n_wmh
all_primary 155        0.207277        0.070538                  0.451613               8.15900          87.0
         CN  79        0.235860        0.074538                  0.481013               7.10018          75.0
        MCI  67        0.162436        0.062194                  0.402985               9.21200          95.0
   Dementia   9        0.290195        0.097536                  0.555556              15.94500         109.0
      other   0             NaN             NaN                       NaN                   NaN           NaN
    unknown   0             NaN             NaN                       NaN                   NaN           NaN

Tests: {'MCA_share_CN_vs_MCI': {'U': 2997.0, 'p': 0.1517396202128571}}

Locked: mean MCA share > mean ACA share remains the primary geography claim;
report subject-level MCA>ACA fractions as secondary descriptive stats.
