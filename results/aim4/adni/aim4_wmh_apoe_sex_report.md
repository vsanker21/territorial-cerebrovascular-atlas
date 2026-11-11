# ADNI Aim4 WMH geography — sex / APOE strata

Primary QC n=155; APOE match rate=0.72
APOE source: `ADNIMERGE2/data/APOERES.rda (RID+GENOTYPE)`

           stratum   n  mean_share_MCA  mean_share_ACA  median_n_wmh  median_TOTAL_WMH_ucd
       all_primary 155        0.207277        0.070538          87.0              8.159000
            female  82        0.203680        0.081592          91.5              8.619005
              male  73        0.211318        0.058120          84.0              7.784100
   APOE_e4_carrier  42        0.210555        0.057258          65.0              6.877195
APOE_e4_noncarrier  69        0.209972        0.070499          92.0              8.998000

Tests: {'MCA_share_female_vs_male': {'U': 3002.5, 'p': 0.9731533461219665}, 'MCA_share_e4_vs_non': {'U': 1436.0, 'p': 0.9366262795599518}}

Locked: descriptive sensitivity only; does not change primary MCA-dominant geography lock.
