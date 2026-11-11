# ADNI WMH → TerritoryVoxels pilot (local segmentation)

**Verdict:** `PILOT_LOCAL_WMH__NO_IDA_MASKS`

- Method: N4 FLAIR + affine MNI + WM-prior z-score hyperintensity (exploratory)
- OK: 3 / 3

```
 n_wmh_voxels  wmh_vox_ACA  tv_vox_ACA  wmh_frac_of_wmh_ACA  wmh_burden_in_ACA  wmh_vox_MCA  tv_vox_MCA  wmh_frac_of_wmh_MCA  wmh_burden_in_MCA  wmh_vox_PCA  tv_vox_PCA  wmh_frac_of_wmh_PCA  wmh_burden_in_PCA  wmh_vox_VB  tv_vox_VB  wmh_frac_of_wmh_VB  wmh_burden_in_VB        ptid   RID  TOTAL_WMH_ucd                    series  z_thr    seg_thr  seg_thr_z  seg_thr_p  seg_n_wm status
           28            0       10047                  0.0                0.0            0       58027                  0.0                0.0            0       11252                  0.0                0.0           0      10725                 0.0               0.0  002_S_0413   413       2.863230         Sagittal_3D_FLAIR    1.5 150.068235 154.416353 150.068235     92505     ok
            0            0       10047                  0.0                0.0            0       58027                  0.0                0.0            0       11252                  0.0                0.0           0      10725                 0.0               0.0 002_S_10814 10814       0.996863 Sagittal_3D_FLAIR__MSV22_    1.5  92.375527  92.375527  93.220525     90356     ok
            0            0       10047                  0.0                0.0            0       58027                  0.0                0.0            0       11252                  0.0                0.0           0      10725                 0.0               0.0  002_S_1155  1155       1.101240         Sagittal_3D_FLAIR    1.5 134.748635 134.748635 139.648750     92391     ok
```

## Locked
- Not using IDA spatial WMH (unavailable).
- Territory fractions are exploratory until scaled + QC'd.
- UCD TOTAL_WMH remains the validated volume endpoint for full-cohort analyses.
