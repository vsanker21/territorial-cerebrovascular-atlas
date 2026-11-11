# Aim 4 exploratory molecular × ADNI WMH nest

**Verdict:** `EXPLORATORY_ECOLOGICAL_n3__NOT_CONFIRMATORY_H3`

## Territory table
territory  n_donors_allen  n_cells_allen  adni_mean_wmh_share  mra_vessel_probability  mra_vessel_radius  mod_BBB_tight_junction  mod_angiogenesis  mod_ECM_basement  mod_hypoxia_stress  mod_inflammation
      ACA               5            165             0.070538                1.313527           0.239551                2.991409          3.679895          1.732916            1.868508          0.801701
      MCA               6            425             0.207277                0.499102           0.192529                3.322516          3.686326          1.530510            1.592489          0.574053
      PCA               5            102             0.039610                0.290947           0.177866                3.866130          3.905240          1.388604            1.485910          0.687042

## Spearman
               feature  spearman_rho_vs_wmh_share        p  n
mod_BBB_tight_junction                       -0.5 0.666667  3
      mod_angiogenesis                       -0.5 0.666667  3
      mod_ECM_basement                        0.5 0.666667  3
    mod_hypoxia_stress                        0.5 0.666667  3
      mod_inflammation                       -0.5 0.666667  3
mra_vessel_probability                        0.5 0.666667  3
     mra_vessel_radius                        0.5 0.666667  3

## Nested R² (descriptive; n=3)
- Geometry: 1.0
- Modules: 1.0
- Geometry+modules: 1.0

Locked: exploratory only; not confirmatory molecular H3.
