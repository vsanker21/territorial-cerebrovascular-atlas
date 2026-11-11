# Aim 4 exploratory molecular × ATLAS nest

**Verdict:** `EXPLORATORY_ECOLOGICAL_n3__NOT_CONFIRMATORY_H3`

## Design
- Molecular: Allen multi-region Endo module means (ACA/MCA/PCA; VB absent).
- Modules: BBB / angiogenesis / ECM / hypoxia / inflammation (Class B panel).
- Outcome: ATLAS primary mean lesion share by territory.
- Null: MRA vessel probability + radius.

## Territory table
territory  n_donors_allen  n_cells_allen  atlas_mean_lesion_share  atlas_mean_lesion_burden  mra_vessel_probability  mra_vessel_radius  mod_BBB_tight_junction  mod_angiogenesis  mod_ECM_basement  mod_hypoxia_stress  mod_inflammation
      ACA               5            165                 0.095858                  0.011028                1.313527           0.239551                2.991409          3.679895          1.732916            1.868508          0.801701
      MCA               6            425                 0.695156                  0.032676                0.499102           0.192529                3.322516          3.686326          1.530510            1.592489          0.574053
      PCA               5            102                 0.091666                  0.006737                0.290947           0.177866                3.866130          3.905240          1.388604            1.485910          0.687042

## Spearman (module/geometry vs ATLAS lesion share)
               feature  spearman_rho_vs_lesion_share        p  n
mod_BBB_tight_junction                          -0.5 0.666667  3
      mod_angiogenesis                          -0.5 0.666667  3
      mod_ECM_basement                           0.5 0.666667  3
    mod_hypoxia_stress                           0.5 0.666667  3
      mod_inflammation                          -0.5 0.666667  3
mra_vessel_probability                           0.5 0.666667  3
     mra_vessel_radius                           0.5 0.666667  3

## Nested R² (descriptive only; n=3)
- Geometry only: 1.0
- Modules only: 1.0
- Geometry + modules: 1.0

## Locked interpretation
- This is an **exploratory ecological** screen, **not** confirmatory molecular H3.
- n_territories=3 precludes stable incremental inference; VB molecular arm missing.
- Enrichment-matched vessel-enriched ACC still required for confirmatory Class B / H3.
- Imaging layer (ATLAS geography + MRA null) remains frozen independently.

Tables: `results\aim4\atlas\aim4_molecular_atlas_ecological_table.csv`, `aim4_molecular_atlas_spearman.csv`
