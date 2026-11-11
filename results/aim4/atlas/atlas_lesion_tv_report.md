# ATLAS lesion x arterial territories (ArterialAtlas_level2)

**Atlas map:** ArterialAtlas_level2 (L/R collapsed to ACA/MCA/PCA/VB)
**Rationale:** TerritoryVoxels_BMM labels 1-4 are effectively unilateral in this release; ArterialAtlas_level2 provides bilateral major territories per ArterialAtlasLabels.txt.

**n_ok / primary / strict:** 1452 / 1445 / 1444
- Median lesion voxels (primary): 7449.0
- Median assigned TV frac (primary): 0.9973753280839895
- Mean lesion share ACA/MCA/PCA/VB (primary): {'ACA': 0.09585793267727814, 'MCA': 0.6951560812113908, 'PCA': 0.09166610860713116, 'VB': 0.09731892338873774}
- Dominant territory counts (primary): {'MCA': 1114, 'VB': 139, 'PCA': 113, 'ACA': 78, 'None': 1}

## Locked interpretation
- Independent stroke outcome arm for Aim 4 (not Liu construction lesions).
- Use bilateral level2 major territories; do not use TerritoryVoxels_BMM for ATLAS.
- Geography is descriptive until nested with geometry (+/- molecular features).
- Do not claim molecular H3 from lesion shares alone.

Table: `results\aim4\atlas\atlas_lesion_tv_fractions.csv`
