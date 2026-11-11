# Yang 10x load + TerritoryVoxels attach

- Cells: **319,449**
- Genes: **33,538**
- Libraries: **25**
- Object: `data/processed/yang/yang_vascular.h5ad`

## Streams

analysis_stream  n_libs  n_cells
Aim3_borderzone       8   115991
     H1_primary      17   203458

## Per-library QC (medians)

geo_accession                  region disease analysis_stream  n_cells  median_genes  median_counts  median_pct_mt
   GSM4982083             hippocampus control      H1_primary    12142         708.0          958.0       0.235849
   GSM4982084             hippocampus      AD      H1_primary    18839         810.0         1170.0       0.285307
   GSM4982085             hippocampus control      H1_primary    11048         249.0          321.0       0.000000
   GSM4982086             hippocampus control      H1_primary    10039        1070.0         1787.0       0.458716
   GSM4982087             hippocampus control      H1_primary    19098         852.0         1342.0       0.184867
   GSM4982088             hippocampus control      H1_primary    23227         848.0         1220.0       0.300075
   GSM4982089             hippocampus      AD      H1_primary    10600         462.0          620.0       0.111359
   GSM4982090             hippocampus      AD      H1_primary     9101         361.0          460.0       0.311042
   GSM4982091             hippocampus control      H1_primary     9190         483.0          629.0       0.512821
   GSM4982092             hippocampus control      H1_primary    10929         852.0         1264.0       0.102881
   GSM4982093             hippocampus      AD      H1_primary    11435         529.0          679.0       0.331675
   GSM4982094             hippocampus      AD      H1_primary    10902         455.0          587.0       0.096829
   GSM4982095             hippocampus      AD      H1_primary     4893         321.0          425.0       0.210084
   GSM4982096             hippocampus control      H1_primary     8884         342.5          433.0       0.282619
   GSM4982097             hippocampus      AD      H1_primary     6279         455.0          632.0       0.000000
   GSM4982098             hippocampus      AD      H1_primary     1899         557.0          865.0       0.191939
   GSM4982099             hippocampus      AD      H1_primary    24953         920.0         1341.0       0.098039
   GSM4982100 superior_frontal_cortex control Aim3_borderzone    12992         626.0          801.0       0.346620
   GSM4982101 superior_frontal_cortex control Aim3_borderzone    11259         725.0         1006.0       0.525394
   GSM4982102 superior_frontal_cortex control Aim3_borderzone    15069         829.0         1147.0       0.085837
   GSM4982103 superior_frontal_cortex      AD Aim3_borderzone    17215        1184.0         1857.0       0.386972
   GSM4982104 superior_frontal_cortex      AD Aim3_borderzone    20144         968.0         1432.0       0.104603
   GSM4982105 superior_frontal_cortex      AD Aim3_borderzone    12105         688.0          922.0       0.476758
   GSM4982106 superior_frontal_cortex control Aim3_borderzone    14774         791.0         1135.0       0.124611
   GSM4982107 superior_frontal_cortex      AD Aim3_borderzone    12433         591.0          789.0       0.329308

Contrast lock: `territory_contrast` = PCA_hippocampus vs ACA_MCA_borderzone_SFG.
