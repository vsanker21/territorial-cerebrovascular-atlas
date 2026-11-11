# Walchli pathology strata - exploratory probe

**Verdict:** `EXPLORATORY_NOT_TERRITORY`

Loaded AVM sorted EC (double-gzip Seurat RDS). Adult temporal control present on disk but large; skipped this pass.

PATH=AVM; A-V EC clusters present. Disease tissue, not ACA, not Class B territory proof.

                    file n_cells PATH n_EC_clusters
1 AVM_sorted_endothelial   20305  AVM            12
                                           top_clusters
1 Capillary;Vein;Angiogenic capillary;Venule;Large vein

           Capillary                 Vein Angiogenic capillary 
                5436                 3440                 2275 
              Venule           Large vein               Artery 
                2160                 1997                 1865 
        Large artery            Arteriole        Mitochondrial 
                1129                 1071                  587 
              EndoMT   Proliferating cell           Stem-to-EC 
                 231                   82                   32 

## Locked use
- H4-adjacent A-V / pathology precedent only
- Not enrichment-matched ACA Class B
