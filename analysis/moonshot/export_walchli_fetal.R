suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

rds <- "F:/Cerebrovascular SCS/data/raw/walchli/geo_GSE256493/suppl/_decomp/fetal_cns_sorted.rds"
out_dir <- "F:/Cerebrovascular SCS/data/processed/walchli/fetal_mtx"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(rds)
message("class: ", paste(class(obj), collapse = ","))
message("dim: ", paste(dim(obj), collapse = " x "))

meta <- slot(obj, "meta.data")
counts <- tryCatch(
  GetAssayData(obj, slot = "counts"),
  error = function(e) GetAssayData(obj, layer = "counts")
)

Matrix::writeMM(counts, file.path(out_dir, "matrix.mtx"))
write.table(rownames(counts), file.path(out_dir, "features.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(counts), file.path(out_dir, "barcodes.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.csv(meta, file.path(out_dir, "meta.csv"))
message("WROTE fetal mtx to ", out_dir)

# Adult temporal sorted for alignment target
adult_rds <- "F:/Cerebrovascular SCS/data/raw/walchli/geo_GSE256493/suppl/_decomp/adult_temporal_sorted.rds"
adult_dir <- "F:/Cerebrovascular SCS/data/processed/walchli/adult_temporal_mtx"
if (file.exists(adult_rds)) {
  dir.create(adult_dir, showWarnings = FALSE, recursive = TRUE)
  aobj <- readRDS(adult_rds)
  message("adult dim: ", paste(dim(aobj), collapse = " x "))
  ameta <- slot(aobj, "meta.data")
  acounts <- tryCatch(
    GetAssayData(aobj, slot = "counts"),
    error = function(e) GetAssayData(aobj, layer = "counts")
  )
  Matrix::writeMM(acounts, file.path(adult_dir, "matrix.mtx"))
  write.table(rownames(acounts), file.path(adult_dir, "features.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(colnames(acounts), file.path(adult_dir, "barcodes.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.csv(ameta, file.path(adult_dir, "meta.csv"))
  message("WROTE adult temporal mtx")
}
