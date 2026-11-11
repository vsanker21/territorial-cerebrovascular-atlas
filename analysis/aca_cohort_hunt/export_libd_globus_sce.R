# Export published LIBD sce_Final Azimuth Endo/VLMC to MTX + metadata CSV
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

root <- "F:/Cerebrovascular SCS"
sce_path <- file.path(root, "data/raw/aca_cohort_hunt/libd_annotations/sce_Final.Rdata")
outdir <- file.path(root, "data/processed/aca_cohort_hunt/libd_globus_vasc_mtx")
resdir <- file.path(root, "results/aca_cohort_hunt")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

e <- new.env()
load(sce_path, envir = e)
sce <- e$sce
stopifnot(inherits(sce, "SingleCellExperiment"))
ct <- as.character(colData(sce)$cellType_azimuth)
keep <- ct %in% c("Endo", "VLMC")
message("n_all=", ncol(sce), " n_Endo=", sum(ct == "Endo"), " n_VLMC=", sum(ct == "VLMC"))
vasc <- sce[, keep]

rn <- rownames(vasc)
if ("gene_name" %in% colnames(rowData(vasc))) {
  sym <- as.character(rowData(vasc)$gene_name)
  ok <- !is.na(sym) & sym != ""
  vasc <- vasc[ok, ]
  rownames(vasc) <- make.unique(sym[ok])
}

counts <- counts(vasc)
Matrix::writeMM(counts, file.path(outdir, "matrix.mtx"))
write.table(rownames(counts), file.path(outdir, "features.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(counts), file.path(outdir, "barcodes.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)

md <- as.data.frame(colData(vasc))
md$cell_barcode <- colnames(vasc)
md$cell_type <- as.character(md$cellType_azimuth)
if ("brain" %in% colnames(md)) {
  md$donor_id <- as.character(md$brain)
} else if ("Sample" %in% colnames(md)) {
  md$donor_id <- as.character(md$Sample)
}
write.csv(md, file.path(outdir, "metadata.csv"), row.names = FALSE)

sink(file.path(resdir, "libd_dacc_globus_sce_summary.txt"))
cat("source\tsce_Final.Rdata (jhpce#dACC_spatial)\n")
cat("n_all\t", ncol(sce), "\n", sep = "")
cat("n_Endo\t", sum(ct == "Endo"), "\n", sep = "")
cat("n_VLMC\t", sum(ct == "VLMC"), "\n", sep = "")
cat("annotation\tcellType_azimuth\n")
print(table(md$cell_type))
print(table(md$donor_id, md$cell_type))
sink()
message("DONE_EXPORT_GLOBUS")
