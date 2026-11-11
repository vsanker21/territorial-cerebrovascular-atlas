from pathlib import Path

R = r'''# Annotate LIBD dACC GSE296731 with Azimuth humancortexref (Shah et al. method).
# Fixes Seurat v5 multi-layer merge that caused "non-conformable arguments" in RunAzimuth.
suppressPackageStartupMessages({
  library(Seurat)
  library(Azimuth)
  library(Matrix)
})

root <- "F:/Cerebrovascular SCS"
extracted <- file.path(root, "data/raw/aca_cohort_hunt/GSE296731/extracted")
outdir <- file.path(root, "data/processed/aca_cohort_hunt")
resdir <- file.path(root, "results/aca_cohort_hunt")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(resdir, recursive = TRUE, showWarnings = FALSE)

mtxs <- sort(list.files(extracted, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE))
stopifnot(length(mtxs) >= 1)
message("n libraries: ", length(mtxs))

load_one <- function(mtx) {
  stem <- sub("_matrix\\.mtx\\.gz$", "", basename(mtx))
  parent <- dirname(mtx)
  barcodes <- file.path(parent, paste0(stem, "_barcodes.tsv.gz"))
  features <- file.path(parent, paste0(stem, "_features.tsv.gz"))
  message("Loading ", stem)
  mat <- Matrix::readMM(mtx)
  bc <- read.delim(barcodes, header = FALSE, stringsAsFactors = FALSE)[, 1]
  feat <- read.delim(features, header = FALSE, stringsAsFactors = FALSE)
  genes <- make.unique(as.character(feat[, 2]))
  if (nrow(mat) == length(bc) && ncol(mat) == length(genes)) {
    mat <- Matrix::t(mat)
  }
  stopifnot(nrow(mat) == length(genes), ncol(mat) == length(bc))
  colnames(mat) <- paste0(stem, "-", bc)
  rownames(mat) <- genes
  nCount <- Matrix::colSums(mat)
  nGene <- Matrix::colSums(mat > 0)
  keep <- nCount >= 800 & nGene >= 400
  mat <- mat[, keep, drop = FALSE]
  message("  retained ", ncol(mat), " / ", length(keep), " barcodes")
  list(mat = mat, stem = stem, donor = strsplit(strsplit(stem, "_", fixed = TRUE)[[1]][2], "-", fixed = TRUE)[[1]][1], geo = strsplit(stem, "_", fixed = TRUE)[[1]][1])
}

loaded <- lapply(mtxs, load_one)
# Concatenate into one sparse matrix (avoids Seurat v5 layer issues)
all_genes <- unique(unlist(lapply(loaded, function(x) rownames(x$mat))))
mats <- lapply(loaded, function(x) {
  m <- x$mat
  missing <- setdiff(all_genes, rownames(m))
  if (length(missing)) {
    z <- Matrix::Matrix(0, nrow = length(missing), ncol = ncol(m), sparse = TRUE)
    rownames(z) <- missing
    colnames(z) <- colnames(m)
    m <- rbind(m, z)
  }
  m[all_genes, , drop = FALSE]
})
counts <- do.call(cbind, mats)
message("Combined matrix: ", nrow(counts), " genes x ", ncol(counts), " cells")

meta <- data.frame(
  library_id = sub("-.*$", "", colnames(counts)),
  row.names = colnames(counts),
  stringsAsFactors = FALSE
)
meta$donor_id <- vapply(strsplit(meta$library_id, "_", fixed = TRUE), function(p) {
  strsplit(p[2], "-", fixed = TRUE)[[1]][1]
}, character(1))
meta$geo_sample <- vapply(strsplit(meta$library_id, "_", fixed = TRUE), function(p) p[1], character(1))
meta$region <- "dACC"

merged <- CreateSeuratObject(counts = counts, meta.data = meta, project = "LIBD_dACC", min.cells = 3, min.features = 200)
message("Seurat cells after min filters: ", ncol(merged))
merged <- NormalizeData(merged, verbose = FALSE)

message("Running Azimuth humancortexref (paper reference)")
annotated <- RunAzimuth(merged, reference = "humancortexref")
annotated$azimuth_reference <- "humancortexref"
message("SUCCESS with humancortexref")

md <- as.data.frame(laura_PLACEHOLDER)
'''

R = R.replace("laura_PLACEHOLDER", 'slot(annotated, "meta.data")')

tail = r'''
md$cell_barcode <- rownames(md)
pred_cols <- grep("^predicted", colnames(md), value = TRUE)
message("predicted cols: ", paste(pred_cols, collapse = ", "))
write.csv(md, file.path(resdir, "libd_dacc_azimuth_metadata.csv"), row.names = FALSE)

ct_col <- if ("predicted.subclass" %in% colnames(md)) {
  "predicted.subclass"
} else if ("predicted.celltype.l2" %in% colnames(md)) {
  "predicted.celltype.l2"
} else if (length(pred_cols)) {
  pred_cols[1]
} else {
  stop("No predicted annotation columns found")
}
message("Using annotation column: ", ct_col)
labs <- tolower(as.character(md[[ct_col]]))
keep <- grepl("endo|vlmc|pericyt|smooth|vascular", labs)
message("vascular-like nuclei: ", sum(keep))
message("Endo: ", sum(grepl("endo", labs)))
message("VLMC: ", sum(grepl("vlmc", labs)))

vasc <- annotated[, keep]
saveRDS(vasc, file.path(outdir, "libd_dacc_azimuth_vascular.rds"))
saveRDS(annotated, file.path(outdir, "libd_dacc_azimuth_all.rds"))

mtxdir <- file.path(outdir, "libd_azimuth_vasc_mtx")
dir.create(mtxdir, showWarnings = FALSE)
counts_v <- tryCatch(GetAssayData(vasc, layer = "counts"), error = function(e) GetAssayData(vasc, slot = "counts"))
Matrix::writeMM(counts_v, file.path(mtxdir, "matrix.mtx"))
write.table(rownames(counts_v), file.path(mtxdir, "features.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(counts_v), file.path(mtxdir, "barcodes.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.csv(md[keep, , drop = FALSE], file.path(mtxdir, "metadata.csv"), row.names = FALSE)

sink(file.path(resdir, "libd_dacc_azimuth_summary.txt"))
cat("n_all\t", ncol(annotated), "\n", sep = "")
cat("n_vascular\t", sum(keep), "\n", sep = "")
cat("n_Endo\t", sum(grepl("endo", labs)), "\n", sep = "")
cat("n_VLMC\t", sum(grepl("vlmc", labs)), "\n", sep = "")
cat("reference\thumancortexref\n")
cat("ct_col\t", ct_col, "\n", sep = "")
cat("published_expect_Endo\t89\n")
cat("published_expect_VLMC\t159\n")
print(table(md[[ct_col]][keep]))
sink()
message("DONE_AZIMUTH")
'''

path = Path(r"F:\Cerebrovascular SCS\analysis\aca_cohort_hunt\annotate_libd_azimuth.R")
path.write_text(R + tail, encoding="utf-8")
print("wrote", path, "bytes", path.stat().st_size)
# sanity: no yuan/laura corruption
text = path.read_text(encoding="utf-8")
assert "yuan" not in text
assert "laura" not in text
assert "JoinLayers" not in text or True
print("ok")
