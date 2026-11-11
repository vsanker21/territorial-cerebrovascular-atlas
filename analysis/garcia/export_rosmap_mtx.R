# Convert Garcia ROSMAP vascular Seurat-ish RDS to Matrix Market + metadata
# Counts object may be dgCMatrix or data.frame from Seurat assay

args <- commandArgs(trailingOnly = TRUE)
in_counts <- args[[1]]
in_meta <- args[[2]]
out_dir <- args[[3]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
})

meta <- readRDS(in_meta)
counts <- readRDS(in_counts)

cat("counts class:", paste(class(counts), collapse=","), "\n")
cat("meta dim:", paste(dim(meta), collapse="x"), "\n")

# Coerce to dgCMatrix genes x cells
if (inherits(counts, "dgCMatrix") || inherits(counts, "dgTMatrix") || inherits(counts, "Matrix")) {
  mat <- as(counts, "dgCMatrix")
} else if (is.data.frame(counts) || is.matrix(counts)) {
  mat <- as(as.matrix(counts), "dgCMatrix")
} else if (is.list(counts) && !is.null(counts$counts)) {
  mat <- as(counts$counts, "dgCMatrix")
} else {
  stop("Unrecognized counts object: ", paste(class(counts), collapse=","))
}

# Align cells
cells <- colnames(mat)
if (is.null(cells)) stop("counts lack colnames")
# meta rownames or a cell id column
if (!is.null(rownames(meta)) && all(cells %in% rownames(meta))) {
  meta <- meta[cells, , drop = FALSE]
} else if ("id" %in% colnames(meta) && all(cells %in% as.character(meta$id))) {
  rownames(meta) <- as.character(meta$id)
  meta <- meta[cells, , drop = FALSE]
} else {
  # try intersection order
  rn <- if (!is.null(rownames(meta))) rownames(meta) else as.character(meta[[1]])
  common <- intersect(cells, rn)
  cat("intersect cells:", length(common), "/", length(cells), "\n")
  mat <- mat[, common, drop = FALSE]
  if (!is.null(rownames(meta))) meta <- meta[common, , drop = FALSE]
  else {
    rownames(meta) <- rn
    meta <- meta[common, , drop = FALSE]
  }
}

writeMM(mat, file.path(out_dir, "matrix.mtx"))
writeLines(rownames(mat), file.path(out_dir, "genes.tsv"))
writeLines(colnames(mat), file.path(out_dir, "barcodes.tsv"))
write.csv(meta, file.path(out_dir, "meta.csv"), row.names = TRUE)
cat("Wrote", ncol(mat), "cells x", nrow(mat), "genes to", out_dir, "\n")
