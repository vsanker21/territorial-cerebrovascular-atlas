#!/usr/bin/env Rscript
# B0 gate: extract mural cell counts + region labels from Sun 2023 Seurat RDS
rds_path <- "F:/Cerebrovascular SCS/data/raw/sun2023_scADbbb/ROSMAP.VascularCells.counts_full.rds"
out_dir <- "F:/Cerebrovascular SCS/data/processed/sun2023"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading", rds_path, "...\n")
obj <- readRDS(rds_path)
cat("Class:", paste(class(obj), collapse = ", "), "\n")

meta <- NULL
if (inherits(obj, "Seurat")) {
  cat("Seurat cells:", ncol(obj), " genes:", nrow(obj), "\n")
  meta <- as.data.frame(obj[[]])
  meta$cell_barcode <- colnames(obj)
} else if (is.data.frame(obj)) {
  meta <- obj
  if (is.null(meta$cell_barcode)) meta$cell_barcode <- rownames(obj)
} else if (is.list(obj)) {
  cat("List names:", paste(names(obj), collapse = ", "), "\n")
  for (nm in names(obj)) {
    cat("  ", nm, "->", paste(class(obj[[nm]]), collapse = ","),
        " dim=", paste(dim(obj[[nm]]), collapse = "x"), "\n")
  }
  if (!is.null(obj$meta.data)) {
    meta <- as.data.frame(obj$meta.data)
    meta$cell_barcode <- rownames(meta)
  } else if (!is.null(obj$metadata)) {
    meta <- as.data.frame(obj$metadata)
    meta$cell_barcode <- rownames(meta)
  }
} else {
  print(str(obj, max.level = 2))
}

if (is.null(meta)) stop("Could not extract metadata from RDS")

write.csv(meta, file.path(out_dir, "sun2023_cell_metadata.csv"), row.names = FALSE)
cat("Wrote metadata n=", nrow(meta), " cols=", ncol(meta), "\n")
cat("Columns:\n")
print(colnames(meta))
for (c in colnames(meta)) {
  u <- unique(as.character(meta[[c]]))
  if (length(u) <= 50) {
    cat(sprintf("  %s (%d): %s\n", c, length(u), paste(head(u, 25), collapse = " | ")))
  } else {
    cat(sprintf("  %s (%d unique)\n", c, length(u)))
  }
}
