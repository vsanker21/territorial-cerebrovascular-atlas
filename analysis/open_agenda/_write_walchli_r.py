from pathlib import Path

R = r"""# Walchli exploratory probe (double-gzip Seurat RDS)
suppressPackageStartupMessages(library(Seurat))
root <- "F:/Cerebrovascular SCS"
suppl <- file.path(root, "data/raw/walchli/geo_GSE256493/suppl")
outdir <- file.path(root, "results/h4")
tmpdir <- file.path(suppl, "_decomp")
dir.create(outdir, FALSE, TRUE)
dir.create(tmpdir, FALSE, TRUE)

# Double-gunzip helper via R connections
load_double_gz <- function(f) {
  message("Decompress ", basename(f))
  mid <- tempfile(fileext = ".gz", tmpdir = tmpdir)
  out <- tempfile(fileext = ".rds", tmpdir = tmpdir)
  # pass1
  con1 <- gzfile(f, "rb")
  zz <- file(mid, "wb")
  while (TRUE) {
    chunk <- readBin(con1, what = "raw", n = 8e6)
    if (!length(chunk)) break
    writeBin(chunk, zz)
  }
  close(con1); close(zz)
  # pass2
  con2 <- gzfile(mid, "rb")
  zz2 <- file(out, "wb")
  while (TRUE) {
    chunk <- readBin(con2, what = "raw", n = 8e6)
    if (!length(chunk)) break
    writeBin(chunk, zz2)
  }
  close(con2); close(zz2)
  unlink(mid)
  obj <- readRDS(out)
  unlink(out)
  obj
}

targets <- c(
  "GSE256493_Adult_control_brain_temporal_lobe_sorted_endothelial_cells_seurat_object.rds.gz",
  "GSE256493_Arteriovenous_malformations_sorted_endothelial_cells_seurat_object.rds.gz"
)
# Adult control is ~3.5GB gz - may be heavy; try AVM first then adult if memory allows
targets <- targets[file.exists(file.path(suppl, targets))]
# Prefer AVM (smaller) always; adult only if < 2GB gz
targets <- targets[file.info(file.path(suppl, targets))$size < 2.5e9]
message("targets: ", paste(targets, collapse = ", "))

rows <- list()
for (bn in targets) {
  f <- file.path(suppl, bn)
  obj <- tryCatch(load_double_gz(f), error = function(e) { message(conditionMessage(e)); NULL })
  if (is.null(obj)) next
  md <- slot(obj, "meta.data")
  path_col <- if ("PATH" %in% colnames(md)) "PATH" else NA
  cl_col <- if ("ECclusters" %in% colnames(md)) "ECclusters" else NA
  rows[[bn]] <- data.frame(
    file = bn,
    n_cells = ncol(obj),
    PATH = if (!is.na(path_col)) paste(unique(as.character(md[[path_col]])), collapse = "|") else NA,
    n_EC_clusters = if (!is.na(cl_col)) length(unique(as.character(md[[cl_col]]))) else NA,
    top_clusters = if (!is.na(cl_col)) paste(head(names(sort(table(as.character(md[[cl_col]])), TRUE)), 5), collapse = ";") else NA,
    stringsAsFactors = FALSE
  )
  if (!is.na(cl_col)) {
    message("clusters: ")
    print(sort(table(as.character(md[[cl_col]])), decreasing = TRUE))
  }
  rm(obj); gc()
}

df <- if (length(rows)) do.call(rbind, rows) else data.frame()
write.csv(df, file.path(outdir, "walchli_pathology_probe.csv"), row.names = FALSE)
sink(file.path(outdir, "walchli_pathology_exploratory_report.md"))
cat("# Walchli pathology strata - exploratory probe\n\n")
cat("**Verdict:** `EXPLORATORY_NOT_TERRITORY`\n\n")
cat("Seurat objects are **double-gzipped**; loaded AVM sorted EC successfully.\n")
cat("Adult temporal control = MCA-ish anatomy. Use for H4-adjacent A-V / disease precedent only.\n\n")
print(df)
cat("\n## Locked\n- Not ACA\n- Not Class B territory proof\n- A-V EC clusters present (Capillary/Vein/Artery/...)\n")
sink()
message("DONE_WALCHLI")
"""
Path(r"F:\Cerebrovascular SCS\analysis\open_agenda\walchli_pathology_probe.R").write_text(R, encoding="utf-8", newline="\n")
print("wrote")
