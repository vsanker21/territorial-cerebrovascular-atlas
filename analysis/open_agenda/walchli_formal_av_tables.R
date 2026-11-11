# Formalize Wälchli adult temporal control + pathology A–V exploratory tables
# Objects on GEO are double-gzipped (.rds.gz that remains gzip after one pass).
suppressPackageStartupMessages({
  library(Seurat)
})

root <- "F:/Cerebrovascular SCS"
suppl <- file.path(root, "data/raw/walchli/geo_GSE256493/suppl")
outdir <- file.path(root, "results/h4")
tmpdir <- file.path(suppl, "_decomp")
dir.create(outdir, FALSE, TRUE)
dir.create(tmpdir, FALSE, TRUE)

double_gunzip_to_rds <- function(f, out_rds) {
  message("Double-gunzip ", basename(f), " -> ", basename(out_rds))
  mid <- file.path(tmpdir, paste0(basename(f), ".pass1.gz"))
  if (file.exists(out_rds) && file.info(out_rds)$size > 1e6) {
    message("  reuse existing ", basename(out_rds))
    return(out_rds)
  }
  # pass 1
  con1 <- gzfile(f, "rb")
  zz <- file(mid, "wb")
  repeat {
    chunk <- readBin(con1, what = "raw", n = 16e6)
    if (!length(chunk)) break
    writeBin(chunk, zz)
  }
  close(con1); close(zz)
  magic <- readBin(file(mid, "rb"), what = "raw", n = 2); close(file(mid, "rb"))
  message("  after pass1 magic: ", paste(as.character(magic), collapse = " "))
  # pass 2
  con2 <- gzfile(mid, "rb")
  zz2 <- file(out_rds, "wb")
  repeat {
    chunk <- readBin(con2, what = "raw", n = 16e6)
    if (!length(chunk)) break
    writeBin(chunk, zz2)
  }
  close(con2); close(zz2)
  unlink(mid)
  message("  wrote ", basename(out_rds), " (", round(file.info(out_rds)$size / 1e9, 2), " GB)")
  out_rds
}

summarize_obj <- function(obj, label) {
  md <- slot(obj, "meta.data")
  cl_col <- if ("ECclusters" %in% colnames(md)) "ECclusters" else if ("ec_cluster" %in% colnames(md)) "ec_cluster" else NA
  path_col <- if ("PATH" %in% colnames(md)) "PATH" else if ("pathology" %in% colnames(md)) "pathology" else NA
  donor_cols <- intersect(c("orig.ident", "sample", "Sample", "donor", "Donor", "patient", "Patient"), colnames(md))
  cluster_tab <- if (!is.na(cl_col)) as.data.frame(sort(table(as.character(md[[cl_col]])), decreasing = TRUE)) else data.frame()
  if (nrow(cluster_tab)) colnames(cluster_tab) <- c("cluster", "n")
  # A–V coarse bins from Wälchli naming if present
  av <- NA_character_
  if (!is.na(cl_col)) {
    cl <- as.character(md[[cl_col]])
    av <- rep("other", length(cl))
    av[grepl("Arter|art", cl, ignore.case = TRUE)] <- "arterial"
    av[grepl("Vein|Venous|ven", cl, ignore.case = TRUE)] <- "venous"
    av[grepl("Capill", cl, ignore.case = TRUE)] <- "capillary"
    av[grepl("Tip|tip|prolif|Angiogen", cl, ignore.case = TRUE)] <- "tip_angiogenic"
  }
  list(
    label = label,
    n_cells = ncol(obj),
    n_genes = nrow(obj),
    path_values = if (!is.na(path_col)) paste(unique(as.character(md[[path_col]])), collapse = "|") else NA,
    cluster_col = cl_col,
    cluster_tab = cluster_tab,
    av_tab = if (!all(is.na(av))) as.data.frame(sort(table(av), decreasing = TRUE)) else data.frame(),
    donor_cols = paste(donor_cols, collapse = ","),
    n_donors_proxy = if (length(donor_cols)) length(unique(as.character(md[[donor_cols[1]]]))) else NA,
    meta_cols = paste(colnames(md), collapse = ";")
  )
}

targets <- list(
  adult_temporal_sorted = "GSE256493_Adult_control_brain_temporal_lobe_sorted_endothelial_cells_seurat_object.rds.gz",
  avm_sorted = "GSE256493_Arteriovenous_malformations_sorted_endothelial_cells_seurat_object.rds.gz",
  lgg_sorted = "GSE256493_Lower_grade_glioma_sorted_endothelial_cells_seurat_object.rds.gz"
)

rows <- list()
av_long <- list()
cl_long <- list()

for (nm in names(targets)) {
  bn <- targets[[nm]]
  f <- file.path(suppl, bn)
  if (!file.exists(f)) {
    message("MISSING ", bn)
    next
  }
  out_rds <- file.path(tmpdir, paste0(nm, ".rds"))
  ok <- tryCatch({
    double_gunzip_to_rds(f, out_rds)
    TRUE
  }, error = function(e) {
    message("gunzip fail ", nm, ": ", conditionMessage(e))
    FALSE
  })
  if (!ok) next
  message("readRDS ", nm)
  obj <- tryCatch(readRDS(out_rds), error = function(e) {
    message("readRDS fail: ", conditionMessage(e)); NULL
  })
  if (is.null(obj)) next
  s <- summarize_obj(obj, nm)
  rows[[nm]] <- data.frame(
    cohort = nm,
    file = bn,
    n_cells = s$n_cells,
    n_genes = s$n_genes,
    path_values = s$path_values,
    cluster_col = s$cluster_col,
    n_donors_proxy = s$n_donors_proxy,
    donor_cols = s$donor_cols,
    stringsAsFactors = FALSE
  )
  if (nrow(s$cluster_tab)) {
    ct <- s$cluster_tab
    ct$cohort <- nm
    cl_long[[nm]] <- ct
  }
  if (nrow(s$av_tab)) {
    at <- s$av_tab
    colnames(at) <- c("av_bin", "n")
    at$cohort <- nm
    av_long[[nm]] <- at
  }
  # free memory aggressively; keep decompressed rds for reuse
  rm(obj); gc()
}

df <- if (length(rows)) do.call(rbind, rows) else data.frame()
write.csv(df, file.path(outdir, "walchli_av_pathology_cohort_table.csv"), row.names = FALSE)
if (length(cl_long)) {
  write.csv(do.call(rbind, cl_long), file.path(outdir, "walchli_ec_cluster_counts.csv"), row.names = FALSE)
}
if (length(av_long)) {
  write.csv(do.call(rbind, av_long), file.path(outdir, "walchli_av_bin_counts.csv"), row.names = FALSE)
}

# Formal report
sink(file.path(outdir, "walchli_av_pathology_formal_report.md"))
cat("# Wälchli adult temporal + pathology — formal exploratory A–V tables\n\n")
cat("**Verdict:** `EXPLORATORY_A-V_AND_PATHOLOGY__NOT_TERRITORY`\n\n")
cat("## Methods locks\n")
cat("- GEO GSE256493 / GSE256490 Seurat objects are **double-gzipped**; one `gunzip` leaves gzip magic `1f 8b`.\n")
cat("- Adult control = **temporal lobe** FACS CD31+/CD45− EC → **MCA anatomy**, not ACA.\n")
cat("- Use for A–V zonation / pathology EC-state precedent only; **not** Class B territory proof.\n")
cat("- Enrichment-matched relative to unsorted ACC cohorts, but wrong territory for ACA claims.\n\n")
cat("## Cohort table\n\n")
print(df)
cat("\n## EC cluster counts\n\n")
if (length(cl_long)) print(do.call(rbind, cl_long)) else cat("(none)\n")
cat("\n## Coarse A–V bins (name-based)\n\n")
if (length(av_long)) print(do.call(rbind, av_long)) else cat("(none — cluster names not matching Arter/Vein/Capill)\n")
cat("\n## What this enables\n")
cat("- Document adult temporal control A–V cluster composition vs AVM / LGG.\n")
cat("- Precedent that vessel-enriched temporal EC carry clear A–V structure (H4-adjacent / Aim2 zonation).\n")
cat("- Does **not** unlock enrichment-matched ACA↔MCA DE.\n")
sink()
message("DONE_WALCHLI_FORMAL")
