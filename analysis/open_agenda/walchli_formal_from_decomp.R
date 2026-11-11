# Formalize Wälchli from already double-gunzipped RDS (see walchli_double_gunzip.py)
suppressPackageStartupMessages(library(Seurat))

root <- "F:/Cerebrovascular SCS"
decomp <- file.path(root, "data/raw/walchli/geo_GSE256493/suppl/_decomp")
outdir <- file.path(root, "results/h4")
dir.create(outdir, FALSE, TRUE)

targets <- c(
  adult_temporal_sorted = "adult_temporal_sorted.rds",
  avm_sorted = "avm_sorted.rds",
  lgg_sorted = "lgg_sorted.rds"
)

av_bin <- function(cl) {
  av <- rep("other", length(cl))
  av[grepl("Arteriole|Large artery|Artery", cl, ignore.case = TRUE)] <- "arterial"
  av[grepl("Venule|Large vein|Vein", cl, ignore.case = TRUE)] <- "venous"
  av[grepl("Capillary", cl, ignore.case = TRUE) & !grepl("Angiogenic", cl, ignore.case = TRUE)] <- "capillary"
  av[grepl("Angiogenic", cl, ignore.case = TRUE)] <- "angiogenic_capillary"
  av[grepl("Prolif|Stem|EndoMT|Mitochon", cl, ignore.case = TRUE)] <- "activated_other"
  av
}

rows <- list(); cl_long <- list(); av_long <- list()

for (nm in names(targets)) {
  f <- file.path(decomp, targets[[nm]])
  message("Loading ", nm, " (", round(file.info(f)$size/1e9, 2), " GB)")
  obj <- readRDS(f)
  md <- slot(obj, "meta.data")
  cl_col <- if ("ECclusters" %in% colnames(md)) "ECclusters" else NA
  path_col <- if ("PATH" %in% colnames(md)) "PATH" else NA
  donor_cols <- intersect(c("orig.ident", "sample", "Sample", "donor", "Donor", "patient", "Patient", "ID"), colnames(md))
  n_don <- if (length(donor_cols)) length(unique(as.character(md[[donor_cols[1]]]))) else NA_integer_

  rows[[nm]] <- data.frame(
    cohort = nm,
    file = targets[[nm]],
    n_cells = ncol(obj),
    n_genes = nrow(obj),
    path_values = if (!is.na(path_col)) paste(sort(unique(as.character(md[[path_col]]))), collapse = "|") else NA_character_,
    cluster_col = cl_col,
    n_ec_clusters = if (!is.na(cl_col)) length(unique(as.character(md[[cl_col]]))) else NA_integer_,
    n_donors_proxy = n_don,
    donor_col = if (length(donor_cols)) donor_cols[1] else NA_character_,
    stringsAsFactors = FALSE
  )

  if (!is.na(cl_col)) {
    tab <- as.data.frame(sort(table(as.character(md[[cl_col]])), decreasing = TRUE), stringsAsFactors = FALSE)
    colnames(tab) <- c("cluster", "n")
    tab$cohort <- nm
    tab$frac <- tab$n / sum(tab$n)
    cl_long[[nm]] <- tab

    av <- av_bin(as.character(md[[cl_col]]))
    at <- as.data.frame(sort(table(av), decreasing = TRUE), stringsAsFactors = FALSE)
    colnames(at) <- c("av_bin", "n")
    at$cohort <- nm
    at$frac <- at$n / sum(at$n)
    av_long[[nm]] <- at
  }
  rm(obj); gc()
}

df <- do.call(rbind, rows)
write.csv(df, file.path(outdir, "walchli_av_pathology_cohort_table.csv"), row.names = FALSE)
cldf <- do.call(rbind, cl_long)
write.csv(cldf, file.path(outdir, "walchli_ec_cluster_counts.csv"), row.names = FALSE)
avdf <- do.call(rbind, av_long)
write.csv(avdf, file.path(outdir, "walchli_av_bin_counts.csv"), row.names = FALSE)

# Wide A–V composition for paper table
wide <- reshape(avdf[, c("cohort", "av_bin", "frac")], idvar = "cohort", timevar = "av_bin", direction = "wide")
write.csv(wide, file.path(outdir, "walchli_av_bin_fractions_wide.csv"), row.names = FALSE)

sink(file.path(outdir, "walchli_av_pathology_formal_report.md"))
cat("# Wälchli adult temporal + pathology — formal exploratory A–V tables\n\n")
cat("**Verdict:** `EXPLORATORY_A-V_AND_PATHOLOGY__NOT_TERRITORY`\n\n")
cat("## Methods locks\n")
cat("- GEO GSE256493 Seurat objects are **double-gzipped**; decompressed to `data/raw/walchli/geo_GSE256493/suppl/_decomp/`.\n")
cat("- Adult control = **temporal lobe** FACS CD31+/CD45− EC → **MCA anatomy**, not ACA.\n")
cat("- Use for A–V zonation / pathology EC-state precedent only; **not** Class B territory proof.\n")
cat("- Enrichment yes relative to unsorted ACC, but wrong territory for ACA marker claims.\n\n")
cat("## Cohort table\n\n```\n"); print(df); cat("```\n\n")
cat("## EC cluster counts\n\n```\n"); print(cldf); cat("```\n\n")
cat("## Coarse A–V bins\n\n```\n"); print(avdf); cat("```\n\n")
cat("## What this enables / does not\n")
cat("- Enables: adult temporal control A–V composition vs AVM / LGG (exploratory H4-adjacent / Aim2 zonation).\n")
cat("- Does **not** unlock enrichment-matched ACA↔MCA DE or territorial H1.\n")
sink()
message("DONE_WALCHLI_FORMAL")
