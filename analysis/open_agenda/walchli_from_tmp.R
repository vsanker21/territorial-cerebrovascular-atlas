suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("F:/Cerebrovascular SCS/data/raw/walchli/geo_GSE256493/suppl/_tmp_avm_sorted.rds")
md <- slot(obj, "meta.data")
df <- data.frame(
  file = "AVM_sorted_endothelial",
  n_cells = ncol(obj),
  PATH = paste(unique(as.character(md$PATH)), collapse = "|"),
  n_EC_clusters = length(unique(as.character(md$ECclusters))),
  top_clusters = paste(head(names(sort(table(as.character(md$ECclusters)), TRUE)), 5), collapse = ";"),
  stringsAsFactors = FALSE
)
dir.create("F:/Cerebrovascular SCS/results/h4", FALSE, TRUE)
write.csv(df, "F:/Cerebrovascular SCS/results/h4/walchli_pathology_probe.csv", row.names = FALSE)
sink("F:/Cerebrovascular SCS/results/h4/walchli_pathology_exploratory_report.md")
cat("# Walchli pathology strata - exploratory probe\n\n")
cat("**Verdict:** `EXPLORATORY_NOT_TERRITORY`\n\n")
cat("Loaded AVM sorted EC (double-gzip Seurat RDS). Adult temporal control present on disk but large; skipped this pass.\n\n")
cat("PATH=AVM; A-V EC clusters present. Disease tissue, not ACA, not Class B territory proof.\n\n")
print(df)
print(sort(table(as.character(md$ECclusters)), decreasing = TRUE))
cat("\n## Locked use\n- H4-adjacent A-V / pathology precedent only\n- Not enrichment-matched ACA Class B\n")
sink()
message("DONE_WALCHLI")
