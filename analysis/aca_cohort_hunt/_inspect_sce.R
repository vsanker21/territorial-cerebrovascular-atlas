options(warn = 1)
cat("start\n")
sce_path <- "F:/Cerebrovascular SCS/data/raw/aca_cohort_hunt/libd_annotations/sce_Final.Rdata"
cat("size", file.info(sce_path)$size, "\n")
ok <- FALSE
tryCatch({
  x <- readRDS(sce_path)
  cat("readRDS class:", paste(class(x), collapse = "/"), "\n")
  assign("obj", x, envir = .GlobalEnv)
  ok <- TRUE
}, error = function(e) cat("readRDS fail:", conditionMessage(e), "\n"))
if (!ok) {
  tryCatch({
    e <- new.env()
    load(sce_path, envir = e)
    cat("load names:", paste(ls(e), collapse = ", "), "\n")
    for (nm in ls(e)) {
      x <- e[[nm]]
      cat(nm, "class", paste(class(x), collapse = "/"),
          "size", format(object.size(x), units = "MB"), "\n")
      assign("obj", x, envir = .GlobalEnv)
      assign("obj_name", nm, envir = .GlobalEnv)
    }
  }, error = function(e) cat("load fail:", conditionMessage(e), "\n"))
}
if (exists("obj")) {
  suppressPackageStartupMessages(library(SingleCellExperiment))
  x <- obj
  cat("dims", paste(dim(x), collapse = "x"), "\n")
  cd <- colnames(colData(x))
  cat("colData:", paste(cd, collapse = ", "), "\n")
  for (c in cd) {
    if (grepl("cell|azimuth|type|subclass|cluster", c, ignore.case = TRUE)) {
      cat("HIT", c, "\n")
      print(utils::head(sort(table(as.character(colData(x)[[c]])), decreasing = TRUE), 30))
    }
  }
  cat("assays:", paste(assayNames(x), collapse = ", "), "\n")
}
cat("DONE_INSPECT\n")
