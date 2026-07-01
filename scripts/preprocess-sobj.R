# p(re)processing data from Trailmaker Seurat object
# CJS 2026-06-30

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name) {
  i <- match(name, args)
  if (is.na(i) || i == length(args)) NULL else args[[i + 1]]
}
normalization <- arg("--normalization")
if (is.null(normalization)) normalization <- "log1p"
stopifnot(normalization %in% c("log1p", "pflog"))
cli_args <- list(input = arg("--input"), normalization = normalization)

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
  library(scclrR)
  pkgload::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

in_path <- if (is.null(cli_args$input)) {
  file.path(INPUT_OBJECT_DIR, "pipseq_processed_matrix_with_egfp.rds")
} else {
  cli_args$input
}
sobj <- readRDS(in_path)
sobj[[c("pca_for_harmony", "harmony", "pca", "umap")]] <- NULL
sobj[["RNA"]] <- as(sobj[["RNA"]], Class = "Assay5")

sobj$sample_id <- paste0(
  "M",
  sobj$Mouse,
  "_",
  gsub("[^A-Za-z0-9]", "", sobj$Condition)
)
