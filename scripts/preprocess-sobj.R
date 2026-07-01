# p(re)processing data from Trailmaker Seurat object
# CJS 2026-06-30
# This script takes a Seurat object and prepares it for downstream analysis.
# TODO: Describe the branches in this script, and the expected inputs and outputs.

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name) {
  i <- match(name, args)
  if (is.na(i) || i == length(args)) NULL else args[[i + 1]]
}
normalization <- arg("--normalization")
if (is.null(normalization)) {
  normalization <- "log1p"
}
stopifnot(normalization %in% c("log1p", "pflog"))
cli_args <- list(input = arg("--input"), normalization = normalization)
# TODO: Add argument parsing for cell-cycle gene filtering from VariableFeatures.

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

# TODO: splot_qc_metrics_violin(sobj) (See preprocess-plots.R)

sobj <- FindVariableFeatures(sobj, nfeatures = 2000)

# TODO: Filter cell-cycle genes from VariableFeatures
# IF cell-cycle filtering is requested (via argument)
#   all_cc_genes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
#   keep features <- setdiff(VariableFeatures(sobj), all_cc_genes)
#   VariableFeatures(sobj) <- keep features

# TODO: splot_hvg_scatter(sobj, n_top = 10) (See preprocess-plots.R)

# TODO: Add normalization and PCA steps here
# Args: (with defaults)
#   normalization method = cli_args$normalization
#   n_pcs = 50
# Output: seurat object with reduction saved to "pca",
#
# IF normalization == "log1p" (or NULL)
#   sobj <- Seurat::NormalizeData() with defaults.
#   sobj <- Seurat::RunPCA(npcs = n_pcs) with defaults.
#   set sobj@misc$preprocessing to record the normalization method and n_pcs used.
# ELIF normalization == "pflog"
#   sobj <- scclrR::pflog() with defaults.
#   scclrR::pca_matrix on variable features
#   set sobj@misc$preprocessing to record the normalization method and n_pcs used.

# TODO:
