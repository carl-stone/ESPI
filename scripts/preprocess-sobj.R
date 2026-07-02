# Preprocess a Trailmaker Seurat object into one normalization branch.
# CJS 2026-06-30
# Inputs: raw `.rds` at `INPUT_OBJECT_DIR` or `--input <path>`.
# Outputs: preprocessed object saved to `CURRENT_OBJECT_DIR/preprocess_<norm>.rds`
# where `<norm>` is `log1p` or `pflog`; QC/HVG/PCA-diagnostic plots saved by
# `splot_*` helpers.
# Branch note: this script produces one normalization branch per run; rerun with
# `--normalization pflog` for the other branch.
# Next step: `scripts/cluster-sobj.R` consumes this output with `--elbow-n <N>`
# after visual inspection of the elbow plot.
# Terms: see CONTEXT.md (normalization branch, candidate clustering, chosen clustering, pseudobulk sample, focused test).

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name) {
  i <- match(name, args)
  if (is.na(i)) {
    return(NULL)
  }
  if (i == length(args) || startsWith(args[[i + 1]], "--")) {
    return(TRUE)
  }
  args[[i + 1]]
}
normalization <- arg("--normalization")
if (is.null(normalization)) {
  normalization <- "log1p"
}
filter_cc_arg <- arg("--filter-cell-cycle")
filter_cc <- identical(filter_cc_arg, TRUE) ||
  (!is.null(filter_cc_arg) && tolower(filter_cc_arg) == "true")
cli_args <- list(
  input = arg("--input"),
  normalization = normalization,
  filter_cc = filter_cc
)
stopifnot(
  normalization %in% c("log1p", "pflog"),
  is.logical(cli_args$filter_cc)
)

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
  library(scclrR)
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

in_path <- if (is.null(cli_args$input)) {
  file.path(INPUT_OBJECT_DIR, "pipseq_processed_matrix_with_egfp.rds")
} else {
  cli_args$input
}
sobj <- readRDS(in_path)
# Drop any pre-existing reductions from upstream Trailmaker export.
sobj@reductions <- list()
sobj[["RNA"]] <- as(sobj[["RNA"]], Class = "Assay5")

sobj$sample_id <- paste0(
  "M",
  sobj$Mouse,
  "_",
  gsub("[^A-Za-z0-9]", "", sobj$Condition)
)

sobj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(
  sobj,
  pattern = "^Rp[sl]"
)

# Plot QC diagnostics.
splot_qc_metrics_violin(sobj)

sobj <- FindVariableFeatures(sobj, nfeatures = 2000)

# Plot gene mean-vs-variance scatter with top 10 HVGs labeled.
splot_hvg_scatter(sobj, n_top = 10)

if (cli_args$filter_cc) {
  VariableFeatures(sobj) <- setdiff(
    VariableFeatures(sobj),
    mouse_cell_cycle_genes
  )
}

sobj <- switch(
  cli_args$normalization,
  log1p = run_log1p_pca(sobj, n_pcs = 50),
  pflog = run_pflog_pca(sobj, n_pcs = 50),
  stop("Unknown normalization: ", cli_args$normalization, call. = FALSE)
)

sobj@misc$preprocessing$filter_cc <- cli_args$filter_cc

splot_viz_dim_loadings(sobj, n_pcs = 30)
splot_elbow(sobj, n_pcs = 50)

dir.create(CURRENT_OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(
  CURRENT_OBJECT_DIR,
  sprintf("preprocess_%s.rds", cli_args$normalization)
)
saveRDS(sobj, out_path)

message(
  "Saved ",
  out_path,
  ". Inspect the elbow plot, choose elbow_n, then run: ",
  "Rscript scripts/cluster-sobj.R --input ",
  out_path,
  " --elbow-n <N>"
)
