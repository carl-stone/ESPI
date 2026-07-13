#!/usr/bin/env Rscript

# Preprocess one Seurat object into one normalization branch.
#
# Usage:
#   Rscript scripts/preprocess-sobj.R \
#     --input <seurat-object.rds> | --input-source <legacy|counts-qc> \
#     --normalization <log1p|pflog> \
#     --filter-cell-cycle <true|false>
#
# Arguments:
#   --input
#     Explicit Seurat object to preprocess. Cannot be combined with
#     --input-source.
#   --input-source
#     Named source object: legacy selects
#     INPUT_OBJECT_DIR/pipseq_processed_matrix_with_egfp.rds; counts-qc selects
#     INPUT_OBJECT_DIR/sobj_qc_filtered.rds. Defaults to legacy.
#   --normalization
#     Normalization branch to run. Must be log1p or pflog. Defaults to log1p.
#   --filter-cell-cycle
#     If true, removes mouse_cell_cycle_genes from the HVG set before PCA.
#     Defaults to false.
#
# Outputs:
#   CURRENT_OBJECT_DIR/preprocess_<normalization>_<filter-cc|no-filter-cc>.rds
#   QC, HVG, DimHeatmap, and elbow figures under FIGURE_DIR/preprocess.
#
# Next step:
#   Run scripts/cluster-sobj.R with --input <preprocess output> and --elbow-n.

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
  library(scclrR)
})
here::i_am("scripts/preprocess-sobj.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

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
arg_value <- function(name, default = NULL, required = FALSE) {
  value <- arg(name)
  if (identical(value, TRUE)) {
    stop("Missing value for ", name, call. = FALSE)
  }
  if (is.null(value)) {
    if (required) {
      stop("Missing required argument ", name, call. = FALSE)
    }
    return(default)
  }
  value
}
arg_flag <- function(name) {
  identical(arg(name), TRUE)
}

input <- arg_value("--input", default = NULL)
input_source <- arg_value("--input-source", default = "legacy")
if (!is.null(input) && !is.null(arg("--input-source"))) {
  stop("Use either --input or --input-source, not both.", call. = FALSE)
}
if (!input_source %in% c("legacy", "counts-qc")) {
  stop(
    "--input-source must be one of legacy or counts-qc; got ",
    input_source,
    call. = FALSE
  )
}
if (is.null(input)) {
  input <- analysis_input_path(input_source)
} else {
  input_source <- "custom"
}
normalization <- arg_value("--normalization", default = "log1p")
filter_cc <- arg_flag("--filter-cell-cycle") ||
  identical(
    tolower(arg_value("--filter-cell-cycle", default = "false")),
    "true"
  )

# ---- validation ----

if (!normalization %in% c("log1p", "pflog")) {
  stop("Unknown normalization: ", normalization, call. = FALSE)
}
if (!file.exists(input)) {
  stop("Input Seurat object does not exist: ", input, call. = FALSE)
}

# ---- work ----

sobj <- readRDS(input)
emit_tripwire_checkpoint(
  "raw_data_available",
  input = input,
  input_source = input_source,
  n_cells = ncol(sobj),
  n_features = nrow(sobj)
)
validate_required_metadata(sobj@meta.data, c("Mouse", "Condition"))

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

sobj@misc$preprocessing$input_source <- input_source
sobj@misc$preprocessing$normalization <- normalization
sobj@misc$preprocessing$filtered_cell_cycle <- filter_cc

# Plot QC diagnostics.
splot_qc_metrics_violin(sobj)

sobj <- FindVariableFeatures(sobj, nfeatures = 2000)

if (filter_cc) {
  utils::data("mouse_cell_cycle_genes", package = "ESPI", envir = environment())
  VariableFeatures(sobj) <- setdiff(
    VariableFeatures(sobj),
    mouse_cell_cycle_genes
  )
}
emit_tripwire_checkpoint(
  "variable_features_selected",
  normalization = normalization,
  filtered_cell_cycle = filter_cc,
  n_variable_features = length(VariableFeatures(sobj))
)

# Plot gene mean-vs-variance scatter with top 20 retained HVGs labeled.
splot_hvg_scatter(sobj, n_top = 20)

sobj <- switch(
  normalization,
  log1p = run_log1p_pca(sobj, n_pcs = 50),
  pflog = run_pflog_pca(sobj, n_pcs = 50)
)
emit_tripwire_checkpoint(
  "pca_ready",
  normalization = normalization,
  filtered_cell_cycle = filter_cc,
  n_pcs = ncol(SeuratObject::Embeddings(sobj, reduction = "pca"))
)

splot_dim_heatmap(sobj)
splot_elbow(sobj, n_pcs = 50)

# ---- output ----

cc_tag <- if (filter_cc) "filter-cc" else "no-filter-cc"
out_path <- file.path(
  CURRENT_OBJECT_DIR,
  sprintf("preprocess_%s_%s.rds", normalization, cc_tag)
)

dir.create(CURRENT_OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(sobj, out_path)
emit_tripwire_checkpoint(
  "preprocess_object_written",
  output = out_path,
  normalization = normalization,
  filtered_cell_cycle = cc_tag
)

message(
  "Saved ",
  out_path,
  ". Inspect the elbow plot, choose elbow_n, then run: ",
  "Rscript scripts/cluster-sobj.R --input ",
  out_path,
  " --elbow-n <N>"
)
