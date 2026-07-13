#!/usr/bin/env Rscript

# Cluster one preprocessed normalization branch across a candidate grid.
#
# Usage:
#   Rscript scripts/04-cluster.R \
#     --input <preprocessed-seurat-object.rds> \
#     --elbow-n <positive integer> \
#     --extra-dims <comma-separated integers> \
#     --resolutions <comma-separated numbers>
#
# Arguments:
#   --input
#     Preprocessed Seurat object. Required.
#   --elbow-n
#     Primary PC count selected from the elbow plot. Required.
#   --extra-dims
#     Additional PC counts to cluster. Defaults to 30,50.
#   --resolutions
#     Leiden resolutions to cluster. Defaults to 0.3,0.5,0.8.
#
#   CURRENT_OBJECT_DIR/cluster_<normalization>_<cc_tag>_elbow<N>.rds
#   CURRENT_OBJECT_DIR/cluster_<normalization>_<dataset_tag>_<cc_tag>_elbow<N>.rds
#   UMAP overlays and clustree figures under FIGURE_DIR/cluster.
#
# Notes:
#   Cluster, UMAP, clustree, and RDS artifact tags use Seurat-safe underscores.

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
})
here::i_am("scripts/04-cluster.R")
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
parse_csv_int <- function(x, default) {
  if (is.null(x)) {
    return(default)
  }
  as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}
parse_csv_num <- function(x, default) {
  if (is.null(x)) {
    return(default)
  }
  as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}
res_tag <- function(x) {
  format(x, trim = TRUE, scientific = FALSE)
}

input <- arg_value("--input", required = TRUE)
elbow_n <- as.integer(arg_value("--elbow-n", required = TRUE))
extra_dims <- parse_csv_int(
  arg_value("--extra-dims", default = NULL),
  default = c(30, 50)
)
resolutions <- parse_csv_num(
  arg_value("--resolutions", default = NULL),
  default = c(0.3, 0.5, 0.8)
)

# ---- validation ----

if (!file.exists(input)) {
  stop("Input Seurat object does not exist: ", input, call. = FALSE)
}
if (
  length(elbow_n) != 1 ||
    is.na(elbow_n) ||
    !is.finite(elbow_n) ||
    elbow_n <= 0
) {
  stop("--elbow-n must be a positive integer.", call. = FALSE)
}
if (
  any(is.na(extra_dims)) || any(!is.finite(extra_dims)) || any(extra_dims <= 0)
) {
  stop("--extra-dims must contain positive integers.", call. = FALSE)
}
if (
  any(is.na(resolutions)) ||
    any(!is.finite(resolutions)) ||
    any(resolutions <= 0)
) {
  stop("--resolutions must contain positive numbers.", call. = FALSE)
}

# ---- work ----

sobj <- readRDS(input)
norm <- sobj@misc$preprocessing$normalization
if (
  !is.character(norm) ||
    length(norm) != 1 ||
    !norm %in% c("log1p", "pflog")
) {
  stop("Missing or invalid preprocessing normalization: ", norm, call. = FALSE)
}
dataset_tag <- sobj@misc$preprocessing$dataset_tag
if (is.null(dataset_tag) || identical(dataset_tag, "all_cells")) {
  dataset_tag <- NULL
} else if (
  !is.character(dataset_tag) ||
    length(dataset_tag) != 1 ||
    is.na(dataset_tag) ||
    !nzchar(dataset_tag) ||
    !grepl("^[A-Za-z0-9_]+$", dataset_tag)
) {
  stop(
    "Missing or invalid preprocessing dataset_tag: ",
    paste(dataset_tag, collapse = ", "),
    call. = FALSE
  )
}
cc_tag <- if (isTRUE(sobj@misc$preprocessing$filtered_cell_cycle)) {
  "filter_cc"
} else {
  "no_filter_cc"
}
branch_tag <- paste(c(norm, dataset_tag, cc_tag), collapse = "_")
if (!grepl("^[A-Za-z0-9_]+$", branch_tag)) {
  stop("Unsafe clustering branch tag: ", branch_tag, call. = FALSE)
}
emit_tripwire_checkpoint(
  "cluster_input_available",
  input = input,
  normalization = norm,
  filtered_cell_cycle = cc_tag,
  n_cells = ncol(sobj),
  n_features = nrow(sobj)
)

dims_grid <- sort(unique(c(elbow_n, extra_dims)))
candidate_names <- character()

for (d in dims_grid) {
  old_idents <- SeuratObject::Idents(sobj)
  sobj <- Seurat::FindNeighbors(sobj, reduction = "pca", dims = 1:d)
  for (r in resolutions) {
    name <- sprintf("cluster_%s_dims%d_res%s", branch_tag, d, res_tag(r))
    sobj <- Seurat::FindClusters(
      sobj,
      algorithm = 4,
      leiden_method = "igraph",
      resolution = r,
      random.seed = SEED
    )
    sobj@meta.data[[name]] <- SeuratObject::Idents(sobj)
    SeuratObject::Idents(sobj) <- old_idents
    candidate_names <- c(candidate_names, name)
  }
}

for (d in dims_grid) {
  reduction_name <- sprintf("umap_%s_dims%d", branch_tag, d)
  sobj <- Seurat::RunUMAP(
    sobj,
    reduction = "pca",
    dims = 1:d,
    reduction.name = reduction_name,
    reduction.key = paste0(gsub("[^A-Za-z0-9]", "", reduction_name), "_"),
    seed.use = SEED
  )
}

for (d in dims_grid) {
  for (r in resolutions) {
    splot_umap_by(
      sobj,
      umap = sprintf("umap_%s_dims%d", branch_tag, d),
      color_by = sprintf("cluster_%s_dims%d_res%s", branch_tag, d, res_tag(r))
    )
  }
}

min_clustree_resolutions <- 2L
clustree_plotted <- length(unique(resolutions)) >= min_clustree_resolutions
if (clustree_plotted) {
  for (d in dims_grid) {
    splot_clustree(
      sobj,
      prefix = sprintf("cluster_%s_dims%d_res", branch_tag, d),
      out_tag = sprintf("%s_dims%d", branch_tag, d)
    )
  }
} else {
  message(
    "Skipping clustree output because fewer than two resolutions were requested."
  )
}

sobj@misc$clustering <- list(
  algorithm = "leiden",
  filtered_cell_cycle = isTRUE(sobj@misc$preprocessing$filtered_cell_cycle),
  branch_tag = branch_tag,
  resolutions = resolutions,
  dims_grid = dims_grid,
  elbow_n = elbow_n,
  candidate_names = candidate_names,
  clustree_plotted = clustree_plotted
)

# ---- output ----

out_path <- file.path(
  CURRENT_OBJECT_DIR,
  sprintf("cluster_%s_elbow%d.rds", branch_tag, elbow_n)
)

dir.create(CURRENT_OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(sobj, out_path)
emit_tripwire_checkpoint(
  "cluster_artifacts_written",
  output = out_path,
  branch_tag = branch_tag,
  n_candidates = length(candidate_names)
)
message("Saved ", out_path)
