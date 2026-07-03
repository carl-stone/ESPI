# Cluster one preprocessed normalization branch across a candidate grid.
# CJS 2026-07-02
# Inputs: preprocessed `.rds` at `--input <path>` (required); `--elbow-n <N>`
# (required; integer chosen from the elbow plot); optional `--extra-dims <csv>`
# (default "30,50"); optional `--resolutions <csv>` (default "0.3,0.5,0.8").
# Outputs: clustered object saved to a Seurat-safe branch-tagged path:
# `CURRENT_OBJECT_DIR/cluster_<normalization>_<cc_tag>_elbow<N>.rds`;
# Terms: see CONTEXT.md (normalization branch, candidate clustering, chosen clustering,
# pseudobulk sample, focused test).

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

elbow_n <- arg("--elbow-n")
stopifnot(!is.null(elbow_n), !identical(elbow_n, TRUE))

cli_args <- list(
  input = arg("--input"),
  elbow_n = as.integer(elbow_n),
  extra_dims = parse_csv_int(arg("--extra-dims"), default = c(30, 50)),
  resolutions = parse_csv_num(arg("--resolutions"), default = c(0.3, 0.5, 0.8))
)
stopifnot(
  !is.null(cli_args$input),
  length(cli_args$elbow_n) == 1,
  is.finite(cli_args$elbow_n),
  all(is.finite(cli_args$extra_dims)),
  all(cli_args$extra_dims > 0),
  all(is.finite(cli_args$resolutions)),
  all(cli_args$resolutions > 0)
)

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

sobj <- readRDS(cli_args$input)
norm <- sobj@misc$preprocessing$normalization
if (
  !is.character(norm) ||
    length(norm) != 1 ||
    !norm %in% c("log1p", "pflog")
) {
  stop("Missing or invalid preprocessing normalization: ", norm, call. = FALSE)
}
cc_tag <- if (isTRUE(sobj@misc$preprocessing$filtered_cell_cycle)) {
  "filter_cc"
} else {
  "no_filter_cc"
}
branch_tag <- sprintf("%s_%s", norm, cc_tag)
if (!grepl("^[A-Za-z0-9_]+$", branch_tag)) {
  stop("Unsafe clustering branch tag: ", branch_tag, call. = FALSE)
}
emit_tripwire_checkpoint(
  "cluster_input_available",
  input = cli_args$input,
  normalization = norm,
  filtered_cell_cycle = cc_tag,
  n_cells = ncol(sobj),
  n_features = nrow(sobj)
)

dims_grid <- sort(unique(c(cli_args$elbow_n, cli_args$extra_dims)))
candidate_names <- character()

for (d in dims_grid) {
  old_idents <- SeuratObject::Idents(sobj)
  sobj <- Seurat::FindNeighbors(sobj, reduction = "pca", dims = 1:d)
  for (r in cli_args$resolutions) {
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
  for (r in cli_args$resolutions) {
    splot_umap_by(
      sobj,
      umap = sprintf("umap_%s_dims%d", branch_tag, d),
      color_by = sprintf("cluster_%s_dims%d_res%s", branch_tag, d, res_tag(r))
    )
  }
}

min_clustree_resolutions <- 2L
clustree_plotted <- length(unique(cli_args$resolutions)) >=
  min_clustree_resolutions
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
  resolutions = cli_args$resolutions,
  dims_grid = dims_grid,
  elbow_n = cli_args$elbow_n,
  candidate_names = candidate_names,
  clustree_plotted = clustree_plotted
)

dir.create(CURRENT_OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(
  CURRENT_OBJECT_DIR,
  sprintf("cluster_%s_elbow%d.rds", branch_tag, cli_args$elbow_n)
)
saveRDS(sobj, out_path)
emit_tripwire_checkpoint(
  "cluster_artifacts_written",
  output = out_path,
  branch_tag = branch_tag,
  n_candidates = length(candidate_names)
)
message("Saved ", out_path)
