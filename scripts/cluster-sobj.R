# Cluster one preprocessed normalization branch across a candidate grid.
# CJS 2026-07-02
# Inputs: preprocessed `.rds` at `--input <path>` (required); `--elbow-n <N>`
# (required; integer chosen from the elbow plot); optional `--extra-dims <csv>`
# (default "30,50"); optional `--resolutions <csv>` (default "0.3,0.5,0.8").
# Outputs: clustered object saved to `CURRENT_OBJECT_DIR/cluster_<norm>_elbow<N>.rds`;
# candidate-clustering UMAP overlays and clustree plots.
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
  is.finite(cli_args$elbow_n)
)

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

sobj <- readRDS(cli_args$input)
norm <- sobj@misc$preprocessing$normalization

dims_grid <- sort(unique(c(cli_args$elbow_n, cli_args$extra_dims)))
candidate_names <- character()

for (d in dims_grid) {
  old_idents <- SeuratObject::Idents(sobj)
  sobj <- Seurat::FindNeighbors(sobj, reduction = "pca", dims = 1:d)
  for (r in cli_args$resolutions) {
    name <- sprintf("cluster_%s_dims%d_res%s", norm, d, res_tag(r))
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
  reduction_name <- sprintf("umap_%s_dims%d", norm, d)
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
      umap = sprintf("umap_%s_dims%d", norm, d),
      color_by = sprintf("cluster_%s_dims%d_res%s", norm, d, res_tag(r))
    )
  }
}

for (d in dims_grid) {
  splot_clustree(
    sobj,
    prefix = sprintf("cluster_%s_dims%d_res", norm, d),
    out_tag = sprintf("%s_dims%d", norm, d)
  )
}

sobj@misc$clustering <- list(
  algorithm = "leiden",
  resolutions = cli_args$resolutions,
  dims_grid = dims_grid,
  elbow_n = cli_args$elbow_n,
  candidate_names = candidate_names
)

dir.create(CURRENT_OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(
  CURRENT_OBJECT_DIR,
  sprintf("cluster_%s_elbow%d.rds", norm, cli_args$elbow_n)
)
saveRDS(sobj, out_path)
message("Saved ", out_path)
