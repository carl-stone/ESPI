# Cluster one preprocessed normalization branch across a candidate grid.
# CJS 2026-07-02
# Inputs: preprocessed `.rds` at `--input <path>` (required); `--elbow-n <N>`
# (required; integer chosen from the elbow plot); optional `--extra-dims <csv>`
# (default "30,50"); optional `--resolutions <csv>` (default "0.3,0.5,0.8").
# Outputs: clustered object saved to `CURRENT_OBJECT_DIR/cluster_<norm>_elbow<N>.rds`;
# candidate-clustering UMAP overlays and clustree plots.
# Terms: see CONTEXT.md (normalization branch, candidate clustering, chosen clustering,
# pseudobulk sample, focused test).

# CLI parsing pseudocode.
# args <- commandArgs(trailingOnly = TRUE)
# arg <- helper as in preprocess-sobj.R
# parse_csv_int <- function(x, default) {
#   If x is NULL, return default.
#   Otherwise split x on ",", trim whitespace, coerce to integer, and return.
# }
# parse_csv_num <- function(x, default) {
#   If x is NULL, return default.
#   Otherwise split x on ",", trim whitespace, coerce to numeric, and return.
# }
# cli_args <- list(
#   input = arg("--input"),                     # required path
#   elbow_n = as.integer(arg("--elbow-n")),     # required integer
#   extra_dims = parse_csv_int(arg("--extra-dims"), default = c(30, 50)),
#   resolutions = parse_csv_num(arg("--resolutions"), default = c(0.3, 0.5, 0.8))
# )
# stopifnot(!is.null(cli_args$input), is.finite(cli_args$elbow_n))

# 1. Load preprocessed branch object.
# sobj <- readRDS(cli_args$input).
# Assert "pca" %in% Reductions(sobj).
# Assert sobj@misc$preprocessing$normalization exists.
# norm <- sobj@misc$preprocessing$normalization.

# 2. Build dims grid.
# dims_grid <- sort(unique(c(cli_args$elbow_n, cli_args$extra_dims))).
# Record dims_grid on sobj@misc$clustering$dims_grid.

# 3. Compute neighbors + Leiden clusters across grid.
# For each d in dims_grid, for each r in cli_args$resolutions:
# name <- sprintf("cluster_%s_dims%d_res%s", norm, d, gsub("\\.", "", format(r))).
# sobj <- find_leiden_clusters(sobj, dims = 1:d, resolution = r, name = name).
# The helper writes labels into sobj@meta.data[[name]] and returns sobj.

# 4. UMAP per (branch × dims).
# For each d in dims_grid:
# sobj <- run_umap_for_dims(
#   sobj,
#   dims = 1:d,
#   reduction_name = sprintf("umap_%s_dims%d", norm, d)
# ).
# The reduction name lands in sobj@reductions.

# 5. Save candidate UMAP overlays.
# For each d in dims_grid, for each r in cli_args$resolutions:
# splot_umap_by(
#   sobj,
#   umap = sprintf("umap_%s_dims%d", norm, d),
#   color_by = sprintf("cluster_%s_dims%d_res%s", norm, d, gsub("\\.", "", format(r)))
# ).
# The helper embeds norm, d, and r in the output filename.

# 6. Save clustree per dims.
# For each d in dims_grid:
# splot_clustree(
#   sobj,
#   prefix = sprintf("cluster_%s_dims%d_res", norm, d),
#   out_tag = sprintf("%s_dims%d", norm, d)
# ).
# clustree varies resolution at fixed dims.

# 7. Record clustering provenance.
# Extend sobj@misc$clustering with fields:
# algorithm = "leiden", resolutions, dims_grid, elbow_n, and candidate_names.
# candidate_names is the vector of all cluster metadata columns written.

# 8. Save clustered object.
# saveRDS(
#   sobj,
#   file.path(CURRENT_OBJECT_DIR, sprintf("cluster_%s_elbow%d.rds", norm, cli_args$elbow_n))
# ).
