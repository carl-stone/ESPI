#!/usr/bin/env Rscript

# Generate per-cluster marker module and p27 enrichment heatmaps.
#
# Usage:
#   Rscript scripts/plot-cluster-marker-heatmaps.R \
#     [--input <clustered-seurat-object.rds>] \
#     [--dims <positive integer>] \
#     [--resolution <resolution string>] \
#     [--layer <assay layer>] \
#     [--slot <module-score slot/layer>] \
#     [--n-perm <positive integer>] \
#     [--out-dir <output directory>]
#
# Arguments:
#   --input
#     Clustered Seurat object to plot. Defaults to the full-dataset clustered
#     RDS in CURRENT_OBJECT_DIR.
#   --dims
#     PC count embedded in the cluster metadata column name. Defaults to 50.
#   --resolution
#     Leiden resolution embedded in the cluster metadata column name. Defaults
#     to 0.3.
#   --layer
#     Assay layer used for p27 expression values. Defaults to pflog because the
#     p27 strip summarizes the PFlog expression signal used for cluster review.
#   --slot
#     Assay slot/layer used for module scoring. Defaults to data to mirror
#     Seurat::AddModuleScore() defaults and the MG-selection module scoring.
#   --n-perm
#     Number of within-sample cluster-label permutations for descriptive
#     cluster-level p27 enrichment. Defaults to 2000.
#   --out-dir
#     Directory for PNG/PDF outputs. Defaults to FIGURE_DIR/annotation.
#
# Outputs:
#   Writes PNG/PDF files under --out-dir with the
#   cell_type_module_p27_heatmap_<layer>_<branch>_dims<dims>_res<resolution>
#   stem.
#   Writes module-score and p27-enrichment TSVs under TABLE_DIR/annotation.
#   Creates or replaces notebook/figures/<png filename> as a symlink to the PNG.
#
# The cluster metadata column is derived from the input object as
# cluster_<sobj@misc$clustering$branch_tag>_dims<dims>_res<resolution>.
# The input object must have been produced by scripts/cluster-sobj.R.

suppressPackageStartupMessages({
  library(here)
})
here::i_am("scripts/plot-cluster-marker-heatmaps.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})
palette_dotplot_pair <- get(
  "palette_dotplot_pair",
  envir = asNamespace("ESPI"),
  inherits = FALSE
)
compute_cluster_module_scores <- get(
  "compute_cluster_module_scores",
  envir = asNamespace("ESPI"),
  inherits = FALSE
)
compute_cluster_p27_enrichment <- get(
  "compute_cluster_p27_enrichment",
  envir = asNamespace("ESPI"),
  inherits = FALSE
)

# ---- parameters ----

get_arg <- function(args, flag, default) {
  match_index <- match(flag, args)
  if (is.na(match_index)) {
    return(default)
  }
  if (match_index == length(args)) {
    stop("Missing value for ", flag, call. = FALSE)
  }
  args[[match_index + 1]]
}

cli_args <- commandArgs(trailingOnly = TRUE)
allowed_flags <- c(
  "--input",
  "--dims",
  "--resolution",
  "--layer",
  "--slot",
  "--n-perm",
  "--out-dir"
)
unknown_flags <- cli_args[
  startsWith(cli_args, "--") & !cli_args %in% allowed_flags
]
if (length(unknown_flags) > 0) {
  stop(
    "Unknown argument(s): ",
    paste(unknown_flags, collapse = ", "),
    call. = FALSE
  )
}

required_packages <- c("ComplexHeatmap", "circlize", "Matrix", "SeuratObject")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

HEATMAP_Z_SCORE_LIMIT <- 2

dims <- as.integer(get_arg(cli_args, "--dims", "50"))
if (is.na(dims) || dims <= 0) {
  stop("--dims must be a positive integer.", call. = FALSE)
}
resolution <- get_arg(cli_args, "--resolution", "0.3")
if (!nzchar(resolution)) {
  stop("--resolution must not be empty.", call. = FALSE)
}
expression_layer <- get_arg(cli_args, "--layer", "pflog")
if (!nzchar(expression_layer)) {
  stop("--layer must not be empty.", call. = FALSE)
}
score_slot <- get_arg(cli_args, "--slot", "data")
if (!nzchar(score_slot)) {
  stop("--slot must not be empty.", call. = FALSE)
}
n_perm <- as.integer(get_arg(cli_args, "--n-perm", "2000"))
if (is.na(n_perm) || n_perm <= 0) {
  stop("--n-perm must be a positive integer.", call. = FALSE)
}

input_path <- get_arg(
  cli_args,
  "--input",
  file.path(CURRENT_OBJECT_DIR, "cluster_pflog_no_filter_cc_elbow20.rds")
)
out_dir <- get_arg(cli_args, "--out-dir", file.path(FIGURE_DIR, "annotation"))

# ---- validation ----

if (!file.exists(input_path)) {
  stop("Input Seurat object does not exist: ", input_path, call. = FALSE)
}

if (!identical(names(cell_type_marker_genes), names(cell_type_marker_labels))) {
  stop(
    "cell_type_marker_genes and cell_type_marker_labels names must match.",
    call. = FALSE
  )
}

# ---- work ----

sobj <- readRDS(input_path)
if (!inherits(sobj, "Seurat")) {
  stop("Input is not a Seurat object: ", input_path, call. = FALSE)
}
if (!"RNA" %in% SeuratObject::Assays(sobj)) {
  stop("Input Seurat object does not contain an RNA assay.", call. = FALSE)
}
SeuratObject::DefaultAssay(sobj) <- "RNA"
if (!inherits(sobj[["RNA"]], "Assay5")) {
  sobj[["RNA"]] <- as(sobj[["RNA"]], Class = "Assay5")
}
assay <- "RNA"

branch_tag <- sobj@misc$clustering$branch_tag
if (
  !is.character(branch_tag) ||
    length(branch_tag) != 1 ||
    is.na(branch_tag) ||
    !nzchar(branch_tag) ||
    !grepl("^[A-Za-z0-9_]+$", branch_tag)
) {
  stop(
    "Input object is missing a valid sobj@misc$clustering$branch_tag; ",
    "run scripts/cluster-sobj.R first.",
    call. = FALSE
  )
}
cluster_column <- sprintf(
  "cluster_%s_dims%s_res%s",
  branch_tag,
  dims,
  resolution
)
out_tag <- sprintf(
  "cell_type_module_p27_heatmap_%s_%s_dims%s_res%s",
  expression_layer,
  branch_tag,
  dims,
  resolution
)
if (!cluster_column %in% colnames(sobj@meta.data)) {
  stop("Missing cluster metadata column: ", cluster_column, call. = FALSE)
}

available_layers <- SeuratObject::Layers(sobj[[assay]])
if (!expression_layer %in% available_layers) {
  stop(
    "Missing expression layer '",
    expression_layer,
    "' in assay ",
    assay,
    ". Available layers: ",
    paste(available_layers, collapse = ", "),
    call. = FALSE
  )
}
if (!score_slot %in% available_layers) {
  stop(
    "Missing score slot/layer '",
    score_slot,
    "' in assay ",
    assay,
    ". Available layers: ",
    paste(available_layers, collapse = ", "),
    call. = FALSE
  )
}
if (!"Cdkn1b" %in% rownames(sobj)) {
  stop("Cdkn1b is missing from the Seurat object.", call. = FALSE)
}

marker_table <- stack(cell_type_marker_genes)
colnames(marker_table) <- c("gene", "cell_type")
duplicated_markers <- marker_table$gene[duplicated(marker_table$gene)]
if (length(duplicated_markers) > 0) {
  stop(
    "Marker gene(s) assigned to more than one cell type: ",
    paste(unique(duplicated_markers), collapse = ", "),
    call. = FALSE
  )
}
marker_table$cell_type_label <- unname(
  cell_type_marker_labels[marker_table$cell_type]
)

missing_markers <- setdiff(marker_table$gene, rownames(sobj))
if (length(missing_markers) > 0) {
  stop(
    "Marker gene(s) missing from the Seurat object: ",
    paste(missing_markers, collapse = ", "),
    call. = FALSE
  )
}

module_matrix <- compute_cluster_module_scores(
  sobj,
  cluster_column,
  cell_type_marker_genes,
  assay = assay,
  slot = score_slot
)
rownames(module_matrix) <- unname(cell_type_marker_labels[rownames(
  module_matrix
)])
module_z <- t(scale(t(module_matrix)))
module_z[is.na(module_z)] <- 0
module_z[module_z > HEATMAP_Z_SCORE_LIMIT] <- HEATMAP_Z_SCORE_LIMIT
module_z[module_z < -HEATMAP_Z_SCORE_LIMIT] <- -HEATMAP_Z_SCORE_LIMIT

p27 <- compute_cluster_p27_enrichment(
  sobj,
  cluster_column,
  layer = expression_layer,
  n_perm = n_perm
)
p27_z <- p27$z_score[match(colnames(module_z), p27$cluster)]
if (
  !identical(
    p27$cluster[match(colnames(module_z), p27$cluster)],
    colnames(module_z)
  )
) {
  stop(
    "p27 enrichment rows do not align to module-score columns.",
    call. = FALSE
  )
}

body_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c(palette_dotplot_pair[1], "white", palette_dotplot_pair[2])
)
zlim <- max(abs(p27_z), na.rm = TRUE)
if (!is.finite(zlim) || zlim == 0) {
  zlim <- 1
}
p27_fun <- circlize::colorRamp2(
  c(-zlim, 0, zlim),
  c(palette_dotplot_pair[1], "white", palette_dotplot_pair[2])
)
p27_legend <- ComplexHeatmap::Legend(
  col_fun = p27_fun,
  title = "p27 z-score"
)
heatmap <- ComplexHeatmap::Heatmap(
  module_z,
  name = "Module z-score",
  col = body_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = ComplexHeatmap::HeatmapAnnotation(
    `p27 z-score` = ComplexHeatmap::anno_simple(
      p27_z,
      col = p27_fun,
      na_col = "grey85"
    ),
    annotation_name_gp = grid::gpar(fontsize = 8),
    show_legend = FALSE
  ),
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 8),
  column_title = NULL,
  use_raster = FALSE
)

# ---- output ----

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(TABLE_DIR, "annotation")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
png_path <- file.path(out_dir, paste0(out_tag, ".png"))
pdf_path <- file.path(out_dir, paste0(out_tag, ".pdf"))
module_scores_path <- file.path(
  table_dir,
  paste0(out_tag, "_module_scores.tsv")
)
p27_path <- file.path(table_dir, paste0(out_tag, "_p27_enrichment.tsv"))

module_scores_out <- data.frame(
  cell_type = rownames(module_matrix),
  module_matrix,
  check.names = FALSE
)
utils::write.table(
  module_scores_out,
  module_scores_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  p27,
  p27_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

grDevices::png(png_path, width = 8, height = 6, units = "in", res = 300)
ComplexHeatmap::draw(
  heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  heatmap_legend_list = list(p27_legend)
)
grDevices::dev.off()

grDevices::pdf(pdf_path, width = 8, height = 6)
ComplexHeatmap::draw(
  heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  heatmap_legend_list = list(p27_legend)
)
grDevices::dev.off()

notebook_figure_dir <- here::here("notebook", "figures")
dir.create(notebook_figure_dir, recursive = TRUE, showWarnings = FALSE)
notebook_png_path <- file.path(notebook_figure_dir, basename(png_path))
if (file.exists(notebook_png_path) || nzchar(Sys.readlink(notebook_png_path))) {
  unlink(notebook_png_path)
}
link_created <- file.symlink(png_path, notebook_png_path)
if (!isTRUE(link_created)) {
  stop("Failed to link notebook figure: ", notebook_png_path, call. = FALSE)
}

message("Wrote module score TSV: ", module_scores_path)
message("Wrote p27 enrichment TSV: ", p27_path)
message("Wrote marker module p27 heatmap PNG: ", png_path)
message("Wrote marker module p27 heatmap PDF: ", pdf_path)
message("Linked notebook figure: ", notebook_png_path)
