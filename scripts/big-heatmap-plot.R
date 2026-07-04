#!/usr/bin/env Rscript

# Generate a per-cell marker heatmap for the manuscript.
#
# Usage:
#   Rscript scripts/big-heatmap-plot.R \
#     --input <clustered-seurat-object.rds> \
#     --dims <positive integer> \
#     --resolution <resolution string> \
#     --layer <assay layer> \
#     --out-dir <output directory>
#
# Arguments:
#   --input
#     Clustered Seurat object to plot. Defaults to
#     CURRENT_OBJECT_DIR/cluster_pflog_filter_cc_elbow20.rds.
#   --dims
#     PC count embedded in the cluster metadata column name. Defaults to 50.
#   --resolution
#     Leiden resolution embedded in the cluster metadata column name. Defaults
#     to 0.3.
#   --layer
#     Assay layer used for marker expression values. Defaults to pflog.
#   --out-dir
#     Directory for PNG/PDF outputs. Defaults to FIGURE_DIR/annotation.
#
# Outputs:
#   Writes PNG and PDF files under --out-dir, named
#   cell_type_marker_heatmap_<layer>_cells_dims<dims>_res<resolution>.(png|pdf).
#   Creates or replaces notebook/figures/<png filename> as a symlink to the PNG.
#
# The cluster metadata column is constructed as
# cluster_pflog_filter_cc_dims<dims>_res<resolution>.

suppressPackageStartupMessages({
  library(here)
})
here::i_am("scripts/big-heatmap-plot.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

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
allowed_flags <- c("--input", "--dims", "--resolution", "--layer", "--out-dir")
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
expression_layer <- get_arg(cli_args, "--layer", "pflog")
cluster_column <- sprintf(
  "cluster_pflog_filter_cc_dims%s_res%s",
  dims,
  resolution
)

input_path <- get_arg(
  cli_args,
  "--input",
  file.path(CURRENT_OBJECT_DIR, "cluster_pflog_filter_cc_elbow20.rds")
)
out_dir <- get_arg(cli_args, "--out-dir", file.path(FIGURE_DIR, "annotation"))
out_tag <- sprintf(
  "cell_type_marker_heatmap_%s_cells_dims%s_res%s",
  expression_layer,
  dims,
  resolution
)

# ---- validation ----

if (!file.exists(input_path)) {
  stop("Input Seurat object does not exist: ", input_path, call. = FALSE)
}

if (!identical(names(cell_type_marker_genes), names(cell_type_marker_labels))) {
  stop(
    "cell_type_marker_genes and cell_type_marker_labels must have identical names.",
    call. = FALSE
  )
}

# ---- work ----

sobj <- readRDS(input_path)
if (!cluster_column %in% colnames(sobj@meta.data)) {
  stop("Missing cluster metadata column: ", cluster_column, call. = FALSE)
}
assay <- SeuratObject::DefaultAssay(sobj)
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

cluster_values <- as.character(sobj@meta.data[[cluster_column]])
unique_cluster_values <- unique(cluster_values)
if (all(grepl("^-?[0-9]+$", unique_cluster_values))) {
  cluster_levels <- as.character(sort(as.integer(unique_cluster_values)))
} else {
  cluster_levels <- sort(unique_cluster_values, method = "radix")
}
cell_cluster_labels <- factor(
  paste("Cluster", cluster_values),
  levels = paste("Cluster", cluster_levels)
)

marker_expression <- SeuratObject::GetAssayData(
  sobj,
  assay = assay,
  layer = expression_layer
)[marker_table$gene, , drop = FALSE]
if (!identical(colnames(marker_expression), rownames(sobj@meta.data))) {
  stop(
    "Expression matrix columns do not match Seurat metadata rows.",
    call. = FALSE
  )
}

scaled_expression <- t(scale(t(as.matrix(marker_expression))))
scaled_expression[is.na(scaled_expression)] <- 0
scaled_expression[scaled_expression > HEATMAP_Z_SCORE_LIMIT] <-
  HEATMAP_Z_SCORE_LIMIT
scaled_expression[scaled_expression < -HEATMAP_Z_SCORE_LIMIT] <-
  -HEATMAP_Z_SCORE_LIMIT

cell_type_groups <- factor(
  marker_table$cell_type_label,
  levels = unname(cell_type_marker_labels)
)

cell_type_colors <- stats::setNames(
  grDevices::hcl.colors(nlevels(cell_type_groups), palette = "Dark 3"),
  levels(cell_type_groups)
)
cluster_colors <- stats::setNames(
  grDevices::hcl.colors(length(cluster_levels), palette = "Temps"),
  levels(cell_cluster_labels)
)

row_annotation <- ComplexHeatmap::rowAnnotation(
  `Cell type` = cell_type_groups,
  col = list(`Cell type` = cell_type_colors),
  show_annotation_name = FALSE
)
column_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Cluster = cell_cluster_labels,
  col = list(Cluster = cluster_colors),
  show_annotation_name = TRUE
)

heatmap <- ComplexHeatmap::Heatmap(
  scaled_expression,
  name = "Row z-score",
  col = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  left_annotation = row_annotation,
  top_annotation = column_annotation,
  row_split = cell_type_groups,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  column_split = cell_cluster_labels,
  show_column_names = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  column_title = sprintf(
    "%s (%s layer; %s cells)",
    cluster_column,
    expression_layer,
    ncol(scaled_expression)
  ),
  row_title = NULL,
  use_raster = TRUE
)

# ---- output ----

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
png_path <- file.path(out_dir, paste0(out_tag, ".png"))
pdf_path <- file.path(out_dir, paste0(out_tag, ".pdf"))

grDevices::png(png_path, width = 9, height = 11, units = "in", res = 300)
ComplexHeatmap::draw(
  heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
grDevices::dev.off()

grDevices::pdf(pdf_path, width = 9, height = 11)
ComplexHeatmap::draw(
  heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
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

message("Wrote marker heatmap PNG: ", png_path)
message("Wrote marker heatmap PDF: ", pdf_path)
message("Linked notebook figure: ", notebook_png_path)
