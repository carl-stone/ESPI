#!/usr/bin/env Rscript

# Generate the publication marker heatmap with the marker list below.

suppressPackageStartupMessages({
  here::i_am("scripts/scratch-custom-curated-marker-heatmap.R")
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(ggview)
  library(tidyverse)
})

# ---- parameters ----

config <- publication_config()
branch_settings <- config$selected$mg
expression_layer <- "pflog"
output_stem <- here::here(
  "analysis",
  "exploratory",
  "figures",
  "marker_heatmap_cells"
)
heatmap_font_size <- 8
heatmap_z_limit <- 2
heatmap_group_label_width <- grid::unit(1.15, "in")
heatmap_colors <- c(
  config$palettes$dotplot[1],
  "white",
  config$palettes$dotplot[2]
)
cluster_palette <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")
marker_dot_plot_stem <- file.path(dirname(output_stem), "marker_dotplot")
dot_plot_width <- 7
dot_plot_height <- 9
dot_plot_dpi <- 300
dot_plot_font_size <- 11
dot_plot_max_size <- 6

# ---- validation ----

heatmap_marker_list <- list(
  "Muller glia" = c(
    "Rlbp1",
    "Glul",
    "Vim",
    "Slc1a3",
    "Sox9",
    "Hes1",
    "Aqp4",
    "Kcnj10"
  ),
  "Proliferative" = c("Pcna", "Mcm2", "Mcm6", "Ccnd1", "Cdk4", "Cdk6"),
  "Neurogenic progenitor" = c("Ascl1", "Hes6", "Hes5", "Dll1", "Neurog2"),
  "Cone bipolar" = c("Otx2", "Cabp5", "Scgn", "Lhx4", "Grik1", "Neurod1")
)

if (
  !is.list(heatmap_marker_list) ||
    is.null(names(heatmap_marker_list)) ||
    any(names(heatmap_marker_list) == "")
) {
  stop("heatmap_marker_list must be a fully named list.", call. = FALSE)
}

marker_table <- stack(heatmap_marker_list)
colnames(marker_table) <- c("gene", "cell_type")
duplicated_markers <- marker_table$gene[duplicated(marker_table$gene)]
if (length(duplicated_markers) > 0L) {
  stop(
    "Marker gene(s) assigned to more than one cell type: ",
    paste(unique(duplicated_markers), collapse = ", "),
    call. = FALSE
  )
}
marker_table$cell_type_label <- as.character(marker_table$cell_type)

# ---- inputs ----

sobj <- readRDS(branch_settings$path)
cluster_column <- branch_settings$column
assay <- SeuratObject::DefaultAssay(sobj)
expected_cluster_levels <- as.character(1:5)
cluster_values <- as.character(sobj[[cluster_column, drop = TRUE]])
observed_cluster_levels <- sort(unique(cluster_values))
if (!identical(observed_cluster_levels, expected_cluster_levels)) {
  stop(
    "Expected clusters 1 through 5; found: ",
    paste(observed_cluster_levels, collapse = ", "),
    call. = FALSE
  )
}

missing_markers <- setdiff(marker_table$gene, rownames(sobj))
if (length(missing_markers) > 0L) {
  stop(
    "Marker gene(s) missing from the Seurat object: ",
    paste(missing_markers, collapse = ", "),
    call. = FALSE
  )
}

# ---- heatmap data ----

marker_expression <- SeuratObject::GetAssayData(
  sobj,
  assay = assay,
  layer = expression_layer
)[marker_table$gene, , drop = FALSE]
scaled_expression <- t(scale(t(as.matrix(marker_expression))))
scaled_expression[is.na(scaled_expression)] <- 0
scaled_expression[scaled_expression > heatmap_z_limit] <- heatmap_z_limit
scaled_expression[scaled_expression < -heatmap_z_limit] <- -heatmap_z_limit

cell_cluster_labels <- factor(
  paste("Cluster", cluster_values),
  levels = paste("Cluster", expected_cluster_levels)
)
cell_type_groups <- factor(
  marker_table$cell_type_label,
  levels = names(heatmap_marker_list)
)
cluster_colors <- stats::setNames(
  cluster_palette,
  paste("Cluster", expected_cluster_levels)
)

row_annotation <- ComplexHeatmap::rowAnnotation(
  `Cell type` = ComplexHeatmap::anno_block(
    gp = grid::gpar(fill = NA, col = NA),
    labels = levels(cell_type_groups),
    labels_gp = grid::gpar(fontsize = heatmap_font_size),
    labels_rot = 0,
    labels_just = "right",
    labels_offset = grid::unit(1, "npc"),
    width = heatmap_group_label_width
  ),
  Gene = ComplexHeatmap::anno_text(
    marker_table$gene,
    gp = grid::gpar(fontsize = heatmap_font_size),
    just = "right",
    location = grid::unit(1, "npc")
  ),
  show_annotation_name = FALSE
)
column_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Cluster = ComplexHeatmap::anno_block(panel_fun = function(index, nm) {
    grid::grid.rect(gp = grid::gpar(fill = cluster_colors[[nm]], col = NA))
    grid::grid.text(
      sub("^Cluster ", "", nm),
      gp = grid::gpar(fontsize = heatmap_font_size)
    )
  }),
  show_annotation_name = FALSE
)

heatmap <- ComplexHeatmap::Heatmap(
  scaled_expression,
  name = "Relative\nexpression",
  col = circlize::colorRamp2(
    c(-heatmap_z_limit, 0, heatmap_z_limit),
    heatmap_colors
  ),
  left_annotation = row_annotation,
  top_annotation = column_annotation,
  row_split = cell_type_groups,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = cell_cluster_labels,
  cluster_column_slices = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_title = NULL,
  row_title = NULL,
  use_raster = TRUE
)

view_heatmap <- function(heatmap, width, height, dpi = 300) {
  heatmap_grob <- grid::grid.grabExpr(
    ComplexHeatmap::draw(heatmap),
    width = width,
    height = height
  )

  ggplotify::as.ggplot(heatmap_grob) +
    ggview::canvas(width, height, units = "in", dpi = dpi)
}

# ---- figure ----

dir.create(dirname(output_stem), recursive = TRUE, showWarnings = FALSE)
png_path <- paste0(output_stem, ".png")
pdf_path <- paste0(output_stem, ".pdf")


view_heatmap(heatmap, 6.5, 3.7)

heatmap_width <- 6.5
heatmap_height <- 3.7
heatmap_dpi <- 300

heatmap_plot <- grid::grid.grabExpr(
  ComplexHeatmap::draw(heatmap),
  width = heatmap_width,
  height = heatmap_height
)

ggplot2::ggsave(
  png_path,
  heatmap_plot,
  width = heatmap_width,
  height = heatmap_height,
  dpi = heatmap_dpi
)
ggplot2::ggsave(
  pdf_path,
  heatmap_plot,
  width = heatmap_width,
  height = heatmap_height
)

# ---- marker dot-plot data ----

cluster_average_expression <- vapply(
  expected_cluster_levels,
  function(cluster_value) {
    Matrix::rowMeans(marker_expression[,
      cluster_values == cluster_value,
      drop = FALSE
    ])
  },
  numeric(nrow(marker_expression))
)
percent_expressing <- vapply(
  expected_cluster_levels,
  function(cluster_value) {
    100 *
      Matrix::rowMeans(
        marker_expression[, cluster_values == cluster_value, drop = FALSE] > 0
      )
  },
  numeric(nrow(marker_expression))
)
rownames(cluster_average_expression) <- marker_table$gene
rownames(percent_expressing) <- marker_table$gene
colnames(cluster_average_expression) <- paste(
  "Cluster",
  expected_cluster_levels
)
colnames(percent_expressing) <- colnames(cluster_average_expression)

scaled_cluster_average_expression <- t(scale(t(cluster_average_expression)))
scaled_cluster_average_expression[is.na(scaled_cluster_average_expression)] <- 0
scaled_cluster_average_expression[scaled_cluster_average_expression > 2] <- 2
scaled_cluster_average_expression[scaled_cluster_average_expression < -2] <- -2

dot_plot_data <- data.frame(
  gene = rep(marker_table$gene, times = length(expected_cluster_levels)),
  cluster = rep(
    colnames(cluster_average_expression),
    each = nrow(marker_table)
  ),
  average_expression = as.vector(cluster_average_expression),
  average_expression_z = as.vector(scaled_cluster_average_expression),
  percent_expressing = as.vector(percent_expressing)
)
dot_plot_data$gene <- factor(
  dot_plot_data$gene,
  levels = rev(marker_table$gene)
)
dot_plot_data$cluster <- factor(
  dot_plot_data$cluster,
  levels = paste("Cluster", expected_cluster_levels)
)

marker_dot_plot <- ggplot2::ggplot(
  dot_plot_data,
  ggplot2::aes(x = cluster, y = gene)
) +
  ggplot2::geom_point(ggplot2::aes(
    color = average_expression_z,
    size = percent_expressing
  )) +
  ggplot2::scale_color_gradient2(
    low = heatmap_colors[1],
    mid = heatmap_colors[2],
    high = heatmap_colors[3],
    limits = c(-heatmap_z_limit, heatmap_z_limit),
    name = "Mean expression\nrow z-score"
  ) +
  ggplot2::scale_size_continuous(
    range = c(0, dot_plot_max_size),
    limits = c(0, 100),
    name = "% expressing"
  ) +
  ggplot2::scale_x_discrete(labels = expected_cluster_levels) +
  ggplot2::labs(x = "Cluster", y = NULL) +
  ggplot2::theme_bw(base_size = dot_plot_font_size) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

# ---- marker dot plot ----

marker_dot_plot_png_path <- paste0(marker_dot_plot_stem, ".png")
marker_dot_plot_pdf_path <- paste0(marker_dot_plot_stem, ".pdf")

ggplot2::ggsave(
  marker_dot_plot_png_path,
  marker_dot_plot,
  width = dot_plot_width,
  height = dot_plot_height,
  dpi = dot_plot_dpi
)
ggplot2::ggsave(
  marker_dot_plot_pdf_path,
  marker_dot_plot,
  width = dot_plot_width,
  height = dot_plot_height
)
