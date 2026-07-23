#!/usr/bin/env Rscript

# Rank significant markers by absolute effect size for the fixed MG-selected
# clustering with Wilcox tests and make a dotplot.

suppressPackageStartupMessages({
  here::i_am("scripts/03-marker-analysis.R")
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(tidyverse)
})

# ---- parameters ----

config <- publication_config()
input_path <- config$selected$mg$path
cluster_column <- config$selected$mg$column
assay <- "RNA"
expression_layer <- "data"
counts_layer <- "counts"
top_n <- 5L
significance_threshold <- 0.01

# ---- paths ----

table_dir <- file.path(config$paths$tables, "mg_selected")
figure_dir <- file.path(config$paths$figures, "mg_selected")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

analysis_tag <- "data_pflog_mg_selected_no_filter_cc_dims20_res0.3"
cluster_tag <- "pflog_mg_selected_no_filter_cc_dims20_res0.3"
dotplot_tag <- paste0(
  "mg_selected_cluster_marker_dotplot_",
  analysis_tag,
  "_top5"
)
full_marker_path <- file.path(
  table_dir,
  paste0("find_all_markers_wilcox_", analysis_tag, ".csv")
)
top_marker_path <- file.path(
  table_dir,
  paste0("find_all_markers_wilcox_top5_", analysis_tag, ".csv")
)
summary_path <- file.path(
  table_dir,
  paste0("find_all_markers_summary_", analysis_tag, ".csv")
)
identity_map_path <- file.path(
  table_dir,
  paste0("find_all_markers_identity_map_", cluster_tag, ".csv")
)
figure_stem <- file.path(figure_dir, dotplot_tag)

assert_output_available(
  c(
    full_marker_path,
    top_marker_path,
    summary_path,
    identity_map_path,
    paste0(figure_stem, ".png"),
    paste0(figure_stem, ".pdf")
  ),
  config$overwrite
)

# ---- marker analysis ----

sobj <- readRDS(input_path)
cluster_values <- as.character(sobj[[cluster_column, drop = TRUE]])
identity_levels <- as.character(sort(as.integer(unique(cluster_values))))
marker_identities <- factor(cluster_values, levels = identity_levels)
SeuratObject::Idents(sobj) <- marker_identities

identity_map <- tibble::tibble(
  source_cluster = identity_levels,
  marker_identity = identity_levels,
  n_cells = as.integer(table(marker_identities)),
  decision_source = "confirmed_no_merge",
  input_path = input_path,
  cluster_column = cluster_column,
  assay = assay,
  expression_layer = expression_layer,
  counts_layer = counts_layer
)

markers <- Seurat::FindAllMarkers(
  object = sobj,
  assay = assay,
  slot = expression_layer,
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.01,
  logfc.threshold = 0.1,
  min.diff.pct = -Inf,
  return.thresh = significance_threshold,
  verbose = FALSE
) |>
  tibble::as_tibble() |>
  dplyr::mutate(
    gene = as.character(gene),
    cluster = factor(as.character(cluster), levels = identity_levels),
    pct_diff = pct.1 - pct.2
  ) |>
  dplyr::filter(p_val_adj <= significance_threshold) |>
  dplyr::arrange(
    cluster,
    dplyr::desc(abs(avg_log2FC)),
    p_val_adj,
    p_val,
    gene
  ) |>
  dplyr::group_by(cluster) |>
  dplyr::mutate(rank_within_cluster = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::relocate(
    gene,
    cluster,
    rank_within_cluster,
    p_val,
    avg_log2FC,
    pct.1,
    pct.2,
    pct_diff
  )

marker_counts <- markers |>
  dplyr::count(cluster) |>
  tidyr::complete(
    cluster = factor(identity_levels, levels = identity_levels),
    fill = list(n = 0L)
  )
if (any(marker_counts$n < top_n)) {
  stop(
    "Fewer than ",
    top_n,
    " significant Wilcox markers are available for at least one cluster.",
    call. = FALSE
  )
}

top_markers <- markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_head(n = top_n) |>
  dplyr::ungroup()

marker_summary <- identity_map |>
  dplyr::select(marker_identity, n_cells) |>
  dplyr::left_join(
    markers |>
      dplyr::count(cluster) |>
      dplyr::rename(n_significant_markers = n),
    by = c("marker_identity" = "cluster")
  ) |>
  dplyr::left_join(
    top_markers |> dplyr::count(cluster) |> dplyr::rename(n_top_markers = n),
    by = c("marker_identity" = "cluster")
  ) |>
  dplyr::mutate(
    dplyr::across(dplyr::starts_with("n_"), ~ tidyr::replace_na(.x, 0L)),
    decision_source = "confirmed_no_merge"
  )

# ---- marker dotplot ----

cell_order <- colnames(sobj)
identity_membership <- Matrix::sparseMatrix(
  i = seq_along(marker_identities),
  j = as.integer(marker_identities),
  x = 1,
  dims = c(length(marker_identities), length(identity_levels)),
  dimnames = list(cell_order, identity_levels)
)
identity_sizes <- Matrix::colSums(identity_membership)
blue_ramp <- grDevices::colorRampPalette(c(
  config$palettes$dotplot[[1]],
  "white"
))(4L)
pink_ramp <- grDevices::colorRampPalette(c(
  "white",
  config$palettes$dotplot[[2]]
))(4L)

make_marker_dotplot <- function(top_marker_data) {
  marker_rows <- paste0("marker_", seq_len(nrow(top_marker_data)))
  marker_labels <- stats::setNames(top_marker_data$gene, marker_rows)
  expression_matrix <- SeuratObject::LayerData(
    sobj[[assay]],
    layer = expression_layer
  )[top_marker_data$gene, cell_order, drop = FALSE]
  counts_matrix <- SeuratObject::LayerData(sobj[[assay]], layer = counts_layer)[
    top_marker_data$gene,
    cell_order,
    drop = FALSE
  ]
  rownames(expression_matrix) <- marker_rows
  rownames(counts_matrix) <- marker_rows

  mean_expression_matrix <- expression_matrix %*%
    identity_membership |>
    sweep(2, identity_sizes, "/")
  pct_detected_matrix <- 100 *
    sweep((counts_matrix > 0) %*% identity_membership, 2, identity_sizes, "/")
  mean_expression_values <- as.vector(t(as.matrix(mean_expression_matrix)))
  pct_detected_values <- as.vector(t(as.matrix(pct_detected_matrix)))

  plot_data <- tidyr::expand_grid(
    marker_row = marker_rows,
    marker_identity = identity_levels
  ) |>
    dplyr::mutate(
      mean_expression = .env$mean_expression_values,
      pct_detected = .env$pct_detected_values
    ) |>
    dplyr::group_by(marker_row) |>
    dplyr::mutate(
      scaled_mean_expression = if (stats::sd(mean_expression) == 0) {
        rep(0, dplyr::n())
      } else {
        as.numeric(scale(mean_expression))
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      marker_row = factor(marker_row, levels = rev(marker_rows)),
      marker_identity = factor(marker_identity, levels = identity_levels)
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = marker_identity,
      y = marker_row,
      size = pct_detected,
      color = scaled_mean_expression
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(labels = marker_labels) +
    ggplot2::scale_size(
      range = c(0.5, 6),
      limits = c(0, 100),
      name = "Detected cells (%)"
    ) +
    ggplot2::scale_color_stepsn(
      colours = c(
        config$palettes$dotplot[[1]],
        blue_ramp[[2]],
        blue_ramp[[3]],
        pink_ramp[[2]],
        pink_ramp[[3]],
        config$palettes$dotplot[[2]]
      ),
      breaks = c(-3, -2, -1, 0, 1, 2, 3),
      labels = c("<= -2", "-2", "-1", "0", "1", "2", ">= 2"),
      limits = c(-3, 3),
      oob = scales::squish,
      name = "Mean data expression\n(row z-score bin)"
    ) +
    ggplot2::labs(
      title = "MG-selected significant Wilcox markers",
      x = "Cluster number",
      y = "Marker gene"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

plot <- make_marker_dotplot(top_markers)

# ---- output ----

utils::write.csv(markers, full_marker_path, row.names = FALSE, na = "")
utils::write.csv(top_markers, top_marker_path, row.names = FALSE, na = "")
utils::write.csv(marker_summary, summary_path, row.names = FALSE, na = "")
utils::write.csv(identity_map, identity_map_path, row.names = FALSE, na = "")

save_publication_plot(
  plot,
  figure_stem,
  width = max(7, 2.5 + 0.45 * length(identity_levels)),
  height = max(5, 2 + 0.18 * nrow(top_markers)),
  notebook_basename = paste0(dotplot_tag, ".png")
)

message("Saved MG-selected Wilcox marker tables and dotplot under ", table_dir)
