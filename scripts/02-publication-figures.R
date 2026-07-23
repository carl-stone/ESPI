#!/usr/bin/env Rscript

# Generate publication figures and tables from the frozen analysis objects.

suppressPackageStartupMessages({
  here::i_am("scripts/02-publication-figures.R")
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(tidyverse)
})

# ---- parameters ----

config <- publication_config()
seed <- config$seed
condition_col <- config$conditions$column
control_label <- config$conditions$control
estim_label <- config$conditions$estim
control_display_label <- config$conditions$control_display
estim_display_label <- config$conditions$estim_display
expression_layer <- "pflog"
module_score_layer <- "data"
n_perm <- 2000L
source_suffix <- paste0(
  "_dims",
  config$selected$source$dimensions,
  "_res",
  config$selected$source$resolution
)
mg_settings <- list(config$selected$mg, config$selected$mg_filter_cc)
branch_tags <- vapply(mg_settings, `[[`, character(1), "branch")

# ---- paths ----

annotation_dir <- file.path(config$paths$figures, "annotation")
mg_figure_dir <- file.path(config$paths$figures, "mg_selected")
annotation_table_dir <- file.path(config$paths$tables, "annotation")
mg_table_dir <- file.path(config$paths$tables, "mg_selected")
notebook_figure_dir <- config$paths$notebook_figures

purrr::walk(
  c(annotation_dir, mg_figure_dir, annotation_table_dir, mg_table_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)

source_marker_stem <- file.path(
  annotation_dir,
  paste0(
    "cell_type_marker_heatmap_pflog_",
    config$selected$source$branch,
    "_cells",
    source_suffix
  )
)
source_module_stem <- file.path(
  annotation_dir,
  paste0(
    "cell_type_module_p27_heatmap_pflog_",
    config$selected$source$branch,
    source_suffix
  )
)

output_paths <- c(
  paste0(source_marker_stem, c(".png", ".pdf")),
  paste0(source_module_stem, c(".png", ".pdf")),
  file.path(
    annotation_table_dir,
    paste0(
      basename(source_module_stem),
      c("_module_scores.tsv", "_p27_enrichment.tsv")
    )
  )
)
for (branch_settings in mg_settings) {
  branch_tag <- branch_settings$branch
  branch_suffix <- paste0(
    "_dims",
    branch_settings$dimensions,
    "_res",
    branch_settings$resolution
  )
  mg_stems <- file.path(
    mg_figure_dir,
    c(
      paste0("mg_selected_cluster_umap_", branch_tag, branch_suffix),
      paste0("mg_selected_condition_umap_", branch_tag, branch_suffix),
      paste0("mg_selected_feature_umap_pflog_", branch_tag, branch_suffix),
      paste0("mg_selected_ascl1_hes6_coexpression_", branch_tag, branch_suffix),
      paste0(
        "mg_selected_cluster_abundance_enrichment_",
        branch_tag,
        branch_suffix
      ),
      paste0(
        "mg_selected_cluster_proportion_by_mouse_",
        branch_tag,
        branch_suffix
      )
    )
  )
  annotation_stems <- file.path(
    annotation_dir,
    c(
      paste0(
        "cell_type_marker_heatmap_pflog_",
        branch_tag,
        "_cells",
        branch_suffix
      ),
      paste0("cell_type_module_p27_heatmap_pflog_", branch_tag, branch_suffix)
    )
  )
  output_paths <- c(
    output_paths,
    paste0(mg_stems, rep(c(".png", ".pdf"), each = length(mg_stems))),
    paste0(
      annotation_stems,
      rep(c(".png", ".pdf"), each = length(annotation_stems))
    ),
    file.path(
      mg_table_dir,
      c(
        paste0(
          "mg_selected_cluster_abundance_enrichment_",
          branch_tag,
          paste0(branch_suffix, ".tsv")
        ),
        paste0(
          "mg_selected_cluster_proportion_randomization_",
          branch_tag,
          paste0(branch_suffix, ".tsv")
        ),
        paste0(
          "mg_selected_sample_cluster_proportions_",
          branch_tag,
          paste0(branch_suffix, ".tsv")
        )
      )
    ),
    file.path(
      annotation_table_dir,
      c(
        paste0(basename(annotation_stems[[2]]), "_module_scores.tsv"),
        paste0(basename(annotation_stems[[2]]), "_p27_enrichment.tsv")
      )
    )
  )
}
assert_output_available(output_paths, config$overwrite)

# ---- inputs ----

source_sobj <- readRDS(config$selected$source$path)
mg_sobj <- readRDS(config$selected$mg$path)
mg_filter_cc_sobj <- readRDS(config$selected$mg_filter_cc$path)
load(here::here("data", "umap_feature_list.rda"))
umap_features <- umap_feature_list

# ---- source annotation ----

source_marker_paths <- write_curated_marker_heatmap(
  source_sobj,
  config$selected$source$column,
  expression_layer,
  source_marker_stem,
  width = 10,
  height = 9
)
ESPI:::.copy_notebook_figure(
  source_marker_paths[["png"]],
  file.path(notebook_figure_dir, basename(source_marker_paths[["png"]]))
)

source_module_scores <- compute_cluster_module_scores(
  source_sobj,
  config$selected$source$column,
  cell_type_marker_genes,
  assay = "RNA",
  slot = module_score_layer,
  seed = seed
)
source_module_scores_out <- source_module_scores
rownames(source_module_scores_out) <- unname(cell_type_marker_labels[rownames(
  source_module_scores_out
)])
source_module_scores_out <- data.frame(
  cell_type = rownames(source_module_scores_out),
  source_module_scores_out,
  check.names = FALSE
)
source_p27 <- compute_cluster_p27_enrichment(
  source_sobj,
  config$selected$source$column,
  layer = expression_layer,
  assay = "RNA",
  mouse_col = "Mouse",
  condition_col = condition_col,
  n_perm = n_perm,
  seed = seed
)
source_module_paths <- write_module_p27_heatmap(
  source_module_scores,
  source_p27,
  source_module_stem,
  width = 8,
  height = 6
)
ESPI:::.copy_notebook_figure(
  source_module_paths[["png"]],
  file.path(notebook_figure_dir, basename(source_module_paths[["png"]]))
)
readr::write_tsv(
  source_module_scores_out,
  file.path(
    annotation_table_dir,
    paste0(basename(source_module_stem), "_module_scores.tsv")
  )
)
readr::write_tsv(
  source_p27,
  file.path(
    annotation_table_dir,
    paste0(basename(source_module_stem), "_p27_enrichment.tsv")
  )
)

# ---- MG figures and tables ----

branches <- list(
  list(
    sobj = mg_sobj,
    tag = branch_tags[[1]],
    settings = config$selected$mg,
    cluster_column = config$selected$mg$column,
    legend_position = c(0.98, 0.98),
    legend_justification = c(1, 1)
  ),
  list(
    sobj = mg_filter_cc_sobj,
    tag = branch_tags[[2]],
    settings = config$selected$mg_filter_cc,
    cluster_column = config$selected$mg_filter_cc$column,
    legend_position = c(0.98, 0.02),
    legend_justification = c(1, 0)
  )
)
condition_colors <- stats::setNames(
  config$palettes$dotplot,
  c(control_label, estim_label)
)
condition_labels <- stats::setNames(
  c(control_display_label, estim_display_label),
  c(control_label, estim_label)
)

for (branch in branches) {
  sobj <- branch$sobj
  branch_tag <- branch$tag
  cluster_column <- branch$cluster_column
  branch_suffix <- paste0(
    "_dims",
    branch$settings$dimensions,
    "_res",
    branch$settings$resolution
  )
  reduction <- paste0("umap_", branch_tag, "_dims", branch$settings$dimensions)
  settings_label <- sprintf(
    "MG-selected PFlog; %d PCs; res %s",
    branch$settings$dimensions,
    format(branch$settings$resolution, trim = TRUE, scientific = FALSE)
  )

  cluster_stem <- file.path(
    mg_figure_dir,
    paste0("mg_selected_cluster_umap_", branch_tag, branch_suffix)
  )
  condition_stem <- file.path(
    mg_figure_dir,
    paste0("mg_selected_condition_umap_", branch_tag, branch_suffix)
  )
  feature_stem <- file.path(
    mg_figure_dir,
    paste0("mg_selected_feature_umap_pflog_", branch_tag, branch_suffix)
  )
  coexpression_stem <- file.path(
    mg_figure_dir,
    paste0("mg_selected_ascl1_hes6_coexpression_", branch_tag, branch_suffix)
  )
  abundance_stem <- file.path(
    mg_figure_dir,
    paste0(
      "mg_selected_cluster_abundance_enrichment_",
      branch_tag,
      branch_suffix
    )
  )
  proportion_stem <- file.path(
    mg_figure_dir,
    paste0(
      "mg_selected_cluster_proportion_by_mouse_",
      branch_tag,
      branch_suffix
    )
  )
  randomization_stem <- paste0(
    "mg_selected_cluster_proportion_randomization_",
    branch_tag,
    branch_suffix
  )
  sample_props_stem <- paste0(
    "mg_selected_sample_cluster_proportions_",
    branch_tag,
    branch_suffix
  )
  marker_stem <- file.path(
    annotation_dir,
    paste0(
      "cell_type_marker_heatmap_pflog_",
      branch_tag,
      "_cells",
      branch_suffix
    )
  )
  module_stem <- file.path(
    annotation_dir,
    paste0("cell_type_module_p27_heatmap_pflog_", branch_tag, branch_suffix)
  )

  cluster_values <- as.character(sobj[[cluster_column, drop = TRUE]])
  cluster_levels <- as.character(sort(as.integer(unique(cluster_values))))
  SeuratObject::Idents(sobj) <- factor(cluster_values, levels = cluster_levels)

  cluster_plot <- Seurat::DimPlot(
    sobj,
    reduction = reduction,
    group.by = cluster_column,
    label = TRUE,
    pt.size = 0.5,
    stroke.size = 0
  ) +
    ggplot2::ggtitle(settings_label) +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2")

  condition_plot <- Seurat::DimPlot(
    sobj,
    reduction = reduction,
    group.by = condition_col,
    label = FALSE,
    pt.size = 0.5,
    stroke.size = 0
  ) +
    ggplot2::scale_color_manual(
      values = condition_colors,
      breaks = c(control_label, estim_label),
      labels = condition_labels,
      drop = FALSE
    ) +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2", color = "Condition") +
    ggplot2::theme(
      legend.position = branch$legend_position,
      legend.justification = branch$legend_justification,
      legend.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  cell_order <- rownames(sobj[[]])
  embeddings <- SeuratObject::Embeddings(sobj, reduction = reduction)[
    cell_order,
    1:2,
    drop = FALSE
  ]
  colnames(embeddings) <- c("UMAP_1", "UMAP_2")
  expression <- SeuratObject::GetAssayData(
    sobj,
    assay = "RNA",
    layer = expression_layer
  )[umap_features, cell_order, drop = FALSE]

  umap_span <- max(
    diff(range(embeddings[, "UMAP_1"])),
    diff(range(embeddings[, "UMAP_2"]))
  )
  umap_x_limits <- mean(range(embeddings[, "UMAP_1"])) +
    c(-0.5, 0.5) * umap_span
  umap_y_limits <- mean(range(embeddings[, "UMAP_2"])) +
    c(-0.5, 0.5) * umap_span

  feature_plots <- lapply(umap_features, function(feature) {
    feature_expression <- as.numeric(expression[feature, ])
    expression_range <- range(feature_expression)
    scaled_expression <- if (diff(expression_range) > 0) {
      (feature_expression - expression_range[[1]]) / diff(expression_range)
    } else {
      rep(0, length(feature_expression))
    }
    plot_data <- tibble::tibble(
      UMAP_1 = embeddings[, "UMAP_1"],
      UMAP_2 = embeddings[, "UMAP_2"],
      scaled_expression
    ) |>
      dplyr::arrange(scaled_expression)

    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = UMAP_1, y = UMAP_2, color = scaled_expression)
    ) +
      ggplot2::geom_point(size = 0.5, stroke = 0) +
      ggplot2::scale_color_gradient(
        low = "grey85",
        high = config$palettes$dotplot[[2]],
        limits = c(0, 1),
        breaks = c(0, 1),
        labels = c("0", "1"),
        name = "Scaled expression"
      ) +
      ggplot2::coord_equal(
        xlim = umap_x_limits,
        ylim = umap_y_limits,
        expand = FALSE
      ) +
      ggplot2::ggtitle(feature) +
      ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 10, hjust = 0.5),
        legend.position = "right"
      )
  })
  feature_plot <- patchwork::wrap_plots(feature_plots, ncol = 3L) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = "pflog layer marker expression") &
    ggplot2::theme(legend.position = "right")

  coexpression_expression <- SeuratObject::GetAssayData(
    sobj,
    assay = "RNA",
    layer = expression_layer
  )[c("Ascl1", "Hes6"), cell_order, drop = FALSE]
  coexpression_data <- tibble::tibble(
    Ascl1 = as.numeric(coexpression_expression["Ascl1", ]),
    Hes6 = as.numeric(coexpression_expression["Hes6", ]),
    cluster = factor(cluster_values, levels = cluster_levels)
  )
  coexpression_plot <- ggplot2::ggplot(
    coexpression_data,
    ggplot2::aes(x = Ascl1, y = Hes6, color = cluster)
  ) +
    ggplot2::geom_point(size = 0.6, alpha = 0.6, stroke = 0) +
    ggplot2::labs(
      x = "Ascl1 (pflog)",
      y = "Hes6 (pflog)",
      color = "Cluster",
      title = "Ascl1 / Hes6 coexpression by cell"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))
    ) +
    theme_stone(base_size = 11)

  abundance_table <- compute_cluster_abundance(
    sobj,
    cluster_column,
    condition_col,
    control_label,
    estim_label
  )
  abundance_plot_data <- abundance_table |>
    dplyr::arrange(dplyr::desc(log2_enrichment)) |>
    dplyr::mutate(
      cluster = factor(cluster, levels = cluster),
      direction = factor(
        direction,
        levels = c(
          "Enriched in E-Stim",
          "Depleted in E-Stim",
          "Not significant"
        )
      ),
      significance = dplyr::case_when(
        padj < 0.001 ~ "***",
        padj < 0.01 ~ "**",
        padj < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  abundance_label_offset <- max(
    0.12,
    0.05 * diff(range(abundance_plot_data$log2_enrichment))
  )
  abundance_plot_data <- abundance_plot_data |>
    dplyr::mutate(
      label_y = log2_enrichment +
        dplyr::if_else(
          log2_enrichment >= 0,
          abundance_label_offset,
          -abundance_label_offset
        ),
      label_vjust = dplyr::if_else(log2_enrichment >= 0, 0, 1)
    )
  abundance_plot <- ggplot2::ggplot(
    abundance_plot_data,
    ggplot2::aes(x = cluster, y = log2_enrichment, fill = direction)
  ) +
    ggplot2::geom_col(color = "black", width = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(y = label_y, label = significance, vjust = label_vjust),
      size = 4,
      show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::scale_fill_manual(
      values = c(
        "Enriched in E-Stim" = config$palettes$analysis[["high"]],
        "Depleted in E-Stim" = config$palettes$analysis[["low"]],
        "Not significant" = config$palettes$analysis[["mid"]]
      ),
      labels = c(
        "Enriched in E-Stim" = paste("Enriched in", estim_display_label),
        "Depleted in E-Stim" = paste("Depleted in", estim_display_label),
        "Not significant" = "Not significant"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.12, 0.18))
    ) +
    ggplot2::labs(
      x = "MG-selected cluster",
      y = paste("CLR log2 enrichment", config$conditions$contrast_display),
      fill = NULL,
      subtitle = "Pooled cell-level Fisher/CLR summary; descriptive relative to Mouse x Condition DE unit."
    ) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(sub("PFlog", "cluster abundance", settings_label))

  sample_props <- compute_sample_cluster_proportions(
    sobj,
    cluster_column,
    "Mouse",
    condition_col,
    control_label,
    estim_label
  )
  randomization_table <- test_cluster_proportion_randomization(
    sample_props,
    control_label,
    estim_label
  )
  proportion_data <- sample_props |>
    dplyr::mutate(
      condition = factor(
        as.character(condition),
        levels = c(control_label, estim_label)
      ),
      cluster = factor(as.character(cluster), levels = cluster_levels),
      mouse_role = factor(
        as.character(mouse_role),
        levels = c("paired", "estim_only", "control_only")
      )
    )
  paired_data <- proportion_data |> dplyr::filter(mouse_role == "paired")
  proportion_plot <- ggplot2::ggplot(
    proportion_data,
    ggplot2::aes(x = condition, y = proportion)
  ) +
    ggplot2::geom_line(
      data = paired_data,
      ggplot2::aes(group = mouse),
      color = "grey60",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = mouse, shape = mouse_role),
      size = 2
    ) +
    ggplot2::facet_wrap(~cluster, scales = "free_y") +
    ggplot2::scale_shape_manual(
      values = c(paired = 16, estim_only = 17, control_only = 15),
      drop = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Cluster proportion of MG-selected cells",
      color = "Mouse",
      shape = "Mouse role",
      subtitle = "Mouse x Condition sample-level proportions; paired mice connected, singletons shown separately."
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)) +
    ggplot2::ggtitle(sub(
      "PFlog",
      "cluster proportion by mouse",
      settings_label
    ))

  module_scores <- compute_cluster_module_scores(
    sobj,
    cluster_column,
    cell_type_marker_genes,
    assay = "RNA",
    slot = module_score_layer,
    seed = seed
  )
  module_scores_out <- module_scores
  rownames(module_scores_out) <- unname(cell_type_marker_labels[rownames(
    module_scores_out
  )])
  module_scores_out <- data.frame(
    cell_type = rownames(module_scores_out),
    module_scores_out,
    check.names = FALSE
  )
  p27_enrichment <- compute_cluster_p27_enrichment(
    sobj,
    cluster_column,
    layer = expression_layer,
    assay = "RNA",
    mouse_col = "Mouse",
    condition_col = condition_col,
    n_perm = n_perm,
    seed = seed
  )

  marker_paths <- write_curated_marker_heatmap(
    sobj,
    cluster_column,
    expression_layer,
    marker_stem,
    width = 10,
    height = 9,
    cluster_cells = TRUE
  )
  module_paths <- write_module_p27_heatmap(
    module_scores,
    p27_enrichment,
    module_stem,
    width = 8,
    height = 6
  )
  ESPI:::.copy_notebook_figure(
    marker_paths[["png"]],
    file.path(notebook_figure_dir, basename(marker_paths[["png"]]))
  )
  ESPI:::.copy_notebook_figure(
    module_paths[["png"]],
    file.path(notebook_figure_dir, basename(module_paths[["png"]]))
  )

  save_publication_plot(
    cluster_plot,
    cluster_stem,
    width = 5.5,
    height = 5,
    notebook_basename = paste0(basename(cluster_stem), ".png")
  )
  save_publication_plot(
    condition_plot,
    condition_stem,
    width = 5.5,
    height = 5,
    notebook_basename = paste0(basename(condition_stem), ".png")
  )
  save_publication_plot(
    feature_plot,
    feature_stem,
    width = 10.5,
    height = 9,
    notebook_basename = paste0(basename(feature_stem), ".png")
  )
  save_publication_plot(
    coexpression_plot,
    coexpression_stem,
    width = 6,
    height = 5,
    notebook_basename = paste0(basename(coexpression_stem), ".png")
  )
  save_publication_plot(
    abundance_plot,
    abundance_stem,
    width = max(6.5, 0.45 * nrow(abundance_table)),
    height = 4.5,
    notebook_basename = paste0(basename(abundance_stem), ".png")
  )
  save_publication_plot(
    proportion_plot,
    proportion_stem,
    width = 8,
    height = 6,
    notebook_basename = paste0(basename(proportion_stem), ".png")
  )

  readr::write_tsv(
    abundance_table,
    file.path(mg_table_dir, paste0(basename(abundance_stem), ".tsv"))
  )
  readr::write_tsv(
    randomization_table,
    file.path(mg_table_dir, paste0(randomization_stem, ".tsv"))
  )
  readr::write_tsv(
    sample_props,
    file.path(mg_table_dir, paste0(sample_props_stem, ".tsv"))
  )
  readr::write_tsv(
    module_scores_out,
    file.path(
      annotation_table_dir,
      paste0(basename(module_stem), "_module_scores.tsv")
    )
  )
  readr::write_tsv(
    p27_enrichment,
    file.path(
      annotation_table_dir,
      paste0(basename(module_stem), "_p27_enrichment.tsv")
    )
  )
}

message("Saved publication figures under ", config$paths$figures)
message("Saved publication tables under ", config$paths$tables)
