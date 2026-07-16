#!/usr/bin/env Rscript

# Generate all publication figures and tables from the three frozen analysis objects.

suppressPackageStartupMessages({
  library(here)
})
here::i_am("scripts/02-publication-figures.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
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
source_branch <- "pflog_no_filter_cc"
mg_branch <- "pflog_mg_selected_no_filter_cc"
mg_filter_cc_branch <- "pflog_mg_selected_filter_cc"
source_dims <- 30L
source_resolution <- "0.3"
mg_dims <- 20L
mg_resolution <- "0.5"

source_spec <- config$selected$source
mg_spec <- config$selected$mg
mg_filter_cc_spec <- config$selected$mg_filter_cc

annotation_dir <- file.path(config$paths$figures, "annotation")
mg_figure_dir <- file.path(config$paths$figures, "mg_selected")
annotation_table_dir <- file.path(config$paths$tables, "annotation")
mg_table_dir <- file.path(config$paths$tables, "mg_selected")
notebook_figure_dir <- config$paths$notebook_figures

# ---- fixed output inventory ----

source_marker_stem <- file.path(
  annotation_dir,
  sprintf(
    "cell_type_marker_heatmap_%s_%s_cells_dims%d_res%s",
    expression_layer,
    source_branch,
    source_dims,
    source_resolution
  )
)
source_module_stem <- file.path(
  annotation_dir,
  sprintf(
    "cell_type_module_p27_heatmap_%s_%s_dims%d_res%s",
    expression_layer,
    source_branch,
    source_dims,
    source_resolution
  )
)

branch_specs <- list(
  list(
    sobj_name = "mg",
    branch = mg_branch,
    cluster_column = mg_spec$column
  ),
  list(
    sobj_name = "mg_filter_cc",
    branch = mg_filter_cc_branch,
    cluster_column = mg_filter_cc_spec$column
  )
)

branch_output_stems <- lapply(branch_specs, function(spec) {
  tags <- list(
    cluster = sprintf(
      "mg_selected_cluster_umap_%s_dims%d_res%s",
      spec$branch,
      mg_dims,
      mg_resolution
    ),
    condition = sprintf(
      "mg_selected_condition_umap_%s_dims%d_res%s",
      spec$branch,
      mg_dims,
      mg_resolution
    ),
    feature = sprintf(
      "mg_selected_feature_umap_%s_%s_dims%d_res%s",
      expression_layer,
      spec$branch,
      mg_dims,
      mg_resolution
    ),
    coexpression = sprintf(
      "mg_selected_ascl1_hes6_coexpression_%s_dims%d_res%s",
      spec$branch,
      mg_dims,
      mg_resolution
    ),
    abundance = sprintf(
      "mg_selected_cluster_abundance_enrichment_%s_dims%d_res%s",
      spec$branch,
      mg_dims,
      mg_resolution
    ),
    proportion = sprintf(
      "mg_selected_cluster_proportion_by_mouse_%s_dims%d_res%s",
      spec$branch,
      mg_dims,
      mg_resolution
    ),
    randomization = sprintf(
      "mg_selected_cluster_proportion_randomization_%s_dims%d_res%s",
      spec$branch,
      mg_dims,
      mg_resolution
    ),
    sample_props = sprintf(
      "mg_selected_sample_cluster_proportions_%s_dims%d_res%s",
      spec$branch,
      mg_dims,
      mg_resolution
    ),
    marker = sprintf(
      "cell_type_marker_heatmap_%s_%s_cells_dims%d_res%s",
      expression_layer,
      spec$branch,
      mg_dims,
      mg_resolution
    ),
    module = sprintf(
      "cell_type_module_p27_heatmap_%s_%s_dims%d_res%s",
      expression_layer,
      spec$branch,
      mg_dims,
      mg_resolution
    )
  )
  stems <- stats::setNames(
    lapply(tags, function(tag) file.path(mg_figure_dir, tag)),
    names(tags)
  )
  stems[c("marker", "module")] <- lapply(
    tags[c("marker", "module")],
    function(tag) file.path(annotation_dir, tag)
  )
  stems
})

output_paths <- c(
  paste0(source_marker_stem, c(".png", ".pdf")),
  paste0(source_module_stem, c(".png", ".pdf")),
  file.path(
    annotation_table_dir,
    paste0(
      basename(source_module_stem),
      c("_module_scores.tsv", "_p27_enrichment.tsv")
    )
  ),
  unlist(
    lapply(branch_output_stems, function(stems) {
      c(
        paste0(stems$marker, c(".png", ".pdf")),
        paste0(stems$module, c(".png", ".pdf")),
        paste0(stems$cluster, c(".png", ".pdf")),
        paste0(stems$condition, c(".png", ".pdf")),
        paste0(stems$feature, c(".png", ".pdf")),
        paste0(stems$coexpression, c(".png", ".pdf")),
        paste0(stems$abundance, c(".png", ".pdf")),
        paste0(stems$proportion, c(".png", ".pdf")),
        file.path(mg_table_dir, paste0(basename(stems$abundance), ".tsv")),
        file.path(mg_table_dir, paste0(basename(stems$randomization), ".tsv")),
        file.path(mg_table_dir, paste0(basename(stems$sample_props), ".tsv")),
        file.path(
          annotation_table_dir,
          paste0(
            basename(stems$module),
            c("_module_scores.tsv", "_p27_enrichment.tsv")
          )
        )
      )
    }),
    use.names = FALSE
  )
)

notebook_png_paths <- c(
  paste0(source_marker_stem, ".png"),
  paste0(source_module_stem, ".png"),
  unlist(
    lapply(branch_output_stems, function(stems) {
      paste0(
        stems[c(
          "marker",
          "module",
          "cluster",
          "condition",
          "feature",
          "coexpression",
          "abundance",
          "proportion"
        )],
        ".png"
      )
    }),
    use.names = FALSE
  )
)

# Fail before loading any scientific input when a fixed primary output is protected.
assert_output_available(output_paths, config$overwrite)

# Runtime notebook copies must already be regular files. Publication writers never
# write through symlinks and the mirror below replaces only existing regular files.
for (png_path in notebook_png_paths) {
  destination <- file.path(notebook_figure_dir, basename(png_path))
  destination_link <- Sys.readlink(destination)
  if (
    length(destination_link) == 1L &&
      !is.na(destination_link) &&
      nzchar(destination_link)
  ) {
    stop(
      "Refusing to replace symlinked notebook figure: ",
      destination,
      call. = FALSE
    )
  }
  destination_info <- file.info(destination)
  if (
    !file.exists(destination) ||
      is.na(destination_info$isdir) ||
      isTRUE(destination_info$isdir)
  ) {
    stop(
      "Notebook figure destination must be an existing regular file: ",
      destination,
      call. = FALSE
    )
  }
}

# ---- frozen inputs ----

source_sobj <- readRDS(source_spec$path)
assert_frozen_input(source_spec$path, source_sobj, config$frozen$source)
mg_sobj <- readRDS(mg_spec$path)
assert_frozen_input(mg_spec$path, mg_sobj, config$frozen$mg)
mg_filter_cc_sobj <- readRDS(mg_filter_cc_spec$path)
assert_frozen_input(
  mg_filter_cc_spec$path,
  mg_filter_cc_sobj,
  config$frozen$mg_filter_cc
)

if (!all(c("RNA") %in% SeuratObject::Assays(source_sobj))) {
  stop("Source Seurat object must contain an RNA assay.", call. = FALSE)
}
if (!all(c("RNA") %in% SeuratObject::Assays(mg_sobj))) {
  stop("MG-selected Seurat object must contain an RNA assay.", call. = FALSE)
}
if (!all(c("RNA") %in% SeuratObject::Assays(mg_filter_cc_sobj))) {
  stop(
    "Filtered-CC MG-selected Seurat object must contain an RNA assay.",
    call. = FALSE
  )
}

feature_env <- new.env(parent = emptyenv())
load(here::here("data", "umap_feature_list.rda"), envir = feature_env)
if (!exists("umap_feature_list", envir = feature_env, inherits = FALSE)) {
  stop(
    "data/umap_feature_list.rda must define umap_feature_list.",
    call. = FALSE
  )
}
umap_features <- get("umap_feature_list", envir = feature_env, inherits = FALSE)
# ANALYSIS_OK[R002]: fixed nine-feature panel cardinality preserves the audited notebook feature list.
if (
  !is.character(umap_features) ||
    length(umap_features) != 9L ||
    anyNA(umap_features) ||
    any(!nzchar(umap_features)) ||
    anyDuplicated(umap_features) > 0L
) {
  stop(
    "umap_feature_list must contain nine unique non-empty gene names.",
    call. = FALSE
  )
}

# ---- source marker and module heatmaps ----

source_marker_paths <- write_curated_marker_heatmap(
  source_sobj,
  source_spec$column,
  expression_layer,
  source_marker_stem,
  width = 10,
  height = 9
)
source_module_scores <- compute_cluster_module_scores(
  source_sobj,
  source_spec$column,
  cell_type_marker_genes,
  assay = "RNA",
  slot = module_score_layer,
  seed = seed
)
source_module_scores_out <- source_module_scores
rownames(source_module_scores_out) <- unname(
  cell_type_marker_labels[rownames(source_module_scores_out)]
)
source_module_scores_out <- data.frame(
  cell_type = rownames(source_module_scores_out),
  source_module_scores_out,
  check.names = FALSE
)
source_p27 <- compute_cluster_p27_enrichment(
  source_sobj,
  source_spec$column,
  layer = expression_layer,
  assay = "RNA",
  mouse_col = "Mouse",
  condition_col = condition_col,
  n_perm = n_perm,
  seed = seed
)
source_module_paths <- write_module_p27_heatmap(
  source_module_scores[,
    setdiff(colnames(source_module_scores), "label"),
    drop = FALSE
  ],
  source_p27,
  source_module_stem,
  width = 8,
  height = 6
)
dir.create(annotation_table_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  source_module_scores_out,
  file.path(
    annotation_table_dir,
    paste0(basename(source_module_stem), "_module_scores.tsv")
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  source_p27,
  file.path(
    annotation_table_dir,
    paste0(basename(source_module_stem), "_p27_enrichment.tsv")
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ---- MG UMAP, feature, coexpression, abundance, and module figures ----

branch_objects <- list(mg = mg_sobj, mg_filter_cc = mg_filter_cc_sobj)
heatmap_png_paths <- c(
  source_marker_paths[["png"]],
  source_module_paths[["png"]]
)

for (branch_index in seq_along(branch_specs)) {
  spec <- branch_specs[[branch_index]]
  sobj <- branch_objects[[spec$sobj_name]]
  stems <- branch_output_stems[[branch_index]]
  cluster_column <- spec$cluster_column
  reduction <- sprintf("umap_%s_dims%d", spec$branch, mg_dims)

  if (!reduction %in% names(sobj@reductions)) {
    stop("Missing UMAP reduction: ", reduction, call. = FALSE)
  }
  if (!cluster_column %in% colnames(sobj@meta.data)) {
    stop("Missing cluster metadata column: ", cluster_column, call. = FALSE)
  }
  if (!condition_col %in% colnames(sobj@meta.data)) {
    stop("Missing condition metadata column: ", condition_col, call. = FALSE)
  }
  if (
    anyNA(sobj@meta.data[[cluster_column]]) ||
      anyNA(sobj@meta.data[[condition_col]])
  ) {
    stop(
      "Selected cluster and condition metadata must not contain NA values.",
      call. = FALSE
    )
  }
  condition_values <- as.character(sobj@meta.data[[condition_col]])
  if (
    !all(unique(condition_values) %in% c(control_label, estim_label)) ||
      !all(c(control_label, estim_label) %in% unique(condition_values))
  ) {
    stop(
      "Selected object must contain both expected condition labels.",
      call. = FALSE
    )
  }
  assay <- "RNA"
  available_layers <- SeuratObject::Layers(sobj[[assay]])
  if (
    !expression_layer %in% available_layers ||
      !module_score_layer %in% available_layers
  ) {
    stop("Selected object lacks required RNA expression layers.", call. = FALSE)
  }
  missing_features <- setdiff(
    c(umap_features, "Ascl1", "Hes6", "Cdkn1b"),
    rownames(sobj)
  )
  if (length(missing_features) > 0L) {
    stop(
      "Selected object is missing required feature(s): ",
      paste(missing_features, collapse = ", "),
      call. = FALSE
    )
  }

  cluster_values <- as.character(sobj@meta.data[[cluster_column]])
  cluster_levels <- unique(cluster_values)
  if (all(grepl("^-?[0-9]+$", cluster_levels))) {
    cluster_levels <- as.character(sort(as.integer(cluster_levels)))
  } else {
    cluster_levels <- sort(cluster_levels, method = "radix")
  }
  SeuratObject::Idents(sobj) <- factor(cluster_values, levels = cluster_levels)

  cluster_plot <- Seurat::DimPlot(
    sobj,
    reduction = reduction,
    group.by = cluster_column,
    label = TRUE,
    pt.size = 0.5,
    stroke.size = 0
  ) +
    ggplot2::ggtitle(sprintf(
      "MG-selected PFlog; %d PCs; res %s",
      mg_dims,
      mg_resolution
    )) +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2")
  condition_colors <- stats::setNames(
    config$palettes$dotplot,
    c(control_label, estim_label)
  )
  condition_labels <- stats::setNames(
    c(control_display_label, estim_display_label),
    c(control_label, estim_label)
  )
  condition_filter_cc <- identical(spec$branch, mg_filter_cc_branch)
  condition_legend_position <- if (condition_filter_cc) {
    c(0.98, 0.02)
  } else {
    c(0.98, 0.98)
  }
  condition_legend_justification <- if (condition_filter_cc) {
    c(1, 0)
  } else {
    c(1, 1)
  }
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
      legend.position = condition_legend_position,
      legend.justification = condition_legend_justification,
      legend.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  embeddings <- SeuratObject::Embeddings(sobj, reduction = reduction)
  # ANALYSIS_OK[R002]: fixed two-coordinate UMAP requirement preserves the plotted reduction contract.
  if (ncol(embeddings) < 2L) {
    stop("UMAP reduction must expose at least two coordinates.", call. = FALSE)
  }
  if (!setequal(rownames(embeddings), rownames(sobj@meta.data))) {
    stop(
      "UMAP reduction cell names do not match Seurat metadata.",
      call. = FALSE
    )
  }
  # ANALYSIS_OK[R001]: select the first two named UMAP coordinates as the audited plotting axes.
  # ANALYSIS_OK[R005]: reorder embeddings to metadata cell order and select UMAP coordinates; no cells are dropped.
  embeddings <- embeddings[
    rownames(sobj@meta.data),
    seq_len(2L),
    drop = FALSE
  ]
  colnames(embeddings) <- c("UMAP_1", "UMAP_2")
  expression <- SeuratObject::GetAssayData(
    sobj,
    assay = assay,
    layer = expression_layer
  )[umap_features, rownames(sobj@meta.data), drop = FALSE]
  if (!identical(colnames(expression), rownames(sobj@meta.data))) {
    stop(
      "Feature-expression cell order does not match Seurat metadata.",
      call. = FALSE
    )
  }
  umap_span <- max(
    diff(range(embeddings[, "UMAP_1"])),
    diff(range(embeddings[, "UMAP_2"]))
  )
  umap_x_center <- mean(range(embeddings[, "UMAP_1"]))
  umap_y_center <- mean(range(embeddings[, "UMAP_2"]))
  umap_x_limits <- umap_x_center + c(-0.5, 0.5) * umap_span
  umap_y_limits <- umap_y_center + c(-0.5, 0.5) * umap_span
  feature_plots <- lapply(umap_features, function(feature) {
    feature_expression <- as.numeric(expression[feature, ])
    expression_range <- range(feature_expression, na.rm = TRUE)
    if (!all(is.finite(expression_range))) {
      stop(
        "Expression values are not finite for feature: ",
        feature,
        call. = FALSE
      )
    }
    expression_min <- min(expression_range)
    expression_max <- max(expression_range)
    scaled_expression <- if (expression_max > expression_min) {
      (feature_expression - expression_min) / (expression_max - expression_min)
    } else {
      rep(0, length(feature_expression))
    }
    # ANALYSIS_OK[R019]: construct marker plot data from the validated feature and cell populations for visualization.
    plot_data <- data.frame(
      UMAP_1 = embeddings[, "UMAP_1"],
      UMAP_2 = embeddings[, "UMAP_2"],
      scaled_expression = scaled_expression,
      stringsAsFactors = FALSE
    )
    # ANALYSIS_OK[R005]: sort marker plot rows for deterministic rendering without dropping observations.
    # ANALYSIS_OK[R019]: row ordering is a visualization-only sort of the validated marker population.
    plot_data <- plot_data[order(plot_data$scaled_expression), , drop = FALSE]
    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data[["UMAP_1"]],
        y = .data[["UMAP_2"]],
        color = .data[["scaled_expression"]]
      )
    ) +
      ggplot2::geom_point(size = 0.5, stroke = 0) +
      ggplot2::scale_color_gradient(
        low = "grey85",
        # ANALYSIS_OK[R001]: the second named palette entry is the audited high-expression color.
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
    patchwork::plot_annotation(
      title = sprintf("%s layer marker expression", expression_layer)
    ) &
    ggplot2::theme(legend.position = "right")

  coexpression_expression <- SeuratObject::GetAssayData(
    sobj,
    assay = assay,
    layer = expression_layer
  )[c("Ascl1", "Hes6"), rownames(sobj@meta.data), drop = FALSE]
  if (
    !identical(
      colnames(coexpression_expression),
      rownames(sobj@meta.data)
    )
  ) {
    stop(
      "Coexpression cell order does not match Seurat metadata.",
      call. = FALSE
    )
  }
  coexpression_data <- data.frame(
    x = as.numeric(coexpression_expression["Ascl1", ]),
    y = as.numeric(coexpression_expression["Hes6", ]),
    cluster = factor(cluster_values, levels = cluster_levels),
    stringsAsFactors = FALSE
  )
  coexpression_plot <- ggplot2::ggplot(
    coexpression_data,
    ggplot2::aes(x = .data[["x"]], y = .data[["y"]], color = .data[["cluster"]])
  ) +
    ggplot2::geom_point(size = 0.6, alpha = 0.6, stroke = 0) +
    ggplot2::labs(
      x = sprintf("Ascl1 (%s)", expression_layer),
      y = sprintf("Hes6 (%s)", expression_layer),
      color = "Cluster",
      title = "Ascl1 / Hes6 coexpression by cell"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))
    ) +
    theme_stone(base_size = 11)

  abundance_table <- compute_cluster_abundance(
    sobj = sobj,
    cluster_col = cluster_column,
    condition_col = condition_col,
    control_label = control_label,
    estim_label = estim_label
  )
  abundance_order <- order(
    -abundance_table$log2_enrichment,
    seq_len(nrow(abundance_table))
  )
  abundance_plot_data <- abundance_table
  # ANALYSIS_OK[R005]: factor levels reorder plotted abundance categories without dropping rows.
  abundance_plot_data$cluster <- factor(
    abundance_plot_data$cluster,
    levels = abundance_plot_data$cluster[abundance_order]
  )
  abundance_plot_data$direction <- factor(
    abundance_plot_data$direction,
    levels = c("Enriched in E-Stim", "Depleted in E-Stim", "Not significant")
  )
  abundance_plot_data$significance <- ifelse(
    # ANALYSIS_OK[R002]: fixed adjusted-p-value cutoff preserves the audited significance annotation.
    abundance_plot_data$padj < 0.05,
    # ANALYSIS_OK[R002]: fixed star cutoffs preserve the audited publication significance labels.
    ifelse(
      abundance_plot_data$padj < 0.001,
      "***",
      ifelse(abundance_plot_data$padj < 0.01, "**", "*")
    ),
    ""
  )
  abundance_y_span <- diff(range(
    abundance_plot_data$log2_enrichment,
    na.rm = TRUE
  ))
  abundance_label_offset <- max(0.12, 0.05 * abundance_y_span)
  abundance_plot_data$label_y <- abundance_plot_data$log2_enrichment +
    ifelse(
      abundance_plot_data$log2_enrichment >= 0,
      abundance_label_offset,
      -abundance_label_offset
    )
  abundance_plot_data$label_vjust <- ifelse(
    abundance_plot_data$log2_enrichment >= 0,
    0,
    1
  )
  abundance_plot <- ggplot2::ggplot(
    abundance_plot_data,
    ggplot2::aes(
      x = .data[["cluster"]],
      y = .data[["log2_enrichment"]],
      fill = .data[["direction"]]
    )
  ) +
    ggplot2::geom_col(color = "black", width = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(
        y = .data[["label_y"]],
        label = .data[["significance"]],
        vjust = .data[["label_vjust"]]
      ),
      size = 4,
      show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::scale_fill_manual(
      values = c(
        "Enriched in E-Stim" = unname(config$palettes$analysis[["high"]]),
        "Depleted in E-Stim" = unname(config$palettes$analysis[["low"]]),
        "Not significant" = unname(config$palettes$analysis[["mid"]])
      ),
      labels = c(
        "Enriched in E-Stim" = sprintf("Enriched in %s", estim_display_label),
        "Depleted in E-Stim" = sprintf("Depleted in %s", estim_display_label),
        "Not significant" = "Not significant"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.12, 0.18))
    ) +
    ggplot2::labs(
      x = "MG-selected cluster",
      y = sprintf("CLR log2 enrichment %s", config$conditions$contrast_display),
      fill = NULL,
      subtitle = "Pooled cell-level Fisher/CLR summary; descriptive relative to Mouse x Condition DE unit."
    ) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(sprintf(
      "MG-selected cluster abundance; %d PCs; res %s",
      mg_dims,
      mg_resolution
    ))

  sample_props <- compute_sample_cluster_proportions(
    sobj = sobj,
    cluster_col = cluster_column,
    mouse_col = "Mouse",
    condition_col = condition_col,
    control_label = control_label,
    estim_label = estim_label
  )
  randomization_table <- test_cluster_proportion_randomization(
    sample_props,
    control_label = control_label,
    estim_label = estim_label
  )
  proportion_data <- sample_props
  proportion_cluster_levels <- unique(as.character(proportion_data$cluster))
  if (all(grepl("^-?[0-9]+$", proportion_cluster_levels))) {
    proportion_cluster_levels <- as.character(sort(as.integer(
      proportion_cluster_levels
    )))
  } else {
    proportion_cluster_levels <- sort(
      proportion_cluster_levels,
      method = "radix"
    )
  }
  proportion_data$condition <- factor(
    as.character(proportion_data$condition),
    levels = c(control_label, estim_label)
  )
  proportion_data$cluster <- factor(
    as.character(proportion_data$cluster),
    levels = proportion_cluster_levels
  )
  proportion_data$mouse_role <- factor(
    as.character(proportion_data$mouse_role),
    levels = c("paired", "estim_only", "control_only")
  )
  paired_data <- proportion_data[
    proportion_data$mouse_role == "paired",
    ,
    drop = FALSE
  ]
  proportion_plot <- ggplot2::ggplot(
    proportion_data,
    ggplot2::aes(x = .data[["condition"]], y = .data[["proportion"]])
  ) +
    ggplot2::geom_line(
      data = paired_data,
      ggplot2::aes(group = .data[["mouse"]]),
      color = "grey60",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[["mouse"]], shape = .data[["mouse_role"]]),
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
    ggplot2::ggtitle(sprintf(
      "MG-selected cluster proportion by mouse; %d PCs; res %s",
      mg_dims,
      mg_resolution
    ))

  module_scores <- compute_cluster_module_scores(
    sobj,
    cluster_column,
    cell_type_marker_genes,
    assay = assay,
    slot = module_score_layer,
    seed = seed
  )
  module_scores_out <- module_scores
  rownames(module_scores_out) <- unname(
    cell_type_marker_labels[rownames(module_scores_out)]
  )
  module_scores_out <- data.frame(
    cell_type = rownames(module_scores_out),
    module_scores_out,
    check.names = FALSE
  )
  p27_enrichment <- compute_cluster_p27_enrichment(
    sobj,
    cluster_column,
    layer = expression_layer,
    assay = assay,
    mouse_col = "Mouse",
    condition_col = condition_col,
    n_perm = n_perm,
    seed = seed
  )

  branch_marker_paths <- write_curated_marker_heatmap(
    sobj,
    cluster_column,
    expression_layer,
    stems$marker,
    width = 10,
    height = 9,
    cluster_cells = TRUE
  )
  branch_module_paths <- write_module_p27_heatmap(
    module_scores,
    p27_enrichment,
    stems$module,
    width = 8,
    height = 6
  )
  heatmap_png_paths <- c(
    heatmap_png_paths,
    branch_marker_paths[["png"]],
    branch_module_paths[["png"]]
  )

  dir.create(mg_figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mg_table_dir, recursive = TRUE, showWarnings = FALSE)
  save_publication_plot(
    cluster_plot,
    stems$cluster,
    width = 5.5,
    height = 5.0,
    notebook_basename = paste0(basename(stems$cluster), ".png")
  )
  save_publication_plot(
    condition_plot,
    stems$condition,
    width = 5.5,
    height = 5.0,
    notebook_basename = paste0(basename(stems$condition), ".png")
  )
  save_publication_plot(
    feature_plot,
    stems$feature,
    width = 10.5,
    height = 9.0,
    notebook_basename = paste0(basename(stems$feature), ".png")
  )
  save_publication_plot(
    coexpression_plot,
    stems$coexpression,
    width = 6.0,
    height = 5.0,
    notebook_basename = paste0(basename(stems$coexpression), ".png")
  )
  save_publication_plot(
    abundance_plot,
    stems$abundance,
    width = max(6.5, 0.45 * nrow(abundance_table)),
    height = 4.5,
    notebook_basename = paste0(basename(stems$abundance), ".png")
  )
  save_publication_plot(
    proportion_plot,
    stems$proportion,
    width = 8,
    height = 6,
    notebook_basename = paste0(basename(stems$proportion), ".png")
  )
  utils::write.table(
    abundance_table,
    file.path(mg_table_dir, paste0(basename(stems$abundance), ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  utils::write.table(
    randomization_table,
    file.path(mg_table_dir, paste0(basename(stems$randomization), ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  utils::write.table(
    sample_props,
    file.path(mg_table_dir, paste0(basename(stems$sample_props), ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  utils::write.table(
    module_scores_out,
    file.path(
      annotation_table_dir,
      paste0(basename(stems$module), "_module_scores.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    p27_enrichment,
    file.path(
      annotation_table_dir,
      paste0(basename(stems$module), "_p27_enrichment.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# Replace each heatmap notebook image with a hash-verified regular-file sibling.
for (png_path in heatmap_png_paths) {
  destination <- file.path(notebook_figure_dir, basename(png_path))
  destination_link <- Sys.readlink(destination)
  if (
    length(destination_link) == 1L &&
      !is.na(destination_link) &&
      nzchar(destination_link)
  ) {
    stop(
      "Refusing to replace symlinked notebook figure: ",
      destination,
      call. = FALSE
    )
  }
  destination_info <- file.info(destination)
  if (
    !file.exists(destination) ||
      is.na(destination_info$isdir) ||
      isTRUE(destination_info$isdir)
  ) {
    stop(
      "Notebook figure destination must be an existing regular file: ",
      destination,
      call. = FALSE
    )
  }
  temporary <- tempfile(
    pattern = paste0(".", basename(destination), "."),
    tmpdir = dirname(destination)
  )
  on.exit(unlink(temporary), add = TRUE)
  if (
    !file.copy(png_path, temporary, overwrite = FALSE, copy.date = TRUE) ||
      !file.exists(temporary) ||
      isTRUE(file.info(temporary)$isdir)
  ) {
    stop(
      "Failed to create temporary notebook figure: ",
      destination,
      call. = FALSE
    )
  }
  source_hash <- digest::digest(png_path, algo = "sha256", file = TRUE)
  temporary_hash <- digest::digest(temporary, algo = "sha256", file = TRUE)
  if (!identical(source_hash, temporary_hash)) {
    stop(
      "Temporary notebook figure hash mismatch: ",
      destination,
      call. = FALSE
    )
  }
  source_dimensions <- magick::image_info(
    magick::image_read(png_path)
  )[1L, c("width", "height")]
  temporary_dimensions <- magick::image_info(
    magick::image_read(temporary)
  )[1L, c("width", "height")]
  if (!identical(source_dimensions, temporary_dimensions)) {
    stop(
      "Temporary notebook figure dimensions mismatch: ",
      destination,
      call. = FALSE
    )
  }
  unlink(destination)
  if (!file.rename(temporary, destination)) {
    stop("Failed to replace notebook figure: ", destination, call. = FALSE)
  }
  destination_hash <- digest::digest(destination, algo = "sha256", file = TRUE)
  if (!identical(source_hash, destination_hash)) {
    stop(
      "Notebook figure hash mismatch after replacement: ",
      destination,
      call. = FALSE
    )
  }
}

invisible(output_paths)
