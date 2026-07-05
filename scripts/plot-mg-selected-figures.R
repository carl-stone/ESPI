#!/usr/bin/env Rscript

# Generate representative mg-selected UMAP figures.
#
# Usage:
#   Rscript scripts/plot-mg-selected-figures.R \
#     [--input <clustered-seurat-object.rds>] \
#     [--branch-tag <branch tag>] \
#     [--elbow-n <positive integer>] \
#     [--dims <positive integer>] \
#     [--resolution <resolution string>] \
#     [--layer <assay layer>] \
#     [--feature-list <data/umap_feature_list.rda>]
#
# Outputs:
#   FIGURE_DIR/mg_selected/mg_selected_cluster_umap_<branch>_dims<dims>_res<resolution>.(png|pdf)
#   FIGURE_DIR/mg_selected/mg_selected_feature_umap_<layer>_<branch>_dims<dims>_res<resolution>.(png|pdf)
#   FIGURE_DIR/mg_selected/mg_selected_cluster_abundance_enrichment_<branch>_dims<dims>_res<resolution>.(png|pdf)
#   TABLE_DIR/mg_selected/mg_selected_cluster_abundance_enrichment_<branch>_dims<dims>_res<resolution>.tsv
#   notebook/figures/<png filename> symlinks for all PNG outputs.

suppressPackageStartupMessages({
  library(here)
})
here::i_am("scripts/plot-mg-selected-figures.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})
palette_dotplot_pair <- get(
  "palette_dotplot_pair",
  envir = asNamespace("ESPI"),
  inherits = FALSE
)
compute_cluster_abundance <- get(
  "compute_cluster_abundance",
  envir = asNamespace("ESPI"),
  inherits = FALSE
)
plot_clr_fisher_enrichment <- get(
  "plot_clr_fisher_enrichment",
  envir = asNamespace("ESPI"),
  inherits = FALSE
)

# ---- parameters ----

get_arg <- function(args, flag, default) {
  match_index <- match(flag, args)
  if (is.na(match_index)) {
    return(default)
  }
  if (
    match_index == length(args) || startsWith(args[[match_index + 1]], "--")
  ) {
    stop("Missing value for ", flag, call. = FALSE)
  }
  args[[match_index + 1]]
}

cli_args <- commandArgs(trailingOnly = TRUE)
allowed_flags <- c(
  "--input",
  "--branch-tag",
  "--elbow-n",
  "--dims",
  "--resolution",
  "--layer",
  "--feature-list"
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

branch_tag <- get_arg(
  cli_args,
  "--branch-tag",
  "pflog_mg_selected_no_filter_cc"
)
if (
  !is.character(branch_tag) ||
    length(branch_tag) != 1L ||
    is.na(branch_tag) ||
    !nzchar(branch_tag) ||
    !grepl("^[A-Za-z0-9_]+$", branch_tag)
) {
  stop("--branch-tag must be a safe non-empty branch tag.", call. = FALSE)
}

elbow_n <- as.integer(get_arg(cli_args, "--elbow-n", "20"))
if (length(elbow_n) != 1L || is.na(elbow_n) || elbow_n <= 0) {
  stop("--elbow-n must be a positive integer.", call. = FALSE)
}

dims <- as.integer(get_arg(cli_args, "--dims", "30"))
if (length(dims) != 1L || is.na(dims) || dims <= 0) {
  stop("--dims must be a positive integer.", call. = FALSE)
}

resolution <- get_arg(cli_args, "--resolution", "0.3")
if (
  !is.character(resolution) ||
    length(resolution) != 1L ||
    is.na(resolution) ||
    !nzchar(resolution)
) {
  stop("--resolution must be a non-empty string.", call. = FALSE)
}
resolution_number <- suppressWarnings(as.numeric(resolution))
if (
  is.na(resolution_number) ||
    !is.finite(resolution_number) ||
    resolution_number <= 0
) {
  stop("--resolution must be a positive number string.", call. = FALSE)
}

expression_layer <- get_arg(cli_args, "--layer", "pflog")
if (
  !is.character(expression_layer) ||
    length(expression_layer) != 1L ||
    is.na(expression_layer) ||
    !nzchar(expression_layer)
) {
  stop("--layer must be a non-empty assay layer name.", call. = FALSE)
}

feature_list_path <- get_arg(
  cli_args,
  "--feature-list",
  here::here("data", "umap_feature_list.rda")
)

input_path <- get_arg(
  cli_args,
  "--input",
  file.path(
    CURRENT_OBJECT_DIR,
    sprintf("cluster_%s_elbow%d.rds", branch_tag, elbow_n)
  )
)

# ---- helpers ----

filename_tag <- function(value) {
  gsub("[^A-Za-z0-9_.-]", "_", value)
}

load_umap_feature_list <- function(path) {
  if (!file.exists(path)) {
    stop("Feature list file does not exist: ", path, call. = FALSE)
  }
  feature_env <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = feature_env)
  if ("umap_feature_list" %in% loaded_names) {
    features <- get("umap_feature_list", envir = feature_env)
  } else if (length(loaded_names) == 1L) {
    features <- get(loaded_names[[1L]], envir = feature_env)
  } else {
    stop(
      "Feature list file must contain umap_feature_list or exactly one object: ",
      path,
      call. = FALSE
    )
  }
  if (
    !is.character(features) || is.matrix(features) || is.data.frame(features)
  ) {
    stop("UMAP feature list must be a character vector.", call. = FALSE)
  }
  features <- unname(features)
  if (length(features) != 9L) {
    stop(
      "UMAP feature list must contain exactly 9 genes; found ",
      length(features),
      call. = FALSE
    )
  }
  if (anyNA(features) || any(!nzchar(features))) {
    stop(
      "UMAP feature list contains missing or empty gene names.",
      call. = FALSE
    )
  }
  duplicated_features <- unique(features[duplicated(features)])
  if (length(duplicated_features) > 0L) {
    stop(
      "UMAP feature list contains duplicated gene(s): ",
      paste(duplicated_features, collapse = ", "),
      call. = FALSE
    )
  }
  features
}

link_notebook_png <- function(png_path) {
  notebook_figure_dir <- here::here("notebook", "figures")
  dir.create(notebook_figure_dir, recursive = TRUE, showWarnings = FALSE)
  notebook_png_path <- file.path(notebook_figure_dir, basename(png_path))
  if (
    file.exists(notebook_png_path) || nzchar(Sys.readlink(notebook_png_path))
  ) {
    unlink(notebook_png_path)
  }
  link_created <- file.symlink(png_path, notebook_png_path)
  if (!isTRUE(link_created)) {
    stop("Failed to link notebook figure: ", notebook_png_path, call. = FALSE)
  }
  notebook_png_path
}

cluster_levels_for_labels <- function(values) {
  unique_values <- unique(as.character(values))
  if (all(grepl("^-?[0-9]+$", unique_values))) {
    as.character(sort(as.integer(unique_values)))
  } else {
    sort(unique_values, method = "radix")
  }
}

feature_umap_plot <- function(sobj, features, reduction, assay, layer) {
  embeddings <- SeuratObject::Embeddings(sobj, reduction = reduction)
  if (ncol(embeddings) < 2L) {
    stop(
      "UMAP reduction has fewer than two dimensions: ",
      reduction,
      call. = FALSE
    )
  }
  if (!setequal(rownames(embeddings), rownames(sobj@meta.data))) {
    stop(
      "UMAP reduction cells do not match Seurat metadata rows: ",
      reduction,
      call. = FALSE
    )
  }
  embeddings <- embeddings[rownames(sobj@meta.data), , drop = FALSE]
  expression <- SeuratObject::GetAssayData(
    sobj,
    assay = assay,
    layer = layer
  )[features, rownames(sobj@meta.data), drop = FALSE]
  if (!identical(colnames(expression), rownames(sobj@meta.data))) {
    stop(
      "Expression matrix columns do not match Seurat metadata rows.",
      call. = FALSE
    )
  }

  umap_x_range <- range(embeddings[, 1L])
  umap_y_range <- range(embeddings[, 2L])
  umap_span <- max(diff(umap_x_range), diff(umap_y_range))
  umap_x_center <- mean(umap_x_range)
  umap_y_center <- mean(umap_y_range)
  umap_x_limits <- umap_x_center + c(-0.5, 0.5) * umap_span
  umap_y_limits <- umap_y_center + c(-0.5, 0.5) * umap_span

  plots <- lapply(features, function(feature) {
    feature_expression <- as.numeric(expression[feature, ])
    expression_range <- range(feature_expression, na.rm = TRUE)
    if (!all(is.finite(expression_range))) {
      stop(
        "Expression values are not finite for feature: ",
        feature,
        call. = FALSE
      )
    }
    if (expression_range[[2L]] > expression_range[[1L]]) {
      scaled_expression <- (feature_expression - expression_range[[1L]]) /
        (expression_range[[2L]] - expression_range[[1L]])
    } else {
      scaled_expression <- rep(0, length(feature_expression))
    }
    plot_data <- data.frame(
      UMAP_1 = embeddings[, 1L],
      UMAP_2 = embeddings[, 2L],
      scaled_expression = scaled_expression,
      stringsAsFactors = FALSE
    )
    plot_data <- plot_data[order(plot_data$scaled_expression), ]
    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data[["UMAP_1"]],
        y = .data[["UMAP_2"]],
        color = .data[["scaled_expression"]]
      )
    ) +
      ggplot2::geom_point(size = 0.25, stroke = 0) +
      ggplot2::scale_color_gradient(
        low = "grey85",
        high = palette_dotplot_pair[[2L]],
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
  patchwork::wrap_plots(plots, ncol = 3L) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")
}

# ---- validation ----

if (!file.exists(input_path)) {
  stop("Input Seurat object does not exist: ", input_path, call. = FALSE)
}
features <- load_umap_feature_list(feature_list_path)

# ---- work ----

sobj <- readRDS(input_path)
object_branch_tag <- sobj@misc$clustering$branch_tag
if (
  is.character(object_branch_tag) &&
    length(object_branch_tag) == 1L &&
    nzchar(object_branch_tag) &&
    !identical(object_branch_tag, branch_tag)
) {
  stop(
    "Input object branch tag mismatch: expected ",
    branch_tag,
    ", found ",
    object_branch_tag,
    call. = FALSE
  )
}

reduction <- sprintf("umap_%s_dims%d", branch_tag, dims)
cluster_column <- sprintf(
  "cluster_%s_dims%d_res%s",
  branch_tag,
  dims,
  resolution
)
if (!reduction %in% names(sobj@reductions)) {
  stop("Missing UMAP reduction: ", reduction, call. = FALSE)
}
if (!cluster_column %in% colnames(sobj@meta.data)) {
  stop("Missing cluster metadata column: ", cluster_column, call. = FALSE)
}
if (anyNA(sobj@meta.data[[cluster_column]])) {
  stop(
    "Cluster metadata column contains NA values: ",
    cluster_column,
    call. = FALSE
  )
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
missing_features <- setdiff(features, rownames(sobj))
if (length(missing_features) > 0L) {
  stop(
    "UMAP feature gene(s) missing from the Seurat object: ",
    paste(missing_features, collapse = ", "),
    call. = FALSE
  )
}

cluster_values <- as.character(sobj@meta.data[[cluster_column]])
cluster_levels <- cluster_levels_for_labels(cluster_values)
SeuratObject::Idents(sobj) <- factor(cluster_values, levels = cluster_levels)

cluster_plot <- Seurat::DimPlot(
  sobj,
  reduction = reduction,
  group.by = cluster_column,
  label = TRUE,
  pt.size = 0.25
) +
  ggplot2::ggtitle(sprintf(
    "MG-selected PFlog; %d PCs; res %s",
    dims,
    resolution
  )) +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2")

abundance_table <- compute_cluster_abundance(
  sobj = sobj,
  cluster_col = cluster_column
)
abundance_plot <- plot_clr_fisher_enrichment(abundance_table) +
  ggplot2::ggtitle(sprintf(
    "MG-selected cluster abundance; %d PCs; res %s",
    dims,
    resolution
  ))

feature_plot <- feature_umap_plot(
  sobj = sobj,
  features = features,
  reduction = reduction,
  assay = assay,
  layer = expression_layer
) +
  patchwork::plot_annotation(
    title = sprintf("%s layer marker expression", expression_layer)
  )

# ---- output ----

out_dir <- file.path(FIGURE_DIR, "mg_selected")
table_dir <- file.path(TABLE_DIR, "mg_selected")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
resolution_tag <- filename_tag(resolution)
cluster_out_tag <- sprintf(
  "mg_selected_cluster_umap_%s_dims%d_res%s",
  branch_tag,
  dims,
  resolution_tag
)
feature_out_tag <- sprintf(
  "mg_selected_feature_umap_%s_%s_dims%d_res%s",
  filename_tag(expression_layer),
  branch_tag,
  dims,
  resolution_tag
)
abundance_out_tag <- sprintf(
  "mg_selected_cluster_abundance_enrichment_%s_dims%d_res%s",
  branch_tag,
  dims,
  resolution_tag
)
cluster_png_path <- file.path(out_dir, sprintf("%s.png", cluster_out_tag))
cluster_pdf_path <- file.path(out_dir, sprintf("%s.pdf", cluster_out_tag))
feature_png_path <- file.path(out_dir, sprintf("%s.png", feature_out_tag))
feature_pdf_path <- file.path(out_dir, sprintf("%s.pdf", feature_out_tag))
abundance_png_path <- file.path(out_dir, sprintf("%s.png", abundance_out_tag))
abundance_pdf_path <- file.path(out_dir, sprintf("%s.pdf", abundance_out_tag))
abundance_tsv_path <- file.path(table_dir, sprintf("%s.tsv", abundance_out_tag))

ggplot2::ggsave(
  cluster_png_path,
  cluster_plot,
  width = 5.5,
  height = 5.0,
  bg = "white"
)
ggplot2::ggsave(
  cluster_pdf_path,
  cluster_plot,
  width = 5.5,
  height = 5.0,
  bg = "white"
)
ggplot2::ggsave(
  feature_png_path,
  feature_plot,
  width = 10.5,
  height = 9.0,
  bg = "white"
)
ggplot2::ggsave(
  feature_pdf_path,
  feature_plot,
  width = 10.5,
  height = 9.0,
  bg = "white"
)
ggplot2::ggsave(
  abundance_png_path,
  abundance_plot,
  width = max(6.5, 0.45 * nrow(abundance_table)),
  height = 4.5,
  bg = "white"
)
ggplot2::ggsave(
  abundance_pdf_path,
  abundance_plot,
  width = max(6.5, 0.45 * nrow(abundance_table)),
  height = 4.5,
  bg = "white"
)
utils::write.table(
  abundance_table,
  file = abundance_tsv_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cluster_notebook_path <- link_notebook_png(cluster_png_path)
feature_notebook_path <- link_notebook_png(feature_png_path)
abundance_notebook_path <- link_notebook_png(abundance_png_path)

message("Wrote mg-selected cluster UMAP PNG: ", cluster_png_path)
message("Wrote mg-selected cluster UMAP PDF: ", cluster_pdf_path)
message("Wrote mg-selected feature UMAP PNG: ", feature_png_path)
message("Wrote mg-selected feature UMAP PDF: ", feature_pdf_path)
message("Wrote mg-selected cluster abundance TSV: ", abundance_tsv_path)
message("Wrote mg-selected cluster abundance PNG: ", abundance_png_path)
message("Wrote mg-selected cluster abundance PDF: ", abundance_pdf_path)
message("Linked notebook figure: ", cluster_notebook_path)
message("Linked notebook figure: ", feature_notebook_path)
message("Linked notebook figure: ", abundance_notebook_path)
