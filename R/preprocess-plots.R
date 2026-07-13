#' Save QC metric violin plots by sample.
#'
#' Draws one violin plot each for `nFeature_RNA`, `nCount_RNA`, `percent.mt`,
#' and `percent.ribo`, groups cells by `sample_id`, arranges the plots in a
#' 2x2 grid, and writes branch-tagged PNG and PDF files to `FIGURE_DIR/preprocess`.
#'
#' @param sobj Seurat object with a `sample_id` metadata column and the QC
#'   metadata columns `nFeature_RNA`, `nCount_RNA`, `percent.mt`, and
#'   `percent.ribo`.
#'
#' @return `NULL`, invisibly.
#' @export
splot_qc_metrics_violin <- function(sobj) {
  required_cols <- c(
    "sample_id",
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.ribo"
  )
  missing_cols <- setdiff(required_cols, colnames(sobj@meta.data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  out_dir <- file.path(FIGURE_DIR, "preprocess")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  norm <- sobj@misc$preprocessing$normalization
  cc_tag <- if (isTRUE(sobj@misc$preprocessing$filtered_cell_cycle)) {
    "filter-cc"
  } else {
    "no-filter-cc"
  }
  branch_tag <- sprintf("%s_%s", norm, cc_tag)

  metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
  plot <- Seurat::VlnPlot(
    sobj,
    features = metrics,
    group.by = "sample_id",
    ncol = 2
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("qc_metrics_violin_%s.png", branch_tag)),
    plot,
    width = 8,
    height = 8
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("qc_metrics_violin_%s.pdf", branch_tag)),
    plot,
    width = 8,
    height = 8
  )

  invisible(NULL)
}

#' Save the HVG mean-vs-variance scatter plot.
#'
#' Labels the first `n_top` genes from `VariableFeatures(sobj)` and writes
#' branch-tagged PNG and PDF files to `FIGURE_DIR/preprocess`.
#'
#' @param sobj Seurat object with `VariableFeatures` populated.
#' @param n_top Integer count of top variable features to label. Default 10.
#'
#' @return `NULL`, invisibly.
#' @export
splot_hvg_scatter <- function(sobj, n_top = 10) {
  stopifnot(length(n_top) == 1, is.numeric(n_top), is.finite(n_top), n_top >= 0)

  variable_features <- SeuratObject::VariableFeatures(sobj)
  if (length(variable_features) == 0) {
    stop("VariableFeatures(sobj) is empty.", call. = FALSE)
  }

  out_dir <- file.path(FIGURE_DIR, "preprocess")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  norm <- sobj@misc$preprocessing$normalization
  cc_tag <- if (isTRUE(sobj@misc$preprocessing$filtered_cell_cycle)) {
    "filter-cc"
  } else {
    "no-filter-cc"
  }
  branch_tag <- sprintf("%s_%s", norm, cc_tag)

  top_features <- head(variable_features, n_top)
  plot <- Seurat::VariableFeaturePlot(sobj)
  if (length(top_features) > 0) {
    plot <- Seurat::LabelPoints(
      plot = plot,
      points = top_features,
      repel = TRUE,
      xnudge = 0,
      ynudge = 0,
      max.overlaps = Inf
    )
  }

  ggplot2::ggsave(
    file.path(out_dir, sprintf("hvg_scatter_%s.png", branch_tag)),
    plot,
    width = 6,
    height = 5
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("hvg_scatter_%s.pdf", branch_tag)),
    plot,
    width = 6,
    height = 5
  )

  invisible(NULL)
}

#' Save PCA DimHeatmap plots.
#'
#' The output filename includes the normalization and cell-cycle filtering branch
#' tag from `sobj@misc$preprocessing`.
#'
#' @param sobj Seurat object with `pca` reduction and preprocessing metadata.
#'
#' @return `invisible(NULL)`.
#' @export
splot_dim_heatmap <- function(sobj) {
  norm <- sobj@misc$preprocessing$normalization
  cc_tag <- if (isTRUE(sobj@misc$preprocessing$filtered_cell_cycle)) {
    "filter-cc"
  } else {
    "no-filter-cc"
  }
  branch_tag <- sprintf("%s_%s", norm, cc_tag)
  out_dir <- file.path(FIGURE_DIR, "preprocess")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pca_source_layer <- sobj@misc$preprocessing$pca_source_layer
  if (
    !is.character(pca_source_layer) ||
      length(pca_source_layer) != 1 ||
      is.na(pca_source_layer) ||
      !nzchar(pca_source_layer)
  ) {
    stop("Missing PCA source layer in preprocessing metadata.", call. = FALSE)
  }

  plot <- Seurat::DimHeatmap(
    sobj,
    dims = 1:6,
    cells = 500,
    balanced = TRUE,
    fast = FALSE,
    combine = TRUE,
    slot = pca_source_layer,
    ncol = 2
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("dim_heatmap_%s.png", branch_tag)),
    plot,
    width = 8,
    height = 12
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("dim_heatmap_%s.pdf", branch_tag)),
    plot,
    width = 8,
    height = 12
  )

  invisible(NULL)
}

#' Save an ElbowPlot for the `pca` reduction.
#'
#' Filename includes the normalization and cell-cycle filtering branch tag.
#'
#' @param sobj Seurat object with `pca` reduction and preprocessing metadata.
#' @param n_pcs Integer number of PCs to plot. Default 50.
#'
#' @return `invisible(NULL)`.
#' @export
splot_elbow <- function(sobj, n_pcs = 50) {
  norm <- sobj@misc$preprocessing$normalization
  cc_tag <- if (isTRUE(sobj@misc$preprocessing$filtered_cell_cycle)) {
    "filter-cc"
  } else {
    "no-filter-cc"
  }
  branch_tag <- sprintf("%s_%s", norm, cc_tag)
  out_dir <- file.path(FIGURE_DIR, "preprocess")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot <- Seurat::ElbowPlot(sobj, ndims = n_pcs, reduction = "pca")
  ggplot2::ggsave(
    file.path(out_dir, sprintf("elbow_%s.png", branch_tag)),
    plot,
    width = 5,
    height = 3
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("elbow_%s.pdf", branch_tag)),
    plot,
    width = 5,
    height = 3
  )

  invisible(NULL)
}
