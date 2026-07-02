#' Save QC metric violin plots by sample.
#'
#' Draws one violin plot each for `nFeature_RNA`, `nCount_RNA`, `percent.mt`,
#' and `percent.ribo`, groups cells by `sample_id`, arranges the plots in a
#' 2x2 grid, and writes PNG and PDF files to `FIGURE_DIR/preprocess`.
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

  metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
  plot <- Seurat::VlnPlot(
    sobj,
    features = metrics,
    group.by = "sample_id",
    ncol = 2
  )
  ggplot2::ggsave(
    file.path(out_dir, "qc_metrics_violin.png"),
    plot,
    width = 8,
    height = 8
  )
  ggplot2::ggsave(
    file.path(out_dir, "qc_metrics_violin.pdf"),
    plot,
    width = 8,
    height = 8
  )

  invisible(NULL)
}

#' Save the HVG mean-vs-variance scatter plot.
#'
#' Labels the first `n_top` genes from `VariableFeatures(sobj)` and writes PNG
#' and PDF files to `FIGURE_DIR/preprocess`.
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

  top_features <- head(variable_features, n_top)
  plot <- Seurat::VariableFeaturePlot(sobj)
  if (length(top_features) > 0) {
    plot <- Seurat::LabelPoints(
      plot = plot,
      points = top_features,
      repel = TRUE,
      xnudge = 0,
      ynudge = 0,
      max.overlaps = 15
    )
  }

  ggplot2::ggsave(
    file.path(out_dir, "hvg_scatter.png"),
    plot,
    width = 6,
    height = 5
  )
  ggplot2::ggsave(
    file.path(out_dir, "hvg_scatter.pdf"),
    plot,
    width = 6,
    height = 5
  )

  invisible(NULL)
}

#' Save VizDimLoadings for the first n PCs of the `pca` reduction.
#'
#' The output filename includes the normalization branch tag from
#' `sobj@misc$preprocessing$normalization`.
#'
#' @param sobj Seurat object with `pca` reduction and preprocessing metadata.
#' @param n_pcs Integer number of PCs to plot. Default 30.
#'
#' @return `invisible(NULL)`.
#' @export
splot_viz_dim_loadings <- function(sobj, n_pcs = 30) {
  norm <- sobj@misc$preprocessing$normalization
  out_dir <- file.path(FIGURE_DIR, "preprocess")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot <- Seurat::VizDimLoadings(sobj, dims = 1:n_pcs, reduction = "pca")
  height <- max(6, ceiling(n_pcs / 2) * 1.2)
  ggplot2::ggsave(
    file.path(out_dir, sprintf("viz_dim_loadings_%s.png", norm)),
    plot,
    width = 8,
    height = height
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("viz_dim_loadings_%s.pdf", norm)),
    plot,
    width = 8,
    height = height
  )

  invisible(NULL)
}

#' Save an ElbowPlot for the `pca` reduction.
#'
#' Filename includes the normalization branch tag.
#'
#' @param sobj Seurat object with `pca` reduction and preprocessing metadata.
#' @param n_pcs Integer number of PCs to plot. Default 50.
#'
#' @return `invisible(NULL)`.
#' @export
splot_elbow <- function(sobj, n_pcs = 50) {
  norm <- sobj@misc$preprocessing$normalization
  out_dir <- file.path(FIGURE_DIR, "preprocess")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot <- Seurat::ElbowPlot(sobj, ndims = n_pcs, reduction = "pca")
  ggplot2::ggsave(
    file.path(out_dir, sprintf("elbow_%s.png", norm)),
    plot,
    width = 5,
    height = 3
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("elbow_%s.pdf", norm)),
    plot,
    width = 5,
    height = 3
  )

  invisible(NULL)
}

#' Save a DimPlot on a specified UMAP reduction colored by a metadata column.
#'
#' Filename encodes the UMAP reduction and coloring column.
#'
#' @param sobj Seurat object with the named UMAP reduction populated.
#' @param umap Character name of the UMAP reduction in `sobj@reductions`.
#' @param color_by Character name of a `sobj@meta.data` column.
#'
#' @return `invisible(NULL)`.
#' @export
splot_umap_by <- function(sobj, umap, color_by) {
  # 1. p <- Seurat::DimPlot(sobj, reduction = umap, group.by = color_by)
  # 2. save p as PNG and PDF at 5 in x 5 in to
  #    FIGURE_DIR/preprocess/<umap>__<color_by>.{png,pdf}
  # 3. invisible(NULL)
  invisible(NULL)
}

#' Save a clustree plot for a family of resolution-varying cluster columns.
#'
#' @param sobj Seurat object with the resolution-family cluster columns
#'   populated in `meta.data`, all sharing a common prefix.
#' @param prefix Character prefix shared by the cluster columns; clustree
#'   requires the shared prefix to walk resolution levels.
#' @param out_tag Character tag embedded in the output filename to distinguish
#'   dims choices.
#'
#' @return `invisible(NULL)`.
#' @export
splot_clustree <- function(sobj, prefix, out_tag) {
  # 1. p <- clustree::clustree(sobj, prefix = prefix)
  # 2. save p as PNG and PDF at 6 in x 6 in to
  #    FIGURE_DIR/preprocess/clustree_<out_tag>.{png,pdf}
  # 3. invisible(NULL)
  # Note for implementation: if clustree is not installed at runtime, wrap the
  # body in requireNamespace("clustree", quietly = TRUE), warn, and return
  # invisible(NULL) when unavailable.
  invisible(NULL)
}
