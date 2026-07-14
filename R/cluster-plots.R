#' Save a DimPlot on a specified UMAP reduction colored by a metadata column.
#'
#' Filename combines the UMAP reduction and coloring column.
#' Saves the PNG to the notebook figure directory as a symlink after writing
#' both plot files.
#'
#' @param sobj Seurat object with the named UMAP reduction populated.
#' @param umap Character name of the UMAP reduction in `sobj@reductions`.
#' @param color_by Character name of a `sobj@meta.data` column.
#'
#' @return `invisible(NULL)`.
#' @export
splot_umap_by <- function(sobj, umap, color_by) {
  out_dir <- file.path(FIGURE_DIR, "cluster")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  umap_tag <- gsub("[^A-Za-z0-9_-]", "_", umap)
  color_tag <- gsub("[^A-Za-z0-9_-]", "_", color_by)

  plot <- Seurat::DimPlot(
    sobj,
    reduction = umap,
    group.by = color_by,
    label = TRUE,
    pt.size = 0.25
  )
  png_path <- file.path(out_dir, sprintf("%s_by_%s.png", umap_tag, color_tag))
  ggplot2::ggsave(
    png_path,
    plot,
    width = 5,
    height = 5
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("%s_by_%s.pdf", umap_tag, color_tag)),
    plot,
    width = 5,
    height = 5
  )
  link_notebook_png(png_path)

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
#' @import ggraph
#' @export
splot_clustree <- function(sobj, prefix, out_tag) {
  out_dir <- file.path(FIGURE_DIR, "cluster")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  matching_cols <- colnames(sobj@meta.data)[startsWith(
    colnames(sobj@meta.data),
    prefix
  )]
  if (length(matching_cols) < MIN_CLUSTREE_RESOLUTION_COLUMNS) {
    stop(
      "Need at least two resolution columns for clustree prefix: ",
      prefix,
      call. = FALSE
    )
  }
  cluster_data <- sobj@meta.data[, matching_cols, drop = FALSE]
  plot <- clustree::clustree(cluster_data, prefix = prefix) +
    ggplot2::guides(edge_colour = "none")
  ggplot2::ggsave(
    file.path(out_dir, sprintf("clustree_%s.png", out_tag)),
    plot,
    width = 6,
    height = 6
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("clustree_%s.pdf", out_tag)),
    plot,
    width = 6,
    height = 6
  )

  invisible(NULL)
}
