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
  out_dir <- file.path(FIGURE_DIR, "preprocess")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot <- Seurat::DimPlot(sobj, reduction = umap, group.by = color_by)
  ggplot2::ggsave(
    file.path(out_dir, sprintf("%s__%s.png", umap, color_by)),
    plot,
    width = 5,
    height = 5
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("%s__%s.pdf", umap, color_by)),
    plot,
    width = 5,
    height = 5
  )

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
  out_dir <- file.path(FIGURE_DIR, "preprocess")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plot <- clustree::clustree(sobj, prefix = prefix) +
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
