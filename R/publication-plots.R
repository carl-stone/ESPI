#' Save a publication plot as PNG and PDF and copy the PNG to the notebook
#' @param plot A ggplot object.
#' @param output_stem Full output path without a file extension.
#' @param width,height Figure size in inches.
#' @param dpi PNG resolution.
#' @return Named PNG, PDF, and notebook paths, invisibly.
#' @export
save_publication_plot <- function(plot, output_stem, width, height, dpi = 300) {
  paths <- output_path(paste0(output_stem, c(".png", ".pdf")))
  ggplot2::ggsave(paths[[1]], plot, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(paths[[2]], plot, width = width, height = height)
  notebook <- copy_notebook_figure(paths[[1]])
  invisible(c(png = paths[[1]], pdf = paths[[2]], notebook = notebook))
}

#' Copy a figure to the notebook
#'
#' The caller chooses what to copy; this does not parse the notebook or hash images.
#' New destinations are allowed. Existing files require overwrite opt-in.
#' @param source Figure to copy.
#' @param destination Destination path, normally the same basename in notebook/figures.
#' @return The destination, invisibly.
#' @export
copy_notebook_figure <- function(
  source,
  destination = here::here("notebook", "figures", basename(source))
) {
  destination <- output_path(destination)
  if (!file.copy(source, destination, overwrite = TRUE)) {
    stop("Could not copy figure to ", destination)
  }
  invisible(destination)
}

#' Write the curated per-cell marker heatmap as PNG and PDF.
#'
#' @param sobj Seurat object containing marker genes and cluster metadata.
#' @param cluster_column Metadata column containing cluster labels.
#' @param layer Assay layer to plot from the default assay.
#' @param output_stem Full output path without a file extension.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param cluster_cells Whether to cluster cells within cluster slices and
#'   hierarchically reorder the slices.
#' @param sequential_clusters Keep cluster blocks in sorted cluster-ID order,
#'   preserving input cell order within blocks. Disables all column clustering
#'   and overrides `cluster_cells` when TRUE.
#'
#' @return Named PNG and PDF paths, invisibly.
#' @export
write_curated_marker_heatmap <- function(
  sobj,
  cluster_column,
  layer,
  output_stem,
  width,
  height,
  cluster_cells = FALSE,
  sequential_clusters = FALSE
) {
  if (isTRUE(sequential_clusters)) {
    cluster_cells <- FALSE
  }
  if (!cluster_column %in% colnames(sobj@meta.data)) {
    stop("Missing cluster metadata column: ", cluster_column, call. = FALSE)
  }
  assay <- SeuratObject::DefaultAssay(sobj)
  if (!layer %in% SeuratObject::Layers(sobj[[assay]])) {
    stop(
      "Missing expression layer '",
      layer,
      "' in assay ",
      assay,
      call. = FALSE
    )
  }

  marker_table <- stack(cell_type_marker_genes)
  colnames(marker_table) <- c("gene", "cell_type")
  marker_table$cell_type_label <- unname(cell_type_marker_labels[
    marker_table$cell_type
  ])
  missing_markers <- setdiff(marker_table$gene, rownames(sobj))
  if (length(missing_markers) > 0L) {
    stop(
      "Marker gene(s) missing from the Seurat object: ",
      paste(missing_markers, collapse = ", "),
      call. = FALSE
    )
  }

  cluster_values <- as.character(sobj@meta.data[[cluster_column]])
  cluster_levels <- .sort_cluster_labels(cluster_values)
  # Assign colors by cluster ID before any marker-dependent reordering.
  cluster_colors <- stats::setNames(
    grDevices::hcl.colors(length(cluster_levels), palette = "Temps"),
    paste("Cluster", cluster_levels)
  )
  marker_expression <- SeuratObject::GetAssayData(
    sobj,
    assay = assay,
    layer = layer
  )[marker_table$gene, , drop = FALSE]
  if (!identical(colnames(marker_expression), rownames(sobj@meta.data))) {
    stop(
      "Expression matrix columns do not match Seurat metadata rows.",
      call. = FALSE
    )
  }

  scaled_expression <- t(scale(t(as.matrix(marker_expression))))
  scaled_expression[is.na(scaled_expression)] <- 0
  z_score_limit <- 2
  scaled_expression[scaled_expression > z_score_limit] <- z_score_limit
  scaled_expression[scaled_expression < -z_score_limit] <- -z_score_limit
  cluster_dendrogram <- NULL
  if (!isTRUE(cluster_cells) && !isTRUE(sequential_clusters)) {
    cluster_means <- vapply(
      cluster_levels,
      function(cluster_value) {
        rowMeans(scaled_expression[,
          cluster_values == cluster_value,
          drop = FALSE
        ])
      },
      numeric(nrow(scaled_expression))
    )
    cluster_dendrogram <- stats::as.dendrogram(stats::hclust(stats::dist(t(
      cluster_means
    ))))
    cluster_levels <- labels(cluster_dendrogram)
  }
  cell_cluster_labels <- factor(
    paste("Cluster", cluster_values),
    levels = paste("Cluster", cluster_levels)
  )
  cell_type_groups <- factor(
    marker_table$cell_type_label,
    levels = unname(cell_type_marker_labels)
  )
  row_annotation <- ComplexHeatmap::rowAnnotation(
    `Cell type` = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = NA, col = NA),
      labels = levels(cell_type_groups),
      labels_gp = grid::gpar(fontsize = 8),
      labels_rot = 0,
      labels_just = "right",
      labels_offset = grid::unit(0.98, "npc"),
      width = grid::unit(1.15, "in")
    ),
    `Cell type marker divider` = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = "black", col = "black", lwd = 0.5),
      width = grid::unit(0.5, "mm")
    ),
    Gene = ComplexHeatmap::anno_text(
      marker_table$gene,
      gp = grid::gpar(fontsize = 8),
      just = "right",
      location = grid::unit(1, "npc"),
      width = grid::unit(0.5, "in")
    ),
    show_annotation_name = FALSE
  )
  column_annotation <- if (is.null(cluster_dendrogram)) {
    ComplexHeatmap::HeatmapAnnotation(
      Cluster = ComplexHeatmap::anno_block(
        height = grid::unit(4, "mm"),
        panel_fun = function(index, nm) {
          cluster_label <- nm
          grid::grid.rect(
            gp = grid::gpar(fill = cluster_colors[[cluster_label]], col = NA)
          )
          grid::grid.text(
            sub("^Cluster ", "", cluster_label),
            gp = grid::gpar(col = "black", fontsize = 8)
          )
        }
      ),
      show_annotation_name = FALSE
    )
  } else {
    ComplexHeatmap::HeatmapAnnotation(
      Cluster_dendrogram = ComplexHeatmap::anno_empty(
        border = FALSE,
        height = grid::unit(12, "mm")
      ),
      Cluster = ComplexHeatmap::anno_block(
        height = grid::unit(4, "mm"),
        panel_fun = function(index, nm) {
          cluster_label <- nm
          grid::grid.rect(
            gp = grid::gpar(fill = cluster_colors[[cluster_label]], col = NA)
          )
          grid::grid.text(
            sub("^Cluster ", "", cluster_label),
            gp = grid::gpar(col = "black", fontsize = 8)
          )
        }
      ),
      show_annotation_name = FALSE
    )
  }
  heatmap_arguments <- list(
    matrix = scaled_expression,
    name = "Row z-score",
    col = circlize::colorRamp2(
      c(-2, 0, 2),
      c(palette_dotplot_pair[1], "white", palette_dotplot_pair[2])
    ),
    left_annotation = row_annotation,
    top_annotation = column_annotation,
    row_split = cell_type_groups,
    cluster_rows = FALSE,
    cluster_columns = isTRUE(cluster_cells),
    column_split = cell_cluster_labels,
    show_column_names = FALSE,
    show_column_dend = isTRUE(cluster_cells),
    cluster_column_slices = isTRUE(cluster_cells),
    show_row_names = FALSE,
    row_title = NULL,
    column_title = NULL,
    use_raster = TRUE
  )
  heatmap <- do.call(ComplexHeatmap::Heatmap, heatmap_arguments)

  paths <- output_path(paste0(output_stem, c(".png", ".pdf")))
  png_path <- paths[[1]]
  pdf_path <- paths[[2]]
  grDevices::png(
    png_path,
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  if (is.null(cluster_dendrogram)) {
    ComplexHeatmap::draw(
      heatmap,
      heatmap_legend_side = "right",
      annotation_legend_side = "right"
    )
  } else {
    .draw_curated_marker_heatmap(heatmap, cluster_dendrogram)
  }
  grDevices::dev.off()
  grDevices::pdf(pdf_path, width = width, height = height)
  if (is.null(cluster_dendrogram)) {
    ComplexHeatmap::draw(
      heatmap,
      heatmap_legend_side = "right",
      annotation_legend_side = "right"
    )
  } else {
    .draw_curated_marker_heatmap(heatmap, cluster_dendrogram)
  }
  grDevices::dev.off()
  invisible(c(png = png_path, pdf = pdf_path))
}

.draw_curated_marker_heatmap <- function(heatmap, cluster_dendrogram) {
  drawn_heatmap <- ComplexHeatmap::draw(
    heatmap,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  slice_count <- length(labels(cluster_dendrogram))
  slice_centers <- vapply(
    seq_len(slice_count),
    function(slice_index) {
      grid::seekViewport(sprintf(
        "annotation_Cluster_dendrogram_%s",
        slice_index
      ))
      center <- grid::deviceLoc(
        x = grid::unit(0.5, "npc"),
        y = grid::unit(0, "npc")
      )
      grid::convertX(center$x, "in", valueOnly = TRUE)
    },
    numeric(1)
  )
  grid::seekViewport("annotation_Cluster_dendrogram_1")
  lower_left <- grid::deviceLoc(
    x = grid::unit(0, "npc"),
    y = grid::unit(0, "npc")
  )
  grid::seekViewport(sprintf("annotation_Cluster_dendrogram_%s", slice_count))
  upper_right <- grid::deviceLoc(
    x = grid::unit(1, "npc"),
    y = grid::unit(1, "npc")
  )
  left <- grid::convertX(lower_left$x, "in", valueOnly = TRUE)
  bottom <- grid::convertY(lower_left$y, "in", valueOnly = TRUE)
  right <- grid::convertX(upper_right$x, "in", valueOnly = TRUE)
  top <- grid::convertY(upper_right$y, "in", valueOnly = TRUE)
  leaf_positions <- (slice_centers - left) / (right - left)
  positioned_dendrogram <- ComplexHeatmap::adjust_dend_by_x(
    cluster_dendrogram,
    leaf_pos = leaf_positions
  )
  grid::seekViewport("global")
  grid::pushViewport(grid::viewport(
    x = grid::unit(left, "in"),
    y = grid::unit(bottom, "in"),
    width = grid::unit(right - left, "in"),
    height = grid::unit(top - bottom, "in"),
    just = c("left", "bottom"),
    xscale = c(0, 1),
    yscale = c(0, attr(cluster_dendrogram, "height"))
  ))
  grid::grid.draw(ComplexHeatmap::dendrogramGrob(
    positioned_dendrogram,
    facing = "bottom",
    gp = grid::gpar(col = "black", lwd = 0.5)
  ))
  grid::popViewport()
  invisible(drawn_heatmap)
}

#' Write the module-score and p27-enrichment heatmap as PNG and PDF.
#'
#' @param module_scores Numeric matrix returned by
#'   `compute_cluster_module_scores()`.
#' @param p27_enrichment Data frame returned by
#'   `compute_cluster_p27_enrichment()`.
#' @param output_stem Full output path without a file extension.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#'
#' @return Named PNG and PDF paths, invisibly.
#' @export
write_module_p27_heatmap <- function(
  module_scores,
  p27_enrichment,
  output_stem,
  width,
  height
) {
  module_matrix <- module_scores
  module_matrix_rownames <- unname(cell_type_marker_labels[rownames(
    module_matrix
  )])
  if (anyNA(module_matrix_rownames)) {
    stop(
      "module_scores row names do not match curated cell-type marker names.",
      call. = FALSE
    )
  }
  rownames(module_matrix) <- module_matrix_rownames
  module_z <- t(scale(t(module_matrix)))
  module_z[is.na(module_z)] <- 0
  z_score_limit <- 2
  module_z[module_z > z_score_limit] <- z_score_limit
  module_z[module_z < -z_score_limit] <- -z_score_limit
  p27_index <- match(colnames(module_z), p27_enrichment$cluster)
  if (!identical(p27_enrichment$cluster[p27_index], colnames(module_z))) {
    stop(
      "p27 enrichment rows do not align to module-score columns.",
      call. = FALSE
    )
  }
  p27_z <- p27_enrichment$z_score[p27_index]

  body_fun <- circlize::colorRamp2(
    c(-2, 0, 2),
    c(palette_dotplot_pair[1], "white", palette_dotplot_pair[2])
  )
  finite_p27 <- p27_z[is.finite(p27_z)]
  zlim <- if (length(finite_p27) > 0L) max(abs(finite_p27)) else NA_real_
  if (!is.finite(zlim) || zlim == 0) {
    zlim <- 1
  }
  p27_fun <- circlize::colorRamp2(
    c(-zlim, 0, zlim),
    c(palette_dotplot_pair[1], "white", palette_dotplot_pair[2])
  )
  p27_legend <- ComplexHeatmap::Legend(col_fun = p27_fun, title = "p27 z-score")
  heatmap <- ComplexHeatmap::Heatmap(
    module_z,
    name = "Module z-score",
    col = body_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      `p27 z-score` = ComplexHeatmap::anno_simple(
        p27_z,
        col = p27_fun,
        na_col = "grey85"
      ),
      annotation_name_gp = grid::gpar(fontsize = 8),
      show_legend = FALSE
    ),
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    column_title = NULL,
    use_raster = FALSE
  )

  paths <- output_path(paste0(output_stem, c(".png", ".pdf")))
  png_path <- paths[[1]]
  pdf_path <- paths[[2]]
  grDevices::png(
    png_path,
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  ComplexHeatmap::draw(
    heatmap,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    heatmap_legend_list = list(p27_legend)
  )
  grDevices::dev.off()
  grDevices::pdf(pdf_path, width = width, height = height)
  ComplexHeatmap::draw(
    heatmap,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    heatmap_legend_list = list(p27_legend)
  )
  grDevices::dev.off()
  invisible(c(png = png_path, pdf = pdf_path))
}
