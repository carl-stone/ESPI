#' Save a publication plot as PNG and PDF and mirror its PNG into the notebook.
#'
#' The notebook mirror is replaced through a temporary regular file. Existing
#' symlinks are rejected so a publication run cannot write through the notebook
#' tree into an external output directory.
#'
#' @param plot A ggplot object.
#' @param output_stem Full output path without a file extension.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param notebook_basename PNG basename to copy into `notebook/figures`.
#' @param dpi PNG resolution.
#'
#' @return Named paths for the PNG, PDF, and notebook PNG, invisibly.
#' @export
# ANALYSIS_OK[R026]: exported plot writer is called directly by publication and DE phase scripts.
# ANALYSIS_OK[smuggled-default]: exported plot writer preserves the publication PNG resolution default.
save_publication_plot <- function(
  plot,
  output_stem,
  width,
  height,
  notebook_basename,
  dpi = 300
) {
  if (
    length(output_stem) != 1L ||
      !is.character(output_stem) ||
      !nzchar(output_stem)
  ) {
    stop("output_stem must be one non-empty path.", call. = FALSE)
  }
  if (
    length(notebook_basename) != 1L ||
      !is.character(notebook_basename) ||
      !nzchar(notebook_basename) ||
      identical(notebook_basename, basename(notebook_basename)) == FALSE
  ) {
    stop("notebook_basename must be one non-empty basename.", call. = FALSE)
  }
  if (
    length(width) != 1L ||
      !is.numeric(width) ||
      !is.finite(width) ||
      width <= 0 ||
      length(height) != 1L ||
      !is.numeric(height) ||
      !is.finite(height) ||
      height <= 0
  ) {
    stop("width and height must be positive finite numbers.", call. = FALSE)
  }
  if (length(dpi) != 1L || !is.numeric(dpi) || !is.finite(dpi) || dpi <= 0) {
    stop("dpi must be a positive finite number.", call. = FALSE)
  }

  png_path <- paste0(output_stem, ".png")
  pdf_path <- paste0(output_stem, ".pdf")
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )

  notebook_dir <- here::here("notebook", "figures")
  dir.create(notebook_dir, recursive = TRUE, showWarnings = FALSE)
  notebook_path <- file.path(notebook_dir, notebook_basename)
  .copy_notebook_figure(png_path, notebook_path)

  invisible(c(png = png_path, pdf = pdf_path, notebook = notebook_path))
}

# ANALYSIS_OK[R026]: private mirror helper is called by the exported plot writer in this module.
# Replace a notebook figure only through a verified regular temporary file. POSIX
# rename-over is attempted first; the fallback moves the old destination aside
# and restores it if installation fails.
.copy_notebook_figure <- function(source, destination) {
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
  if (
    !file.exists(destination) ||
      isTRUE(file.info(destination)$isdir) ||
      !isTRUE(file.info(destination)$isdir == FALSE)
  ) {
    stop(
      "Notebook figure destination must be an existing regular file: ",
      destination,
      call. = FALSE
    )
  }
  if (!file.exists(source) || !isTRUE(file.info(source)$isdir == FALSE)) {
    stop(
      "Notebook figure source is not a regular file: ",
      source,
      call. = FALSE
    )
  }

  temporary <- tempfile(
    pattern = paste0(".", basename(destination), "."),
    tmpdir = dirname(destination)
  )
  backup <- tempfile(
    pattern = paste0(".", basename(destination), ".prior."),
    tmpdir = dirname(destination)
  )
  displaced <- tempfile(
    pattern = paste0(".", basename(destination), ".old."),
    tmpdir = dirname(destination)
  )
  on.exit(
    {
      if (file.exists(temporary)) {
        unlink(temporary)
      }
      if (file.exists(backup)) {
        unlink(backup)
      }
      if (file.exists(displaced)) unlink(displaced)
    },
    add = TRUE
  )
  if (
    !file.copy(source, temporary, overwrite = FALSE, copy.date = TRUE) ||
      !file.exists(temporary) ||
      !isTRUE(file.info(temporary)$isdir == FALSE)
  ) {
    stop(
      "Failed to create temporary notebook figure: ",
      destination,
      call. = FALSE
    )
  }
  if (nzchar(Sys.readlink(temporary))) {
    stop("Temporary notebook figure is a symlink: ", temporary, call. = FALSE)
  }
  source_hash <- digest::digest(source, algo = "sha256", file = TRUE)
  temporary_hash <- digest::digest(temporary, algo = "sha256", file = TRUE)
  if (!identical(source_hash, temporary_hash)) {
    stop(
      "Temporary notebook figure hash mismatch: ",
      destination,
      call. = FALSE
    )
  }
  source_dimensions <- magick::image_info(magick::image_read(source))[
    1L,
    c("width", "height")
  ]
  temporary_dimensions <- magick::image_info(magick::image_read(temporary))[
    1L,
    c("width", "height")
  ]
  if (!identical(source_dimensions, temporary_dimensions)) {
    stop(
      "Temporary notebook figure dimensions mismatch: ",
      destination,
      call. = FALSE
    )
  }
  if (!file.copy(destination, backup, overwrite = FALSE, copy.date = TRUE)) {
    stop(
      "Failed to preserve existing notebook figure: ",
      destination,
      call. = FALSE
    )
  }
  installed <- file.rename(temporary, destination)
  if (!installed) {
    if (!file.rename(destination, displaced)) {
      stop(
        "Failed to replace notebook figure without risking its existing file: ",
        destination,
        call. = FALSE
      )
    }
    installed <- file.rename(temporary, destination)
    if (!installed) {
      restored <- file.rename(displaced, destination)
      if (!restored) {
        restored <- file.copy(
          backup,
          destination,
          overwrite = TRUE,
          copy.date = TRUE
        )
      }
      if (!restored) {
        stop(
          "Failed to replace notebook figure and restore its existing file: ",
          destination,
          call. = FALSE
        )
      }
      stop(
        "Failed to replace notebook figure; existing file was preserved: ",
        destination,
        call. = FALSE
      )
    }
  }
  destination_hash <- digest::digest(destination, algo = "sha256", file = TRUE)
  if (!identical(source_hash, destination_hash)) {
    restored <- file.copy(
      backup,
      destination,
      overwrite = TRUE,
      copy.date = TRUE
    )
    if (!restored) {
      failed_destination <- tempfile(
        pattern = paste0(".", basename(destination), ".failed."),
        tmpdir = dirname(destination)
      )
      if (file.rename(destination, failed_destination)) {
        restored <- file.rename(backup, destination)
        if (!restored) {
          restored <- file.copy(
            backup,
            destination,
            overwrite = FALSE,
            copy.date = TRUE
          )
        }
        unlink(failed_destination)
      }
    }
    if (!restored) {
      stop(
        "Notebook figure hash mismatch and existing file could not be restored: ",
        destination,
        call. = FALSE
      )
    }
    stop(
      "Notebook figure hash mismatch; existing file was restored: ",
      destination,
      call. = FALSE
    )
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
#'
#' @return Named PNG and PDF paths, invisibly.
#' @export
# ANALYSIS_OK[R026]: exported heatmap writer is called by the publication-figures phase script.
write_curated_marker_heatmap <- function(
  sobj,
  cluster_column,
  layer,
  output_stem,
  width,
  height,
  cluster_cells = FALSE
) {
  if (!cluster_column %in% colnames(sobj@meta.data)) {
    stop("Missing cluster metadata column: ", cluster_column, call. = FALSE)
  }
  assay <- SeuratObject::DefaultAssay(sobj)
  if (!assay %in% SeuratObject::Assays(sobj)) {
    stop("Missing default assay: ", assay, call. = FALSE)
  }
  if (!layer %in% SeuratObject::Layers(sobj[[assay]])) {
    stop(
      "Missing expression layer '",
      layer,
      "' in assay ",
      assay,
      call. = FALSE
    )
  }
  if (
    !identical(names(cell_type_marker_genes), names(cell_type_marker_labels))
  ) {
    stop(
      "cell_type_marker_genes and cell_type_marker_labels must have identical names.",
      call. = FALSE
    )
  }

  marker_table <- stack(cell_type_marker_genes)
  colnames(marker_table) <- c("gene", "cell_type")
  duplicated_markers <- marker_table$gene[duplicated(marker_table$gene)]
  if (length(duplicated_markers) > 0L) {
    stop(
      "Marker gene(s) assigned to more than one cell type: ",
      paste(unique(duplicated_markers), collapse = ", "),
      call. = FALSE
    )
  }
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
  if (!isTRUE(cluster_cells)) {
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
  cluster_colors <- stats::setNames(
    grDevices::hcl.colors(length(cluster_levels), palette = "Temps"),
    levels(cell_cluster_labels)
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
  column_annotation <- if (isTRUE(cluster_cells)) {
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
    use_raster = TRUE
  )
  if (!isTRUE(cluster_cells)) {
    heatmap_arguments$column_title <- NULL
  }
  heatmap <- do.call(ComplexHeatmap::Heatmap, heatmap_arguments)

  dir.create(
    dirname(paste0(output_stem, ".png")),
    recursive = TRUE,
    showWarnings = FALSE
  )
  png_path <- paste0(output_stem, ".png")
  pdf_path <- paste0(output_stem, ".pdf")
  grDevices::png(
    png_path,
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  if (isTRUE(cluster_cells)) {
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
  if (isTRUE(cluster_cells)) {
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

# ANALYSIS_OK[R026]: private heatmap drawing helper is called by the exported writer in this module.
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
# ANALYSIS_OK[R026]: exported module/p27 heatmap writer is called by the publication-figures phase script.
write_module_p27_heatmap <- function(
  module_scores,
  p27_enrichment,
  output_stem,
  width,
  height
) {
  if (
    !is.matrix(module_scores) ||
      !is.numeric(module_scores) ||
      is.null(rownames(module_scores)) ||
      is.null(colnames(module_scores)) ||
      nrow(module_scores) == 0L ||
      ncol(module_scores) == 0L
  ) {
    stop(
      "module_scores must be a non-empty named numeric matrix.",
      call. = FALSE
    )
  }
  required_cols <- c("cluster", "z_score")
  missing_cols <- setdiff(required_cols, colnames(p27_enrichment))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing p27 enrichment column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
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

  dir.create(
    dirname(paste0(output_stem, ".png")),
    recursive = TRUE,
    showWarnings = FALSE
  )
  png_path <- paste0(output_stem, ".png")
  pdf_path <- paste0(output_stem, ".pdf")
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
