#!/usr/bin/env Rscript

# Find positive markers for the fixed final MG-selected clustering and write
# the publication marker tables and dotplot.

suppressPackageStartupMessages({
  library(here)
})
here::i_am("scripts/03-marker-analysis.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

config <- publication_config()
seed <- config$seed
assay <- "RNA"
expression_layer <- "data"
counts_layer <- "counts"
top_n <- 5L
min_pct <- 0.10
logfc_threshold <- 0.25
min_diff_pct <- 0
min_cells_group <- 3L
return_threshold <- 1
branch_tag <- "pflog_mg_selected_no_filter_cc"
dimensions <- 20L
resolution <- 0.5

input_path <- config$selected$mg$path
cluster_column <- config$selected$mg$column
table_dir <- file.path(config$paths$tables, "mg_selected")
figure_dir <- file.path(config$paths$figures, "mg_selected")

layer_tag <- expression_layer
resolution_tag <- "0.5"
full_marker_path <- file.path(
  table_dir,
  sprintf(
    "find_all_markers_%s_%s_dims%d_res%s.csv",
    layer_tag,
    branch_tag,
    dimensions,
    resolution_tag
  )
)
top_marker_path <- file.path(
  table_dir,
  sprintf(
    "find_all_markers_top%d_%s_%s_dims%d_res%s.csv",
    top_n,
    layer_tag,
    branch_tag,
    dimensions,
    resolution_tag
  )
)
summary_path <- file.path(
  table_dir,
  sprintf(
    "find_all_markers_summary_%s_%s_dims%d_res%s.csv",
    layer_tag,
    branch_tag,
    dimensions,
    resolution_tag
  )
)
identity_map_path <- file.path(
  table_dir,
  sprintf(
    "find_all_markers_identity_map_%s_dims%d_res%s.csv",
    branch_tag,
    dimensions,
    resolution_tag
  )
)
dotplot_tag <- sprintf(
  "mg_selected_cluster_marker_dotplot_%s_%s_dims%d_res%s_top%d",
  layer_tag,
  branch_tag,
  dimensions,
  resolution_tag,
  top_n
)
png_path <- file.path(figure_dir, paste0(dotplot_tag, ".png"))
pdf_path <- file.path(figure_dir, paste0(dotplot_tag, ".pdf"))
notebook_png_path <- file.path(
  config$paths$notebook_figures,
  basename(png_path)
)

output_paths <- c(
  full_marker_path,
  top_marker_path,
  summary_path,
  identity_map_path,
  png_path,
  pdf_path
)
assert_output_available(output_paths, config$overwrite)

# ---- validation and work ----

if (!identical(seed, 1312L)) {
  stop("Unexpected publication seed.", call. = FALSE)
}
if (!identical(input_path, config$selected$mg$path)) {
  stop("Marker input path is not the fixed final MG object.", call. = FALSE)
}
if (!identical(cluster_column, config$selected$mg$column)) {
  stop("Marker cluster column is not the fixed final MG column.", call. = FALSE)
}
if (!file.exists(input_path)) {
  stop("Input Seurat object does not exist: ", input_path, call. = FALSE)
}

sobj <- readRDS(input_path)
assert_frozen_input(input_path, sobj, config$frozen$mg)
if (!inherits(sobj, "Seurat")) {
  stop(
    "Input RDS does not contain a Seurat object: ",
    input_path,
    call. = FALSE
  )
}
object_branch_tag <- sobj@misc$clustering$branch_tag
if (!is.null(object_branch_tag) && !identical(object_branch_tag, branch_tag)) {
  stop(
    "Input object branch tag mismatch: expected ",
    branch_tag,
    ", found ",
    object_branch_tag,
    call. = FALSE
  )
}
if (!cluster_column %in% colnames(sobj[[]])) {
  stop("Missing cluster metadata column: ", cluster_column, call. = FALSE)
}
cluster_values <- trimws(as.character(sobj@meta.data[[cluster_column]]))
if (anyNA(cluster_values) || any(!nzchar(cluster_values))) {
  stop(
    "Cluster metadata column contains missing or empty values: ",
    cluster_column,
    call. = FALSE
  )
}
cluster_levels <- unique(cluster_values)
if (all(grepl("^-?[0-9]+$", cluster_levels))) {
  cluster_levels <- as.character(sort(as.integer(cluster_levels)))
} else {
  cluster_levels <- sort(cluster_levels, method = "radix")
}

if (!assay %in% names(sobj@assays)) {
  stop("Missing fixed assay '", assay, "'.", call. = FALSE)
}
available_layers <- SeuratObject::Layers(sobj[[assay]])
missing_layers <- setdiff(c(expression_layer, counts_layer), available_layers)
if (length(missing_layers) > 0L) {
  stop(
    "Missing fixed assay layer(s): ",
    paste(missing_layers, collapse = ", "),
    call. = FALSE
  )
}

# The selected no-merge identity map is deliberately encoded here: each
# observed Leiden label remains its own marker identity.
identity_map <- data.frame(
  source_cluster = cluster_levels,
  marker_identity = cluster_levels,
  n_cells = as.integer(table(factor(cluster_values, levels = cluster_levels))),
  stringsAsFactors = FALSE
)
identity_levels <- cluster_levels
# ANALYSIS_OK[R002]: fixed minimum identity count preserves the audited marker-analysis contract.
if (length(identity_levels) < 2L) {
  stop("At least two marker identities are required.", call. = FALSE)
}
marker_identities <- cluster_values
identity_map$decision_source <- "confirmed_no_merge"
identity_map$input_path <- input_path
identity_map$cluster_column <- cluster_column
identity_map$assay <- assay
identity_map$expression_layer <- expression_layer
identity_map$counts_layer <- counts_layer

identity_counts <- table(factor(marker_identities, levels = identity_levels))
if (any(identity_counts < min_cells_group)) {
  stop("A fixed marker identity has fewer than three cells.", call. = FALSE)
}
SeuratObject::Idents(sobj) <- factor(
  marker_identities,
  levels = identity_levels
)

markers <- Seurat::FindAllMarkers(
  object = sobj,
  assay = assay,
  slot = expression_layer,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = min_pct,
  logfc.threshold = logfc_threshold,
  min.diff.pct = min_diff_pct,
  min.cells.group = min_cells_group,
  return.thresh = return_threshold,
  verbose = FALSE
)
if (nrow(markers) == 0L) {
  stop("FindAllMarkers returned zero marker rows.", call. = FALSE)
}
if (!"gene" %in% colnames(markers)) {
  markers$gene <- rownames(markers)
}
markers$gene <- as.character(markers$gene)
if (anyNA(markers$gene) || any(!nzchar(markers$gene))) {
  stop("FindAllMarkers output contains missing gene names.", call. = FALSE)
}
required_marker_columns <- c(
  "p_val",
  "pct.1",
  "pct.2",
  "p_val_adj",
  "cluster"
)
missing_marker_columns <- setdiff(required_marker_columns, colnames(markers))
if (length(missing_marker_columns) > 0L) {
  stop(
    "FindAllMarkers output is missing required column(s): ",
    paste(missing_marker_columns, collapse = ", "),
    call. = FALSE
  )
}
fold_change_column <- if ("avg_log2FC" %in% colnames(markers)) {
  "avg_log2FC"
} else if ("avg_logFC" %in% colnames(markers)) {
  "avg_logFC"
} else {
  stop(
    "FindAllMarkers output must contain avg_log2FC or avg_logFC.",
    call. = FALSE
  )
}
markers$cluster <- as.character(markers$cluster)
unknown_clusters <- setdiff(unique(markers$cluster), identity_levels)
if (length(unknown_clusters) > 0L) {
  stop(
    "FindAllMarkers returned unexpected marker identity value(s): ",
    paste(unknown_clusters, collapse = ", "),
    call. = FALSE
  )
}

markers$pct_diff <- markers$pct.1 - markers$pct.2
# Preserve the existing positive-marker guard: anti-markers are excluded.
ranked_markers <- markers[
  !is.na(markers$pct_diff) & markers$pct_diff > min_diff_pct,
  ,
  drop = FALSE
]
marker_groups <- split(
  seq_len(nrow(ranked_markers)),
  factor(ranked_markers$cluster, levels = identity_levels),
  drop = FALSE
)
ranked_groups <- lapply(marker_groups, function(indices) {
  if (length(indices) == 0L) {
    group_markers <- ranked_markers[integer(0), , drop = FALSE]
    group_markers$rank_within_cluster <- integer(0)
    return(group_markers)
  }
  group_markers <- ranked_markers[indices, , drop = FALSE]
  group_order <- order(
    is.na(group_markers$p_val_adj),
    group_markers$p_val_adj,
    is.na(group_markers$p_val),
    group_markers$p_val,
    is.na(group_markers[[fold_change_column]]),
    -group_markers[[fold_change_column]],
    is.na(group_markers$pct_diff),
    -group_markers$pct_diff,
    is.na(group_markers$pct.1),
    -group_markers$pct.1,
    group_markers$gene,
    method = "radix"
  )
  # ANALYSIS_OK[R005]: reorder each cluster's markers for deterministic ranking without dropping markers.
  group_markers <- group_markers[group_order, , drop = FALSE]
  group_markers$rank_within_cluster <- seq_len(nrow(group_markers))
  group_markers
})
ranked_markers <- do.call(rbind, ranked_groups)
rownames(ranked_markers) <- NULL
preferred_columns <- c(
  "gene",
  "cluster",
  "rank_within_cluster",
  "p_val",
  fold_change_column,
  "pct.1",
  "pct.2",
  "pct_diff",
  "p_val_adj"
)
# ANALYSIS_OK[R005]: reorder output columns to the fixed marker-table schema without dropping rows.
ranked_markers <- ranked_markers[,
  c(preferred_columns, setdiff(colnames(ranked_markers), preferred_columns)),
  drop = FALSE
]

ranked_marker_groups <- split(
  seq_len(nrow(ranked_markers)),
  factor(ranked_markers$cluster, levels = identity_levels),
  drop = FALSE
)
top_groups <- lapply(ranked_marker_groups, function(indices) {
  ranked_markers[utils::head(indices, top_n), , drop = FALSE]
})
top_markers <- do.call(rbind, top_groups)
rownames(top_markers) <- NULL
if (any(!is.na(top_markers$pct_diff) & top_markers$pct_diff <= min_diff_pct)) {
  stop("Top marker selection retained an anti-marker.", call. = FALSE)
}

n_cells <- stats::aggregate(
  n_cells ~ marker_identity,
  data = identity_map,
  FUN = sum
)
n_cells <- stats::setNames(n_cells$n_cells, n_cells$marker_identity)
retained_counts <- table(factor(
  ranked_markers$cluster,
  levels = identity_levels
))
top_counts <- table(factor(top_markers$cluster, levels = identity_levels))
marker_summary <- data.frame(
  marker_identity = identity_levels,
  n_cells = as.integer(n_cells[identity_levels]),
  n_retained_markers = as.integer(retained_counts[identity_levels]),
  n_top_markers = as.integer(top_counts[identity_levels]),
  decision_source = "confirmed_no_merge",
  stringsAsFactors = FALSE
)

# ---- dotplot ----

gene_order <- unique(top_markers$gene)
expression_matrix <- SeuratObject::LayerData(
  sobj[[assay]],
  layer = expression_layer
)
counts_matrix <- SeuratObject::LayerData(sobj[[assay]], layer = counts_layer)
if (length(setdiff(gene_order, rownames(expression_matrix))) > 0L) {
  stop("Top marker gene is missing from the RNA data layer.", call. = FALSE)
}
if (length(setdiff(gene_order, rownames(counts_matrix))) > 0L) {
  stop("Top marker gene is missing from the RNA counts layer.", call. = FALSE)
}
metadata_cells <- rownames(sobj@meta.data)
if (!setequal(colnames(expression_matrix), metadata_cells)) {
  stop(
    "RNA data-layer columns do not match Seurat metadata rows.",
    call. = FALSE
  )
}
if (!setequal(colnames(counts_matrix), metadata_cells)) {
  stop(
    "RNA counts-layer columns do not match Seurat metadata rows.",
    call. = FALSE
  )
}
# ANALYSIS_OK[R005]: reorder expression rows and cells to the validated marker plotting order.
expression_matrix <- expression_matrix[gene_order, metadata_cells, drop = FALSE]
# ANALYSIS_OK[R005]: reorder count rows and cells to the validated marker plotting order.
counts_matrix <- counts_matrix[gene_order, metadata_cells, drop = FALSE]
marker_identities_factor <- factor(marker_identities, levels = identity_levels)
names(marker_identities_factor) <- metadata_cells
# ANALYSIS_OK[R019]: construct dotplot data from the validated marker genes and identities for visualization.
plot_data <- do.call(
  rbind,
  lapply(gene_order, function(gene) {
    do.call(
      rbind,
      lapply(identity_levels, function(identity) {
        cells <- names(marker_identities_factor)[
          marker_identities_factor == identity
        ]
        data.frame(
          gene = gene,
          marker_identity = identity,
          mean_expression = mean(as.numeric(expression_matrix[gene, cells])),
          pct_detected = mean(as.numeric(counts_matrix[gene, cells]) > 0) * 100,
          stringsAsFactors = FALSE
        )
      })
    )
  })
)
scaled_expression <- numeric(nrow(plot_data))
for (indices in split(seq_len(nrow(plot_data)), plot_data$gene)) {
  values <- plot_data$mean_expression[indices]
  value_sd <- stats::sd(values)
  # ANALYSIS_OK[R002]: zero-variance handling preserves the audited marker z-score definition.
  if (length(values) < 2L || is.na(value_sd) || value_sd == 0) {
    scaled_expression[indices] <- 0
  } else {
    scaled_expression[indices] <- (values - mean(values)) / value_sd
  }
}
scaled_expression[is.na(scaled_expression)] <- 0
plot_data$scaled_mean_expression <- scaled_expression
plot_data$gene <- factor(plot_data$gene, levels = rev(gene_order))
plot_data$marker_identity <- factor(
  plot_data$marker_identity,
  levels = identity_levels
)

palette_dotplot <- stats::setNames(
  config$palettes$dotplot,
  c("negative", "positive")
)
blue_ramp <- stats::setNames(
  grDevices::colorRampPalette(c(palette_dotplot[["negative"]], "white"))(4L),
  c("low", "mid_low", "mid_high", "high")
)
pink_ramp <- stats::setNames(
  grDevices::colorRampPalette(c("white", palette_dotplot[["positive"]]))(4L),
  c("low", "mid_low", "mid_high", "high")
)
dotplot_colour_breaks <- c(-3, -2, -1, 0, 1, 2, 3)
dotplot_colour_limits <- c(lower = -3, upper = 3)
# ANALYSIS_OK[R026]: local dotplot-label helper is called by the marker plotting code below.
dotplot_colour_labels <- function(breaks) {
  labels <- as.character(breaks)
  # ANALYSIS_OK[R002]: fixed endpoint tolerance preserves the audited dotplot label mapping.
  labels[abs(breaks - dotplot_colour_limits[["lower"]]) < 1e-8] <- "<= -2"
  # ANALYSIS_OK[R002]: fixed endpoint tolerance preserves the audited dotplot label mapping.
  labels[abs(breaks - dotplot_colour_limits[["upper"]]) < 1e-8] <- ">= 2"
  labels
}
plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = marker_identity,
    y = gene,
    size = pct_detected,
    color = scaled_mean_expression
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_size(
    range = c(0.5, 6),
    limits = c(0, 100),
    name = "Detected cells (%)"
  ) +
  ggplot2::scale_color_stepsn(
    colours = c(
      palette_dotplot[["negative"]],
      blue_ramp[["mid_low"]],
      blue_ramp[["mid_high"]],
      pink_ramp[["mid_low"]],
      pink_ramp[["mid_high"]],
      palette_dotplot[["positive"]]
    ),
    breaks = dotplot_colour_breaks,
    labels = dotplot_colour_labels,
    limits = dotplot_colour_limits,
    oob = function(x, range, ...) {
      range_lower <- min(range)
      range_upper <- max(range)
      pmin(pmax(x, range_lower), range_upper)
    },
    guide = ggplot2::guide_coloursteps(show.limits = TRUE),
    name = sprintf("Mean %s expression\n(row z-score bin)", expression_layer)
  ) +
  ggplot2::labs(
    title = "MG-selected cluster markers",
    x = "Marker identity",
    y = "Top marker gene"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = ggplot2::element_line(linewidth = 0.2),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(hjust = 0.5)
  )
width_in <- max(7.0, 2.5 + 0.45 * length(identity_levels))
height_in <- max(5.0, 2.0 + 0.18 * length(gene_order))

# ---- output ----

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  ranked_markers,
  file = full_marker_path,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  top_markers,
  file = top_marker_path,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  marker_summary,
  file = summary_path,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  identity_map,
  file = identity_map_path,
  row.names = FALSE,
  na = ""
)
ggplot2::ggsave(
  png_path,
  plot,
  width = width_in,
  height = height_in,
  dpi = 300,
  bg = "white"
)
ggplot2::ggsave(
  pdf_path,
  plot,
  width = width_in,
  height = height_in,
  bg = "white"
)

# Replace only an existing regular notebook destination. Never write through a
# symlink: copy and hash-verify a regular temporary sibling before renaming.
if (nzchar(Sys.readlink(notebook_png_path))) {
  stop(
    "Notebook figure destination is a symlink: ",
    notebook_png_path,
    call. = FALSE
  )
}
if (
  !file.exists(notebook_png_path) || isTRUE(file.info(notebook_png_path)$isdir)
) {
  stop(
    "Notebook figure destination must be an existing regular file: ",
    notebook_png_path,
    call. = FALSE
  )
}
temporary_notebook_png <- tempfile(
  pattern = paste0(".", basename(notebook_png_path), "."),
  tmpdir = dirname(notebook_png_path),
  fileext = ".tmp"
)
tryCatch(
  {
    if (
      !isTRUE(file.copy(png_path, temporary_notebook_png, overwrite = FALSE))
    ) {
      stop("Failed to copy notebook figure to temporary file.", call. = FALSE)
    }
    if (
      nzchar(Sys.readlink(temporary_notebook_png)) ||
        !file.exists(temporary_notebook_png) ||
        isTRUE(file.info(temporary_notebook_png)$isdir) ||
        !identical(
          digest::digest(png_path, algo = "sha256", file = TRUE),
          digest::digest(temporary_notebook_png, algo = "sha256", file = TRUE)
        )
    ) {
      stop(
        "Temporary notebook figure failed regular-file/hash verification.",
        call. = FALSE
      )
    }
    source_dimensions <- magick::image_info(
      magick::image_read(png_path)
    )[1L, c("width", "height")]
    temporary_dimensions <- magick::image_info(
      magick::image_read(temporary_notebook_png)
    )[1L, c("width", "height")]
    if (!identical(source_dimensions, temporary_dimensions)) {
      stop(
        "Temporary notebook figure dimensions do not match source.",
        call. = FALSE
      )
    }
    unlink(notebook_png_path)
    if (!isTRUE(file.rename(temporary_notebook_png, notebook_png_path))) {
      stop(
        "Failed to replace notebook figure: ",
        notebook_png_path,
        call. = FALSE
      )
    }
    if (
      !identical(
        digest::digest(png_path, algo = "sha256", file = TRUE),
        digest::digest(notebook_png_path, algo = "sha256", file = TRUE)
      )
    ) {
      stop(
        "Replaced notebook figure hash does not match source.",
        call. = FALSE
      )
    }
  },
  finally = {
    if (file.exists(temporary_notebook_png)) {
      unlink(temporary_notebook_png)
    }
  }
)

message("Wrote full marker CSV: ", full_marker_path)
message("Wrote top marker CSV: ", top_marker_path)
message("Wrote marker summary CSV: ", summary_path)
message("Wrote marker identity map CSV: ", identity_map_path)
message("Wrote marker dotplot PNG: ", png_path)
message("Wrote marker dotplot PDF: ", pdf_path)
message("Replaced notebook figure: ", notebook_png_path)
