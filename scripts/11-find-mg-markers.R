#!/usr/bin/env Rscript

# Run descriptive Seurat FindAllMarkers ranking for the chosen MG-selected
# clustering.
#
# Usage:
#   Rscript scripts/11-find-mg-markers.R \
#     [--input <clustered-seurat-object.rds>] \
#     [--branch-tag <branch tag>] \
#     [--elbow-n <positive integer>] \
#     [--dims <positive integer>] \
#     [--resolution <resolution string>] \
#     [--assay <assay>] \
#     [--layer <assay layer>] \
#     [--counts-layer <raw counts layer>] \
#     [--top-n <positive integer>] \
#     [--min-pct <number>] \
#     [--logfc-threshold <number>] \
#     [--min-diff-pct <number>] \
#     [--min-cells-group <positive integer>] \
#     (--confirm-no-merge | --cluster-map <csv>) \
#     [--table-dir <output directory>] \
#     [--figure-dir <output directory>] \
#     [--overwrite]
#
# Defaults target the MG-selected data-layer marker run on the PFlog
# no-filter-cc branch:
#   CURRENT_OBJECT_DIR/cluster_pflog_mg_selected_no_filter_cc_elbow20.rds
#   cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5
# Marker testing defaults to the Seurat `data` layer even though the chosen
# clustering is PFlog-derived. Seurat's default FoldChange math is not valid for
# the PFlog layer; this script enforces pct.1 > pct.2 for retained positive
# markers to prevent anti-markers from being reported as cluster markers.
#
# The cluster-map CSV must contain source_cluster and marker_identity columns.
# Duplicate marker_identity values are the explicit cluster-merge mechanism.
# Use --confirm-no-merge to record that the Leiden labels should be used as-is.
#
# Outputs:
#   TABLE_DIR/mg_selected/find_all_markers_<layer>_<branch>_dims<dims>_res<resolution>.csv
#   TABLE_DIR/mg_selected/find_all_markers_top<N>_<layer>_<branch>_dims<dims>_res<resolution>.csv
#   TABLE_DIR/mg_selected/find_all_markers_summary_<layer>_<branch>_dims<dims>_res<resolution>.csv
#   TABLE_DIR/mg_selected/find_all_markers_identity_map_<branch>_dims<dims>_res<resolution>.csv
#   FIGURE_DIR/mg_selected/mg_selected_cluster_marker_dotplot_<layer>_<branch>_dims<dims>_res<resolution>_top<N>.(png|pdf)
#   notebook/figures/<dotplot png filename> symlink.

suppressPackageStartupMessages({
  library(here)
})
here::i_am("scripts/11-find-mg-markers.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

palette_dotplot_pair <- get(
  "palette_dotplot_pair",
  envir = asNamespace("ESPI"),
  inherits = FALSE
)
palette_dotplot_pair <- stats::setNames(
  palette_dotplot_pair,
  c("negative", "positive")
)

# ---- parameters ----

# ANALYSIS_OK[R026]: standalone CLI entrypoint helper is intentionally local to this script.
get_mg_marker_arg <- function(args, flag, default = NULL) {
  match_index <- match(flag, args)
  if (is.na(match_index)) {
    return(default)
  }
  if (
    match_index == length(args) || startsWith(args[[match_index + 1L]], "--")
  ) {
    stop("Missing value for ", flag, call. = FALSE)
  }
  args[[match_index + 1L]]
}


cli_args <- commandArgs(trailingOnly = TRUE)
value_flags <- c(
  "--input",
  "--branch-tag",
  "--elbow-n",
  "--dims",
  "--resolution",
  "--assay",
  "--layer",
  "--counts-layer",
  "--top-n",
  "--min-pct",
  "--logfc-threshold",
  "--min-diff-pct",
  "--min-cells-group",
  "--cluster-map",
  "--table-dir",
  "--figure-dir"
)
boolean_flags <- c("--confirm-no-merge", "--overwrite")
allowed_flags <- c(value_flags, boolean_flags)
seen_flags <- cli_args[startsWith(cli_args, "--")]
unknown_flags <- seen_flags[!seen_flags %in% allowed_flags]
if (length(unknown_flags) > 0L) {
  stop(
    "Unknown argument(s): ",
    paste(unknown_flags, collapse = ", "),
    call. = FALSE
  )
}
repeated_flags <- unique(seen_flags[duplicated(seen_flags)])
if (length(repeated_flags) > 0L) {
  stop(
    "Repeated argument(s): ",
    paste(repeated_flags, collapse = ", "),
    call. = FALSE
  )
}
flag_positions <- which(startsWith(cli_args, "--"))
value_positions <- flag_positions[cli_args[flag_positions] %in% value_flags] +
  1L
unexpected_values <- cli_args[
  !startsWith(cli_args, "--") & !(seq_along(cli_args) %in% value_positions)
]
if (length(unexpected_values) > 0L) {
  stop(
    "Unexpected positional argument(s): ",
    paste(unexpected_values, collapse = ", "),
    call. = FALSE
  )
}

required_packages <- c("Seurat", "SeuratObject", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

branch_tag <- get_mg_marker_arg(
  cli_args,
  "--branch-tag",
  "pflog_mg_selected_no_filter_cc"
)
elbow_n <- as.integer(get_mg_marker_arg(cli_args, "--elbow-n", "20"))
dims <- as.integer(get_mg_marker_arg(cli_args, "--dims", "20"))
resolution <- get_mg_marker_arg(cli_args, "--resolution", "0.5")
assay_arg <- get_mg_marker_arg(cli_args, "--assay", NULL)
expression_layer <- get_mg_marker_arg(cli_args, "--layer", "data")
counts_layer <- get_mg_marker_arg(cli_args, "--counts-layer", "counts")
top_n <- as.integer(get_mg_marker_arg(cli_args, "--top-n", "5"))
min_pct <- as.numeric(get_mg_marker_arg(cli_args, "--min-pct", "0.10"))
logfc_threshold <- as.numeric(get_mg_marker_arg(
  cli_args,
  "--logfc-threshold",
  "0.25"
))
min_diff_pct <- as.numeric(get_mg_marker_arg(cli_args, "--min-diff-pct", "0"))
min_cells_group <- as.integer(get_mg_marker_arg(
  cli_args,
  "--min-cells-group",
  "3"
))
cluster_map_path <- get_mg_marker_arg(cli_args, "--cluster-map", NULL)
confirm_no_merge <- "--confirm-no-merge" %in% cli_args
overwrite_outputs <- "--overwrite" %in% cli_args
input_path <- get_mg_marker_arg(
  cli_args,
  "--input",
  file.path(
    CURRENT_OBJECT_DIR,
    sprintf("cluster_%s_elbow%d.rds", branch_tag, elbow_n)
  )
)
table_dir <- get_mg_marker_arg(
  cli_args,
  "--table-dir",
  file.path(TABLE_DIR, "mg_selected")
)
figure_dir <- get_mg_marker_arg(
  cli_args,
  "--figure-dir",
  file.path(FIGURE_DIR, "mg_selected")
)

# ---- helpers ----
MIN_MARKER_IDENTITIES <- 2L
ZERO_STANDARD_DEVIATION <- 0
DOTPLOT_ZSCORE_EPSILON <- 1e-8
DOTPLOT_ZSCORE_LIMITS <- c(lower = -3, upper = 3)

# ANALYSIS_OK[R026]: standalone validation helper is intentionally local to this script.
assert_positive_integer <- function(x, name) {
  if (length(x) != 1L || is.na(x) || x <= 0L) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
}

# ANALYSIS_OK[R026]: standalone validation helper is intentionally local to this script.
assert_fraction <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0 || x > 1) {
    stop(name, " must be a finite number from 0 to 1.", call. = FALSE)
  }
}

# ANALYSIS_OK[R026]: standalone validation helper is intentionally local to this script.
assert_nonnegative_number <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    stop(name, " must be a finite non-negative number.", call. = FALSE)
  }
}


# ANALYSIS_OK[R026]: standalone validation helper is intentionally local to this script.
assert_absent_or_overwrite <- function(paths, overwrite) {
  existing_paths <- paths[
    file.exists(paths) | nzchar(Sys.readlink(paths))
  ]
  if (length(existing_paths) > 0L && !isTRUE(overwrite)) {
    stop(
      "Output file(s) already exist; rerun with --overwrite to replace: ",
      paste(existing_paths, collapse = ", "),
      call. = FALSE
    )
  }
}


# ANALYSIS_OK[R026]: standalone mapping helper is intentionally local to this script.
read_identity_map <- function(path, cluster_levels, cluster_values) {
  if (!file.exists(path)) {
    stop("Cluster map CSV does not exist: ", path, call. = FALSE)
  }
  identity_map <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_columns <- c("source_cluster", "marker_identity")
  missing_columns <- setdiff(required_columns, colnames(identity_map))
  if (length(missing_columns) > 0L) {
    stop(
      "Cluster map CSV is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  # ANALYSIS_OK[filter-ledger]: retain only the two validated cluster-map columns.
  identity_map <- identity_map[, required_columns, drop = FALSE]
  identity_map$source_cluster <- trimws(as.character(
    identity_map$source_cluster
  ))
  identity_map$marker_identity <- trimws(as.character(
    identity_map$marker_identity
  ))
  if (
    anyNA(identity_map$source_cluster) ||
      any(!nzchar(identity_map$source_cluster))
  ) {
    stop("Cluster map source_cluster values must be non-empty.", call. = FALSE)
  }
  if (
    anyNA(identity_map$marker_identity) ||
      any(!nzchar(identity_map$marker_identity))
  ) {
    stop("Cluster map marker_identity values must be non-empty.", call. = FALSE)
  }
  duplicated_sources <- unique(
    identity_map$source_cluster[duplicated(identity_map$source_cluster)]
  )
  if (length(duplicated_sources) > 0L) {
    stop(
      "Cluster map source_cluster values must be unique; duplicated: ",
      paste(duplicated_sources, collapse = ", "),
      call. = FALSE
    )
  }
  unknown_sources <- setdiff(identity_map$source_cluster, cluster_levels)
  if (length(unknown_sources) > 0L) {
    stop(
      "Cluster map contains unknown source_cluster value(s): ",
      paste(unknown_sources, collapse = ", "),
      call. = FALSE
    )
  }
  missing_sources <- setdiff(cluster_levels, identity_map$source_cluster)
  if (length(missing_sources) > 0L) {
    stop(
      "Cluster map is missing observed source_cluster value(s): ",
      paste(missing_sources, collapse = ", "),
      call. = FALSE
    )
  }
  # ANALYSIS_OK[filter-ledger]: reorder the validated map to observed cluster order.
  identity_map <- identity_map[
    match(cluster_levels, identity_map$source_cluster),
  ]
  identity_map$n_cells <- as.integer(
    table(factor(cluster_values, levels = cluster_levels))
  )
  rownames(identity_map) <- NULL
  identity_map
}

# ANALYSIS_OK[R026]: standalone identity-map helper is intentionally local to this script.
confirmed_identity_map <- function(cluster_levels, cluster_values) {
  data.frame(
    source_cluster = cluster_levels,
    marker_identity = cluster_levels,
    n_cells = as.integer(table(factor(
      cluster_values,
      levels = cluster_levels
    ))),
    stringsAsFactors = FALSE
  )
}

# ANALYSIS_OK[R026]: standalone identity-map helper is intentionally local to this script.
marker_identity_values <- function(cluster_values, identity_map) {
  marker_by_cluster <- stats::setNames(
    identity_map$marker_identity,
    identity_map$source_cluster
  )
  unname(marker_by_cluster[cluster_values])
}

# ANALYSIS_OK[R026]: standalone identity validation helper is intentionally local to this script.
validate_marker_identities <- function(
  marker_identities,
  identity_levels,
  min_cells
) {
  if (anyNA(marker_identities) || any(!nzchar(marker_identities))) {
    stop("Marker identity mapping produced missing labels.", call. = FALSE)
  }
  if (length(identity_levels) < MIN_MARKER_IDENTITIES) {
    stop("At least two marker identities are required.", call. = FALSE)
  }
  identity_counts <- table(factor(marker_identities, levels = identity_levels))
  too_small <- names(identity_counts)[identity_counts < min_cells]
  if (length(too_small) > 0L) {
    stop(
      "Marker identity group(s) have fewer than --min-cells-group cells: ",
      paste(
        sprintf("%s=%d", too_small, identity_counts[too_small]),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
}

# ANALYSIS_OK[R026]: standalone marker helper is intentionally local to this script.
find_fold_change_column <- function(markers) {
  if ("avg_log2FC" %in% colnames(markers)) {
    return("avg_log2FC")
  }
  if ("avg_logFC" %in% colnames(markers)) {
    return("avg_logFC")
  }
  stop(
    "FindAllMarkers output must contain avg_log2FC or avg_logFC.",
    call. = FALSE
  )
}

# ANALYSIS_OK[R026]: standalone marker helper is intentionally local to this script.
prepare_marker_table <- function(markers, identity_levels, min_diff_pct) {
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
  required_columns <- c("p_val", "pct.1", "pct.2", "p_val_adj", "cluster")
  missing_columns <- setdiff(required_columns, colnames(markers))
  if (length(missing_columns) > 0L) {
    stop(
      "FindAllMarkers output is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  fold_change_column <- find_fold_change_column(markers)
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
  # ANALYSIS_OK[filter-ledger]: retain markers meeting the predeclared pct-difference rule.
  markers <- markers[
    !is.na(markers$pct_diff) & markers$pct_diff > min_diff_pct,
    ,
    drop = FALSE
  ]
  marker_groups <- split(
    seq_len(nrow(markers)),
    factor(markers$cluster, levels = identity_levels),
    drop = FALSE
  )
  ranked_groups <- lapply(marker_groups, function(indices) {
    if (length(indices) == 0L) {
      group_markers <- markers[integer(0), , drop = FALSE]
      group_markers$rank_within_cluster <- integer(0)
      return(group_markers)
    }
    group_markers <- markers[indices, , drop = FALSE]
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
    # ANALYSIS_OK[filter-ledger]: sort each cluster's retained markers by deterministic rank.
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
  ranked_markers[,
    c(
      preferred_columns,
      setdiff(colnames(ranked_markers), preferred_columns)
    ),
    drop = FALSE
  ]
}

# ANALYSIS_OK[R026]: standalone marker helper is intentionally local to this script.
select_top_markers <- function(markers, identity_levels, top_n, min_diff_pct) {
  marker_groups <- split(
    seq_len(nrow(markers)),
    factor(markers$cluster, levels = identity_levels),
    drop = FALSE
  )
  top_groups <- lapply(marker_groups, function(indices) {
    markers[utils::head(indices, top_n), , drop = FALSE]
  })
  top_markers <- do.call(rbind, top_groups)
  rownames(top_markers) <- NULL
  invalid_top_markers <- !is.na(top_markers$pct_diff) &
    top_markers$pct_diff <= min_diff_pct
  if (any(invalid_top_markers)) {
    stop(
      "Top marker selection retained row(s) with pct_diff <= --min-diff-pct.",
      call. = FALSE
    )
  }
  top_markers
}

# ANALYSIS_OK[R026]: standalone marker helper is intentionally local to this script.
summarise_markers <- function(
  identity_map,
  identity_levels,
  ranked_markers,
  top_markers,
  decision_source
) {
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
  data.frame(
    marker_identity = identity_levels,
    n_cells = as.integer(n_cells[identity_levels]),
    n_retained_markers = as.integer(retained_counts[identity_levels]),
    n_top_markers = as.integer(top_counts[identity_levels]),
    decision_source = decision_source,
    stringsAsFactors = FALSE
  )
}

# ANALYSIS_OK[R026]: standalone plotting helper is intentionally local to this script.
build_dotplot_data <- function(
  sobj,
  assay,
  expression_layer,
  counts_layer,
  genes,
  marker_identities,
  identity_levels
) {
  expression_matrix <- SeuratObject::LayerData(
    sobj[[assay]],
    layer = expression_layer
  )
  counts_matrix <- SeuratObject::LayerData(sobj[[assay]], layer = counts_layer)
  missing_expression_genes <- setdiff(genes, rownames(expression_matrix))
  if (length(missing_expression_genes) > 0L) {
    stop(
      "Top marker gene(s) missing from expression layer: ",
      paste(missing_expression_genes, collapse = ", "),
      call. = FALSE
    )
  }
  missing_count_genes <- setdiff(genes, rownames(counts_matrix))
  if (length(missing_count_genes) > 0L) {
    stop(
      "Top marker gene(s) missing from counts layer: ",
      paste(missing_count_genes, collapse = ", "),
      call. = FALSE
    )
  }
  metadata_cells <- rownames(sobj@meta.data)
  if (!setequal(colnames(expression_matrix), metadata_cells)) {
    stop(
      "Expression layer columns do not match Seurat metadata rows.",
      call. = FALSE
    )
  }
  if (!setequal(colnames(counts_matrix), metadata_cells)) {
    stop(
      "Counts layer columns do not match Seurat metadata rows.",
      call. = FALSE
    )
  }
  # ANALYSIS_OK[alignment-check]: align expression and counts layers to validated metadata cells.
  expression_matrix <- expression_matrix[genes, metadata_cells, drop = FALSE]
  # ANALYSIS_OK[alignment-check]: align expression and counts layers to validated metadata cells.
  counts_matrix <- counts_matrix[genes, metadata_cells, drop = FALSE]
  marker_identities <- factor(marker_identities, levels = identity_levels)
  names(marker_identities) <- metadata_cells
  # ANALYSIS_OK[plot-filter]: build plot-only summaries from validated marker identities.
  plot_data <- do.call(
    rbind,
    lapply(genes, function(gene) {
      do.call(
        rbind,
        lapply(identity_levels, function(identity) {
          cells <- names(marker_identities)[marker_identities == identity]
          data.frame(
            gene = gene,
            marker_identity = identity,
            mean_expression = mean(as.numeric(expression_matrix[gene, cells])),
            pct_detected = mean(as.numeric(counts_matrix[gene, cells]) > 0) *
              100,
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )
  scaled_expression <- numeric(nrow(plot_data))
  gene_indices <- split(seq_len(nrow(plot_data)), plot_data$gene)
  for (indices in gene_indices) {
    values <- plot_data$mean_expression[indices]
    value_sd <- stats::sd(values)
    if (
      length(values) < MIN_MARKER_IDENTITIES ||
        is.na(value_sd) ||
        value_sd == ZERO_STANDARD_DEVIATION
    ) {
      scaled_expression[indices] <- 0
    } else {
      scaled_expression[indices] <- (values - mean(values)) / value_sd
    }
  }
  scaled_expression[is.na(scaled_expression)] <- 0
  plot_data$scaled_mean_expression <- scaled_expression
  plot_data$gene <- factor(plot_data$gene, levels = rev(genes))
  plot_data$marker_identity <- factor(
    plot_data$marker_identity,
    levels = identity_levels
  )
  plot_data
}

# ANALYSIS_OK[R026]: standalone plotting helper is intentionally local to this script.
marker_dotplot <- function(plot_data, expression_layer) {
  blue_ramp <- stats::setNames(
    grDevices::colorRampPalette(c(
      palette_dotplot_pair[["negative"]],
      "white"
    ))(4L),
    c("low", "mid_low", "mid_high", "high")
  )
  pink_ramp <- stats::setNames(
    grDevices::colorRampPalette(c(
      "white",
      palette_dotplot_pair[["positive"]]
    ))(4L),
    c("low", "mid_low", "mid_high", "high")
  )
  dotplot_colour_breaks <- c(-3, -2, -1, 0, 1, 2, 3)
  dotplot_colour_limits <- DOTPLOT_ZSCORE_LIMITS
  # ANALYSIS_OK[R026]: nested plotting label helper is intentionally local to this script.
  dotplot_colour_labels <- function(breaks) {
    labels <- as.character(breaks)
    labels[
      abs(breaks - dotplot_colour_limits[["lower"]]) < DOTPLOT_ZSCORE_EPSILON
    ] <- "<= -2"
    labels[
      abs(breaks - dotplot_colour_limits[["upper"]]) < DOTPLOT_ZSCORE_EPSILON
    ] <- ">= 2"
    labels
  }
  ggplot2::ggplot(
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
        palette_dotplot_pair[["negative"]],
        blue_ramp[["mid_low"]],
        blue_ramp[["mid_high"]],
        pink_ramp[["mid_low"]],
        pink_ramp[["mid_high"]],
        palette_dotplot_pair[["positive"]]
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
}

# ---- validation ----

assert_scalar_character(branch_tag, "--branch-tag")
if (!grepl("^[A-Za-z0-9_]+$", branch_tag)) {
  stop(
    "--branch-tag must contain only letters, numbers, and underscores.",
    call. = FALSE
  )
}
assert_positive_integer(elbow_n, "--elbow-n")
assert_positive_integer(dims, "--dims")
assert_scalar_character(resolution, "--resolution")
# ANALYSIS_OK[warning-suppression]: invalid numeric CLI text is rejected by explicit validation below.
resolution_number <- suppressWarnings(as.numeric(resolution))
if (
  is.na(resolution_number) ||
    !is.finite(resolution_number) ||
    resolution_number <= 0
) {
  stop("--resolution must be a positive number string.", call. = FALSE)
}
if (!is.null(assay_arg)) {
  assert_scalar_character(assay_arg, "--assay")
}
assert_scalar_character(expression_layer, "--layer")
assert_scalar_character(counts_layer, "--counts-layer")
assert_positive_integer(top_n, "--top-n")
assert_fraction(min_pct, "--min-pct")
assert_nonnegative_number(logfc_threshold, "--logfc-threshold")
assert_fraction(min_diff_pct, "--min-diff-pct")
assert_positive_integer(min_cells_group, "--min-cells-group")
assert_scalar_character(input_path, "--input")
assert_scalar_character(table_dir, "--table-dir")
assert_scalar_character(figure_dir, "--figure-dir")
if (!xor(confirm_no_merge, !is.null(cluster_map_path))) {
  stop(
    "Provide exactly one of --confirm-no-merge or --cluster-map before ",
    "running FindAllMarkers.",
    call. = FALSE
  )
}
if (!file.exists(input_path)) {
  stop("Input Seurat object does not exist: ", input_path, call. = FALSE)
}
if (!is.null(cluster_map_path)) {
  assert_scalar_character(cluster_map_path, "--cluster-map")
}

resolution_tag <- gsub("[^A-Za-z0-9_.-]", "_", resolution)
layer_tag <- gsub("[^A-Za-z0-9_.-]", "_", expression_layer)
full_marker_path <- file.path(
  table_dir,
  sprintf(
    "find_all_markers_%s_%s_dims%d_res%s.csv",
    layer_tag,
    branch_tag,
    dims,
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
    dims,
    resolution_tag
  )
)
summary_path <- file.path(
  table_dir,
  sprintf(
    "find_all_markers_summary_%s_%s_dims%d_res%s.csv",
    layer_tag,
    branch_tag,
    dims,
    resolution_tag
  )
)
identity_map_path <- file.path(
  table_dir,
  sprintf(
    "find_all_markers_identity_map_%s_dims%d_res%s.csv",
    branch_tag,
    dims,
    resolution_tag
  )
)
dotplot_tag <- sprintf(
  "mg_selected_cluster_marker_dotplot_%s_%s_dims%d_res%s_top%d",
  layer_tag,
  branch_tag,
  dims,
  resolution_tag,
  top_n
)
png_path <- file.path(figure_dir, sprintf("%s.png", dotplot_tag))
pdf_path <- file.path(figure_dir, sprintf("%s.pdf", dotplot_tag))
assert_absent_or_overwrite(
  c(
    full_marker_path,
    top_marker_path,
    summary_path,
    identity_map_path,
    png_path,
    pdf_path
  ),
  overwrite_outputs
)

# ---- work ----

sobj <- readRDS(input_path)
if (!inherits(sobj, "Seurat")) {
  stop(
    "Input RDS does not contain a Seurat object: ",
    input_path,
    call. = FALSE
  )
}
object_branch_tag <- sobj@misc$clustering$branch_tag
if (!is.null(object_branch_tag)) {
  assert_scalar_character(object_branch_tag, "sobj@misc$clustering$branch_tag")
  if (!identical(object_branch_tag, branch_tag)) {
    stop(
      "Input object branch tag mismatch: expected ",
      branch_tag,
      ", found ",
      object_branch_tag,
      call. = FALSE
    )
  }
}
cluster_column <- sprintf(
  "cluster_%s_dims%d_res%s",
  branch_tag,
  dims,
  resolution
)
if (!cluster_column %in% colnames(sobj@meta.data)) {
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
cluster_levels <- cluster_levels_for_labels(cluster_values)

available_assays <- names(sobj@assays)
if (is.null(assay_arg)) {
  assay <- if ("RNA" %in% available_assays) {
    "RNA"
  } else {
    SeuratObject::DefaultAssay(sobj)
  }
} else {
  assay <- assay_arg
}
if (!assay %in% available_assays) {
  stop(
    "Missing assay '",
    assay,
    "'. Available assays: ",
    paste(available_assays, collapse = ", "),
    call. = FALSE
  )
}
available_layers <- SeuratObject::Layers(sobj[[assay]])
missing_layers <- setdiff(c(expression_layer, counts_layer), available_layers)
if (length(missing_layers) > 0L) {
  stop(
    "Missing assay layer(s) in assay ",
    assay,
    ": ",
    paste(missing_layers, collapse = ", "),
    ". Available layers: ",
    paste(available_layers, collapse = ", "),
    call. = FALSE
  )
}

if (confirm_no_merge) {
  identity_map <- confirmed_identity_map(cluster_levels, cluster_values)
  decision_source <- "confirmed_no_merge"
} else {
  identity_map <- read_identity_map(
    cluster_map_path,
    cluster_levels,
    cluster_values
  )
  decision_source <- "cluster_map"
}
identity_levels <- unique(identity_map$marker_identity)
marker_identities <- marker_identity_values(cluster_values, identity_map)
validate_marker_identities(marker_identities, identity_levels, min_cells_group)
identity_map$decision_source <- decision_source
identity_map$input_path <- input_path
identity_map$cluster_column <- cluster_column
identity_map$assay <- assay
identity_map$expression_layer <- expression_layer
identity_map$counts_layer <- counts_layer

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
  return.thresh = 1,
  verbose = FALSE
)
ranked_markers <- prepare_marker_table(markers, identity_levels, min_diff_pct)
top_markers <- select_top_markers(
  ranked_markers,
  identity_levels,
  top_n,
  min_diff_pct
)
marker_summary <- summarise_markers(
  identity_map,
  identity_levels,
  ranked_markers,
  top_markers,
  decision_source
)
gene_order <- unique(top_markers$gene)
plot_data <- build_dotplot_data(
  sobj = sobj,
  assay = assay,
  expression_layer = expression_layer,
  counts_layer = counts_layer,
  genes = gene_order,
  marker_identities = marker_identities,
  identity_levels = identity_levels
)
plot <- marker_dotplot(plot_data, expression_layer)
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
notebook_png_path <- link_notebook_png(png_path)

message("Wrote full marker CSV: ", full_marker_path)
message("Wrote top marker CSV: ", top_marker_path)
message("Wrote marker summary CSV: ", summary_path)
message("Wrote marker identity map CSV: ", identity_map_path)
message("Wrote marker dotplot PNG: ", png_path)
message("Wrote marker dotplot PDF: ", pdf_path)
message("Linked notebook figure: ", notebook_png_path)
