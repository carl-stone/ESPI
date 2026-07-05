# Cluster grid summaries and supplemental plots.

DEFAULT_CLUSTER_ELBOW_N <- 20L
DEFAULT_CLUSTER_DIMS <- c(20L, 30L, 50L)
DEFAULT_CLUSTER_RESOLUTIONS <- c(0.3, 0.5, 0.8)
DEFAULT_CLUSTER_NORMALIZATIONS <- c("log1p", "pflog")
DEFAULT_CLUSTER_FILTER_STATES <- c(FALSE, TRUE)
DEFAULT_CLUSTER_REFERENCE_COLUMN <- "cluster_pflog_filter_cc_dims50_res0.3"
SMALL_CLUSTER_CELL_THRESHOLD <- 50L
CLUSTREE_GRID_WIDTH <- 16
CLUSTREE_GRID_HEIGHT <- 12
UMAP_SWEEP_PANEL_WIDTH <- 5
UMAP_SWEEP_HEIGHT <- 5
CLUSTER_COLUMN_DIMS_CAPTURE <- 2L
CLUSTER_COLUMN_RES_CAPTURE <- 3L
MIN_CLUSTREE_RESOLUTION_COLUMNS <- 2L
CLUSTREE_GRID_OUT_TAG <- "cluster_grid_clustree_12_panel"
REPRESENTATIVE_NORMALIZATION <- "pflog"
REPRESENTATIVE_FILTERED_CC <- TRUE
REPRESENTATIVE_DIMS <- 50L
UMAP_SWEEP_OUT_TAG <- "umap_resolution_sweep_pflog_filter_cc_dims50"
CLUSTER_GRID_STABILITY_SUMMARY_FILENAME <- "cluster_grid_stability_summary.tsv"
CLUSTER_GRID_PAIRWISE_STABILITY_FILENAME <- "cluster_grid_pairwise_stability.tsv"
ENTROPY_LOG_BASE <- 2
CLUSTREE_STABLE_CHILD_FRACTION <- 0.8


cluster_grid_branch_tag <- function(normalization, filtered_cell_cycle) {
  cc_tag <- if (isTRUE(filtered_cell_cycle)) {
    "filter_cc"
  } else {
    "no_filter_cc"
  }
  branch_tag <- sprintf("%s_%s", normalization, cc_tag)
  if (!grepl("^[A-Za-z0-9_]+$", branch_tag)) {
    stop("Unsafe clustering branch tag: ", branch_tag, call. = FALSE)
  }
  branch_tag
}

cluster_grid_branch_label <- function(normalization, filtered_cell_cycle) {
  cc_label <- if (isTRUE(filtered_cell_cycle)) {
    "CC-HVG filtered"
  } else {
    "CC-HVG retained"
  }
  norm_label <- if (identical(normalization, "pflog")) {
    "PFlog"
  } else {
    normalization
  }
  sprintf("%s, %s", norm_label, cc_label)
}

cluster_grid_res_tag <- function(resolution) {
  format(resolution, trim = TRUE, scientific = FALSE)
}

cluster_grid_branches <- function(
  normalizations = DEFAULT_CLUSTER_NORMALIZATIONS,
  filter_states = DEFAULT_CLUSTER_FILTER_STATES
) {
  rows <- list()
  idx <- 1L
  for (normalization in normalizations) {
    for (filtered_cell_cycle in filter_states) {
      rows[[idx]] <- data.frame(
        normalization = normalization,
        filtered_cell_cycle = filtered_cell_cycle,
        branch_tag = cluster_grid_branch_tag(
          normalization,
          filtered_cell_cycle
        ),
        branch_label = cluster_grid_branch_label(
          normalization,
          filtered_cell_cycle
        ),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

cluster_grid_object_path <- function(
  branch_tag,
  elbow_n = DEFAULT_CLUSTER_ELBOW_N
) {
  file.path(
    CURRENT_OBJECT_DIR,
    sprintf("cluster_%s_elbow%d.rds", branch_tag, elbow_n)
  )
}

load_cluster_grid_objects <- function(
  elbow_n = DEFAULT_CLUSTER_ELBOW_N,
  normalizations = DEFAULT_CLUSTER_NORMALIZATIONS,
  filter_states = DEFAULT_CLUSTER_FILTER_STATES
) {
  branches <- cluster_grid_branches(normalizations, filter_states)
  paths <- vapply(
    branches$branch_tag,
    cluster_grid_object_path,
    character(1),
    elbow_n = elbow_n
  )
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0) {
    stop(
      "Missing clustered object(s): ",
      paste(missing_paths, collapse = ", "),
      call. = FALSE
    )
  }
  objects <- lapply(paths, readRDS)
  names(objects) <- branches$branch_tag
  list(branches = branches, paths = paths, objects = objects)
}

parse_cluster_grid_column <- function(column, branch_tag) {
  pattern <- sprintf("^cluster_%s_dims([0-9]+)_res(.+)$", branch_tag)
  hit <- regexec(pattern, column, perl = TRUE)
  parts <- regmatches(column, hit)[[1]]
  expected_parts <- 3L
  if (length(parts) != expected_parts) {
    stop("Cannot parse cluster column: ", column, call. = FALSE)
  }
  data.frame(
    cluster_column = column,
    dims = as.integer(parts[[CLUSTER_COLUMN_DIMS_CAPTURE]]),
    resolution = as.numeric(parts[[CLUSTER_COLUMN_RES_CAPTURE]]),
    stringsAsFactors = FALSE
  )
}

cluster_grid_candidate_columns <- function(sobj, branch_tag) {
  candidate_names <- sobj@misc$clustering$candidate_names
  if (is.null(candidate_names) || length(candidate_names) == 0) {
    candidate_names <- grep(
      sprintf("^cluster_%s_dims[0-9]+_res", branch_tag),
      colnames(sobj@meta.data),
      value = TRUE
    )
  }
  missing_columns <- setdiff(candidate_names, colnames(sobj@meta.data))
  if (length(missing_columns) > 0) {
    stop(
      "Missing candidate cluster column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(candidate_names) == 0) {
    stop(
      "No candidate cluster columns found for branch: ",
      branch_tag,
      call. = FALSE
    )
  }
  parsed <- lapply(
    candidate_names,
    parse_cluster_grid_column,
    branch_tag = branch_tag
  )
  do.call(rbind, parsed)
}

cluster_grid_labels <- function(sobj, column) {
  labels <- sobj@meta.data[[column]]
  names(labels) <- rownames(sobj@meta.data)
  labels
}

best_jaccard_to_reference <- function(labels, reference_labels) {
  contingency <- table(labels, reference_labels)
  cluster_sizes <- rowSums(contingency)
  reference_sizes <- colSums(contingency)
  best <- vapply(
    seq_len(nrow(contingency)),
    function(i) {
      intersections <- contingency[i, ]
      unions <- cluster_sizes[[i]] + reference_sizes - intersections
      max(intersections / unions)
    },
    numeric(1)
  )
  list(
    mean_best_jaccard = stats::weighted.mean(best, cluster_sizes),
    min_best_jaccard = min(best)
  )
}

cluster_grid_size_summary <- function(
  labels,
  small_cluster_cell_threshold = SMALL_CLUSTER_CELL_THRESHOLD
) {
  cluster_sizes <- table(labels)
  small_clusters <- cluster_sizes[cluster_sizes < small_cluster_cell_threshold]
  data.frame(
    n_cells = length(labels),
    n_clusters = length(cluster_sizes),
    min_cluster_n = min(cluster_sizes),
    q25_cluster_n = as.numeric(stats::quantile(
      cluster_sizes,
      0.25,
      names = FALSE
    )),
    median_cluster_n = as.numeric(stats::median(cluster_sizes)),
    max_cluster_n = max(cluster_sizes),
    n_small_clusters = length(small_clusters),
    n_cells_in_small_clusters = sum(small_clusters),
    fraction_cells_in_small_clusters = sum(small_clusters) /
      length(labels),
    stringsAsFactors = FALSE
  )
}

collect_cluster_grid_candidates <- function(grid) {
  rows <- list()
  labels <- list()
  idx <- 1L
  for (branch_idx in seq_len(nrow(grid$branches))) {
    branch_info <- grid$branches[branch_idx, ]
    sobj <- grid$objects[[branch_info$branch_tag]]
    candidates <- cluster_grid_candidate_columns(sobj, branch_info$branch_tag)
    for (candidate_idx in seq_len(nrow(candidates))) {
      candidate <- candidates[candidate_idx, ]
      rows[[idx]] <- data.frame(
        normalization = branch_info$normalization,
        filtered_cell_cycle = branch_info$filtered_cell_cycle,
        branch_tag = branch_info$branch_tag,
        dims = candidate$dims,
        resolution = candidate$resolution,
        cluster_column = candidate$cluster_column,
        stringsAsFactors = FALSE
      )
      labels[[candidate$cluster_column]] <- cluster_grid_labels(
        sobj,
        candidate$cluster_column
      )
      idx <- idx + 1L
    }
  }
  metadata <- do.call(rbind, rows)
  ordered <- order(
    metadata$normalization,
    metadata$filtered_cell_cycle,
    metadata$dims,
    metadata$resolution
  )
  ordered_metadata <- metadata[ordered, ]
  rownames(ordered_metadata) <- NULL
  list(
    metadata = ordered_metadata,
    labels = labels[ordered_metadata$cluster_column]
  )
}

cluster_grid_entropy <- function(labels) {
  proportions <- as.numeric(table(labels)) / length(labels)
  -sum(proportions * log(proportions, base = ENTROPY_LOG_BASE))
}

cluster_grid_pairwise_metrics <- function(labels_a, labels_b) {
  if (!setequal(names(labels_a), names(labels_b))) {
    stop("Cell names differ between pairwise clusterings.", call. = FALSE)
  }
  ordered_a <- labels_a[sort(names(labels_a))]
  ordered_b <- labels_b[names(ordered_a)]
  contingency <- table(
    as.character(ordered_a),
    as.character(ordered_b)
  )
  total <- sum(contingency)
  joint <- contingency / total
  row_probability <- rowSums(joint)
  col_probability <- colSums(joint)
  independent <- outer(row_probability, col_probability)
  nonzero <- joint > 0
  mutual_information <- sum(
    joint[nonzero] *
      log(joint[nonzero] / independent[nonzero], base = ENTROPY_LOG_BASE)
  )
  entropy_a <- cluster_grid_entropy(ordered_a)
  entropy_b <- cluster_grid_entropy(ordered_b)
  entropy_sum <- entropy_a + entropy_b
  variation_of_information <- entropy_sum - 2 * mutual_information
  normalized_mutual_information <- if (entropy_sum > 0) {
    2 * mutual_information / entropy_sum
  } else {
    NA_real_
  }
  jaccard_ab <- best_jaccard_to_reference(ordered_a, ordered_b)
  jaccard_ba <- best_jaccard_to_reference(ordered_b, ordered_a)

  data.frame(
    ari = mclust::adjustedRandIndex(
      as.character(ordered_a),
      as.character(ordered_b)
    ),
    normalized_mutual_information = max(
      0,
      min(1, normalized_mutual_information)
    ),
    variation_of_information = max(0, variation_of_information),
    mean_best_jaccard_a_to_b = jaccard_ab$mean_best_jaccard,
    min_best_jaccard_a_to_b = jaccard_ab$min_best_jaccard,
    mean_best_jaccard_b_to_a = jaccard_ba$mean_best_jaccard,
    min_best_jaccard_b_to_a = jaccard_ba$min_best_jaccard,
    mean_best_jaccard_bidirectional = mean(c(
      jaccard_ab$mean_best_jaccard,
      jaccard_ba$mean_best_jaccard
    )),
    min_best_jaccard_bidirectional = min(
      jaccard_ab$min_best_jaccard,
      jaccard_ba$min_best_jaccard
    ),
    stringsAsFactors = FALSE
  )
}

cluster_grid_neighbor_axis <- function(candidate_a, candidate_b) {
  changed_axes <- c(
    normalization = candidate_a$normalization != candidate_b$normalization,
    filtered_cell_cycle = candidate_a$filtered_cell_cycle !=
      candidate_b$filtered_cell_cycle,
    dims = candidate_a$dims != candidate_b$dims,
    resolution = candidate_a$resolution != candidate_b$resolution
  )
  if (sum(changed_axes) != 1L) {
    return(NA_character_)
  }
  names(changed_axes)[changed_axes]
}

cluster_grid_is_local_neighbor <- function(candidate_a, candidate_b) {
  axis <- cluster_grid_neighbor_axis(candidate_a, candidate_b)
  if (is.na(axis)) {
    return(FALSE)
  }
  if (axis == "dims") {
    dim_positions <- match(
      c(candidate_a$dims, candidate_b$dims),
      DEFAULT_CLUSTER_DIMS
    )
    return(abs(diff(dim_positions)) == 1L)
  }
  if (axis == "resolution") {
    resolution_positions <- match(
      c(candidate_a$resolution, candidate_b$resolution),
      DEFAULT_CLUSTER_RESOLUTIONS
    )
    return(abs(diff(resolution_positions)) == 1L)
  }
  TRUE
}

cluster_grid_pairwise_stability <- function(candidate_map) {
  metadata <- candidate_map$metadata
  labels <- candidate_map$labels
  rows <- list()
  idx <- 1L
  for (candidate_a_idx in seq_len(nrow(metadata) - 1L)) {
    candidate_a <- metadata[candidate_a_idx, ]
    for (candidate_b_idx in seq(
      from = candidate_a_idx + 1L,
      to = nrow(metadata)
    )) {
      candidate_b <- metadata[candidate_b_idx, ]
      metrics <- cluster_grid_pairwise_metrics(
        labels[[candidate_a$cluster_column]],
        labels[[candidate_b$cluster_column]]
      )
      neighbor_axis <- cluster_grid_neighbor_axis(candidate_a, candidate_b)
      rows[[idx]] <- data.frame(
        cluster_column_a = candidate_a$cluster_column,
        normalization_a = candidate_a$normalization,
        filtered_cell_cycle_a = candidate_a$filtered_cell_cycle,
        dims_a = candidate_a$dims,
        resolution_a = candidate_a$resolution,
        cluster_column_b = candidate_b$cluster_column,
        normalization_b = candidate_b$normalization,
        filtered_cell_cycle_b = candidate_b$filtered_cell_cycle,
        dims_b = candidate_b$dims,
        resolution_b = candidate_b$resolution,
        neighbor_axis = neighbor_axis,
        is_local_neighbor = cluster_grid_is_local_neighbor(
          candidate_a,
          candidate_b
        ),
        metrics,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

summarize_cluster_grid_pairwise <- function(metadata, pairwise) {
  rows <- lapply(metadata$cluster_column, function(column) {
    local_pairs <- pairwise[
      pairwise$is_local_neighbor &
        (pairwise$cluster_column_a == column |
          pairwise$cluster_column_b == column),
    ]
    data.frame(
      cluster_column = column,
      local_pairwise_n = nrow(local_pairs),
      local_mean_ari = if (nrow(local_pairs) > 0) {
        mean(local_pairs$ari)
      } else {
        NA_real_
      },
      local_min_ari = if (nrow(local_pairs) > 0) {
        min(local_pairs$ari)
      } else {
        NA_real_
      },
      local_mean_nmi = if (nrow(local_pairs) > 0) {
        mean(local_pairs$normalized_mutual_information)
      } else {
        NA_real_
      },
      local_min_nmi = if (nrow(local_pairs) > 0) {
        min(local_pairs$normalized_mutual_information)
      } else {
        NA_real_
      },
      local_mean_vi = if (nrow(local_pairs) > 0) {
        mean(local_pairs$variation_of_information)
      } else {
        NA_real_
      },
      local_max_vi = if (nrow(local_pairs) > 0) {
        max(local_pairs$variation_of_information)
      } else {
        NA_real_
      },
      local_mean_best_jaccard = if (nrow(local_pairs) > 0) {
        mean(local_pairs$mean_best_jaccard_bidirectional)
      } else {
        NA_real_
      },
      local_min_best_jaccard = if (nrow(local_pairs) > 0) {
        min(local_pairs$min_best_jaccard_bidirectional)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

cluster_grid_clustree_metrics <- function(labels, next_labels) {
  if (!setequal(names(labels), names(next_labels))) {
    stop("Cell names differ between clustree clusterings.", call. = FALSE)
  }
  ordered_labels <- labels[sort(names(labels))]
  ordered_next_labels <- next_labels[names(ordered_labels)]
  contingency <- table(
    as.character(ordered_labels),
    as.character(ordered_next_labels)
  )
  cluster_sizes <- rowSums(contingency)
  largest_child_fraction <- apply(contingency, 1, max) / cluster_sizes
  child_counts <- rowSums(contingency > 0)
  child_entropy <- apply(contingency, 1, function(counts) {
    proportions <- counts[counts > 0] / sum(counts)
    -sum(proportions * log(proportions, base = ENTROPY_LOG_BASE))
  })

  data.frame(
    clustree_weighted_largest_child_fraction = stats::weighted.mean(
      largest_child_fraction,
      cluster_sizes
    ),
    clustree_min_largest_child_fraction = min(largest_child_fraction),
    clustree_weighted_child_entropy = stats::weighted.mean(
      child_entropy,
      cluster_sizes
    ),
    clustree_mean_child_count = mean(child_counts),
    clustree_max_child_count = max(child_counts),
    clustree_splitting_clusters = sum(
      largest_child_fraction < CLUSTREE_STABLE_CHILD_FRACTION
    ),
    stringsAsFactors = FALSE
  )
}

summarize_cluster_grid_clustree <- function(candidate_map) {
  metadata <- candidate_map$metadata
  labels <- candidate_map$labels
  rows <- vector("list", nrow(metadata))
  for (candidate_idx in seq_len(nrow(metadata))) {
    candidate <- metadata[candidate_idx, ]
    next_resolutions <- DEFAULT_CLUSTER_RESOLUTIONS[
      DEFAULT_CLUSTER_RESOLUTIONS > candidate$resolution
    ]
    next_resolution <- if (length(next_resolutions) > 0) {
      min(next_resolutions)
    } else {
      NA_real_
    }
    next_row <- metadata[
      metadata$branch_tag == candidate$branch_tag &
        metadata$dims == candidate$dims &
        metadata$resolution == next_resolution,
    ]
    if (is.na(next_resolution) || nrow(next_row) != 1L) {
      metrics <- data.frame(
        clustree_next_resolution = NA_real_,
        clustree_weighted_largest_child_fraction = NA_real_,
        clustree_min_largest_child_fraction = NA_real_,
        clustree_weighted_child_entropy = NA_real_,
        clustree_mean_child_count = NA_real_,
        clustree_max_child_count = NA_integer_,
        clustree_splitting_clusters = NA_integer_,
        stringsAsFactors = FALSE
      )
    } else {
      metrics <- cbind(
        data.frame(
          clustree_next_resolution = next_resolution,
          stringsAsFactors = FALSE
        ),
        cluster_grid_clustree_metrics(
          labels[[candidate$cluster_column]],
          labels[[next_row$cluster_column]]
        )
      )
    }
    rows[[candidate_idx]] <- data.frame(
      cluster_column = candidate$cluster_column,
      metrics,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_cluster_grid_column <- function(
  sobj,
  branch_info,
  column,
  dims,
  resolution,
  reference_labels,
  reference_column,
  small_cluster_cell_threshold = SMALL_CLUSTER_CELL_THRESHOLD
) {
  labels <- cluster_grid_labels(sobj, column)
  if (!setequal(names(labels), names(reference_labels))) {
    stop(
      "Cell names differ between ",
      column,
      " and ",
      reference_column,
      call. = FALSE
    )
  }
  ordered_labels <- labels[names(reference_labels)]
  size_summary <- cluster_grid_size_summary(
    ordered_labels,
    small_cluster_cell_threshold = small_cluster_cell_threshold
  )
  jaccard <- best_jaccard_to_reference(ordered_labels, reference_labels)

  data.frame(
    normalization = branch_info$normalization,
    filtered_cell_cycle = branch_info$filtered_cell_cycle,
    branch_tag = branch_info$branch_tag,
    dims = dims,
    resolution = resolution,
    cluster_column = column,
    size_summary,
    ari_vs_reference = mclust::adjustedRandIndex(
      as.character(ordered_labels),
      as.character(reference_labels)
    ),
    mean_best_jaccard_to_reference = jaccard$mean_best_jaccard,
    min_best_jaccard_to_reference = jaccard$min_best_jaccard,
    reference_column = reference_column,
    stringsAsFactors = FALSE
  )
}

#' Write a supplemental clustering grid summary table.
#'
#' Computes cluster counts, cluster-size summaries, adjusted Rand index (ARI),
#' and best-overlap Jaccard summaries for every clustered branch in the current
#' grid. The default reference is the real cluster metadata column for PFlog,
#' cell-cycle-filtered HVGs, 50 PCs, and resolution 0.3:
#' `cluster_pflog_filter_cc_dims50_res0.3`.
#'
#' @param reference_column Character cluster metadata column used as the ARI
#'   reference.
#' @param elbow_n Integer elbow value embedded in clustered object filenames.
#' @param out_path Output TSV path. Defaults to
#'   `TABLE_DIR/cluster/cluster_grid_summary.tsv`.
#' @param small_cluster_cell_threshold Integer threshold for counting small
#'   clusters in the supplemental table.
#'
#' @return The summary `data.frame`, invisibly.
#' @export
write_cluster_grid_summary <- function(
  reference_column = DEFAULT_CLUSTER_REFERENCE_COLUMN,
  elbow_n = DEFAULT_CLUSTER_ELBOW_N,
  out_path = file.path(TABLE_DIR, "cluster", "cluster_grid_summary.tsv"),
  small_cluster_cell_threshold = SMALL_CLUSTER_CELL_THRESHOLD
) {
  grid <- load_cluster_grid_objects(elbow_n = elbow_n)
  reference_hits <- vapply(
    grid$objects,
    function(sobj) {
      reference_column %in% colnames(sobj@meta.data)
    },
    logical(1)
  )
  if (sum(reference_hits) != 1L) {
    stop(
      "Expected exactly one reference column match for ",
      reference_column,
      "; found ",
      sum(reference_hits),
      call. = FALSE
    )
  }
  reference_sobj <- grid$objects[[which(reference_hits)]]
  reference_labels <- cluster_grid_labels(reference_sobj, reference_column)

  rows <- list()
  idx <- 1L
  for (branch_idx in seq_len(nrow(grid$branches))) {
    branch_info <- grid$branches[branch_idx, ]
    sobj <- grid$objects[[branch_info$branch_tag]]
    candidates <- cluster_grid_candidate_columns(sobj, branch_info$branch_tag)
    for (candidate_idx in seq_len(nrow(candidates))) {
      candidate <- candidates[candidate_idx, ]
      rows[[idx]] <- summarize_cluster_grid_column(
        sobj = sobj,
        branch_info = branch_info,
        column = candidate$cluster_column,
        dims = candidate$dims,
        resolution = candidate$resolution,
        reference_labels = reference_labels,
        reference_column = reference_column,
        small_cluster_cell_threshold = small_cluster_cell_threshold
      )
      idx <- idx + 1L
    }
  }

  summary <- do.call(rbind, rows)
  ordered_summary <- summary[
    order(
      summary$normalization,
      summary$filtered_cell_cycle,
      summary$dims,
      summary$resolution
    ),
  ]
  rownames(ordered_summary) <- NULL
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    ordered_summary,
    out_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
  invisible(ordered_summary)
}

#' Write label-blind clustering stability tables.
#'
#' Computes two supplemental stability outputs from the clustering grid without
#' using condition labels. The summary table has one row per candidate clustering
#' with cluster-size summaries, local pairwise stability to neighboring grid
#' settings, and clustree-style split stability to the next resolution. The
#' pairwise table has one row per candidate pair with ARI, normalized mutual
#' information, variation of information, and bidirectional best-overlap Jaccard.
#'
#' @param elbow_n Integer elbow value embedded in clustered object filenames.
#' @param summary_out_path Output path for the one-row-per-candidate stability
#'   summary.
#' @param pairwise_out_path Output path for the one-row-per-pair stability table.
#' @param small_cluster_cell_threshold Integer threshold for counting small
#'   clusters in the stability summary.
#'
#' @return A list with `summary` and `pairwise` data frames, invisibly.
#' @export
write_cluster_grid_stability_tables <- function(
  elbow_n = DEFAULT_CLUSTER_ELBOW_N,
  summary_out_path = file.path(
    TABLE_DIR,
    "cluster",
    CLUSTER_GRID_STABILITY_SUMMARY_FILENAME
  ),
  pairwise_out_path = file.path(
    TABLE_DIR,
    "cluster",
    CLUSTER_GRID_PAIRWISE_STABILITY_FILENAME
  ),
  small_cluster_cell_threshold = SMALL_CLUSTER_CELL_THRESHOLD
) {
  grid <- load_cluster_grid_objects(elbow_n = elbow_n)
  candidate_map <- collect_cluster_grid_candidates(grid)
  metadata <- candidate_map$metadata
  size_rows <- lapply(metadata$cluster_column, function(column) {
    cluster_grid_size_summary(
      candidate_map$labels[[column]],
      small_cluster_cell_threshold = small_cluster_cell_threshold
    )
  })
  size_summary <- do.call(rbind, size_rows)
  pairwise <- cluster_grid_pairwise_stability(candidate_map)
  pairwise_summary <- summarize_cluster_grid_pairwise(metadata, pairwise)
  clustree_summary <- summarize_cluster_grid_clustree(candidate_map)

  if (!identical(metadata$cluster_column, pairwise_summary$cluster_column)) {
    stop("Pairwise stability summary order mismatch.", call. = FALSE)
  }
  if (!identical(metadata$cluster_column, clustree_summary$cluster_column)) {
    stop("Clustree stability summary order mismatch.", call. = FALSE)
  }

  summary <- cbind(
    metadata,
    size_summary,
    pairwise_summary[
      setdiff(colnames(pairwise_summary), "cluster_column")
    ],
    clustree_summary[
      setdiff(colnames(clustree_summary), "cluster_column")
    ]
  )
  rownames(summary) <- NULL

  dir.create(dirname(summary_out_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    summary,
    summary_out_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
  dir.create(dirname(pairwise_out_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    pairwise,
    pairwise_out_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )

  invisible(list(summary = summary, pairwise = pairwise))
}

#' Save a 12-panel clustree grid across clustering branches.
#'
#' The default grid shows two normalization methods, two cell-cycle-HVG policies,
#' and three PC counts. Each panel contains the resolution sweep for one branch
#' and PC count.
#'
#' @param elbow_n Integer elbow value embedded in clustered object filenames.
#' @param dims Integer PC counts to plot.
#' @param out_tag Filename tag used for PNG and PDF outputs.
#'
#' @return `invisible(NULL)`.
#' @export
splot_cluster_grid_clustree <- function(
  elbow_n = DEFAULT_CLUSTER_ELBOW_N,
  dims = DEFAULT_CLUSTER_DIMS,
  out_tag = CLUSTREE_GRID_OUT_TAG
) {
  grid <- load_cluster_grid_objects(elbow_n = elbow_n)
  plots <- list()
  idx <- 1L
  for (d in dims) {
    for (branch_idx in seq_len(nrow(grid$branches))) {
      branch_info <- grid$branches[branch_idx, ]
      sobj <- grid$objects[[branch_info$branch_tag]]
      prefix <- sprintf("cluster_%s_dims%d_res", branch_info$branch_tag, d)
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
      plots[[idx]] <- clustree::clustree(
        cluster_data,
        prefix = prefix,
        node_text_size = 2
      ) +
        ggplot2::guides(edge_colour = "none") +
        ggplot2::ggtitle(sprintf("%s; %d PCs", branch_info$branch_label, d)) +
        ggplot2::theme(
          legend.position = "none",
          plot.title = ggplot2::element_text(size = 10)
        )
      idx <- idx + 1L
    }
  }

  out_dir <- file.path(FIGURE_DIR, "cluster")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plot <- patchwork::wrap_plots(plots, ncol = nrow(grid$branches))
  ggplot2::ggsave(
    file.path(out_dir, sprintf("%s.png", out_tag)),
    plot,
    width = CLUSTREE_GRID_WIDTH,
    height = CLUSTREE_GRID_HEIGHT
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("%s.pdf", out_tag)),
    plot,
    width = CLUSTREE_GRID_WIDTH,
    height = CLUSTREE_GRID_HEIGHT
  )

  invisible(NULL)
}

#' Save a representative UMAP resolution sweep.
#'
#' The default representative branch is PFlog with cell-cycle-filtered HVGs at 30
#' PCs, matching the current reference clustering branch.
#'
#' @param normalization Normalization branch.
#' @param filtered_cell_cycle Logical indicating whether cell-cycle genes were
#'   removed from HVGs.
#' @param dims Integer PC count whose UMAP reduction should be plotted.
#' @param elbow_n Integer elbow value embedded in clustered object filenames.
#' @param resolutions Numeric resolutions to plot.
#' @param out_tag Filename tag used for PNG and PDF outputs.
#'
#' @return `invisible(NULL)`.
#' @export
splot_umap_resolution_sweep <- function(
  normalization = REPRESENTATIVE_NORMALIZATION,
  filtered_cell_cycle = REPRESENTATIVE_FILTERED_CC,
  dims = REPRESENTATIVE_DIMS,
  elbow_n = DEFAULT_CLUSTER_ELBOW_N,
  resolutions = DEFAULT_CLUSTER_RESOLUTIONS,
  out_tag = UMAP_SWEEP_OUT_TAG
) {
  branch_tag <- cluster_grid_branch_tag(normalization, filtered_cell_cycle)
  path <- cluster_grid_object_path(branch_tag, elbow_n = elbow_n)
  if (!file.exists(path)) {
    stop("Missing clustered object: ", path, call. = FALSE)
  }
  sobj <- readRDS(path)
  reduction <- sprintf("umap_%s_dims%d", branch_tag, dims)
  if (!reduction %in% names(sobj@reductions)) {
    stop("Missing UMAP reduction: ", reduction, call. = FALSE)
  }

  plots <- lapply(resolutions, function(resolution) {
    column <- sprintf(
      "cluster_%s_dims%d_res%s",
      branch_tag,
      dims,
      cluster_grid_res_tag(resolution)
    )
    if (!column %in% colnames(sobj@meta.data)) {
      stop("Missing cluster metadata column: ", column, call. = FALSE)
    }
    Seurat::DimPlot(
      sobj,
      reduction = reduction,
      group.by = column,
      label = TRUE,
      pt.size = 0.25
    ) +
      ggplot2::ggtitle(sprintf("res %s", resolution)) +
      ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 12)
      )
  })

  out_dir <- file.path(FIGURE_DIR, "cluster")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plot <- patchwork::wrap_plots(plots, ncol = length(resolutions)) +
    patchwork::plot_annotation(
      title = sprintf(
        "%s, %s, %d PCs",
        if (identical(normalization, "pflog")) "PFlog" else normalization,
        if (isTRUE(filtered_cell_cycle)) {
          "CC-HVG filtered"
        } else {
          "CC-HVG retained"
        },
        dims
      )
    )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("%s.png", out_tag)),
    plot,
    width = UMAP_SWEEP_PANEL_WIDTH * length(resolutions),
    height = UMAP_SWEEP_HEIGHT
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("%s.pdf", out_tag)),
    plot,
    width = UMAP_SWEEP_PANEL_WIDTH * length(resolutions),
    height = UMAP_SWEEP_HEIGHT
  )

  invisible(NULL)
}
