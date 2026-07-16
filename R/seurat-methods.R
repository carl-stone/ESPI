# Nonstandard Seurat preprocessing and clustering-grid methods.

#' Normalize with log1p and compute PCA for one Seurat object.
#'
#' Writes results to `sobj[["pca"]]` and records preprocessing provenance.
#'
#' @param sobj Seurat object with `counts` layer populated and `VariableFeatures`
#'   already selected.
#' @param n_pcs Integer number of PCs to compute. Default 50.
#'
#' @return `sobj` with `pca` reduction and `misc$preprocessing` populated.
#' @export
# ANALYSIS_OK[smuggled-default]: intentional package API default for PCA dimensionality.
run_log1p_pca <- function(sobj, n_pcs = 50) {
  stopifnot(length(n_pcs) == 1, is.numeric(n_pcs), is.finite(n_pcs), n_pcs > 0)

  features <- SeuratObject::VariableFeatures(sobj)
  if (length(features) == 0) {
    stop("VariableFeatures(sobj) is empty.", call. = FALSE)
  }
  filtered_cell_cycle <- isTRUE(sobj@misc$preprocessing$filtered_cell_cycle)

  sobj <- Seurat::NormalizeData(sobj)
  sobj <- Seurat::ScaleData(sobj, features = features)
  sobj <- Seurat::RunPCA(sobj, features = features, npcs = n_pcs)

  sobj@misc$active.reduction <- "pca"
  sobj@misc$preprocessing <- list(
    normalization = "log1p",
    pca_method = "Seurat::RunPCA",
    pca_source_layer = "scale.data",
    hvg_method = "Seurat::FindVariableFeatures(selection.method = 'vst')",
    n_variable_features = length(features),
    n_pcs = n_pcs,
    filtered_cell_cycle = filtered_cell_cycle
  )

  sobj
}

#' Normalize with PFlog (scclrR) and compute PCA on variable features.
#'
#' PFlog normalizes over the full retained feature set so that its per-cell
#' center is well defined; PCA is then computed on `VariableFeatures(sobj)` via
#' `scclrR::pca_matrix` using the full PFlog center.
#'
#' @param sobj Seurat object with `counts` layer populated and `VariableFeatures`
#'   already selected.
#' @param n_pcs Integer number of PCs to compute. Default 50.
#'
#' @return `sobj` with `pca` reduction and `misc$preprocessing` populated.
#' @export
# ANALYSIS_OK[smuggled-default]: intentional package API default for PCA dimensionality.
run_pflog_pca <- function(sobj, n_pcs = 50) {
  stopifnot(length(n_pcs) == 1, is.numeric(n_pcs), is.finite(n_pcs), n_pcs > 0)

  assay <- SeuratObject::DefaultAssay(sobj)
  if (!inherits(sobj[[assay]], "Assay5")) {
    sobj[[assay]] <- as(sobj[[assay]], Class = "Assay5")
  }
  features <- SeuratObject::VariableFeatures(sobj)
  if (length(features) == 0) {
    stop("VariableFeatures(sobj) is empty.", call. = FALSE)
  }
  filtered_cell_cycle <- isTRUE(sobj@misc$preprocessing$filtered_cell_cycle)
  sobj <- Seurat::NormalizeData(sobj)

  sobj <- scclrR::pflog(sobj)
  pflog_layer <- SeuratObject::LayerData(sobj[[assay]], layer = "pflog")
  center <- sobj[["pflog_center", drop = TRUE]]
  pca <- scclrR::pca_matrix(
    sparse = pflog_layer[features, , drop = FALSE],
    center = center,
    n.components = n_pcs,
    seed = SEED
  )

  pc_names <- paste0("PC_", seq_len(ncol(pca$scores)))
  colnames(pca$scores) <- pc_names
  colnames(pca$loadings) <- pc_names
  if (is.null(rownames(pca$scores))) {
    rownames(pca$scores) <- colnames(sobj)
  }
  if (is.null(rownames(pca$loadings))) {
    rownames(pca$loadings) <- features
  }

  sobj[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = pca$scores,
    loadings = pca$loadings,
    stdev = sqrt(pca$explained_variance),
    key = "PC_",
    assay = assay
  )

  sobj@misc$active.reduction <- "pca"
  sobj@misc$preprocessing <- list(
    normalization = "pflog",
    pca_method = "scclrR::pca_matrix",
    pca_source_layer = "pflog",
    pca_center_key = "pflog_center",
    hvg_method = "Seurat::FindVariableFeatures(selection.method = 'vst')",
    n_variable_features = length(features),
    n_pcs = n_pcs,
    filtered_cell_cycle = filtered_cell_cycle
  )

  sobj
}

#' Write clustering-grid summary and pairwise-stability tables.
#'
#' The supplied objects and branch table are already loaded by the caller. The
#' function writes the cluster summary, stability summary, and pairwise table
#' using the supplied output prefix, preserving the existing column order.
#'
#' @param grid_objects Named list of clustered Seurat objects, one per branch.
#' @param branch_table Data frame describing the branches. It must contain
#'   `normalization`, `filtered_cell_cycle`, and `branch_tag` columns.
#' @param output_prefix Filename prefix for the three TSV outputs.
#' @param table_dir Directory receiving the TSV outputs.
#'
#' @return A list with `summary` and `pairwise` data frames, invisibly.
#' @export
# ANALYSIS_OK[R026]: exported grid-table entrypoint is called by the frozen-regeneration phase script.
write_cluster_grid_tables <- function(
  grid_objects,
  branch_table,
  output_prefix,
  table_dir
) {
  write_grid_tables_impl(
    grid_objects = grid_objects,
    branch_table = branch_table,
    summary_path = file.path(table_dir, paste0(output_prefix, "_summary.tsv")),
    stability_path = file.path(
      table_dir,
      paste0(output_prefix, "_stability_summary.tsv")
    ),
    pairwise_path = file.path(
      table_dir,
      paste0(output_prefix, "_pairwise_stability.tsv")
    )
  )
}

# Fixed-interface implementation; phase scripts own object loading and paths.
# ANALYSIS_OK[R026]: private grid implementation is called by its exported entrypoint in this module.
write_grid_tables_impl <- function(
  grid_objects,
  branch_table,
  summary_path = NULL,
  stability_path = NULL,
  pairwise_path = NULL,
  reference_column = DEFAULT_CLUSTER_REFERENCE_COLUMN,
  small_cluster_cell_threshold = SMALL_CLUSTER_CELL_THRESHOLD
) {
  # ANALYSIS_OK[R026]: local parser helper is called by the grid implementation below.
  candidate_columns <- function(sobj, branch_tag) {
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
    parsed <- lapply(candidate_names, function(column) {
      pattern <- sprintf("^cluster_%s_dims([0-9]+)_res(.+)$", branch_tag)
      parts <- regmatches(column, regexec(pattern, column, perl = TRUE))[[1]]
      capture_group_count <- 3L
      if (length(parts) != capture_group_count) {
        stop("Cannot parse cluster column: ", column, call. = FALSE)
      }
      dims_capture_index <- 2L
      resolution_capture_index <- 3L
      data.frame(
        cluster_column = column,
        dims = as.integer(parts[[dims_capture_index]]),
        resolution = as.numeric(parts[[resolution_capture_index]]),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, parsed)
  }

  # ANALYSIS_OK[R026]: local metadata-label helper is called by the grid implementation below.
  labels_for <- function(sobj, column) {
    labels <- sobj@meta.data[[column]]
    names(labels) <- rownames(sobj@meta.data)
    labels
  }

  # ANALYSIS_OK[R026]: local similarity helper is called by the grid implementation below.
  best_jaccard <- function(labels, reference_labels) {
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

  # ANALYSIS_OK[R026]: local size-summary helper is called by the grid implementation below.
  size_summary <- function(labels) {
    cluster_sizes <- table(labels)
    small_clusters <- cluster_sizes[
      cluster_sizes < small_cluster_cell_threshold
    ]
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

  rows <- list()
  labels <- list()
  idx <- 1L
  for (branch_idx in seq_len(nrow(branch_table))) {
    branch_info <- branch_table[branch_idx, ]
    sobj <- grid_objects[[branch_info$branch_tag]]
    candidates <- candidate_columns(sobj, branch_info$branch_tag)
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
      labels[[candidate$cluster_column]] <- labels_for(
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
  # ANALYSIS_OK[R005]: reorder rows to the validated grid order without dropping any metadata rows.
  metadata <- metadata[ordered, ]
  rownames(metadata) <- NULL
  # ANALYSIS_OK[R005]: reorder label vectors to the same validated grid order without dropping cells.
  labels <- labels[metadata$cluster_column]
  candidate_map <- list(metadata = metadata, labels = labels)

  reference_hits <- vapply(
    grid_objects,
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
  reference_sobj <- grid_objects[[which(reference_hits)]]
  reference_labels <- labels_for(reference_sobj, reference_column)

  summary_rows <- lapply(seq_len(nrow(metadata)), function(candidate_idx) {
    candidate <- metadata[candidate_idx, ]
    labels <- candidate_map$labels[[candidate$cluster_column]]
    if (!setequal(names(labels), names(reference_labels))) {
      stop(
        "Cell names differ between ",
        candidate$cluster_column,
        " and ",
        reference_column,
        call. = FALSE
      )
    }
    ordered_labels <- labels[names(reference_labels)]
    size <- size_summary(ordered_labels)
    jaccard <- best_jaccard(ordered_labels, reference_labels)
    data.frame(
      normalization = candidate$normalization,
      filtered_cell_cycle = candidate$filtered_cell_cycle,
      branch_tag = candidate$branch_tag,
      dims = candidate$dims,
      resolution = candidate$resolution,
      cluster_column = candidate$cluster_column,
      size,
      ari_vs_reference = mclust::adjustedRandIndex(
        as.character(ordered_labels),
        as.character(reference_labels)
      ),
      mean_best_jaccard_to_reference = jaccard$mean_best_jaccard,
      min_best_jaccard_to_reference = jaccard$min_best_jaccard,
      reference_column = reference_column,
      stringsAsFactors = FALSE
    )
  })
  grid_summary <- do.call(rbind, summary_rows)
  rownames(grid_summary) <- NULL

  # ANALYSIS_OK[R026]: local entropy helper is called by pairwise grid calculations below.
  entropy <- function(labels) {
    proportions <- as.numeric(table(labels)) / length(labels)
    -sum(proportions * log(proportions, base = ENTROPY_LOG_BASE))
  }
  # ANALYSIS_OK[R026]: local pairwise-metrics helper is called by grid calculations below.
  pairwise_metrics <- function(labels_a, labels_b) {
    if (!setequal(names(labels_a), names(labels_b))) {
      stop("Cell names differ between pairwise clusterings.", call. = FALSE)
    }
    ordered_a <- labels_a[sort(names(labels_a))]
    ordered_b <- labels_b[names(ordered_a)]
    contingency <- table(as.character(ordered_a), as.character(ordered_b))
    total <- sum(contingency)
    joint <- contingency / total
    row_probability <- rowSums(joint)
    col_probability <- colSums(joint)
    independent <- outer(row_probability, col_probability)
    nonzero <- joint > 0
    mutual_information <- sum(
      joint[nonzero] *
        log(
          joint[nonzero] / independent[nonzero],
          base = ENTROPY_LOG_BASE
        )
    )
    entropy_a <- entropy(ordered_a)
    entropy_b <- entropy(ordered_b)
    entropy_sum <- entropy_a + entropy_b
    variation_of_information <- entropy_sum - 2 * mutual_information
    normalized_mutual_information <- if (entropy_sum > 0) {
      2 * mutual_information / entropy_sum
    } else {
      NA_real_
    }
    jaccard_ab <- best_jaccard(ordered_a, ordered_b)
    jaccard_ba <- best_jaccard(ordered_b, ordered_a)
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
  # ANALYSIS_OK[R026]: local neighbor-axis helper is called by pairwise grid calculations below.
  neighbor_axis <- function(candidate_a, candidate_b) {
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
  # ANALYSIS_OK[R026]: local neighbor predicate is called by pairwise grid calculations below.
  is_local_neighbor <- function(candidate_a, candidate_b) {
    axis <- neighbor_axis(candidate_a, candidate_b)
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

  pairwise_rows <- list()
  idx <- 1L
  if (nrow(metadata) > 1L) {
    for (candidate_a_idx in seq_len(nrow(metadata) - 1L)) {
      candidate_a <- metadata[candidate_a_idx, ]
      for (candidate_b_idx in seq(
        from = candidate_a_idx + 1L,
        to = nrow(metadata)
      )) {
        candidate_b <- metadata[candidate_b_idx, ]
        metrics <- pairwise_metrics(
          labels[[candidate_a$cluster_column]],
          labels[[candidate_b$cluster_column]]
        )
        axis <- neighbor_axis(candidate_a, candidate_b)
        pairwise_rows[[idx]] <- data.frame(
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
          neighbor_axis = axis,
          is_local_neighbor = is_local_neighbor(candidate_a, candidate_b),
          metrics,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  pairwise <- if (length(pairwise_rows) > 0L) {
    do.call(rbind, pairwise_rows)
  } else {
    data.frame()
  }

  pairwise_summary_rows <- lapply(metadata$cluster_column, function(column) {
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
  pairwise_summary <- do.call(rbind, pairwise_summary_rows)

  # ANALYSIS_OK[R026]: local clustree helper is called by the grid summary below.
  clustree_metrics <- function(labels_a, labels_b) {
    if (!setequal(names(labels_a), names(labels_b))) {
      stop("Cell names differ between clustree clusterings.", call. = FALSE)
    }
    ordered_a <- labels_a[sort(names(labels_a))]
    ordered_b <- labels_b[names(ordered_a)]
    contingency <- table(as.character(ordered_a), as.character(ordered_b))
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
  clustree_rows <- lapply(seq_len(nrow(metadata)), function(candidate_idx) {
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
        clustree_metrics(
          labels[[candidate$cluster_column]],
          labels[[next_row$cluster_column]]
        )
      )
    }
    data.frame(cluster_column = candidate$cluster_column, metrics)
  })
  clustree_summary <- do.call(rbind, clustree_rows)

  if (!identical(metadata$cluster_column, pairwise_summary$cluster_column)) {
    stop("Pairwise stability summary order mismatch.", call. = FALSE)
  }
  if (!identical(metadata$cluster_column, clustree_summary$cluster_column)) {
    stop("Clustree stability summary order mismatch.", call. = FALSE)
  }
  size_rows <- lapply(metadata$cluster_column, function(column) {
    size_summary(labels[[column]])
  })
  stability_summary <- cbind(
    metadata,
    do.call(rbind, size_rows),
    pairwise_summary[setdiff(colnames(pairwise_summary), "cluster_column")],
    clustree_summary[setdiff(colnames(clustree_summary), "cluster_column")]
  )
  rownames(stability_summary) <- NULL

  # ANALYSIS_OK[R026]: local TSV writer is called by the grid-table implementation below.
  write_tsv <- function(data, path) {
    if (is.null(path)) {
      return(invisible(NULL))
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      data,
      path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = ""
    )
  }
  write_tsv(grid_summary, summary_path)
  write_tsv(stability_summary, stability_path)
  write_tsv(pairwise, pairwise_path)
  invisible(list(summary = stability_summary, pairwise = pairwise))
}

DEFAULT_CLUSTER_DIMS <- c(20L, 30L, 50L)
DEFAULT_CLUSTER_RESOLUTIONS <- c(0.3, 0.5, 0.8)
DEFAULT_CLUSTER_REFERENCE_COLUMN <- "cluster_pflog_filter_cc_dims50_res0.3"
SMALL_CLUSTER_CELL_THRESHOLD <- 50L
ENTROPY_LOG_BASE <- 2
CLUSTREE_STABLE_CHILD_FRACTION <- 0.8
