# Cluster summaries and sample-level inference.

#' Compute pooled cluster abundance enrichment by condition
#'
#' Fisher tests use pooled cell counts, with Holm adjustment across clusters.
#' CLR effects use a 0.5 pseudocount. These are descriptive cell-level summaries,
#' not condition-level inference with Mouse x Condition replication.
#' @param sobj Seurat object containing cluster and condition metadata.
#' @param cluster_col Metadata column containing cluster labels.
#' @param condition_col Metadata column containing condition labels.
#' @param control_label,estim_label Condition labels.
#' @return A data frame of counts, CLR effects, p-values, and directions.
#' @export
compute_cluster_abundance <- function(
  sobj,
  cluster_col,
  condition_col = CONDITION_COL,
  control_label = CTRL_LABEL,
  estim_label = ESTIM_LABEL
) {
  meta <- sobj[[]]
  cluster <- as.character(meta[[cluster_col]])
  condition <- as.character(meta[[condition_col]])
  if (is.null(meta[[cluster_col]]) || anyNA(cluster) || any(!nzchar(cluster))) {
    stop("Missing or empty cluster labels: ", cluster_col)
  }
  expected_conditions <- c(control_label, estim_label)
  if (anyNA(condition) || !setequal(unique(condition), expected_conditions)) {
    stop("Condition labels must contain exactly the configured comparison.")
  }
  counts_mat <- table(
    cluster = factor(cluster, levels = .sort_cluster_labels(cluster)),
    condition = factor(condition, levels = expected_conditions)
  ) |>
    as.matrix()
  storage.mode(counts_mat) <- "integer"
  condition_totals <- colSums(counts_mat)

  clr_mat <- apply(counts_mat + 0.5, 2L, function(x) {
    log_x <- log(x)
    log_x - mean(log_x)
  })
  if (is.null(dim(clr_mat))) {
    clr_mat <- matrix(
      clr_mat,
      nrow = nrow(counts_mat),
      dimnames = dimnames(counts_mat)
    )
  }
  fisher_p <- vapply(
    seq_len(nrow(counts_mat)),
    function(row_idx) {
      estim_count <- counts_mat[row_idx, estim_label]
      control_count <- counts_mat[row_idx, control_label]
      fisher_mat <- matrix(
        c(
          estim_count,
          condition_totals[[estim_label]] - estim_count,
          control_count,
          condition_totals[[control_label]] - control_count
        ),
        nrow = 2L,
        byrow = TRUE
      )
      stats::fisher.test(fisher_mat)$p.value
    },
    numeric(1)
  )
  clr_diff <- clr_mat[, estim_label] - clr_mat[, control_label]
  log2_enrichment <- as.numeric(clr_diff) / log(2)
  padj <- stats::p.adjust(fisher_p, method = "holm")
  direction <- rep("Not significant", length(log2_enrichment))
  direction[padj < 0.05 & log2_enrichment > 0] <- "Enriched in E-Stim"
  direction[padj < 0.05 & log2_enrichment < 0] <- "Depleted in E-Stim"

  data.frame(
    cluster = rownames(counts_mat),
    control_count = as.integer(counts_mat[, control_label]),
    estim_count = as.integer(counts_mat[, estim_label]),
    control_clr = as.numeric(clr_mat[, control_label]),
    estim_clr = as.numeric(clr_mat[, estim_label]),
    clr_diff = as.numeric(clr_diff),
    log2_enrichment = log2_enrichment,
    p_value = fisher_p,
    padj = padj,
    direction = direction,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Compute cluster proportions for each Mouse x Condition sample
#'
#' Includes zero-count clusters within observed samples, but does not invent
#' samples for missing Mouse x Condition combinations.
#' @param sobj Seurat object containing cluster, mouse, and condition metadata.
#' @param cluster_col,mouse_col,condition_col Metadata column names.
#' @param control_label,estim_label Condition labels.
#' @return A data frame with one row per observed sample and cluster.
#' @export
compute_sample_cluster_proportions <- function(
  sobj,
  cluster_col,
  mouse_col = "Mouse",
  condition_col = CONDITION_COL,
  control_label = CTRL_LABEL,
  estim_label = ESTIM_LABEL
) {
  meta <- sobj[[]]
  columns <- c(cluster_col, mouse_col, condition_col)
  if (!all(columns %in% names(meta)) || anyNA(meta[columns])) {
    stop(
      "Cluster, mouse, and condition metadata must be present and non-missing."
    )
  }
  cells <- tibble::tibble(
    mouse = as.character(meta[[mouse_col]]),
    condition = as.character(meta[[condition_col]]),
    cluster = as.character(meta[[cluster_col]])
  )
  if (any(as.matrix(cells) == "")) {
    stop("Cell metadata contains empty labels.")
  }
  expected_conditions <- c(control_label, estim_label)
  if (!setequal(unique(cells$condition), expected_conditions)) {
    stop("Condition labels must contain exactly the configured comparison.")
  }
  cluster_levels <- .sort_cluster_labels(cells$cluster)
  mouse_levels <- .sort_cluster_labels(cells$mouse)
  samples <- cells |> dplyr::count(mouse, condition, name = "sample_total")
  paired_mice <- samples |>
    dplyr::count(mouse) |>
    dplyr::filter(n == 2L) |>
    dplyr::pull(mouse)
  samples <- samples |>
    dplyr::mutate(
      mouse_role = dplyr::case_when(
        mouse %in% paired_mice ~ "paired",
        condition == control_label ~ "control_only",
        TRUE ~ "estim_only"
      )
    )
  cells |>
    dplyr::count(mouse, condition, cluster, name = "cluster_n") |>
    tidyr::complete(
      tidyr::nesting(mouse, condition),
      cluster = cluster_levels,
      fill = list(cluster_n = 0L)
    ) |>
    dplyr::left_join(
      samples,
      by = c("mouse", "condition"),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(proportion = cluster_n / sample_total) |>
    dplyr::arrange(
      match(condition, expected_conditions),
      match(mouse, mouse_levels),
      match(cluster, cluster_levels)
    ) |>
    dplyr::select(
      mouse,
      condition,
      mouse_role,
      cluster,
      cluster_n,
      sample_total,
      proportion
    ) |>
    as.data.frame()
}

# All sign assignments for the paired contrasts.
.sign_vectors <- function(k) {
  as.matrix(expand.grid(rep(list(c(-1, 1)), k)))
}

#' Test sample-level cluster-proportion shifts by exact randomization
#'
#' Enumerates signs of paired-mouse contrasts on stabilized logits:
#' `qlogis((cluster_n + 0.5) / (sample_total + 1))`. A sensitivity analysis adds
#' one singleton contrast when exactly one control-only and one E-Stim-only
#' mouse are present. Both tests are two-sided, with BH adjustment across clusters.
#' @param sample_props Output of [compute_sample_cluster_proportions()].
#' @param control_label,estim_label Condition labels.
#' @return A data frame with effects, p-values, q-values, and sensitivity status.
#' @export
test_cluster_proportion_randomization <- function(
  sample_props,
  control_label = CTRL_LABEL,
  estim_label = ESTIM_LABEL
) {
  sample_props <- as.data.frame(sample_props)
  if (anyNA(sample_props) || nrow(sample_props) == 0L) {
    stop("Sample proportions must be non-empty and non-missing.")
  }
  if (
    any(
      sample_props$sample_total <= 0 |
        sample_props$cluster_n < 0 |
        sample_props$cluster_n > sample_props$sample_total
    )
  ) {
    stop("Invalid cluster counts or sample totals.")
  }
  sample_props$logit_effect <- stats::qlogis(
    (sample_props$cluster_n + 0.5) / (sample_props$sample_total + 1)
  )
  cluster_levels <- .sort_cluster_labels(sample_props$cluster)
  paired_mice <- sort(
    unique(sample_props$mouse[sample_props$mouse_role == "paired"]),
    method = "radix"
  )
  if (length(paired_mice) == 0L) {
    stop("Paired randomization requires paired mice.")
  }
  control_only_mice <- sort(
    unique(sample_props$mouse[sample_props$mouse_role == "control_only"]),
    method = "radix"
  )
  estim_only_mice <- sort(
    unique(sample_props$mouse[sample_props$mouse_role == "estim_only"]),
    method = "radix"
  )
  singleton_runs <- length(control_only_mice) == 1L &&
    length(estim_only_mice) == 1L
  singleton_status <- sprintf(
    "skipped: expected 1 control-only and 1 E-Stim-only mouse; found %d and %d",
    length(control_only_mice),
    length(estim_only_mice)
  )

  result_rows <- lapply(cluster_levels, function(cluster_label) {
    cluster_data <- sample_props[
      sample_props$cluster == cluster_label,
      ,
      drop = FALSE
    ]
    get_logit <- function(mouse_id, condition_label) {
      row_idx <- cluster_data$mouse == mouse_id &
        cluster_data$condition == condition_label
      if (sum(row_idx) != 1L) {
        stop(
          "Expected one row for cluster ",
          cluster_label,
          ", mouse ",
          mouse_id,
          ", condition ",
          condition_label
        )
      }
      cluster_data$logit_effect[row_idx]
    }
    diffs <- vapply(
      paired_mice,
      function(mouse_id) {
        get_logit(mouse_id, estim_label) - get_logit(mouse_id, control_label)
      },
      numeric(1)
    )
    effect_paired <- mean(diffs)
    paired_signs <- .sign_vectors(length(diffs))
    paired_null <- as.numeric(paired_signs %*% diffs) / length(diffs)
    p_value_paired <- mean(abs(paired_null) >= abs(effect_paired) - 1e-8)
    control_prop_idx <- cluster_data$mouse %in%
      paired_mice &
      cluster_data$condition == control_label
    estim_prop_idx <- cluster_data$mouse %in%
      paired_mice &
      cluster_data$condition == estim_label

    if (singleton_runs) {
      singleton_diff <- get_logit(estim_only_mice[[1L]], estim_label) -
        get_logit(control_only_mice[[1L]], control_label)
      contrasts <- c(diffs, singleton_diff)
      signs <- .sign_vectors(length(contrasts))
      effect_paired_singleton <- mean(contrasts)
      null <- as.numeric(signs %*% contrasts) / length(contrasts)
      p_value_paired_singleton <- mean(
        abs(null) >= abs(effect_paired_singleton) - 1e-8
      )
      n_perm_paired_singleton <- nrow(signs)
      paired_singleton_status <- "ok"
    } else {
      effect_paired_singleton <- NA_real_
      p_value_paired_singleton <- NA_real_
      n_perm_paired_singleton <- NA_integer_
      paired_singleton_status <- singleton_status
    }
    direction_paired <- if (effect_paired > 0) {
      "Higher in E-Stim"
    } else if (effect_paired < 0) {
      "Lower in E-Stim"
    } else {
      "No difference"
    }
    data.frame(
      cluster = cluster_label,
      mean_proportion_control_paired = mean(cluster_data$proportion[
        control_prop_idx
      ]),
      mean_proportion_estim_paired = mean(cluster_data$proportion[
        estim_prop_idx
      ]),
      effect_paired = effect_paired,
      p_value_paired = p_value_paired,
      q_value_paired = NA_real_,
      n_perm_paired = nrow(paired_signs),
      direction_paired = direction_paired,
      effect_paired_singleton = effect_paired_singleton,
      p_value_paired_singleton = p_value_paired_singleton,
      q_value_paired_singleton = NA_real_,
      n_perm_paired_singleton = n_perm_paired_singleton,
      paired_singleton_status = paired_singleton_status,
      stringsAsFactors = FALSE
    )
  })
  results <- do.call(rbind, result_rows)
  results$q_value_paired <- stats::p.adjust(
    results$p_value_paired,
    method = "BH"
  )
  results$q_value_paired_singleton <- stats::p.adjust(
    results$p_value_paired_singleton,
    method = "BH"
  )
  rownames(results) <- NULL
  results
}

#' Compute cluster-level cell-type marker module scores
#'
#' Scores cells with [Seurat::AddModuleScore()] and averages scores within clusters.
#' @param sobj Seurat object containing expression and cluster metadata.
#' @param cluster_col Cluster metadata column.
#' @param marker_genes Named list of marker-gene vectors.
#' @param assay,slot Assay and expression layer used by AddModuleScore.
#' @param seed Random seed for scoring.
#' @return A numeric matrix with modules in rows and sorted clusters in columns.
#' @export
compute_cluster_module_scores <- function(
  sobj,
  cluster_col,
  marker_genes,
  assay = "RNA",
  slot = "data",
  seed = SEED
) {
  cluster <- as.character(sobj[[cluster_col, drop = TRUE]])
  missing_genes <- setdiff(
    unlist(marker_genes, use.names = FALSE),
    rownames(sobj[[assay]])
  )
  if (length(missing_genes) > 0L) {
    stop(
      "Marker genes missing from the assay: ",
      paste(missing_genes, collapse = ", ")
    )
  }
  if (anyNA(cluster) || any(!nzchar(cluster))) {
    stop("Missing or empty cluster labels.")
  }
  cluster_levels <- .sort_cluster_labels(cluster)
  module_prefix <- "celltype_module_score"
  module_score_cols <- paste0(module_prefix, seq_along(marker_genes))
  sobj@meta.data[intersect(module_score_cols, colnames(sobj@meta.data))] <- NULL
  sobj <- Seurat::AddModuleScore(
    object = sobj,
    features = marker_genes,
    assay = assay,
    name = module_prefix,
    seed = seed,
    search = FALSE,
    slot = slot
  )
  module_scores <- sobj@meta.data[module_score_cols]
  colnames(module_scores) <- names(marker_genes)
  cluster_marker_scores <- stats::aggregate(
    module_scores,
    by = list(cluster = cluster),
    FUN = mean
  )
  cluster_marker_scores <- cluster_marker_scores[
    match(cluster_levels, cluster_marker_scores$cluster),
    ,
    drop = FALSE
  ]
  score_matrix <- t(as.matrix(cluster_marker_scores[names(marker_genes)]))
  storage.mode(score_matrix) <- "numeric"
  rownames(score_matrix) <- names(marker_genes)
  colnames(score_matrix) <- cluster_levels
  score_matrix
}

#' Compute sample-aware cluster p27 enrichment z-scores
#'
#' Permutes cluster labels within each Mouse x Condition sample, preserving
#' sample-specific expression and cluster sizes. Reports NA for a cluster with
#' zero or non-finite null standard deviation. Restores the caller's RNG state.
#' @param sobj Seurat object containing expression and metadata.
#' @param cluster_col Cluster metadata column.
#' @param gene Gene to summarize.
#' @param layer,assay Expression layer and assay.
#' @param mouse_col,condition_col Sample metadata columns.
#' @param n_perm Positive integer number of within-sample permutations.
#' @param seed Random seed for permutations.
#' @return A data frame with cluster sizes, observed/null means, null SDs, and z-scores.
#' @export
compute_cluster_p27_enrichment <- function(
  sobj,
  cluster_col,
  gene = "Cdkn1b",
  layer = "pflog",
  assay = "RNA",
  mouse_col = "Mouse",
  condition_col = CONDITION_COL,
  n_perm = 2000L,
  seed = SEED
) {
  meta <- sobj[[]]
  columns <- c(cluster_col, mouse_col, condition_col)
  if (!all(columns %in% names(meta)) || anyNA(meta[columns])) {
    stop(
      "Cluster, mouse, and condition metadata must be present and non-missing."
    )
  }
  if (any(as.matrix(meta[columns]) == "")) {
    stop("Cell metadata contains empty labels.")
  }
  stopifnot(
    length(n_perm) == 1L,
    is.finite(n_perm),
    n_perm >= 1L,
    n_perm == as.integer(n_perm)
  )
  cluster <- as.character(meta[[cluster_col]])
  mouse <- as.character(meta[[mouse_col]])
  condition <- as.character(meta[[condition_col]])
  expression_matrix <- SeuratObject::LayerData(sobj[[assay]], layer = layer)
  if (!identical(colnames(expression_matrix), rownames(meta))) {
    stop(
      "Expression columns and metadata rows must refer to the same cells in order."
    )
  }
  expr <- as.numeric(expression_matrix[gene, , drop = TRUE])
  if (any(!is.finite(expr))) {
    stop("Expression values must be finite.")
  }
  cluster_levels <- .sort_cluster_labels(cluster)
  cluster_factor <- factor(cluster, levels = cluster_levels)
  observed <- as.numeric(tapply(expr, cluster_factor, mean))
  n_cells <- as.integer(table(cluster_factor))
  sample_id <- paste(mouse, condition, sep = "__")
  sample_indices <- split(seq_along(sample_id), sample_id)
  null_mat <- matrix(
    NA_real_,
    nrow = n_perm,
    ncol = length(cluster_levels),
    dimnames = list(NULL, cluster_levels)
  )

  withr::local_seed(seed)
  for (perm_idx in seq_len(n_perm)) {
    permuted_cluster <- cluster
    for (idx in sample_indices) {
      permuted_cluster[idx] <- sample(
        cluster[idx],
        length(idx),
        replace = FALSE
      )
    }
    null_mat[perm_idx, ] <- as.numeric(tapply(
      expr,
      factor(permuted_cluster, levels = cluster_levels),
      mean
    ))
  }
  null_mean <- colMeans(null_mat)
  null_sd <- apply(null_mat, 2L, stats::sd)
  z_score <- (observed - null_mean) / null_sd
  z_score[!is.finite(null_sd) | null_sd == 0] <- NA_real_
  data.frame(
    cluster = cluster_levels,
    n_cells = n_cells,
    observed_mean = observed,
    null_mean = null_mean,
    null_sd = null_sd,
    z_score = z_score,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.sort_cluster_labels <- function(x) {
  x <- unique(as.character(x))
  if (all(grepl("^-?[0-9]+$", x))) {
    return(as.character(sort(as.integer(x))))
  }
  sort(x, method = "radix")
}
