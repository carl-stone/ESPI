# Expected design cardinality and significance threshold for two-condition abundance tests.
EXPECTED_CONDITION_COUNT <- 2L
FDR_THRESHOLD <- 0.05

#' Compute pooled cluster abundance enrichment by condition.
#'
#' Uses pooled cell-level cluster x condition counts. Fisher tests use raw counts;
#' CLR effects use counts plus a 0.5 pseudocount and are descriptive relative to
#' the Mouse x Condition statistical unit used by differential analyses.
#'
#' @param sobj Seurat object containing cluster and condition metadata.
#' @param cluster_col Metadata column containing cluster labels.
#' @param condition_col Metadata column containing condition labels.
#' @param control_label Control condition label.
#' @param estim_label E-Stim condition label.
#'
#' @return Data frame with per-cluster counts, CLR log2 enrichment, Fisher
#'   p-values, Holm-adjusted p-values, and enrichment direction.
#' @export
compute_cluster_abundance <- function(
  sobj,
  cluster_col,
  condition_col = CONDITION_COL,
  control_label = CTRL_LABEL,
  estim_label = ESTIM_LABEL
) {
  meta <- sobj@meta.data
  required_cols <- c(cluster_col, condition_col)
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing required metadata column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  cluster <- as.character(meta[[cluster_col]])
  condition <- as.character(meta[[condition_col]])
  if (length(cluster) != nrow(meta) || length(condition) != nrow(meta)) {
    stop(
      "Cluster and condition metadata are not aligned to cells.",
      call. = FALSE
    )
  }
  if (anyNA(cluster) || any(!nzchar(cluster))) {
    stop(
      "Cluster metadata contains missing or empty labels: ",
      cluster_col,
      call. = FALSE
    )
  }
  if (anyNA(condition) || any(!nzchar(condition))) {
    stop(
      "Condition metadata contains missing or empty labels: ",
      condition_col,
      call. = FALSE
    )
  }

  expected_conditions <- c(control_label, estim_label)
  unexpected_conditions <- setdiff(unique(condition), expected_conditions)
  if (length(unexpected_conditions) > 0L) {
    stop(
      "Unexpected condition label(s) in ",
      condition_col,
      ": ",
      paste(unexpected_conditions, collapse = ", "),
      call. = FALSE
    )
  }
  missing_conditions <- setdiff(expected_conditions, unique(condition))
  if (length(missing_conditions) > 0L) {
    stop(
      "Missing expected condition label(s) in ",
      condition_col,
      ": ",
      paste(missing_conditions, collapse = ", "),
      call. = FALSE
    )
  }

  cluster_levels <- unique(cluster)
  if (all(grepl("^-?[0-9]+$", cluster_levels))) {
    cluster_levels <- as.character(sort(as.integer(cluster_levels)))
  } else {
    cluster_levels <- sort(cluster_levels, method = "radix")
  }
  counts_mat <- table(
    cluster = factor(cluster, levels = cluster_levels),
    condition = factor(condition, levels = expected_conditions)
  )
  counts_mat <- as.matrix(counts_mat)
  storage.mode(counts_mat) <- "integer"
  if (
    nrow(counts_mat) == 0L ||
      ncol(counts_mat) != EXPECTED_CONDITION_COUNT ||
      sum(counts_mat) == 0L
  ) {
    stop("Cluster x condition count table is empty.", call. = FALSE)
  }
  if (!identical(colnames(counts_mat), expected_conditions)) {
    stop(
      "Condition count columns are not aligned to expected labels.",
      call. = FALSE
    )
  }

  condition_totals <- colSums(counts_mat)
  if (any(condition_totals <= 0L)) {
    stop(
      "Each expected condition must contain at least one cell.",
      call. = FALSE
    )
  }

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
  if (
    !identical(rownames(clr_mat), rownames(counts_mat)) ||
      !identical(colnames(clr_mat), colnames(counts_mat))
  ) {
    stop("CLR matrix names do not align to count matrix names.", call. = FALSE)
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
    numeric(1L)
  )

  clr_diff <- clr_mat[, estim_label] - clr_mat[, control_label]
  log2_enrichment <- as.numeric(clr_diff) / log(2)
  if (any(!is.finite(log2_enrichment))) {
    stop("Non-finite CLR log2 enrichment value(s) produced.", call. = FALSE)
  }

  padj <- stats::p.adjust(fisher_p, method = "holm")

  direction <- rep("Not significant", length(log2_enrichment))
  direction[padj < FDR_THRESHOLD & log2_enrichment > 0] <- "Enriched in E-Stim"
  direction[padj < FDR_THRESHOLD & log2_enrichment < 0] <- "Depleted in E-Stim"

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

#' Compute sample-level cluster proportions by Mouse x Condition.
#'
#' Emits one row for every observed Mouse x Condition sample and cluster,
#' including zero-count clusters within observed samples.
#'
#' @param sobj Seurat object containing cluster, mouse, and condition metadata.
#' @param cluster_col Metadata column containing cluster labels.
#' @param mouse_col Metadata column containing mouse labels.
#' @param condition_col Metadata column containing condition labels.
#' @param control_label Control condition label.
#' @param estim_label E-Stim condition label.
#'
#' @return Data frame with one row per observed Mouse x Condition sample and
#'   cluster.
#' @export
# ANALYSIS_OK[smuggled-default]: intentional package API default for Mouse metadata.
compute_sample_cluster_proportions <- function(
  sobj,
  cluster_col,
  mouse_col = "Mouse",
  condition_col = CONDITION_COL,
  control_label = CTRL_LABEL,
  estim_label = ESTIM_LABEL
) {
  meta <- sobj@meta.data
  required_cols <- c(cluster_col, mouse_col, condition_col)
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing required metadata column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  cluster <- as.character(meta[[cluster_col]])
  mouse <- as.character(meta[[mouse_col]])
  condition <- as.character(meta[[condition_col]])
  if (
    length(cluster) != nrow(meta) ||
      length(mouse) != nrow(meta) ||
      length(condition) != nrow(meta)
  ) {
    stop(
      "Cluster, mouse, and condition metadata are not aligned to cells.",
      call. = FALSE
    )
  }
  if (anyNA(cluster) || any(!nzchar(cluster))) {
    stop(
      "Cluster metadata contains missing or empty labels: ",
      cluster_col,
      call. = FALSE
    )
  }
  if (anyNA(mouse) || any(!nzchar(mouse))) {
    stop(
      "Mouse metadata contains missing or empty labels: ",
      mouse_col,
      call. = FALSE
    )
  }
  if (anyNA(condition) || any(!nzchar(condition))) {
    stop(
      "Condition metadata contains missing or empty labels: ",
      condition_col,
      call. = FALSE
    )
  }

  expected_conditions <- c(control_label, estim_label)
  unexpected_conditions <- setdiff(unique(condition), expected_conditions)
  if (length(unexpected_conditions) > 0L) {
    stop(
      "Unexpected condition label(s) in ",
      condition_col,
      ": ",
      paste(unexpected_conditions, collapse = ", "),
      call. = FALSE
    )
  }
  missing_conditions <- setdiff(expected_conditions, unique(condition))
  if (length(missing_conditions) > 0L) {
    stop(
      "Missing expected condition label(s) in ",
      condition_col,
      ": ",
      paste(missing_conditions, collapse = ", "),
      call. = FALSE
    )
  }

  cluster_levels <- unique(cluster)
  if (all(grepl("^-?[0-9]+$", cluster_levels))) {
    cluster_levels <- as.character(sort(as.integer(cluster_levels)))
  } else {
    cluster_levels <- sort(cluster_levels, method = "radix")
  }

  mouse_levels <- unique(mouse)
  if (all(grepl("^-?[0-9]+$", mouse_levels))) {
    mouse_levels <- as.character(sort(as.integer(mouse_levels)))
  } else {
    mouse_levels <- sort(mouse_levels, method = "radix")
  }
  sample_tab <- table(
    mouse = factor(mouse, levels = mouse_levels),
    condition = factor(condition, levels = expected_conditions)
  )
  sample_tab <- as.matrix(sample_tab)
  storage.mode(sample_tab) <- "integer"
  if (
    nrow(sample_tab) == 0L ||
      ncol(sample_tab) != EXPECTED_CONDITION_COUNT ||
      sum(sample_tab) == 0L
  ) {
    stop("Mouse x condition count table is empty.", call. = FALSE)
  }
  if (!identical(colnames(sample_tab), expected_conditions)) {
    stop(
      "Condition count columns are not aligned to expected labels.",
      call. = FALSE
    )
  }

  present_mat <- sample_tab > 0L
  paired_mouse <- present_mat[, control_label] & present_mat[, estim_label]
  if (!any(paired_mouse)) {
    stop(
      "No mice have both conditions; paired randomization is undefined.",
      call. = FALSE
    )
  }

  mouse_roles <- rep(NA_character_, length(mouse_levels))
  names(mouse_roles) <- mouse_levels
  mouse_roles[paired_mouse] <- "paired"
  mouse_roles[
    !paired_mouse &
      present_mat[, estim_label] &
      !present_mat[, control_label]
  ] <- "estim_only"
  mouse_roles[
    !paired_mouse &
      present_mat[, control_label] &
      !present_mat[, estim_label]
  ] <- "control_only"
  if (anyNA(mouse_roles)) {
    stop("Unable to assign mouse roles from observed samples.", call. = FALSE)
  }

  counts_array <- table(
    mouse = factor(mouse, levels = mouse_levels),
    condition = factor(condition, levels = expected_conditions),
    cluster = factor(cluster, levels = cluster_levels)
  )
  storage.mode(counts_array) <- "integer"

  output_rows <- vector(
    "list",
    sum(present_mat) * length(cluster_levels)
  )
  row_idx <- 1L
  for (condition_label in expected_conditions) {
    for (mouse_id in mouse_levels) {
      sample_total <- as.integer(sample_tab[mouse_id, condition_label])
      if (sample_total == 0L) {
        next
      }
      for (cluster_label in cluster_levels) {
        cluster_n <- as.integer(
          counts_array[mouse_id, condition_label, cluster_label]
        )
        output_rows[[row_idx]] <- data.frame(
          mouse = mouse_id,
          condition = condition_label,
          mouse_role = mouse_roles[[mouse_id]],
          cluster = cluster_label,
          cluster_n = cluster_n,
          sample_total = sample_total,
          proportion = cluster_n / sample_total,
          stringsAsFactors = FALSE
        )
        row_idx <- row_idx + 1L
      }
    }
  }

  sample_props <- do.call(rbind, output_rows)
  if (any(sample_props$sample_total == 0L)) {
    stop("Observed sample has zero cells.", call. = FALSE)
  }
  row.names(sample_props) <- NULL
  sample_props
}

# ANALYSIS_OK[R026]: package helper is loaded by devtools::load_all and called by the same-file permutation test.
.sign_vectors <- function(k) {
  if (
    length(k) != 1L ||
      is.na(k) ||
      k < 1L ||
      k != as.integer(k)
  ) {
    stop("k must be a positive integer.", call. = FALSE)
  }
  sign_grid <- expand.grid(rep(list(c(-1, 1)), as.integer(k)))
  as.matrix(sign_grid)
}

#' Test sample-level cluster-proportion shifts by randomization.
#'
#' Uses exact sign enumeration over paired Mouse x Condition contrasts, with an
#' optional paired-plus-singleton sensitivity when exactly one E-Stim-only and
#' one control-only mouse are present.
#'
#' @param sample_props Data frame returned by
#'   `compute_sample_cluster_proportions()`.
#' @param control_label Control condition label.
#' @param estim_label E-Stim condition label.
#'
#' @return Data frame with per-cluster effect estimates, exact randomization
#'   p-values, BH-adjusted q-values, and singleton sensitivity status.
#' @export
# ANALYSIS_OK[R026]: package export is loaded by devtools::load_all and invoked by executable analysis scripts.
test_cluster_proportion_randomization <- function(
  sample_props,
  control_label = CTRL_LABEL,
  estim_label = ESTIM_LABEL
) {
  if (!is.data.frame(sample_props)) {
    stop("sample_props must be a data frame.", call. = FALSE)
  }
  required_cols <- c(
    "mouse",
    "condition",
    "mouse_role",
    "cluster",
    "cluster_n",
    "sample_total",
    "proportion"
  )
  missing_cols <- setdiff(required_cols, colnames(sample_props))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing sample proportion column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(sample_props) == 0L) {
    stop("Sample proportion table is empty.", call. = FALSE)
  }

  sample_props$mouse <- as.character(sample_props$mouse)
  sample_props$condition <- as.character(sample_props$condition)
  sample_props$mouse_role <- as.character(sample_props$mouse_role)
  sample_props$cluster <- as.character(sample_props$cluster)
  if (anyNA(sample_props$mouse) || any(!nzchar(sample_props$mouse))) {
    stop(
      "Sample proportion table contains missing mouse labels.",
      call. = FALSE
    )
  }
  if (anyNA(sample_props$condition) || any(!nzchar(sample_props$condition))) {
    stop(
      "Sample proportion table contains missing condition labels.",
      call. = FALSE
    )
  }
  if (anyNA(sample_props$cluster) || any(!nzchar(sample_props$cluster))) {
    stop(
      "Sample proportion table contains missing cluster labels.",
      call. = FALSE
    )
  }

  expected_conditions <- c(control_label, estim_label)
  unexpected_conditions <- setdiff(
    unique(sample_props$condition),
    expected_conditions
  )
  if (length(unexpected_conditions) > 0L) {
    stop(
      "Unexpected condition label(s): ",
      paste(unexpected_conditions, collapse = ", "),
      call. = FALSE
    )
  }
  missing_conditions <- setdiff(
    expected_conditions,
    unique(sample_props$condition)
  )
  if (length(missing_conditions) > 0L) {
    stop(
      "Missing expected condition label(s): ",
      paste(missing_conditions, collapse = ", "),
      call. = FALSE
    )
  }

  expected_roles <- c("paired", "estim_only", "control_only")
  unexpected_roles <- setdiff(unique(sample_props$mouse_role), expected_roles)
  if (length(unexpected_roles) > 0L) {
    stop(
      "Unexpected mouse role(s): ",
      paste(unexpected_roles, collapse = ", "),
      call. = FALSE
    )
  }

  if (
    anyNA(sample_props$cluster_n) ||
      anyNA(sample_props$sample_total) ||
      anyNA(sample_props$proportion)
  ) {
    stop(
      "Sample proportion table contains missing count or proportion values.",
      call. = FALSE
    )
  }
  if (
    any(sample_props$sample_total <= 0L) ||
      any(sample_props$cluster_n < 0L) ||
      any(sample_props$cluster_n > sample_props$sample_total)
  ) {
    stop(
      "Sample proportion table contains invalid cluster or sample counts.",
      call. = FALSE
    )
  }
  if (
    any(!is.finite(sample_props$proportion)) ||
      any(sample_props$proportion < 0) ||
      any(sample_props$proportion > 1)
  ) {
    stop(
      "Sample proportion values must be finite and in [0, 1].",
      call. = FALSE
    )
  }

  sample_props$logit_effect <- stats::qlogis(
    (sample_props$cluster_n + 0.5) / (sample_props$sample_total + 1)
  )
  if (any(!is.finite(sample_props$logit_effect))) {
    stop("Non-finite stabilized logit value(s) produced.", call. = FALSE)
  }

  cluster_levels <- unique(sample_props$cluster)
  if (all(grepl("^-?[0-9]+$", cluster_levels))) {
    cluster_levels <- as.character(sort(as.integer(cluster_levels)))
  } else {
    cluster_levels <- sort(cluster_levels, method = "radix")
  }
  paired_mice <- sort(
    unique(sample_props$mouse[sample_props$mouse_role == "paired"]),
    method = "radix"
  )
  if (length(paired_mice) == 0L) {
    stop(
      "No mice have both conditions; paired randomization is undefined.",
      call. = FALSE
    )
  }
  control_only_mice <- sort(
    unique(sample_props$mouse[sample_props$mouse_role == "control_only"]),
    method = "radix"
  )
  estim_only_mice <- sort(
    unique(sample_props$mouse[sample_props$mouse_role == "estim_only"]),
    method = "radix"
  )
  singleton_status <- sprintf(
    "skipped: expected 1 control-only and 1 E-Stim-only mouse; found %d and %d",
    length(control_only_mice),
    length(estim_only_mice)
  )
  singleton_runs <- length(control_only_mice) == 1L &&
    length(estim_only_mice) == 1L

  result_rows <- vector("list", length(cluster_levels))
  for (cluster_idx in seq_along(cluster_levels)) {
    cluster_label <- cluster_levels[[cluster_idx]]
    cluster_data <- sample_props[
      sample_props$cluster == cluster_label,
      ,
      drop = FALSE
    ]
    # ANALYSIS_OK[R026]: package helper is loaded by devtools::load_all and called by same-file logit calculations.
    get_logit <- function(mouse_id, condition_label) {
      row_idx <- cluster_data$mouse == mouse_id &
        cluster_data$condition == condition_label
      if (sum(row_idx) != 1L) {
        stop(
          "Expected exactly one row for cluster ",
          cluster_label,
          ", mouse ",
          mouse_id,
          ", condition ",
          condition_label,
          ".",
          call. = FALSE
        )
      }
      cluster_data$logit_effect[row_idx]
    }

    diffs <- vapply(
      paired_mice,
      function(mouse_id) {
        get_logit(mouse_id, estim_label) - get_logit(mouse_id, control_label)
      },
      numeric(1L)
    )
    effect_paired <- mean(diffs)
    paired_signs <- .sign_vectors(length(diffs))
    paired_null <- as.numeric(paired_signs %*% diffs) / length(diffs)
    p_value_paired <- mean(
      abs(paired_null) >= abs(effect_paired) - 1e-8
    )

    control_prop_idx <- cluster_data$mouse %in%
      paired_mice &
      cluster_data$condition == control_label
    estim_prop_idx <- cluster_data$mouse %in%
      paired_mice &
      cluster_data$condition == estim_label
    if (
      sum(control_prop_idx) != length(paired_mice) ||
        sum(estim_prop_idx) != length(paired_mice)
    ) {
      stop(
        "Paired sample rows are incomplete for cluster ",
        cluster_label,
        ".",
        call. = FALSE
      )
    }
    mean_proportion_control_paired <- mean(
      cluster_data$proportion[control_prop_idx]
    )
    mean_proportion_estim_paired <- mean(
      cluster_data$proportion[estim_prop_idx]
    )

    if (singleton_runs) {
      singleton_diff <- get_logit(estim_only_mice[[1L]], estim_label) -
        get_logit(control_only_mice[[1L]], control_label)
      paired_singleton_contrasts <- c(diffs, singleton_diff)
      paired_singleton_signs <- .sign_vectors(
        length(paired_singleton_contrasts)
      )
      effect_paired_singleton <- mean(paired_singleton_contrasts)
      paired_singleton_null <- as.numeric(
        paired_singleton_signs %*% paired_singleton_contrasts
      ) /
        length(paired_singleton_contrasts)
      p_value_paired_singleton <- mean(
        abs(paired_singleton_null) >= abs(effect_paired_singleton) - 1e-8
      )
      n_perm_paired_singleton <- nrow(paired_singleton_signs)
      paired_singleton_status <- "ok"
    } else {
      effect_paired_singleton <- NA_real_
      p_value_paired_singleton <- NA_real_
      n_perm_paired_singleton <- NA_integer_
      paired_singleton_status <- singleton_status
    }

    direction_paired <- "No difference"
    if (effect_paired > 0) {
      direction_paired <- "Higher in E-Stim"
    } else if (effect_paired < 0) {
      direction_paired <- "Lower in E-Stim"
    }

    result_rows[[cluster_idx]] <- data.frame(
      cluster = cluster_label,
      mean_proportion_control_paired = mean_proportion_control_paired,
      mean_proportion_estim_paired = mean_proportion_estim_paired,
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
  }

  results <- do.call(rbind, result_rows)
  # ANALYSIS_OK[R026]: package helper is loaded by devtools::load_all and called by same-file q-value adjustment.
  adjust_non_missing <- function(p_value) {
    q_value <- rep(NA_real_, length(p_value))
    observed <- !is.na(p_value)
    q_value[observed] <- stats::p.adjust(p_value[observed], method = "BH")
    q_value
  }
  results$q_value_paired <- adjust_non_missing(results$p_value_paired)
  results$q_value_paired_singleton <- adjust_non_missing(
    results$p_value_paired_singleton
  )
  # ANALYSIS_OK[result-schema]: explicit column selection only reorders/enforces the output schema; no rows are dropped.
  results <- results[,
    c(
      "cluster",
      "mean_proportion_control_paired",
      "mean_proportion_estim_paired",
      "effect_paired",
      "p_value_paired",
      "q_value_paired",
      "n_perm_paired",
      "direction_paired",
      "effect_paired_singleton",
      "p_value_paired_singleton",
      "q_value_paired_singleton",
      "n_perm_paired_singleton",
      "paired_singleton_status"
    )
  ]
  row.names(results) <- NULL
  results
}
