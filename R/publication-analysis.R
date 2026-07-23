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
# ANALYSIS_OK[R026]: exported abundance computation is called by the publication-figures phase script.
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
    !paired_mouse & present_mat[, estim_label] & !present_mat[, control_label]
  ] <- "estim_only"
  mouse_roles[
    !paired_mouse & present_mat[, control_label] & !present_mat[, estim_label]
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

  output_rows <- vector("list", sum(present_mat) * length(cluster_levels))
  row_idx <- 1L
  for (condition_label in expected_conditions) {
    for (mouse_id in mouse_levels) {
      sample_total <- as.integer(sample_tab[mouse_id, condition_label])
      if (sample_total == 0L) {
        next
      }
      for (cluster_label in cluster_levels) {
        cluster_n <- as.integer(counts_array[
          mouse_id,
          condition_label,
          cluster_label
        ])
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
  if (length(k) != 1L || is.na(k) || k < 1L || k != as.integer(k)) {
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
    p_value_paired <- mean(abs(paired_null) >= abs(effect_paired) - 1e-8)

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
    mean_proportion_control_paired <- mean(cluster_data$proportion[
      control_prop_idx
    ])
    mean_proportion_estim_paired <- mean(cluster_data$proportion[
      estim_prop_idx
    ])

    if (singleton_runs) {
      singleton_diff <- get_logit(estim_only_mice[[1L]], estim_label) -
        get_logit(control_only_mice[[1L]], control_label)
      paired_singleton_contrasts <- c(diffs, singleton_diff)
      paired_singleton_signs <- .sign_vectors(length(
        paired_singleton_contrasts
      ))
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
  results <- results[, c(
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
  )]
  row.names(results) <- NULL
  results
}

#' Compute cluster-level cell-type marker module scores.
#'
#' Scores each cell for each named marker-gene set with
#' [Seurat::AddModuleScore()], then averages the module scores within clusters.
#'
#' @param sobj Seurat object containing expression and cluster metadata.
#' @param cluster_col Metadata column containing cluster labels.
#' @param marker_genes Named non-empty list of marker-gene character vectors.
#' @param assay Assay used for module scoring.
#' @param slot Assay slot/layer used by [Seurat::AddModuleScore()].
#' @param seed Random seed forwarded to [Seurat::AddModuleScore()].
#'
#' @return Numeric matrix with marker-gene sets in rows, sorted clusters in
#'   columns, and mean module scores as values.
#' @export
compute_cluster_module_scores <- function(
  sobj,
  cluster_col,
  marker_genes,
  # ANALYSIS_OK[smuggled-default]: package API default selects the RNA assay.
  assay = "RNA",
  # ANALYSIS_OK[smuggled-default]: package API default selects the data layer.
  slot = "data",
  seed = SEED
) {
  meta <- sobj@meta.data
  if (!cluster_col %in% colnames(meta)) {
    stop("Missing cluster metadata column: ", cluster_col, call. = FALSE)
  }
  if (!assay %in% SeuratObject::Assays(sobj)) {
    stop("Missing assay: ", assay, call. = FALSE)
  }
  if (!slot %in% SeuratObject::Layers(sobj[[assay]])) {
    stop(
      "Missing score slot/layer '",
      slot,
      "' in assay ",
      assay,
      ". Available layers: ",
      paste(SeuratObject::Layers(sobj[[assay]]), collapse = ", "),
      call. = FALSE
    )
  }
  if (
    !is.list(marker_genes) ||
      length(marker_genes) == 0L ||
      is.null(names(marker_genes)) ||
      anyNA(names(marker_genes)) ||
      any(!nzchar(names(marker_genes)))
  ) {
    stop("marker_genes must be a named non-empty list.", call. = FALSE)
  }
  if (any(duplicated(names(marker_genes)))) {
    stop("marker_genes names must be unique.", call. = FALSE)
  }
  invalid_sets <- names(marker_genes)[
    !vapply(
      marker_genes,
      function(x) {
        is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x))
      },
      logical(1L)
    )
  ]
  if (length(invalid_sets) > 0L) {
    stop(
      "Marker gene set(s) must be non-empty character vectors: ",
      paste(invalid_sets, collapse = ", "),
      call. = FALSE
    )
  }

  cluster <- as.character(meta[[cluster_col]])
  if (length(cluster) != nrow(meta)) {
    stop("Cluster metadata is not aligned to cells.", call. = FALSE)
  }
  if (anyNA(cluster) || any(!nzchar(cluster))) {
    stop(
      "Cluster metadata contains missing or empty labels: ",
      cluster_col,
      call. = FALSE
    )
  }

  marker_vector <- unlist(marker_genes, use.names = FALSE)
  missing_genes <- setdiff(unique(marker_vector), rownames(sobj))
  if (length(missing_genes) > 0L) {
    stop(
      "Marker gene(s) missing from the Seurat object: ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  cluster_levels <- .sort_cluster_labels(cluster)
  module_prefix <- "celltype_module_score"
  module_score_cols <- paste0(module_prefix, seq_along(marker_genes))
  existing_module_cols <- intersect(module_score_cols, colnames(sobj@meta.data))
  if (length(existing_module_cols) > 0L) {
    sobj@meta.data[existing_module_cols] <- NULL
  }
  sobj <- Seurat::AddModuleScore(
    object = sobj,
    features = marker_genes,
    assay = assay,
    name = module_prefix,
    seed = seed,
    search = FALSE,
    slot = slot
  )
  missing_module_cols <- setdiff(module_score_cols, colnames(sobj@meta.data))
  if (length(missing_module_cols) > 0L) {
    stop(
      "Seurat::AddModuleScore did not create expected column(s): ",
      paste(missing_module_cols, collapse = ", "),
      call. = FALSE
    )
  }

  module_scores <- sobj@meta.data[module_score_cols]
  colnames(module_scores) <- names(marker_genes)
  module_scores$cluster <- cluster
  cluster_marker_scores <- stats::aggregate(
    module_scores[names(marker_genes)],
    by = list(cluster = module_scores$cluster),
    FUN = mean
  )
  # ANALYSIS_OK[row-order]: reorders aggregate rows to the validated cluster_levels; no rows are dropped and the output schema is preserved.
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

#' Compute sample-aware cluster p27 enrichment z-scores.
#'
#' Computes observed per-cluster p27 expression, then estimates a null
#' distribution by permuting cluster labels within each Mouse x Condition
#' sample. The permutation preserves each sample's p27 values and cluster-size
#' composition. A cluster receives `NA` when its permutation null has zero or
#' non-finite standard deviation, which can occur when p27 is constant within
#' every sample.
#'
#' @param sobj Seurat object containing expression and metadata.
#' @param cluster_col Metadata column containing cluster labels.
#' @param gene Gene used for p27 expression. Defaults to `Cdkn1b`.
#' @param layer Assay layer used for expression values.
#' @param assay Assay containing `layer`.
#' @param mouse_col Metadata column containing mouse labels.
#' @param condition_col Metadata column containing condition labels.
#' @param n_perm Positive integer number of within-sample permutations.
#' @param seed Random seed used for permutations.
#'
#' @return Data frame with one row per sorted cluster and columns `cluster`,
#'   `n_cells`, `observed_mean`, `null_mean`, `null_sd`, and `z_score`.
#' @export
compute_cluster_p27_enrichment <- function(
  sobj,
  cluster_col,
  # ANALYSIS_OK[smuggled-default]: package API default selects the p27 gene.
  gene = "Cdkn1b",
  # ANALYSIS_OK[smuggled-default]: package API default selects the pflog layer.
  layer = "pflog",
  # ANALYSIS_OK[smuggled-default]: package API default selects the RNA assay.
  assay = "RNA",
  # ANALYSIS_OK[smuggled-default]: package API default selects the Mouse metadata column.
  mouse_col = "Mouse",
  condition_col = CONDITION_COL,
  # ANALYSIS_OK[smuggled-default]: package API default sets permutation count to 2000.
  n_perm = 2000L,
  seed = SEED
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
  if (!assay %in% SeuratObject::Assays(sobj)) {
    stop("Missing assay: ", assay, call. = FALSE)
  }
  if (!gene %in% rownames(sobj)) {
    stop("Gene is missing from the Seurat object: ", gene, call. = FALSE)
  }
  available_layers <- SeuratObject::Layers(sobj[[assay]])
  if (!layer %in% available_layers) {
    stop(
      "Missing layer '",
      layer,
      "' in assay ",
      assay,
      ". Available layers: ",
      paste(available_layers, collapse = ", "),
      call. = FALSE
    )
  }
  if (
    length(n_perm) != 1L ||
      is.na(n_perm) ||
      n_perm < 1L ||
      n_perm != as.integer(n_perm)
  ) {
    stop("n_perm must be a positive integer.", call. = FALSE)
  }
  n_perm <- as.integer(n_perm)

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

  expression_matrix <- SeuratObject::LayerData(sobj[[assay]], layer = layer)
  expr <- as.numeric(expression_matrix[gene, colnames(sobj), drop = TRUE])
  if (length(expr) != nrow(meta)) {
    stop("Expression values are not aligned to cells.", call. = FALSE)
  }
  if (anyNA(expr) || any(!is.finite(expr))) {
    stop(
      "Expression values contain missing or non-finite values.",
      call. = FALSE
    )
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

  had_random_seed <- exists(
    ".Random.seed",
    envir = globalenv(),
    inherits = FALSE
  )
  if (had_random_seed) {
    saved_random_seed <- get(
      ".Random.seed",
      envir = globalenv(),
      inherits = FALSE
    )
  }
  on.exit(
    {
      if (had_random_seed) {
        assign(".Random.seed", saved_random_seed, envir = globalenv())
      } else if (
        exists(".Random.seed", envir = globalenv(), inherits = FALSE)
      ) {
        rm(".Random.seed", envir = globalenv())
      }
    },
    add = TRUE
  )

  # ANALYSIS_OK[random-seed-only]: RNG is scoped to permutation null generation; on.exit restores prior state.
  set.seed(seed)
  for (perm_idx in seq_len(n_perm)) {
    permuted_cluster <- cluster
    for (idx in sample_indices) {
      # ANALYSIS_OK[random-seed-only]: RNG is scoped to permutation null generation; on.exit restores prior state.
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
  degenerate_null <- !is.finite(null_sd) | null_sd == 0
  z_score[degenerate_null] <- NA_real_

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

# ANALYSIS_OK[R026]: package helper is loaded by devtools::load_all and called by same-file marker-score entrypoints.
.sort_cluster_labels <- function(x) {
  x <- unique(as.character(x))
  if (all(grepl("^-?[0-9]+$", x))) {
    return(as.character(sort(as.integer(x))))
  }
  sort(x, method = "radix")
}
