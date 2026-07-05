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
    nrow(counts_mat) == 0L || ncol(counts_mat) != 2L || sum(counts_mat) == 0L
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

cluster_abundance_significance_label <- function(p_value) {
  labels <- rep("", length(p_value))
  labels[p_value < 0.05] <- "*"
  labels[p_value < 0.01] <- "**"
  labels[p_value < 0.001] <- "***"
  labels[is.na(p_value)] <- ""
  labels
}

#' Plot CLR/Fisher cluster enrichment.
#'
#' @param enrichment Data frame returned by `compute_cluster_abundance()`.
#'
#' @return A ggplot object.
#' @export
plot_clr_fisher_enrichment <- function(enrichment) {
  required_cols <- c("cluster", "log2_enrichment", "padj", "direction")
  missing_cols <- setdiff(required_cols, colnames(enrichment))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing enrichment column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(enrichment) == 0L) {
    stop("Enrichment table is empty.", call. = FALSE)
  }
  if (any(!is.finite(enrichment$log2_enrichment))) {
    stop("Enrichment table contains non-finite log2 effects.", call. = FALSE)
  }

  plot_data <- enrichment
  plot_data$cluster <- factor(plot_data$cluster, levels = plot_data$cluster)
  plot_data$direction <- factor(
    plot_data$direction,
    levels = c("Enriched in E-Stim", "Depleted in E-Stim", "Not significant")
  )
  plot_data$significance <- ifelse(
    plot_data$padj < 0.05,
    cluster_abundance_significance_label(plot_data$padj),
    ""
  )
  y_span <- diff(range(plot_data$log2_enrichment, na.rm = TRUE))
  label_offset <- max(0.12, 0.05 * y_span)
  plot_data$label_y <- plot_data$log2_enrichment +
    ifelse(plot_data$log2_enrichment >= 0, label_offset, -label_offset)
  plot_data$label_vjust <- ifelse(plot_data$log2_enrichment >= 0, 0, 1)

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[["cluster"]],
      y = .data[["log2_enrichment"]],
      fill = .data[["direction"]]
    )
  ) +
    ggplot2::geom_col(color = "black", width = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(
        y = .data[["label_y"]],
        label = .data[["significance"]],
        vjust = .data[["label_vjust"]]
      ),
      size = 4,
      show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::scale_fill_manual(
      values = c(
        "Enriched in E-Stim" = "#e31a8c",
        "Depleted in E-Stim" = "#2166ac",
        "Not significant" = "grey75"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.12, 0.18))
    ) +
    ggplot2::labs(
      x = "MG-selected cluster",
      y = "CLR log2 enrichment (E-Stim vs control)",
      fill = NULL,
      subtitle = "Pooled cell-level Fisher/CLR summary; descriptive relative to Mouse x Condition DE unit."
    ) +
    ggplot2::theme_bw()
}
