#' Plot sample-level cluster proportions by mouse.
#'
#' @param sample_props Data frame returned by
#'   `compute_sample_cluster_proportions()`.
#' @param control_label Control condition label.
#' @param estim_label E-Stim condition label.
#'
#' @return A ggplot object.
#' @export
plot_cluster_proportion_by_mouse <- function(
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

  expected_conditions <- c(control_label, estim_label)
  cluster_levels <- unique(as.character(sample_props$cluster))
  if (all(grepl("^-?[0-9]+$", cluster_levels))) {
    cluster_levels <- as.character(sort(as.integer(cluster_levels)))
  } else {
    cluster_levels <- sort(cluster_levels, method = "radix")
  }

  plot_data <- sample_props
  plot_data$mouse <- as.character(plot_data$mouse)
  plot_data$condition <- factor(
    as.character(plot_data$condition),
    levels = expected_conditions
  )
  plot_data$cluster <- factor(
    as.character(plot_data$cluster),
    levels = cluster_levels
  )
  plot_data$mouse_role <- factor(
    as.character(plot_data$mouse_role),
    levels = c("paired", "estim_only", "control_only")
  )
  if (anyNA(plot_data$condition)) {
    stop(
      "Sample proportion table contains unexpected conditions.",
      call. = FALSE
    )
  }
  if (anyNA(plot_data$cluster)) {
    stop(
      "Sample proportion table contains missing cluster labels.",
      call. = FALSE
    )
  }
  if (anyNA(plot_data$mouse_role)) {
    stop(
      "Sample proportion table contains unexpected mouse roles.",
      call. = FALSE
    )
  }
  if (
    anyNA(plot_data$proportion) ||
      any(!is.finite(plot_data$proportion)) ||
      any(plot_data$proportion < 0) ||
      any(plot_data$proportion > 1)
  ) {
    stop(
      "Sample proportion values must be finite and in [0, 1].",
      call. = FALSE
    )
  }

  paired_data <- plot_data[plot_data$mouse_role == "paired", , drop = FALSE]
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[["condition"]],
      y = .data[["proportion"]]
    )
  ) +
    ggplot2::geom_line(
      data = paired_data,
      ggplot2::aes(group = .data[["mouse"]]),
      color = "grey60",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        color = .data[["mouse"]],
        shape = .data[["mouse_role"]]
      ),
      size = 2
    ) +
    ggplot2::facet_wrap(~cluster, scales = "free_y") +
    ggplot2::scale_shape_manual(
      values = c(paired = 16, estim_only = 17, control_only = 15),
      drop = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Cluster proportion of MG-selected cells",
      color = "Mouse",
      shape = "Mouse role",
      subtitle = paste(
        "Mouse x Condition sample-level proportions;",
        "paired mice connected, singletons shown separately."
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)
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
  enrichment_order <- order(
    -plot_data$log2_enrichment,
    seq_len(nrow(plot_data))
  )
  plot_data$cluster <- factor(
    plot_data$cluster,
    levels = plot_data$cluster[enrichment_order]
  )
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

  direction_labels <- c(
    "Enriched in E-Stim" = sprintf("Enriched in %s", ESTIM_DISPLAY_LABEL),
    "Depleted in E-Stim" = sprintf("Depleted in %s", ESTIM_DISPLAY_LABEL),
    "Not significant" = "Not significant"
  )
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
        "Enriched in E-Stim" = unname(palette_analysis_three[["high"]]),
        "Depleted in E-Stim" = unname(palette_analysis_three[["low"]]),
        "Not significant" = unname(palette_analysis_three[["mid"]])
      ),
      labels = direction_labels,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.12, 0.18))
    ) +
    ggplot2::labs(
      x = "MG-selected cluster",
      y = sprintf("CLR log2 enrichment %s", CONTRAST_DISPLAY_LABEL),
      fill = NULL,
      subtitle = "Pooled cell-level Fisher/CLR summary; descriptive relative to Mouse x Condition DE unit."
    ) +
    ggplot2::theme_bw()
}
