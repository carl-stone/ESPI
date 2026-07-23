# Interactive sandbox for last-mile manuscript plots.

devtools::load_all()

library(Seurat)
library(tidyverse)
library(lme4)
library(ggview)

# ---- inputs ----

config <- publication_config()
cluster_column <- config$selected$mg$column

seurat_path <- config$selected$mg$path
deg_dir <- file.path(config$paths$degs, "mg_selected")
de_path <- file.path(deg_dir, "deseq2_full_results.tsv")
de_significant_path <- file.path(deg_dir, "deseq2_significant_degs.tsv")

input_paths <- c(seurat_path, de_path, de_significant_path)
missing_paths <- input_paths[!file.exists(input_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing input file(s):\n", paste(missing_paths, collapse = "\n"))
}

# ---- load ----

sobj <- readRDS(seurat_path)
de_results <- readr::read_tsv(de_path, show_col_types = FALSE)
de_significant <- readr::read_tsv(de_significant_path, show_col_types = FALSE)

if (!cluster_column %in% colnames(sobj[[]])) {
  stop("Missing cluster column: ", cluster_column)
}
Seurat::Idents(sobj) <- cluster_column

# ---- plots ----

condition_levels <- c(config$conditions$control, config$conditions$estim)

plot_gene_violin <- function(object, genes, title = NULL, clusters = NULL) {
  features <- base::unname(genes)
  gene_labels <- base::names(genes)

  if (length(features) == 0L) {
    stop("genes must contain at least one gene")
  }

  if (is.null(gene_labels) || anyNA(gene_labels) || any(gene_labels == "")) {
    gene_labels <- features
  }

  cells <- if (is.null(clusters)) {
    NULL
  } else {
    Seurat::WhichCells(object, idents = clusters)
  }

  plot_data <- Seurat::FetchData(
    object,
    vars = c("Condition", features),
    cells = cells,
    layer = "pflog"
  ) |>
    tibble::rownames_to_column("cell") |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(features),
      names_to = "feature",
      values_to = "expression"
    ) |>
    dplyr::mutate(
      Condition = factor(Condition, levels = condition_levels),
      gene = factor(feature, levels = features, labels = gene_labels)
    )

  plot <- plot_data |>
    ggplot(aes(x = Condition, y = expression)) +
    geom_violin(trim = TRUE, scale = "area", fill = "#d4d4d4") +
    geom_boxplot(width = 0.12, outlier.shape = NA)

  if (length(features) > 1L) {
    plot <- plot + facet_wrap(vars(gene), scales = "free_y", ncol = 3)
  }

  plot +
    labs(title = title, x = NULL, y = "Log normalized expression") +
    ESPI::theme_stone() +
    theme(axis.text.x = element_text(face = "bold"))
}

plot_gene_pair <- function(
  object,
  gene_pairs,
  title = NULL,
  clusters = NULL,
  positive_only = TRUE
) {
  pairs <- if (is.list(gene_pairs)) {
    gene_pairs
  } else {
    list(gene_pairs)
  }

  if (length(pairs) == 0L || any(lengths(pairs) != 2L)) {
    stop("gene_pairs must contain one or more two-gene pairs")
  }

  pair_labels <- names(pairs)
  if (is.null(pair_labels) || anyNA(pair_labels) || any(pair_labels == "")) {
    pair_labels <- vapply(
      pairs,
      paste,
      collapse = " + ",
      FUN.VALUE = character(1)
    )
  }

  cells <- if (is.null(clusters)) {
    NULL
  } else {
    Seurat::WhichCells(object, idents = clusters)
  }

  plot_data <- purrr::map2(pairs, pair_labels, function(features, pair_label) {
    Seurat::FetchData(
      object,
      vars = c("Condition", features),
      cells = cells,
      layer = "pflog"
    ) |>
      tibble::rownames_to_column("cell") |>
      dplyr::transmute(
        cell,
        Condition = factor(Condition, levels = condition_levels),
        pair = pair_label,
        gene_1 = .data[[features[[1L]]]],
        gene_2 = .data[[features[[2L]]]]
      )
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(pair = factor(pair, levels = pair_labels))

  if (positive_only) {
    plot_data <- plot_data |> dplyr::filter(gene_1 > 0, gene_2 > 0)
  }

  plot <- plot_data |>
    ggplot(aes(x = gene_1, y = gene_2)) +
    geom_point(alpha = 0.25, size = 0.7)

  if (length(pairs) > 1L) {
    plot <- plot +
      facet_grid(rows = vars(pair), cols = vars(Condition), scales = "free")
  } else {
    plot <- plot + facet_grid(cols = vars(Condition), scales = "free")
  }

  if (length(pairs) == 1L) {
    features <- pairs[[1L]]
    x_label <- paste(features[[1L]], "log normalized expression")
    y_label <- paste(features[[2L]], "log normalized expression")
  } else {
    x_label <- "First gene log normalized expression"
    y_label <- "Second gene log normalized expression"
  }

  plot + labs(title = title, x = x_label, y = y_label) + ESPI::theme_stone()
}

plot_gene_pair_binary <- function(
  object,
  gene_pairs,
  title = NULL,
  clusters = NULL
) {
  pairs <- if (is.list(gene_pairs)) {
    gene_pairs
  } else {
    list(gene_pairs)
  }

  if (length(pairs) == 0L || any(lengths(pairs) != 2L)) {
    stop("gene_pairs must contain one or more two-gene pairs")
  }

  pair_labels <- names(pairs)
  if (is.null(pair_labels) || anyNA(pair_labels) || any(pair_labels == "")) {
    pair_labels <- vapply(
      pairs,
      paste,
      collapse = " + ",
      FUN.VALUE = character(1)
    )
  }

  cells <- if (is.null(clusters)) {
    NULL
  } else {
    Seurat::WhichCells(object, idents = clusters)
  }

  state_levels <- c("Neither", "Gene 1 only", "Gene 2 only", "Both")

  plot_data <- purrr::map2(pairs, pair_labels, function(features, pair_label) {
    Seurat::FetchData(
      object,
      vars = c("Sample", "Condition", features),
      cells = cells,
      layer = "counts"
    ) |>
      tibble::rownames_to_column("cell") |>
      dplyr::transmute(
        cell,
        Sample,
        Condition = factor(Condition, levels = condition_levels),
        pair = pair_label,
        state = dplyr::case_when(
          .data[[features[[1L]]]] > 0 & .data[[features[[2L]]]] > 0 ~ "Both",
          .data[[features[[1L]]]] > 0 ~ "Gene 1 only",
          .data[[features[[2L]]]] > 0 ~ "Gene 2 only",
          TRUE ~ "Neither"
        )
      )
  }) |>
    purrr::list_rbind() |>
    dplyr::count(Sample, Condition, pair, state, name = "n") |>
    tidyr::complete(
      tidyr::nesting(Sample, Condition, pair),
      state = state_levels,
      fill = list(n = 0L)
    ) |>
    dplyr::group_by(Sample, Condition, pair) |>
    dplyr::mutate(
      state = factor(state, levels = state_levels),
      proportion = n / sum(n)
    ) |>
    dplyr::ungroup()

  plot <- ggplot(
    plot_data,
    aes(x = state, y = proportion, color = Condition, group = Condition)
  ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65),
      size = 2,
      alpha = 0.75
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      position = position_dodge(width = 0.65),
      size = 3,
      shape = 18
    ) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      position = position_dodge(width = 0.65),
      width = 0.2
    )

  if (length(pairs) > 1L) {
    plot <- plot + facet_wrap(vars(pair), ncol = 1)
  }

  plot +
    scale_y_continuous(labels = scales::label_percent()) +
    labs(
      title = title,
      x = NULL,
      y = "Proportion of cells",
      color = "Condition"
    ) +
    ESPI::theme_stone() +
    theme(axis.text.x = element_text(face = "bold"))
}

plot_gene_pair_expression <- function(
  object,
  genes,
  title = NULL,
  clusters = NULL,
  smooth = list()
) {
  if (!is.list(smooth)) {
    stop("smooth must be a list of geom_smooth() arguments")
  }

  smooth <- utils::modifyList(
    list(method = MASS::rlm, se = FALSE, linewidth = 0.6, na.rm = TRUE),
    smooth
  )

  features <- base::unname(genes)
  gene_labels <- base::names(genes)

  if (length(features) < 2L) {
    stop("genes must contain at least two genes")
  }

  if (is.null(gene_labels) || anyNA(gene_labels) || any(gene_labels == "")) {
    gene_labels <- features
  }

  cells <- if (is.null(clusters)) {
    NULL
  } else {
    Seurat::WhichCells(object, idents = clusters)
  }

  expression_data <- Seurat::FetchData(
    object,
    vars = c("Condition", features),
    cells = cells,
    layer = "pflog"
  ) |>
    tibble::rownames_to_column("cell") |>
    dplyr::mutate(Condition = factor(Condition, levels = condition_levels))

  pair_grid <- tidyr::expand_grid(x_feature = features, y_feature = features)

  plot_data <- purrr::map2(
    pair_grid$x_feature,
    pair_grid$y_feature,
    function(x_feature, y_feature) {
      expression_data |>
        dplyr::transmute(
          cell,
          Condition,
          x_feature,
          y_feature,
          x_expression = .data[[x_feature]],
          y_expression = .data[[y_feature]]
        )
    }
  ) |>
    purrr::list_rbind() |>
    dplyr::mutate(
      x_gene = factor(x_feature, levels = features, labels = gene_labels),
      y_gene = factor(y_feature, levels = features, labels = gene_labels)
    )

  plot_data <- plot_data |>
    dplyr::mutate(
      x_index = match(x_feature, features),
      y_index = match(y_feature, features)
    )

  diagonal <- plot_data |> dplyr::filter(x_index == y_index)
  upper_triangle <- plot_data |> dplyr::filter(y_index < x_index)
  lower_triangle <- plot_data |> dplyr::filter(y_index > x_index)

  correlations <- lower_triangle |>
    dplyr::group_by(x_feature, y_feature, x_gene, y_gene, Condition) |>
    dplyr::summarise(
      rho = if (dplyr::n() > 1L) {
        suppressWarnings(stats::cor(
          x_expression,
          y_expression,
          method = "spearman",
          use = "complete.obs"
        ))
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      label = paste0("ρ = ", formatC(rho, format = "f", digits = 2)),
      vjust = 1.2 + (as.integer(Condition) - 1L) * 1.5
    )

  condition_colors <- stats::setNames(
    config$palettes$dotplot[seq_along(condition_levels)],
    condition_levels
  )

  smooth_layer <- do.call(
    ggplot2::geom_smooth,
    c(
      list(
        data = upper_triangle,
        mapping = ggplot2::aes(
          x = x_expression,
          y = y_expression,
          color = Condition,
          group = Condition
        )
      ),
      smooth
    )
  )

  ggplot(plot_data) +
    geom_histogram(
      data = diagonal,
      aes(
        x = x_expression,
        y = after_stat(density),
        fill = Condition,
        color = Condition
      ),
      bins = 30,
      position = "identity",
      alpha = 0.3,
      na.rm = TRUE
    ) +
    geom_point(
      data = upper_triangle,
      aes(x = x_expression, y = y_expression, color = Condition),
      alpha = 0.15,
      size = 0.5,
      na.rm = TRUE
    ) +
    smooth_layer +
    geom_text(
      data = correlations,
      aes(x = Inf, y = Inf, label = label, color = Condition, vjust = vjust),
      hjust = 1.1,
      na.rm = TRUE
    ) +
    scale_color_manual(values = condition_colors) +
    scale_fill_manual(values = condition_colors) +
    facet_grid(rows = vars(y_gene), cols = vars(x_gene), scales = "free") +
    labs(
      title = title,
      x = "Log normalized expression",
      y = "Log normalized expression",
      color = "Condition",
      fill = "Condition"
    ) +
    ESPI::theme_stone()
}

{
  all_cell_genes <- c(
    OTX2 = "Otx2",
    Ascl1 = "Ascl1",
    Hes6 = "Hes6",
    SCGN = "Scgn",
    LHX4 = "Lhx4",
    VSX2 = "Vsx2",
    Grik1 = "Grik1",
    CABP5 = "Cabp5",
    Dll1 = "Dll1",
    Neurod1 = "Neurod1",
    Neurog2 = "Neurog2"
  )
  all_cell_pairs <- list(
    `ASCL1 + OTX2` = c("Ascl1", "Otx2"),
    `OTX2 + CABP5` = c("Otx2", "Cabp5")
  )

  cluster_5_cells <- Seurat::WhichCells(sobj, idents = "5")
  cluster_5_genes <- all_cell_genes[c(
    "OTX2",
    "SCGN",
    "LHX4",
    "VSX2",
    "Grik1",
    "CABP5",
    "Dll1",
    "Neurod1"
  )]
  cluster_5_pairs <- list(
    `OTX2 + CABP5` = c("Otx2", "Cabp5"),
    `OTX2 + Hes6` = c("Otx2", "Hes6")
  )

  cluster_4_cells <- Seurat::WhichCells(sobj, idents = "4")
  cluster_4_genes <- all_cell_genes[c("Ascl1", "Hes6", "Neurog2", "OTX2")]
  cluster_4_pairs <- list(
    `ASCL1 + OTX2` = c("Ascl1", "Otx2"),
    `ASCL1 + Hes6` = c("Ascl1", "Hes6")
  )

  all_cells_single_gene_plot <- plot_gene_violin(
    sobj,
    all_cell_genes,
    "All cells"
  )
  all_cells_gene_pair_plot <- plot_gene_pair(
    sobj,
    all_cell_pairs,
    "All cells: positive paired gene expression"
  )
  all_cells_gene_pair_binary_plot <- plot_gene_pair_binary(
    sobj,
    all_cell_pairs,
    "All cells: paired gene detection"
  )
  all_cells_gene_pair_expression_plot <- plot_gene_pair_expression(
    sobj,
    genes = all_cell_genes,
    title = "All cells: pairwise gene expression",
    smooth = list(method = "loess")
  )
  cluster_5_single_gene_plot <- plot_gene_violin(
    sobj,
    cluster_5_genes,
    "Cluster 5",
    clusters = "5"
  )
  cluster_5_gene_pair_plot <- plot_gene_pair(
    sobj,
    cluster_5_pairs,
    "Cluster 5: positive paired gene expression",
    clusters = "5"
  )
  cluster_5_gene_pair_binary_plot <- plot_gene_pair_binary(
    sobj,
    cluster_5_pairs,
    "Cluster 5: paired gene detection",
    clusters = "5"
  )
  cluster_5_gene_pair_expression_plot <- plot_gene_pair_expression(
    sobj,
    cluster_5_genes,
    "Cluster 5: pairwise gene expression",
    clusters = "5"
  )
  cluster_4_single_gene_plot <- plot_gene_violin(
    sobj,
    cluster_4_genes,
    "Cluster 4",
    clusters = "4"
  )
  cluster_4_gene_pair_plot <- plot_gene_pair(
    sobj,
    cluster_4_pairs,
    "Cluster 4: positive paired gene expression",
    clusters = "4"
  )
  cluster_4_gene_pair_binary_plot <- plot_gene_pair_binary(
    sobj,
    cluster_4_pairs,
    "Cluster 4: paired gene detection",
    clusters = "4"
  )
  cluster_4_gene_pair_expression_plot <- plot_gene_pair_expression(
    sobj,
    cluster_4_genes,
    "Cluster 4: pairwise gene expression",
    clusters = "4"
  )
}

ascl1_violin <- plot_gene_violin(sobj, "Ascl1")

c5_ascl1_otx2_plot <- plot_gene_pair(
  sobj,
  c("Ascl1", "Otx2"),
  title = "Cluster 5: positive ASCL1-OTX2 expression",
  clusters = "5"
) +
  geom_smooth()

c5_ascl1_otx2_binary_plot <- plot_gene_pair_binary(
  sobj,
  c("Ascl1", "Otx2"),
  title = "Cluster 5: ASCL1-OTX2 detection",
  clusters = "5"
)

c5_ascl1_otx2_plot

c4_ascl1_otx2_plot <- plot_gene_pair(
  sobj,
  c("Ascl1", "Otx2"),
  title = "Cluster 4: positive ASCL1-OTX2 expression",
  clusters = "4"
) +
  geom_smooth()

c4_ascl1_otx2_binary_plot <- plot_gene_pair_binary(
  sobj,
  c("Ascl1", "Otx2"),
  title = "Cluster 4: ASCL1-OTX2 detection",
  clusters = "4"
)

all_cells_single_gene_plot + canvas(7, 7)
all_cells_gene_pair_plot + geom_smooth(method = MASS::rlm) + canvas(4.5, 4)
all_cells_gene_pair_binary_plot + canvas(6, 5)
all_cells_gene_pair_expression_plot + canvas(10, 10)
cluster_5_single_gene_plot + canvas(7, 5)
cluster_5_gene_pair_plot + geom_smooth(method = MASS::rlm) + canvas(4.5, 4)
cluster_5_gene_pair_binary_plot + canvas(6, 5)
cluster_5_gene_pair_expression_plot + canvas(8.5, 8)
cluster_4_single_gene_plot + canvas(6.5, 5)
cluster_4_gene_pair_plot + geom_smooth(method = MASS::rlm) + canvas(4.5, 4)
cluster_4_gene_pair_binary_plot + canvas(6, 5)
cluster_4_gene_pair_expression_plot + canvas(6, 5)
c5_ascl1_otx2_binary_plot + canvas(6, 4)
c4_ascl1_otx2_binary_plot + canvas(6, 4)

fig_path <- file.path(config$paths$figures, "random_pairs")

# Save plots
ggsave(
  file.path(fig_path, "all_gene_violins.pdf"),
  all_cells_single_gene_plot,
  width = 7,
  height = 7
)

ggsave(
  file.path(fig_path, "all_ascl1_otx2_capb5_scatter.pdf"),
  all_cells_gene_pair_plot + geom_smooth(method = MASS::rlm),
  width = 4.5,
  height = 4
)

ggsave(
  file.path(fig_path, "all_ascl1_otx2_capb5_detection.pdf"),
  all_cells_gene_pair_binary_plot,
  width = 6,
  height = 5
)

ggsave(
  file.path(fig_path, "all_pairs_plot.pdf"),
  all_cells_gene_pair_expression_plot,
  width = 10,
  height = 10
)

ggsave(
  file.path(fig_path, "c5_gene_violins.pdf"),
  cluster_5_single_gene_plot,
  width = 7,
  height = 5
)

ggsave(
  file.path(fig_path, "c5_otx2_capb5_otx2_hes6_scatter.pdf"),
  cluster_5_gene_pair_plot + geom_smooth(method = MASS::rlm),
  width = 4.5,
  height = 4
)

ggsave(
  file.path(fig_path, "c5_otx2_capb5_otx2_hes6_detection.pdf"),
  cluster_5_gene_pair_binary_plot,
  width = 6,
  height = 5
)

ggsave(
  file.path(fig_path, "c5_pairs_plot.pdf"),
  cluster_5_gene_pair_expression_plot,
  width = 8.5,
  height = 8
)

ggsave(
  file.path(fig_path, "c4_gene_violins.pdf"),
  cluster_4_single_gene_plot,
  width = 6.5,
  height = 5
)

ggsave(
  file.path(fig_path, "c4_ascl1_otx2_ascl1_hes6_scatter.pdf"),
  cluster_4_gene_pair_plot + geom_smooth(method = MASS::rlm),
  width = 4.5,
  height = 4
)

ggsave(
  file.path(fig_path, "c4_ascl1_otx2_ascl1_hes6_detection.pdf"),
  cluster_4_gene_pair_binary_plot,
  width = 6,
  height = 5
)

ggsave(
  file.path(fig_path, "c4_pairs_plot.pdf"),
  cluster_4_gene_pair_expression_plot,
  width = 6,
  height = 5
)

ggsave(
  file.path(fig_path, "c5_ascl1_otx2_detection.pdf"),
  c5_ascl1_otx2_binary_plot,
  width = 6,
  height = 4
)

ggsave(
  file.path(fig_path, "c4_ascl1_otx2_detection.pdf"),
  c4_ascl1_otx2_binary_plot,
  width = 6,
  height = 4
)

# Do stats on cluster 5 expression of Cabp5, Scgn, Dll1, Otx2, Neurod1

# Ascl1 + Otx2 co-expression test by cluster
ascl1_otx2_coexp_dat <- Seurat::FetchData(
  sobj,
  vars = c("Sample", "Mouse", "ident", "Condition", "Ascl1", "Otx2"),
  layer = "counts",
  cells = Seurat::WhichCells(sobj, idents = c(1, 2, 3, 4, 5, 6, 7))
) |>
  # dplyr::filter(ident != "8") |>
  mutate(both = (Ascl1 > 0) & (Otx2 > 0)) |>
  group_by(Sample, Mouse, Condition, ident) |>
  summarize(n_cluster = n(), n_both = sum(both))

fit_ascl1_otx2_additive <- glm(
  formula = cbind(n_both, n_cluster - n_both) ~ Condition + ident + Sample,
  family = quasibinomial(),
  data = ascl1_otx2_coexp_dat
)

fit_ascl1_otx2_bycluster <- glm(
  formula = cbind(n_both, n_cluster - n_both) ~ Condition * ident + Sample,
  family = quasibinomial(),
  data = ascl1_otx2_coexp_dat
)

fit_ascl1_otx2_condition_bycluster_test <- anova(
  fit_ascl1_otx2_additive,
  fit_ascl1_otx2_bycluster,
  test = "F"
)

# Ascl1 + Hes6 co-expression test by cluster
ascl1_hes6_coexp_dat <- Seurat::FetchData(
  sobj,
  vars = c("Sample", "Mouse", "ident", "Condition", "Ascl1", "Hes6"),
  layer = "counts",
  cells = Seurat::WhichCells(sobj, idents = c(1, 2, 3, 4, 5, 6, 7))
) |>
  # dplyr::filter(ident != "8") |>
  mutate(both = (Ascl1 > 0) & (Hes6 > 0)) |>
  group_by(Sample, Mouse, Condition, ident) |>
  summarize(n_cluster = n(), n_both = sum(both))

fit_ascl1_hes6_additive <- glmer(
  formula = cbind(n_both, n_cluster - n_both) ~ Condition +
    ident +
    (1 | Sample),
  family = binomial(),
  data = ascl1_hes6_coexp_dat
)

fit_ascl1_hes6_bycluster <- glmer(
  formula = cbind(n_both, n_cluster - n_both) ~ Condition *
    ident +
    (1 | Sample),
  family = binomial(),
  data = ascl1_hes6_coexp_dat
)

fit_ascl1_hes6_condition_bycluster_test <- anova(
  fit_ascl1_hes6_additive,
  fit_ascl1_hes6_bycluster,
  test = "F"
)


# Combined detection per cluster
# cluster 4, Ascl1 + Otx2, Ascl1 + Hes6
c4_ascl1_otx2_dat <- Seurat::FetchData(
  sobj,
  vars = c("Sample", "Mouse", "Condition", "Ascl1", "Otx2"),
  cells = WhichCells(sobj, idents = "4"),
  layer = "counts"
) |>
  mutate(both = (Ascl1 > 0) & (Otx2 > 0))

c4_ascl1_otx2_glm_summary <- c4_ascl1_otx2_dat |>
  group_by(Sample, Mouse, Condition) |>
  summarize(n_cluster = n(), n_both = sum(both))

fit_c4_ascl1_otx2_glm <- glm(
  formula = cbind(n_both, n_cluster - n_both) ~ Condition,
  family = quasibinomial(),
  data = c4_ascl1_otx2_glm_summary
)

binomial_glm_clust_multigenes <- function(
  cluster = seq_along(unique(Idents(sobj))),
  genes
) {
  Seurat::FetchData(
    sobj,
    vars = c("Sample", "Mouse", "Condition", genes),
    cells = WhichCells(sobj, idents = cluster),
    layer = "counts"
  ) |>
    mutate(
      across(all_of(genes), \(x) x > 0),
      all_genes = if_all(all_of(genes))
    ) |>
    group_by(Sample, Mouse, Condition) |>
    summarize(n_cluster = n(), n_all = sum(all_genes), .groups = "drop") |>
    glm(
      formula = cbind(n_all, n_cluster - n_all) ~ Condition,
      family = quasibinomial(),
      data = _
    )
}

# focused stats
# c4 & c5 scatter
# all ascl1, otx2, cabp5 scatter
# c4 and c5 violins
# for pairs plots and violins, possibly remove cells with 0's

# c4 ascl1 and otx2 binary expression
c4_ascl1_otx2_glm <- binomial_glm_clust_multigenes("4", c("Ascl1", "Otx2"))
c4_ascl1_hes6_glm <- binomial_glm_clust_multigenes("4", c("Ascl1", "Hes6"))

all_ascl1_otx2_glm <- binomial_glm_clust_multigenes(genes = c("Ascl1", "Otx2"))

# c5
c5_otx2_cabp5_glm <- binomial_glm_clust_multigenes("5", c("Otx2", "Cabp5"))
c5_otx2_hes6_glm <- binomial_glm_clust_multigenes("5", c("Otx2", "Hes6"))
