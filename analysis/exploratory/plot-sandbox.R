# Interactive sandbox for last-mile manuscript plots.

devtools::load_all()

library(Seurat)
library(tidyverse)
library(lme4)

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

plot_gene_violin_single <- function(
  object,
  gene,
  title = NULL,
  cells = NULL
) {
  Seurat::FetchData(
    object,
    vars = c("Condition", unname(gene)),
    cells = cells,
    layer = "pflog"
  ) |>
    tibble::rownames_to_column("cell") |>
    dplyr::rename(expression = dplyr::all_of(gene)) |>
    dplyr::mutate(Condition = factor(Condition, levels = condition_levels)) |>
    ggplot(aes(x = Condition, y = expression)) +
    geom_violin(trim = TRUE, scale = "area", fill = "#d4d4d4") +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    labs(
      title = title,
      x = NULL,
      y = "Log normalized expression"
    ) +
    ESPI::theme_stone() +
    theme(axis.text.x = element_text(face = "bold"))
}

plot_gene_violin_facet <- function(object, genes, title, cells = NULL) {
  Seurat::FetchData(
    object,
    vars = c("Condition", unname(genes)),
    cells = cells,
    layer = "pflog"
  ) |>
    tibble::rownames_to_column("cell") |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(unname(genes)),
      names_to = "feature",
      values_to = "expression"
    ) |>
    dplyr::mutate(
      Condition = factor(Condition, levels = condition_levels),
      gene = factor(
        feature,
        levels = unname(genes),
        labels = names(genes)
      )
    ) |>
    ggplot(aes(x = Condition, y = expression)) +
    geom_violin(trim = TRUE, scale = "area", fill = "#d4d4d4") +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    facet_wrap(vars(gene), scales = "free_y", ncol = 4) +
    labs(
      title = title,
      x = NULL,
      y = "Normalized expression"
    ) +
    ESPI::theme_stone() +
    theme(axis.text.x = element_text(face = "bold"))
}


plot_gene_pairs <- function(object, gene_pairs, title, cells = NULL) {
  purrr::imap(
    gene_pairs,
    function(features, pair_label) {
      Seurat::FetchData(
        object,
        vars = c("Condition", features),
        cells = cells,
        layer = "data"
      ) |>
        tibble::rownames_to_column("cell") |>
        dplyr::transmute(
          cell,
          Condition = factor(Condition, levels = condition_levels),
          pair = pair_label,
          gene_1 = .data[[features[[1L]]]],
          gene_2 = .data[[features[[2L]]]]
        )
    }
  ) |>
    purrr::list_rbind() |>
    dplyr::mutate(pair = factor(pair, levels = names(gene_pairs))) |>
    ggplot(aes(x = gene_1, y = gene_2)) +
    geom_point(alpha = 0.25, size = 0.7) +
    facet_grid(
      rows = vars(pair),
      cols = vars(Condition),
      scales = "free"
    ) +
    labs(
      title = title,
      x = "First gene normalized expression",
      y = "Second gene normalized expression"
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

  all_cells_single_gene_plot <- plot_gene_expression(
    sobj,
    all_cell_genes,
    "All cells"
  )
  all_cells_gene_pair_plot <- plot_gene_pairs(
    sobj,
    all_cell_pairs,
    "All cells: paired gene expression"
  )
  cluster_5_single_gene_plot <- plot_gene_expression(
    sobj,
    cluster_5_genes,
    "Cluster 5",
    cells = cluster_5_cells
  )
  cluster_5_gene_pair_plot <- plot_gene_pairs(
    sobj,
    cluster_5_pairs,
    "Cluster 5: paired gene expression",
    cells = cluster_5_cells
  )
  cluster_4_single_gene_plot <- plot_gene_expression(
    sobj,
    cluster_4_genes,
    "Cluster 4",
    cells = cluster_4_cells
  )
  cluster_4_gene_pair_plot <- plot_gene_pairs(
    sobj,
    cluster_4_pairs,
    "Cluster 4: paired gene expression",
    cells = cluster_4_cells
  )
}

ascl1_violin <- plot_gene_violin_single(sobj, "Ascl1")

all_cells_single_gene_plot
all_cells_gene_pair_plot
cluster_5_single_gene_plot
cluster_5_gene_pair_plot
cluster_4_single_gene_plot
cluster_4_gene_pair_plot

# Combined detection per cluster
# cluster 4, Ascl1 + Otx2, Ascl1 + Hes6
c4_ascl1_otx2_dat <- Seurat::FetchData(
  sobj,
  vars = c("Sample", "Mouse", "Condition", "Ascl1", "Otx2"),
  cells = WhichCells(sobj, idents = "4"),
  layer = "counts"
) |>
  mutate(
    both = (Ascl1 > 0) & (Otx2 > 0)
  )

c4_ascl1_otx2_glm_summary <- c4_ascl1_otx2_dat |>
  group_by(Sample, Mouse, Condition) |>
  summarize(
    n_cluster = n(),
    n_both = sum(both)
  )

fit_c4_ascl1_otx2_glm <- glm(
  formula = cbind(n_both, n_cluster - n_both) ~ Condition,
  family = quasibinomial(),
  data = c4_ascl1_otx2_glm_summary
)

binomial_glm_clust_multigenes <- function(cluster, genes) {
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
    summarize(
      n_cluster = n(),
      n_all = sum(all_genes),
      .groups = "drop"
    ) |>
    glm(
      formula = cbind(n_all, n_cluster - n_all) ~ Condition,
      family = quasibinomial(),
      data = _
    )
}

c4_ascl1_hes6_glm <- binomial_glm_clust_multigenes("4", c("Ascl1", "Hes6"))

x <- binomial_glm_clust_multigenes("4", "Ascl1")
