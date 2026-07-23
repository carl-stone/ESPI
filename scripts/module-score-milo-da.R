#!/usr/bin/env Rscript

# Milo neighborhood differential-abundance analysis in four-module score space.
#
# Cells define local neighborhoods using standardized progenitor, cone bipolar,
# Müller glia, and proliferation scores. Mouse × Condition samples remain the
# units used for differential-abundance inference.
#
# Usage: ESPI_OVERWRITE=false ESPI_MODULE_SCORE_MILO_K=50 \
#   ESPI_MODULE_SCORE_MILO_PROP=0.1 Rscript scripts/module-score-milo-da.R

suppressPackageStartupMessages({
  library(here)
  here::i_am("scripts/module-score-milo-da.R")
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(tidyverse)
  library(Seurat)
  library(SingleCellExperiment)
  library(miloR)
})

# ---- parameters ----

config <- publication_config()
input_path <- config$selected$mg$path
condition_col <- config$conditions$column
control_label <- config$conditions$control
estim_label <- config$conditions$estim
seed <- config$seed

score_reduction <- "MODULE_SCORE_4D"
score_dimensions <- c(
  "progenitor_score",
  "cone_bipolar_score",
  "muller_score",
  "proliferation_score"
)
k_neighbors <- as.integer(Sys.getenv("ESPI_MODULE_SCORE_MILO_K", unset = "60"))
nhood_proportion <- as.numeric(Sys.getenv(
  "ESPI_MODULE_SCORE_MILO_PROP",
  unset = "0.04"
))
spatial_fdr_threshold <- 0.1

if (length(k_neighbors) != 1L || is.na(k_neighbors) || k_neighbors < 2L) {
  cli::cli_abort("ESPI_MODULE_SCORE_MILO_K must be one integer of at least 2.")
}
if (
  length(nhood_proportion) != 1L ||
    is.na(nhood_proportion) ||
    nhood_proportion <= 0 ||
    nhood_proportion >= 1
) {
  cli::cli_abort("ESPI_MODULE_SCORE_MILO_PROP must be between 0 and 1.")
}

marker_modules <- list(
  progenitor = cell_type_marker_genes$neurogenic_progenitor,
  cone_bipolar = cell_type_marker_genes$cone_bipolar,
  muller = cell_type_marker_genes$muller_glia,
  proliferation = cell_type_marker_genes$proliferative
)

# ---- paths ----

default_output_dir <- file.path(
  config$paths$degs,
  "mg_selected",
  "module_score_milo_da",
  paste0("k_", k_neighbors, "__prop_", nhood_proportion)
)
output_dir <- Sys.getenv(
  "ESPI_MODULE_SCORE_MILO_OUTPUT_DIR",
  unset = default_output_dir
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_paths <- file.path(
  output_dir,
  c(
    "sample_table.tsv",
    "module_scores.tsv",
    "marker_features.tsv",
    "nhood_da_results.tsv",
    "nhood_score_profiles.tsv",
    "nhood_cell_membership.tsv",
    "sample_nhood_abundance.tsv",
    "score_pca_loadings.tsv",
    "parameters.tsv",
    "nhood_da_score_space.png",
    "nhood_da_score_space.pdf",
    "nhood_component_profiles.png",
    "nhood_component_profiles.pdf"
  )
)
assert_output_available(output_paths, config$overwrite)

# ---- input and module scores ----

sobj <- readRDS(input_path)
assert_frozen_input(input_path, sobj, config$frozen$mg)

if (!condition_col %in% colnames(sobj[[]])) {
  cli::cli_abort(
    "The configured condition column is absent from cell metadata."
  )
}

rna_features <- rownames(sobj[["RNA"]])
marker_features <- purrr::imap(marker_modules, \(genes, module) {
  tibble::tibble(module, gene = genes, present_in_rna = genes %in% rna_features)
}) |>
  purrr::list_rbind()

missing_modules <- marker_features |>
  dplyr::group_by(module) |>
  dplyr::summarise(n_present = sum(present_in_rna), .groups = "drop") |>
  dplyr::filter(n_present == 0L) |>
  dplyr::pull(module)

if (length(missing_modules) > 0L) {
  cli::cli_abort(c(
    "No marker genes are present for one or more modules.",
    "x" = "Missing modules: {paste(missing_modules, collapse = ', ')}"
  ))
}

score_features <- purrr::map(marker_modules, \(genes) {
  base::intersect(genes, rna_features)
})

set.seed(seed)
sobj <- Seurat::AddModuleScore(
  object = sobj,
  features = score_features,
  name = paste0("score_space_", names(score_features)),
  assay = "RNA",
  search = FALSE,
  seed = seed
)

generated_score_columns <- stats::setNames(
  paste0("score_space_", names(score_features), seq_along(score_features)),
  names(score_features)
)
if (!all(generated_score_columns %in% colnames(sobj[[]]))) {
  cli::cli_abort("AddModuleScore did not create the expected score columns.")
}

# ANALYSIS_OK[R026]: local helper standardizes all four module-score columns.
standardize_score <- function(x) {
  score_sd <- stats::sd(x)
  if (!is.finite(score_sd) || score_sd == 0) {
    cli::cli_abort("A module score has zero or non-finite standard deviation.")
  }
  as.numeric((x - mean(x)) / score_sd)
}

module_scores <- sobj[[]] |>
  tibble::rownames_to_column("cell") |>
  dplyr::transmute(
    cell,
    Mouse = as.character(Mouse),
    condition_label = as.character(.data[[condition_col]]),
    condition = dplyr::case_when(
      condition_label == control_label ~ "control",
      condition_label == estim_label ~ "estim"
    ),
    sample_id = paste0("Mouse_", Mouse, "__", condition),
    progenitor_score = standardize_score(.data[[generated_score_columns[[
      "progenitor"
    ]]]]),
    cone_bipolar_score = standardize_score(.data[[generated_score_columns[[
      "cone_bipolar"
    ]]]]),
    muller_score = standardize_score(.data[[generated_score_columns[[
      "muller"
    ]]]]),
    proliferation_score = standardize_score(.data[[generated_score_columns[[
      "proliferation"
    ]]]])
  )

if (anyNA(module_scores$condition)) {
  cli::cli_abort("Cells have a condition outside the configured comparison.")
}
if (
  anyNA(module_scores) ||
    any(!is.finite(as.matrix(module_scores[score_dimensions])))
) {
  cli::cli_abort("Module-score data contain missing or non-finite values.")
}
if (anyDuplicated(module_scores$cell)) {
  cli::cli_abort("Cell identifiers must be unique.")
}

sample_table <- module_scores |>
  dplyr::count(
    sample_id,
    Mouse,
    condition_label,
    condition,
    name = "n_cells"
  ) |>
  dplyr::arrange(condition, Mouse) |>
  dplyr::mutate(condition = factor(condition, levels = c("control", "estim")))

model_matrix <- stats::model.matrix(~condition, data = sample_table)
if (qr(model_matrix)$rank < ncol(model_matrix)) {
  cli::cli_abort("The sample-level condition design matrix is not full rank.")
}

score_matrix <- as.matrix(module_scores[score_dimensions])
rownames(score_matrix) <- module_scores$cell

# ---- neighborhoods and differential abundance ----

sce <- Seurat::as.SingleCellExperiment(sobj, assay = "RNA")
if (!identical(colnames(sce), module_scores$cell)) {
  cli::cli_abort(
    "Seurat-to-SingleCellExperiment conversion changed cell order."
  )
}

SingleCellExperiment::reducedDim(sce, score_reduction) <- score_matrix
SummarizedExperiment::colData(sce)$sample_id <- module_scores$sample_id
SummarizedExperiment::colData(sce)$condition <- module_scores$condition

set.seed(seed)
milo <- miloR::Milo(sce)
milo <- miloR::buildGraph(
  milo,
  k = k_neighbors,
  d = length(score_dimensions),
  reduced.dim = score_reduction
)
milo <- miloR::makeNhoods(
  milo,
  prop = nhood_proportion,
  k = k_neighbors,
  d = length(score_dimensions),
  refined = TRUE,
  reduced_dims = score_reduction
)
milo <- miloR::countCells(
  milo,
  meta.data = as.data.frame(SummarizedExperiment::colData(sce)),
  samples = "sample_id"
)
milo <- miloR::calcNhoodDistance(
  milo,
  d = length(score_dimensions),
  reduced.dim = score_reduction
)

nhood_counts <- miloR::nhoodCounts(milo)
sample_table_model <- sample_table |> tibble::column_to_rownames("sample_id")
# ANALYSIS_OK[sample-order]: countCells determines the required design row order.
sample_table_model <- sample_table_model[colnames(nhood_counts), , drop = FALSE]

if (anyNA(rownames(sample_table_model))) {
  cli::cli_abort("Neighborhood count columns do not match the sample design.")
}

# ANALYSIS_OK[contrast-definition]: the six Mouse × Condition samples define the primary contrast.
da_results <- miloR::testNhoods(
  milo,
  design = ~condition,
  design.df = sample_table_model,
  reduced.dim = score_reduction
)

# ---- neighborhood profiles ----

nhood_membership <- miloR::nhoods(milo)
if (!identical(rownames(nhood_membership), rownames(score_matrix))) {
  cli::cli_abort(
    "Neighborhood membership and module scores have different cells."
  )
}
if (ncol(nhood_membership) != nrow(nhood_counts)) {
  cli::cli_abort(
    "Neighborhood membership and count matrices have different sizes."
  )
}

nhood_ids <- rownames(nhood_counts)
nhood_index_cells <- colnames(nhood_membership)
nhood_sizes <- as.numeric(Matrix::colSums(nhood_membership))
score_centroids <- Matrix::crossprod(nhood_membership, score_matrix) |>
  as.matrix() |>
  sweep(1, nhood_sizes, "/")

score_pca <- stats::prcomp(score_matrix, center = TRUE, scale. = FALSE)
plot_pcs <- paste0("PC", seq_len(2L))
pca_centroids <- Matrix::crossprod(
  nhood_membership,
  score_pca$x[, plot_pcs, drop = FALSE]
) |>
  as.matrix() |>
  sweep(1, nhood_sizes, "/")

nhood_profiles <- tibble::as_tibble(score_centroids) |>
  dplyr::mutate(
    nhood = nhood_ids,
    index_cell = nhood_index_cells,
    score_pc1 = pca_centroids[, "PC1"],
    score_pc2 = pca_centroids[, "PC2"],
    n_cells = as.integer(nhood_sizes),
    .before = 1
  )

nhood_summary <- tibble::tibble(
  nhood = rownames(nhood_counts),
  mean_cells_per_sample = Matrix::rowMeans(nhood_counts),
  samples_with_cells = as.integer(Matrix::rowSums(nhood_counts > 0))
)

# ANALYSIS_OK[nhood-summary-join]: all inputs have one row per unique neighborhood.
da_table <- da_results |>
  as.data.frame() |>
  tibble::rownames_to_column("nhood") |>
  dplyr::left_join(nhood_profiles, by = "nhood") |>
  dplyr::left_join(nhood_summary, by = "nhood") |>
  dplyr::mutate(
    significant = !is.na(SpatialFDR) & SpatialFDR <= spatial_fdr_threshold
  ) |>
  dplyr::arrange(SpatialFDR, PValue)

membership_table <- Matrix::summary(nhood_membership) |>
  tibble::as_tibble() |>
  dplyr::transmute(
    nhood = nhood_ids[j],
    index_cell = nhood_index_cells[j],
    cell = rownames(nhood_membership)[i]
  )

sample_nhood_abundance <- Matrix::summary(nhood_counts) |>
  tibble::as_tibble() |>
  dplyr::transmute(
    nhood = rownames(nhood_counts)[i],
    sample_id = colnames(nhood_counts)[j],
    n_cells = x
  ) |>
  tidyr::complete(
    nhood = rownames(nhood_counts),
    sample_id = colnames(nhood_counts),
    fill = list(n_cells = 0)
  ) |>
  dplyr::left_join(
    sample_table |> dplyr::select(sample_id, sample_total = n_cells),
    by = "sample_id"
  ) |>
  dplyr::mutate(sample_fraction = n_cells / sample_total)

score_pca_loadings <- score_pca$rotation |>
  as.data.frame() |>
  tibble::rownames_to_column("score")

# ---- figures ----

plot_nhood_da <- ggplot2::ggplot(
  da_table,
  ggplot2::aes(
    x = score_pc1,
    y = score_pc2,
    color = logFC,
    size = n_cells,
    shape = significant
  )
) +
  ggplot2::geom_point(alpha = 0.85) +
  ggplot2::scale_color_gradient2(
    low = scales::muted("blue"),
    mid = "grey90",
    high = scales::muted("red"),
    midpoint = 0,
    name = "log2 fold change\n(EStim / control)"
  ) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
  ggplot2::labs(
    x = "Module-score PC1",
    y = "Module-score PC2",
    size = "Neighborhood cells",
    shape = paste0("Spatial FDR <= ", spatial_fdr_threshold)
  ) +
  theme_stone()

component_profile_data <- da_table |>
  dplyr::select(
    nhood,
    logFC,
    significant,
    n_cells,
    dplyr::all_of(score_dimensions)
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(score_dimensions),
    names_to = "module",
    values_to = "mean_score"
  )

plot_component_profiles <- ggplot2::ggplot(
  component_profile_data,
  ggplot2::aes(x = logFC, y = mean_score, size = n_cells, alpha = significant)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~module, scales = "free_y") +
  ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.25, `TRUE` = 1)) +
  ggplot2::labs(
    x = "Neighborhood log2 fold change (EStim / control)",
    y = "Mean standardized module score",
    size = "Neighborhood cells",
    alpha = paste0("Spatial FDR <= ", spatial_fdr_threshold)
  ) +
  theme_stone()

ggplot2::ggsave(
  file.path(output_dir, "nhood_da_score_space.png"),
  plot_nhood_da,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)
ggplot2::ggsave(
  file.path(output_dir, "nhood_da_score_space.pdf"),
  plot_nhood_da,
  width = 7,
  height = 5,
  units = "in"
)
ggplot2::ggsave(
  file.path(output_dir, "nhood_component_profiles.png"),
  plot_component_profiles,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)
ggplot2::ggsave(
  file.path(output_dir, "nhood_component_profiles.pdf"),
  plot_component_profiles,
  width = 8,
  height = 6,
  units = "in"
)

# ---- output tables ----

parameter_table <- tibble::tibble(
  parameter = c(
    "input_path",
    "output_dir",
    "score_reduction",
    "score_dimensions",
    "score_scaling",
    "k_neighbors",
    "nhood_proportion",
    "spatial_fdr_threshold",
    "seed",
    "design",
    "contrast",
    "normalization",
    "fdr_weighting",
    "miloR_version"
  ),
  value = c(
    input_path,
    output_dir,
    score_reduction,
    paste(score_dimensions, collapse = ","),
    "global cell-level z score per module",
    as.character(k_neighbors),
    as.character(nhood_proportion),
    as.character(spatial_fdr_threshold),
    as.character(seed),
    "~ condition",
    "estim_vs_control",
    "TMM with Mouse x Condition cell totals",
    "k-distance",
    as.character(utils::packageVersion("miloR"))
  )
)

readr::write_tsv(sample_table, file.path(output_dir, "sample_table.tsv"))
readr::write_tsv(module_scores, file.path(output_dir, "module_scores.tsv"))
readr::write_tsv(marker_features, file.path(output_dir, "marker_features.tsv"))
readr::write_tsv(da_table, file.path(output_dir, "nhood_da_results.tsv"))
readr::write_tsv(
  nhood_profiles,
  file.path(output_dir, "nhood_score_profiles.tsv")
)
readr::write_tsv(
  membership_table,
  file.path(output_dir, "nhood_cell_membership.tsv")
)
readr::write_tsv(
  sample_nhood_abundance,
  file.path(output_dir, "sample_nhood_abundance.tsv")
)
readr::write_tsv(
  score_pca_loadings,
  file.path(output_dir, "score_pca_loadings.tsv")
)
readr::write_tsv(parameter_table, file.path(output_dir, "parameters.tsv"))
