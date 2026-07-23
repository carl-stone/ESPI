#!/usr/bin/env Rscript

# Milo neighborhood differential-abundance analysis of the frozen MG-selected object.
#
# Usage: ESPI_OVERWRITE=false Rscript scripts/miloR-da.R

suppressPackageStartupMessages({
  library(here)
  here::i_am("scripts/miloR-da.R")
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

n_pcs <- 20
k_neighbors <- 50L
nhood_proportion <- 0.1

# ---- paths ----

output_dir <- file.path(config$paths$degs, "mg_selected", "miloR_da")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_paths <- file.path(
  output_dir,
  c(
    "sample_table.tsv",
    "nhood_da_results.tsv",
    "nhood_cell_membership.tsv",
    "parameters.tsv"
  )
)
assert_output_available(output_paths, config$overwrite)

# ---- input and sample design ----

sobj <- readRDS(input_path)
assert_frozen_input(input_path, sobj, config$frozen$mg)

meta <- sobj[[]] |>
  tibble::rownames_to_column("cell") |>
  dplyr::transmute(
    cell,
    Mouse = as.character(Mouse),
    condition_label = as.character(.data[[condition_col]]),
    condition = dplyr::case_when(
      condition_label == control_label ~ "control",
      condition_label == estim_label ~ "estim"
    ),
    sample_id = paste0("Mouse_", Mouse, "__", condition)
  )

if (anyNA(meta$condition)) {
  cli::cli_abort("Cells have a condition outside the configured comparison.")
}

sample_table <- meta |>
  dplyr::count(
    sample_id,
    Mouse,
    condition_label,
    condition,
    name = "n_cells"
  ) |>
  dplyr::arrange(condition, Mouse) |>
  dplyr::mutate(condition = factor(condition, levels = c("control", "estim")))

# ---- neighborhoods and differential abundance ----

sce <- Seurat::as.SingleCellExperiment(
  sobj,
  assay = SeuratObject::DefaultAssay(sobj)
)
if (!identical(colnames(sce), meta$cell)) {
  cli::cli_abort(
    "Seurat-to-SingleCellExperiment conversion changed cell order."
  )
}

SummarizedExperiment::colData(sce)$sample_id <- meta$sample_id
SummarizedExperiment::colData(sce)$condition <- meta$condition

set.seed(seed)

milo <- miloR::Milo(sce)
milo <- miloR::buildGraph(milo, k = k_neighbors, d = n_pcs, reduced.dim = "PCA")
milo <- miloR::makeNhoods(
  milo,
  prop = nhood_proportion,
  k = k_neighbors,
  d = n_pcs,
  refined = TRUE,
  reduced_dims = "PCA"
)

plotNhoodSizeHist(milo)

milo <- miloR::countCells(
  milo,
  meta.data = as.data.frame(SummarizedExperiment::colData(sce)),
  samples = "sample_id"
)
milo <- miloR::calcNhoodDistance(milo, d = n_pcs, reduced.dim = "PCA")

nhood_counts <- miloR::nhoodCounts(milo)
sample_table <- sample_table |> tibble::column_to_rownames("sample_id")
# ANALYSIS_OK[sample-order]: countCells determines the required design row order.
sample_table <- sample_table[colnames(nhood_counts), , drop = FALSE]

if (anyNA(rownames(sample_table))) {
  cli::cli_abort("Neighborhood count columns do not match the sample design.")
}

# ANALYSIS_OK[contrast-definition]: condition is the prespecified primary six-sample contrast.
da_results <- miloR::testNhoods(
  milo,
  design = ~condition,
  design.df = sample_table,
  reduced.dim = "PCA"
)

milo <- buildNhoodGraph(milo)
plotNhoodGraphDA(
  milo,
  da_results,
  alpha = 0.1,
  layout = "UMAP_PFLOG_MG_SELECTED_NO_FILTER_CC_DIMS20"
)

# ---- output tables ----

nhood_membership <- miloR::nhoods(milo)
membership_table <- Matrix::summary(nhood_membership) |>
  tibble::as_tibble() |>
  dplyr::transmute(
    nhood = colnames(nhood_membership)[j],
    cell = rownames(nhood_membership)[i]
  )

nhood_summary <- tibble::tibble(
  nhood = rownames(nhood_counts),
  n_cells = as.integer(Matrix::colSums(nhood_membership)),
  mean_cells_per_sample = Matrix::rowMeans(nhood_counts),
  samples_with_cells = as.integer(Matrix::rowSums(nhood_counts > 0))
)

# ANALYSIS_OK[nhood-summary-join]: both tables have one row per unique neighborhood.
da_table <- da_results |>
  as.data.frame() |>
  tibble::rownames_to_column("nhood") |>
  dplyr::left_join(nhood_summary, by = "nhood") |>
  dplyr::arrange(SpatialFDR, PValue)

parameter_table <- tibble::tibble(
  parameter = c(
    "input_path",
    "reduction",
    "n_pcs",
    "k_neighbors",
    "nhood_proportion",
    "seed",
    "design",
    "contrast",
    "normalization",
    "fdr_weighting",
    "miloR_version"
  ),
  value = c(
    input_path,
    "PCA",
    as.character(n_pcs),
    as.character(k_neighbors),
    as.character(nhood_proportion),
    as.character(seed),
    "~ condition",
    "estim_vs_control",
    "TMM with Mouse x Condition cell totals",
    "k-distance",
    as.character(utils::packageVersion("miloR"))
  )
)

readr::write_tsv(
  tibble::rownames_to_column(sample_table, "sample_id"),
  file.path(output_dir, "sample_table.tsv")
)
readr::write_tsv(da_table, file.path(output_dir, "nhood_da_results.tsv"))
readr::write_tsv(
  membership_table,
  file.path(output_dir, "nhood_cell_membership.tsv")
)
readr::write_tsv(parameter_table, file.path(output_dir, "parameters.tsv"))
