#!/usr/bin/env Rscript

# Summarize sample-level QC metrics after cell filtering and before MG selection.
# Run with `Rscript scripts/qc_table.R`; set ESPI_OVERWRITE=true to replace output.

suppressPackageStartupMessages({
  here::i_am("scripts/qc_table.R")
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(tidyverse)
})

# ---- parameters and paths ----

config <- publication_config()
input_path <- file.path(config$paths$input_objects, "sobj_qc_filtered.rds")
output_dir <- file.path(config$paths$tables, "qc")
output_path <- file.path(output_dir, "post_filtering_qc_stats.tsv")

if (!file.exists(input_path)) {
  stop("Filtered Seurat object does not exist: ", input_path, call. = FALSE)
}
assert_output_available(output_path, config$overwrite)

# ---- summarize filtered cells ----

sobj <- readRDS(input_path)
metadata <- sobj[[]]
required_columns <- c(
  "Sample",
  "Mouse",
  "Condition",
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "pass_qc"
)
missing_columns <- setdiff(required_columns, colnames(metadata))
if (length(missing_columns) > 0L) {
  stop(
    "Filtered Seurat metadata lacks required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}
if (any(is.na(metadata$pass_qc)) || any(!metadata$pass_qc)) {
  stop(
    "Filtered Seurat object contains cells that did not pass QC.",
    call. = FALSE
  )
}
if (anyNA(metadata[, setdiff(required_columns, "pass_qc"), drop = FALSE])) {
  stop(
    "Required filtered-cell metadata contains missing values.",
    call. = FALSE
  )
}

qc_table <- metadata |>
  dplyr::group_by(Sample, Mouse, Condition) |>
  dplyr::summarize(
    n_cells = dplyr::n(),
    mean_umi_per_cell = mean(nCount_RNA),
    mean_features_per_cell = mean(nFeature_RNA),
    mean_percent_mt = mean(percent.mt),
    .groups = "drop"
  ) |>
  dplyr::arrange(Sample)

if (nrow(qc_table) != dplyr::n_distinct(metadata$Sample)) {
  stop(
    "A sample maps to more than one Mouse or Condition value.",
    call. = FALSE
  )
}

# ---- write table ----

readr::write_tsv(qc_table, output_path)
