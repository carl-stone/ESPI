#!/usr/bin/env Rscript

# Preview the analysis pipeline configuration without running any stages.
#
# Usage:
#   Rscript scripts/run-pipeline.R --dry-run \
#     [--input-source counts-qc|legacy | --input <seurat.rds>] [--overwrite]
#
# Arguments:
#   --dry-run       Print the pipeline contract without creating outputs.
#   --input-source  Use the named counts-qc or legacy source object.
#   --input         Use an explicit Seurat RDS input path.
#   --overwrite     Permit replacement of pipeline outputs when execution lands.

suppressPackageStartupMessages({
  library(here)
})
suppressMessages(here::i_am("scripts/run-pipeline.R"))
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

arguments <- commandArgs(trailingOnly = TRUE)
input_source <- "counts-qc"
input_path <- NULL
input_source_supplied <- FALSE
dry_run <- FALSE
overwrite <- FALSE

argument_index <- 1L
while (argument_index <= length(arguments)) {
  argument <- arguments[[argument_index]]

  if (identical(argument, "--dry-run")) {
    dry_run <- TRUE
    argument_index <- argument_index + 1L
    next
  }
  if (identical(argument, "--overwrite")) {
    overwrite <- TRUE
    argument_index <- argument_index + 1L
    next
  }
  if (argument %in% c("--input-source", "--input")) {
    if (
      argument_index == length(arguments) ||
        startsWith(arguments[[argument_index + 1L]], "--")
    ) {
      stop("Missing value for ", argument, ".", call. = FALSE)
    }

    value <- arguments[[argument_index + 1L]]
    if (identical(argument, "--input-source")) {
      input_source <- value
      input_source_supplied <- TRUE
    } else {
      input_path <- value
    }
    argument_index <- argument_index + 2L
    next
  }

  stop("Unknown argument: ", argument, ".", call. = FALSE)
}

if (!is.null(input_path) && input_source_supplied) {
  stop("Use either --input or --input-source, not both.", call. = FALSE)
}
if (!input_source %in% c("counts-qc", "legacy")) {
  stop(
    "--input-source must be one of counts-qc or legacy.",
    call. = FALSE
  )
}
if (!is.null(input_path)) {
  input_source <- "explicit"
} else {
  input_path <- if (identical(input_source, "counts-qc")) {
    file.path(INPUT_OBJECT_DIR, "sobj_qc_filtered.rds")
  } else {
    file.path(INPUT_OBJECT_DIR, "pipseq_processed_matrix_with_egfp.rds")
  }
}

run_spec <- list(
  input_source = input_source,
  input_path = input_path,
  overwrite = overwrite,
  normalization = "pflog",
  filter_cell_cycle_hvgs = FALSE,
  chosen_dims = 20L,
  sensitivity_dims = c(30L, 50L),
  resolution = 0.3,
  mg_pca_dims = 50L,
  counts_layer = "counts"
)

cell_cycle_tag <- if (isTRUE(run_spec$filter_cell_cycle_hvgs)) {
  "filter_cc"
} else {
  "no_filter_cc"
}
run_spec$source_branch_tag <- paste(
  run_spec$normalization,
  cell_cycle_tag,
  sep = "_"
)
run_spec$mg_branch_tag <- paste(
  run_spec$normalization,
  "mg_selected",
  cell_cycle_tag,
  sep = "_"
)
run_spec$source_cluster_column <- paste0(
  "cluster_",
  run_spec$source_branch_tag,
  "_dims",
  run_spec$chosen_dims,
  "_res",
  run_spec$resolution
)
run_spec$mg_cluster_column <- paste0(
  "cluster_",
  run_spec$mg_branch_tag,
  "_dims",
  run_spec$chosen_dims,
  "_res",
  run_spec$resolution
)

# ---- validation ----

if (!dry_run) {
  stop(
    "Only --dry-run is supported; refusing to execute the pipeline.",
    call. = FALSE
  )
}

# ---- dry run ----

contract_lines <- c(
  "mode: dry-run",
  paste0("input_source: ", run_spec$input_source),
  paste0("overwrite: ", tolower(as.character(run_spec$overwrite))),
  paste0("source_cluster_column: ", run_spec$source_cluster_column),
  paste0("mg_cluster_column: ", run_spec$mg_cluster_column),
  paste0("mg_pca_dims: ", run_spec$mg_pca_dims),
  "first_stage: process-counts",
  "final_stage: tripwires"
)

base::cat(contract_lines, sep = "\n")
