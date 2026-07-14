#!/usr/bin/env Rscript

# Preview the analysis pipeline configuration without running any stages.
#
# Usage:
#   Rscript scripts/run-pipeline.R --dry-run
#
# Arguments:
#   --dry-run  Print the pipeline contract without creating outputs.

suppressPackageStartupMessages({
  library(here)
})
suppressMessages(here::i_am("scripts/run-pipeline.R"))
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

source_cluster_column <- "cluster_pflog_no_filter_cc_dims20_res0.3"
mg_cluster_column <- "cluster_pflog_mg_selected_no_filter_cc_dims20_res0.3"
mg_pca_dims <- 50L

# ---- validation ----

arguments <- commandArgs(trailingOnly = TRUE)
unknown_arguments <- base::setdiff(arguments, "--dry-run")

if (length(unknown_arguments) > 0) {
  stop(
    "Unknown argument(s): ",
    paste(unknown_arguments, collapse = ", "),
    call. = FALSE
  )
}
if (!("--dry-run" %in% arguments)) {
  stop(
    "Only --dry-run is supported; refusing to execute the pipeline.",
    call. = FALSE
  )
}

# ---- dry run ----

contract_lines <- c(
  "mode: dry-run",
  "input_source: counts-qc",
  "overwrite: false",
  paste0("source_cluster_column: ", source_cluster_column),
  paste0("mg_cluster_column: ", mg_cluster_column),
  paste0("mg_pca_dims: ", mg_pca_dims),
  "first_stage: process-counts",
  "final_stage: tripwires"
)

base::cat(contract_lines, sep = "\n")
