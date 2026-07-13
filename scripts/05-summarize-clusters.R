#!/usr/bin/env Rscript

# Generate supplemental clustering grid summaries and figures.
#
# Usage:
#   Rscript scripts/05-summarize-clusters.R
#
# Arguments:
#   None.
#
# Inputs:
#   Current clustered objects in CURRENT_OBJECT_DIR for the full clustering grid.
#
# Outputs:
#   TABLE_DIR/cluster/cluster_grid_summary.tsv
#   TABLE_DIR/cluster/cluster_grid_stability_summary.tsv
#   TABLE_DIR/cluster/cluster_grid_pairwise_stability.tsv
#   FIGURE_DIR/cluster/cluster_grid_clustree_12_panel.{png,pdf}
#   FIGURE_DIR/cluster/umap_resolution_sweep_pflog_filter_cc_dims50.{png,pdf}

suppressPackageStartupMessages({
  library(here)
})
here::i_am("scripts/05-summarize-clusters.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  stop("This script does not take command-line arguments.", call. = FALSE)
}

# ---- work ----

summary <- write_cluster_grid_summary()
stability <- write_cluster_grid_stability_tables()
splot_cluster_grid_clustree()
splot_umap_resolution_sweep()

# ---- output ----

message(
  "Wrote cluster grid summary with ",
  nrow(summary),
  " rows to ",
  file.path(TABLE_DIR, "cluster", "cluster_grid_summary.tsv")
)
message(
  "Wrote cluster stability summary with ",
  nrow(stability$summary),
  " rows and pairwise stability table with ",
  nrow(stability$pairwise),
  " rows to ",
  file.path(TABLE_DIR, "cluster")
)
message(
  "Wrote cluster supplemental figures to ",
  file.path(FIGURE_DIR, "cluster")
)
