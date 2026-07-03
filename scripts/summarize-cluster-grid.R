#!/usr/bin/env Rscript

# Generate supplemental clustering grid summaries and figures.

suppressPackageStartupMessages({
  library(here)
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

summary <- write_cluster_grid_summary()
stability <- write_cluster_grid_stability_tables()
splot_cluster_grid_clustree()
splot_umap_resolution_sweep()

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
