#!/usr/bin/env Rscript

# Compare saved pseudotime ranks with the configured MG-selected clustering.
# Run from the project root: source("exploratory/diffusion-pseudotime/plot-clusters.R")

suppressPackageStartupMessages({
  library(here)
  here::i_am("exploratory/diffusion-pseudotime/plot-clusters.R")
  library(devtools)
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggbeeswarm)
})
ggplot2::theme_set(theme_stone())

# ---- parameters ----

pt_cluster_config <- publication_config()
pt_cluster_column <- pt_cluster_config$selected$mg$column
pt_cluster_dir <- Sys.getenv(
  "ESPI_PSEUDOTIME_OUTPUT_DIR",
  unset = here::here(
    "exploratory",
    "diffusion-pseudotime",
    "outputs",
    "mg-minus-neuron"
  )
)

# ---- saved fit and matching cell metadata ----

pt_cluster_fit <- readRDS(file.path(pt_cluster_dir, "diffusion_pseudotime.rds"))
pt_cluster_sobj <- readRDS(pt_cluster_fit$input_path)
pt_cluster_metadata <- pt_cluster_sobj[[]]
pt_cluster_reduction <- paste0(
  "umap_",
  pt_cluster_config$selected$mg$branch,
  "_dims",
  pt_cluster_fit$settings$dimensions
)
stopifnot(
  pt_cluster_column %in% names(pt_cluster_metadata),
  pt_cluster_reduction %in% SeuratObject::Reductions(pt_cluster_sobj),
  !anyDuplicated(pt_cluster_fit$cells$cell),
  setequal(pt_cluster_fit$cells$cell, rownames(pt_cluster_metadata))
)
pt_cluster_coordinates <- SeuratObject::Embeddings(
  pt_cluster_sobj,
  pt_cluster_reduction
)
stopifnot(setequal(rownames(pt_cluster_coordinates), pt_cluster_fit$cells$cell))
pt_cluster_data <- pt_cluster_fit$cells |>
  dplyr::select(cell, pseudotime_percentile) |>
  dplyr::mutate(
    cluster_id = as.factor(pt_cluster_metadata[cell, pt_cluster_column]),
    UMAP1 = pt_cluster_coordinates[cell, 1],
    UMAP2 = pt_cluster_coordinates[cell, 2]
  )
stopifnot(
  !anyNA(pt_cluster_data),
  all(is.finite(as.matrix(pt_cluster_data[, c(
    "pseudotime_percentile",
    "UMAP1",
    "UMAP2"
  )]))),
  all(
    pt_cluster_data$pseudotime_percentile >= 0 &
      pt_cluster_data$pseudotime_percentile <= 1
  )
)

# ---- plots ----

ggsave(
  output_path(
    pt_cluster_dir,
    paste0(
      "pseudotime_percentile_by_cluster_res",
      pt_cluster_config$selected$mg$resolution,
      "_beeswarm.png"
    )
  ),
  plot = ggplot(pt_cluster_data, aes(cluster_id, pseudotime_percentile)) +
    geom_beeswarm(size = 0.7, cex = 0.7) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::label_percent()
    ) +
    labs(x = "Cluster ID", y = "Pseudotime percentile"),
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  output_path(pt_cluster_dir, "umap_pseudotime_percentile.png"),
  plot = ggplot(
    pt_cluster_data,
    aes(UMAP1, UMAP2, color = pseudotime_percentile)
  ) +
    geom_point(size = 0.8) +
    scale_color_viridis_c(limits = c(0, 1), labels = scales::label_percent()) +
    coord_equal() +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Pseudotime\nrank"),
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  output_path(pt_cluster_dir, "umap_clusters.png"),
  plot = ggplot(pt_cluster_data, aes(UMAP1, UMAP2, color = cluster_id)) +
    geom_point(size = 0.8) +
    scale_color_viridis_d() +
    coord_equal() +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Cluster ID"),
  width = 8,
  height = 6,
  dpi = 300
)
