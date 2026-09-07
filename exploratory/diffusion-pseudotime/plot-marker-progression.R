#!/usr/bin/env Rscript

# Replot saved MG-selected pseudotime without refitting or filtering cells.
# Each gene's y-axis runs from zero to its maximum bin mean, not its cell maximum.
# Run from the project root: source("exploratory/diffusion-pseudotime/plot-marker-progression.R")

suppressPackageStartupMessages({
  library(here)
  here::i_am("exploratory/diffusion-pseudotime/plot-marker-progression.R")
  library(devtools)
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
ggplot2::theme_set(theme_stone())

# ---- parameters ----

pt_plot_bins <- 100L
pt_plot_scale_factor <- 10000
pt_plot_genes <- c("Rlbp1", "Hes1", "Hes6", "Ascl1", "Otx2")
pt_plot_dir <- Sys.getenv(
  "ESPI_PSEUDOTIME_OUTPUT_DIR",
  unset = here::here(
    "exploratory",
    "diffusion-pseudotime",
    "outputs",
    "mg-minus-neuron"
  )
)

# ---- saved ordering and normalized counts ----

pt_plot_fit <- readRDS(file.path(pt_plot_dir, "diffusion_pseudotime.rds"))
pt_plot_sobj <- readRDS(pt_plot_fit$input_path)
pt_plot_counts <- SeuratObject::LayerData(
  pt_plot_sobj,
  assay = "RNA",
  layer = "counts"
)
stopifnot(
  identical(colnames(pt_plot_counts), pt_plot_fit$cells$cell),
  !anyDuplicated(pt_plot_fit$cells$cell),
  all(pt_plot_genes %in% rownames(pt_plot_counts)),
  all(is.finite(pt_plot_fit$cells$pseudotime)),
  all(is.finite(pt_plot_counts@x)),
  all(pt_plot_counts@x >= 0),
  nrow(pt_plot_fit$cells) >= pt_plot_bins
)
pt_plot_library_sizes <- Matrix::colSums(pt_plot_counts)
stopifnot(all(is.finite(pt_plot_library_sizes)), all(pt_plot_library_sizes > 0))

# Normalize each cell by its total RNA count before averaging; retain zeros.
pt_plot_normalized <- sweep(
  as.matrix(pt_plot_counts[pt_plot_genes, , drop = FALSE]),
  2L,
  pt_plot_library_sizes,
  "/"
) *
  pt_plot_scale_factor
pt_marker_trends_100 <- as.data.frame(t(pt_plot_normalized)) |>
  tibble::rownames_to_column("cell") |>
  dplyr::left_join(
    dplyr::select(pt_plot_fit$cells, cell, pseudotime, pseudotime_percentile),
    by = "cell"
  ) |>
  dplyr::mutate(bin = dplyr::ntile(pseudotime, pt_plot_bins)) |>
  tidyr::pivot_longer(
    dplyr::all_of(pt_plot_genes),
    names_to = "gene",
    values_to = "normalized_counts"
  ) |>
  dplyr::group_by(gene, bin) |>
  dplyr::summarise(
    mean_normalized_counts = mean(normalized_counts),
    fraction_detected = mean(normalized_counts > 0),
    n_cells = dplyr::n(),
    pseudotime = median(pseudotime),
    pseudotime_percentile = median(pseudotime_percentile),
    .groups = "drop"
  ) |>
  dplyr::mutate(gene = factor(gene, levels = pt_plot_genes))
stopifnot(
  nrow(pt_marker_trends_100) == length(pt_plot_genes) * pt_plot_bins,
  all(is.finite(pt_marker_trends_100$mean_normalized_counts)),
  all(pt_marker_trends_100$mean_normalized_counts >= 0)
)

# ---- plot and exports ----

print(
  ggplot(
    pt_marker_trends_100,
    aes(pseudotime_percentile, mean_normalized_counts)
  ) +
    geom_point(size = 1.2) +
    facet_wrap(~gene, scales = "free_y", ncol = 1) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = 0)) +
    coord_cartesian(clip = "off") +
    labs(x = "Pseudotime percentile", y = "Mean counts per 10,000")
)
ggsave(
  output_path(pt_plot_dir, "marker_progression_100bins_normalized_counts.png"),
  plot = ggplot2::last_plot(),
  width = 7,
  height = 9,
  dpi = 200
)
readr::write_csv(
  pt_marker_trends_100,
  output_path(pt_plot_dir, "marker_trends_100bins_normalized_counts.csv")
)
