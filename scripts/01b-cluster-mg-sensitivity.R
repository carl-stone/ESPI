#!/usr/bin/env Rscript

# Rebuild MG-selected clustering sensitivity grids from frozen preprocessing.
# Loads both MG preprocessing branches written by 01-regenerate-frozen.R,
# runs the fixed Leiden/UMAP grid, and writes selected frozen MG objects.

suppressPackageStartupMessages({
  here::i_am("scripts/01b-cluster-mg-sensitivity.R")
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(tidyverse)
})

# ---- parameters and inputs ----

config <- publication_config()
seed <- config$seed
selected_seed <- config$selected$mg$seed
current_object_dir <- config$paths$current_objects
figure_dir <- config$paths$figures
table_dir <- config$paths$tables
cluster_figure_dir <- file.path(figure_dir, "cluster")
mg_figure_dir <- file.path(figure_dir, "mg_selected")
mg_table_dir <- file.path(table_dir, "mg_selected")
dims_grid <- config$grid$dimensions
resolutions <- config$grid$resolutions
small_cluster_threshold <- 50L
mg_branch_tags <- c(
  "pflog_mg_selected_no_filter_cc",
  "pflog_mg_selected_filter_cc"
)
mg_preprocess_tags <- c(
  "pflog_mg_selected_no-filter-cc",
  "pflog_mg_selected_filter-cc"
)
mg_preprocess_paths <- file.path(
  current_object_dir,
  paste0("preprocess_", mg_preprocess_tags, ".rds")
)

if (
  !dir.exists(current_object_dir) || file.access(current_object_dir, 2L) != 0L
) {
  stop(
    "Refusing MG clustering regeneration; object directory is not writable: ",
    current_object_dir,
    call. = FALSE
  )
}
if (any(!file.exists(mg_preprocess_paths))) {
  stop(
    "MG clustering regeneration requires both preprocessing objects: ",
    paste(mg_preprocess_paths, collapse = ", "),
    call. = FALSE
  )
}


# ---- load frozen preprocessing branches ----

mg_preprocessed <- mg_preprocess_paths |>
  purrr::map(readRDS) |>
  stats::setNames(mg_branch_tags)
# ---- both MG Leiden/UMAP grids and summaries ----

mg_branches <- data.frame(
  filtered_cell_cycle = c(FALSE, TRUE),
  branch_tag = mg_branch_tags,
  stringsAsFactors = FALSE
)
mg_clustered <- list()
for (branch_index in seq_len(nrow(mg_branches))) {
  branch_info <- mg_branches[branch_index, ]
  branch_sobj <- mg_preprocessed[[branch_info$branch_tag]]
  candidate_names <- character()
  candidate_seeds <- integer()
  umap_seeds <- integer()
  is_selected_branch <- identical(
    branch_info$branch_tag,
    config$selected$mg$branch
  )
  for (dims in dims_grid) {
    old_idents <- SeuratObject::Idents(branch_sobj)
    branch_sobj <- Seurat::FindNeighbors(
      branch_sobj,
      reduction = "pca",
      dims = 1:dims
    )
    for (resolution in resolutions) {
      column <- paste0(
        "cluster_",
        branch_info$branch_tag,
        "_dims",
        dims,
        "_res",
        resolution
      )
      candidate_seed <- if (
        is_selected_branch &&
          dims == config$selected$mg$dimensions &&
          resolution == config$selected$mg$resolution
      ) {
        selected_seed
      } else {
        seed
      }
      branch_sobj <- Seurat::FindClusters(
        branch_sobj,
        algorithm = 4,
        leiden_method = "igraph",
        resolution = resolution,
        random.seed = candidate_seed
      )
      branch_sobj@meta.data[[column]] <- SeuratObject::Idents(branch_sobj)
      SeuratObject::Idents(branch_sobj) <- old_idents
      candidate_names <- c(candidate_names, column)
      candidate_seeds[[column]] <- candidate_seed
    }
    reduction_name <- sprintf("umap_%s_dims%d", branch_info$branch_tag, dims)
    umap_seed <- if (
      is_selected_branch && dims == config$selected$mg$dimensions
    ) {
      selected_seed
    } else {
      seed
    }
    branch_sobj <- Seurat::RunUMAP(
      branch_sobj,
      reduction = "pca",
      dims = 1:dims,
      reduction.name = reduction_name,
      reduction.key = paste0(gsub("[^A-Za-z0-9]", "", reduction_name), "_"),
      seed.use = umap_seed
    )
    umap_seeds[[reduction_name]] <- umap_seed
    for (resolution in resolutions) {
      column <- paste0(
        "cluster_",
        branch_info$branch_tag,
        "_dims",
        dims,
        "_res",
        resolution
      )
      plot <- Seurat::DimPlot(
        branch_sobj,
        reduction = reduction_name,
        group.by = column,
        label = TRUE,
        pt.size = 0.25
      )
      png_path <- file.path(
        cluster_figure_dir,
        paste0(reduction_name, "_by_", column, ".png")
      )
      ggplot2::ggsave(output_path(png_path), plot, width = 5, height = 5)
      ggplot2::ggsave(
        output_path(sub("\\.png$", ".pdf", png_path)),
        plot,
        width = 5,
        height = 5
      )
      if (
        column %in%
          c(config$selected$mg$column, config$selected$mg_filter_cc$column)
      ) {
        copy_notebook_figure(png_path)
      }
    }
    prefix <- sprintf("cluster_%s_dims%d_res", branch_info$branch_tag, dims)
    cluster_data <- branch_sobj@meta.data[,
      startsWith(colnames(branch_sobj@meta.data), prefix),
      drop = FALSE
    ]
    clustree_plot <- clustree::clustree(cluster_data, prefix = prefix) +
      ggplot2::guides(edge_colour = "none")
    clustree_png <- file.path(
      cluster_figure_dir,
      sprintf("clustree_%s_dims%d.png", branch_info$branch_tag, dims)
    )
    ggplot2::ggsave(
      output_path(clustree_png),
      clustree_plot,
      width = 6,
      height = 6
    )
    ggplot2::ggsave(
      output_path(sub("\\.png$", ".pdf", clustree_png)),
      clustree_plot,
      width = 6,
      height = 6
    )
  }
  branch_sobj@misc$clustering <- list(
    algorithm = "leiden",
    filtered_cell_cycle = branch_info$filtered_cell_cycle,
    branch_tag = branch_info$branch_tag,
    resolutions = resolutions,
    dims_grid = dims_grid,
    elbow_n = 20L,
    candidate_names = candidate_names,
    candidate_seeds = candidate_seeds,
    umap_seeds = umap_seeds
  )
  saveRDS(
    branch_sobj,
    output_path(
      current_object_dir,
      paste0("cluster_", branch_info$branch_tag, "_elbow20.rds")
    )
  )
  mg_clustered[[branch_info$branch_tag]] <- branch_sobj
}
mg_summary_rows <- list()
summary_index <- 1L
for (branch_index in seq_len(nrow(mg_branches))) {
  branch_info <- mg_branches[branch_index, ]
  sobj <- mg_clustered[[branch_info$branch_tag]]
  candidates <- tidyr::crossing(dims = dims_grid, resolution = resolutions) |>
    dplyr::mutate(
      cluster_column = paste0(
        "cluster_",
        branch_info$branch_tag,
        "_dims",
        dims,
        "_res",
        resolution
      )
    )
  for (candidate_index in seq_len(nrow(candidates))) {
    candidate <- candidates[candidate_index, ]
    labels <- sobj@meta.data[[candidate$cluster_column]]
    cluster_sizes <- table(as.character(labels))
    small_clusters <- cluster_sizes[cluster_sizes < small_cluster_threshold]
    mg_summary_rows[[summary_index]] <- data.frame(
      branch_tag = branch_info$branch_tag,
      filtered_cell_cycle = branch_info$filtered_cell_cycle,
      dims = candidate$dims,
      resolution = candidate$resolution,
      cluster_column = candidate$cluster_column,
      n_cells = length(labels),
      n_clusters = length(cluster_sizes),
      min_cluster_size = as.integer(min(cluster_sizes)),
      q25_cluster_size = as.numeric(stats::quantile(
        cluster_sizes,
        0.25,
        names = FALSE
      )),
      median_cluster_size = as.numeric(stats::median(cluster_sizes)),
      max_cluster_size = as.integer(max(cluster_sizes)),
      n_clusters_lt50 = length(small_clusters),
      fraction_clusters_lt50 = length(small_clusters) / length(cluster_sizes),
      stringsAsFactors = FALSE
    )
    summary_index <- summary_index + 1L
  }
}
mg_summary <- do.call(rbind, mg_summary_rows)
mg_summary_ordered <- mg_summary[
  order(mg_summary$branch_tag, mg_summary$dims, mg_summary$resolution),
]
utils::write.table(
  mg_summary_ordered,
  output_path(mg_table_dir, "mg_selected_cluster_grid_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
for (branch_index in seq_len(nrow(mg_branches))) {
  branch_info <- mg_branches[branch_index, ]
  sobj <- mg_clustered[[branch_info$branch_tag]]
  candidates <- tidyr::crossing(dims = dims_grid, resolution = resolutions) |>
    dplyr::mutate(
      cluster_column = paste0(
        "cluster_",
        branch_info$branch_tag,
        "_dims",
        dims,
        "_res",
        resolution
      )
    )
  for (dims in sort(unique(candidates$dims))) {
    reduction <- sprintf("umap_%s_dims%d", branch_info$branch_tag, dims)
    dim_candidates <- candidates |>
      dplyr::filter(dims == .env$dims) |>
      dplyr::arrange(resolution)
    plots <- lapply(seq_len(nrow(dim_candidates)), function(index) {
      candidate <- dim_candidates[index, ]
      Seurat::DimPlot(
        sobj,
        reduction = reduction,
        group.by = candidate$cluster_column,
        label = TRUE,
        pt.size = 0.25
      ) +
        ggplot2::ggtitle(sprintf(
          "res %s",
          format(candidate$resolution, trim = TRUE, scientific = FALSE)
        )) +
        ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
        ggplot2::theme(
          legend.position = "none",
          plot.title = ggplot2::element_text(size = 11)
        )
    })
    sweep <- patchwork::wrap_plots(plots, ncol = nrow(dim_candidates)) +
      patchwork::plot_annotation(
        title = sprintf("%s; %d PCs", branch_info$branch_tag, dims)
      )
    stem <- paste0(
      "mg_selected_umap_resolution_sweep_",
      branch_info$branch_tag,
      "_dims",
      dims
    )
    png_path <- file.path(mg_figure_dir, paste0(stem, ".png"))
    ggplot2::ggsave(
      output_path(png_path),
      sweep,
      width = 3.6 * nrow(dim_candidates),
      height = 4.2
    )
    ggplot2::ggsave(
      output_path(mg_figure_dir, paste0(stem, ".pdf")),
      sweep,
      width = 3.6 * nrow(dim_candidates),
      height = 4.2
    )
    if (dims == config$selected$mg$dimensions) copy_notebook_figure(png_path)
  }
}

message(
  "MG clustering complete. Selected column: ",
  config$selected$mg$column,
  "; seed: ",
  selected_seed
)
