#!/usr/bin/env Rscript

# Summarize and plot the mg-selected PFlog clustering grid.
#
# Usage:
#   Rscript scripts/summarize-mg-selected-grid.R [--elbow-n <positive integer>]
#
# Inputs:
#   CURRENT_OBJECT_DIR/cluster_pflog_mg_selected_no_filter_cc_elbow<N>.rds
#   CURRENT_OBJECT_DIR/cluster_pflog_mg_selected_filter_cc_elbow<N>.rds
#
# Outputs:
#   TABLE_DIR/mg_selected/mg_selected_cluster_grid_summary.tsv
#   FIGURE_DIR/mg_selected/mg_selected_umap_resolution_sweep_<branch>_dims<dims>.(png|pdf)
#   notebook/figures/<sweep png filename> symlinks.

suppressPackageStartupMessages({
  library(here)
})
here::i_am("scripts/summarize-mg-selected-grid.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

get_arg <- function(args, flag, default) {
  match_index <- match(flag, args)
  if (is.na(match_index)) {
    return(default)
  }
  if (
    match_index == length(args) || startsWith(args[[match_index + 1]], "--")
  ) {
    stop("Missing value for ", flag, call. = FALSE)
  }
  args[[match_index + 1]]
}

cli_args <- commandArgs(trailingOnly = TRUE)
allowed_flags <- c("--elbow-n")
unknown_flags <- cli_args[
  startsWith(cli_args, "--") & !cli_args %in% allowed_flags
]
if (length(unknown_flags) > 0) {
  stop(
    "Unknown argument(s): ",
    paste(unknown_flags, collapse = ", "),
    call. = FALSE
  )
}

elbow_n <- as.integer(get_arg(cli_args, "--elbow-n", "20"))
if (length(elbow_n) != 1L || is.na(elbow_n) || elbow_n <= 0) {
  stop("--elbow-n must be a positive integer.", call. = FALSE)
}

MG_SELECTED_BRANCHES <- data.frame(
  branch_tag = c(
    "pflog_mg_selected_no_filter_cc",
    "pflog_mg_selected_filter_cc"
  ),
  filtered_cell_cycle = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)
MG_SELECTED_SMALL_CLUSTER_THRESHOLD <- 50L
MG_SELECTED_SWEEP_PANEL_WIDTH <- 3.6
MG_SELECTED_SWEEP_HEIGHT <- 4.2

# ---- helpers ----

mg_selected_object_path <- function(branch_tag, elbow_n) {
  file.path(
    CURRENT_OBJECT_DIR,
    sprintf("cluster_%s_elbow%d.rds", branch_tag, elbow_n)
  )
}

mg_selected_res_tag <- function(resolution) {
  format(resolution, trim = TRUE, scientific = FALSE)
}

parse_mg_selected_cluster_column <- function(column, branch_tag) {
  pattern <- sprintf("^cluster_%s_dims([0-9]+)_res(.+)$", branch_tag)
  hit <- regexec(pattern, column, perl = TRUE)
  parts <- regmatches(column, hit)[[1]]
  if (length(parts) != 3L) {
    stop("Cannot parse mg-selected cluster column: ", column, call. = FALSE)
  }
  resolution <- suppressWarnings(as.numeric(parts[[3L]]))
  if (is.na(resolution) || !is.finite(resolution) || resolution <= 0) {
    stop("Invalid resolution in cluster column: ", column, call. = FALSE)
  }
  data.frame(
    cluster_column = column,
    dims = as.integer(parts[[2L]]),
    resolution = resolution,
    stringsAsFactors = FALSE
  )
}

mg_selected_candidate_columns <- function(sobj, branch_tag) {
  candidate_names <- sobj@misc$clustering$candidate_names
  if (!is.character(candidate_names) || length(candidate_names) == 0L) {
    stop(
      "Missing sobj@misc$clustering$candidate_names for branch: ",
      branch_tag,
      call. = FALSE
    )
  }
  if (anyNA(candidate_names) || any(!nzchar(candidate_names))) {
    stop(
      "Invalid empty candidate cluster column name for branch: ",
      branch_tag,
      call. = FALSE
    )
  }
  missing_columns <- setdiff(candidate_names, colnames(sobj@meta.data))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing candidate cluster column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  parsed <- lapply(
    candidate_names,
    parse_mg_selected_cluster_column,
    branch_tag = branch_tag
  )
  candidates <- do.call(rbind, parsed)
  candidates <- candidates[order(candidates$dims, candidates$resolution), ]
  rownames(candidates) <- NULL
  candidates
}

mg_selected_size_summary <- function(labels) {
  if (anyNA(labels)) {
    stop("Cluster labels contain NA values.", call. = FALSE)
  }
  cluster_sizes <- table(as.character(labels))
  small_clusters <- cluster_sizes[
    cluster_sizes < MG_SELECTED_SMALL_CLUSTER_THRESHOLD
  ]
  data.frame(
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
}

load_mg_selected_objects <- function(branches, elbow_n) {
  objects <- vector("list", nrow(branches))
  names(objects) <- branches$branch_tag
  for (i in seq_len(nrow(branches))) {
    branch_tag <- branches$branch_tag[[i]]
    object_path <- mg_selected_object_path(branch_tag, elbow_n)
    if (!file.exists(object_path)) {
      stop("Missing clustered mg-selected object: ", object_path, call. = FALSE)
    }
    sobj <- readRDS(object_path)
    object_branch_tag <- sobj@misc$clustering$branch_tag
    if (
      is.character(object_branch_tag) &&
        length(object_branch_tag) == 1L &&
        nzchar(object_branch_tag) &&
        !identical(object_branch_tag, branch_tag)
    ) {
      stop(
        "Object branch tag mismatch for ",
        object_path,
        ": expected ",
        branch_tag,
        ", found ",
        object_branch_tag,
        call. = FALSE
      )
    }
    objects[[branch_tag]] <- sobj
  }
  objects
}

collect_mg_selected_summary <- function(objects, branches) {
  rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(branches))) {
    branch_info <- branches[i, ]
    sobj <- objects[[branch_info$branch_tag]]
    candidates <- mg_selected_candidate_columns(sobj, branch_info$branch_tag)
    for (j in seq_len(nrow(candidates))) {
      candidate <- candidates[j, ]
      labels <- sobj@meta.data[[candidate$cluster_column]]
      rows[[idx]] <- cbind(
        data.frame(
          branch_tag = branch_info$branch_tag,
          filtered_cell_cycle = branch_info$filtered_cell_cycle,
          dims = candidate$dims,
          resolution = candidate$resolution,
          cluster_column = candidate$cluster_column,
          stringsAsFactors = FALSE
        ),
        mg_selected_size_summary(labels)
      )
      idx <- idx + 1L
    }
  }
  summary <- do.call(rbind, rows)
  summary <- summary[
    order(
      summary$branch_tag,
      summary$dims,
      summary$resolution
    ),
  ]
  rownames(summary) <- NULL
  summary
}

link_notebook_png <- function(png_path) {
  notebook_figure_dir <- here::here("notebook", "figures")
  dir.create(notebook_figure_dir, recursive = TRUE, showWarnings = FALSE)
  notebook_png_path <- file.path(notebook_figure_dir, basename(png_path))
  if (
    file.exists(notebook_png_path) || nzchar(Sys.readlink(notebook_png_path))
  ) {
    unlink(notebook_png_path)
  }
  link_created <- file.symlink(png_path, notebook_png_path)
  if (!isTRUE(link_created)) {
    stop("Failed to link notebook figure: ", notebook_png_path, call. = FALSE)
  }
  notebook_png_path
}

save_mg_selected_resolution_sweeps <- function(objects, branches) {
  out_dir <- file.path(FIGURE_DIR, "mg_selected")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  linked_paths <- character()

  for (i in seq_len(nrow(branches))) {
    branch_info <- branches[i, ]
    branch_tag <- branch_info$branch_tag
    sobj <- objects[[branch_tag]]
    candidates <- mg_selected_candidate_columns(sobj, branch_tag)
    for (dims in sort(unique(candidates$dims))) {
      reduction <- sprintf("umap_%s_dims%d", branch_tag, dims)
      if (!reduction %in% names(sobj@reductions)) {
        stop("Missing UMAP reduction: ", reduction, call. = FALSE)
      }
      dim_candidates <- candidates[candidates$dims == dims, ]
      dim_candidates <- dim_candidates[order(dim_candidates$resolution), ]
      plots <- lapply(seq_len(nrow(dim_candidates)), function(row_idx) {
        candidate <- dim_candidates[row_idx, ]
        column <- candidate$cluster_column
        if (!column %in% colnames(sobj@meta.data)) {
          stop("Missing cluster metadata column: ", column, call. = FALSE)
        }
        Seurat::DimPlot(
          sobj,
          reduction = reduction,
          group.by = column,
          label = TRUE
        ) +
          ggplot2::ggtitle(sprintf(
            "res %s",
            mg_selected_res_tag(candidate$resolution)
          )) +
          ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
          ggplot2::theme(
            legend.position = "none",
            plot.title = ggplot2::element_text(size = 11)
          )
      })
      plot <- patchwork::wrap_plots(plots, ncol = nrow(dim_candidates)) +
        patchwork::plot_annotation(
          title = sprintf("%s; %d PCs", branch_tag, dims)
        )
      out_tag <- sprintf(
        "mg_selected_umap_resolution_sweep_%s_dims%d",
        branch_tag,
        dims
      )
      png_path <- file.path(out_dir, sprintf("%s.png", out_tag))
      pdf_path <- file.path(out_dir, sprintf("%s.pdf", out_tag))
      ggplot2::ggsave(
        png_path,
        plot,
        width = MG_SELECTED_SWEEP_PANEL_WIDTH * nrow(dim_candidates),
        height = MG_SELECTED_SWEEP_HEIGHT
      )
      ggplot2::ggsave(
        pdf_path,
        plot,
        width = MG_SELECTED_SWEEP_PANEL_WIDTH * nrow(dim_candidates),
        height = MG_SELECTED_SWEEP_HEIGHT
      )
      linked_paths <- c(linked_paths, link_notebook_png(png_path))
    }
  }

  linked_paths
}

# ---- work ----

objects <- load_mg_selected_objects(MG_SELECTED_BRANCHES, elbow_n = elbow_n)
summary <- collect_mg_selected_summary(objects, MG_SELECTED_BRANCHES)

# ---- output ----

summary_path <- file.path(
  TABLE_DIR,
  "mg_selected",
  "mg_selected_cluster_grid_summary.tsv"
)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  summary,
  summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
linked_paths <- save_mg_selected_resolution_sweeps(
  objects,
  MG_SELECTED_BRANCHES
)

message(
  "Wrote mg-selected cluster grid summary with ",
  nrow(summary),
  " rows to ",
  summary_path
)
message(
  "Wrote mg-selected UMAP resolution sweeps to ",
  file.path(FIGURE_DIR, "mg_selected")
)
message("Linked notebook figure(s): ", paste(linked_paths, collapse = ", "))
