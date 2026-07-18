#!/usr/bin/env Rscript

# Regenerate all frozen preprocessing, clustering, and MG-selection artifacts.
# The top-to-bottom sections cover counts, QC, four source preprocessing
# branches, fixed source and MG clustering grids, MG selection, summaries,
# supplemental plots, and frozen writes.
#
# This deliberate maintenance script has no command-line parser or stage
# discovery. It uses the data root resolved by R/config.R and refuses to run
# unless both frozen-object directories already exist and are writable.

suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- fixed parameters and immediate writable-root refusal ----

config <- publication_config()
seed <- config$seed
project_root <- config$paths$project
data_root <- config$paths$data
object_dir <- config$paths$objects
input_object_dir <- config$paths$input_objects
current_object_dir <- config$paths$current_objects
figure_dir <- config$paths$figures
table_dir <- config$paths$tables
notebook_figure_dir <- config$paths$notebook_figures

if (
  !dir.exists(input_object_dir) ||
    !dir.exists(current_object_dir) ||
    file.access(input_object_dir, 2L) != 0L ||
    file.access(current_object_dir, 2L) != 0L
) {
  stop(
    paste0(
      "Refusing frozen regeneration: seurat_objects/input and ",
      "seurat_objects/current must both exist and be writable."
    ),
    call. = FALSE
  )
}

raw_counts_dir <- file.path(data_root, "data", "input", "Raw Matrices")
metadata_path <- file.path(raw_counts_dir, "Sample_Metadata_MS1.txt")
raw_object_path <- file.path(data_root, "data", "input", "sobj_raw.rds")
qc_input_dir <- input_object_dir
qc_figure_dir <- file.path(figure_dir, "qc")
qc_table_dir <- file.path(table_dir, "qc")
preprocess_figure_dir <- file.path(figure_dir, "preprocess")
cluster_figure_dir <- file.path(figure_dir, "cluster")
cluster_table_dir <- file.path(table_dir, "cluster")
mg_figure_dir <- file.path(figure_dir, "mg_selected")
mg_table_dir <- file.path(table_dir, "mg_selected")

normalizations <- c("log1p", "pflog")
filter_states <- c(FALSE, TRUE)
dims_grid <- c(20L, 30L, 50L)
resolutions <- c(0.3, 0.5, 0.8)
source_branch_tags <- c(
  "log1p_no_filter_cc",
  "log1p_filter_cc",
  "pflog_no_filter_cc",
  "pflog_filter_cc"
)
source_preprocess_tags <- c(
  "log1p_no-filter-cc",
  "log1p_filter-cc",
  "pflog_no-filter-cc",
  "pflog_filter-cc"
)
mg_branch_tags <- c(
  "pflog_mg_selected_no_filter_cc",
  "pflog_mg_selected_filter_cc"
)
mg_preprocess_tags <- c(
  "pflog_mg_selected_no-filter-cc",
  "pflog_mg_selected_filter-cc"
)
# ANALYSIS_OK[R026]: local path helper is called by this executable script.
cluster_column <- function(branch, dims, resolution) {
  sprintf(
    "cluster_%s_dims%d_res%s",
    branch,
    dims,
    format(resolution, trim = TRUE, scientific = FALSE)
  )
}
# ANALYSIS_OK[R026]: local cluster-object path helper is called by this executable script.
cluster_object_path <- function(branch) {
  file.path(current_object_dir, sprintf("cluster_%s_elbow20.rds", branch))
}
# ANALYSIS_OK[R026]: local preprocessing path helper is called by this executable script.
preprocess_object_path <- function(tag) {
  file.path(current_object_dir, sprintf("preprocess_%s.rds", tag))
}

# Every primary path is fixed before any input read. Notebook mirrors are
# intentionally excluded: they are existing regular files replaced only after
# the corresponding primary PNG is complete.
phase_output_paths <- c(
  raw_object_path,
  file.path(qc_input_dir, c("sobj_raw_with_qc.rds", "sobj_qc_filtered.rds")),
  file.path(
    qc_table_dir,
    c(
      "sample_cell_call_summary.tsv",
      "sample_qc_summary.tsv",
      "sample_qc_mad_thresholds.tsv",
      "sobj_qc_summary_by_sample.tsv",
      "sample_qc_attrition.tsv"
    )
  ),
  file.path(
    qc_figure_dir,
    c(
      "knee_plot.png",
      "cell_call_plot.png",
      "sample_cell_call_fraction.png",
      "sample_qc_summary.png",
      "count_sample_qc_plot.png",
      "feature_sample_qc_plot.png",
      "mt_sample_qc_plot.png",
      "count_feature_mt_sample_scatter.png",
      "feature_vs_mt.png"
    )
  ),
  file.path(
    cluster_table_dir,
    c(
      "cluster_grid_summary.tsv",
      "cluster_grid_stability_summary.tsv",
      "cluster_grid_pairwise_stability.tsv"
    )
  ),
  file.path(
    cluster_figure_dir,
    c(
      "cluster_grid_clustree_12_panel.png",
      "cluster_grid_clustree_12_panel.pdf",
      "umap_resolution_sweep_pflog_filter_cc_dims50.png",
      "umap_resolution_sweep_pflog_filter_cc_dims50.pdf"
    )
  ),
  file.path(mg_table_dir, "mg_selected_cluster_selection.tsv"),
  file.path(
    mg_figure_dir,
    c(
      "mg_selected_cluster_selection_diagnostics.png",
      "mg_selected_cluster_selection_diagnostics.pdf",
      "elbow_pflog_mg_selected_no_filter_cc.png",
      "elbow_pflog_mg_selected_no_filter_cc.pdf",
      "elbow_pflog_mg_selected_filter_cc.png",
      "elbow_pflog_mg_selected_filter_cc.pdf"
    )
  ),
  file.path(mg_table_dir, "mg_selected_cluster_grid_summary.tsv"),
  file.path(table_dir, "frozen_object_numbers.tsv")
)
for (tag in source_preprocess_tags) {
  phase_output_paths <- c(
    phase_output_paths,
    preprocess_object_path(tag),
    file.path(
      preprocess_figure_dir,
      paste0(
        c("qc_metrics_violin_", "hvg_scatter_", "dim_heatmap_", "elbow_"),
        tag,
        ".png"
      )
    ),
    file.path(
      preprocess_figure_dir,
      paste0(
        c("qc_metrics_violin_", "hvg_scatter_", "dim_heatmap_", "elbow_"),
        tag,
        ".pdf"
      )
    )
  )
}
for (tag in mg_preprocess_tags) {
  phase_output_paths <- c(
    phase_output_paths,
    preprocess_object_path(tag)
  )
}
for (branch in c(source_branch_tags, mg_branch_tags)) {
  phase_output_paths <- c(phase_output_paths, cluster_object_path(branch))
  for (dims in dims_grid) {
    phase_output_paths <- c(
      phase_output_paths,
      file.path(
        cluster_figure_dir,
        paste0("clustree_", branch, "_dims", dims, c(".png", ".pdf"))
      )
    )
    for (resolution in resolutions) {
      column <- cluster_column(branch, dims, resolution)
      umap <- sprintf("umap_%s_dims%d", branch, dims)
      phase_output_paths <- c(
        phase_output_paths,
        file.path(
          cluster_figure_dir,
          paste0(umap, "_by_", column, c(".png", ".pdf"))
        )
      )
    }
  }
}
for (branch in mg_branch_tags) {
  for (dims in dims_grid) {
    phase_output_paths <- c(
      phase_output_paths,
      file.path(
        mg_figure_dir,
        paste0(
          "mg_selected_umap_resolution_sweep_",
          branch,
          "_dims",
          dims,
          c(".png", ".pdf")
        )
      )
    )
  }
}
phase_output_paths <- unique(phase_output_paths)
assert_output_available(phase_output_paths, config$overwrite)
dir.create(preprocess_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cluster_figure_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(preprocess_figure_dir) || !dir.exists(cluster_figure_dir)) {
  stop("Failed to create phase-01 figure directories.", call. = FALSE)
}

notebook_preprocess_pngs <- c(
  "elbow_log1p_no-filter-cc.png",
  "elbow_pflog_no-filter-cc.png",
  "dim_heatmap_log1p_no-filter-cc.png",
  "dim_heatmap_log1p_filter-cc.png",
  "dim_heatmap_pflog_no-filter-cc.png",
  "dim_heatmap_pflog_filter-cc.png",
  "hvg_scatter_log1p_no-filter-cc.png"
)

# Safe runtime notebook mirror: reject symlinks, copy to a regular temporary
# sibling, hash-verify, replace only an existing regular destination. POSIX
# rename-over is attempted first; the fallback moves the old destination aside
# and restores it if installation fails.
# ANALYSIS_OK[R026]: local notebook mirror helper is called by this executable script.
mirror_notebook_png <- function(source, destination) {
  if (!file.exists(source) || !isTRUE(file.info(source)$isdir == FALSE)) {
    stop(
      "Notebook mirror source is not a regular file: ",
      source,
      call. = FALSE
    )
  }
  if (nzchar(Sys.readlink(destination))) {
    stop("Refusing symlink notebook destination: ", destination, call. = FALSE)
  }
  destination_info <- file.info(destination)
  if (!file.exists(destination) || isTRUE(destination_info$isdir)) {
    stop(
      "Notebook destination must be an existing regular file: ",
      destination,
      call. = FALSE
    )
  }
  source_hash <- digest::digest(source, algo = "sha256", file = TRUE)
  temporary <- tempfile(
    pattern = paste0(".", basename(destination), "."),
    tmpdir = dirname(destination)
  )
  backup <- tempfile(
    pattern = paste0(".", basename(destination), ".prior."),
    tmpdir = dirname(destination)
  )
  displaced <- tempfile(
    pattern = paste0(".", basename(destination), ".old."),
    tmpdir = dirname(destination)
  )
  on.exit(
    {
      if (file.exists(temporary)) {
        unlink(temporary)
      }
      if (file.exists(backup)) {
        unlink(backup)
      }
      if (file.exists(displaced)) unlink(displaced)
    },
    add = TRUE
  )
  if (!file.copy(source, temporary, overwrite = FALSE)) {
    stop("Failed to copy notebook mirror: ", source, call. = FALSE)
  }
  if (nzchar(Sys.readlink(temporary))) {
    stop("Temporary notebook mirror is a symlink: ", temporary, call. = FALSE)
  }
  if (
    !identical(
      digest::digest(temporary, algo = "sha256", file = TRUE),
      source_hash
    )
  ) {
    stop("Notebook mirror hash mismatch: ", destination, call. = FALSE)
  }
  source_dimensions <- magick::image_info(
    magick::image_read(source)
  )[1L, c("width", "height")]
  temporary_dimensions <- magick::image_info(
    magick::image_read(temporary)
  )[1L, c("width", "height")]
  if (!identical(source_dimensions, temporary_dimensions)) {
    stop(
      "Temporary notebook mirror dimensions mismatch: ",
      destination,
      call. = FALSE
    )
  }
  if (!file.copy(destination, backup, overwrite = FALSE, copy.date = TRUE)) {
    stop(
      "Failed to preserve existing notebook mirror: ",
      destination,
      call. = FALSE
    )
  }
  installed <- file.rename(temporary, destination)
  if (!installed) {
    if (!file.rename(destination, displaced)) {
      stop(
        "Failed to replace notebook mirror without risking its existing file: ",
        destination,
        call. = FALSE
      )
    }
    installed <- file.rename(temporary, destination)
    if (!installed) {
      restored <- file.rename(displaced, destination)
      if (!restored) {
        restored <- file.copy(
          backup,
          destination,
          overwrite = TRUE,
          copy.date = TRUE
        )
      }
      if (!restored) {
        stop(
          "Failed to replace notebook mirror and restore its existing file: ",
          destination,
          call. = FALSE
        )
      }
      stop(
        "Failed to replace notebook mirror; existing file was preserved: ",
        destination,
        call. = FALSE
      )
    }
  }
  installed_hash <- digest::digest(destination, algo = "sha256", file = TRUE)
  if (!identical(installed_hash, source_hash)) {
    restored <- file.copy(
      backup,
      destination,
      overwrite = TRUE,
      copy.date = TRUE
    )
    if (!restored) {
      failed_destination <- tempfile(
        pattern = paste0(".", basename(destination), ".failed."),
        tmpdir = dirname(destination)
      )
      if (file.rename(destination, failed_destination)) {
        restored <- file.rename(backup, destination)
        if (!restored) {
          restored <- file.copy(
            backup,
            destination,
            overwrite = FALSE,
            copy.date = TRUE
          )
        }
        unlink(failed_destination)
      }
    }
    if (!restored) {
      stop(
        "Notebook mirror hash mismatch and existing file could not be restored: ",
        destination,
        call. = FALSE
      )
    }
    stop(
      "Installed notebook mirror hash mismatch; existing file was restored: ",
      destination,
      call. = FALSE
    )
  }
  invisible(destination)
}


# ---- six 10X folders and metadata reconciliation ----

required_metadata_columns <- c("Sample", "Mouse", "Condition")
if (!dir.exists(raw_counts_dir)) {
  stop(
    "Raw 10X count directory does not exist: ",
    raw_counts_dir,
    call. = FALSE
  )
}
if (!file.exists(metadata_path)) {
  stop("Sample metadata does not exist: ", metadata_path, call. = FALSE)
}
metadata <- utils::read.delim(
  metadata_path,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!all(required_metadata_columns %in% colnames(metadata))) {
  stop("Sample metadata lacks required columns.", call. = FALSE)
}
metadata$Condition <- sub("\\+\\s+EStim$", "+EStim", trimws(metadata$Condition))
# ANALYSIS_OK[R005]: select only the validated metadata schema; no metadata rows are dropped.
metadata <- metadata[, required_metadata_columns, drop = FALSE]
if (anyDuplicated(metadata$Sample)) {
  stop("Sample metadata contains duplicate Sample IDs.", call. = FALSE)
}
count_directories <- list.dirs(
  raw_counts_dir,
  full.names = FALSE,
  recursive = FALSE
)
# ANALYSIS_OK[R005]: remove empty directory names before the explicit sample-set equality check.
count_directories <- basename(count_directories[nzchar(count_directories)])
if (!setequal(metadata$Sample, count_directories)) {
  stop(
    "10X count folder names do not match Sample metadata IDs.",
    call. = FALSE
  )
}
sample_dirs <- stats::setNames(
  file.path(raw_counts_dir, metadata$Sample),
  metadata$Sample
)
tenx_counts <- Seurat::Read10X(
  data.dir = sample_dirs,
  gene.column = 2,
  unique.features = TRUE
)
if (is.list(tenx_counts)) {
  tenx_counts <- tenx_counts[["Gene Expression"]]
}
cell_names <- colnames(tenx_counts)
cell_sample_ids <- sub("_.*$", "", cell_names)
cell_metadata <- metadata[
  match(cell_sample_ids, metadata$Sample),
  required_metadata_columns,
  drop = FALSE
]
row.names(cell_metadata) <- cell_names
sobj <- Seurat::CreateSeuratObject(
  counts = tenx_counts,
  project = "ESPI",
  meta.data = cell_metadata,
  min.cells = 1,
  min.features = 1
)
dir.create(dirname(raw_object_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(sobj, raw_object_path)

# ---- current QC and cell calling ----

is_cell_FDR <- 0.01
mad_multiplier <- 3
MIN_CELL_COUNTS <- 108L
MIN_CELL_FEATURES <- 99L
mt_features <- rownames(sobj)[
  # ANALYSIS_OK[sample-exclusion]: fixed mitochondrial reference gene set is the audited QC contract.
  rownames(sobj) %in%
    c(
      "mt-Rnr1",
      "mt-Rnr2",
      "ND1",
      "ND2",
      "COX1",
      "COX2",
      "ATP8",
      "ATP6",
      "COX3",
      "ND3",
      "ND4L",
      "ND4",
      "ND5",
      "ND6",
      "CYTB",
      "TrnF",
      "TrnV",
      "TrnL1",
      "TrnI",
      "TrnQ",
      "TrnM",
      "TrnW",
      "TrnA",
      "TrnN",
      "TrnC",
      "TrnY",
      "TrnS1",
      "TrnD",
      "TrnK",
      "TrnG",
      "TrnR",
      "TrnH",
      "TrnS2",
      "TrnL2",
      "TrnE",
      "TrnT",
      "TrnP"
    )
]
ribo_features <- grep("^Rp[sl]", rownames(sobj), value = TRUE)
sobj[["percent.mt"]] <- Seurat::PercentageFeatureSet(
  sobj,
  features = mt_features
)
sobj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(
  sobj,
  features = ribo_features
)
set.seed(seed)
br.sobj <- DropletUtils::barcodeRanks(sobj[["RNA"]]$counts)
set.seed(seed)
e.out <- DropletUtils::emptyDrops(sobj[["RNA"]]$counts)
sobj[["cellcall_LogProb"]] <- e.out$LogProb
sobj[["cellcall_FDR"]] <- e.out$FDR
sobj@misc$knee <- S4Vectors::metadata(br.sobj)$knee
sobj[["is_cell"]] <- sobj$cellcall_FDR < is_cell_FDR
sobj_cells <- subset(
  sobj,
  subset = is_cell &
    nCount_RNA >= MIN_CELL_COUNTS &
    nFeature_RNA >= MIN_CELL_FEATURES
)
sce <- Seurat::as.SingleCellExperiment(sobj_cells, assay = "RNA")
sce <- scDblFinder::scDblFinder(
  sce,
  samples = "Sample",
  BPPARAM = BiocParallel::SerialParam(RNGseed = seed),
  dbr = 0.01
)
sobj$doublet_score <- NA_real_
sobj$doublet_call <- NA_character_
sobj$doublet_score[colnames(sce)] <- sce$scDblFinder.score
sobj$doublet_call[colnames(sce)] <- as.character(sce$scDblFinder.class)
sobj$is_singlet <- sobj$doublet_call == "singlet"

sample_cell_call_summary <- with(sobj[[]], table(Sample, is_cell))
sample_summary_table <- sobj[[]] |>
  dplyr::filter(is_cell & is_singlet) |>
  dplyr::group_by(Sample) |>
  dplyr::summarize(
    n_cells = dplyr::n(),
    median_counts = median(nCount_RNA),
    median_features = median(nFeature_RNA),
    median_percent_mt = median(percent.mt),
    median_percent_ribo = median(percent.ribo)
  )
readr::write_tsv(
  sample_summary_table,
  file.path(qc_table_dir, "sample_cell_call_summary.tsv")
)
qc_md <- sobj[[]] |>
  tibble::rownames_to_column("barcode") |>
  dplyr::filter(is_cell & is_singlet)
qc_summary_table <- qc_md |>
  dplyr::group_by(Sample) |>
  dplyr::summarize(
    n_cells = dplyr::n(),
    counts_q01 = quantile(nCount_RNA, 0.01),
    counts_q05 = quantile(nCount_RNA, 0.05),
    counts_median = median(nCount_RNA),
    counts_q95 = quantile(nCount_RNA, 0.95),
    counts_q99 = quantile(nCount_RNA, 0.99),
    features_q01 = quantile(nFeature_RNA, 0.01),
    features_q05 = quantile(nFeature_RNA, 0.05),
    features_median = median(nFeature_RNA),
    features_q95 = quantile(nFeature_RNA, 0.95),
    features_q99 = quantile(nFeature_RNA, 0.99),
    mt_median = median(percent.mt),
    mt_q90 = quantile(percent.mt, 0.90),
    mt_q95 = quantile(percent.mt, 0.95),
    mt_q99 = quantile(percent.mt, 0.99),
    ribo_median = median(percent.ribo),
    ribo_q05 = quantile(percent.ribo, 0.05),
    ribo_q95 = quantile(percent.ribo, 0.95),
    .groups = "drop"
  )
utils::write.table(
  qc_summary_table,
  file.path(qc_table_dir, "sample_qc_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
qc_thresholds <- qc_md |>
  dplyr::group_by(Sample) |>
  dplyr::summarize(
    min_count_mad = 10^(median(log10(nCount_RNA)) -
      mad_multiplier * mad(log10(nCount_RNA))),
    min_feature_mad = 10^(median(log10(nFeature_RNA)) -
      mad_multiplier * mad(log10(nFeature_RNA))),
    max_percent_mt_mad = median(percent.mt) + mad_multiplier * mad(percent.mt),
    .groups = "drop"
  )
readr::write_tsv(
  qc_thresholds,
  file.path(qc_table_dir, "sample_qc_mad_thresholds.tsv")
)
n_qc_md_before_threshold_join <- nrow(qc_md)
qc_md <- qc_md |>
  dplyr::left_join(qc_thresholds, by = "Sample")
stopifnot(nrow(qc_md) == n_qc_md_before_threshold_join)
n_sobj_metadata_before_threshold_join <- nrow(sobj[[]])
sobj[[]] <- sobj[[]] |>
  dplyr::left_join(qc_thresholds, by = "Sample")
stopifnot(nrow(sobj[[]]) == n_sobj_metadata_before_threshold_join)
sobj$fail_low_counts <- sobj$nCount_RNA < sobj$min_count_mad
sobj$fail_low_features <- sobj$nFeature_RNA < sobj$min_feature_mad
sobj$fail_high_mt <- sobj$percent.mt > sobj$max_percent_mt_mad

sobj$pass_qc <- !sobj$fail_low_counts &
  !sobj$fail_low_features &
  !sobj$fail_high_mt &
  !is.na(sobj$is_cell) &
  sobj$is_cell &
  !is.na(sobj$is_singlet) &
  sobj$is_singlet
sobj_filtered <- subset(sobj, subset = pass_qc)
dir.create(qc_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_table_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(sobj, file.path(qc_input_dir, "sobj_raw_with_qc.rds"))
saveRDS(sobj_filtered, file.path(qc_input_dir, "sobj_qc_filtered.rds"))
sobj_qc_summary_table <- sobj[[]] |>
  dplyr::filter(is_cell) |>
  dplyr::group_by(Sample) |>
  dplyr::summarize(
    called_cells = dplyr::n(),
    low_counts = sum(fail_low_counts),
    low_features = sum(fail_low_features),
    high_mt = sum(fail_high_mt),
    failed_any = sum(!pass_qc),
    retained = sum(pass_qc)
  )
readr::write_tsv(
  sobj_qc_summary_table,
  file.path(qc_table_dir, "sobj_qc_summary_by_sample.tsv")
)
sample_qc_attrition <- sobj[[]] |>
  dplyr::group_by(Sample) |>
  dplyr::summarize(
    n_barcodes = dplyr::n(),
    n_is_cell_true = sum(is_cell, na.rm = TRUE),
    n_is_cell_false = sum(!is_cell, na.rm = TRUE),
    n_is_cell_na = sum(is.na(is_cell)),
    n_scored = sum(!is.na(doublet_call)),
    n_singlet = sum(doublet_call == "singlet", na.rm = TRUE),
    n_doublet = sum(doublet_call == "doublet", na.rm = TRUE),
    n_missing_singlet = sum(is_cell & is.na(is_singlet), na.rm = TRUE),
    n_fail_low_counts = sum(fail_low_counts),
    n_fail_low_features = sum(fail_low_features),
    n_fail_high_mt = sum(fail_high_mt),
    n_pass_qc = sum(pass_qc),
    .groups = "drop"
  )
sample_qc_attrition_total <- sobj[[]] |>
  dplyr::summarize(
    Sample = "TOTAL",
    n_barcodes = dplyr::n(),
    n_is_cell_true = sum(is_cell, na.rm = TRUE),
    n_is_cell_false = sum(!is_cell, na.rm = TRUE),
    n_is_cell_na = sum(is.na(is_cell)),
    n_scored = sum(!is.na(doublet_call)),
    n_singlet = sum(doublet_call == "singlet", na.rm = TRUE),
    n_doublet = sum(doublet_call == "doublet", na.rm = TRUE),
    n_missing_singlet = sum(is_cell & is.na(is_singlet), na.rm = TRUE),
    n_fail_low_counts = sum(fail_low_counts),
    n_fail_low_features = sum(fail_low_features),
    n_fail_high_mt = sum(fail_high_mt),
    n_pass_qc = sum(pass_qc)
  )
sample_qc_attrition <- dplyr::bind_rows(
  sample_qc_attrition,
  sample_qc_attrition_total
)
readr::write_tsv(
  sample_qc_attrition,
  file.path(qc_table_dir, "sample_qc_attrition.tsv")
)

# QC plots are kept inline; they are the standard ggplot calls from the legacy
# stage and intentionally keep their existing filenames and dimensions.
knee_plot <- ggplot2::ggplot(br.sobj, ggplot2::aes(rank, total)) +
  ggplot2::geom_point(alpha = 0.5, size = 0.5) +
  ggplot2::scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  ggplot2::scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  ggplot2::geom_hline(
    yintercept = S4Vectors::metadata(br.sobj)$knee,
    linetype = "dashed",
    color = "blue"
  ) +
  ggplot2::geom_hline(
    yintercept = S4Vectors::metadata(br.sobj)$inflection,
    linetype = "dashed",
    color = "red"
  ) +
  theme_stone()
ggplot2::ggsave(
  file.path(qc_figure_dir, "knee_plot.png"),
  knee_plot,
  width = 7,
  height = 5,
  dpi = 300
)
cell_call_plot <- ggplot2::ggplot(
  e.out,
  ggplot2::aes(
    Total,
    -LogProb,
    color = dplyr::if_else(FDR <= is_cell_FDR, "Cell", "Empty", "Empty")
  )
) +
  ggplot2::geom_point(alpha = 0.7, size = 1) +
  ggplot2::scale_x_log10(
    limits = c(100, NA),
    labels = scales::label_number(big.mark = ",")
  ) +
  ggplot2::geom_vline(
    xintercept = S4Vectors::metadata(br.sobj)$knee,
    linetype = "dashed",
    color = "blue"
  ) +
  theme_stone() +
  ggplot2::labs(color = "Cell call")
ggplot2::ggsave(
  file.path(qc_figure_dir, "cell_call_plot.png"),
  cell_call_plot,
  width = 7,
  height = 5,
  dpi = 300
)
sample_cell_frac_plot <- sobj[[]] |>
  dplyr::filter(is.finite(is_cell)) |>
  dplyr::count(Sample, is_cell) |>
  dplyr::group_by(Sample) |>
  dplyr::mutate(fraction = n / sum(n)) |>
  ggplot2::ggplot(ggplot2::aes(Sample, fraction, fill = is_cell)) +
  ggplot2::geom_col() +
  ggplot2::scale_y_continuous(labels = scales::label_percent()) +
  ggplot2::labs(
    x = "Sample",
    y = "Fraction of barcodes called as cells",
    fill = "Cell call"
  ) +
  theme_stone() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggplot2::ggsave(
  file.path(qc_figure_dir, "sample_cell_call_fraction.png"),
  sample_cell_frac_plot,
  width = 7,
  height = 5,
  dpi = 300
)
qc_summary_plot <- qc_md |>
  dplyr::group_by(Sample) |>
  dplyr::select(
    barcode,
    Sample,
    nCount_RNA,
    nFeature_RNA,
    percent.mt,
    percent.ribo
  ) |>
  tidyr::pivot_longer(
    c(nCount_RNA, nFeature_RNA, percent.mt, percent.ribo),
    names_to = "metric",
    values_to = "value"
  ) |>
  ggplot2::ggplot(ggplot2::aes(Sample, value, fill = Sample)) +
  gghalves::geom_half_violin(side = "r") +
  gghalves::geom_half_boxplot(side = "l") +
  ggplot2::facet_wrap(~metric, scales = "free", ncol = 2) +
  ggplot2::labs(
    x = "Sample",
    y = "Value",
    title = "QC metrics for cells called by emptyDrops"
  ) +
  theme_stone() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
ggplot2::ggsave(
  file.path(qc_figure_dir, "sample_qc_summary.png"),
  qc_summary_plot,
  width = 7,
  height = 5,
  dpi = 300
)
count_sample_qc_plot <- ggplot2::ggplot(qc_md, ggplot2::aes(nCount_RNA)) +
  ggplot2::geom_histogram() +
  ggplot2::geom_vline(ggplot2::aes(xintercept = min_count_mad)) +
  ggplot2::scale_x_log10() +
  ggplot2::facet_wrap(~Sample, scales = "free_y") +
  theme_stone()
feature_sample_qc_plot <- ggplot2::ggplot(qc_md, ggplot2::aes(nFeature_RNA)) +
  ggplot2::geom_histogram() +
  ggplot2::geom_vline(ggplot2::aes(xintercept = min_feature_mad)) +
  ggplot2::scale_x_log10() +
  ggplot2::facet_wrap(~Sample, scales = "free_y") +
  theme_stone()
mt_sample_qc_plot <- ggplot2::ggplot(qc_md, ggplot2::aes(percent.mt)) +
  ggplot2::geom_histogram() +
  ggplot2::geom_vline(ggplot2::aes(xintercept = max_percent_mt_mad)) +
  ggplot2::facet_wrap(~Sample, scales = "free_y") +
  theme_stone()
count_feature_mt_sample_scatter <- ggplot2::ggplot(
  qc_md,
  ggplot2::aes(
    nCount_RNA,
    nFeature_RNA,
    color = percent.mt < max_percent_mt_mad
  )
) +
  ggplot2::geom_point(alpha = 0.5, size = 0.5) +
  ggplot2::scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  ggplot2::scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  ggplot2::geom_vline(ggplot2::aes(xintercept = min_count_mad)) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = min_feature_mad)) +
  ggplot2::facet_wrap(~Sample) +
  theme_stone()
feature_vs_mt <- qc_md |>
  dplyr::mutate(
    mt_pass = percent.mt < max_percent_mt_mad,
    feature_pass = nFeature_RNA > min_feature_mad,
    pass_both = mt_pass & feature_pass
  ) |>
  dplyr::group_by(Sample) |>
  dplyr::mutate(
    frac_pass = mean(pass_both),
    strip_lab = paste0(
      Sample,
      ", ",
      round(unique(frac_pass) * 100),
      "% pass both"
    )
  ) |>
  dplyr::ungroup() |>
  ggplot2::ggplot(ggplot2::aes(nFeature_RNA, percent.mt)) +
  ggplot2::geom_point(size = 0.5) +
  ggplot2::scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  ggplot2::geom_vline(ggplot2::aes(xintercept = min_feature_mad)) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = max_percent_mt_mad)) +
  ggplot2::facet_wrap(~strip_lab) +
  theme_stone()
ggplot2::ggsave(
  file.path(qc_figure_dir, "count_sample_qc_plot.png"),
  count_sample_qc_plot,
  width = 10,
  height = 6,
  dpi = 300
)
ggplot2::ggsave(
  file.path(qc_figure_dir, "feature_sample_qc_plot.png"),
  feature_sample_qc_plot,
  width = 10,
  height = 6,
  dpi = 300
)
ggplot2::ggsave(
  file.path(qc_figure_dir, "mt_sample_qc_plot.png"),
  mt_sample_qc_plot,
  width = 10,
  height = 6,
  dpi = 300
)
ggplot2::ggsave(
  file.path(qc_figure_dir, "count_feature_mt_sample_scatter.png"),
  count_feature_mt_sample_scatter,
  width = 10,
  height = 6,
  dpi = 300
)
ggplot2::ggsave(
  file.path(qc_figure_dir, "feature_vs_mt.png"),
  feature_vs_mt,
  width = 10,
  height = 6,
  dpi = 300
)

# ---- four independent source preprocessing branches ----

utils::data("mouse_cell_cycle_genes", package = "ESPI", envir = environment())
source_preprocessed <- list()
branch_numbers <- list()
for (normalization in normalizations) {
  for (filter_cc in filter_states) {
    branch_tag <- sprintf(
      "%s_%s",
      normalization,
      if (filter_cc) "filter_cc" else "no_filter_cc"
    )
    output_tag <- sprintf(
      "%s_%s",
      normalization,
      if (filter_cc) "filter-cc" else "no-filter-cc"
    )
    branch_sobj <- readRDS(file.path(qc_input_dir, "sobj_qc_filtered.rds"))
    branch_sobj@reductions <- list()
    branch_sobj[["RNA"]] <- as(branch_sobj[["RNA"]], Class = "Assay5")
    branch_sobj$sample_id <- paste0(
      "M",
      branch_sobj$Mouse,
      "_",
      gsub("[^A-Za-z0-9]", "", branch_sobj$Condition)
    )
    branch_sobj@misc$preprocessing <- list(
      normalization = normalization,
      filtered_cell_cycle = filter_cc
    )
    branch_sobj <- Seurat::FindVariableFeatures(branch_sobj, nfeatures = 2000)
    if (filter_cc) {
      SeuratObject::VariableFeatures(branch_sobj) <- setdiff(
        SeuratObject::VariableFeatures(branch_sobj),
        mouse_cell_cycle_genes
      )
    }
    branch_sobj <- if (identical(normalization, "log1p")) {
      run_log1p_pca(branch_sobj, n_pcs = 50)
    } else {
      run_pflog_pca(branch_sobj, n_pcs = 50)
    }
    branch_numbers[[output_tag]] <- tibble::tibble(
      object = output_tag,
      path = preprocess_object_path(output_tag),
      n_cells = ncol(branch_sobj),
      n_genes = nrow(branch_sobj),
      n_hvg = length(SeuratObject::VariableFeatures(branch_sobj)),
      selected_column = NA_character_,
      n_clusters = NA_integer_
    )
    pca_source_layer <- branch_sobj@misc$preprocessing$pca_source_layer
    qc_plot <- Seurat::VlnPlot(
      branch_sobj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
      group.by = "sample_id",
      ncol = 2
    )
    hvg_plot <- Seurat::LabelPoints(
      Seurat::VariableFeaturePlot(branch_sobj),
      points = head(SeuratObject::VariableFeatures(branch_sobj), 10),
      repel = TRUE,
      xnudge = 0,
      ynudge = 0,
      max.overlaps = Inf
    )
    dim_plot <- Seurat::DimHeatmap(
      branch_sobj,
      dims = 1:6,
      cells = 500,
      balanced = TRUE,
      fast = FALSE,
      combine = TRUE,
      slot = pca_source_layer,
      ncol = 2
    )
    elbow_plot <- Seurat::ElbowPlot(branch_sobj, ndims = 50, reduction = "pca")
    ggplot2::ggsave(
      file.path(
        preprocess_figure_dir,
        paste0("qc_metrics_violin_", output_tag, ".png")
      ),
      qc_plot,
      width = 8,
      height = 8
    )
    ggplot2::ggsave(
      file.path(
        preprocess_figure_dir,
        paste0("qc_metrics_violin_", output_tag, ".pdf")
      ),
      qc_plot,
      width = 8,
      height = 8
    )
    ggplot2::ggsave(
      file.path(
        preprocess_figure_dir,
        paste0("hvg_scatter_", output_tag, ".png")
      ),
      hvg_plot,
      width = 6,
      height = 5
    )
    ggplot2::ggsave(
      file.path(
        preprocess_figure_dir,
        paste0("hvg_scatter_", output_tag, ".pdf")
      ),
      hvg_plot,
      width = 6,
      height = 5
    )
    ggplot2::ggsave(
      file.path(
        preprocess_figure_dir,
        paste0("dim_heatmap_", output_tag, ".png")
      ),
      dim_plot,
      width = 8,
      height = 12
    )
    ggplot2::ggsave(
      file.path(
        preprocess_figure_dir,
        paste0("dim_heatmap_", output_tag, ".pdf")
      ),
      dim_plot,
      width = 8,
      height = 12
    )
    ggplot2::ggsave(
      file.path(preprocess_figure_dir, paste0("elbow_", output_tag, ".png")),
      elbow_plot,
      width = 5,
      height = 3
    )
    ggplot2::ggsave(
      file.path(preprocess_figure_dir, paste0("elbow_", output_tag, ".pdf")),
      elbow_plot,
      width = 5,
      height = 3
    )
    preprocess_png_paths <- file.path(
      preprocess_figure_dir,
      c(
        paste0("qc_metrics_violin_", output_tag, ".png"),
        paste0("hvg_scatter_", output_tag, ".png"),
        paste0("dim_heatmap_", output_tag, ".png"),
        paste0("elbow_", output_tag, ".png")
      )
    )
    for (preprocess_png_path in preprocess_png_paths) {
      if (basename(preprocess_png_path) %in% notebook_preprocess_pngs) {
        mirror_notebook_png(
          preprocess_png_path,
          file.path(notebook_figure_dir, basename(preprocess_png_path))
        )
      }
    }
    saveRDS(branch_sobj, preprocess_object_path(output_tag))
    source_preprocessed[[branch_tag]] <- branch_sobj
  }
}

# ---- fixed source Leiden/UMAP grids and source summaries ----

source_branches <- data.frame(
  normalization = c("log1p", "log1p", "pflog", "pflog"),
  filtered_cell_cycle = c(FALSE, TRUE, FALSE, TRUE),
  branch_tag = source_branch_tags,
  branch_label = c(
    "log1p, CC-HVG retained",
    "log1p, CC-HVG filtered",
    "PFlog, CC-HVG retained",
    "PFlog, CC-HVG filtered"
  ),
  stringsAsFactors = FALSE
)
source_clustered <- list()
for (branch_index in seq_len(nrow(source_branches))) {
  branch_info <- source_branches[branch_index, ]
  branch_sobj <- source_preprocessed[[branch_info$branch_tag]]
  candidate_names <- character()
  for (dims in dims_grid) {
    old_idents <- SeuratObject::Idents(branch_sobj)
    branch_sobj <- Seurat::FindNeighbors(
      branch_sobj,
      reduction = "pca",
      dims = 1:dims
    )
    for (resolution in resolutions) {
      column <- cluster_column(branch_info$branch_tag, dims, resolution)
      branch_sobj <- Seurat::FindClusters(
        branch_sobj,
        algorithm = 4,
        leiden_method = "igraph",
        resolution = resolution,
        random.seed = seed
      )
      branch_sobj@meta.data[[column]] <- SeuratObject::Idents(branch_sobj)
      SeuratObject::Idents(branch_sobj) <- old_idents
      candidate_names <- c(candidate_names, column)
    }
    reduction_name <- sprintf("umap_%s_dims%d", branch_info$branch_tag, dims)
    branch_sobj <- Seurat::RunUMAP(
      branch_sobj,
      reduction = "pca",
      dims = 1:dims,
      reduction.name = reduction_name,
      reduction.key = paste0(gsub("[^A-Za-z0-9]", "", reduction_name), "_"),
      seed.use = seed
    )
    for (resolution in resolutions) {
      column <- cluster_column(branch_info$branch_tag, dims, resolution)
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
      ggplot2::ggsave(png_path, plot, width = 5, height = 5)
      ggplot2::ggsave(
        sub("\\.png$", ".pdf", png_path),
        plot,
        width = 5,
        height = 5
      )
      mirror_notebook_png(
        png_path,
        file.path(notebook_figure_dir, basename(png_path))
      )
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
    ggplot2::ggsave(clustree_png, clustree_plot, width = 6, height = 6)
    ggplot2::ggsave(
      sub("\\.png$", ".pdf", clustree_png),
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
    clustree_plotted = TRUE
  )
  saveRDS(branch_sobj, cluster_object_path(branch_info$branch_tag))
  source_clustered[[branch_info$branch_tag]] <- branch_sobj
}
dir.create(cluster_table_dir, recursive = TRUE, showWarnings = FALSE)
write_cluster_grid_tables(
  source_clustered,
  source_branches[, c("normalization", "filtered_cell_cycle", "branch_tag")],
  "cluster_grid",
  cluster_table_dir
)
clustree_plots <- list()
plot_index <- 1L
for (dims in dims_grid) {
  for (branch_index in seq_len(nrow(source_branches))) {
    branch_info <- source_branches[branch_index, ]
    sobj <- source_clustered[[branch_info$branch_tag]]
    prefix <- sprintf("cluster_%s_dims%d_res", branch_info$branch_tag, dims)
    clustree_plots[[plot_index]] <- clustree::clustree(
      sobj@meta.data[,
        startsWith(colnames(sobj@meta.data), prefix),
        drop = FALSE
      ],
      prefix = prefix,
      node_text_size = 2
    ) +
      ggplot2::guides(edge_colour = "none") +
      ggplot2::ggtitle(sprintf("%s; %d PCs", branch_info$branch_label, dims)) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 10)
      )
    plot_index <- plot_index + 1L
  }
}
cluster_grid_plot <- patchwork::wrap_plots(
  clustree_plots,
  ncol = nrow(source_branches)
)
ggplot2::ggsave(
  file.path(cluster_figure_dir, "cluster_grid_clustree_12_panel.png"),
  cluster_grid_plot,
  width = 16,
  height = 12
)
mirror_notebook_png(
  file.path(cluster_figure_dir, "cluster_grid_clustree_12_panel.png"),
  file.path(notebook_figure_dir, "cluster_grid_clustree_12_panel.png")
)
ggplot2::ggsave(
  file.path(cluster_figure_dir, "cluster_grid_clustree_12_panel.pdf"),
  cluster_grid_plot,
  width = 16,
  height = 12
)
representative <- source_clustered[["pflog_filter_cc"]]
representative_plots <- lapply(resolutions, function(resolution) {
  column <- cluster_column("pflog_filter_cc", 50L, resolution)
  Seurat::DimPlot(
    representative,
    reduction = "umap_pflog_filter_cc_dims50",
    group.by = column,
    label = TRUE,
    pt.size = 0.25
  ) +
    ggplot2::ggtitle(sprintf("res %s", resolution)) +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 12)
    )
})
representative_plot <- patchwork::wrap_plots(
  representative_plots,
  ncol = length(resolutions)
) +
  patchwork::plot_annotation(title = "PFlog, CC-HVG filtered, 50 PCs")
ggplot2::ggsave(
  file.path(
    cluster_figure_dir,
    "umap_resolution_sweep_pflog_filter_cc_dims50.png"
  ),
  representative_plot,
  width = 15,
  height = 5
)
ggplot2::ggsave(
  file.path(
    cluster_figure_dir,
    "umap_resolution_sweep_pflog_filter_cc_dims50.pdf"
  ),
  representative_plot,
  width = 15,
  height = 5
)

# ---- MG selection: fixed source branch, AddModuleScore, and exclusion tests ----

utils::data("cell_type_marker_genes", package = "ESPI", envir = environment())
utils::data("cell_type_marker_labels", package = "ESPI", envir = environment())
mg_source_column <- "cluster_pflog_no_filter_cc_dims30_res0.3"
mg_source <- source_clustered[["pflog_no_filter_cc"]]
SeuratObject::DefaultAssay(mg_source) <- "RNA"
if (!inherits(mg_source[["RNA"]], "Assay5")) {
  mg_source[["RNA"]] <- as(mg_source[["RNA"]], Class = "Assay5")
}
marker_table <- stack(cell_type_marker_genes)
colnames(marker_table) <- c("gene", "cell_type")
marker_score_prefix <- "mg_selection_marker_score"
module_score_cols <- paste0(
  marker_score_prefix,
  seq_along(cell_type_marker_genes)
)
mg_source@meta.data[intersect(
  module_score_cols,
  colnames(mg_source@meta.data)
)] <- NULL
mg_source <- Seurat::AddModuleScore(
  mg_source,
  features = cell_type_marker_genes,
  assay = "RNA",
  name = marker_score_prefix,
  seed = seed,
  search = FALSE,
  slot = "data"
)
marker_score_cols <- paste0("marker_score_", names(cell_type_marker_genes))
module_scores <- mg_source@meta.data[module_score_cols]
colnames(module_scores) <- marker_score_cols
cluster_values <- as.character(mg_source@meta.data[[mg_source_column]])
cluster_levels <- sort(unique(cluster_values))
module_scores$cluster <- cluster_values
cluster_marker_scores <- stats::aggregate(
  module_scores[marker_score_cols],
  by = list(cluster = module_scores$cluster),
  FUN = mean
)
# ANALYSIS_OK[R005]: reorder cluster summaries to validated cluster order without dropping clusters.
cluster_marker_scores <- cluster_marker_scores[
  match(cluster_levels, cluster_marker_scores$cluster),
]
counts <- SeuratObject::LayerData(mg_source[["RNA"]], layer = "counts")
cdkn1b_counts <- as.numeric(counts["Cdkn1b", colnames(mg_source), drop = TRUE])
cdkn1b_detected <- as.integer(cdkn1b_counts > 0)
expression_layer <- if ("pflog" %in% SeuratObject::Layers(mg_source[["RNA"]])) {
  "pflog"
} else if ("data" %in% SeuratObject::Layers(mg_source[["RNA"]])) {
  "data"
} else {
  "counts"
}
cdkn1b_expression <- as.numeric(SeuratObject::LayerData(
  mg_source[["RNA"]],
  layer = expression_layer
)["Cdkn1b", colnames(mg_source), drop = TRUE])
mg_source$Cdkn1b_selection_expression <- cdkn1b_expression
marker_decisions <- lapply(
  seq_len(nrow(cluster_marker_scores)),
  function(index) {
    # ANALYSIS_OK[R001]: fixed rank/capture indexing preserves the audited marker-selection ordering.
    scores <- as.numeric(cluster_marker_scores[index, marker_score_cols])
    ord <- order(scores, decreasing = TRUE)
    top <- ord[[1L]]
    second <- ord[[2L]]
    # ANALYSIS_OK[R002]: fixed marker score and margin cutoffs preserve the audited MG exclusion rule.
    data.frame(
      cluster = cluster_marker_scores$cluster[[index]],
      top_marker_class = names(cell_type_marker_genes)[[top]],
      top_marker_label = unname(cell_type_marker_labels[[names(
        cell_type_marker_genes
      )[[top]]]]),
      top_marker_score = scores[[top]],
      second_marker_class = names(cell_type_marker_genes)[[second]],
      second_marker_label = unname(cell_type_marker_labels[[names(
        cell_type_marker_genes
      )[[second]]]]),
      second_marker_score = scores[[second]],
      marker_score_margin = scores[[top]] - scores[[second]],
      marker_exclude = names(cell_type_marker_genes)[[top]] %in%
        c("microglia", "photoreceptor") &&
        scores[[top]] >= 0.5 &&
        scores[[top]] - scores[[second]] >= 0.25,
      stringsAsFactors = FALSE
    )
  }
)
marker_decisions <- do.call(rbind, marker_decisions)
# ANALYSIS_OK[R026]: local Wilcoxon helper is called by the Cdkn1b selection calculations.
wilcox_greater <- function(in_values, out_values) {
  if (length(unique(c(in_values, out_values))) <= 1L) {
    return(1)
  }
  stats::wilcox.test(
    in_values,
    out_values,
    alternative = "greater",
    exact = FALSE
  )$p.value
}
cdkn1b_stats <- lapply(cluster_levels, function(cluster) {
  in_cluster <- cluster_values == cluster
  data.frame(
    cluster = cluster,
    n_cells = sum(in_cluster),
    cdkn1b_detection_fraction = mean(cdkn1b_detected[in_cluster]),
    cdkn1b_mean_expression = mean(cdkn1b_expression[in_cluster]),
    cdkn1b_expression_p = wilcox_greater(
      cdkn1b_expression[in_cluster],
      cdkn1b_expression[!in_cluster]
    ),
    cdkn1b_detection_p = wilcox_greater(
      cdkn1b_detected[in_cluster],
      cdkn1b_detected[!in_cluster]
    ),
    stringsAsFactors = FALSE
  )
})
cdkn1b_stats <- do.call(rbind, cdkn1b_stats)
cdkn1b_stats$cdkn1b_expression_q <- stats::p.adjust(
  cdkn1b_stats$cdkn1b_expression_p,
  method = "BH"
)
cdkn1b_stats$cdkn1b_detection_q <- stats::p.adjust(
  cdkn1b_stats$cdkn1b_detection_p,
  method = "BH"
)
# ANALYSIS_OK[R002]: fixed FDR and detection cutoffs preserve the audited Cdkn1b exclusion rule.
cdkn1b_stats$cdkn1b_exclude <- cdkn1b_stats$cdkn1b_expression_q < 0.05 &
  cdkn1b_stats$cdkn1b_detection_q < 0.05 &
  cdkn1b_stats$cdkn1b_detection_fraction >= 0.20
decision_table <- merge(
  merge(cluster_marker_scores, marker_decisions, by = "cluster", sort = FALSE),
  cdkn1b_stats,
  by = "cluster",
  sort = FALSE
)
# ANALYSIS_OK[R005]: reorder the decision table to validated cluster order without dropping clusters.
decision_table <- decision_table[
  match(cluster_levels, decision_table$cluster),
]
decision_table$exclude <- decision_table$marker_exclude |
  decision_table$cdkn1b_exclude
decision_table$exclusion_reasons <- vapply(
  seq_len(nrow(decision_table)),
  function(index) {
    paste(
      c(
        if (decision_table$marker_exclude[[index]]) {
          paste0("marker_", decision_table$top_marker_class[[index]])
        },
        if (decision_table$cdkn1b_exclude[[index]]) "Cdkn1b_high"
      ),
      collapse = ";"
    )
  },
  character(1)
)
decision_table$marker_exclude_classes <- "microglia;photoreceptor"
decision_table$marker_min_top_score <- 0.5
decision_table$marker_min_score_margin <- 0.25
decision_table$cdkn1b_expression_max_q <- 0.05
decision_table$cdkn1b_detection_max_q <- 0.05
decision_table$cdkn1b_min_detection_fraction <- 0.20
excluded_clusters <- decision_table$cluster[decision_table$exclude]
if (length(excluded_clusters) == 0L) {
  stop("No clusters meet exclusion criteria.", call. = FALSE)
}
dir.create(mg_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mg_table_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  decision_table,
  file.path(mg_table_dir, "mg_selected_cluster_selection.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
selection_plot <- Seurat::VlnPlot(
  mg_source,
  features = c(
    module_score_cols[match("microglia", names(cell_type_marker_genes))],
    module_score_cols[match("photoreceptor", names(cell_type_marker_genes))],
    "Cdkn1b_selection_expression"
  ),
  group.by = mg_source_column,
  assay = "RNA",
  layer = "data",
  pt.size = 0
)
ggplot2::ggsave(
  file.path(mg_figure_dir, "mg_selected_cluster_selection_diagnostics.png"),
  selection_plot,
  width = 10,
  height = 6,
  dpi = 300
)
ggplot2::ggsave(
  file.path(mg_figure_dir, "mg_selected_cluster_selection_diagnostics.pdf"),
  selection_plot,
  width = 10,
  height = 6
)
retained_cells <- colnames(mg_source)[!cluster_values %in% excluded_clusters]
base_subset <- subset(mg_source, cells = retained_cells)
base_subset@reductions <- list()
base_subset@graphs <- list()
base_subset@neighbors <- list()
base_subset@misc$clustering <- NULL
base_subset@meta.data[setdiff(
  grep(
    "^cluster_|^seurat_clusters$",
    colnames(base_subset@meta.data),
    value = TRUE
  ),
  mg_source_column
)] <- NULL
base_subset@meta.data[intersect(
  c(module_score_cols, "Cdkn1b_selection_expression"),
  colnames(base_subset@meta.data)
)] <- NULL
SeuratObject::DefaultAssay(base_subset) <- "RNA"
if (!inherits(base_subset[["RNA"]], "Assay5")) {
  base_subset[["RNA"]] <- as(base_subset[["RNA"]], Class = "Assay5")
}
# ANALYSIS_OK[R005]: retain exclusion reasons only for explicitly excluded clusters.
exclusion_reasons <- stats::setNames(
  decision_table$exclusion_reasons[decision_table$exclude],
  excluded_clusters
)
for (filter_cc in filter_states) {
  branch_tag <- paste0(
    "pflog_mg_selected_",
    if (filter_cc) "filter_cc" else "no_filter_cc"
  )
  output_tag <- paste0(
    "pflog_mg_selected_",
    if (filter_cc) "filter-cc" else "no-filter-cc"
  )
  branch_sobj <- base_subset
  branch_sobj@misc$preprocessing <- list(filtered_cell_cycle = filter_cc)
  branch_sobj <- Seurat::FindVariableFeatures(branch_sobj, nfeatures = 2000)
  if (filter_cc) {
    SeuratObject::VariableFeatures(branch_sobj) <- setdiff(
      SeuratObject::VariableFeatures(branch_sobj),
      mouse_cell_cycle_genes
    )
  }
  branch_sobj <- run_pflog_pca(branch_sobj, n_pcs = 50)
  branch_numbers[[output_tag]] <- tibble::tibble(
    object = output_tag,
    path = preprocess_object_path(output_tag),
    n_cells = ncol(branch_sobj),
    n_genes = nrow(branch_sobj),
    n_hvg = length(SeuratObject::VariableFeatures(branch_sobj)),
    selected_column = NA_character_,
    n_clusters = NA_integer_
  )
  branch_sobj@misc$preprocessing$dataset_tag <- "mg_selected"
  branch_sobj@misc$preprocessing$source_cluster_column <- mg_source_column
  branch_sobj@misc$preprocessing$source_input <- cluster_object_path(
    "pflog_no_filter_cc"
  )
  branch_sobj@misc$preprocessing$source_cluster_selection_table <- file.path(
    mg_table_dir,
    "mg_selected_cluster_selection.tsv"
  )
  branch_sobj@misc$preprocessing$source_cluster_selection_figure <- file.path(
    mg_figure_dir,
    "mg_selected_cluster_selection_diagnostics.png"
  )
  branch_sobj@misc$preprocessing$source_cluster_excluded <- excluded_clusters
  branch_sobj@misc$preprocessing$source_cluster_exclusion_reasons <- exclusion_reasons
  branch_sobj@misc$preprocessing$cdkn1b_expression_layer <- expression_layer
  branch_sobj@misc$preprocessing$marker_score_slot <- "data"
  branch_sobj@misc$preprocessing$marker_exclude_classes <- c(
    "microglia",
    "photoreceptor"
  )
  branch_sobj@misc$preprocessing$marker_min_top_score <- 0.5
  branch_sobj@misc$preprocessing$marker_min_score_margin <- 0.25
  branch_sobj@misc$preprocessing$cdkn1b_expression_max_q <- 0.05
  branch_sobj@misc$preprocessing$cdkn1b_detection_max_q <- 0.05
  branch_sobj@misc$preprocessing$cdkn1b_min_detection_fraction <- 0.20
  elbow_plot <- Seurat::ElbowPlot(branch_sobj, ndims = 50, reduction = "pca")
  elbow_stem <- paste0(
    "elbow_pflog_mg_selected_",
    if (filter_cc) "filter_cc" else "no_filter_cc"
  )
  ggplot2::ggsave(
    file.path(mg_figure_dir, paste0(elbow_stem, ".png")),
    elbow_plot,
    width = 5,
    height = 3,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(mg_figure_dir, paste0(elbow_stem, ".pdf")),
    elbow_plot,
    width = 5,
    height = 3,
    bg = "white"
  )
  saveRDS(branch_sobj, preprocess_object_path(output_tag))
  assign(branch_tag, branch_sobj)
}

# ---- both MG Leiden/UMAP grids and summaries ----

mg_branches <- data.frame(
  normalization = c("pflog", "pflog"),
  filtered_cell_cycle = c(FALSE, TRUE),
  branch_tag = mg_branch_tags,
  branch_label = c(
    "PFlog MG-selected, CC-HVG retained",
    "PFlog MG-selected, CC-HVG filtered"
  ),
  stringsAsFactors = FALSE
)
mg_clustered <- list()
for (branch_index in seq_len(nrow(mg_branches))) {
  branch_info <- mg_branches[branch_index, ]
  branch_sobj <- get(branch_info$branch_tag)
  candidate_names <- character()
  for (dims in dims_grid) {
    old_idents <- SeuratObject::Idents(branch_sobj)
    branch_sobj <- Seurat::FindNeighbors(
      branch_sobj,
      reduction = "pca",
      dims = 1:dims
    )
    for (resolution in resolutions) {
      column <- cluster_column(branch_info$branch_tag, dims, resolution)
      branch_sobj <- Seurat::FindClusters(
        branch_sobj,
        algorithm = 4,
        leiden_method = "igraph",
        resolution = resolution,
        random.seed = seed
      )
      branch_sobj@meta.data[[column]] <- SeuratObject::Idents(branch_sobj)
      SeuratObject::Idents(branch_sobj) <- old_idents
      candidate_names <- c(candidate_names, column)
    }
    reduction_name <- sprintf("umap_%s_dims%d", branch_info$branch_tag, dims)
    branch_sobj <- Seurat::RunUMAP(
      branch_sobj,
      reduction = "pca",
      dims = 1:dims,
      reduction.name = reduction_name,
      reduction.key = paste0(gsub("[^A-Za-z0-9]", "", reduction_name), "_"),
      seed.use = seed
    )
    for (resolution in resolutions) {
      column <- cluster_column(branch_info$branch_tag, dims, resolution)
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
      ggplot2::ggsave(png_path, plot, width = 5, height = 5)
      ggplot2::ggsave(
        sub("\\.png$", ".pdf", png_path),
        plot,
        width = 5,
        height = 5
      )
      mirror_notebook_png(
        png_path,
        file.path(notebook_figure_dir, basename(png_path))
      )
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
    ggplot2::ggsave(clustree_png, clustree_plot, width = 6, height = 6)
    ggplot2::ggsave(
      sub("\\.png$", ".pdf", clustree_png),
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
    clustree_plotted = TRUE
  )
  saveRDS(branch_sobj, cluster_object_path(branch_info$branch_tag))
  mg_clustered[[branch_info$branch_tag]] <- branch_sobj
}
mg_summary_rows <- list()
summary_index <- 1L
# ANALYSIS_OK[R026]: local candidate-column parser is called by the MG summary and plotting loops.
mg_selected_candidate_columns <- function(sobj, branch_tag) {
  candidate_names <- sobj@misc$clustering$candidate_names
  parsed <- lapply(candidate_names, function(column) {
    parts <- regmatches(
      column,
      regexec(
        sprintf("^cluster_%s_dims([0-9]+)_res(.+)$", branch_tag),
        column,
        perl = TRUE
      )
    )[[1]]
    # ANALYSIS_OK[R001]: fixed regex-capture indexing preserves the audited candidate-column parser.
    data.frame(
      cluster_column = column,
      dims = as.integer(parts[[2L]]),
      resolution = as.numeric(parts[[3L]]),
      stringsAsFactors = FALSE
    )
  })
  candidates <- do.call(rbind, parsed)
  # ANALYSIS_OK[R005]: order candidate columns for deterministic grid summaries without dropping candidates.
  candidates <- candidates[order(candidates$dims, candidates$resolution), ]
  rownames(candidates) <- NULL
  candidates
}
for (branch_index in seq_len(nrow(mg_branches))) {
  branch_info <- mg_branches[branch_index, ]
  sobj <- mg_clustered[[branch_info$branch_tag]]
  candidates <- mg_selected_candidate_columns(sobj, branch_info$branch_tag)
  for (candidate_index in seq_len(nrow(candidates))) {
    candidate <- candidates[candidate_index, ]
    labels <- sobj@meta.data[[candidate$cluster_column]]
    # ANALYSIS_OK[R002]: fixed small-cluster threshold preserves the audited MG summary definition.
    cluster_sizes <- table(as.character(labels))
    small_clusters <- cluster_sizes[cluster_sizes < 50L]
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
# ANALYSIS_OK[R005]: order the complete MG summary for deterministic output; no rows are dropped.
mg_summary <- mg_summary[
  order(mg_summary$branch_tag, mg_summary$dims, mg_summary$resolution),
]
utils::write.table(
  mg_summary,
  file.path(mg_table_dir, "mg_selected_cluster_grid_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
for (branch_index in seq_len(nrow(mg_branches))) {
  branch_info <- mg_branches[branch_index, ]
  sobj <- mg_clustered[[branch_info$branch_tag]]
  candidates <- mg_selected_candidate_columns(sobj, branch_info$branch_tag)
  for (dims in sort(unique(candidates$dims))) {
    # ANALYSIS_OK[R005]: select one dimensions slice for its corresponding resolution-sweep plot.
    reduction <- sprintf("umap_%s_dims%d", branch_info$branch_tag, dims)
    dim_candidates <- candidates[candidates$dims == dims, ]
    dim_candidates <- dim_candidates[order(dim_candidates$resolution), ]
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
      png_path,
      sweep,
      width = 3.6 * nrow(dim_candidates),
      height = 4.2
    )
    ggplot2::ggsave(
      file.path(mg_figure_dir, paste0(stem, ".pdf")),
      sweep,
      width = 3.6 * nrow(dim_candidates),
      height = 4.2
    )
    mirror_notebook_png(
      png_path,
      file.path(notebook_figure_dir, basename(png_path))
    )
  }
}

# ---- frozen-object numbers and structural assertions ----

frozen_object_numbers <- dplyr::bind_rows(
  tibble::tibble(
    object = "pflog_no_filter_cc",
    path = cluster_object_path("pflog_no_filter_cc"),
    n_cells = ncol(source_clustered[["pflog_no_filter_cc"]]),
    n_genes = nrow(source_clustered[["pflog_no_filter_cc"]]),
    n_hvg = NA_integer_,
    selected_column = config$selected$source$column,
    n_clusters = dplyr::n_distinct(
      source_clustered[["pflog_no_filter_cc"]][[]][[
        config$selected$source$column
      ]]
    )
  ),
  tibble::tibble(
    object = "pflog_mg_selected_no_filter_cc",
    path = cluster_object_path("pflog_mg_selected_no_filter_cc"),
    n_cells = ncol(mg_clustered[["pflog_mg_selected_no_filter_cc"]]),
    n_genes = nrow(mg_clustered[["pflog_mg_selected_no_filter_cc"]]),
    n_hvg = NA_integer_,
    selected_column = config$selected$mg$column,
    n_clusters = dplyr::n_distinct(
      mg_clustered[["pflog_mg_selected_no_filter_cc"]][[]][[
        config$selected$mg$column
      ]]
    )
  ),
  tibble::tibble(
    object = "pflog_mg_selected_filter_cc",
    path = cluster_object_path("pflog_mg_selected_filter_cc"),
    n_cells = ncol(mg_clustered[["pflog_mg_selected_filter_cc"]]),
    n_genes = nrow(mg_clustered[["pflog_mg_selected_filter_cc"]]),
    n_hvg = NA_integer_,
    selected_column = config$selected$mg_filter_cc$column,
    n_clusters = dplyr::n_distinct(
      mg_clustered[["pflog_mg_selected_filter_cc"]][[]][[
        config$selected$mg_filter_cc$column
      ]]
    )
  ),
  dplyr::bind_rows(branch_numbers)
)
readr::write_tsv(
  frozen_object_numbers,
  file.path(table_dir, "frozen_object_numbers.tsv")
)

source_cells <- ncol(source_clustered[["pflog_no_filter_cc"]])
mg_cells <- ncol(mg_clustered[["pflog_mg_selected_no_filter_cc"]])
stopifnot(
  ncol(source_clustered[["log1p_no_filter_cc"]]) == source_cells,
  ncol(source_clustered[["log1p_filter_cc"]]) == source_cells,
  ncol(source_clustered[["pflog_filter_cc"]]) == source_cells,
  mg_cells < source_cells,
  ncol(mg_clustered[["pflog_mg_selected_filter_cc"]]) == mg_cells,
  "cluster_pflog_no_filter_cc_dims30_res0.3" %in%
    colnames(source_clustered[["pflog_no_filter_cc"]][[]]),
  "cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5" %in%
    colnames(mg_clustered[["pflog_mg_selected_no_filter_cc"]][[]]),
  "cluster_pflog_mg_selected_filter_cc_dims20_res0.5" %in%
    colnames(mg_clustered[["pflog_mg_selected_filter_cc"]][[]])
)
message(
  "Frozen regeneration complete: source cells=",
  source_cells,
  "; MG cells=",
  mg_cells,
  "; QC removals=",
  sample_qc_attrition_total$n_barcodes - sample_qc_attrition_total$n_pass_qc
)
