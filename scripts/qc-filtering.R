#!/usr/bin/env Rscript

# Calculate complete mitochondrial and diagnostic ribosomal QC metrics, retain cells
# meeting conservative complexity and mitochondrial-quality thresholds, and save the filtered Seurat v5 object.
#
# Usage:
#   Rscript scripts/qc-filtering.R
#
# Outputs:
#   FIGURE_DIR/qc/*.png
#   TABLE_DIR/qc/*.tsv
#   INPUT_OBJECT_DIR/sobj_qc_filtered.rds
#
# Next step:
#   Run scripts/preprocess-sobj.R with --input-source counts-qc.

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
})
here::i_am("scripts/qc-filtering.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

input_path <- file.path(DATA_ROOT_DIR, "data", "input", "sobj_raw.rds")
figure_dir <- file.path(FIGURE_DIR, "qc")
table_dir <- file.path(TABLE_DIR, "qc")
output_path <- file.path(INPUT_OBJECT_DIR, "sobj_qc_filtered.rds")

# ---- work ----

sobj <- readRDS(input_path)
# The custom reference uses complete mouse mitochondrial features with mixed labels.
mt_features <- rownames(sobj)[
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

keep_cell <- sobj$nFeature_RNA >= 50 &
  sobj$nCount_RNA >= 100 &
  sobj$percent.mt <= 20
sobj_filtered <- sobj[, keep_cell]

# Rank every cell within sample before excluding zero-UMI cells from log-scale plots.
barcode_rank_before <- sobj[[]] |>
  dplyr::select(Sample, nCount_RNA) |>
  dplyr::group_by(Sample) |>
  dplyr::mutate(barcode_rank = dplyr::row_number(dplyr::desc(nCount_RNA))) |>
  dplyr::ungroup() |>
  dplyr::filter(nCount_RNA > 0)
barcode_rank_after <- sobj_filtered[[]] |>
  dplyr::select(Sample, nCount_RNA) |>
  dplyr::group_by(Sample) |>
  dplyr::mutate(barcode_rank = dplyr::row_number(dplyr::desc(nCount_RNA))) |>
  dplyr::ungroup() |>
  dplyr::filter(nCount_RNA > 0)

# ---- output ----

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(INPUT_OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)

filtering_decisions <- data.frame(
  metric = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.ribo",
    "high-count/high-feature cells",
    "sample membership"
  ),
  criterion = c(
    ">= 50",
    ">= 100",
    "<= 20%",
    "no threshold",
    "no threshold",
    "no sample-level exclusion"
  ),
  decision = c(
    "Required by keep_cell.",
    "Required by keep_cell.",
    "Required by keep_cell; calculated from the complete 37-feature mixed-label mitochondrial set. The >20% sparse extreme tail is excluded (Q97.5 = 19.313%).",
    "Diagnostic only; no separated high-ribosomal population (median 8.62%) and biology-dependent.",
    "Retained and visualized as possible multiplets.",
    "All samples remain in scope."
  )
)
utils::write.table(
  filtering_decisions,
  file.path(table_dir, "filtering_thresholds_and_decisions.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

global_retention <- data.frame(
  raw_cells = ncol(sobj),
  retained_cells = ncol(sobj_filtered),
  removed_cells = ncol(sobj) - ncol(sobj_filtered),
  retained_percent = 100 * ncol(sobj_filtered) / ncol(sobj),
  raw_samples = length(unique(sobj$Sample)),
  retained_samples = length(unique(sobj_filtered$Sample))
)
utils::write.table(
  global_retention,
  file.path(table_dir, "global_retention_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

sample_retention <- sobj[[]] |>
  dplyr::group_by(Sample) |>
  dplyr::summarise(
    raw_cells = dplyr::n(),
    median_nCount_RNA_before = stats::median(nCount_RNA),
    median_nFeature_RNA_before = stats::median(nFeature_RNA),
    median_percent_mt_before = stats::median(percent.mt, na.rm = TRUE),
    median_percent_ribo_before = stats::median(percent.ribo, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    sobj_filtered[[]] |>
      dplyr::group_by(Sample) |>
      dplyr::summarise(
        retained_cells = dplyr::n(),
        median_nCount_RNA_after = stats::median(nCount_RNA),
        median_nFeature_RNA_after = stats::median(nFeature_RNA),
        median_percent_mt_after = stats::median(percent.mt, na.rm = TRUE),
        median_percent_ribo_after = stats::median(percent.ribo, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Sample"
  ) |>
  dplyr::mutate(retained_percent = 100 * retained_cells / raw_cells)
utils::write.table(
  sample_retention,
  file.path(table_dir, "sample_retention_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

qc_violin_before <- Seurat::VlnPlot(
  sobj,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
  group.by = "Sample",
  pt.size = 0,
  ncol = 2
)
ggplot2::ggsave(
  file.path(figure_dir, "qc_metrics_before_filter.png"),
  qc_violin_before,
  width = 10,
  height = 8,
  dpi = 300
)

barcode_rank_plot_before <- ggplot2::ggplot(
  barcode_rank_before,
  ggplot2::aes(
    x = nCount_RNA,
    y = barcode_rank,
    colour = Sample,
    group = Sample
  )
) +
  ggplot2::geom_line(alpha = 0.7, linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 100, linetype = "dashed") +
  ggplot2::scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  ggplot2::annotation_logticks(sides = "b") +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title = "Barcode-rank QC before filtering",
    x = "UMI counts per cell",
    y = "Barcode rank within sample",
    caption = "Dashed line marks the 100-UMI retained-cell cutoff. Ranks are calculated independently within each sample; zero-UMI cells are omitted only for log-scale plotting."
  )
ggplot2::ggsave(
  file.path(figure_dir, "barcode_rank_before_filter.png"),
  barcode_rank_plot_before,
  width = 7,
  height = 5,
  dpi = 300
)

count_feature_before <- ggplot2::ggplot(
  sobj[[]],
  ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, colour = Sample)
) +
  ggplot2::geom_point(alpha = 0.05, size = 0.15) +
  ggplot2::geom_vline(xintercept = 100, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 50, linetype = "dashed") +
  ggplot2::labs(
    title = "QC before filtering: counts versus features",
    caption = "Cells above both dashed lines are retained; high values remain as possible multiplets."
  )
ggplot2::ggsave(
  file.path(figure_dir, "count_vs_feature_before_filter.png"),
  count_feature_before,
  width = 7,
  height = 5,
  dpi = 300
)

count_mito_before <- ggplot2::ggplot(
  sobj[[]],
  ggplot2::aes(x = nCount_RNA, y = percent.mt, colour = Sample)
) +
  ggplot2::geom_point(alpha = 0.05, size = 0.15) +
  ggplot2::geom_vline(xintercept = 100, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 20, linetype = "dashed") +
  ggplot2::labs(
    title = "QC before filtering: counts versus mitochondrial proportion",
    caption = "Complete 37-feature mixed-label mitochondrial proportion; the >20% sparse extreme tail is excluded."
  )
ggplot2::ggsave(
  file.path(figure_dir, "count_vs_mito_before_filter.png"),
  count_mito_before,
  width = 7,
  height = 5,
  dpi = 300
)

qc_violin_after <- Seurat::VlnPlot(
  sobj_filtered,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
  group.by = "Sample",
  pt.size = 0,
  ncol = 2
)
ggplot2::ggsave(
  file.path(figure_dir, "qc_metrics_after_filter.png"),
  qc_violin_after,
  width = 10,
  height = 8,
  dpi = 300
)

barcode_rank_plot_after <- ggplot2::ggplot(
  barcode_rank_after,
  ggplot2::aes(
    x = nCount_RNA,
    y = barcode_rank,
    colour = Sample,
    group = Sample
  )
) +
  ggplot2::geom_line(alpha = 0.7, linewidth = 0.3) +
  ggplot2::scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  ggplot2::annotation_logticks(sides = "b") +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title = "Barcode-rank QC after filtering",
    x = "UMI counts per cell",
    y = "Barcode rank within sample",
    caption = "Ranks are calculated independently within each sample; zero-UMI cells are omitted only for log-scale plotting."
  )
ggplot2::ggsave(
  file.path(figure_dir, "barcode_rank_after_filter.png"),
  barcode_rank_plot_after,
  width = 7,
  height = 5,
  dpi = 300
)

count_feature_after <- ggplot2::ggplot(
  sobj_filtered[[]],
  ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, colour = Sample)
) +
  ggplot2::geom_point(alpha = 0.2, size = 0.3) +
  ggplot2::geom_vline(xintercept = 100, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 50, linetype = "dashed") +
  ggplot2::labs(
    title = "QC after filtering: counts versus features",
    caption = "All displayed cells meet the complexity and <=20% mitochondrial thresholds."
  )
ggplot2::ggsave(
  file.path(figure_dir, "count_vs_feature_after_filter.png"),
  count_feature_after,
  width = 7,
  height = 5,
  dpi = 300
)

count_mito_after <- ggplot2::ggplot(
  sobj_filtered[[]],
  ggplot2::aes(x = nCount_RNA, y = percent.mt, colour = Sample)
) +
  ggplot2::geom_point(alpha = 0.2, size = 0.3) +
  ggplot2::geom_vline(xintercept = 100, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 20, linetype = "dashed") +
  ggplot2::labs(
    title = "QC after filtering: counts versus mitochondrial proportion",
    caption = "Complete 37-feature mixed-label mitochondrial proportion; all displayed cells meet the <=20% threshold."
  )
ggplot2::ggsave(
  file.path(figure_dir, "count_vs_mito_after_filter.png"),
  count_mito_after,
  width = 7,
  height = 5,
  dpi = 300
)

saveRDS(sobj_filtered, output_path)
message(
  "QC retained ",
  ncol(sobj_filtered),
  "/",
  ncol(sobj),
  " cells across ",
  length(unique(sobj_filtered$Sample)),
  " samples. Saved filtered object to ",
  output_path,
  ". Next step: Rscript scripts/preprocess-sobj.R --input ",
  output_path,
  " --normalization <log1p|pflog>."
)
