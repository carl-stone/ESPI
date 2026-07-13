#!/usr/bin/env Rscript

# Call cells with emptyDrops, calculate sample-specific MAD QC thresholds,
# save annotated raw and filtered Seurat v5 objects, and write QC diagnostics.
#
# Usage:
#   Rscript scripts/02-qc-filtering.R
#
# Outputs:
#   FIGURE_DIR/qc/*.png
#   TABLE_DIR/qc/*.tsv
#   INPUT_OBJECT_DIR/sobj_qc_filtered.rds
#
# Next step:
#   Run scripts/03-preprocess.R with --input-source counts-qc.

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
  library(DropletUtils)
  library(ggplot2)
  library(gghalves) # erocoar github
  library(dplyr)
  library(scDblFinder)
})
here::i_am("scripts/02-qc-filtering.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

SEED <- 1312

input_path <- file.path(DATA_ROOT_DIR, "data", "input", "sobj_raw.rds")
figure_dir <- file.path(FIGURE_DIR, "qc")
table_dir <- file.path(TABLE_DIR, "qc")
output_path <- file.path(INPUT_OBJECT_DIR, "sobj_qc_filtered.rds")

is_cell_FDR <- 0.01
mad_multiplier <- 3

# ---- work ----

sobj <- readRDS(input_path)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(INPUT_OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)

# Empty droplet ID before any filtering
set.seed(ESPI::SEED)
br.sobj <- DropletUtils::barcodeRanks(sobj[["RNA"]]$counts)

# Plot knee plot
knee_plot <- ggplot2::ggplot(
  br.sobj,
  ggplot2::aes(x = rank, y = total)
) +
  ggplot2::geom_point(alpha = 0.5, size = 0.5) +
  ggplot2::scale_x_log10(
    labels = scales::label_number(big.mark = ",")
  ) +
  ggplot2::scale_y_log10(
    labels = scales::label_number(big.mark = ",")
  ) +
  ggplot2::geom_hline(
    yintercept = metadata(br.sobj)$knee,
    linetype = "dashed",
    color = "blue"
  ) +
  ggplot2::geom_hline(
    yintercept = metadata(br.sobj)$inflection,
    linetype = "dashed",
    color = "red"
  ) +
  ESPI::theme_stone()
# save plot
ggplot2::ggsave(
  file.path(figure_dir, "knee_plot.png"),
  knee_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# Cell calling with DropletUtils permutation test
set.seed(ESPI::SEED)
e.out <- DropletUtils::emptyDrops(sobj[["RNA"]]$counts)

# Plot cell calling results
cell_call_plot <- ggplot2::ggplot(
  e.out,
  ggplot2::aes(
    x = Total,
    y = -LogProb,
    color = dplyr::if_else(FDR <= 0.01, "Cell", "Empty", "Empty")
  )
) +
  ggplot2::geom_point(alpha = 0.7, size = 1) +
  ggplot2::scale_x_log10(
    limits = c(100, NA),
    labels = scales::label_number(big.mark = ",")
  ) +
  ggplot2::geom_vline(
    xintercept = metadata(br.sobj)$knee,
    linetype = "dashed",
    color = "blue"
  ) +
  ESPI::theme_stone() +
  ggplot2::labs(
    color = "Cell call"
  )
# save plot
ggplot2::ggsave(
  file.path(figure_dir, "cell_call_plot.png"),
  cell_call_plot,
  width = 7,
  height = 5,
  dpi = 300
)


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

# Add cell calling metadata to sobj
# br.sobj has knee threshold
# e.out has emptyDrops LogProb and FDR

sobj[["cellcall_LogProb"]] <- e.out$LogProb
sobj[["cellcall_FDR"]] <- e.out$FDR
sobj@misc$knee <- metadata(br.sobj)$knee

# Set is_cell for cell call threshold based off FDR
sobj[["is_cell"]] <- sobj$cellcall_FDR < is_cell_FDR

# Identify floors for UMI and features before doublet detection
sobj[[]] |>
  dplyr::filter(is_cell) |>
  dplyr::summarize(
    n = dplyr::n(),
    min_counts = min(nCount_RNA),
    counts_q01 = quantile(nCount_RNA, 0.01),
    counts_q05 = quantile(nCount_RNA, 0.05),
    min_features = min(nFeature_RNA),
    features_q01 = quantile(nFeature_RNA, 0.01),
    features_q05 = quantile(nFeature_RNA, 0.05)
  )

# Doublet filtering
sobj_cells <- subset(
  sobj,
  subset = is_cell &
    nCount_RNA >= 108 &
    nFeature_RNA >= 99
)
sce <- as.SingleCellExperiment(sobj_cells, assay = "RNA")

bp <- BiocParallel::SerialParam(RNGseed = ESPI::SEED)

sce <- scDblFinder(
  sce,
  samples = "Sample",
  BPPARAM = bp,
  dbr = 0.01
)

sobj$doublet_score <- NA_real_
sobj$doublet_call <- NA_character_

sobj$doublet_score[colnames(sce)] <- sce$scDblFinder.score
sobj$doublet_call[colnames(sce)] <- as.character(sce$scDblFinder.class)
sobj$is_singlet <- sobj$doublet_call == "singlet"

sample_cell_call_summary <- with(
  sobj[[]],
  table(Sample, is_cell)
)

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
  file.path(table_dir, "sample_cell_call_summary.tsv")
)

qc_md <- sobj[[]] |>
  tibble::rownames_to_column("barcode") |>
  dplyr::filter(is_cell & is_singlet)

sample_cell_frac_plot <- sobj[[]] |>
  dplyr::filter(is.finite(is_cell)) |>
  dplyr::count(Sample, is_cell) |>
  dplyr::group_by(Sample) |>
  dplyr::mutate(fraction = n / sum(n)) |>
  ggplot2::ggplot(
    ggplot2::aes(x = Sample, y = fraction, fill = is_cell)
  ) +
  ggplot2::geom_col() +
  ggplot2::scale_y_continuous(labels = scales::label_percent()) +
  ggplot2::labs(
    x = "Sample",
    y = "Fraction of barcodes called as cells",
    fill = "Cell call"
  ) +
  ESPI::theme_stone() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggsave(
  file.path(figure_dir, "sample_cell_call_fraction.png"),
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
    cols = c(nCount_RNA, nFeature_RNA, percent.mt, percent.ribo),
    names_to = "metric",
    values_to = "value"
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(x = Sample, y = value, fill = Sample)
  ) +
  gghalves::geom_half_violin(side = "r") +
  gghalves::geom_half_boxplot(side = "l") +
  ggplot2::facet_wrap(~metric, scales = "free", ncol = 2) +
  ggplot2::labs(
    x = "Sample",
    y = "Value",
    title = "QC metrics for cells called by emptyDrops"
  ) +
  ESPI::theme_stone() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(
  file.path(figure_dir, "sample_qc_summary.png"),
  qc_summary_plot,
  width = 7,
  height = 5,
  dpi = 300
)

qc_summary_table <- qc_md |>
  group_by(Sample) |>
  summarize(
    n_cells = n(),

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

write.table(
  qc_summary_table,
  file.path(table_dir, "sample_qc_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

qc_thresholds <- qc_md |>
  group_by(Sample) |>
  summarize(
    min_count_mad = 10^(median(log10(nCount_RNA)) -
      mad_multiplier * mad(log10(nCount_RNA))),
    min_feature_mad = 10^(median(log10(nFeature_RNA)) -
      mad_multiplier * mad(log10(nFeature_RNA))),
    max_percent_mt_mad = median(percent.mt) +
      mad_multiplier * mad(percent.mt),
    .groups = "drop"
  )
readr::write_tsv(
  qc_thresholds,
  file.path(table_dir, "sample_qc_mad_thresholds.tsv")
)

qc_md <- qc_md |>
  left_join(qc_thresholds, by = "Sample")

count_sample_qc_plot <- ggplot(qc_md, aes(x = nCount_RNA)) +
  geom_histogram() +
  geom_vline(aes(xintercept = min_count_mad), linetype = "dashed") +
  scale_x_log10() +
  facet_wrap(~Sample, scales = "free_y") +
  theme_stone()

ggsave(
  plot = count_sample_qc_plot,
  filename = file.path(figure_dir, "count_sample_qc_plot.png"),
  width = 10,
  height = 6,
  dpi = 300
)

feature_sample_qc_plot <- ggplot(qc_md, aes(x = nFeature_RNA)) +
  geom_histogram() +
  geom_vline(aes(xintercept = min_feature_mad), linetype = "dashed") +
  scale_x_log10() +
  facet_wrap(~Sample, scales = "free_y") +
  theme_stone()

ggsave(
  plot = feature_sample_qc_plot,
  filename = file.path(figure_dir, "feature_sample_qc_plot.png"),
  width = 10,
  height = 6,
  dpi = 300
)

mt_sample_qc_plot <- ggplot(qc_md, aes(x = percent.mt)) +
  geom_histogram() +
  geom_vline(aes(xintercept = max_percent_mt_mad), linetype = "dashed") +
  facet_wrap(~Sample, scales = "free_y") +
  theme_stone()

ggsave(
  plot = mt_sample_qc_plot,
  filename = file.path(figure_dir, "mt_sample_qc_plot.png"),
  width = 10,
  height = 6,
  dpi = 300
)

# ggplot(qc_md, aes(x = nCount_RNA, color = Sample)) +
#   geom_density() +
#   scale_x_continuous(
#     transform = scales::log10_trans(),
#     guide = "axis_logticks",
#     labels = scales::label_number(big.mark = ",")
#   ) +
#   theme_stone() +
#   labs(
#     x = "UMI count",
#     y = "Density",
#     color = "Sample"
#   )

# ggplot(qc_md, aes(x = nFeature_RNA, color = Sample)) +
#   geom_density() +
#   scale_x_continuous(
#     transform = scales::log10_trans(),
#     guide = "axis_logticks",
#     labels = scales::label_number(big.mark = ",")
#   ) +
#   theme_stone() +
#   labs(
#     x = "Feature count",
#     y = "Density",
#     color = "Sample"
#   )

count_feature_mt_sample_scatter <- ggplot(
  qc_md,
  aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt < max_percent_mt_mad)
) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  geom_vline(aes(xintercept = min_count_mad), linetype = "dashed") +
  geom_hline(aes(yintercept = min_feature_mad), linetype = "dashed") +
  facet_wrap(~Sample) +
  # scale_color_viridis_b(
  #   breaks = c(5, 20),
  #   oob = scales::squish_infinite
  # ) +
  theme_stone()

ggsave(
  plot = count_feature_mt_sample_scatter,
  filename = file.path(figure_dir, "count_feature_mt_sample_scatter.png"),
  width = 10,
  height = 6,
  dpi = 300
)

feature_vs_mt <- qc_md |>
  mutate(
    mt_pass = percent.mt < max_percent_mt_mad,
    feature_pass = nFeature_RNA > min_feature_mad,
    pass_both = mt_pass & feature_pass
  ) |>
  group_by(Sample) |>
  mutate(
    frac_pass = mean(pass_both),
    strip_lab = paste0(
      Sample,
      ", ",
      round(unique(frac_pass) * 100),
      "% pass both"
    )
  ) |>
  ungroup() |>
  ggplot(aes(x = nFeature_RNA, y = percent.mt)) +
  geom_point(size = 0.5) +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  geom_vline(aes(xintercept = min_feature_mad), linetype = "dashed") +
  geom_hline(aes(yintercept = max_percent_mt_mad), linetype = "dashed") +
  facet_wrap(~strip_lab) +
  theme_stone()

ggsave(
  plot = feature_vs_mt,
  filename = file.path(figure_dir, "feature_vs_mt.png"),
  width = 10,
  height = 6,
  dpi = 300
)

## Now set flags for each threshold
sobj[[]] <- sobj[[]] |>
  left_join(qc_thresholds, by = "Sample")

sobj$fail_low_counts <- sobj$nCount_RNA < sobj$min_count_mad
sobj$fail_low_features <- sobj$nFeature_RNA < sobj$min_feature_mad
sobj$fail_high_mt <- sobj$percent.mt > sobj$max_percent_mt_mad

sobj$pass_qc <- TRUE
sobj$pass_qc <- sobj$pass_qc &
  !sobj$fail_low_counts &
  !sobj$fail_low_features &
  !sobj$fail_high_mt

sobj_filtered <- subset(sobj, subset = pass_qc)

saveRDS(sobj, file.path(INPUT_OBJECT_DIR, "sobj_raw_with_qc.rds"))
saveRDS(sobj_filtered, file.path(INPUT_OBJECT_DIR, "sobj_qc_filtered.rds"))

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
  file.path(table_dir, "sobj_qc_summary_by_sample.tsv")
)
