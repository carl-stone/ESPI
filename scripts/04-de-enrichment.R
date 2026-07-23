#!/usr/bin/env Rscript

# MG-selected pseudobulk differential expression and enrichment analysis.

suppressPackageStartupMessages({
  library(here)
  here::i_am("scripts/04-de-enrichment.R")
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(tidyverse)
  library(Seurat)
  library(DESeq2)
  library(ggview)
})

# ---- parameters ----

config <- publication_config()
input_path <- config$selected$mg$path
condition_col <- config$conditions$column
control_label <- config$conditions$control
estim_label <- config$conditions$estim
seed <- config$seed

min_gene_count <- 10L
min_samples <- 3L
min_paired_mice <- 2L
padj_cutoff <- 0.05

# ---- paths ----

deg_dir <- file.path(config$paths$degs, "mg_selected")
enrichment_dir <- file.path(config$paths$enrichment, "mg_selected")
figure_dir <- file.path(config$paths$figures, "mg_selected")
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

output_paths <- c(
  file.path(
    deg_dir,
    c(
      "pseudobulk_sample_summary.tsv",
      "design_summary.tsv",
      "deseq2_full_results.tsv",
      "deseq2_significant_degs.tsv",
      "deseq2_marker_overlap.tsv",
      "deseq2_paired_sensitivity_full_results.tsv",
      "deseq2_paired_sensitivity_significant_degs.tsv",
      "deseq2_paired_sensitivity_marker_overlap.tsv",
      "numbers.json"
    )
  ),
  file.path(
    enrichment_dir,
    c(
      "go_bp_ora_up.tsv",
      "go_bp_ora_down.tsv",
      "go_bp_gsea.tsv",
      "go_bp_gsea_symbol_entrez_mapping.tsv",
      "go_bp_ora_up_simplified.tsv",
      "go_bp_ora_down_simplified.tsv",
      "go_bp_gsea_simplified.tsv",
      "go_bp_ora_up_bayes_simplified.tsv",
      "go_bp_ora_down_bayes_simplified.tsv"
    )
  ),
  file.path(
    figure_dir,
    paste0(
      rep(
        c(
          "mg_selected_de_volcano",
          "mg_selected_go_ora_up_dotplot",
          "mg_selected_go_ora_down_dotplot",
          "mg_selected_go_gsea_dotplot",
          "mg_selected_go_ora_up_bayes_dotplot",
          "mg_selected_go_ora_down_bayes_dotplot"
        ),
        each = 2L
      ),
      c(".png", ".pdf")
    )
  )
)
assert_output_available(output_paths, config$overwrite)

# ---- pseudobulk counts ----

sobj <- readRDS(input_path)
meta <- sobj[[]] |>
  tibble::rownames_to_column("cell") |>
  dplyr::transmute(
    cell,
    Mouse = as.character(Mouse),
    condition_label = as.character(.data[[condition_col]]),
    condition = dplyr::case_when(
      condition_label == control_label ~ "control",
      condition_label == estim_label ~ "estim"
    ),
    sample_id = paste0("Mouse_", Mouse, "__", condition)
  )

sample_table <- meta |>
  tibble::as_tibble() |>
  dplyr::distinct(sample_id, Mouse, condition_label, condition) |>
  dplyr::arrange(condition, Mouse) |>
  dplyr::mutate(
    condition = factor(condition, levels = c("control", "estim")),
    mouse = factor(Mouse),
    pseudobulk_group = paste0("sample", dplyr::row_number())
  )

# ANALYSIS_OK[sample-map-join]: sample_id is unique per Mouse x Condition.
cell_groups <- meta |>
  dplyr::select("cell", "sample_id") |>
  dplyr::left_join(
    dplyr::select(sample_table, "sample_id", "pseudobulk_group"),
    by = "sample_id"
  )
if (
  nrow(cell_groups) != nrow(meta) ||
    anyDuplicated(cell_groups$cell) ||
    anyNA(cell_groups$pseudobulk_group)
) {
  cli::cli_abort(
    "Each cell must map to one Mouse x Condition pseudobulk sample."
  )
}
# ANALYSIS_OK[pseudobulk-group-alignment]: match preserves the original cell order.
sobj$pseudobulk_group <- cell_groups$pseudobulk_group[match(
  colnames(sobj),
  cell_groups$cell
)]

pseudobulk_counts <- Seurat::AggregateExpression(
  sobj,
  assays = SeuratObject::DefaultAssay(sobj),
  group.by = "pseudobulk_group"
)$RNA
pseudobulk_counts <- pseudobulk_counts[,
  sample_table$pseudobulk_group,
  drop = FALSE
]
colnames(pseudobulk_counts) <- sample_table$sample_id

sample_table$n_cells <- as.integer(table(factor(
  meta$sample_id,
  levels = sample_table$sample_id
)))
sample_table$total_counts <- Matrix::colSums(pseudobulk_counts)

mouse_conditions <- table(sample_table$Mouse, sample_table$condition)
paired_mice <- rownames(mouse_conditions)[
  mouse_conditions[, "control"] > 0 & mouse_conditions[, "estim"] > 0
]
unmatched_mice <- base::setdiff(unique(sample_table$Mouse), paired_mice)
sample_table$paired_mouse <- sample_table$Mouse %in% paired_mice

readr::write_tsv(
  sample_table,
  file.path(deg_dir, "pseudobulk_sample_table.tsv")
)

# ---- differential expression ----

keep_genes <- rowSums(pseudobulk_counts >= min_gene_count) >= min_samples
pseudobulk_counts_filtered <- pseudobulk_counts[keep_genes, , drop = FALSE]

deseq_sample_table <- sample_table |>
  dplyr::select(-"pseudobulk_group") |>
  tibble::column_to_rownames("sample_id")

dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts_filtered,
  colData = deseq_sample_table,
  design = ~condition
)

dds <- DESeq(dds)

pseudobulk_res <- results(dds, contrast = c("condition", "estim", "control"))
pseudobulk_res_shrunk <- lfcShrink(
  dds,
  coef = "condition_estim_vs_control",
  res = pseudobulk_res,
  type = "apeglm"
)

pseudobulk_res <- as.data.frame(pseudobulk_res)
pseudobulk_res_shrunk <- as.data.frame(pseudobulk_res_shrunk)

pseudobulk_res$gene <- rownames(pseudobulk_res)
pseudobulk_res_shrunk$gene <- rownames(pseudobulk_res_shrunk)

pseudobulk_res <- pseudobulk_res |>
  dplyr::relocate(gene, .before = everything())
pseudobulk_res_shrunk <- pseudobulk_res_shrunk |>
  dplyr::relocate(gene, .before = everything())


control_samples <- sample_table$sample_id[sample_table$condition == "control"]
estim_samples <- sample_table$sample_id[sample_table$condition == "estim"]

full_de <- tibble::tibble(
  gene = rownames(pseudobulk_res),
  baseMean = pseudobulk_res$baseMean,
  log2FoldChange = pseudobulk_res_shrunk$log2FoldChange,
  lfcSE = pseudobulk_res_shrunk$lfcSE,
  stat = pseudobulk_res$stat,
  pvalue = pseudobulk_res$pvalue,
  padj = pseudobulk_res$padj,
  unshrunkenLog2FoldChange = pseudobulk_res$log2FoldChange,
  mean_count_control = Matrix::rowMeans(pseudobulk_counts_filtered[,
    control_samples,
    drop = FALSE
  ]),
  mean_count_estim = Matrix::rowMeans(pseudobulk_counts_filtered[,
    estim_samples,
    drop = FALSE
  ]),
  contrast = "estim_vs_control",
  design = "unpaired_condition",
  lfc_shrink_type = "apeglm"
)
sig_de <- full_de |> dplyr::filter(!is.na(padj), padj < padj_cutoff)

utils::data("cell_type_marker_genes", package = "ESPI", envir = environment())
utils::data("cell_type_marker_labels", package = "ESPI", envir = environment())

marker_table <- stack(cell_type_marker_genes) |>
  tibble::as_tibble() |>
  dplyr::rename(gene = values, cell_type = ind) |>
  dplyr::mutate(
    cell_type = as.character(cell_type),
    cell_type_label = unname(cell_type_marker_labels[cell_type])
  ) |>
  dplyr::bind_rows(tibble::tibble(
    gene = "Cdkn1b",
    cell_type = "cdkn1b_standalone",
    cell_type_label = "Cdkn1b"
  ))
de_marker_overlap <- marker_table |>
  dplyr::filter(gene %in% rownames(pseudobulk_counts_filtered)) |>
  dplyr::left_join(full_de, by = "gene") |>
  dplyr::mutate(significant = !is.na(padj) & padj < padj_cutoff)

write_tsv(full_de, file.path(deg_dir, "deseq2_full_results.tsv"))
write_tsv(sig_de, file.path(deg_dir, "deseq2_significant_degs.tsv"))
write_tsv(de_marker_overlap, file.path(deg_dir, "deseq2_marker_overlap.tsv"))

# ---- paired-sensitivity differential expression ----

paired_status <- "skipped"
paired_reason <- "fewer than two mice have both control and E-Stim samples"
paired_de_n_degs <- NA_integer_

if (length(paired_mice) >= min_paired_mice) {
  paired_sample_table <- sample_table |>
    dplyr::filter(Mouse %in% paired_mice) |>
    dplyr::arrange(condition, Mouse) |>
    dplyr::mutate(mouse = droplevels(mouse), condition = droplevels(condition))
  paired_deseq_sample_table <- paired_sample_table |>
    dplyr::select(-"pseudobulk_group") |>
    tibble::column_to_rownames("sample_id")
  paired_counts <- pseudobulk_counts[,
    rownames(paired_deseq_sample_table),
    drop = FALSE
  ]
  # ANALYSIS_OK[pseudobulk-prefilter]: preserves the established paired-analysis count filter.
  paired_counts <- paired_counts[
    Matrix::rowSums(paired_counts) >= min_gene_count,
    ,
    drop = FALSE
  ]
  paired_design <- stats::model.matrix(
    ~ mouse + condition,
    data = paired_deseq_sample_table
  )

  if (qr(paired_design)$rank == ncol(paired_design)) {
    # ANALYSIS_OK[contrast-definition]: paired sensitivity uses the prespecified Mouse-adjusted design.
    paired_dds <- DESeqDataSetFromMatrix(
      countData = paired_counts,
      colData = paired_deseq_sample_table,
      design = ~ mouse + condition
    )
    paired_dds <- DESeq(paired_dds)
    paired_res <- results(
      paired_dds,
      contrast = c("condition", "estim", "control")
    )
    paired_res_shrunk <- lfcShrink(
      paired_dds,
      coef = "condition_estim_vs_control",
      res = paired_res,
      type = "apeglm"
    )
    paired_res <- as.data.frame(paired_res)
    paired_res_shrunk <- as.data.frame(paired_res_shrunk)
    paired_full_de <- tibble::tibble(
      gene = rownames(paired_res),
      baseMean = paired_res$baseMean,
      log2FoldChange = paired_res_shrunk$log2FoldChange,
      lfcSE = paired_res_shrunk$lfcSE,
      stat = paired_res$stat,
      pvalue = paired_res$pvalue,
      padj = paired_res$padj,
      unshrunkenLog2FoldChange = paired_res$log2FoldChange,
      mean_count_control = Matrix::rowMeans(paired_counts[,
        paired_deseq_sample_table$condition == "control",
        drop = FALSE
      ]),
      mean_count_estim = Matrix::rowMeans(paired_counts[,
        paired_deseq_sample_table$condition == "estim",
        drop = FALSE
      ]),
      contrast = "estim_vs_control",
      design = "paired_mouse_condition_sensitivity",
      lfc_shrink_type = "apeglm"
    )
    paired_sig_de <- paired_full_de |>
      dplyr::filter(!is.na(padj), padj < padj_cutoff)
    paired_de_marker_overlap <- marker_table |>
      dplyr::filter(gene %in% rownames(paired_counts)) |>
      dplyr::left_join(paired_full_de, by = "gene") |>
      dplyr::mutate(significant = !is.na(padj) & padj < padj_cutoff)

    write_tsv(
      paired_full_de,
      file.path(deg_dir, "deseq2_paired_sensitivity_full_results.tsv")
    )
    write_tsv(
      paired_sig_de,
      file.path(deg_dir, "deseq2_paired_sensitivity_significant_degs.tsv")
    )
    write_tsv(
      paired_de_marker_overlap,
      file.path(deg_dir, "deseq2_paired_sensitivity_marker_overlap.tsv")
    )

    paired_status <- "run"
    paired_reason <- "paired sensitivity used mice with both control and E-Stim samples"
    paired_de_n_degs <- nrow(paired_sig_de)
  } else {
    paired_reason <- paste0(
      "paired sensitivity design was not full rank; columns=",
      paste(colnames(paired_design), collapse = ",")
    )
  }
}

if (!identical(paired_status, "run")) {
  paired_skip <- tibble::tibble(status = "skipped", reason = paired_reason)
  write_tsv(
    paired_skip,
    file.path(deg_dir, "deseq2_paired_sensitivity_full_results.tsv")
  )
  write_tsv(
    paired_skip,
    file.path(deg_dir, "deseq2_paired_sensitivity_significant_degs.tsv")
  )
  write_tsv(
    paired_skip,
    file.path(deg_dir, "deseq2_paired_sensitivity_marker_overlap.tsv")
  )
}

# ---- volcano plot ----

load(here("data", "volcano_genes.rda"))

volcano_data <- full_de |>
  drop_na(log2FoldChange, padj) |>
  dplyr::mutate(
    neg_log10_padj = -log10(padj),
    significance = dplyr::case_when(
      padj < padj_cutoff & log2FoldChange > 0 ~ "Increased",
      padj < padj_cutoff & log2FoldChange < 0 ~ "Decreased",
      TRUE ~ "Not significant"
    ),
    significance = factor(
      significance,
      levels = c("Increased", "Decreased", "Not significant")
    ),
    label = case_when(gene %in% volcano_genes ~ gene, .default = "")
  )

volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_hline(yintercept = -log10(padj_cutoff), color = "grey75") +
  geom_vline(xintercept = c(-0.58, 0.58), linewidth = 0.25, color = "grey70") +
  geom_point(
    aes(color = significance),
    # alpha = 0.55,
    size = 0.8
  ) +
  ggrepel::geom_text_repel(
    aes(label = label),
    seed = seed,
    show.legend = FALSE,
    max.overlaps = Inf,
    max.time = 5,
    force = 2
  ) +
  scale_color_manual(
    values = c(
      "Increased" = config$palettes$analysis[["high"]],
      "Decreased" = config$palettes$analysis[["low"]],
      "Not significant" = config$palettes$analysis[["mid"]]
    ),
    drop = FALSE
  ) +
  labs(
    x = sprintf("Shrunken log2FC %s", config$conditions$contrast_display),
    y = expression(-log[10]("adjusted p-value")),
    color = NULL
  ) +
  theme_stone(base_size = 12) +
  theme(panel.grid.minor = element_blank())

save_publication_plot(
  volcano,
  file.path(figure_dir, "mg_selected_de_volcano"),
  width = 5,
  height = 4.5,
  notebook_basename = "mg_selected_de_volcano.png"
)

# ---- design summary ----

design_summary <- tibble::tibble(
  analysis = c("deseq2_primary", "deseq2_paired_sensitivity"),
  design = c("~ condition", "~ mouse + condition"),
  method = "deseq2_wald",
  status = c("run", paired_status),
  contrast = "estim_vs_control",
  limitation = c(
    "uses all Mouse × Condition samples; mouse pairing is not modeled",
    paired_reason
  ),
  included_mice = c(
    paste(unique(sample_table$Mouse), collapse = ","),
    paste(paired_mice, collapse = ",")
  ),
  unmatched_mice = paste(unique(unmatched_mice), collapse = ","),
  n_samples = c(nrow(sample_table), sum(sample_table$Mouse %in% paired_mice)),
  source_input = input_path,
  cluster_column = config$selected$mg$column,
  counts_layer = "counts",
  lfc_shrink_type = "apeglm",
  gsea_seed = seed
)
readr::write_tsv(design_summary, file.path(deg_dir, "design_summary.tsv"))

# ---- GO enrichment ----

background_map <- clusterProfiler::bitr(
  rownames(pseudobulk_counts_filtered),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db::org.Mm.eg.db
) |>
  dplyr::distinct(SYMBOL, ENTREZID)
background_entrez <- unique(background_map$ENTREZID)

up_entrez <- background_map |>
  dplyr::filter(SYMBOL %in% sig_de$gene[sig_de$log2FoldChange > 0]) |>
  dplyr::pull(ENTREZID) |>
  unique()
ora_up <- clusterProfiler::enrichGO(
  gene = up_entrez,
  universe = background_entrez,
  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)
readr::write_tsv(
  as.data.frame(ora_up) |> dplyr::mutate(direction = "up_estim_vs_control"),
  file.path(enrichment_dir, "go_bp_ora_up.tsv")
)

down_entrez <- background_map |>
  dplyr::filter(SYMBOL %in% sig_de$gene[sig_de$log2FoldChange < 0]) |>
  dplyr::pull(ENTREZID) |>
  unique()
ora_down <- clusterProfiler::enrichGO(
  gene = down_entrez,
  universe = background_entrez,
  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)
readr::write_tsv(
  as.data.frame(ora_down) |> dplyr::mutate(direction = "down_estim_vs_control"),
  file.path(enrichment_dir, "go_bp_ora_down.tsv")
)

ranked_genes <- full_de |>
  dplyr::filter(!is.na(stat)) |>
  dplyr::inner_join(background_map, by = c("gene" = "SYMBOL")) |>
  dplyr::mutate(abs_stat = abs(stat)) |>
  dplyr::arrange(ENTREZID, dplyr::desc(abs_stat), gene) |>
  dplyr::mutate(selected_for_gsea = !duplicated(ENTREZID))
readr::write_tsv(
  dplyr::select(ranked_genes, gene, ENTREZID, stat, selected_for_gsea),
  file.path(enrichment_dir, "go_bp_gsea_symbol_entrez_mapping.tsv")
)

gene_list <- ranked_genes |>
  dplyr::filter(selected_for_gsea) |>
  dplyr::arrange(dplyr::desc(stat)) |>
  dplyr::select(ENTREZID, stat) |>
  tibble::deframe()
set.seed(seed)
gsea <- clusterProfiler::gseGO(
  geneList = gene_list,
  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE
)
readr::write_tsv(
  as.data.frame(gsea),
  file.path(enrichment_dir, "go_bp_gsea.tsv")
)

s_up <- ora_up
s_up@result <- s_up@result |>
  dplyr::filter(!is.na(p.adjust), p.adjust < padj_cutoff)
s_up <- clusterProfiler::simplify(s_up)
readr::write_tsv(
  as.data.frame(s_up),
  file.path(enrichment_dir, "go_bp_ora_up_simplified.tsv")
)

s_down <- ora_down
s_down@result <- s_down@result |>
  dplyr::filter(!is.na(p.adjust), p.adjust < padj_cutoff)
s_down <- clusterProfiler::simplify(s_down)
readr::write_tsv(
  as.data.frame(s_down),
  file.path(enrichment_dir, "go_bp_ora_down_simplified.tsv")
)

s_gsea <- gsea
s_gsea@result <- s_gsea@result |>
  dplyr::filter(!is.na(p.adjust), p.adjust < padj_cutoff)
s_gsea <- clusterProfiler::simplify(s_gsea)
readr::write_tsv(
  as.data.frame(s_gsea),
  file.path(enrichment_dir, "go_bp_gsea_simplified.tsv")
)

b_up <- ora_up
b_up@result <- b_up@result |>
  dplyr::filter(!is.na(p.adjust), p.adjust < padj_cutoff)
b_up <- enrichit::bayes_enrich(b_up, seed = seed) |> clusterProfiler::simplify()
readr::write_tsv(
  as.data.frame(b_up),
  file.path(enrichment_dir, "go_bp_ora_up_bayes_simplified.tsv")
)

b_down <- ora_down
b_down@result <- b_down@result |>
  dplyr::filter(!is.na(p.adjust), p.adjust < padj_cutoff)
b_down <- enrichit::bayes_enrich(b_down, seed = seed) |>
  clusterProfiler::simplify()
readr::write_tsv(
  as.data.frame(b_down),
  file.path(enrichment_dir, "go_bp_ora_down_bayes_simplified.tsv")
)

# ---- enrichment plots ----

plot <- enrichplot::dotplot(s_up, showCategory = 15) +
  ggtitle("GO BP ORA up (simplified)")
save_publication_plot(
  plot,
  file.path(figure_dir, "mg_selected_go_ora_up_dotplot"),
  width = 8,
  height = 7,
  notebook_basename = "mg_selected_go_ora_up_dotplot.png"
)

plot <- enrichplot::dotplot(s_down, showCategory = 15) +
  ggtitle("GO BP ORA down (simplified)")
save_publication_plot(
  plot,
  file.path(figure_dir, "mg_selected_go_ora_down_dotplot"),
  width = 8,
  height = 7,
  notebook_basename = "mg_selected_go_ora_down_dotplot.png"
)

plot <- enrichplot::dotplot(s_gsea, showCategory = 15, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("GO BP GSEA (simplified)")
save_publication_plot(
  plot,
  file.path(figure_dir, "mg_selected_go_gsea_dotplot"),
  width = 8,
  height = 7,
  notebook_basename = "mg_selected_go_gsea_dotplot.png"
)

plot <- enrichplot::dotplot(b_up, showCategory = 15) +
  ggtitle("GO BP ORA up (bayes + simplified)")
save_publication_plot(
  plot,
  file.path(figure_dir, "mg_selected_go_ora_up_bayes_dotplot"),
  width = 8,
  height = 7,
  notebook_basename = "mg_selected_go_ora_up_bayes_dotplot.png"
)

plot <- enrichplot::dotplot(b_down, showCategory = 15) +
  ggtitle("GO BP ORA down (bayes + simplified)")
save_publication_plot(
  plot,
  file.path(figure_dir, "mg_selected_go_ora_down_bayes_dotplot"),
  width = 8,
  height = 7,
  notebook_basename = "mg_selected_go_ora_down_bayes_dotplot.png"
)

# ---- reportable values ----

jsonlite::write_json(
  list(
    n_samples = nrow(sample_table),
    n_cells = ncol(sobj),
    n_tested_genes = nrow(pseudobulk_counts_filtered),
    n_degs = nrow(sig_de),
    n_marker_degs = sum(de_marker_overlap$significant),
    n_unmatched_mice = length(unique(unmatched_mice)),
    source_input = input_path,
    cluster_column = config$selected$mg$column,
    counts_layer = "counts",
    gsea_seed = seed,
    lfc_shrink_type = "apeglm",
    primary_design = "~ condition",
    paired_sensitivity_design = "~ mouse + condition",
    paired_sensitivity_included_mice = paste(paired_mice, collapse = ","),
    paired_sensitivity_status = paired_status,
    paired_sensitivity_n_degs = paired_de_n_degs
  ),
  file.path(deg_dir, "numbers.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  na = "null"
)

message("Saved MG-selected DE outputs under ", deg_dir)
message("Saved MG-selected enrichment outputs under ", enrichment_dir)
