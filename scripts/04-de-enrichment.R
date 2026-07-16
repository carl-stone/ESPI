#!/usr/bin/env Rscript

# Fixed MG-selected pseudobulk differential expression and enrichment analysis.

suppressPackageStartupMessages({
  here::i_am("scripts/04-de-enrichment.R")
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

config <- publication_config()
seed <- config$seed
conditions <- config$conditions
selected_mg <- config$selected$mg
input_path <- selected_mg$path
cluster_column <- selected_mg$column
condition_col <- conditions$column
control_label <- conditions$control
estim_label <- conditions$estim
counts_layer <- "counts"
lfc_shrink_type <- "apeglm"

control_level <- "control"
estim_level <- "estim"
contrast_direction <- "estim_vs_control"
cdkn1b_gene <- "Cdkn1b"
min_pseudobulk_gene_count <- 10L
min_ora_genes <- 5L
min_gsea_genes <- 10L
enrichment_significance_cutoff <- 0.05
de_significance_cutoff <- 0.05
gsea_min_gene_set_size <- 10L
gsea_max_gene_set_size <- 500L
volcano_label_limit <- 20L
min_paired_mice <- 2L
enrichment_dotplot_show_n <- 15L

# ---- fixed paths and overwrite guard ----

deg_dir <- file.path(config$paths$degs, "mg_selected")
enrichment_dir <- file.path(config$paths$enrichment, "mg_selected")
figure_dir <- file.path(config$paths$figures, "mg_selected")

primary_paths <- c(
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
      c(
        "mg_selected_de_volcano",
        "mg_selected_go_ora_up_dotplot",
        "mg_selected_go_ora_down_dotplot",
        "mg_selected_go_gsea_dotplot",
        "mg_selected_go_ora_up_bayes_dotplot",
        "mg_selected_go_ora_down_bayes_dotplot"
      ),
      rep(c(".png", ".pdf"), each = 6L)
    )
  )
)
assert_output_available(primary_paths, config$overwrite)

# ---- narrow serialization and identity checks ----

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA"
  )
}

# ANALYSIS_OK[R026]: local skipped-status writer is called by the DE/enrichment branches below.
write_reason_tsv <- function(
  path,
  reason,
  n_input_genes,
  n_mapped_genes,
  min_required
) {
  write_tsv(
    data.frame(
      status = "skipped",
      reason = reason,
      n_input_genes = as.integer(n_input_genes),
      n_mapped_genes = as.integer(n_mapped_genes),
      min_required = as.integer(min_required),
      stringsAsFactors = FALSE
    ),
    path
  )
}

# ANALYSIS_OK[R026]: local path-component helper is called by this executable script.
safe_component <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "missing")
}

# ANALYSIS_OK[R026]: local metadata validator is called before DE model construction.
assert_metadata <- function(meta, columns) {
  missing_columns <- setdiff(columns, colnames(meta))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required metadata column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  bad_columns <- columns[vapply(
    columns,
    function(column) {
      values <- meta[[column]]
      any(is.na(values) | trimws(as.character(values)) == "")
    },
    logical(1)
  )]
  if (length(bad_columns) > 0L) {
    stop(
      "Missing values in required metadata column(s): ",
      paste(bad_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

# ANALYSIS_OK[R026]: local design-rank validator is called by each DE model branch.
assert_full_rank <- function(design, label) {
  design_rank <- qr(design)$rank
  if (design_rank < ncol(design)) {
    stop(
      label,
      " design matrix is not full rank (rank ",
      design_rank,
      " < ",
      ncol(design),
      "). Columns: ",
      paste(colnames(design), collapse = ", "),
      call. = FALSE
    )
  }
}

# ---- frozen MG object and Mouse x Condition pseudobulk samples ----

if (!file.exists(input_path)) {
  stop("Final MG Seurat object does not exist: ", input_path, call. = FALSE)
}
sobj <- readRDS(input_path)
assert_frozen_input(input_path, sobj, config$frozen$mg)
meta <- sobj@meta.data
if (!is.data.frame(meta) || nrow(meta) == 0L) {
  stop("Final MG object has missing or empty metadata.", call. = FALSE)
}
assert_metadata(meta, c("Mouse", condition_col, cluster_column))

condition_labels <- trimws(as.character(meta[[condition_col]]))
expected_conditions <- sort(c(control_label, estim_label))
observed_conditions <- sort(unique(condition_labels))
if (!identical(observed_conditions, expected_conditions)) {
  stop(
    "Bad contrast labels in metadata column '",
    condition_col,
    "'. Expected exactly: ",
    paste(expected_conditions, collapse = ", "),
    "; observed: ",
    paste(observed_conditions, collapse = ", "),
    call. = FALSE
  )
}
condition <- ifelse(
  condition_labels == control_label,
  control_level,
  ifelse(condition_labels == estim_label, estim_level, NA_character_)
)
if (anyNA(condition)) {
  stop("Internal condition recoding failed.", call. = FALSE)
}
condition <- factor(condition, levels = c(control_level, estim_level))
mouse <- trimws(as.character(meta$Mouse))
constructed_sample_id <- paste(
  paste0("Mouse_", safe_component(mouse)),
  as.character(condition),
  sep = "__"
)

if ("sample_id" %in% colnames(meta)) {
  assert_metadata(meta, "sample_id")
  sample_id <- trimws(as.character(meta$sample_id))
  sample_map <- unique(data.frame(
    sample_id = sample_id,
    Mouse = mouse,
    condition = as.character(condition),
    constructed_sample_id = constructed_sample_id,
    stringsAsFactors = FALSE
  ))
  sample_to_units <- stats::aggregate(
    constructed_sample_id ~ sample_id,
    sample_map,
    function(x) length(unique(x))
  )
  if (any(sample_to_units$constructed_sample_id != 1L)) {
    stop(
      "Metadata sample_id values are not consistent with Mouse x Condition units: ",
      paste(
        sample_to_units$sample_id[sample_to_units$constructed_sample_id != 1L],
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  unit_to_samples <- stats::aggregate(
    sample_id ~ constructed_sample_id,
    sample_map,
    function(x) length(unique(x))
  )
  if (any(unit_to_samples$sample_id != 1L)) {
    stop(
      "Multiple sample_id values found within Mouse x Condition unit(s): ",
      paste(
        unit_to_samples$constructed_sample_id[unit_to_samples$sample_id != 1L],
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  pseudobulk_sample_id <- sample_id
  sample_id_source <- "metadata_sample_id"
} else {
  pseudobulk_sample_id <- constructed_sample_id
  sample_id_source <- "constructed_mouse_condition"
}

sample_table <- unique(data.frame(
  sample_id = pseudobulk_sample_id,
  Mouse = mouse,
  mouse = factor(mouse),
  condition_label = condition_labels,
  condition = condition,
  constructed_sample_id = constructed_sample_id,
  stringsAsFactors = FALSE
))
if (anyDuplicated(sample_table$sample_id)) {
  stop(
    "Pseudobulk sample IDs are duplicated across distinct Mouse x Condition units.",
    call. = FALSE
  )
}
# ANALYSIS_OK[R002]: fixed six-sample cardinality preserves the audited Mouse x Condition design.
if (nrow(sample_table) != 6L) {
  stop(
    "Expected exactly six Mouse x Condition pseudobulk samples; found ",
    nrow(sample_table),
    call. = FALSE
  )
}
# ANALYSIS_OK[R005]: order the complete pseudobulk sample table deterministically without dropping samples.
sample_table <- sample_table[
  order(sample_table$condition, sample_table$Mouse),
  ,
  drop = FALSE
]
rownames(sample_table) <- sample_table$sample_id
sample_table$condition <- factor(
  as.character(sample_table$condition),
  levels = c(control_level, estim_level)
)
sample_table$mouse <- factor(sample_table$Mouse)
condition_counts <- table(sample_table$condition)
if (any(condition_counts[c(control_level, estim_level)] < 1L)) {
  stop(
    "Both conditions must have at least one Mouse x Condition sample.",
    call. = FALSE
  )
}

assay <- SeuratObject::DefaultAssay(sobj)
available_layers <- SeuratObject::Layers(sobj[[assay]])
if (!counts_layer %in% available_layers) {
  stop(
    "Missing raw counts layer '",
    counts_layer,
    "' in assay ",
    assay,
    ". Available layers: ",
    paste(available_layers, collapse = ", "),
    call. = FALSE
  )
}
counts <- SeuratObject::GetAssayData(sobj, assay = assay, layer = counts_layer)
if (!inherits(counts, "Matrix")) {
  counts <- Matrix::Matrix(counts, sparse = TRUE)
}
if (is.null(rownames(counts)) || any(!nzchar(rownames(counts)))) {
  stop("Counts matrix must have non-empty gene row names.", call. = FALSE)
}
if (is.null(colnames(counts)) || !identical(colnames(counts), rownames(meta))) {
  stop(
    "Counts matrix columns do not exactly match Seurat metadata row names.",
    call. = FALSE
  )
}
count_values <- if (inherits(counts, "sparseMatrix")) {
  counts@x
} else {
  as.numeric(counts)
}
if (length(count_values) > 0L && any(count_values < 0)) {
  stop("Raw counts layer contains negative values.", call. = FALSE)
}
if (
  length(count_values) > 0L &&
    any(abs(count_values - round(count_values)) > sqrt(.Machine$double.eps))
) {
  stop("Raw counts layer contains non-integer values.", call. = FALSE)
}

sample_ids <- rownames(sample_table)
pseudobulk_counts <- do.call(
  cbind,
  lapply(sample_ids, function(id) {
    Matrix::rowSums(counts[, pseudobulk_sample_id == id, drop = FALSE])
  })
)
rownames(pseudobulk_counts) <- rownames(counts)
colnames(pseudobulk_counts) <- sample_ids
storage.mode(pseudobulk_counts) <- "integer"

sample_table$n_cells <- as.integer(tabulate(
  match(pseudobulk_sample_id, sample_ids),
  nbins = length(sample_ids)
))
sample_table$total_counts <- as.numeric(Matrix::colSums(pseudobulk_counts))
if (any(sample_table$n_cells < 1L)) {
  stop(
    "At least one Mouse x Condition pseudobulk sample has zero cells.",
    call. = FALSE
  )
}
mouse_condition_table <- table(sample_table$Mouse, sample_table$condition)
paired_mice <- rownames(mouse_condition_table)[
  mouse_condition_table[, control_level] > 0L &
    mouse_condition_table[, estim_level] > 0L
]
unmatched_mice <- setdiff(as.character(sample_table$Mouse), paired_mice)
sample_table$paired_mouse <- sample_table$Mouse %in% paired_mice
sample_table$sample_id_source <- sample_id_source
sample_table$analysis_unit <- "Mouse_x_Condition"

# ---- DESeq2 ----

# ANALYSIS_OK[R026]: local DE result serializer is called by both DE model branches.
results_to_table <- function(
  unshrunk_result,
  shrunk_result,
  count_matrix,
  sample_data,
  shrink_type,
  design_label
) {
  unshrunk <- as.data.frame(unshrunk_result)
  shrunk <- as.data.frame(shrunk_result)
  genes <- rownames(unshrunk)
  if (!identical(genes, rownames(shrunk))) {
    stop(
      "DESeq2 unshrunk and shrunken result rows do not match.",
      call. = FALSE
    )
  }
  control_samples <- rownames(sample_data)[
    sample_data$condition == control_level
  ]
  estim_samples <- rownames(sample_data)[sample_data$condition == estim_level]
  data.frame(
    gene = genes,
    baseMean = unshrunk$baseMean,
    log2FoldChange = shrunk$log2FoldChange,
    lfcSE = shrunk$lfcSE,
    stat = unshrunk$stat,
    pvalue = unshrunk$pvalue,
    padj = unshrunk$padj,
    unshrunkenLog2FoldChange = unshrunk$log2FoldChange,
    mean_count_control = rowMeans(as.matrix(count_matrix[,
      control_samples,
      drop = FALSE
    ])),
    mean_count_estim = rowMeans(as.matrix(count_matrix[,
      estim_samples,
      drop = FALSE
    ])),
    contrast = contrast_direction,
    design = design_label,
    lfc_shrink_type = shrink_type,
    stringsAsFactors = FALSE
  )
}

# ANALYSIS_OK[R026]: local DESeq2 runner is called by the primary and paired model branches.
run_deseq2 <- function(
  count_matrix,
  sample_data,
  design_formula,
  design_label
) {
  design <- stats::model.matrix(design_formula, data = sample_data)
  assert_full_rank(design, paste("DESeq2", design_label))
  keep <- Matrix::rowSums(count_matrix) >= min_pseudobulk_gene_count
  if (!any(keep)) {
    stop(
      "No genes passed DESeq2 pre-filter rowSums >= ",
      min_pseudobulk_gene_count,
      ".",
      call. = FALSE
    )
  }
  filtered_counts <- count_matrix[keep, , drop = FALSE]
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = sample_data,
    design = design_formula
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  coefficient <- "condition_estim_vs_control"
  available_results <- DESeq2::resultsNames(dds)
  if (!coefficient %in% available_results) {
    stop(
      "Expected DESeq2 coefficient not found for estim vs control: ",
      coefficient,
      ". Available coefficient(s): ",
      paste(available_results, collapse = ", "),
      call. = FALSE
    )
  }
  unshrunk_result <- DESeq2::results(
    dds,
    contrast = c("condition", estim_level, control_level)
  )
  # ANALYSIS_OK[warning-suppression]: DESeq2 shrinkage messages are intentionally suppressed after model validation.
  shrunk_result <- suppressMessages(DESeq2::lfcShrink(
    dds,
    coef = coefficient,
    res = unshrunk_result,
    type = lfc_shrink_type,
    quiet = TRUE
  ))
  list(
    result_table = results_to_table(
      unshrunk_result,
      shrunk_result,
      filtered_counts,
      sample_data,
      lfc_shrink_type,
      design_label
    ),
    tested_genes = rownames(filtered_counts),
    shrink_type = lfc_shrink_type
  )
}

# ---- curated marker overlap ----

utils::data("cell_type_marker_genes", package = "ESPI", envir = environment())
utils::data("cell_type_marker_labels", package = "ESPI", envir = environment())
# ANALYSIS_OK[R026]: local curated marker-table builder is called by marker-overlap output.
make_marker_table <- function(gene_universe) {
  if (
    !identical(names(cell_type_marker_genes), names(cell_type_marker_labels))
  ) {
    stop("Curated marker gene and label names do not match.", call. = FALSE)
  }
  marker_table <- stack(cell_type_marker_genes)
  colnames(marker_table) <- c("gene", "cell_type")
  if (anyDuplicated(marker_table$gene)) {
    stop(
      "Curated marker genes must be unique across cell types.",
      call. = FALSE
    )
  }
  marker_table$cell_type_label <- unname(cell_type_marker_labels[
    marker_table$cell_type
  ])
  marker_table$marker_source <- "cell_type_marker_genes"
  # ANALYSIS_OK[R005]: restrict curated markers to the validated gene universe before overlap reporting.
  marker_table <- marker_table[
    marker_table$gene %in% gene_universe,
    ,
    drop = FALSE
  ]
  if (cdkn1b_gene %in% gene_universe) {
    marker_table <- rbind(
      marker_table,
      data.frame(
        gene = cdkn1b_gene,
        cell_type = "cdkn1b_standalone",
        cell_type_label = cdkn1b_gene,
        marker_source = "standalone_gene",
        stringsAsFactors = FALSE
      )
    )
  }
  marker_table[
    order(
      marker_table$marker_source,
      marker_table$cell_type,
      marker_table$gene
    ),
    ,
    drop = FALSE
  ]
}

# ANALYSIS_OK[R026]: local marker-overlap writer is called by the enrichment workflow.
write_marker_overlap <- function(result_table, gene_universe, path) {
  marker_table <- make_marker_table(gene_universe)
  if (anyDuplicated(marker_table$gene) || anyDuplicated(result_table$gene)) {
    stop(
      "Marker and DE gene keys must be unique before overlap.",
      call. = FALSE
    )
  }
  overlap <- merge(
    marker_table,
    result_table,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )
  if (nrow(overlap) != nrow(marker_table)) {
    stop(
      "Marker overlap changed curated marker row cardinality.",
      call. = FALSE
    )
  }
  overlap$significant <- !is.na(overlap$padj) &
    overlap$padj < de_significance_cutoff
  write_tsv(overlap, path)
  invisible(overlap)
}

# ---- primary and paired models ----

sample_summary <- sample_table[, c(
  "sample_id",
  "Mouse",
  "condition_label",
  "condition",
  "constructed_sample_id",
  "sample_id_source",
  "analysis_unit",
  "paired_mouse",
  "n_cells",
  "total_counts"
)]
write_tsv(sample_summary, file.path(deg_dir, "pseudobulk_sample_summary.tsv"))

primary_de <- run_deseq2(
  pseudobulk_counts,
  sample_table,
  ~condition,
  "primary_unpaired_condition"
)
full_de <- primary_de$result_table
sig_de <- full_de[
  !is.na(full_de$padj) & full_de$padj < de_significance_cutoff,
  ,
  drop = FALSE
]
write_tsv(full_de, file.path(deg_dir, "deseq2_full_results.tsv"))
write_tsv(sig_de, file.path(deg_dir, "deseq2_significant_degs.tsv"))
de_marker_overlap <- write_marker_overlap(
  full_de,
  primary_de$tested_genes,
  file.path(deg_dir, "deseq2_marker_overlap.tsv")
)

paired_status <- "skipped"
paired_reason <- "fewer than two mice have both control and estim samples"
paired_de_n_degs <- NA_integer_
if (length(paired_mice) >= min_paired_mice) {
  paired_sample_table <- sample_table[
    sample_table$Mouse %in% paired_mice,
    ,
    drop = FALSE
  ]
  paired_sample_table$mouse <- droplevels(factor(paired_sample_table$Mouse))
  paired_sample_table$condition <- droplevels(paired_sample_table$condition)
  paired_counts <- pseudobulk_counts[,
    rownames(paired_sample_table),
    drop = FALSE
  ]
  paired_design <- stats::model.matrix(
    ~ mouse + condition,
    data = paired_sample_table
  )
  if (qr(paired_design)$rank == ncol(paired_design)) {
    paired_de <- run_deseq2(
      paired_counts,
      paired_sample_table,
      ~ mouse + condition,
      "paired_mouse_condition_sensitivity"
    )
    paired_full_de <- paired_de$result_table
    paired_sig_de <- paired_full_de[
      !is.na(paired_full_de$padj) &
        paired_full_de$padj < de_significance_cutoff,
      ,
      drop = FALSE
    ]
    write_tsv(
      paired_full_de,
      file.path(deg_dir, "deseq2_paired_sensitivity_full_results.tsv")
    )
    write_tsv(
      paired_sig_de,
      file.path(deg_dir, "deseq2_paired_sensitivity_significant_degs.tsv")
    )
    write_marker_overlap(
      paired_full_de,
      paired_de$tested_genes,
      file.path(deg_dir, "deseq2_paired_sensitivity_marker_overlap.tsv")
    )
    paired_status <- "run"
    paired_reason <- "paired sensitivity used mice with both control and estim samples"
    paired_de_n_degs <- nrow(paired_sig_de)
  } else {
    paired_reason <- paste0(
      "paired sensitivity design was not full rank; columns=",
      paste(colnames(paired_design), collapse = ",")
    )
  }
}
if (!identical(paired_status, "run")) {
  write_reason_tsv(
    file.path(deg_dir, "deseq2_paired_sensitivity_full_results.tsv"),
    paired_reason,
    nrow(pseudobulk_counts),
    0L,
    0L
  )
  write_reason_tsv(
    file.path(deg_dir, "deseq2_paired_sensitivity_significant_degs.tsv"),
    paired_reason,
    nrow(pseudobulk_counts),
    0L,
    0L
  )
  write_reason_tsv(
    file.path(deg_dir, "deseq2_paired_sensitivity_marker_overlap.tsv"),
    paired_reason,
    nrow(pseudobulk_counts),
    0L,
    0L
  )
}

# ---- DE volcano ----

volcano_data <- full_de[,
  c("gene", "log2FoldChange", "pvalue", "padj"),
  drop = FALSE
]
# ANALYSIS_OK[R005]: restrict volcano rows to finite reportable statistics for visualization only.
volcano_data <- volcano_data[
  !is.na(volcano_data$gene) &
    nzchar(as.character(volcano_data$gene)) &
    is.finite(volcano_data$log2FoldChange) &
    is.finite(volcano_data$pvalue) &
    is.finite(volcano_data$padj),
  ,
  drop = FALSE
]
if (nrow(volcano_data) == 0L) {
  stop("No finite DE rows available for volcano plot.", call. = FALSE)
}
volcano_data$neg_log10_padj <- -log10(pmax(
  volcano_data$padj,
  .Machine$double.xmin
))
volcano_data$significance <- "Not significant"
significant <- volcano_data$padj < de_significance_cutoff
volcano_data$significance[
  significant & volcano_data$log2FoldChange > 0
] <- "Increased"
volcano_data$significance[
  significant & volcano_data$log2FoldChange < 0
] <- "Decreased"
volcano_data$significance <- factor(
  volcano_data$significance,
  levels = c("Not significant", "Increased", "Decreased")
)
volcano_data$label <- NA_character_
for (indices in list(
  which(volcano_data$log2FoldChange > 0),
  which(volcano_data$log2FoldChange < 0)
)) {
  if (length(indices) == 0L) {
    next
  }
  label_indices <- utils::head(
    order(
      volcano_data$pvalue[indices],
      -abs(volcano_data$log2FoldChange[indices]),
      as.character(volcano_data$gene[indices])
    ),
    volcano_label_limit
  )
  volcano_data$label[indices[
    label_indices
  ]] <- as.character(volcano_data$gene[indices[label_indices]])
}
label_data <- volcano_data[!is.na(volcano_data$label), , drop = FALSE]
volcano <- ggplot2::ggplot(
  volcano_data,
  ggplot2::aes(x = .data[["log2FoldChange"]], y = .data[["neg_log10_padj"]])
) +
  ggplot2::geom_hline(
    yintercept = -log10(de_significance_cutoff),
    linewidth = 0.25,
    color = "grey70"
  ) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "grey70") +
  ggplot2::geom_point(
    ggplot2::aes(color = .data[["significance"]]),
    alpha = 0.55,
    size = 0.8
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Not significant" = config$palettes$analysis[["mid"]],
      "Increased" = config$palettes$analysis[["high"]],
      "Decreased" = config$palettes$analysis[["low"]]
    ),
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "MG-selected differential expression",
    x = sprintf("Shrunken log2 FC %s", conditions$contrast_display),
    y = expression(-log[10]("adjusted p-value")),
    color = NULL
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
if (nrow(label_data) > 0L) {
  volcano <- volcano +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = .data[["label"]]),
      size = 2.2,
      seed = 275,
      box.padding = 0.4,
      point.padding = 0.2,
      force = 2,
      max.time = 5,
      max.overlaps = Inf,
      show.legend = FALSE
    )
}
save_publication_plot(
  volcano,
  file.path(figure_dir, "mg_selected_de_volcano"),
  width = 7,
  height = 6.5,
  notebook_basename = "mg_selected_de_volcano.png",
  dpi = 300
)

# ---- design summary ----

design_summary <- data.frame(
  analysis = c("deseq2_primary", "deseq2_paired_sensitivity"),
  design = c("~ condition", "~ mouse + condition"),
  method = c("deseq2_wald", "deseq2_wald"),
  status = c("run", paired_status),
  contrast = contrast_direction,
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
  cluster_column = cluster_column,
  counts_layer = counts_layer,
  lfc_shrink_type = lfc_shrink_type,
  gsea_seed = seed,
  stringsAsFactors = FALSE
)
write_tsv(design_summary, file.path(deg_dir, "design_summary.tsv"))

# ---- GO BP ORA, GSEA, Bayesian ORA, and simplification ----

# ANALYSIS_OK[R026]: local gene-mapping helper is called by all ORA and GSEA branches.
map_genes_to_entrez <- function(symbols) {
  # ANALYSIS_OK[R005]: remove missing and duplicate gene symbols before deterministic mapping.
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (length(symbols) == 0L) {
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  }
  # ANALYSIS_OK[warning-suppression]: expected unmapped-symbol messages are suppressed after input validation.
  mapped <- suppressMessages(clusterProfiler::bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db::org.Mm.eg.db
  ))
  mapped[
    !is.na(mapped$ENTREZID) & nzchar(mapped$ENTREZID),
    ,
    drop = FALSE
  ] |>
    unique()
}

# ANALYSIS_OK[R026]: local GO ORA helper is called for both DE directions.
run_go_ora <- function(significant_symbols, background_map, direction, path) {
  significant_map <- background_map[
    background_map$SYMBOL %in% significant_symbols,
    ,
    drop = FALSE
  ]
  significant_entrez <- unique(significant_map$ENTREZID)
  background_entrez <- unique(background_map$ENTREZID)
  if (length(significant_entrez) < min_ora_genes) {
    write_reason_tsv(
      path,
      paste0(
        "fewer than ",
        min_ora_genes,
        " mapped significant ",
        direction,
        " DE genes"
      ),
      length(significant_symbols),
      length(significant_entrez),
      min_ora_genes
    )
    return(NULL)
  }
  # ANALYSIS_OK[warning-suppression]: expected GO ORA diagnostics are suppressed after validated inputs.
  enrichment <- suppressMessages(clusterProfiler::enrichGO(
    gene = significant_entrez,
    universe = background_entrez,
    OrgDb = org.Mm.eg.db::org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  ))
  enrichment_table <- as.data.frame(enrichment)
  if (nrow(enrichment_table) == 0L) {
    write_reason_tsv(
      path,
      paste0("no GO BP ORA terms returned for ", direction, " DE genes"),
      length(significant_symbols),
      length(significant_entrez),
      min_ora_genes
    )
    return(NULL)
  }
  enrichment_table$direction <- direction
  write_tsv(enrichment_table, path)
  enrichment
}

# ANALYSIS_OK[R026]: local GO GSEA helper is called by the enrichment workflow.
run_go_gsea <- function(result_table, background_map, path, mapping_path) {
  ranked <- result_table[
    !is.na(result_table$stat),
    c("gene", "stat"),
    drop = FALSE
  ]
  if (anyDuplicated(background_map$SYMBOL)) {
    stop(
      "GO GSEA background map contains duplicate SYMBOL keys.",
      call. = FALSE
    )
  }
  ranked_before_merge <- nrow(ranked)
  ranked <- merge(
    ranked,
    background_map,
    by.x = "gene",
    by.y = "SYMBOL",
    all = FALSE,
    sort = FALSE
  )
  if (nrow(ranked) > ranked_before_merge) {
    stop("GO GSEA symbol mapping changed row cardinality.", call. = FALSE)
  }
  if (nrow(ranked) == 0L) {
    write_reason_tsv(
      path,
      "no DESeq2 ranked statistics mapped to Entrez IDs",
      nrow(result_table),
      0L,
      min_gsea_genes
    )
    return(NULL)
  }
  ranked$abs_stat <- abs(ranked$stat)
  # ANALYSIS_OK[R005]: retain one deterministic Entrez representative per ranked gene identifier.
  ranked <- ranked[
    order(ranked$ENTREZID, -ranked$abs_stat, ranked$gene),
    ,
    drop = FALSE
  ]
  ranked$selected_for_gsea <- !duplicated(ranked$ENTREZID)
  write_tsv(
    ranked[, c("gene", "ENTREZID", "stat", "selected_for_gsea")],
    mapping_path
  )
  # ANALYSIS_OK[R005]: retain only the selected unique Entrez representatives for GSEA.
  ranked <- ranked[ranked$selected_for_gsea, , drop = FALSE]
  gene_list <- ranked$stat
  names(gene_list) <- ranked$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  if (length(gene_list) < min_gsea_genes) {
    write_reason_tsv(
      path,
      paste0("fewer than ", min_gsea_genes, " mapped ranked genes"),
      nrow(result_table),
      length(gene_list),
      min_gsea_genes
    )
    return(NULL)
  }
  # ANALYSIS_OK[warning-suppression]: expected GO GSEA diagnostics are suppressed after validated inputs.
  gsea <- suppressMessages(clusterProfiler::gseGO(
    geneList = gene_list,
    OrgDb = org.Mm.eg.db::org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = gsea_min_gene_set_size,
    maxGSSize = gsea_max_gene_set_size,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  ))
  gsea_table <- as.data.frame(gsea)
  if (nrow(gsea_table) == 0L) {
    write_reason_tsv(
      path,
      "no GO BP GSEA terms returned",
      nrow(result_table),
      length(gene_list),
      min_gsea_genes
    )
    return(NULL)
  }
  write_tsv(gsea_table, path)
  gsea
}

# ANALYSIS_OK[R026]: local enrichment-result filter is called by simplification branches.
significant_enrichresult <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  result <- x@result
  keep <- !is.na(result$p.adjust) &
    result$p.adjust < enrichment_significance_cutoff
  if (!any(keep)) {
    return(NULL)
  }
  x@result <- result[keep, , drop = FALSE]
  x
}

# ANALYSIS_OK[R026]: local GO simplification writer is called by each enrichment result branch.
write_simplified <- function(x, path, label) {
  if (is.null(x)) {
    write_reason_tsv(
      path,
      paste0("no significant terms to simplify for ", label),
      0L,
      0L,
      1L
    )
    return(NULL)
  }
  simplified <- clusterProfiler::simplify(x)
  table <- as.data.frame(simplified)
  if (nrow(table) == 0L) {
    write_reason_tsv(
      path,
      paste0("no simplified terms for ", label),
      nrow(x@result),
      0L,
      1L
    )
    return(NULL)
  }
  write_tsv(table, path)
  simplified
}

# ANALYSIS_OK[R026]: local Bayesian simplification writer is called by Bayesian branches.
write_bayes_simplified <- function(x, path, label) {
  if (is.null(x)) {
    write_reason_tsv(
      path,
      paste0("no significant terms for bayes+simplify: ", label),
      0L,
      0L,
      1L
    )
    return(NULL)
  }
  bayes <- enrichit::bayes_enrich(x, seed = seed)
  simplified <- clusterProfiler::simplify(bayes)
  table <- as.data.frame(simplified)
  if (nrow(table) == 0L) {
    write_reason_tsv(
      path,
      paste0("no bayes+simplified terms for ", label),
      nrow(x@result),
      0L,
      1L
    )
    return(NULL)
  }
  write_tsv(table, path)
  simplified
}

background_map <- map_genes_to_entrez(primary_de$tested_genes)
if (nrow(background_map) == 0L) {
  stop(
    "No tested DESeq2 genes mapped to Entrez IDs for enrichment background.",
    call. = FALSE
  )
}
up_genes <- sig_de$gene[
  !is.na(sig_de$log2FoldChange) & sig_de$log2FoldChange > 0
]
down_genes <- sig_de$gene[
  !is.na(sig_de$log2FoldChange) & sig_de$log2FoldChange < 0
]
ora_up <- run_go_ora(
  up_genes,
  background_map,
  "up_estim_vs_control",
  file.path(enrichment_dir, "go_bp_ora_up.tsv")
)
ora_down <- run_go_ora(
  down_genes,
  background_map,
  "down_estim_vs_control",
  file.path(enrichment_dir, "go_bp_ora_down.tsv")
)
set.seed(seed)
gsea <- run_go_gsea(
  full_de,
  background_map,
  file.path(enrichment_dir, "go_bp_gsea.tsv"),
  file.path(enrichment_dir, "go_bp_gsea_symbol_entrez_mapping.tsv")
)

s_up <- write_simplified(
  significant_enrichresult(ora_up),
  file.path(enrichment_dir, "go_bp_ora_up_simplified.tsv"),
  "GO ORA up"
)
s_down <- write_simplified(
  significant_enrichresult(ora_down),
  file.path(enrichment_dir, "go_bp_ora_down_simplified.tsv"),
  "GO ORA down"
)
s_gsea <- write_simplified(
  significant_enrichresult(gsea),
  file.path(enrichment_dir, "go_bp_gsea_simplified.tsv"),
  "GO GSEA"
)
b_up <- write_bayes_simplified(
  significant_enrichresult(ora_up),
  file.path(enrichment_dir, "go_bp_ora_up_bayes_simplified.tsv"),
  "GO ORA up"
)
b_down <- write_bayes_simplified(
  significant_enrichresult(ora_down),
  file.path(enrichment_dir, "go_bp_ora_down_bayes_simplified.tsv"),
  "GO ORA down"
)

# ---- enrichment dotplots ----

if (is.null(s_up) || nrow(as.data.frame(s_up)) == 0L) {
  plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = "No enrichment terms: GO BP ORA up (simplified)"
    ) +
    ggplot2::theme_void()
} else {
  plot <- enrichplot::dotplot(s_up, showCategory = enrichment_dotplot_show_n) +
    ggplot2::ggtitle("GO BP ORA up (simplified)")
}
save_publication_plot(
  plot,
  file.path(figure_dir, "mg_selected_go_ora_up_dotplot"),
  8,
  7,
  "mg_selected_go_ora_up_dotplot.png"
)

if (is.null(s_down) || nrow(as.data.frame(s_down)) == 0L) {
  plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = "No enrichment terms: GO BP ORA down (simplified)"
    ) +
    ggplot2::theme_void()
} else {
  plot <- enrichplot::dotplot(
    s_down,
    showCategory = enrichment_dotplot_show_n
  ) +
    ggplot2::ggtitle("GO BP ORA down (simplified)")
}
save_publication_plot(
  plot,
  file.path(figure_dir, "mg_selected_go_ora_down_dotplot"),
  8,
  7,
  "mg_selected_go_ora_down_dotplot.png"
)

if (is.null(s_gsea) || nrow(as.data.frame(s_gsea)) == 0L) {
  plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = "No enrichment terms: GO BP GSEA (simplified)"
    ) +
    ggplot2::theme_void()
} else {
  plot <- enrichplot::dotplot(
    s_gsea,
    showCategory = enrichment_dotplot_show_n,
    split = ".sign"
  ) +
    ggplot2::facet_grid(. ~ .sign) +
    ggplot2::ggtitle("GO BP GSEA (simplified)")
}
save_publication_plot(
  plot,
  file.path(figure_dir, "mg_selected_go_gsea_dotplot"),
  8,
  7,
  "mg_selected_go_gsea_dotplot.png"
)

if (is.null(b_up) || nrow(as.data.frame(b_up)) == 0L) {
  plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = "No enrichment terms: GO BP ORA up (bayes + simplified)"
    ) +
    ggplot2::theme_void()
} else {
  plot <- enrichplot::dotplot(b_up, showCategory = enrichment_dotplot_show_n) +
    ggplot2::ggtitle("GO BP ORA up (bayes + simplified)")
}
save_publication_plot(
  plot,
  file.path(figure_dir, "mg_selected_go_ora_up_bayes_dotplot"),
  8,
  7,
  "mg_selected_go_ora_up_bayes_dotplot.png"
)

if (is.null(b_down) || nrow(as.data.frame(b_down)) == 0L) {
  plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = "No enrichment terms: GO BP ORA down (bayes + simplified)"
    ) +
    ggplot2::theme_void()
} else {
  plot <- enrichplot::dotplot(
    b_down,
    showCategory = enrichment_dotplot_show_n
  ) +
    ggplot2::ggtitle("GO BP ORA down (bayes + simplified)")
}
save_publication_plot(
  plot,
  file.path(figure_dir, "mg_selected_go_ora_down_bayes_dotplot"),
  8,
  7,
  "mg_selected_go_ora_down_bayes_dotplot.png"
)

# ---- reportable values ----

# ANALYSIS_OK[R026]: local report-number serializer is called for the final numbers artifact.
write_numbers_json <- function(values, path) {
  keys <- names(values)
  encoded <- vapply(
    seq_along(values),
    function(i) {
      value <- values[[i]]
      value_text <- if (length(value) != 1L || is.na(value)) {
        "null"
      } else if (is.numeric(value) || is.integer(value)) {
        as.character(value)
      } else if (is.logical(value)) {
        if (isTRUE(value)) "true" else "false"
      } else {
        value <- gsub("\\\\", "\\\\\\\\", as.character(value))
        value <- gsub("\"", "\\\\\"", value)
        value <- gsub("\n", "\\\\n", value)
        paste0("\"", value, "\"")
      }
      comma <- if (i < length(values)) "," else ""
      paste0("  \"", keys[[i]], "\": ", value_text, comma)
    },
    character(1)
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c("{", encoded, "}"), con = path)
}

write_numbers_json(
  list(
    n_samples = nrow(sample_table),
    n_cells = ncol(counts),
    n_tested_genes = length(primary_de$tested_genes),
    n_degs = nrow(sig_de),
    n_marker_degs = sum(de_marker_overlap$significant, na.rm = TRUE),
    n_unmatched_mice = length(unique(unmatched_mice)),
    source_input = input_path,
    cluster_column = cluster_column,
    counts_layer = counts_layer,
    gsea_seed = seed,
    lfc_shrink_type = primary_de$shrink_type,
    primary_design = "~ condition",
    paired_sensitivity_design = "~ mouse + condition",
    paired_sensitivity_included_mice = paste(paired_mice, collapse = ","),
    paired_sensitivity_status = paired_status,
    paired_sensitivity_n_degs = paired_de_n_degs
  ),
  file.path(deg_dir, "numbers.json")
)

message("Saved MG-selected DE outputs under ", deg_dir)
message("Saved MG-selected enrichment outputs under ", enrichment_dir)
