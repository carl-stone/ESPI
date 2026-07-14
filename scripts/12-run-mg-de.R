#!/usr/bin/env Rscript

# Run mg-selected pseudobulk differential expression, muscat edgeR_NB_optim
# differential detection, and enrichment follow-up analyses.
#
# Usage:
#   Rscript scripts/12-run-mg-de.R \
#     --input <mg-selected-clustered-seurat-object.rds> \
#     --cluster-column <cluster metadata column> \
#     --condition-col <condition metadata column> \
#     --counts-layer <raw counts layer> \
#     --lfc-shrink-type <normal|apeglm> \
#     --overwrite
# Defaults target the mg-selected no-filter-cc pFlog branch:
#   CURRENT_OBJECT_DIR/cluster_pflog_mg_selected_no_filter_cc_elbow20.rds
#   cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5
#
#
# Differential detection uses `muscat::pbDS(method = "DD")`, equivalent to
# `muscat::pbDD()`, on pseudobulk detected-cell counts with the published
# edgeR_NB_optim workflow.
# Outputs:
#   DEG_DIR/mg_selected/pseudobulk_sample_summary.tsv
#   DEG_DIR/mg_selected/design_summary.tsv
#   DEG_DIR/mg_selected/deseq2_full_results.tsv
#   DEG_DIR/mg_selected/deseq2_significant_degs.tsv
#   DEG_DIR/mg_selected/deseq2_marker_overlap.tsv
#   DEG_DIR/mg_selected/detection_full_results.tsv
#   DEG_DIR/mg_selected/detection_marker_overlap.tsv
#   DEG_DIR/mg_selected/*paired_sensitivity*.tsv when feasible, or explicit
#     skipped TSVs when not feasible
#   ENRICHMENT_DIR/mg_selected/go_bp_ora_{up,down}.tsv
#   ENRICHMENT_DIR/mg_selected/go_bp_gsea.tsv
#   ENRICHMENT_DIR/mg_selected/go_bp_gsea_symbol_entrez_mapping.tsv
#   FIGURE_DIR/mg_selected/mg_selected_de_dd_effect_scatter.(png|pdf)
#   notebook/figures/mg_selected_de_dd_effect_scatter.png symlink

suppressPackageStartupMessages({
  library(here)
})
here::i_am("scripts/12-run-mg-de.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})
utils::data("cell_type_marker_genes", package = "ESPI", envir = environment())
utils::data("cell_type_marker_labels", package = "ESPI", envir = environment())
cell_type_marker_genes <- get("cell_type_marker_genes", envir = environment())
cell_type_marker_labels <- get("cell_type_marker_labels", envir = environment())
palette_analysis_three <- get(
  "palette_analysis_three",
  envir = asNamespace("ESPI"),
  inherits = FALSE
)
CONTRAST_DISPLAY_LABEL <- get(
  "CONTRAST_DISPLAY_LABEL",
  envir = asNamespace("ESPI"),
  inherits = FALSE
)

# ---- parameters ----

get_arg <- function(args, flag, default = NULL) {
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
allowed_flags <- c(
  "--input",
  "--cluster-column",
  "--condition-col",
  "--control-label",
  "--estim-label",
  "--counts-layer",
  "--deg-dir",
  "--enrichment-dir",
  "--lfc-shrink-type",
  "--overwrite"
)
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

required_packages <- c(
  "DESeq2",
  "muscat",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "clusterProfiler",
  "org.Mm.eg.db",
  "ggplot2",
  "Matrix",
  "SeuratObject"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

input_path <- get_arg(
  cli_args,
  "--input",
  file.path(
    CURRENT_OBJECT_DIR,
    "cluster_pflog_mg_selected_no_filter_cc_elbow20.rds"
  )
)
cluster_column <- get_arg(
  cli_args,
  "--cluster-column",
  "cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5"
)
condition_col <- get_arg(cli_args, "--condition-col", CONDITION_COL)
control_label <- get_arg(cli_args, "--control-label", CTRL_LABEL)
estim_label <- get_arg(cli_args, "--estim-label", ESTIM_LABEL)
counts_layer <- get_arg(cli_args, "--counts-layer", "counts")
deg_dir <- get_arg(cli_args, "--deg-dir", file.path(DEG_DIR, "mg_selected"))
enrichment_dir <- get_arg(
  cli_args,
  "--enrichment-dir",
  file.path(ENRICHMENT_DIR, "mg_selected")
)
figure_dir <- file.path(FIGURE_DIR, "mg_selected")
de_dd_effect_scatter_png_path <- file.path(
  figure_dir,
  "mg_selected_de_dd_effect_scatter.png"
)
de_dd_effect_scatter_pdf_path <- file.path(
  figure_dir,
  "mg_selected_de_dd_effect_scatter.pdf"
)
de_dd_effect_scatter_notebook_png_path <- here::here(
  "notebook",
  "figures",
  basename(de_dd_effect_scatter_png_path)
)
lfc_shrink_type <- get_arg(cli_args, "--lfc-shrink-type", "normal")
overwrite_outputs <- "--overwrite" %in% cli_args
if (!lfc_shrink_type %in% c("normal", "apeglm")) {
  stop(
    "--lfc-shrink-type must be one of normal or apeglm; got ",
    lfc_shrink_type,
    call. = FALSE
  )
}
if (
  identical(lfc_shrink_type, "apeglm") &&
    !requireNamespace("apeglm", quietly = TRUE)
) {
  stop("--lfc-shrink-type apeglm requires package apeglm.", call. = FALSE)
}

MIN_PSEUDOBULK_GENE_COUNT <- 10L
MIN_ORA_GENES <- 5L
MIN_GSEA_GENES <- 10L
CONTROL_LEVEL <- "control"
ESTIM_LEVEL <- "estim"
CONTRAST_DIRECTION <- "estim_vs_control"
CDKN1B_GENE <- "Cdkn1b"
GSEA_SEED <- SEED

# ---- helpers ----

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

safe_component <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "missing")
}


assert_no_missing_metadata <- function(meta, columns) {
  missing_columns <- setdiff(columns, colnames(meta))
  if (length(missing_columns) > 0) {
    stop(
      "Missing required metadata column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  columns_with_missing <- columns[vapply(
    columns,
    function(column) {
      values <- meta[[column]]
      any(is.na(values) | trimws(as.character(values)) == "")
    },
    logical(1)
  )]
  if (length(columns_with_missing) > 0) {
    stop(
      "Missing values in required metadata column(s): ",
      paste(columns_with_missing, collapse = ", "),
      call. = FALSE
    )
  }
}

assert_full_rank <- function(design, label) {
  design_rank <- qr(design)$rank
  design_cols <- ncol(design)
  if (design_rank < design_cols) {
    stop(
      label,
      " design matrix is not full rank (rank ",
      design_rank,
      " < ",
      design_cols,
      "). Columns: ",
      paste(colnames(design), collapse = ", "),
      call. = FALSE
    )
  }
}

assert_output_paths_clear <- function(paths, overwrite) {
  link_targets <- Sys.readlink(paths)
  existing_paths <- paths[
    file.exists(paths) | (!is.na(link_targets) & nzchar(link_targets))
  ]
  if (length(existing_paths) > 0L && !isTRUE(overwrite)) {
    stop(
      "Refusing to overwrite existing output file(s) without --overwrite: ",
      paste(existing_paths, collapse = ", "),
      call. = FALSE
    )
  }
}


json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"", "\\\\\"", x)
  x <- gsub("\n", "\\\\n", x)
  x
}

write_numbers_json <- function(values, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  keys <- names(values)
  body <- vapply(
    seq_along(values),
    function(i) {
      value <- values[[i]]
      encoded_value <- if (length(value) != 1L || is.na(value)) {
        "null"
      } else if (is.numeric(value) || is.integer(value)) {
        as.character(value)
      } else if (is.logical(value)) {
        if (isTRUE(value)) "true" else "false"
      } else {
        paste0("\"", json_escape(as.character(value)), "\"")
      }
      comma <- if (i < length(values)) "," else ""
      paste0("  \"", json_escape(keys[[i]]), "\": ", encoded_value, comma)
    },
    character(1)
  )
  writeLines(c("{", body, "}"), con = path)
}

register_or_write_numbers <- function(values, out_dir) {
  write_numbers_json(values, file.path(out_dir, "numbers.json"))
  register_script <- file.path(
    here::here(),
    "skills",
    "core",
    "scripts",
    "register_value.py"
  )
  if (file.exists(register_script)) {
    for (key in names(values)) {
      output <- system2(
        "python3",
        args = c(
          register_script,
          "--name",
          paste0("mg_selected_de.", key),
          "--value",
          as.character(values[[key]])
        ),
        stdout = TRUE,
        stderr = TRUE
      )
      status <- attr(output, "status")
      if (!is.null(status) && status != 0) {
        stop(
          "Failed to register reportable value '",
          key,
          "' with ",
          register_script,
          ": ",
          paste(output, collapse = "\n"),
          call. = FALSE
        )
      }
    }
  }
}

make_marker_table <- function(gene_universe) {
  if (
    !identical(names(cell_type_marker_genes), names(cell_type_marker_labels))
  ) {
    stop(
      "cell_type_marker_genes and cell_type_marker_labels must have identical names.",
      call. = FALSE
    )
  }
  marker_table <- stack(cell_type_marker_genes)
  colnames(marker_table) <- c("gene", "cell_type")
  duplicated_markers <- marker_table$gene[duplicated(marker_table$gene)]
  if (length(duplicated_markers) > 0) {
    stop(
      "Marker gene(s) assigned to more than one cell type: ",
      paste(unique(duplicated_markers), collapse = ", "),
      call. = FALSE
    )
  }
  marker_table$cell_type_label <- unname(
    cell_type_marker_labels[marker_table$cell_type]
  )
  marker_table$marker_source <- "cell_type_marker_genes"
  marker_table <- marker_table[
    marker_table$gene %in% gene_universe,
    ,
    drop = FALSE
  ]
  if (CDKN1B_GENE %in% gene_universe) {
    marker_table <- rbind(
      marker_table,
      data.frame(
        gene = CDKN1B_GENE,
        cell_type = "cdkn1b_standalone",
        cell_type_label = CDKN1B_GENE,
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
  ]
}

write_marker_overlap <- function(
  result_table,
  gene_universe,
  path,
  padj_column
) {
  marker_table <- make_marker_table(gene_universe)
  overlap <- merge(
    marker_table,
    result_table,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )
  overlap$significant <- !is.na(overlap[[padj_column]]) &
    overlap[[padj_column]] < 0.05
  write_tsv(overlap, path)
  invisible(overlap)
}

build_de_dd_effect_data <- function(de_table, detection_table, design_label) {
  de_effects <- de_table[, c("gene", "log2FoldChange", "padj"), drop = FALSE]
  colnames(de_effects) <- c("gene", "de_log2_fold_change", "de_padj")
  detection_effects <- detection_table[,
    c("gene", "logFC", "padj"),
    drop = FALSE
  ]
  colnames(detection_effects) <- c("gene", "dd_log2_fold_change", "dd_padj")

  joined <- merge(
    de_effects,
    detection_effects,
    by = "gene",
    all = FALSE,
    sort = FALSE
  )
  joined <- joined[
    !is.na(joined$de_log2_fold_change) &
      !is.na(joined$dd_log2_fold_change),
    ,
    drop = FALSE
  ]

  marker_genes <- unique(unlist(cell_type_marker_genes, use.names = FALSE))
  joined$curated_marker <- joined$gene %in% marker_genes

  de_significant <- !is.na(joined$de_padj) & joined$de_padj < 0.05
  dd_significant <- !is.na(joined$dd_padj) & joined$dd_padj < 0.05
  joined$fdr_category <- ifelse(
    de_significant & dd_significant,
    "both",
    ifelse(
      de_significant,
      "DE only",
      ifelse(dd_significant, "DD only", "neither")
    )
  )
  joined$fdr_category <- factor(
    joined$fdr_category,
    levels = c("neither", "DE only", "DD only", "both")
  )
  joined$design <- factor(
    design_label,
    levels = c("All samples by condition", "Paired mice by condition")
  )
  joined
}

plot_de_dd_effect_scatter <- function(plot_data) {
  if (nrow(plot_data) == 0L) {
    stop(
      "No genes had both DESeq2 and muscat DD effects after joining on gene.",
      call. = FALSE
    )
  }

  plot_data <- plot_data[
    order(plot_data$fdr_category != "neither", plot_data$fdr_category),
    ,
    drop = FALSE
  ]
  label_data <- plot_data[
    plot_data$curated_marker & plot_data$fdr_category != "neither",
    ,
    drop = FALSE
  ]

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[["de_log2_fold_change"]],
      y = .data[["dd_log2_fold_change"]]
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "grey70") +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[["fdr_category"]]),
      alpha = 0.55,
      size = 0.8
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "neither" = unname(palette_analysis_three[["mid"]]),
        "DE only" = unname(palette_analysis_three[["low"]]),
        "DD only" = unname(palette_analysis_three[["high"]]),
        "both" = "#4daf4a"
      ),
      name = "FDR < 0.05",
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "MG-selected DE and differential detection effects",
      subtitle = "Inner join on gene; genes missing either effect are omitted.",
      x = sprintf("Differential expression log2 FC %s", CONTRAST_DISPLAY_LABEL),
      y = sprintf("Differential detection log2 FC %s", CONTRAST_DISPLAY_LABEL)
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )

  if (nrow(label_data) > 0L) {
    plot <- plot +
      ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(label = .data[["gene"]]),
        size = 2.2,
        seed = 275,
        show.legend = FALSE
      )
  }

  if (length(unique(plot_data$design)) > 1L) {
    plot <- plot + ggplot2::facet_wrap(stats::as.formula("~ design"))
  }

  plot
}

results_to_table <- function(
  unshrunk_result,
  shrunk_result,
  count_matrix,
  sample_table,
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
  control_samples <- rownames(sample_table)[
    sample_table$condition == CONTROL_LEVEL
  ]
  estim_samples <- rownames(sample_table)[sample_table$condition == ESTIM_LEVEL]
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
    contrast = CONTRAST_DIRECTION,
    design = design_label,
    lfc_shrink_type = shrink_type,
    stringsAsFactors = FALSE
  )
}

run_deseq2 <- function(
  count_matrix,
  sample_table,
  design_formula,
  design_label,
  shrink_type
) {
  design <- stats::model.matrix(design_formula, data = sample_table)
  assert_full_rank(design, paste("DESeq2", design_label))

  keep <- Matrix::rowSums(count_matrix) >= MIN_PSEUDOBULK_GENE_COUNT
  if (!any(keep)) {
    stop(
      "No genes passed DESeq2 pre-filter rowSums >= ",
      MIN_PSEUDOBULK_GENE_COUNT,
      ".",
      call. = FALSE
    )
  }
  filtered_counts <- count_matrix[keep, , drop = FALSE]
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = sample_table,
    design = design_formula
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  coef_name <- "condition_estim_vs_control"
  available_results <- DESeq2::resultsNames(dds)
  if (!coef_name %in% available_results) {
    stop(
      "Expected DESeq2 coefficient not found for estim vs control: ",
      coef_name,
      ". Available coefficient(s): ",
      paste(available_results, collapse = ", "),
      call. = FALSE
    )
  }
  unshrunk_result <- DESeq2::results(
    dds,
    contrast = c("condition", ESTIM_LEVEL, CONTROL_LEVEL)
  )
  if (
    identical(shrink_type, "apeglm") &&
      !requireNamespace("apeglm", quietly = TRUE)
  ) {
    stop("DESeq2 shrinkage type apeglm requires package apeglm.", call. = FALSE)
  }
  shrunk_result <- suppressMessages(DESeq2::lfcShrink(
    dds,
    coef = coef_name,
    res = unshrunk_result,
    type = shrink_type,
    quiet = TRUE
  ))
  result_table <- results_to_table(
    unshrunk_result = unshrunk_result,
    shrunk_result = shrunk_result,
    count_matrix = filtered_counts,
    sample_table = sample_table,
    shrink_type = shrink_type,
    design_label = design_label
  )
  list(
    result_table = result_table,
    tested_genes = rownames(filtered_counts),
    shrink_type = shrink_type,
    design_matrix = design
  )
}

run_detection_muscat_dd <- function(
  counts,
  cell_metadata,
  sample_table,
  detection_fraction,
  design_formula,
  design_label
) {
  # counts: raw gene x cell dgCMatrix / Matrix, colnames match rownames(cell_metadata).
  # cell_metadata: data.frame with rownames = cell barcodes; must carry the columns
  #   `pseudobulk_sample_id` (character; identifies the Mouse x Condition sample) and
  #   `condition` (factor with levels c(CONTROL_LEVEL, ESTIM_LEVEL)); may carry `mouse`.
  # sample_table: rownames = pseudobulk sample IDs used as columns of DD outputs;
  #   columns must include `condition` (factor), and for the paired design `mouse`.
  # detection_fraction: gene x sample matrix used to fill mean detected fractions.
  # design_formula: RHS-only formula, e.g. ~condition or ~mouse + condition.
  # design_label: character label written into the output `design` column.

  sample_ids <- rownames(sample_table)
  keep_cells <- cell_metadata$pseudobulk_sample_id %in% sample_ids
  cell_metadata <- cell_metadata[keep_cells, , drop = FALSE]
  counts <- counts[, keep_cells, drop = FALSE]

  cluster_id <- factor(rep("mg_selected", ncol(counts)))
  sample_id <- factor(cell_metadata$pseudobulk_sample_id, levels = sample_ids)
  group_id <- factor(
    as.character(cell_metadata$condition),
    levels = c(CONTROL_LEVEL, ESTIM_LEVEL)
  )

  col_data <- S4Vectors::DataFrame(
    cluster_id = cluster_id,
    sample_id = sample_id,
    group_id = group_id
  )
  rownames(col_data) <- colnames(counts)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = col_data
  )

  sce <- muscat::prepSCE(
    sce,
    kid = "cluster_id",
    sid = "sample_id",
    gid = "group_id",
    drop = FALSE
  )

  pb <- muscat::aggregateData(
    sce,
    assay = "counts",
    fun = "num.detected",
    by = c("cluster_id", "sample_id")
  )

  design <- stats::model.matrix(design_formula, data = sample_table)
  assert_full_rank(design, paste("muscat DD", design_label))
  coef_name <- tail(colnames(design), 1L)
  if (!identical(coef_name, "conditionestim")) {
    stop(
      "Expected muscat DD coefficient 'conditionestim' as the last design column; got: ",
      coef_name,
      call. = FALSE
    )
  }
  contrast <- limma::makeContrasts(contrasts = coef_name, levels = design)

  res <- muscat::pbDS(
    pb,
    method = "DD",
    design = design,
    contrast = contrast,
    min_cells = 0L,
    filter = "none",
    verbose = FALSE
  )

  tbl <- res$table[[coef_name]][["mg_selected"]]
  genes <- tbl$gene
  control_samples <- rownames(sample_table)[
    sample_table$condition == CONTROL_LEVEL
  ]
  estim_samples <- rownames(sample_table)[sample_table$condition == ESTIM_LEVEL]

  data.frame(
    gene = genes,
    logFC = tbl$logFC,
    logCPM = tbl$logCPM,
    F = tbl[["F"]],
    pvalue = tbl$p_val,
    padj = tbl$p_adj.loc,
    mean_detected_fraction_control = rowMeans(
      as.matrix(detection_fraction[genes, control_samples, drop = FALSE])
    ),
    mean_detected_fraction_estim = rowMeans(
      as.matrix(detection_fraction[genes, estim_samples, drop = FALSE])
    ),
    contrast = CONTRAST_DIRECTION,
    design = design_label,
    method = "muscat_edgeR_NB_optim",
    stringsAsFactors = FALSE
  )
}

map_genes_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (length(symbols) == 0) {
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  }
  orgdb <- get("org.Mm.eg.db", envir = asNamespace("org.Mm.eg.db"))
  mapped <- suppressMessages(clusterProfiler::bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = orgdb
  ))
  mapped <- mapped[
    !is.na(mapped$ENTREZID) & nzchar(mapped$ENTREZID),
    ,
    drop = FALSE
  ]
  unique(mapped)
}

write_go_ora <- function(significant_symbols, background_map, direction, path) {
  significant_map <- background_map[
    background_map$SYMBOL %in% significant_symbols,
    ,
    drop = FALSE
  ]
  significant_entrez <- unique(significant_map$ENTREZID)
  background_entrez <- unique(background_map$ENTREZID)
  if (length(significant_entrez) < MIN_ORA_GENES) {
    write_reason_tsv(
      path,
      paste0(
        "fewer than ",
        MIN_ORA_GENES,
        " mapped significant ",
        direction,
        " DE genes"
      ),
      n_input_genes = length(significant_symbols),
      n_mapped_genes = length(significant_entrez),
      min_required = MIN_ORA_GENES
    )
    return(invisible(NULL))
  }
  orgdb <- get("org.Mm.eg.db", envir = asNamespace("org.Mm.eg.db"))
  enrichment <- suppressMessages(clusterProfiler::enrichGO(
    gene = significant_entrez,
    universe = background_entrez,
    OrgDb = orgdb,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  ))
  enrichment_table <- as.data.frame(enrichment)
  if (nrow(enrichment_table) == 0) {
    write_reason_tsv(
      path,
      paste0("no GO BP ORA terms returned for ", direction, " DE genes"),
      n_input_genes = length(significant_symbols),
      n_mapped_genes = length(significant_entrez),
      min_required = MIN_ORA_GENES
    )
    return(invisible(NULL))
  }
  enrichment_table$direction <- direction
  write_tsv(enrichment_table, path)
}

write_go_gsea <- function(result_table, background_map, path, mapping_path) {
  ranked <- result_table[
    !is.na(result_table$stat),
    c("gene", "stat"),
    drop = FALSE
  ]
  ranked <- merge(
    ranked,
    background_map,
    by.x = "gene",
    by.y = "SYMBOL",
    all = FALSE,
    sort = FALSE
  )
  if (nrow(ranked) == 0) {
    write_reason_tsv(
      path,
      "no DESeq2 ranked statistics mapped to Entrez IDs",
      n_input_genes = nrow(result_table),
      n_mapped_genes = 0,
      min_required = MIN_GSEA_GENES
    )
    return(invisible(NULL))
  }
  ranked$abs_stat <- abs(ranked$stat)
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
  ranked <- ranked[ranked$selected_for_gsea, , drop = FALSE]
  ranked$abs_stat <- NULL
  ranked$selected_for_gsea <- NULL
  gene_list <- ranked$stat
  names(gene_list) <- ranked$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  if (length(gene_list) < MIN_GSEA_GENES) {
    write_reason_tsv(
      path,
      paste0("fewer than ", MIN_GSEA_GENES, " mapped ranked genes"),
      n_input_genes = nrow(result_table),
      n_mapped_genes = length(gene_list),
      min_required = MIN_GSEA_GENES
    )
    return(invisible(NULL))
  }
  orgdb <- get("org.Mm.eg.db", envir = asNamespace("org.Mm.eg.db"))
  set.seed(GSEA_SEED)
  gsea <- suppressMessages(clusterProfiler::gseGO(
    geneList = gene_list,
    OrgDb = orgdb,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  ))
  gsea_table <- as.data.frame(gsea)
  if (nrow(gsea_table) == 0) {
    write_reason_tsv(
      path,
      "no GO BP GSEA terms returned",
      n_input_genes = nrow(result_table),
      n_mapped_genes = length(gene_list),
      min_required = MIN_GSEA_GENES
    )
    return(invisible(NULL))
  }
  write_tsv(gsea_table, path)
}

# ---- validation ----

assert_scalar_character(input_path, "--input")
assert_scalar_character(cluster_column, "--cluster-column")
assert_scalar_character(condition_col, "--condition-col")
assert_scalar_character(control_label, "--control-label")
assert_scalar_character(estim_label, "--estim-label")
assert_scalar_character(counts_layer, "--counts-layer")
if (identical(control_label, estim_label)) {
  stop("--control-label and --estim-label must differ.", call. = FALSE)
}
if (!file.exists(input_path)) {
  stop("Input Seurat object does not exist: ", input_path, call. = FALSE)
}

# ---- load object and metadata ----

sobj <- readRDS(input_path)
meta <- sobj@meta.data
if (!is.data.frame(meta) || nrow(meta) == 0) {
  stop("Input object has missing or empty Seurat metadata.", call. = FALSE)
}
assert_no_missing_metadata(meta, c("Mouse", condition_col, cluster_column))

condition_labels <- trimws(as.character(meta[[condition_col]]))
observed_conditions <- sort(unique(condition_labels))
expected_conditions <- sort(c(control_label, estim_label))
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
  CONTROL_LEVEL,
  ifelse(condition_labels == estim_label, ESTIM_LEVEL, NA_character_)
)
if (any(is.na(condition))) {
  stop("Internal condition recoding failed.", call. = FALSE)
}
condition <- factor(condition, levels = c(CONTROL_LEVEL, ESTIM_LEVEL))

mouse <- trimws(as.character(meta$Mouse))
constructed_sample_id <- paste(
  paste0("Mouse_", safe_component(mouse)),
  as.character(condition),
  sep = "__"
)

if ("sample_id" %in% colnames(meta)) {
  assert_no_missing_metadata(meta, "sample_id")
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
    bad_samples <- sample_to_units$sample_id[
      sample_to_units$constructed_sample_id != 1L
    ]
    stop(
      "Metadata sample_id values are not consistent with Mouse × Condition units: ",
      paste(bad_samples, collapse = ", "),
      call. = FALSE
    )
  }
  unit_to_samples <- stats::aggregate(
    sample_id ~ constructed_sample_id,
    sample_map,
    function(x) length(unique(x))
  )
  if (any(unit_to_samples$sample_id != 1L)) {
    bad_units <- unit_to_samples$constructed_sample_id[
      unit_to_samples$sample_id != 1L
    ]
    stop(
      "Multiple sample_id values found within Mouse × Condition unit(s): ",
      paste(bad_units, collapse = ", "),
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
if (any(duplicated(sample_table$sample_id))) {
  stop(
    "Pseudobulk sample IDs are duplicated across distinct Mouse × Condition units.",
    call. = FALSE
  )
}
sample_table <- sample_table[
  order(sample_table$condition, sample_table$Mouse),
  ,
  drop = FALSE
]
rownames(sample_table) <- sample_table$sample_id
sample_table$condition <- factor(
  as.character(sample_table$condition),
  levels = c(CONTROL_LEVEL, ESTIM_LEVEL)
)
sample_table$mouse <- factor(sample_table$Mouse)

sample_counts_by_condition <- table(sample_table$condition)
if (any(sample_counts_by_condition[c(CONTROL_LEVEL, ESTIM_LEVEL)] < 1L)) {
  stop(
    "Need at least one Mouse × Condition pseudobulk sample per condition; observed ",
    paste(
      names(sample_counts_by_condition),
      sample_counts_by_condition,
      sep = "=",
      collapse = ", "
    ),
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
if (is.null(rownames(counts)) || any(rownames(counts) == "")) {
  stop("Counts matrix must have gene row names.", call. = FALSE)
}
if (is.null(colnames(counts)) || !identical(colnames(counts), rownames(meta))) {
  stop(
    "Counts matrix columns do not exactly match Seurat metadata row names.",
    call. = FALSE
  )
}
count_values <- counts@x
if (length(count_values) > 0L && any(count_values < 0)) {
  stop("Raw counts layer contains negative values.", call. = FALSE)
}
if (
  length(count_values) > 0L &&
    any(abs(count_values - round(count_values)) > sqrt(.Machine$double.eps))
) {
  stop("Raw counts layer contains non-integer values.", call. = FALSE)
}

# ---- pseudobulk aggregation ----

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

n_cells <- as.integer(tabulate(
  match(pseudobulk_sample_id, sample_ids),
  nbins = length(sample_ids)
))
sample_table$n_cells <- n_cells
sample_table$total_counts <- as.numeric(Matrix::colSums(pseudobulk_counts))
if (any(sample_table$n_cells < 1L)) {
  stop("At least one pseudobulk sample has zero cells.", call. = FALSE)
}

mouse_condition_table <- table(sample_table$Mouse, sample_table$condition)
paired_mice <- rownames(mouse_condition_table)[
  mouse_condition_table[, CONTROL_LEVEL] > 0L &
    mouse_condition_table[, ESTIM_LEVEL] > 0L
]
unmatched_mice <- setdiff(as.character(sample_table$Mouse), paired_mice)
sample_table$paired_mouse <- sample_table$Mouse %in% paired_mice
sample_table$sample_id_source <- sample_id_source
sample_table$analysis_unit <- "Mouse_x_Condition"

output_paths <- c(
  file.path(deg_dir, "pseudobulk_sample_summary.tsv"),
  file.path(deg_dir, "design_summary.tsv"),
  file.path(deg_dir, "deseq2_full_results.tsv"),
  file.path(deg_dir, "deseq2_significant_degs.tsv"),
  file.path(deg_dir, "deseq2_marker_overlap.tsv"),
  file.path(deg_dir, "detection_full_results.tsv"),
  file.path(deg_dir, "detection_marker_overlap.tsv"),
  file.path(deg_dir, "deseq2_paired_sensitivity_full_results.tsv"),
  file.path(deg_dir, "deseq2_paired_sensitivity_significant_degs.tsv"),
  file.path(deg_dir, "deseq2_paired_sensitivity_marker_overlap.tsv"),
  file.path(deg_dir, "detection_paired_sensitivity_full_results.tsv"),
  file.path(deg_dir, "detection_paired_sensitivity_marker_overlap.tsv"),
  file.path(deg_dir, "numbers.json"),
  file.path(enrichment_dir, "go_bp_ora_up.tsv"),
  file.path(enrichment_dir, "go_bp_ora_down.tsv"),
  file.path(enrichment_dir, "go_bp_gsea.tsv"),
  file.path(enrichment_dir, "go_bp_gsea_symbol_entrez_mapping.tsv"),
  de_dd_effect_scatter_png_path,
  de_dd_effect_scatter_pdf_path,
  de_dd_effect_scatter_notebook_png_path
)
assert_output_paths_clear(output_paths, overwrite_outputs)

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

message(
  "Running mg-selected DESeq2 contrast ",
  CONTRAST_DIRECTION,
  " using primary design ~ condition; paired Mouse design is sensitivity-only."
)

# ---- primary DESeq2 ----

primary_de <- run_deseq2(
  count_matrix = pseudobulk_counts,
  sample_table = sample_table,
  design_formula = ~condition,
  design_label = "primary_unpaired_condition",
  shrink_type = lfc_shrink_type
)
full_de <- primary_de$result_table
sig_de <- full_de[!is.na(full_de$padj) & full_de$padj < 0.05, , drop = FALSE]
write_tsv(full_de, file.path(deg_dir, "deseq2_full_results.tsv"))
write_tsv(sig_de, file.path(deg_dir, "deseq2_significant_degs.tsv"))
de_marker_overlap <- write_marker_overlap(
  full_de,
  primary_de$tested_genes,
  file.path(deg_dir, "deseq2_marker_overlap.tsv"),
  padj_column = "padj"
)

# ---- differential detection ----

detected_counts <- do.call(
  cbind,
  lapply(sample_ids, function(id) {
    Matrix::rowSums(counts[, pseudobulk_sample_id == id, drop = FALSE] > 0)
  })
)
rownames(detected_counts) <- rownames(counts)
colnames(detected_counts) <- sample_ids
detection_fraction <- sweep(detected_counts, 2, sample_table$n_cells, "/")

cell_metadata <- meta
cell_metadata$pseudobulk_sample_id <- pseudobulk_sample_id
cell_metadata$condition <- condition
cell_metadata$mouse <- factor(mouse)

full_detection <- run_detection_muscat_dd(
  counts = counts,
  cell_metadata = cell_metadata,
  sample_table = sample_table,
  detection_fraction = detection_fraction,
  design_formula = ~condition,
  design_label = "primary_unpaired_condition"
)
write_tsv(full_detection, file.path(deg_dir, "detection_full_results.tsv"))
detection_marker_overlap <- write_marker_overlap(
  full_detection,
  full_detection$gene,
  file.path(deg_dir, "detection_marker_overlap.tsv"),
  padj_column = "padj"
)

# ---- paired sensitivity analyses ----

paired_status <- "skipped"
paired_reason <- "fewer than two mice have both control and estim samples"
paired_de_n_degs <- NA_integer_
paired_detection_n_hits <- NA_integer_
paired_detection_n_tested <- NA_integer_

if (length(paired_mice) >= 2L) {
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
      count_matrix = paired_counts,
      sample_table = paired_sample_table,
      design_formula = ~ mouse + condition,
      design_label = "paired_mouse_condition_sensitivity",
      shrink_type = lfc_shrink_type
    )
    paired_full_de <- paired_de$result_table
    paired_sig_de <- paired_full_de[
      !is.na(paired_full_de$padj) & paired_full_de$padj < 0.05,
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
      file.path(deg_dir, "deseq2_paired_sensitivity_marker_overlap.tsv"),
      padj_column = "padj"
    )

    paired_cell_metadata <- cell_metadata[
      cell_metadata$mouse %in% paired_mice,
      ,
      drop = FALSE
    ]
    paired_detection_fraction <- detection_fraction[,
      rownames(paired_sample_table),
      drop = FALSE
    ]

    paired_full_detection <- run_detection_muscat_dd(
      counts = counts[, rownames(paired_cell_metadata), drop = FALSE],
      cell_metadata = paired_cell_metadata,
      sample_table = paired_sample_table,
      detection_fraction = paired_detection_fraction,
      design_formula = ~ mouse + condition,
      design_label = "paired_mouse_condition_sensitivity"
    )
    write_tsv(
      paired_full_detection,
      file.path(deg_dir, "detection_paired_sensitivity_full_results.tsv")
    )
    write_marker_overlap(
      paired_full_detection,
      paired_full_detection$gene,
      file.path(deg_dir, "detection_paired_sensitivity_marker_overlap.tsv"),
      padj_column = "padj"
    )
    paired_status <- "run"
    paired_reason <- "paired sensitivity used mice with both control and estim samples"
    paired_de_n_degs <- nrow(paired_sig_de)
    paired_detection_n_tested <- nrow(paired_full_detection)
    paired_detection_n_hits <- sum(
      !is.na(paired_full_detection$padj) & paired_full_detection$padj < 0.05
    )
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
    n_input_genes = nrow(pseudobulk_counts),
    n_mapped_genes = 0L,
    min_required = 0L
  )
  write_reason_tsv(
    file.path(deg_dir, "deseq2_paired_sensitivity_significant_degs.tsv"),
    paired_reason,
    n_input_genes = nrow(pseudobulk_counts),
    n_mapped_genes = 0L,
    min_required = 0L
  )
  write_reason_tsv(
    file.path(deg_dir, "deseq2_paired_sensitivity_marker_overlap.tsv"),
    paired_reason,
    n_input_genes = nrow(pseudobulk_counts),
    n_mapped_genes = 0L,
    min_required = 0L
  )
  write_reason_tsv(
    file.path(deg_dir, "detection_paired_sensitivity_full_results.tsv"),
    paired_reason,
    n_input_genes = nrow(counts),
    n_mapped_genes = 0L,
    min_required = 0L
  )
  write_reason_tsv(
    file.path(deg_dir, "detection_paired_sensitivity_marker_overlap.tsv"),
    paired_reason,
    n_input_genes = nrow(counts),
    n_mapped_genes = 0L,
    min_required = 0L
  )
}

# ---- DE vs differential detection effect scatter ----

de_dd_plot_data <- build_de_dd_effect_data(
  full_de,
  full_detection,
  "All samples by condition"
)
if (identical(paired_status, "run")) {
  paired_de_dd_plot_data <- build_de_dd_effect_data(
    paired_full_de,
    paired_full_detection,
    "Paired mice by condition"
  )
  de_dd_plot_data <- rbind(de_dd_plot_data, paired_de_dd_plot_data)
}

de_dd_effect_scatter <- plot_de_dd_effect_scatter(de_dd_plot_data)
de_dd_effect_scatter_width <- if (identical(paired_status, "run")) 9 else 5.5
dir.create(
  dirname(de_dd_effect_scatter_png_path),
  recursive = TRUE,
  showWarnings = FALSE
)
ggplot2::ggsave(
  filename = de_dd_effect_scatter_png_path,
  plot = de_dd_effect_scatter,
  width = de_dd_effect_scatter_width,
  height = 4.8,
  dpi = 300
)
ggplot2::ggsave(
  filename = de_dd_effect_scatter_pdf_path,
  plot = de_dd_effect_scatter,
  width = de_dd_effect_scatter_width,
  height = 4.8
)

dir.create(
  dirname(de_dd_effect_scatter_notebook_png_path),
  recursive = TRUE,
  showWarnings = FALSE
)
if (
  file.exists(de_dd_effect_scatter_notebook_png_path) ||
    nzchar(Sys.readlink(de_dd_effect_scatter_notebook_png_path))
) {
  unlink(de_dd_effect_scatter_notebook_png_path)
}
link_created <- file.symlink(
  de_dd_effect_scatter_png_path,
  de_dd_effect_scatter_notebook_png_path
)
if (!isTRUE(link_created)) {
  stop(
    "Failed to link notebook figure: ",
    de_dd_effect_scatter_notebook_png_path,
    call. = FALSE
  )
}
message(
  "Wrote MG-selected DE/DD effect scatter PNG: ",
  de_dd_effect_scatter_png_path
)
message(
  "Wrote MG-selected DE/DD effect scatter PDF: ",
  de_dd_effect_scatter_pdf_path
)
message("Linked notebook figure: ", de_dd_effect_scatter_notebook_png_path)

design_summary <- data.frame(
  analysis = c(
    "deseq2_primary",
    "detection_primary",
    "deseq2_paired_sensitivity",
    "detection_paired_sensitivity"
  ),
  design = c(
    "~ condition",
    "~ condition",
    "~ mouse + condition",
    "~ mouse + condition"
  ),
  method = c(
    "deseq2_wald",
    "muscat_edgeR_NB_optim",
    "deseq2_wald",
    "muscat_edgeR_NB_optim"
  ),
  status = c("run", "run", paired_status, paired_status),
  contrast = CONTRAST_DIRECTION,
  limitation = c(
    "uses all Mouse × Condition samples; mouse pairing is not modeled",
    "muscat edgeR_NB_optim: internal 90% detection filter and CDR-style offset applied per Mouse x Condition sample",
    paired_reason,
    paired_reason
  ),
  included_mice = c(
    paste(unique(sample_table$Mouse), collapse = ","),
    paste(unique(sample_table$Mouse), collapse = ","),
    paste(paired_mice, collapse = ","),
    paste(paired_mice, collapse = ",")
  ),
  unmatched_mice = paste(unique(unmatched_mice), collapse = ","),
  n_samples = c(
    nrow(sample_table),
    nrow(sample_table),
    sum(sample_table$Mouse %in% paired_mice),
    sum(sample_table$Mouse %in% paired_mice)
  ),
  source_input = input_path,
  cluster_column = cluster_column,
  counts_layer = counts_layer,
  lfc_shrink_type = lfc_shrink_type,
  gsea_seed = GSEA_SEED,
  stringsAsFactors = FALSE
)
write_tsv(design_summary, file.path(deg_dir, "design_summary.tsv"))

# ---- enrichment follow-ups ----

background_map <- map_genes_to_entrez(primary_de$tested_genes)
if (nrow(background_map) == 0) {
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
write_go_ora(
  up_genes,
  background_map,
  direction = "up_estim_vs_control",
  path = file.path(enrichment_dir, "go_bp_ora_up.tsv")
)
write_go_ora(
  down_genes,
  background_map,
  direction = "down_estim_vs_control",
  path = file.path(enrichment_dir, "go_bp_ora_down.tsv")
)
write_go_gsea(
  full_de,
  background_map,
  path = file.path(enrichment_dir, "go_bp_gsea.tsv"),
  mapping_path = file.path(
    enrichment_dir,
    "go_bp_gsea_symbol_entrez_mapping.tsv"
  )
)

# ---- reportable values ----

reportable_values <- list(
  n_samples = nrow(sample_table),
  n_cells = ncol(counts),
  n_tested_genes = length(primary_de$tested_genes),
  n_detection_tested_genes = nrow(full_detection),
  n_degs = nrow(sig_de),
  n_marker_degs = sum(de_marker_overlap$significant, na.rm = TRUE),
  n_detection_marker_hits = sum(
    detection_marker_overlap$significant,
    na.rm = TRUE
  ),
  n_unmatched_mice = length(unique(unmatched_mice)),
  source_input = input_path,
  cluster_column = cluster_column,
  counts_layer = counts_layer,
  gsea_seed = GSEA_SEED,
  lfc_shrink_type = primary_de$shrink_type,
  dd_method = "muscat_edgeR_NB_optim",
  primary_design = "~ condition",
  paired_sensitivity_design = "~ mouse + condition",
  paired_sensitivity_included_mice = paste(paired_mice, collapse = ","),
  paired_sensitivity_status = paired_status,
  paired_sensitivity_n_degs = paired_de_n_degs,
  paired_sensitivity_detection_tested_genes = paired_detection_n_tested,
  paired_sensitivity_detection_hits = paired_detection_n_hits
)
register_or_write_numbers(reportable_values, deg_dir)

message("Saved mg-selected DE outputs under ", deg_dir)
message("Saved mg-selected enrichment outputs under ", enrichment_dir)
