# Create data-only supplementary tables from existing analysis exports.
# source("scripts/prepare-supp-tables.R"); ESPI_OVERWRITE=true permits replacement.

library(here)
library(dplyr)
library(purrr)
library(readr)
library(readxl)
library(openxlsx)
devtools::load_all(here::here(), quiet = TRUE)

# ---- parameters ----
config <- publication_config()
base_dir <- config$paths$data
share_dir <- file.path(base_dir, "Carl Megan Share Sep5")
output_dir <- Sys.getenv(
  "ESPI_SUPP_TABLE_DIR",
  file.path(share_dir, "supp_tables")
)
overwrite <- identical(tolower(Sys.getenv("ESPI_OVERWRITE", "false")), "true")
suffix <- sub("^cluster_", "", config$selected$mg$column)
cluster_labels <- c(
  `1` = "Müller glia",
  `2` = "Neurogenic progenitor",
  `3` = "Cone bipolar-like",
  `4` = "Other",
  `5` = "Proliferative Müller glia"
)
condition_codes <- stats::setNames(
  c("control", "estim"),
  c(config$conditions$control, config$conditions$estim)
)
condition_labels <- c(
  control = config$conditions$control_display,
  estim = config$conditions$estim_display
)
marker_references <- character()

# ---- inputs ----
files <- c(
  markers = file.path(
    base_dir,
    "tables/mg_selected",
    paste0("find_all_markers_wilcox_data_", suffix, ".csv")
  ),
  top5 = file.path(
    base_dir,
    "tables/mg_selected",
    paste0("find_all_markers_wilcox_top5_data_", suffix, ".csv")
  ),
  marker_summary = file.path(
    base_dir,
    "tables/mg_selected",
    paste0("find_all_markers_summary_data_", suffix, ".csv")
  ),
  selection = file.path(
    base_dir,
    "tables/mg_selected/mg_selected_cluster_selection.tsv"
  ),
  qc = file.path(base_dir, "tables/qc/post_filtering_qc_stats.tsv")
)
de_files <- c(
  primary = "deseq2_full_results.tsv",
  primary_sig = "deseq2_significant_degs.tsv",
  paired = "deseq2_paired_sensitivity_full_results.tsv",
  paired_sig = "deseq2_paired_sensitivity_significant_degs.tsv",
  samples = "pseudobulk_sample_table.tsv",
  design = "design_summary.tsv"
)
go_files <- c(
  ora_up = "go_bp_ora_up.tsv",
  ora_down = "go_bp_ora_down.tsv",
  gsea = "go_bp_gsea.tsv",
  ora_up_reduced = "go_bp_ora_up_simplified.tsv",
  ora_down_reduced = "go_bp_ora_down_simplified.tsv",
  gsea_reduced = "go_bp_gsea_simplified.tsv",
  mapping = "go_bp_gsea_symbol_entrez_mapping.tsv"
)
files <- c(
  files,
  stats::setNames(
    file.path(base_dir, "degs/mg_selected", de_files),
    names(de_files)
  ),
  stats::setNames(
    file.path(base_dir, "enrichment/mg_selected", go_files),
    names(go_files)
  )
)
stopifnot(all(file.exists(files)))
data <- purrr::map(files, function(path) {
  readr::read_delim(
    path,
    delim = if (grepl("\\.csv$", path)) "," else "\t",
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA", "NaN"),
    show_col_types = FALSE,
    progress = FALSE
  )
})

text_fields <- c(
  "gene",
  "ENTREZID",
  "cluster",
  "marker_identity",
  "Mouse",
  "mouse",
  "Sample",
  "ID",
  "GeneRatio",
  "BgRatio",
  "geneID",
  "core_enrichment",
  "References"
)
excel_column <- function(x, key) {
  if (!is.character(x)) {
    return(x)
  }
  x[x %in% c("", "NA", "NaN")] <- NA_character_
  if (key %in% text_fields || all(is.na(x))) {
    return(x)
  }
  values <- x[!is.na(x)]
  if (all(values %in% c("TRUE", "FALSE"))) {
    return(ifelse(is.na(x), NA_character_, ifelse(x == "TRUE", "Yes", "No")))
  }
  if (all(grepl("^-?[0-9]+$", values))) {
    numbers <- as.numeric(x)
    return(
      if (all(abs(numbers[!is.na(numbers)]) <= .Machine$integer.max)) {
        as.integer(numbers)
      } else {
        numbers
      }
    )
  }
  if (
    all(grepl("^[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?$", values))
  ) {
    numbers <- as.numeric(x)
    # Preserve tiny P values affected by scientific-notation parsing in this R build.
    scientific <- which(!is.na(x) & grepl("[eE]", x))
    numbers[scientific] <- vapply(
      strsplit(x[scientific], "[eE]"),
      function(parts) {
        exponent <- as.integer(parts[2])
        as.numeric(parts[1]) *
          10^max(exponent, -300) *
          10^min(exponent + 300, 0)
      },
      numeric(1)
    )
    if (all(is.finite(numbers[!is.na(numbers)]))) return(numbers)
  }
  x
}

# ---- checks ----
stopifnot(
  setequal(data$markers$cluster, names(cluster_labels)),
  all(as.numeric(data$markers$p_val_adj) <= 0.01),
  all(as.numeric(data$markers$avg_log2FC) > 0),
  !anyDuplicated(data$markers[c("cluster", "gene")]),
  identical(names(data$top5), names(data$markers)),
  all(purrr::map2_lgl(
    data$top5,
    dplyr::filter(data$markers, as.integer(rank_within_cluster) <= 5L),
    identical
  )),
  all(data$design$cluster_column == config$selected$mg$column)
)
stopifnot(all(
  as.integer(table(data$markers$cluster)[
    data$marker_summary$marker_identity
  ]) ==
    as.integer(data$marker_summary$n_significant_markers)
))
for (analysis in c("primary", "paired")) {
  results <- data[[analysis]]
  stopifnot(setequal(
    results$gene[!is.na(results$padj) & as.numeric(results$padj) < 0.05],
    data[[paste0(analysis, "_sig")]]$gene
  ))
}
samples <- data$samples
stopifnot(
  !anyNA(samples[c("Mouse", "condition")]),
  !anyDuplicated(samples[c("Mouse", "condition")])
)
stopifnot(
  sum(as.integer(data$selection$n_cells[data$selection$exclude == "FALSE"])) ==
    sum(as.integer(samples$n_cells)),
  sum(as.integer(data$marker_summary$n_cells)) ==
    sum(as.integer(samples$n_cells))
)

# ---- sample QC ----
qc_original <- readxl::read_excel(
  file.path(share_dir, "supp_table_sample_qc.xlsx"),
  col_types = "text",
  .name_repair = "minimal"
) |>
  dplyr::filter(grepl("^S[0-9]+$", Sample))
stopifnot(!anyDuplicated(qc_original$Sample), !anyDuplicated(data$qc$Sample))
qc_table <- purrr::map(seq_len(nrow(qc_original)), function(i) {
  row <- qc_original[i, ]
  qc <- data$qc |> dplyr::filter(Sample == row$Sample)
  stopifnot(nrow(qc) == 1L, row$Mouse == qc$Mouse)
  code <- unname(condition_codes[qc$Condition])
  selected <- samples |> dplyr::filter(Mouse == row$Mouse, condition == code)
  stopifnot(nrow(selected) == 1L)
  passed <- as.integer(row[["Final retained"]])
  stopifnot(
    passed == as.integer(qc$n_cells),
    as.integer(row$Singlets) - as.integer(row[["Failed QC*"]]) == passed,
    as.integer(row[["Called cells"]]) >= as.integer(row$Singlets),
    as.integer(selected$n_cells) <= passed
  )
  metrics <- as.numeric(row[c(
    "Mean UMI/cell",
    "Mean Features/cell",
    "Mean mitochondrial RNA (%)"
  )])
  stopifnot(isTRUE(all.equal(
    metrics,
    as.numeric(qc[c(
      "mean_umi_per_cell",
      "mean_features_per_cell",
      "mean_percent_mt"
    )]),
    tolerance = 1e-12
  )))
  tibble::tibble(
    Sample = row$Sample,
    Mouse = row$Mouse,
    Condition = unname(condition_labels[code]),
    `Called cells` = as.integer(row[["Called cells"]]),
    Singlets = as.integer(row$Singlets),
    `Failed QC` = as.integer(row[["Failed QC*"]]),
    `QC-passing cells` = passed,
    `MG-selected cells` = as.integer(selected$n_cells),
    `Mean UMI/cell` = metrics[1],
    `Mean genes/cell` = metrics[2],
    `Mean mitochondrial RNA (%)` = metrics[3]
  )
}) |>
  purrr::list_rbind()
stopifnot(
  nrow(qc_table) == nrow(samples),
  nrow(qc_table) == nrow(data$qc),
  sum(qc_table[["QC-passing cells"]]) ==
    sum(as.integer(data$selection$n_cells)),
  sum(qc_table[["MG-selected cells"]]) == sum(as.integer(samples$n_cells))
)
qc_total <- qc_table[NA_integer_, ]
qc_total$Sample <- "All"
qc_total[4:8] <- purrr::map(qc_table[4:8], sum)
qc_total[9:11] <- purrr::map(qc_table[9:11], mean)
qc_table <- dplyr::bind_rows(qc_table, qc_total)

# ---- cell-type markers ----
read_marker_list <- function(path, variable) {
  expressions <- purrr::keep(as.list(parse(path)), function(x) {
    is.call(x) &&
      identical(x[[1]], as.name("<-")) &&
      identical(x[[2]], as.name(variable)) &&
      is.call(x[[3]]) &&
      identical(x[[3]][[1]], as.name("list"))
  })
  stopifnot(length(expressions) == 1L)
  eval(expressions[[1]][[3]], envir = new.env(parent = baseenv()))
}
marker_sets <- read_marker_list(
  file.path(config$paths$project, "data-raw/cell-type-marker-genes.R"),
  "cell_type_marker_genes"
)
main_sets <- read_marker_list(
  file.path(config$paths$project, "scripts/02-publication-figures.R"),
  "main_violin_markers"
)
main_types <- c(
  MG = "muller_glia",
  `Activated MG` = "activated_muller_glia",
  Proliferative = "proliferative",
  `Neurogenic progenitor` = "neurogenic_progenitor",
  `Cone bipolar` = "cone_bipolar"
)
# This is a marker catalogue, not a specification of individual scoring modules.
for (name in names(main_sets)) {
  type <- unname(main_types[name])
  stopifnot(!is.na(type), type %in% names(marker_sets))
  marker_sets[[type]] <- union(marker_sets[[type]], main_sets[[name]])
}
marker_table <- purrr::imap(marker_sets, function(genes, type) {
  label <- unname(cell_type_marker_labels[type])
  stopifnot(!is.na(label), !anyNA(genes), length(genes) > 0L)
  tibble::tibble(
    `Cell type` = label,
    Genes = paste(unique(genes), collapse = ", "),
    References = unname(marker_references[label])
  )
}) |>
  purrr::list_rbind()

# ---- result tables ----
annotate <- function(x) {
  key <- if ("cluster" %in% names(x)) "cluster" else "marker_identity"
  stopifnot(all(x[[key]] %in% names(cluster_labels)))
  x |>
    dplyr::mutate(cluster_label = unname(cluster_labels[x[[key]]])) |>
    dplyr::relocate(dplyr::all_of(key), cluster_label)
}
sample_key <- data$qc |>
  dplyr::transmute(
    Sample,
    Mouse,
    condition = unname(condition_codes[Condition])
  )
stopifnot(!anyDuplicated(sample_key[c("Mouse", "condition")]))
sample_table <- samples |>
  dplyr::left_join(
    sample_key,
    by = c("Mouse", "condition"),
    relationship = "one-to-one"
  ) |>
  dplyr::transmute(
    Sample,
    Mouse,
    Condition = unname(condition_labels[condition]),
    n_cells,
    total_counts,
    paired_mouse
  )
stopifnot(!anyNA(sample_table$Sample))
mapping <- data$mapping |> dplyr::filter(selected_for_gsea == "TRUE")
stopifnot(
  !anyDuplicated(mapping$ENTREZID),
  !anyNA(mapping[c("ENTREZID", "gene")])
)
symbols <- stats::setNames(mapping$gene, mapping$ENTREZID)
add_symbols <- function(x) {
  x |>
    dplyr::mutate(
      core_gene_symbols = purrr::map_chr(core_enrichment, function(value) {
        if (is.na(value)) {
          return(NA_character_)
        }
        ids <- strsplit(value, "/", fixed = TRUE)[[1]]
        stopifnot(all(ids %in% names(symbols)))
        paste(unname(symbols[ids]), collapse = "; ")
      })
    )
}
workbooks <- list(
  supp_table_cluster_annotation_markers.xlsx = list(
    Cluster_annotations = annotate(data$marker_summary),
    All_positive_markers = annotate(data$markers),
    Top5_markers = annotate(data$top5)
  ),
  supp_table_pseudobulk_de.xlsx = list(
    Primary_significant = data$primary_sig,
    Primary_all_genes = data$primary,
    Paired_significant = data$paired_sig,
    Paired_all_genes = data$paired,
    Samples = sample_table
  ),
  supp_table_go_enrichment.xlsx = list(
    GSEA_reduced = add_symbols(data$gsea_reduced),
    ORA_down_reduced = data$ora_down_reduced,
    ORA_up_reduced = data$ora_up_reduced,
    GSEA_all = add_symbols(data$gsea),
    ORA_down_all = data$ora_down,
    ORA_up_all = data$ora_up
  ),
  supp_table_sample_qc.xlsx = list(Sample_QC = qc_table),
  supp_table_curated_marker_module_sets.xlsx = list(
    Cell_type_markers = marker_table
  )
)
remove_columns <- c(
  "decision_source",
  "n_top_markers",
  "contrast",
  "design",
  "lfc_shrink_type",
  "core_enrichment",
  "leading_edge",
  "log2err"
)
headers <- c(
  gene = "Gene",
  cluster = "Cluster",
  marker_identity = "Cluster",
  cluster_label = "Cell type",
  rank_within_cluster = "Rank",
  avg_log2FC = "Mean log2 FC",
  p_val = "Wilcoxon P",
  p_val_adj = "P (Bonferroni)",
  pct.1 = "Detected: cluster (%)",
  pct.2 = "Detected: other cells (%)",
  pct_diff = "Detection difference (%)",
  n_cells = "Cells",
  n_significant_markers = "Positive markers",
  baseMean = "Mean normalized count",
  log2FoldChange = "Log2 FC (shrunken)",
  lfcSE = "LFC uncertainty",
  stat = "Wald statistic",
  pvalue = "P value",
  padj = "Adjusted P",
  unshrunkenLog2FoldChange = "Log2 FC (unshrunken)",
  mean_count_control = "Mean count: p27CKO",
  mean_count_estim = "Mean count: E-Stim",
  total_counts = "Total counts",
  paired_mouse = "Paired mouse",
  mouse = "Mouse",
  condition = "Condition",
  ID = "GO ID",
  Description = "GO biological process",
  GeneRatio = "Gene ratio",
  BgRatio = "Background ratio",
  RichFactor = "Rich factor",
  FoldEnrichment = "Fold enrichment",
  zScore = "Z score",
  p.adjust = "BH-adjusted P",
  qvalue = "Q value",
  geneID = "Genes",
  Count = "Gene count",
  setSize = "Gene-set size",
  enrichmentScore = "Enrichment score",
  NES = "NES",
  rank = "Rank",
  core_gene_symbols = "Leading-edge genes"
)

# ---- write and check ----
write_book <- function(sheets, filename) {
  wb <- openxlsx::createWorkbook()
  prepared <- purrr::imap(sheets, function(raw, sheet) {
    raw <- dplyr::select(raw, -dplyr::any_of(remove_columns))
    x <- purrr::imap(raw, excel_column) |> tibble::as_tibble()
    fallback <- gsub("_", " ", names(x), fixed = TRUE)
    names(x) <- unname(dplyr::coalesce(headers[names(raw)], fallback))
    stopifnot(
      !anyDuplicated(names(x)),
      all(vapply(
        x,
        function(z) all(nchar(as.character(z)) <= 32767, na.rm = TRUE),
        logical(1)
      ))
    )
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, x, headerStyle = NULL, keepNA = FALSE)
    x
  })
  path <- output_path(output_dir, filename, overwrite = overwrite)
  openxlsx::saveWorkbook(wb, path, overwrite = overwrite)
  stopifnot(identical(readxl::excel_sheets(path), names(prepared)))
  purrr::iwalk(prepared, function(expected, sheet) {
    types <- vapply(
      expected,
      function(x) if (is.numeric(x)) "numeric" else "text",
      character(1)
    )
    actual <- readxl::read_excel(
      path,
      sheet = sheet,
      col_types = unname(types),
      range = readxl::cell_limits(
        c(1L, 1L),
        c(nrow(expected) + 1L, ncol(expected))
      ),
      .name_repair = "minimal",
      trim_ws = FALSE
    )
    stopifnot(
      identical(names(actual), names(expected)),
      nrow(actual) == nrow(expected)
    )
    for (key in names(expected)) {
      wanted <- expected[[key]]
      observed <- actual[[key]]
      stopifnot(identical(is.na(wanted), is.na(observed)))
      keep <- !is.na(wanted)
      if (is.numeric(wanted)) {
        stopifnot(all(
          abs(observed[keep] - wanted[keep]) <= 1e-14 * abs(wanted[keep])
        ))
      } else {
        stopifnot(identical(observed[keep], wanted[keep]))
      }
    }
  })
  invisible(path)
}
if (!overwrite && any(file.exists(file.path(output_dir, names(workbooks))))) {
  stop(
    "Output exists; set ESPI_OVERWRITE=true only when replacing publication tables."
  )
}
purrr::iwalk(workbooks, write_book)
