#!/usr/bin/env Rscript

# Select and preprocess the Müller glia-enriched cell subset from a clustered
# source object. The output objects are ready for scripts/cluster-sobj.R.

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
  library(scclrR)
})
here::i_am("scripts/filter-mg-subset.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})
utils::data("cell_type_marker_genes", package = "ESPI", envir = environment())
utils::data("cell_type_marker_labels", package = "ESPI", envir = environment())

# ---- parameters ----

args <- commandArgs(trailingOnly = TRUE)
allowed_flags <- c(
  "--input",
  "--cluster-column",
  "--dims",
  "--score-slot",
  "--score-layer",
  "--dataset-tag"
)
unknown_flags <- args[startsWith(args, "--") & !args %in% allowed_flags]
if (length(unknown_flags) > 0) {
  stop(
    "Unknown argument(s): ",
    paste(unknown_flags, collapse = ", "),
    call. = FALSE
  )
}
arg <- function(name) {
  i <- match(name, args)
  if (is.na(i)) {
    return(NULL)
  }
  if (i == length(args) || startsWith(args[[i + 1]], "--")) {
    return(TRUE)
  }
  args[[i + 1]]
}
arg_value <- function(name, default = NULL, required = FALSE) {
  value <- arg(name)
  if (identical(value, TRUE)) {
    stop("Missing value for ", name, call. = FALSE)
  }
  if (is.null(value)) {
    if (required) {
      stop("Missing required argument ", name, call. = FALSE)
    }
    return(default)
  }
  value
}
parse_positive_int <- function(value, name) {
  parsed <- suppressWarnings(as.integer(value))
  if (
    length(parsed) != 1 ||
      is.na(parsed) ||
      !is.finite(parsed) ||
      parsed <= 0 ||
      !identical(as.character(parsed), as.character(value))
  ) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  parsed
}
validate_tag <- function(tag, name) {
  if (
    !is.character(tag) ||
      length(tag) != 1 ||
      is.na(tag) ||
      !nzchar(tag) ||
      !grepl("^[A-Za-z0-9_]+$", tag)
  ) {
    stop(name, " must match ^[A-Za-z0-9_]+$.", call. = FALSE)
  }
  tag
}
sort_cluster_levels <- function(x) {
  x <- unique(as.character(x))
  if (all(grepl("^-?[0-9]+$", x))) {
    return(as.character(sort(as.integer(x))))
  }
  sort(x, method = "radix")
}
get_assay_layer <- function(sobj, assay, layer) {
  available_layers <- SeuratObject::Layers(sobj[[assay]])
  if (!layer %in% available_layers) {
    stop(
      "Missing layer '",
      layer,
      "' in assay ",
      assay,
      ". Available layers: ",
      paste(available_layers, collapse = ", "),
      call. = FALSE
    )
  }
  SeuratObject::LayerData(sobj[[assay]], layer = layer)
}
wilcox_greater <- function(in_values, out_values) {
  values <- c(in_values, out_values)
  if (length(unique(values)) <= 1) {
    return(1)
  }
  stats::wilcox.test(
    in_values,
    out_values,
    alternative = "greater",
    exact = FALSE
  )$p.value
}
strip_previous_cluster_state <- function(
  sobj,
  source_cluster_column,
  extra_drop
) {
  sobj@reductions <- list()
  sobj@graphs <- list()
  sobj@neighbors <- list()
  sobj@misc$clustering <- NULL
  drop_cols <- union(
    grep("^cluster_|^seurat_clusters$", colnames(sobj@meta.data), value = TRUE),
    extra_drop
  )
  drop_cols <- setdiff(drop_cols, source_cluster_column)
  if (length(drop_cols) > 0) {
    sobj@meta.data[drop_cols] <- NULL
  }
  sobj
}

input_path <- arg_value(
  "--input",
  default = file.path(
    CURRENT_OBJECT_DIR,
    "cluster_pflog_no_filter_cc_elbow20.rds"
  )
)
cluster_column <- arg_value(
  "--cluster-column",
  default = "cluster_pflog_no_filter_cc_dims50_res0.3"
)
dims <- parse_positive_int(arg_value("--dims", default = "50"), "--dims")
score_slot <- arg_value(
  "--score-slot",
  default = arg_value("--score-layer", default = "data")
)
dataset_tag <- validate_tag(
  arg_value("--dataset-tag", default = "mg_selected"),
  "--dataset-tag"
)
MARKER_EXCLUDE_CLASSES <- c("microglia", "photoreceptor")
MARKER_MIN_TOP_SCORE <- 0.5
MARKER_MIN_SCORE_MARGIN <- 0.25
CDKN1B_EXPRESSION_MAX_Q <- 0.05
CDKN1B_DETECTION_MAX_Q <- 0.05
CDKN1B_MIN_DETECTION_FRACTION <- 0.20

# ---- validation ----

if (!file.exists(input_path)) {
  stop("Input Seurat object does not exist: ", input_path, call. = FALSE)
}
if (!score_slot %in% c("data", "counts", "scale.data")) {
  stop(
    "--score-slot must be one of data, counts, or scale.data; got ",
    score_slot,
    call. = FALSE
  )
}
if (!identical(names(cell_type_marker_genes), names(cell_type_marker_labels))) {
  stop(
    "cell_type_marker_genes and cell_type_marker_labels must have identical names.",
    call. = FALSE
  )
}
if (!all(c("microglia", "photoreceptor") %in% names(cell_type_marker_genes))) {
  stop(
    "cell_type_marker_genes must include microglia and photoreceptor marker sets.",
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

# ---- work ----

sobj <- readRDS(input_path)
if (!inherits(sobj, "Seurat")) {
  stop("Input is not a Seurat object: ", input_path, call. = FALSE)
}
if (!"RNA" %in% SeuratObject::Assays(sobj)) {
  stop("Input Seurat object does not contain an RNA assay.", call. = FALSE)
}
DefaultAssay(sobj) <- "RNA"
if (!inherits(sobj[["RNA"]], "Assay5")) {
  sobj[["RNA"]] <- as(sobj[["RNA"]], Class = "Assay5")
}
assay <- "RNA"
available_layers <- SeuratObject::Layers(sobj[[assay]])
if (!score_slot %in% available_layers) {
  stop(
    "Missing score slot/layer '",
    score_slot,
    "' in assay ",
    assay,
    ". Available layers: ",
    paste(available_layers, collapse = ", "),
    call. = FALSE
  )
}
if (!"counts" %in% available_layers) {
  stop("Missing counts layer in RNA assay.", call. = FALSE)
}
if (!cluster_column %in% colnames(sobj@meta.data)) {
  stop("Missing cluster metadata column: ", cluster_column, call. = FALSE)
}
cluster_values <- as.character(sobj@meta.data[[cluster_column]])
if (any(is.na(cluster_values)) || any(!nzchar(cluster_values))) {
  stop(
    "Cluster metadata column contains missing or empty values.",
    call. = FALSE
  )
}
cluster_levels <- sort_cluster_levels(cluster_values)
if (length(cluster_levels) < 2) {
  stop(
    "Cluster metadata column must contain at least two clusters.",
    call. = FALSE
  )
}
sobj@meta.data[[cluster_column]] <- factor(
  cluster_values,
  levels = cluster_levels
)

missing_markers <- setdiff(marker_table$gene, rownames(sobj))
if (length(missing_markers) > 0) {
  stop(
    "Marker gene(s) missing from the Seurat object: ",
    paste(missing_markers, collapse = ", "),
    call. = FALSE
  )
}
if (!"Cdkn1b" %in% rownames(sobj)) {
  stop("Cdkn1b is missing from the Seurat object.", call. = FALSE)
}

module_prefix <- "mg_selection_marker_score"
module_score_cols <- paste0(module_prefix, seq_along(cell_type_marker_genes))
existing_module_cols <- intersect(module_score_cols, colnames(sobj@meta.data))
if (length(existing_module_cols) > 0) {
  sobj@meta.data[existing_module_cols] <- NULL
}
sobj <- Seurat::AddModuleScore(
  object = sobj,
  features = cell_type_marker_genes,
  assay = assay,
  name = module_prefix,
  seed = SEED,
  search = FALSE,
  slot = score_slot
)
missing_module_cols <- setdiff(module_score_cols, colnames(sobj@meta.data))
if (length(missing_module_cols) > 0) {
  stop(
    "Seurat::AddModuleScore did not create expected column(s): ",
    paste(missing_module_cols, collapse = ", "),
    call. = FALSE
  )
}
marker_score_cols <- paste0("marker_score_", names(cell_type_marker_genes))
module_scores <- sobj@meta.data[module_score_cols]
colnames(module_scores) <- marker_score_cols
module_scores$cluster <- as.character(sobj@meta.data[[cluster_column]])
cluster_marker_scores <- stats::aggregate(
  module_scores[marker_score_cols],
  by = list(cluster = module_scores$cluster),
  FUN = mean
)
cluster_marker_scores <- cluster_marker_scores[
  match(cluster_levels, cluster_marker_scores$cluster),
]

counts <- get_assay_layer(sobj, assay, "counts")
cdkn1b_counts <- as.numeric(counts["Cdkn1b", colnames(sobj), drop = TRUE])
cdkn1b_detected <- as.integer(cdkn1b_counts > 0)
cdkn1b_expression_layer <- if ("pflog" %in% available_layers) {
  "pflog"
} else if ("data" %in% available_layers) {
  "data"
} else {
  "counts"
}
expression_matrix <- get_assay_layer(sobj, assay, cdkn1b_expression_layer)
cdkn1b_expression <- as.numeric(
  expression_matrix["Cdkn1b", colnames(sobj), drop = TRUE]
)
sobj$Cdkn1b_selection_expression <- cdkn1b_expression

marker_decisions <- lapply(seq_len(nrow(cluster_marker_scores)), function(i) {
  scores <- as.numeric(cluster_marker_scores[i, marker_score_cols])
  ord <- order(scores, decreasing = TRUE)
  top_index <- ord[[1]]
  second_index <- ord[[2]]
  top_class <- names(cell_type_marker_genes)[[top_index]]
  top_score <- scores[[top_index]]
  second_class <- names(cell_type_marker_genes)[[second_index]]
  second_score <- scores[[second_index]]
  margin <- top_score - second_score
  data.frame(
    cluster = cluster_marker_scores$cluster[[i]],
    top_marker_class = top_class,
    top_marker_label = unname(cell_type_marker_labels[[top_class]]),
    top_marker_score = top_score,
    second_marker_class = second_class,
    second_marker_label = unname(cell_type_marker_labels[[second_class]]),
    second_marker_score = second_score,
    marker_score_margin = margin,
    marker_exclude = top_class %in%
      MARKER_EXCLUDE_CLASSES &&
      top_score >= MARKER_MIN_TOP_SCORE &&
      margin >= MARKER_MIN_SCORE_MARGIN,
    stringsAsFactors = FALSE
  )
})
marker_decisions <- do.call(rbind, marker_decisions)

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
cdkn1b_stats$cdkn1b_exclude <-
  cdkn1b_stats$cdkn1b_expression_q < CDKN1B_EXPRESSION_MAX_Q &
  cdkn1b_stats$cdkn1b_detection_q < CDKN1B_DETECTION_MAX_Q &
  cdkn1b_stats$cdkn1b_detection_fraction >= CDKN1B_MIN_DETECTION_FRACTION

decision_table <- merge(
  cluster_marker_scores,
  marker_decisions,
  by = "cluster",
  sort = FALSE
)
decision_table <- merge(
  decision_table,
  cdkn1b_stats,
  by = "cluster",
  sort = FALSE
)
decision_table <- decision_table[
  match(cluster_levels, decision_table$cluster),
]
decision_table$exclude <- decision_table$marker_exclude |
  decision_table$cdkn1b_exclude
decision_table$exclusion_reasons <- vapply(
  seq_len(nrow(decision_table)),
  function(i) {
    reasons <- character()
    if (isTRUE(decision_table$marker_exclude[[i]])) {
      reasons <- c(
        reasons,
        paste0("marker_", decision_table$top_marker_class[[i]])
      )
    }
    if (isTRUE(decision_table$cdkn1b_exclude[[i]])) {
      reasons <- c(reasons, "Cdkn1b_high")
    }
    paste(reasons, collapse = ";")
  },
  character(1)
)
decision_table$marker_exclude_classes <- paste(
  MARKER_EXCLUDE_CLASSES,
  collapse = ";"
)
decision_table$marker_min_top_score <- MARKER_MIN_TOP_SCORE
decision_table$marker_min_score_margin <- MARKER_MIN_SCORE_MARGIN
decision_table$cdkn1b_expression_max_q <- CDKN1B_EXPRESSION_MAX_Q
decision_table$cdkn1b_detection_max_q <- CDKN1B_DETECTION_MAX_Q
decision_table$cdkn1b_min_detection_fraction <- CDKN1B_MIN_DETECTION_FRACTION
excluded_clusters <- decision_table$cluster[decision_table$exclude]
if (length(excluded_clusters) == 0) {
  stop("No clusters meet exclusion criteria.", call. = FALSE)
}

selection_dir <- file.path(TABLE_DIR, dataset_tag)
dir.create(selection_dir, recursive = TRUE, showWarnings = FALSE)
selection_path <- file.path(selection_dir, "mg_selected_cluster_selection.tsv")
utils::write.table(
  decision_table,
  file = selection_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_dir <- file.path(FIGURE_DIR, dataset_tag)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
selection_plot_features <- c(
  module_score_cols[[match("microglia", names(cell_type_marker_genes))]],
  module_score_cols[[match("photoreceptor", names(cell_type_marker_genes))]],
  "Cdkn1b_selection_expression"
)
selection_plot <- Seurat::VlnPlot(
  sobj,
  features = selection_plot_features,
  group.by = cluster_column,
  assay = assay,
  layer = score_slot,
  pt.size = 0
)
ggplot2::ggsave(
  filename = file.path(
    plot_dir,
    "mg_selected_cluster_selection_diagnostics.png"
  ),
  plot = selection_plot,
  width = 10,
  height = 6,
  dpi = 300
)
ggplot2::ggsave(
  filename = file.path(
    plot_dir,
    "mg_selected_cluster_selection_diagnostics.pdf"
  ),
  plot = selection_plot,
  width = 10,
  height = 6
)

retained_mask <- !as.character(sobj@meta.data[[cluster_column]]) %in%
  excluded_clusters
retained_cells <- colnames(sobj)[retained_mask]
if (length(retained_cells) == 0) {
  stop("No retained cells remain after cluster exclusion.", call. = FALSE)
}

base_subset <- subset(sobj, cells = retained_cells)
base_subset <- strip_previous_cluster_state(
  base_subset,
  source_cluster_column = cluster_column,
  extra_drop = c(module_score_cols, "Cdkn1b_selection_expression")
)
DefaultAssay(base_subset) <- "RNA"
if (!inherits(base_subset[["RNA"]], "Assay5")) {
  base_subset[["RNA"]] <- as(base_subset[["RNA"]], Class = "Assay5")
}

utils::data("mouse_cell_cycle_genes", package = "ESPI", envir = environment())
cell_cycle_genes <- mouse_cell_cycle_genes
exclusion_reasons <- stats::setNames(
  decision_table$exclusion_reasons[decision_table$exclude],
  excluded_clusters
)
preprocess_subset <- function(
  sobj,
  filter_cc,
  excluded_clusters,
  exclusion_reasons
) {
  sobj@misc$preprocessing$filtered_cell_cycle <- filter_cc
  sobj <- Seurat::FindVariableFeatures(sobj, nfeatures = 2000)
  if (filter_cc) {
    VariableFeatures(sobj) <- setdiff(
      VariableFeatures(sobj),
      cell_cycle_genes
    )
  }
  sobj <- run_pflog_pca(sobj, n_pcs = dims)
  sobj@misc$preprocessing$dataset_tag <- dataset_tag
  sobj@misc$preprocessing$source_cluster_column <- cluster_column
  sobj@misc$preprocessing$source_input <- input_path
  sobj@misc$preprocessing$source_cluster_selection_table <- selection_path
  sobj@misc$preprocessing$source_cluster_selection_figure <- file.path(
    plot_dir,
    "mg_selected_cluster_selection_diagnostics.png"
  )
  sobj@misc$preprocessing$source_cluster_excluded <- excluded_clusters
  sobj@misc$preprocessing$source_cluster_exclusion_reasons <- exclusion_reasons
  sobj@misc$preprocessing$cdkn1b_expression_layer <- cdkn1b_expression_layer
  sobj@misc$preprocessing$marker_score_slot <- score_slot
  sobj@misc$preprocessing$marker_exclude_classes <- MARKER_EXCLUDE_CLASSES
  sobj@misc$preprocessing$marker_min_top_score <- MARKER_MIN_TOP_SCORE
  sobj@misc$preprocessing$marker_min_score_margin <- MARKER_MIN_SCORE_MARGIN
  sobj@misc$preprocessing$cdkn1b_expression_max_q <- CDKN1B_EXPRESSION_MAX_Q
  sobj@misc$preprocessing$cdkn1b_detection_max_q <- CDKN1B_DETECTION_MAX_Q
  sobj@misc$preprocessing$cdkn1b_min_detection_fraction <-
    CDKN1B_MIN_DETECTION_FRACTION
  sobj@misc$preprocessing$filtered_cell_cycle <- filter_cc
  cc_tag <- if (filter_cc) "filter_cc" else "no_filter_cc"
  elbow_plot <- Seurat::ElbowPlot(sobj, ndims = dims, reduction = "pca")
  elbow_png_path <- file.path(
    plot_dir,
    sprintf("elbow_pflog_%s_%s.png", dataset_tag, cc_tag)
  )
  elbow_pdf_path <- file.path(
    plot_dir,
    sprintf("elbow_pflog_%s_%s.pdf", dataset_tag, cc_tag)
  )
  ggplot2::ggsave(
    elbow_png_path,
    elbow_plot,
    width = 5,
    height = 3,
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    elbow_pdf_path,
    elbow_plot,
    width = 5,
    height = 3,
    bg = "white"
  )
  sobj@misc$preprocessing$elbow_plot <- elbow_png_path
  sobj
}

no_filter_subset <- preprocess_subset(
  base_subset,
  filter_cc = FALSE,
  excluded_clusters = excluded_clusters,
  exclusion_reasons = exclusion_reasons
)
filter_cc_subset <- preprocess_subset(
  base_subset,
  filter_cc = TRUE,
  excluded_clusters = excluded_clusters,
  exclusion_reasons = exclusion_reasons
)

dir.create(CURRENT_OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
no_filter_path <- file.path(
  CURRENT_OBJECT_DIR,
  sprintf("preprocess_pflog_%s_no-filter-cc.rds", dataset_tag)
)
filter_cc_path <- file.path(
  CURRENT_OBJECT_DIR,
  sprintf("preprocess_pflog_%s_filter-cc.rds", dataset_tag)
)
saveRDS(no_filter_subset, no_filter_path)
saveRDS(filter_cc_subset, filter_cc_path)

emit_tripwire_checkpoint(
  "mg_selected_subset_written",
  input = input_path,
  cluster_column = cluster_column,
  dataset_tag = dataset_tag,
  n_source_cells = ncol(sobj),
  n_retained_cells = length(retained_cells),
  excluded_clusters = paste(excluded_clusters, collapse = ","),
  no_filter_output = no_filter_path,
  filter_cc_output = filter_cc_path,
  selection_table = selection_path
)
message("Saved ", no_filter_path)
message("Saved ", filter_cc_path)
message("Wrote ", selection_path)
