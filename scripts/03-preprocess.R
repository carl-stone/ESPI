#!/usr/bin/env Rscript

# Preprocess one Seurat object into one normalization branch.
#
# Usage:
#   Rscript scripts/03-preprocess.R \
#     --input <seurat-object.rds> | --input-source <legacy|counts-qc> \
#     --normalization <log1p|pflog> \
#     --filter-cell-cycle <true|false>
# Test-only label permutation:
#   ESPI_PERMUTE_CONDITION_SEED=<integer>
#     Requires ESPI_TRIPWIRE_MODE=true, CHECKPOINT_LOG, and
#     STOP_AFTER_CHECKPOINT=blind_qc_complete. The permutation is applied only
#     after metadata validation and sample_id derivation, and no persistent
#     output is written in this mode.
# Arguments:
#   --input
#     Explicit Seurat object to preprocess. Cannot be combined with
#     --input-source.
#   --input-source
#     Named source object: legacy selects
#     INPUT_OBJECT_DIR/pipseq_processed_matrix_with_egfp.rds; counts-qc selects
#     INPUT_OBJECT_DIR/sobj_qc_filtered.rds. Defaults to legacy.
#   --normalization
#     Normalization branch to run. Must be log1p or pflog. Defaults to log1p.
#   --filter-cell-cycle
#     If true, removes mouse_cell_cycle_genes from the HVG set before PCA.
#     Defaults to false.
#
# Outputs:
#   CURRENT_OBJECT_DIR/preprocess_<normalization>_<filter-cc|no-filter-cc>.rds
#   QC, HVG, DimHeatmap, and elbow figures under FIGURE_DIR/preprocess.
#
# Next step:
#   Run scripts/04-cluster.R with --input <preprocess output> and --elbow-n.

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
  library(scclrR)
})
here::i_am("scripts/03-preprocess.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

args <- commandArgs(trailingOnly = TRUE)
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
arg_flag <- function(name) {
  identical(arg(name), TRUE)
}

input <- arg_value("--input", default = NULL)
input_source <- arg_value("--input-source", default = "legacy")
if (!is.null(input) && !is.null(arg("--input-source"))) {
  stop("Use either --input or --input-source, not both.", call. = FALSE)
}
if (!input_source %in% c("legacy", "counts-qc")) {
  stop(
    "--input-source must be one of legacy or counts-qc; got ",
    input_source,
    call. = FALSE
  )
}
if (is.null(input)) {
  input <- if (identical(input_source, "legacy")) {
    file.path(INPUT_OBJECT_DIR, "pipseq_processed_matrix_with_egfp.rds")
  } else {
    file.path(INPUT_OBJECT_DIR, "sobj_qc_filtered.rds")
  }
} else {
  input_source <- "custom"
}
normalization <- arg_value("--normalization", default = "log1p")
filter_cc <- arg_flag("--filter-cell-cycle") ||
  identical(
    tolower(arg_value("--filter-cell-cycle", default = "false")),
    "true"
  )
permute_condition_seed_text <- Sys.getenv(
  "ESPI_PERMUTE_CONDITION_SEED",
  unset = ""
)
tripwire_mode <- identical(
  Sys.getenv("ESPI_TRIPWIRE_MODE", unset = ""),
  "true"
)
permute_condition <- !identical(permute_condition_seed_text, "")
permute_condition_seed <- NULL
if (permute_condition) {
  if (!grepl("^[+]?[0-9]+$", permute_condition_seed_text)) {
    stop(
      "ESPI_PERMUTE_CONDITION_SEED must be a non-empty non-negative integer.",
      call. = FALSE
    )
  }
  parsed_seed <- suppressWarnings(as.numeric(permute_condition_seed_text))
  if (
    length(parsed_seed) != 1L ||
      !is.finite(parsed_seed) ||
      parsed_seed != floor(parsed_seed) ||
      parsed_seed > .Machine$integer.max
  ) {
    stop(
      "ESPI_PERMUTE_CONDITION_SEED must be a non-empty non-negative integer.",
      call. = FALSE
    )
  }
  permute_condition_seed <- as.integer(parsed_seed)
  if (
    !tripwire_mode ||
      identical(Sys.getenv("CHECKPOINT_LOG", unset = ""), "") ||
      !identical(
        Sys.getenv("STOP_AFTER_CHECKPOINT", unset = ""),
        "blind_qc_complete"
      )
  ) {
    stop(
      paste(
        "ESPI_PERMUTE_CONDITION_SEED requires",
        "ESPI_TRIPWIRE_MODE=true, a non-empty CHECKPOINT_LOG, and",
        "STOP_AFTER_CHECKPOINT=blind_qc_complete.",
        "Refusing to mutate production labels."
      ),
      call. = FALSE
    )
  }
}

# ---- validation ----

if (!normalization %in% c("log1p", "pflog")) {
  stop("Unknown normalization: ", normalization, call. = FALSE)
}
if (!file.exists(input)) {
  stop("Input Seurat object does not exist: ", input, call. = FALSE)
}

# ---- work ----

sobj <- readRDS(input)
emit_tripwire_checkpoint(
  "raw_data_available",
  input = input,
  input_source = input_source,
  n_cells = ncol(sobj),
  n_features = nrow(sobj)
)
validate_required_metadata(sobj@meta.data, c("Mouse", "Condition"))

# Drop any pre-existing reductions from upstream Trailmaker export.
sobj@reductions <- list()
sobj[["RNA"]] <- as(sobj[["RNA"]], Class = "Assay5")

sobj$sample_id <- paste0(
  "M",
  sobj$Mouse,
  "_",
  gsub("[^A-Za-z0-9]", "", sobj$Condition)
)
if (permute_condition) {
  metadata_before <- sobj[[]]
  cell_sample_ids <- as.character(metadata_before$sample_id)
  sample_table <- unique(data.frame(
    sample_id = cell_sample_ids,
    Mouse = as.character(metadata_before$Mouse),
    Condition = as.character(metadata_before$Condition),
    stringsAsFactors = FALSE
  ))
  if (nrow(sample_table) != length(unique(cell_sample_ids))) {
    stop(
      "Cannot permute Condition: sample_id maps to multiple Mouse × Condition samples.",
      call. = FALSE
    )
  }

  if (length(unique(sample_table$Condition)) < 2L) {
    stop(
      "Cannot permute Condition: at least two labels are required.",
      call. = FALSE
    )
  }
  set.seed(permute_condition_seed)
  permuted_by_sample <- sample(sample_table$Condition)
  if (identical(permuted_by_sample, sample_table$Condition)) {
    rotation <- c(seq.int(2L, nrow(sample_table)), 1L)
    permuted_by_sample <- sample_table$Condition[rotation]
  }
  if (identical(permuted_by_sample, sample_table$Condition)) {
    stop(
      "Condition permutation did not change any sample-level labels.",
      call. = FALSE
    )
  }
  if (
    !setequal(
      unique(permuted_by_sample),
      unique(sample_table$Condition)
    )
  ) {
    stop(
      "Condition permutation changed the represented label set.",
      call. = FALSE
    )
  }

  permuted_condition <- as.character(metadata_before$Condition)
  for (sample_row in seq_len(nrow(sample_table))) {
    cell_rows <- which(
      cell_sample_ids == sample_table$sample_id[[sample_row]]
    )
    permuted_condition[cell_rows] <- permuted_by_sample[[sample_row]]
  }
  if (is.factor(metadata_before$Condition)) {
    permuted_condition <- factor(
      permuted_condition,
      levels = levels(metadata_before$Condition)
    )
  }
  sobj$Condition <- permuted_condition

  metadata_after <- sobj[[]]
  other_metadata_columns <- setdiff(
    colnames(metadata_before),
    "Condition"
  )
  other_metadata_unchanged <- all(vapply(
    other_metadata_columns,
    function(column) {
      identical(metadata_before[[column]], metadata_after[[column]])
    },
    logical(1)
  ))
  if (
    !identical(rownames(metadata_before), rownames(metadata_after)) ||
      !identical(
        colnames(metadata_before),
        colnames(metadata_after)
      ) ||
      !other_metadata_unchanged
  ) {
    stop(
      "Condition permutation changed cell identities or non-Condition metadata.",
      call. = FALSE
    )
  }
  if (
    !identical(
      as.character(metadata_before$sample_id),
      as.character(metadata_after$sample_id)
    ) ||
      !identical(
        as.character(metadata_before$Mouse),
        as.character(metadata_after$Mouse)
      )
  ) {
    stop(
      "Condition permutation changed sample_id or Mouse metadata.",
      call. = FALSE
    )
  }
  condition_counts <- tapply(
    as.character(metadata_after$Condition),
    as.character(metadata_after$sample_id),
    function(values) length(unique(values))
  )
  if (any(condition_counts != 1L)) {
    stop(
      "Condition permutation did not preserve sample-level labels.",
      call. = FALSE
    )
  }
}

sobj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(
  sobj,
  pattern = "^Rp[sl]"
)

sobj@misc$preprocessing$input_source <- input_source
sobj@misc$preprocessing$normalization <- normalization
sobj@misc$preprocessing$filtered_cell_cycle <- filter_cc

if (!tripwire_mode) {
  splot_qc_metrics_violin(sobj)
}
if (tripwire_mode) {
  # Keep baseline and permutation tripwire runs on the same deterministic RNG
  # stream without changing the production path when test mode is unset.
  set.seed(SEED)
}


sobj <- FindVariableFeatures(sobj, nfeatures = 2000)

if (filter_cc) {
  utils::data("mouse_cell_cycle_genes", package = "ESPI", envir = environment())
  VariableFeatures(sobj) <- setdiff(
    VariableFeatures(sobj),
    mouse_cell_cycle_genes
  )
}
emit_tripwire_checkpoint(
  "variable_features_selected",
  normalization = normalization,
  filtered_cell_cycle = filter_cc,
  n_variable_features = length(VariableFeatures(sobj))
)

if (!tripwire_mode) {
  splot_hvg_scatter(sobj, n_top = 20)
}

sobj <- switch(
  normalization,
  log1p = run_log1p_pca(sobj, n_pcs = 50),
  pflog = run_pflog_pca(sobj, n_pcs = 50)
)
emit_tripwire_checkpoint(
  "pca_ready",
  normalization = normalization,
  filtered_cell_cycle = filter_cc,
  n_pcs = ncol(SeuratObject::Embeddings(sobj, reduction = "pca"))
)
emit_tripwire_checkpoint(
  "blind_qc_complete",
  fingerprint_algorithm = "exact-v1",
  hvg_fingerprint = paste(
    sort(SeuratObject::VariableFeatures(sobj)),
    collapse = ","
  ),
  pca_sdev_fingerprint = paste(
    as.character(
      round(
        SeuratObject::Stdev(sobj[["pca"]]),
        digits = 10
      )
    ),
    collapse = ","
  ),
  normalization = normalization,
  filtered_cell_cycle = filter_cc,
  plots_skipped = tripwire_mode
)

splot_dim_heatmap(sobj)
splot_elbow(sobj, n_pcs = 50)

# ---- output ----

cc_tag <- if (filter_cc) "filter-cc" else "no-filter-cc"
out_path <- file.path(
  CURRENT_OBJECT_DIR,
  sprintf("preprocess_%s_%s.rds", normalization, cc_tag)
)

dir.create(CURRENT_OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(sobj, out_path)
emit_tripwire_checkpoint(
  "preprocess_object_written",
  output = out_path,
  normalization = normalization,
  filtered_cell_cycle = cc_tag
)

message(
  "Saved ",
  out_path,
  ". Inspect the elbow plot, choose elbow_n, then run: ",
  "Rscript scripts/04-cluster.R --input ",
  out_path,
  " --elbow-n <N>"
)
