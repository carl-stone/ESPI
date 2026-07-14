#!/usr/bin/env Rscript

# Read six 10X count directories and sample metadata into one combined Seurat object.
#
# Usage:
#   Rscript scripts/01-process-counts.R
#   Rscript scripts/01-process-counts.R \
#     [--raw-counts-dir <path>] [--metadata <path>] [--output <path>]
#
# Arguments:
#   --raw-counts-dir
#     Optional expert/test-only override for the directory containing the
#     per-sample 10X count directories.
#   --metadata
#     Optional expert/test-only override for the sample metadata file.
#   --output
#     Optional expert/test-only override for the output Seurat RDS path.
#
# Outputs:
#   DATA_ROOT_DIR/data/input/sobj_raw.rds (or the explicit --output path)
#
# Next step:
#   Run scripts/02-qc-filtering.R.

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
})
here::i_am("scripts/01-process-counts.R")
suppressPackageStartupMessages({
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
})

# ---- parameters ----

arguments <- commandArgs(trailingOnly = TRUE)
value_flags <- c("--raw-counts-dir", "--metadata", "--output")
parsed_arguments <- list()
argument_index <- 1L
while (argument_index <= length(arguments)) {
  argument <- arguments[[argument_index]]
  if (!argument %in% value_flags) {
    stop("Unknown argument: ", argument, ".", call. = FALSE)
  }
  if (!is.null(parsed_arguments[[argument]])) {
    stop("Duplicate argument: ", argument, ".", call. = FALSE)
  }
  if (
    argument_index == length(arguments) ||
      startsWith(arguments[[argument_index + 1L]], "--")
  ) {
    stop("Missing value for ", argument, ".", call. = FALSE)
  }
  parsed_arguments[[argument]] <- arguments[[argument_index + 1L]]
  if (!nzchar(arguments[[argument_index + 1L]])) {
    stop("Empty value for ", argument, ".", call. = FALSE)
  }
  argument_index <- argument_index + 2L
}

default_raw_counts_dir <- file.path(
  DATA_ROOT_DIR,
  "data",
  "input",
  "Raw Matrices"
)
raw_counts_dir <- if (is.null(parsed_arguments[["--raw-counts-dir"]])) {
  default_raw_counts_dir
} else {
  parsed_arguments[["--raw-counts-dir"]]
}
metadata_path <- if (is.null(parsed_arguments[["--metadata"]])) {
  file.path(raw_counts_dir, "Sample_Metadata_MS1.txt")
} else {
  parsed_arguments[["--metadata"]]
}
output_path <- if (is.null(parsed_arguments[["--output"]])) {
  file.path(DATA_ROOT_DIR, "data", "input", "sobj_raw.rds")
} else {
  parsed_arguments[["--output"]]
}
output_supplied <- !is.null(parsed_arguments[["--output"]])
required_metadata_columns <- c("Sample", "Mouse", "Condition")

# ---- validation ----

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

# Validate required columns and values before normalizing treatment labels.
validate_required_metadata(metadata, required_metadata_columns)

# Raw metadata contains an optional space between `+` and `EStim`; use the
# declared treatment label consistently in all downstream contrasts.
metadata$Condition <- sub("\\+\\s+EStim$", "+EStim", trimws(metadata$Condition))
metadata <- metadata[, required_metadata_columns, drop = FALSE]

if (anyDuplicated(metadata$Sample)) {
  stop("Sample metadata contains duplicate Sample IDs.", call. = FALSE)
}

emit_tripwire_checkpoint(
  "raw_data_available",
  raw_counts_dir = raw_counts_dir,
  metadata_path = metadata_path,
  n_metadata_rows = nrow(metadata)
)

count_directories <- list.files(
  raw_counts_dir,
  full.names = FALSE,
  recursive = FALSE,
  include.dirs = TRUE,
  no.. = TRUE
)
count_directories <- count_directories[
  dir.exists(file.path(raw_counts_dir, count_directories))
]
if (!base::setequal(metadata$Sample, count_directories)) {
  mismatched_samples <- sort(unique(c(
    base::setdiff(metadata$Sample, count_directories),
    base::setdiff(count_directories, metadata$Sample)
  )))
  write_tripwire_drop_ledger(
    sample_ids = mismatched_samples,
    stage = "samples_reconciled",
    reason = "10X count folder names do not match Sample metadata IDs",
    allowed = FALSE
  )
  stop(
    "10X count folder names do not match Sample metadata IDs.",
    call. = FALSE
  )
}
emit_tripwire_checkpoint(
  "samples_reconciled",
  n_samples = nrow(metadata),
  sample_ids = paste(sort(metadata$Sample), collapse = ",")
)

# ---- work ----

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

seurat_object <- Seurat::CreateSeuratObject(
  counts = tenx_counts,
  project = "ESPI",
  meta.data = cell_metadata,
  min.cells = 1,
  min.features = 1
)

# ---- output ----

message(
  "Created a Seurat object with ",
  ncol(seurat_object),
  " cells and ",
  nrow(seurat_object),
  " Gene Expression features across ",
  nrow(metadata),
  " samples."
)

if (output_supplied) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
}
saveRDS(seurat_object, output_path)

message(
  "Saved raw Seurat object to ",
  output_path,
  ". Next step: preprocess this file."
)
