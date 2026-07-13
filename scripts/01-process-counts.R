#!/usr/bin/env Rscript

# Read six 10X count directories and sample metadata into one combined Seurat object.
#
# Usage:
#   Rscript scripts/01-process-counts.R
#
# Arguments:
#   None. Input paths derive from DATA_ROOT_DIR.
#
# Outputs:
#   DATA_ROOT_DIR/data/input/sobj_raw.rds
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

raw_counts_dir <- file.path(DATA_ROOT_DIR, "data", "input", "Raw Matrices")
metadata_path <- file.path(raw_counts_dir, "Sample_Metadata_MS1.txt")
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

missing_columns <- setdiff(required_metadata_columns, colnames(metadata))
if (length(missing_columns) > 0) {
  stop(
    "Sample metadata is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}
metadata <- metadata[, required_metadata_columns, drop = FALSE]

# Raw metadata contains an optional space between `+` and `EStim`; use the
# declared treatment label consistently in all downstream contrasts.
metadata$Condition <- sub("\\+\\s+EStim$", "+EStim", trimws(metadata$Condition))

if (any(is.na(metadata$Sample) | !nzchar(trimws(metadata$Sample)))) {
  stop("Sample metadata contains missing or empty Sample IDs.", call. = FALSE)
}
if (anyDuplicated(metadata$Sample)) {
  stop("Sample metadata contains duplicate Sample IDs.", call. = FALSE)
}

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
if (!setequal(metadata$Sample, count_directories)) {
  stop(
    "10X count folder names do not match Sample metadata IDs.",
    call. = FALSE
  )
}

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
  min.cells = 0,
  min.features = 0
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

saveRDS(
  seurat_object,
  file.path(DATA_ROOT_DIR, "data", "input", "sobj_raw.rds")
)

message(
  "Saved raw Seurat object to ",
  file.path(DATA_ROOT_DIR, "data", "input", "sobj_raw.rds"),
  ". Next step: preprocess this file."
)
