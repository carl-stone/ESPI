#!/usr/bin/env Rscript

# Load Megan's raw 10X counts and sample metadata into one in-memory Seurat object.
#
# Usage:
#   Rscript scripts/process-counts.R
#
# Arguments:
#   None. Input paths derive from DATA_ROOT_DIR.
#
# Outputs:
#   None. This stage validates and prints an in-memory Seurat object only.
#
# Next step:
#   Save or preprocess the validated object in a later pipeline stage.

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
})
here::i_am("scripts/process-counts.R")
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
  "Loaded ",
  ncol(seurat_object),
  " cells and ",
  nrow(seurat_object),
  " Gene Expression features across ",
  nrow(metadata),
  " samples into a Seurat object."
)

saveRDS(
  seurat_object,
  file.path(DATA_ROOT_DIR, "data", "input", "sobj_raw.rds")
)
