# paths.R — Global constants and artifact paths

#' Global random seed used by analysis scripts.
#' @export
SEED <- 1312

# Machine-local artifact root.
#
# Create `config.local.R` from `config.local.example.R` to set either:
#   BOX_PATH <- "/path/to/Box-Box"
# or:
#   MEGAN_SC_DATA_DIR <- "/path/to/megan_sc_data"
#
# If neither is set, use Carl's local Box Drive path. If the resolved data root
# does not exist, fail immediately rather than writing to a different directory.
PROJECT_ROOT <- here::here()
LOCAL_CONFIG <- file.path(PROJECT_ROOT, "config.local.R")
if (file.exists(LOCAL_CONFIG)) {
  source(LOCAL_CONFIG, local = FALSE)
}

MEGAN_SC_DATA_DIR_VALUE <- NULL
if (exists("MEGAN_SC_DATA_DIR", inherits = FALSE)) {
  MEGAN_SC_DATA_DIR_VALUE <- get("MEGAN_SC_DATA_DIR", inherits = FALSE)
} else if (exists("MEGAN_SC_DATA_DIR", envir = .GlobalEnv, inherits = FALSE)) {
  MEGAN_SC_DATA_DIR_VALUE <- get("MEGAN_SC_DATA_DIR", envir = .GlobalEnv)
}

BOX_PATH_VALUE <- NULL
if (exists("BOX_PATH", inherits = FALSE)) {
  BOX_PATH_VALUE <- get("BOX_PATH", inherits = FALSE)
} else if (exists("BOX_PATH", envir = .GlobalEnv, inherits = FALSE)) {
  BOX_PATH_VALUE <- get("BOX_PATH", envir = .GlobalEnv)
}

DEFAULT_BOX_PATH <- path.expand("~/Library/CloudStorage/Box-Box")

if (!is.null(MEGAN_SC_DATA_DIR_VALUE)) {
  DATA_ROOT_DIR <- MEGAN_SC_DATA_DIR_VALUE
} else if (!is.null(BOX_PATH_VALUE)) {
  DATA_ROOT_DIR <- file.path(BOX_PATH_VALUE, "megan_sc_data")
} else {
  DATA_ROOT_DIR <- file.path(DEFAULT_BOX_PATH, "megan_sc_data")
}

if (!dir.exists(DATA_ROOT_DIR)) {
  stop("Data root does not exist: ", DATA_ROOT_DIR, call. = FALSE)
}

#' Root directory for generated project artifacts.
#' @export
DATA_ROOT_DIR <- normalizePath(DATA_ROOT_DIR, mustWork = FALSE)
#' Legacy output root; identical to `DATA_ROOT_DIR`.
#' @export
OUTPUT_DIR <- DATA_ROOT_DIR

#' Legacy object artifact root.
#' @export
OBJECT_DIR <- file.path(OUTPUT_DIR, "seurat_objects")
#' Legacy data artifact root.
#' @export
DATA_DIR <- file.path(OUTPUT_DIR, "data")
#' Figure artifact root.
#' @export
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
#' Table artifact root.
#' @export
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
#' DEG artifact root.
#' @export
DEG_DIR <- file.path(OUTPUT_DIR, "degs")
#' Enrichment artifact root.
#' @export
ENRICHMENT_DIR <- file.path(OUTPUT_DIR, "enrichment")
#' Log artifact root.
#' @export
LOG_DIR <- file.path(OUTPUT_DIR, "logs")
#' Rendered notebook artifact root.
#' @export
RENDERED_NOTEBOOK_DIR <- file.path(OUTPUT_DIR, "rendered_notebooks")

#' Input Seurat object directory.
#' @export
INPUT_OBJECT_DIR <- file.path(OBJECT_DIR, "input")
#' Current Seurat object directory.
#' @export
CURRENT_OBJECT_DIR <- file.path(OBJECT_DIR, "current")


#' Preprocessing data artifact directory.
#' @export
PREPROCESSING_DATA_DIR <- file.path(DATA_DIR, "preprocessing")
#' Analysis data artifact directory.
#' @export
ANALYSIS_DATA_DIR <- file.path(DATA_DIR, "analysis")
