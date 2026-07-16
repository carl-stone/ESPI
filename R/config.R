# Fixed paths, labels, palettes, and package-data documentation.

SEED <- 1312L
PROJECT_ROOT <- here::here()
LOCAL_CONFIG <- file.path(PROJECT_ROOT, "config.local.R")
local_config <- new.env(parent = baseenv())
if (file.exists(LOCAL_CONFIG)) {
  sys.source(LOCAL_CONFIG, envir = local_config)
}

if (exists("MEGAN_SC_DATA_DIR", envir = local_config, inherits = FALSE)) {
  DATA_ROOT_DIR <- get("MEGAN_SC_DATA_DIR", envir = local_config)
} else if (exists("BOX_PATH", envir = local_config, inherits = FALSE)) {
  DATA_ROOT_DIR <- file.path(
    get("BOX_PATH", envir = local_config),
    "megan_sc_data"
  )
} else {
  DATA_ROOT_DIR <- path.expand(
    "~/Library/CloudStorage/Box-Box/megan_sc_data"
  )
}
if (!dir.exists(DATA_ROOT_DIR)) {
  stop("Data root does not exist: ", DATA_ROOT_DIR, call. = FALSE)
}
DATA_ROOT_DIR <- normalizePath(DATA_ROOT_DIR, mustWork = FALSE)
OUTPUT_DIR <- DATA_ROOT_DIR
OBJECT_DIR <- file.path(DATA_ROOT_DIR, "seurat_objects")
DATA_DIR <- file.path(DATA_ROOT_DIR, "data")
FIGURE_DIR <- file.path(DATA_ROOT_DIR, "figures")
TABLE_DIR <- file.path(DATA_ROOT_DIR, "tables")
DEG_DIR <- file.path(DATA_ROOT_DIR, "degs")
ENRICHMENT_DIR <- file.path(DATA_ROOT_DIR, "enrichment")
LOG_DIR <- file.path(DATA_ROOT_DIR, "logs")
RENDERED_NOTEBOOK_DIR <- file.path(DATA_ROOT_DIR, "rendered_notebooks")
INPUT_OBJECT_DIR <- file.path(OBJECT_DIR, "input")
CURRENT_OBJECT_DIR <- file.path(OBJECT_DIR, "current")
PREPROCESSING_DATA_DIR <- file.path(DATA_DIR, "preprocessing")
ANALYSIS_DATA_DIR <- file.path(DATA_DIR, "analysis")

CONDITION_COL <- "Condition"
ESTIM_LABEL <- "p27CKO +EStim"
CTRL_LABEL <- "p27CKO"
ESTIM_DISPLAY_LABEL <- "p27CKO + E-Stim"
CTRL_DISPLAY_LABEL <- "p27CKO"
CONTRAST_DISPLAY_LABEL <- sprintf(
  "(%s vs. %s)",
  ESTIM_DISPLAY_LABEL,
  CTRL_DISPLAY_LABEL
)

palette_analysis_three <- c(
  low = "#2166ac",
  mid = "grey75",
  high = "#e31a8c"
)
palette_dotplot_pair <- unname(
  palette_analysis_three[c("low", "high")]
)

#' Theme for publication figures
#'
#' @param base_size Base font size.
#'
#' @return A ggplot2 theme.
#' @export
# ANALYSIS_OK[smuggled-default]: exported theme API keeps the publication base font size fixed.
# ANALYSIS_OK[R026]: exported theme helper is called by executable figure and regeneration scripts.
theme_stone <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", color = "black"),
      axis.text = ggplot2::element_text(face = "bold", color = "black"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )
}

#' Mouse cell-cycle genes
#'
#' Mouse orthologs of Seurat's human cell-cycle genes, generated with
#' `data-raw/mouse-cell-cycle-genes.R` using `biomaRt`.
#'
#' @docType data
#' @keywords datasets
#' @format A character vector of mouse gene symbols.
"mouse_cell_cycle_genes"

#' Cell type marker genes
#'
#' Curated marker genes for broad retinal cell type annotation.
#'
#' @docType data
#' @keywords datasets
#' @format A named list of mouse gene symbols.
"cell_type_marker_genes"

#' Cell type marker labels
#'
#' Display labels for `cell_type_marker_genes`.
#'
#' @docType data
#' @keywords datasets
#' @format A named character vector.
"cell_type_marker_labels"

#' Fixed publication-analysis configuration
#'
#' Returns the paths, labels, selected clustering contracts, palettes, and
#' execution policy for the publication pipeline.
#'
#' @return A nested list of publication-analysis configuration values.
#' @export
# ANALYSIS_OK[R026]: exported configuration entrypoint is called by all executable phase scripts.
publication_config <- function() {
  overwrite_value <- tolower(Sys.getenv("ESPI_OVERWRITE", unset = "false"))
  if (!overwrite_value %in% c("true", "false")) {
    stop("ESPI_OVERWRITE must be 'true' or 'false'.", call. = FALSE)
  }

  final_source <- file.path(
    CURRENT_OBJECT_DIR,
    "cluster_pflog_no_filter_cc_elbow20.rds"
  )
  final_mg <- file.path(
    CURRENT_OBJECT_DIR,
    "cluster_pflog_mg_selected_no_filter_cc_elbow20.rds"
  )
  sensitivity_mg <- file.path(
    CURRENT_OBJECT_DIR,
    "cluster_pflog_mg_selected_filter_cc_elbow20.rds"
  )

  list(
    seed = as.integer(SEED),
    overwrite = identical(overwrite_value, "true"),
    paths = list(
      project = PROJECT_ROOT,
      data = DATA_ROOT_DIR,
      objects = OBJECT_DIR,
      input_objects = INPUT_OBJECT_DIR,
      current_objects = CURRENT_OBJECT_DIR,
      figures = FIGURE_DIR,
      tables = TABLE_DIR,
      degs = DEG_DIR,
      enrichment = ENRICHMENT_DIR,
      notebook = file.path(PROJECT_ROOT, "notebook"),
      notebook_figures = file.path(PROJECT_ROOT, "notebook", "figures")
    ),
    conditions = list(
      column = CONDITION_COL,
      control = CTRL_LABEL,
      estim = ESTIM_LABEL,
      control_display = CTRL_DISPLAY_LABEL,
      estim_display = ESTIM_DISPLAY_LABEL,
      contrast_display = CONTRAST_DISPLAY_LABEL
    ),
    selected = list(
      source = list(
        path = final_source,
        column = "cluster_pflog_no_filter_cc_dims30_res0.3",
        dimensions = 30L,
        resolution = 0.3
      ),
      mg = list(
        path = final_mg,
        column = "cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5",
        dimensions = 20L,
        resolution = 0.5
      ),
      mg_filter_cc = list(
        path = sensitivity_mg,
        column = "cluster_pflog_mg_selected_filter_cc_dims20_res0.5",
        dimensions = 20L,
        resolution = 0.5
      )
    ),
    frozen = list(
      source = list(
        sha256 = "0c5079e7174a68eed9b62a02f4ec550fcc776ffa14e183e718799cbf9d9bf74a",
        cells = 4146L,
        column = "cluster_pflog_no_filter_cc_dims30_res0.3"
      ),
      mg = list(
        sha256 = "ed320802451f9f426b1ca1e72b5f8bcd6bbd6f1095693fec17253903e62edc67",
        cells = 3456L,
        column = "cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5"
      ),
      mg_filter_cc = list(
        sha256 = "093faaa3ab28e7dfb0dd4121cede133a4f24399f1b8c4a1a920043f689af2c62",
        cells = 3456L,
        column = "cluster_pflog_mg_selected_filter_cc_dims20_res0.5"
      )
    ),
    palettes = list(
      analysis = palette_analysis_three,
      dotplot = palette_dotplot_pair
    )
  )
}

#' Validate a frozen Seurat input
#'
#' @param path Path to the frozen RDS file.
#' @param sobj Loaded Seurat object from `path`.
#' @param contract Frozen contract from [publication_config()].
#'
#' @return `sobj`, invisibly.
#' @export
# ANALYSIS_OK[R026]: exported frozen-input guard is called by publication and marker/DE phase scripts.
assert_frozen_input <- function(path, sobj, contract) {
  actual_hash <- digest::digest(path, algo = "sha256", file = TRUE)
  if (!identical(actual_hash, contract$sha256)) {
    stop("Frozen input hash mismatch: ", path, call. = FALSE)
  }
  if (ncol(sobj) != contract$cells) {
    stop("Frozen input cell-count mismatch: ", path, call. = FALSE)
  }
  if (!contract$column %in% colnames(sobj[[]])) {
    stop("Frozen input lacks cluster column: ", contract$column, call. = FALSE)
  }
  invisible(sobj)
}

#' Refuse to replace fixed outputs without permission
#'
#' @param paths Fixed output paths owned by one phase.
#' @param overwrite Whether existing regular outputs may be replaced.
#'
#' @return `paths`, invisibly.
#' @export
# ANALYSIS_OK[R026]: exported output guard is called by every executable phase script.
assert_output_available <- function(paths, overwrite) {
  link_targets <- Sys.readlink(paths)
  linked <- paths[!is.na(link_targets) & nzchar(link_targets)]
  if (length(linked) > 0L) {
    stop(
      "Refusing to replace symlinked output: ",
      paste(linked, collapse = ", "),
      call. = FALSE
    )
  }
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0L && !isTRUE(overwrite)) {
    stop(
      "Output exists; set ESPI_OVERWRITE=true to replace: ",
      paste(existing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(paths)
}
