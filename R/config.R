# Shared labels, plotting defaults, and the analysis settings to edit.

SEED <- 1312L
CONDITION_COL <- "Condition"
CTRL_LABEL <- "p27CKO"
ESTIM_LABEL <- "p27CKO +EStim"
CTRL_DISPLAY_LABEL <- "p27CKO"
ESTIM_DISPLAY_LABEL <- "p27CKO + E-Stim"
CONTRAST_DISPLAY_LABEL <- sprintf(
  "(%s vs. %s)",
  ESTIM_DISPLAY_LABEL,
  CTRL_DISPLAY_LABEL
)
palette_analysis_three <- c(low = "#2166ac", mid = "grey75", high = "#e31a8c")
palette_dotplot_pair <- unname(palette_analysis_three[c("low", "high")])

#' Theme for publication figures
#' @param base_size Base font size.
#' @return A ggplot2 theme.
#' @export
theme_stone <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", color = "black"),
      text = ggplot2::element_text(color = "black")
    )
}

#' Analysis paths and selected clustering settings
#'
#' Edit the `selected` list here to change clustering choices. Object paths and
#' metadata column names are derived from those settings, not maintained separately.
#' The source and cell-cycle-filtered sensitivity retain seed 1312; the selected
#' MG clustering uses seed 2847.
#'
#' @return A list of paths, conditions, selected clusterings, grid, and palettes.
#' @export
publication_config <- function() {
  project_root <- here::here()
  local_config <- new.env(parent = baseenv())
  local_path <- file.path(project_root, "config.local.R")
  if (file.exists(local_path)) {
    sys.source(local_path, envir = local_config)
  }
  data_root <- if (!is.null(local_config$MEGAN_SC_DATA_DIR)) {
    local_config$MEGAN_SC_DATA_DIR
  } else if (!is.null(local_config$BOX_PATH)) {
    file.path(local_config$BOX_PATH, "megan_sc_data")
  } else {
    "~/Library/CloudStorage/Box-Box/megan_sc_data"
  }
  data_root <- path.expand(data_root)
  object_dir <- file.path(data_root, "seurat_objects")

  # These are the settings used for downstream analysis and plot labels.
  selected <- list(
    source = list(
      branch = "pflog_no_filter_cc",
      dimensions = 20L,
      resolution = 0.3
    ),
    mg = list(
      branch = "pflog_mg_selected_no_filter_cc",
      dimensions = 20L,
      resolution = 0.5,
      seed = 2847L
    ),
    mg_filter_cc = list(
      branch = "pflog_mg_selected_filter_cc",
      dimensions = 20L,
      resolution = 0.5
    )
  )
  for (name in names(selected)) {
    setting <- selected[[name]]
    selected[[name]]$column <- paste0(
      "cluster_",
      setting$branch,
      "_dims",
      setting$dimensions,
      "_res",
      setting$resolution
    )
    selected[[name]]$path <- file.path(
      object_dir,
      "current",
      paste0("cluster_", setting$branch, "_elbow20.rds")
    )
  }

  list(
    seed = SEED,
    paths = list(
      project = project_root,
      data = data_root,
      input_objects = file.path(object_dir, "input"),
      current_objects = file.path(object_dir, "current"),
      figures = file.path(data_root, "figures"),
      tables = file.path(data_root, "tables"),
      degs = file.path(data_root, "degs"),
      enrichment = file.path(data_root, "enrichment")
    ),
    conditions = list(
      column = CONDITION_COL,
      control = CTRL_LABEL,
      estim = ESTIM_LABEL,
      control_display = CTRL_DISPLAY_LABEL,
      estim_display = ESTIM_DISPLAY_LABEL,
      contrast_display = CONTRAST_DISPLAY_LABEL
    ),
    selected = selected,
    grid = list(dimensions = c(20L, 30L, 50L), resolutions = c(0.3, 0.5, 0.8)),
    palettes = list(
      analysis = palette_analysis_three,
      dotplot = palette_dotplot_pair
    )
  )
}

#' Prepare an output path without silently overwriting an existing file
#'
#' Use this at the write, rather than maintaining a separate list of outputs.
#' Existing files require `ESPI_OVERWRITE=true` or an explicit `overwrite = TRUE`.
#' Symlink destinations are rejected. Parent directories are created as needed.
#'
#' @param ... Path components, passed to [file.path()].
#' @param overwrite Whether existing files may be replaced.
#' @return The output path(s).
#' @export
output_path <- function(
  ...,
  overwrite = identical(tolower(Sys.getenv("ESPI_OVERWRITE", "false")), "true")
) {
  paths <- file.path(...)
  if (any(nzchar(Sys.readlink(paths), keepNA = TRUE), na.rm = TRUE)) {
    stop("Refusing to write through a symlink: ", paste(paths, collapse = ", "))
  }
  if (!overwrite && any(file.exists(paths))) {
    stop(
      "Output exists; set ESPI_OVERWRITE=true to replace: ",
      paste(paths, collapse = ", ")
    )
  }
  for (directory in unique(dirname(paths))) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  }
  paths
}

#' Mouse cell-cycle genes
#' Mouse orthologs of Seurat's human cell-cycle genes, generated with
#' `data-raw/mouse-cell-cycle-genes.R` using `biomaRt`.
#' @docType data
#' @keywords datasets
#' @format A character vector of mouse gene symbols.
"mouse_cell_cycle_genes"

#' Cell type marker genes
#' Curated marker genes for broad retinal cell type annotation.
#' @docType data
#' @keywords datasets
#' @format A named list of mouse gene symbols.
"cell_type_marker_genes"

#' Cell type marker labels
#' Display labels for `cell_type_marker_genes`.
#' @docType data
#' @keywords datasets
#' @format A named character vector.
"cell_type_marker_labels"
