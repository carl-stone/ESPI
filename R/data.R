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
#' Curated marker genes for broad retinal cell type annotation in the ESPI
#' single-cell analysis. Ed and Megan selected these markers from domain
#' knowledge and literature review; the repo does not currently record the
#' marker-by-marker selection rationale.
#'
#' Generated from `data-raw/cell-type-marker-genes.R`.
#'
#' @docType data
#' @keywords datasets
#' @format A named list. Each element name is an ASCII cell type key, and each
#'   element is a character vector of mouse gene symbols.
"cell_type_marker_genes"

#' Cell type marker labels
#'
#' Display labels for `cell_type_marker_genes`. Names match the ASCII marker
#' list keys; values are human-readable labels for plots and reports.
#'
#' Generated from `data-raw/cell-type-marker-genes.R`.
#'
#' @docType data
#' @keywords datasets
#' @format A named character vector mapping marker-list keys to display labels.
"cell_type_marker_labels"
