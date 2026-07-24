suppressPackageStartupMessages({
  here::i_am("scripts/findmarkers.R")
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(tidyverse)
  library(Seurat)
  library(ggview)
})

# ---- parameters ----

config <- publication_config()
input_path <- config$selected$mg$path
cluster_column <- config$selected$mg$column
assay <- "RNA"
expression_layer <- "data"
counts_layer <- "counts"

# --- inputs ---
sobj <- readRDS(input_path)
cluster_values <- as.character(sobj[[cluster_column, drop = TRUE]])
identity_levels <- as.character(sort(as.integer(unique(cluster_values))))
marker_identities <- factor(cluster_values, levels = identity_levels)
SeuratObject::Idents(sobj) <- marker_identities

# --- cluster markers ---

Idents(sobj, cells = WhichCells(sobj, idents = "2")) <- "1"

cluster_markers <- FindAllMarkers(
  sobj,
  assay = assay,
  test.use = "wilcox",
  logfc.threshold = 0.1
)
