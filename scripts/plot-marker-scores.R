suppressPackageStartupMessages({
  here::i_am("scripts/plot-marker-scores.R")
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(tidyverse)
  library(Seurat)
  library(ggview)
})

# ---- parameters ----

config <- publication_config()
input_path <- config$selected$mg$path
cluster_column <- config$selected$mg$column
source_cluster_column <- "cluster_pflog_no_filter_cc_dims20_res0.3"
assay <- "RNA"
expression_layer <- "data"
counts_layer <- "counts"

# ---- inputs ----

source_sobj <- readRDS(config$selected$source$path)
sobj <- readRDS(input_path)
cluster_values <- as.character(sobj[[cluster_column, drop = TRUE]])
identity_levels <- as.character(sort(as.integer(unique(cluster_values))))
marker_identities <- factor(cluster_values, levels = identity_levels)
SeuratObject::Idents(sobj) <- marker_identities
Idents(source_sobj) <- source_sobj$cluster_pflog_no_filter_cc_dims20_res0.3

marker_table <- stack(cell_type_marker_genes) |>
  tibble::as_tibble() |>
  dplyr::rename(gene = values, cell_type = ind) |>
  dplyr::mutate(
    cell_type = as.character(cell_type),
    cell_type_label = unname(cell_type_marker_labels[cell_type])
  ) |>
  dplyr::bind_rows(tibble::tibble(
    gene = "Cdkn1b",
    cell_type = "cdkn1b_standalone",
    cell_type_label = "Cdkn1b"
  ))

# ---- plot UMAPS ---
#
# DimPlot(sobj)
# DimPlot(source_sobj)

# --- plot marker scores ---

scoreMarkerList <- function(sobj) {
  sobj <- AddModuleScore(
    object = sobj,
    features = cell_type_marker_genes,
    name = cell_type_marker_labels,
    assay = "RNA",
    seed = config$seed
  )
}

source_sobj <- scoreMarkerList(source_sobj)

full_marker_violin_plot <- VlnPlot(
  source_sobj,
  features = c(paste0(cell_type_marker_labels, 1:11), "Cdkn1b"),
  combine = FALSE
) |>
  lapply(\(x) {
    x +
      labs(x = "Cluster identity") +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  })

full_marker_violin_plot[[12]] + canvas(4, 4)

purrr::iwalk(full_marker_violin_plot, \(p, i) {
  ggsave(
    filename = file.path(
      config$paths$figures,
      "full_marker_score_violins",
      paste0(names(cell_type_marker_labels)[i], ".pdf")
    ),
    plot = p,
    width = 4,
    height = 4
  )
})

ggsave(
  filename = file.path(
    config$paths$figures,
    "full_marker_score_violins",
    "p27.pdf"
  ),
  plot = full_marker_violin_plot[[12]],
  width = 4,
  height = 4
)

sobj <- scoreMarkerList(sobj)

mg_marker_scores <- FeaturePlot(
  sobj,
  features = c(paste0(cell_type_marker_labels, 1:11), "Cdkn1b"),
  stroke.size = NULL,
  pt.size = 0.2,
  order = TRUE,
  combine = FALSE
) |>
  lapply(\(x) {
    x +
      labs(x = "UMAP1", y = "UMAP2") +
      theme(
        legend.position = "inside",
        legend.position.inside = c(1, 0.5),
        legend.justification = c("right")
      )
  })

mg_marker_scores[[12]] + canvas(4, 4)

purrr::iwalk(mg_marker_scores, \(p, i) {
  ggsave(
    filename = file.path(
      config$paths$figures,
      "mg_marker_score_umaps",
      paste0(c(names(cell_type_marker_labels), "p27")[i], ".pdf")
    ),
    plot = p,
    width = 4,
    height = 4
  )
})
