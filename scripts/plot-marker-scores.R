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

# ---- inputs ----

source_sobj <- readRDS(config$selected$source$path)
sobj <- readRDS(input_path)
cluster_values <- as.character(sobj[[cluster_column, drop = TRUE]])
identity_levels <- as.character(sort(as.integer(unique(cluster_values))))
marker_identities <- factor(cluster_values, levels = identity_levels)
SeuratObject::Idents(sobj) <- marker_identities
Idents(source_sobj) <- config$selected$source$column


# --- plot marker scores ---

scoreMarkerList <- function(sobj) {
  AddModuleScore(
    object = sobj,
    features = cell_type_marker_genes,
    name = cell_type_marker_labels,
    assay = "RNA",
    seed = config$seed
  )
}

source_sobj <- scoreMarkerList(source_sobj)
score_features <- c(
  paste0(cell_type_marker_labels, seq_along(cell_type_marker_labels)),
  "Cdkn1b"
)
marker_names <- c(names(cell_type_marker_labels), "p27")

full_marker_violin_plot <- VlnPlot(
  source_sobj,
  features = score_features,
  combine = FALSE
) |>
  lapply(\(x) {
    x +
      labs(x = "Cluster identity") +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  })

full_marker_violin_plot[[length(score_features)]] + canvas(4, 4)

purrr::iwalk(full_marker_violin_plot, \(p, i) {
  ggsave(
    filename = output_path(
      config$paths$figures,
      "full_marker_score_violins",
      paste0(marker_names[i], ".pdf")
    ),
    plot = p,
    width = 4,
    height = 4
  )
})


sobj <- scoreMarkerList(sobj)

mg_marker_scores <- FeaturePlot(
  sobj,
  features = score_features,
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

mg_marker_scores[[length(score_features)]] + canvas(4, 4)

purrr::iwalk(mg_marker_scores, \(p, i) {
  ggsave(
    filename = output_path(
      config$paths$figures,
      "mg_marker_score_umaps",
      paste0(marker_names[i], ".pdf")
    ),
    plot = p,
    width = 4,
    height = 4
  )
})
