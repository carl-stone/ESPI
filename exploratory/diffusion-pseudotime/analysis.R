#!/usr/bin/env Rscript

# Cluster-independent diffusion pseudotime rooted by MG minus neuronal identity.
# Preserve the MG-selected input's existing Cdkn1b-based contamination filtering.
# No new cell filtering or Cdkn1b scoring is performed here.
# Run from the project root: source("exploratory/diffusion-pseudotime/analysis.R")

suppressPackageStartupMessages({
  library(here)
  here::i_am("exploratory/diffusion-pseudotime/analysis.R")
  library(devtools)
  devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(viridisLite)
  library(Matrix)
  library(BiocNeighbors)
  library(igraph)
  library(destiny)
})
ggplot2::theme_set(theme_stone())

# ---- parameters ----

pt_config <- publication_config()
pt_input_path <- pt_config$selected$mg$path
pt_dimensions <- pt_config$selected$mg$dimensions
pt_seed <- 5938L
pt_k <- 60L
pt_n_eigs <- 30L
pt_output_dir <- Sys.getenv(
  "ESPI_PSEUDOTIME_OUTPUT_DIR",
  unset = here::here(
    "exploratory",
    "diffusion-pseudotime",
    "outputs",
    "mg-minus-neuron"
  )
)
pt_neuron_types <- c(
  "cone_bipolar",
  "rod_bipolar",
  "photoreceptor",
  "retinal_ganglion",
  "horizontal",
  "amacrine"
)
pt_signatures <- c(
  cell_type_marker_genes,
  list(neuron_union = unique(unlist(cell_type_marker_genes[pt_neuron_types])))
)

# ---- input validation and scores ----

pt_sobj <- readRDS(pt_input_path)
pt_metadata <- pt_sobj[[]]
stopifnot(
  all(
    c(
      "Mouse",
      "Condition",
      "Sample",
      "nCount_RNA",
      "nFeature_RNA",
      "pass_qc",
      "is_singlet"
    ) %in%
      names(pt_metadata)
  ),
  !anyNA(pt_metadata[, c(
    "Mouse",
    "Condition",
    "Sample",
    "nCount_RNA",
    "nFeature_RNA"
  )]),
  all(pt_metadata$pass_qc),
  all(pt_metadata$is_singlet),
  all(
    pt_metadata$Condition %in%
      c(pt_config$conditions$control, pt_config$conditions$estim)
  )
)
pt_expression <- SeuratObject::LayerData(pt_sobj, assay = "RNA", layer = "data")
pt_pca <- SeuratObject::Embeddings(pt_sobj, "pca")[,
  seq_len(pt_dimensions),
  drop = FALSE
]
stopifnot(
  identical(colnames(pt_expression), rownames(pt_metadata)),
  identical(rownames(pt_pca), rownames(pt_metadata)),
  !anyDuplicated(rownames(pt_metadata)),
  all(is.finite(pt_pca)),
  nrow(pt_pca) > max(pt_k, pt_n_eigs + 1L)
)
pt_marker_coverage <- purrr::imap(pt_signatures, function(genes, signature) {
  tibble::tibble(
    signature,
    gene = genes,
    present = genes %in% rownames(pt_expression)
  )
}) |>
  purrr::list_rbind()
if (!all(pt_marker_coverage$present)) {
  stop("Some curated markers are missing; review coverage before scoring.")
}

# Reuse the project's expression-matched control-gene module scoring.
set.seed(pt_seed)
pt_sobj <- Seurat::AddModuleScore(
  pt_sobj,
  features = pt_signatures,
  assay = "RNA",
  slot = "data",
  name = "trajectory_signature_",
  seed = pt_seed,
  search = FALSE
)
pt_scores <- as.matrix(pt_sobj[[]][, paste0(
  "trajectory_signature_",
  seq_along(pt_signatures)
)])
colnames(pt_scores) <- names(pt_signatures)
stopifnot(all(is.finite(pt_scores)), all(apply(pt_scores, 2L, stats::sd) > 0))
pt_scores_z <- scale(pt_scores)
pt_cells <- pt_metadata |>
  tibble::rownames_to_column("cell") |>
  dplyr::transmute(
    cell,
    Mouse = as.character(Mouse),
    Condition,
    Sample,
    sample_id = paste(Mouse, Condition, sep = "__"),
    nCount_RNA,
    nFeature_RNA
  ) |>
  dplyr::bind_cols(tibble::as_tibble(pt_scores))

# ---- overlapping PCA neighborhoods ----

# Each neighborhood contains its center plus k nearest cells, not a cluster.
pt_knn <- BiocNeighbors::findKNN(pt_pca, k = pt_k)
pt_neighbors <- cbind(seq_len(nrow(pt_pca)), pt_knn$index)
pt_weights <- Matrix::sparseMatrix(
  i = rep(seq_len(nrow(pt_pca)), each = ncol(pt_neighbors)),
  j = as.vector(t(pt_neighbors)),
  x = 1 / ncol(pt_neighbors),
  dims = c(nrow(pt_pca), nrow(pt_pca))
)
pt_nhood_scores <- as.matrix(pt_weights %*% pt_scores)
pt_nhood_scores_z <- as.matrix(pt_weights %*% pt_scores_z)
pt_sample_membership <- Matrix::sparseMatrix(
  i = seq_len(nrow(pt_cells)),
  j = as.integer(factor(pt_cells$sample_id)),
  x = 1,
  dims = c(nrow(pt_cells), nlevels(factor(pt_cells$sample_id))),
  dimnames = list(NULL, levels(factor(pt_cells$sample_id)))
)
pt_sample_counts <- round(
  as.matrix(pt_weights %*% pt_sample_membership) * ncol(pt_neighbors)
)
pt_nhoods <- tibble::tibble(
  center_cell = pt_cells$cell,
  n_cells = ncol(pt_neighbors),
  n_samples = rowSums(pt_sample_counts > 0),
  max_sample_fraction = apply(pt_sample_counts, 1L, max) / ncol(pt_neighbors),
  mean_RNA_counts = as.numeric(pt_weights %*% pt_cells$nCount_RNA),
  mean_RNA_features = as.numeric(pt_weights %*% pt_cells$nFeature_RNA)
) |>
  dplyr::bind_cols(tibble::as_tibble(pt_nhood_scores)) |>
  dplyr::mutate(
    mg_vs_neuron = pt_nhood_scores_z[, "muller_glia"] -
      pt_nhood_scores_z[, "neuron_union"]
  )

# Neuronal identity is the union of the curated retinal-neuron subtype markers.
# Late-state references are descriptive: they do not constrain the DPT fit.
pt_anchor_centers <- c(
  MG = which.max(pt_nhoods$mg_vs_neuron),
  progenitor = which.max(pt_nhoods$neurogenic_progenitor),
  neuronal = which.max(pt_nhoods$neuron_union)
)
pt_anchor_members <- lapply(pt_anchor_centers, function(i) pt_neighbors[i, ])
pt_anchor_medoids <- vapply(
  pt_anchor_members,
  function(ids) {
    ids[which.min(rowSums(as.matrix(stats::dist(pt_pca[ids, , drop = FALSE]))))]
  },
  integer(1)
)
pt_root <- unname(pt_anchor_medoids[["MG"]])
pt_anchor_summary <- pt_nhoods[pt_anchor_centers, ] |>
  dplyr::mutate(
    anchor = names(pt_anchor_centers),
    medoid_cell = pt_cells$cell[pt_anchor_medoids]
  ) |>
  dplyr::select(anchor, center_cell, medoid_cell, dplyr::everything())
pt_anchor_samples <- purrr::imap(pt_anchor_members, function(ids, anchor) {
  pt_cells[ids, ] |>
    dplyr::count(Mouse, Condition, Sample) |>
    dplyr::mutate(anchor)
}) |>
  purrr::list_rbind()

pt_graph <- igraph::graph_from_adjacency_matrix(
  pt_weights > 0,
  mode = "max",
  diag = FALSE
)
pt_components <- igraph::components(pt_graph)
if (pt_components$no != 1L) {
  stop(
    "PCA neighborhoods are disconnected; review components before fitting one trajectory."
  )
}

# ---- diffusion map and rooted pseudotime ----

set.seed(pt_seed)
pt_dm <- destiny::DiffusionMap(
  pt_pca,
  k = pt_k,
  n_eigs = pt_n_eigs,
  sigma = "local",
  density_norm = TRUE,
  n_pcs = NA
)
pt_dpt <- destiny::DPT(pt_dm, tips = pt_root)
# Numeric DPT matrix columns are distances from the indexed input cell.
pt_cells$dpt_raw <- as.numeric(pt_dpt[, pt_root])
stopifnot(
  identical(rownames(destiny::eigenvectors(pt_dm)), pt_cells$cell),
  length(pt_cells$dpt_raw) == nrow(pt_cells),
  all(is.finite(pt_cells$dpt_raw)),
  all(pt_cells$dpt_raw >= 0),
  pt_cells$dpt_raw[pt_root] == 0,
  max(pt_cells$dpt_raw) > 0
)
pt_cells$pseudotime <- pt_cells$dpt_raw / max(pt_cells$dpt_raw)
pt_cells$pseudotime_percentile <- dplyr::percent_rank(pt_cells$dpt_raw)
pt_cells$DC1 <- destiny::eigenvectors(pt_dm)[, 1]
pt_cells$DC2 <- destiny::eigenvectors(pt_dm)[, 2]
pt_cells$DC3 <- destiny::eigenvectors(pt_dm)[, 3]
pt_nhoods$mean_pseudotime <- as.numeric(pt_weights %*% pt_cells$pseudotime)
pt_anchor_times <- purrr::imap(pt_anchor_members, function(ids, anchor) {
  tibble::tibble(
    anchor,
    min_pseudotime = min(pt_cells$pseudotime[ids]),
    median_pseudotime = median(pt_cells$pseudotime[ids]),
    max_pseudotime = max(pt_cells$pseudotime[ids])
  )
}) |>
  purrr::list_rbind()
pt_anchor_summary <- dplyr::left_join(
  pt_anchor_summary,
  pt_anchor_times,
  by = "anchor"
)

# ---- descriptive checks, not cell-level hypothesis tests ----

pt_score_associations <- purrr::map(
  c(names(pt_signatures), "nCount_RNA", "nFeature_RNA"),
  function(variable) {
    tibble::tibble(
      variable,
      spearman_rho = stats::cor(
        pt_cells[[variable]],
        pt_cells$pseudotime,
        method = "spearman"
      )
    )
  }
) |>
  purrr::list_rbind()
pt_sample_summary <- pt_cells |>
  dplyr::group_by(Mouse, Condition, Sample) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    median_pseudotime = median(pseudotime),
    q10_pseudotime = quantile(pseudotime, 0.1),
    q90_pseudotime = quantile(pseudotime, 0.9),
    .groups = "drop"
  )
pt_root_sensitivity <- purrr::map(pt_anchor_members[["MG"]], function(root) {
  tibble::tibble(
    root_cell = pt_cells$cell[root],
    spearman_rho = stats::cor(
      pt_cells$dpt_raw,
      pt_dpt[, root],
      method = "spearman"
    )
  )
}) |>
  purrr::list_rbind()

# Vary graph k and reselect the MG neighborhood using the same biological rule.
pt_k_sensitivity <- purrr::map(c(30L, 90L), function(k) {
  neighbors <- cbind(
    seq_len(nrow(pt_pca)),
    BiocNeighbors::findKNN(pt_pca, k = k)$index
  )
  contrast <- pt_scores_z[, "muller_glia"] - pt_scores_z[, "neuron_union"]
  nhood_contrast <- rowMeans(matrix(
    contrast[neighbors],
    nrow = nrow(neighbors)
  ))
  ids <- neighbors[which.max(nhood_contrast), ]
  root <- ids[which.min(rowSums(as.matrix(stats::dist(pt_pca[
    ids,
    ,
    drop = FALSE
  ]))))]
  set.seed(pt_seed)
  dm <- destiny::DiffusionMap(
    pt_pca,
    k = k,
    n_eigs = pt_n_eigs,
    sigma = "local",
    density_norm = TRUE,
    n_pcs = NA
  )
  dpt <- destiny::DPT(dm, tips = root)
  values <- as.numeric(dpt[, root])
  stopifnot(all(is.finite(values)), values[root] == 0)
  tibble::tibble(
    k,
    root_cell = pt_cells$cell[root],
    spearman_rho = stats::cor(pt_cells$dpt_raw, values, method = "spearman")
  )
}) |>
  purrr::list_rbind()

pt_landmarks <- c("Rlbp1", "Hes1", "Hes6", "Ascl1", "Otx2")
stopifnot(all(pt_landmarks %in% rownames(pt_expression)))

# ---- outputs: exploratory only; existing files require explicit opt-in ----

saveRDS(
  list(
    input_path = pt_input_path,
    settings = list(
      seed = pt_seed,
      dimensions = pt_dimensions,
      k = pt_k,
      n_eigs = pt_n_eigs,
      root_criterion = "mean(z(muller_glia) - z(neuron_union))"
    ),
    signatures = pt_signatures,
    marker_coverage = pt_marker_coverage,
    score_centers = attr(pt_scores_z, "scaled:center"),
    score_scales = attr(pt_scores_z, "scaled:scale"),
    cells = pt_cells,
    neighborhoods = pt_nhoods,
    neighborhood_members = pt_neighbors,
    neighborhood_sample_counts = pt_sample_counts,
    anchors = pt_anchor_summary,
    anchor_members = pt_anchor_members,
    anchor_samples = pt_anchor_samples,
    root_index = pt_root,
    diffusion_map = pt_dm,
    dpt = pt_dpt,
    root_sensitivity = pt_root_sensitivity,
    k_sensitivity = pt_k_sensitivity,
    score_associations = pt_score_associations,
    sample_summary = pt_sample_summary,
    session_info = utils::sessionInfo()
  ),
  output_path(pt_output_dir, "diffusion_pseudotime.rds")
)
readr::write_csv(pt_cells, output_path(pt_output_dir, "cell_pseudotime.csv"))
readr::write_csv(
  pt_nhoods,
  output_path(pt_output_dir, "neighborhood_scores.csv")
)
readr::write_csv(
  pt_anchor_summary,
  output_path(pt_output_dir, "anchor_neighborhoods.csv")
)
readr::write_csv(
  pt_root_sensitivity,
  output_path(pt_output_dir, "root_sensitivity.csv")
)
readr::write_csv(
  pt_k_sensitivity,
  output_path(pt_output_dir, "k_sensitivity.csv")
)
readr::write_csv(
  pt_score_associations,
  output_path(pt_output_dir, "score_associations.csv")
)
readr::write_csv(
  pt_sample_summary,
  output_path(pt_output_dir, "sample_summary.csv")
)
ggsave(
  output_path(pt_output_dir, "diffusion_pseudotime.png"),
  plot = ggplot(pt_cells, aes(DC1, DC2, color = pseudotime)) +
    geom_point(size = 0.7) +
    scale_color_viridis_c() +
    labs(color = "Pseudotime"),
  width = 7,
  height = 5,
  dpi = 200
)

pt_umap_reduction <- paste0(
  "umap_",
  pt_config$selected$mg$branch,
  "_dims",
  pt_dimensions
)
stopifnot(pt_umap_reduction %in% SeuratObject::Reductions(pt_sobj))
pt_umap_coordinates <- SeuratObject::Embeddings(pt_sobj, pt_umap_reduction)
stopifnot(identical(rownames(pt_umap_coordinates), pt_cells$cell))
pt_umap_cells <- tibble::tibble(
  cell = pt_cells$cell,
  UMAP1 = pt_umap_coordinates[, 1],
  UMAP2 = pt_umap_coordinates[, 2],
  pseudotime = pt_cells$pseudotime
)
ggsave(
  output_path(pt_output_dir, "umap_diffusion_pseudotime.png"),
  plot = ggplot(pt_umap_cells, aes(UMAP1, UMAP2, color = pseudotime)) +
    geom_point(size = 0.8) +
    scale_color_viridis_c(limits = c(0, 1)) +
    coord_equal() +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Pseudotime"),
  width = 8,
  height = 6,
  dpi = 300
)

pt_marker_counts <- SeuratObject::LayerData(
  pt_sobj,
  assay = "RNA",
  layer = "counts"
)
pt_marker_library_sizes <- Matrix::colSums(pt_marker_counts)
stopifnot(
  identical(colnames(pt_marker_counts), pt_cells$cell),
  all(pt_landmarks %in% rownames(pt_marker_counts)),
  all(is.finite(pt_marker_library_sizes)),
  all(pt_marker_library_sizes > 0)
)
pt_plot_normalized <- sweep(
  as.matrix(pt_marker_counts[pt_landmarks, , drop = FALSE]),
  2L,
  pt_marker_library_sizes,
  "/"
) *
  10000
pt_marker_cells <- as.data.frame(t(pt_plot_normalized)) |>
  tibble::rownames_to_column("cell") |>
  dplyr::left_join(
    dplyr::select(pt_cells, cell, pseudotime_percentile),
    by = "cell"
  ) |>
  tidyr::pivot_longer(
    dplyr::all_of(pt_landmarks),
    names_to = "gene",
    values_to = "normalized_counts"
  ) |>
  dplyr::mutate(
    gene = factor(gene, levels = pt_landmarks),
    detected = as.numeric(normalized_counts > 0)
  )
pt_bar_width <- 1 / (nrow(pt_cells) - 1)
stopifnot(all(
  abs(diff(sort(unique(pt_cells$pseudotime_percentile))) - pt_bar_width) < 1e-12
))
pt_marker_cols <- pt_marker_cells |>
  dplyr::filter(pseudotime_percentile >= 0.75) |>
  dplyr::arrange(gene, pseudotime_percentile)
ggsave(
  output_path(pt_output_dir, "marker_expression_top25_cell_bars_viridis.png"),
  plot = ggplot(pt_marker_cols, aes(pseudotime_percentile, normalized_counts)) +
    geom_col(
      width = pt_bar_width,
      position = "identity",
      colour = NA,
      fill = viridisLite::viridis(2)[2]
    ) +
    facet_wrap(~gene, scales = "free_y", ncol = 1) +
    scale_x_continuous(
      breaks = seq(0.75, 1, by = 0.05),
      labels = scales::label_percent(accuracy = 1)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = 0)) +
    coord_cartesian(
      xlim = c(
        min(pt_marker_cols$pseudotime_percentile) - pt_bar_width / 2,
        1 + pt_bar_width / 2
      ),
      expand = FALSE
    ) +
    theme(
      panel.background = element_rect(
        fill = viridisLite::viridis(2)[1],
        colour = NA
      ),
      axis.line = element_blank()
    ) +
    labs(x = "Pseudotime percentile", y = "Counts per 10,000"),
  width = 8,
  height = 10,
  dpi = 300
)
