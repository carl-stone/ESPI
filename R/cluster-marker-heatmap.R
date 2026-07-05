#' Compute cluster-level cell-type marker module scores.
#'
#' Scores each cell for each named marker-gene set with
#' [Seurat::AddModuleScore()], then averages the module scores within clusters.
#'
#' @param sobj Seurat object containing expression and cluster metadata.
#' @param cluster_col Metadata column containing cluster labels.
#' @param marker_genes Named non-empty list of marker-gene character vectors.
#' @param assay Assay used for module scoring.
#' @param slot Assay slot/layer used by [Seurat::AddModuleScore()].
#' @param seed Random seed forwarded to [Seurat::AddModuleScore()].
#'
#' @return Numeric matrix with marker-gene sets in rows, sorted clusters in
#'   columns, and mean module scores as values.
#' @export
compute_cluster_module_scores <- function(
  sobj,
  cluster_col,
  marker_genes,
  assay = "RNA",
  slot = "data",
  seed = SEED
) {
  meta <- sobj@meta.data
  if (!cluster_col %in% colnames(meta)) {
    stop("Missing cluster metadata column: ", cluster_col, call. = FALSE)
  }
  if (!assay %in% SeuratObject::Assays(sobj)) {
    stop("Missing assay: ", assay, call. = FALSE)
  }
  if (!slot %in% SeuratObject::Layers(sobj[[assay]])) {
    stop(
      "Missing score slot/layer '",
      slot,
      "' in assay ",
      assay,
      ". Available layers: ",
      paste(SeuratObject::Layers(sobj[[assay]]), collapse = ", "),
      call. = FALSE
    )
  }
  if (
    !is.list(marker_genes) ||
      length(marker_genes) == 0L ||
      is.null(names(marker_genes)) ||
      anyNA(names(marker_genes)) ||
      any(!nzchar(names(marker_genes)))
  ) {
    stop("marker_genes must be a named non-empty list.", call. = FALSE)
  }
  if (any(duplicated(names(marker_genes)))) {
    stop("marker_genes names must be unique.", call. = FALSE)
  }
  invalid_sets <- names(marker_genes)[
    !vapply(
      marker_genes,
      function(x) {
        is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x))
      },
      logical(1L)
    )
  ]
  if (length(invalid_sets) > 0L) {
    stop(
      "Marker gene set(s) must be non-empty character vectors: ",
      paste(invalid_sets, collapse = ", "),
      call. = FALSE
    )
  }

  cluster <- as.character(meta[[cluster_col]])
  if (length(cluster) != nrow(meta)) {
    stop("Cluster metadata is not aligned to cells.", call. = FALSE)
  }
  if (anyNA(cluster) || any(!nzchar(cluster))) {
    stop(
      "Cluster metadata contains missing or empty labels: ",
      cluster_col,
      call. = FALSE
    )
  }

  marker_vector <- unlist(marker_genes, use.names = FALSE)
  missing_genes <- setdiff(unique(marker_vector), rownames(sobj))
  if (length(missing_genes) > 0L) {
    stop(
      "Marker gene(s) missing from the Seurat object: ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  cluster_levels <- .sort_cluster_labels(cluster)
  module_prefix <- "celltype_module_score"
  module_score_cols <- paste0(module_prefix, seq_along(marker_genes))
  existing_module_cols <- intersect(module_score_cols, colnames(sobj@meta.data))
  if (length(existing_module_cols) > 0L) {
    sobj@meta.data[existing_module_cols] <- NULL
  }
  sobj <- Seurat::AddModuleScore(
    object = sobj,
    features = marker_genes,
    assay = assay,
    name = module_prefix,
    seed = seed,
    search = FALSE,
    slot = slot
  )
  missing_module_cols <- setdiff(module_score_cols, colnames(sobj@meta.data))
  if (length(missing_module_cols) > 0L) {
    stop(
      "Seurat::AddModuleScore did not create expected column(s): ",
      paste(missing_module_cols, collapse = ", "),
      call. = FALSE
    )
  }

  module_scores <- sobj@meta.data[module_score_cols]
  colnames(module_scores) <- names(marker_genes)
  module_scores$cluster <- cluster
  cluster_marker_scores <- stats::aggregate(
    module_scores[names(marker_genes)],
    by = list(cluster = module_scores$cluster),
    FUN = mean
  )
  cluster_marker_scores <- cluster_marker_scores[
    match(cluster_levels, cluster_marker_scores$cluster),
    ,
    drop = FALSE
  ]
  score_matrix <- t(as.matrix(cluster_marker_scores[names(marker_genes)]))
  storage.mode(score_matrix) <- "numeric"
  rownames(score_matrix) <- names(marker_genes)
  colnames(score_matrix) <- cluster_levels
  score_matrix
}

#' Compute sample-aware cluster p27 enrichment z-scores.
#'
#' Computes observed per-cluster p27 expression, then estimates a null
#' distribution by permuting cluster labels within each Mouse x Condition
#' sample. The permutation preserves each sample's p27 values and cluster-size
#' composition. A cluster receives `NA` when its permutation null has zero or
#' non-finite standard deviation, which can occur when p27 is constant within
#' every sample.
#'
#' @param sobj Seurat object containing expression and metadata.
#' @param cluster_col Metadata column containing cluster labels.
#' @param gene Gene used for p27 expression. Defaults to `Cdkn1b`.
#' @param layer Assay layer used for expression values.
#' @param assay Assay containing `layer`.
#' @param mouse_col Metadata column containing mouse labels.
#' @param condition_col Metadata column containing condition labels.
#' @param n_perm Positive integer number of within-sample permutations.
#' @param seed Random seed used for permutations.
#'
#' @return Data frame with one row per sorted cluster and columns `cluster`,
#'   `n_cells`, `observed_mean`, `null_mean`, `null_sd`, and `z_score`.
#' @export
compute_cluster_p27_enrichment <- function(
  sobj,
  cluster_col,
  gene = "Cdkn1b",
  layer = "pflog",
  assay = "RNA",
  mouse_col = "Mouse",
  condition_col = CONDITION_COL,
  n_perm = 2000L,
  seed = SEED
) {
  meta <- sobj@meta.data
  required_cols <- c(cluster_col, mouse_col, condition_col)
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing required metadata column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (!assay %in% SeuratObject::Assays(sobj)) {
    stop("Missing assay: ", assay, call. = FALSE)
  }
  if (!gene %in% rownames(sobj)) {
    stop("Gene is missing from the Seurat object: ", gene, call. = FALSE)
  }
  available_layers <- SeuratObject::Layers(sobj[[assay]])
  if (!layer %in% available_layers) {
    stop(
      "Missing layer '",
      layer,
      "' in assay ",
      assay,
      ". Available layers: ",
      paste(available_layers, collapse = ", "),
      call. = FALSE
    )
  }
  if (
    length(n_perm) != 1L ||
      is.na(n_perm) ||
      n_perm < 1L ||
      n_perm != as.integer(n_perm)
  ) {
    stop("n_perm must be a positive integer.", call. = FALSE)
  }
  n_perm <- as.integer(n_perm)

  cluster <- as.character(meta[[cluster_col]])
  mouse <- as.character(meta[[mouse_col]])
  condition <- as.character(meta[[condition_col]])
  if (
    length(cluster) != nrow(meta) ||
      length(mouse) != nrow(meta) ||
      length(condition) != nrow(meta)
  ) {
    stop(
      "Cluster, mouse, and condition metadata are not aligned to cells.",
      call. = FALSE
    )
  }
  if (anyNA(cluster) || any(!nzchar(cluster))) {
    stop(
      "Cluster metadata contains missing or empty labels: ",
      cluster_col,
      call. = FALSE
    )
  }
  if (anyNA(mouse) || any(!nzchar(mouse))) {
    stop(
      "Mouse metadata contains missing or empty labels: ",
      mouse_col,
      call. = FALSE
    )
  }
  if (anyNA(condition) || any(!nzchar(condition))) {
    stop(
      "Condition metadata contains missing or empty labels: ",
      condition_col,
      call. = FALSE
    )
  }

  expression_matrix <- SeuratObject::LayerData(sobj[[assay]], layer = layer)
  expr <- as.numeric(expression_matrix[gene, colnames(sobj), drop = TRUE])
  if (length(expr) != nrow(meta)) {
    stop("Expression values are not aligned to cells.", call. = FALSE)
  }
  if (anyNA(expr) || any(!is.finite(expr))) {
    stop(
      "Expression values contain missing or non-finite values.",
      call. = FALSE
    )
  }

  cluster_levels <- .sort_cluster_labels(cluster)
  cluster_factor <- factor(cluster, levels = cluster_levels)
  observed <- as.numeric(tapply(expr, cluster_factor, mean))
  n_cells <- as.integer(table(cluster_factor))

  sample_id <- paste(mouse, condition, sep = "__")
  sample_indices <- split(seq_along(sample_id), sample_id)
  null_mat <- matrix(
    NA_real_,
    nrow = n_perm,
    ncol = length(cluster_levels),
    dimnames = list(NULL, cluster_levels)
  )

  had_random_seed <- exists(
    ".Random.seed",
    envir = globalenv(),
    inherits = FALSE
  )
  if (had_random_seed) {
    saved_random_seed <- get(
      ".Random.seed",
      envir = globalenv(),
      inherits = FALSE
    )
  }
  on.exit(
    {
      if (had_random_seed) {
        assign(".Random.seed", saved_random_seed, envir = globalenv())
      } else if (
        exists(".Random.seed", envir = globalenv(), inherits = FALSE)
      ) {
        rm(".Random.seed", envir = globalenv())
      }
    },
    add = TRUE
  )

  set.seed(seed)
  for (perm_idx in seq_len(n_perm)) {
    permuted_cluster <- cluster
    for (idx in sample_indices) {
      permuted_cluster[idx] <- sample(
        cluster[idx],
        length(idx),
        replace = FALSE
      )
    }
    null_mat[perm_idx, ] <- as.numeric(tapply(
      expr,
      factor(permuted_cluster, levels = cluster_levels),
      mean
    ))
  }

  null_mean <- colMeans(null_mat)
  null_sd <- apply(null_mat, 2L, stats::sd)
  z_score <- (observed - null_mean) / null_sd
  degenerate_null <- !is.finite(null_sd) | null_sd == 0
  z_score[degenerate_null] <- NA_real_

  data.frame(
    cluster = cluster_levels,
    n_cells = n_cells,
    observed_mean = observed,
    null_mean = null_mean,
    null_sd = null_sd,
    z_score = z_score,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.sort_cluster_labels <- function(x) {
  x <- unique(as.character(x))
  if (all(grepl("^-?[0-9]+$", x))) {
    return(as.character(sort(as.integer(x))))
  }
  sort(x, method = "radix")
}
