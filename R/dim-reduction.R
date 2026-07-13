#' Normalize with log1p and compute PCA for one Seurat object.
#'
#' Writes results to `sobj[["pca"]]`.
#' Records provenance under `sobj@misc$preprocessing` matching the fields listed
#' in the ESPI preprocessing contract.
#'
#' @param sobj Seurat object with `counts` layer populated and `VariableFeatures`
#'   already selected.
#' @param n_pcs Integer number of PCs to compute. Default 50.
#'
#' @return `sobj` with `pca` reduction and `misc$preprocessing` populated.
#' @export
run_log1p_pca <- function(sobj, n_pcs = 50) {
  stopifnot(length(n_pcs) == 1, is.numeric(n_pcs), is.finite(n_pcs), n_pcs > 0)

  features <- SeuratObject::VariableFeatures(sobj)
  if (length(features) == 0) {
    stop("VariableFeatures(sobj) is empty.", call. = FALSE)
  }
  filtered_cell_cycle <- isTRUE(sobj@misc$preprocessing$filtered_cell_cycle)

  sobj <- Seurat::NormalizeData(sobj)
  sobj <- Seurat::ScaleData(sobj, features = features)
  sobj <- Seurat::RunPCA(sobj, features = features, npcs = n_pcs)

  sobj@misc$active.reduction <- "pca"
  sobj@misc$preprocessing <- list(
    normalization = "log1p",
    pca_method = "Seurat::RunPCA",
    pca_source_layer = "scale.data",
    hvg_method = "Seurat::FindVariableFeatures(selection.method = 'vst')",
    n_variable_features = length(features),
    n_pcs = n_pcs,
    filtered_cell_cycle = filtered_cell_cycle
  )

  sobj
}

#' Normalize with PFlog (scclrR) and compute PCA on variable features.
#'
#' Writes results to `sobj[["pca"]]`.
#' PFlog normalizes over the full retained feature set so that its per-cell
#' center is well defined; PCA is then computed on `VariableFeatures(sobj)`
#' via `scclrR::pca_matrix` using the full PFlog center.
#'
#' @param sobj Seurat object with `counts` layer populated and `VariableFeatures`
#'   already selected.
#' @param n_pcs Integer number of PCs to compute. Default 50.
#'
#' @return `sobj` with `pca` reduction and `misc$preprocessing` populated.
#' @export
run_pflog_pca <- function(sobj, n_pcs = 50) {
  stopifnot(length(n_pcs) == 1, is.numeric(n_pcs), is.finite(n_pcs), n_pcs > 0)

  assay <- SeuratObject::DefaultAssay(sobj)
  if (!inherits(sobj[[assay]], "Assay5")) {
    sobj[[assay]] <- as(sobj[[assay]], Class = "Assay5")
  }
  features <- SeuratObject::VariableFeatures(sobj)
  if (length(features) == 0) {
    stop("VariableFeatures(sobj) is empty.", call. = FALSE)
  }
  filtered_cell_cycle <- isTRUE(sobj@misc$preprocessing$filtered_cell_cycle)
  sobj <- Seurat::NormalizeData(sobj)

  sobj <- scclrR::pflog(sobj)
  pflog_layer <- SeuratObject::LayerData(sobj[[assay]], layer = "pflog")
  center <- sobj[["pflog_center", drop = TRUE]]
  pca <- scclrR::pca_matrix(
    sparse = pflog_layer[features, , drop = FALSE],
    center = center,
    n.components = n_pcs,
    seed = SEED
  )

  pc_names <- paste0("PC_", seq_len(ncol(pca$scores)))
  colnames(pca$scores) <- pc_names
  colnames(pca$loadings) <- pc_names
  if (is.null(rownames(pca$scores))) {
    rownames(pca$scores) <- colnames(sobj)
  }
  if (is.null(rownames(pca$loadings))) {
    rownames(pca$loadings) <- features
  }

  sobj[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = pca$scores,
    loadings = pca$loadings,
    stdev = sqrt(pca$explained_variance),
    key = "PC_",
    assay = assay
  )

  sobj@misc$active.reduction <- "pca"
  sobj@misc$preprocessing <- list(
    normalization = "pflog",
    pca_method = "scclrR::pca_matrix",
    pca_source_layer = "pflog",
    pca_center_key = "pflog_center",
    hvg_method = "Seurat::FindVariableFeatures(selection.method = 'vst')",
    n_variable_features = length(features),
    n_pcs = n_pcs,
    filtered_cell_cycle = filtered_cell_cycle
  )

  sobj
}
