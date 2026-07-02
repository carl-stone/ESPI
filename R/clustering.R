#' Run FindNeighbors + Leiden FindClusters and store the labels under a name.
#'
#' Uses the `pca` reduction. Writes the resulting cluster labels into
#' `sobj@meta.data[[name]]` without overwriting existing identities. The
#' active identity is left unchanged.
#'
#' @param sobj Seurat object with `pca` reduction populated.
#' @param dims Integer vector of PCs to use, e.g. `1:30`.
#' @param resolution Numeric Leiden resolution.
#' @param name Character column name to store cluster labels.
#'
#' @return `sobj` with `meta.data[[name]]` populated.
#' @export
find_leiden_clusters <- function(sobj, dims, resolution, name) {
  # 1. sobj <- Seurat::FindNeighbors(sobj, reduction = "pca", dims = dims)
  # 2. sobj <- Seurat::FindClusters(
  #      sobj,
  #      algorithm  = 4,                 # Leiden; requires the Leiden runtime
  #      resolution = resolution,
  #      random.seed = SEED
  #    )
  # 3. sobj@meta.data[[name]] <- SeuratObject::Idents(sobj)
  # 4. Restore prior Idents so the active identity is unchanged.
  # 5. return(sobj)
  stop("not implemented")
}

#' Compute a UMAP embedding for a specific PC-dims choice on `pca`.
#'
#' @param sobj Seurat object with `pca` reduction populated.
#' @param dims Integer vector of PCs to use, e.g. `1:30`.
#' @param reduction_name Character name for the stored reduction.
#'
#' @return `sobj` with `reductions[[reduction_name]]` populated.
#' @export
run_umap_for_dims <- function(sobj, dims, reduction_name) {
  # 1. sobj <- Seurat::RunUMAP(
  #      sobj,
  #      reduction        = "pca",
  #      dims             = dims,
  #      reduction.name   = reduction_name,
  #      reduction.key    = paste0(reduction_name, "_"),
  #      seed.use         = SEED
  #    )
  # 2. return(sobj)
  stop("not implemented")
}
