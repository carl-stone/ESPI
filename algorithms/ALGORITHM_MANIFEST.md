# Algorithm Manifest

Reusable ESPI methods currently live in the R package source. Keep implementation in `R/`; use this manifest as a map, not as a relocation target.

| Entry | Location | Type | Status | Notes |
|-------|----------|------|--------|-------|
| PCA branch helpers | `R/dim-reduction.R` | R package code | active | Runs log1p and PFlog PCA using retained HVGs. |
| Preprocessing diagnostic plots | `R/preprocess-plots.R` | R package code | active | Saves QC violin, HVG scatter, DimHeatmap, and elbow plots. |
| Clustering diagnostic plots | `R/cluster-plots.R` | R package code | active | Saves UMAP overlays and clustree plots. |

## Rule

Do not create tiny helper functions for a few commands unless they are called repeatedly or hide a real conceptual operation.
