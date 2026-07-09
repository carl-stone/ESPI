# Algorithm Manifest

Reusable ESPI methods currently live in the R package source. Keep implementation in `R/`; use this manifest as a map, not as a relocation target.

| Entry | Location | Type | Status | Notes |
|-------|----------|------|--------|-------|
| QC filtering | `scripts/qc-filtering.R` | R analysis script | active | Retains cells with `nFeature_RNA >= 50`, `nCount_RNA >= 100`, and `percent.mt <= 20`; percent.ribo is diagnostic only. Writes QC artifacts and `INPUT_OBJECT_DIR/sobj_qc_filtered.rds`. |
| PCA branch helpers | `R/dim-reduction.R` | R package code | active | Runs log1p and PFlog PCA using retained HVGs. |
| Preprocessing diagnostic plots | `R/preprocess-plots.R` | R package code | active | Saves QC violin, HVG scatter, DimHeatmap, and elbow plots. |
| Clustering diagnostic plots | `R/cluster-plots.R` | R package code | active | Saves UMAP overlays and clustree plots. |
| Cluster grid summaries | `R/cluster-grid-summary.R` | R package code | active | Writes supplemental clustering grid summaries against the current PFlog filtered 50-PC resolution-0.3 reference, plus pairwise local-neighbor stability tables, clustree-style split stability tables, and multi-panel clustering diagnostics. |

## Rule

Do not create tiny helper functions for a few commands unless they are called repeatedly or hide a real conceptual operation.
