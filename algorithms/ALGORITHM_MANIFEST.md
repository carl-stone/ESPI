# Algorithm Manifest

Reusable ESPI methods live in the four focused R modules. The phase scripts
keep standard Seurat and plotting calls visible; this manifest maps the
nonstandard shared operations.

| Entry | Location | Type | Status | Notes |
|-------|----------|------|--------|-------|
| Fixed configuration and invariants | `R/config.R` | R package code | active/current | Owns paths, labels, seed, palettes, frozen-object contracts, output guards, and package-data documentation. |
| PCA branch methods | `R/seurat-methods.R` | R package code | active/current | Runs log1p and PFlog PCA with retained HVGs. |
| Cluster-grid summaries | `R/seurat-methods.R` | R package code | active/current | Writes the nonstandard grid summary and pairwise-stability calculations used by frozen regeneration. |
| Publication statistics | `R/publication-analysis.R` | R package code | active/current | Computes cluster abundance, sample cluster proportions, exact randomization tests, module scores, and p27 enrichment. |
| Publication plot writers | `R/publication-plots.R` | R package code | active/current | Provides the shared theme, PNG/PDF writer, curated marker heatmap, and module/p27 heatmap writers with safe notebook mirroring. |
| Frozen preprocessing and MG selection | `scripts/01-regenerate-frozen.R` | R analysis script | active/current | Applies fixed counts, QC, source preprocessing/clustering, and MG-selection methods. `ESPI_REGENERATION_START=mg-selection` resumes from the configured source RDS. |
| MG clustering sensitivity | `scripts/01b-cluster-mg-sensitivity.R` | R analysis script | active/current | Loads saved MG preprocessing objects and rebuilds the fixed Leiden/UMAP grid; the selected 20-PC/resolution-0.3 clustering and UMAP use seed 2847. |
| Publication phases | `scripts/02-publication-figures.R`, `scripts/03-marker-analysis.R`, `scripts/04-de-enrichment.R` | R analysis scripts | active/current | Keep standard Seurat, ggplot2, DESeq2, and enrichment calls visible in scientific order. |

## Rule

Do not create tiny helper functions for a few commands unless they are called
repeatedly or hide a real conceptual operation.
