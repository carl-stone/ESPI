# Analysis Manifest

ESPI's executable analysis pipeline currently lives in `scripts/` and `notebook/`, not under this Mycelium `analysis/` directory. Keep existing pipeline files in place unless the user explicitly asks to reorganize.

| Entry | Location | Type | Status | Notes |
|-------|----------|------|--------|-------|
| Preprocess Seurat object branches | `scripts/preprocess-sobj.R`, `scripts/preprocess-all.R` | R pipeline | active | Produces normalization and cell-cycle-filter branches plus QC/HVG/PCA diagnostics. |
| Candidate clustering | `scripts/cluster-sobj.R`, `scripts/cluster-all.R`, `scripts/summarize-cluster-grid.R` | R pipeline | active | Produces UMAP/clustree candidate clustering outputs, a 36-row supplemental grid summary table, a 12-panel clustree grid, and a representative UMAP resolution sweep from preprocessed objects. |
| Clustering criteria ideation | `analysis/ideas/2026-07-03-clustering-criteria-brainstorm/` | Mycelium ideation session | active | Persona-generated criteria ideas for label-blind selection of normalization, PC count, and Leiden resolution. |
| Single-cell analysis notebook | `notebook/sc_analysis.qmd` | Quarto notebook | active | Uses notebook-relative figure paths; rerender after source figure updates. |
| Cell type marker references | `data-raw/cell-type-marker-genes.R`, `data/cell_type_marker_genes.rda`, `data/cell_type_marker_labels.rda` | R package data | active | Curated broad retinal cell type marker lists and display labels for downstream single-cell annotation. |
| Tripwire checks | `tools/run-tripwires.R`, `analysis_labels.yml` | R QA runner | active | Checks cluster wrapper execution, branch artifact separation, report freshness, missing-input failure, metadata contract, label firewall, and future contrast direction. |

## For Future Mycelium Analyses

Create a subdirectory under `analysis/` only for new standalone analysis reports. Each subdirectory should include an UPPER_SNAKE_CASE `.md` note describing inputs, commands, outputs, and validation.
