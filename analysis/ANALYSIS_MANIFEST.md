# Analysis Manifest

ESPI's executable analysis pipeline currently lives in `scripts/` and `notebook/`, not under this Mycelium `analysis/` directory. Keep existing pipeline files in place unless the user explicitly asks to reorganize.

| Entry | Location | Type | Status | Notes |
|-------|----------|------|--------|-------|
| Preprocess Seurat object branches | `scripts/preprocess-sobj.R`, `scripts/preprocess-all.R` | R pipeline | active | Produces normalization and cell-cycle-filter branches plus QC/HVG/PCA diagnostics. |
| Candidate clustering | `scripts/cluster-sobj.R`, `scripts/cluster-all.R`, `scripts/summarize-cluster-grid.R` | R pipeline | active | Produces UMAP/clustree candidate clustering outputs, a 36-row supplemental grid summary table, a 12-panel clustree grid, and a representative UMAP resolution sweep from preprocessed objects. |
| Clustering criteria ideation | `analysis/ideas/2026-07-03-clustering-criteria-brainstorm/` | Mycelium ideation session | active | Persona-generated criteria ideas for label-blind selection of normalization, PC count, and Leiden resolution. |
| Single-cell analysis notebook | `notebook/sc_analysis.qmd` | Quarto notebook | active | Uses notebook-relative figure paths; rerender after source figure updates. |
| Cell type marker references | `data-raw/cell-type-marker-genes.R`, `data/cell_type_marker_genes.rda`, `data/cell_type_marker_labels.rda` | R package data | active | Curated broad retinal cell type marker lists and display labels for downstream single-cell annotation. |
| Cell type marker heatmap | `scripts/big-heatmap-plot.R` | R figure script | active | Generates per-cell marker heatmap PNG/PDF outputs from any clustered Seurat object (branch tag derived from `sobj@misc$clustering$branch_tag`) and symlinks the PNG into `notebook/figures/`. |
| MG-selected analysis branch | `scripts/filter-mg-subset.R`, `scripts/summarize-mg-selected-grid.R`, `scripts/plot-mg-selected-figures.R`, `scripts/run-mg-selected-de.R`, `scripts/find-markers-mg-selected.R`, `R/cluster-abundance.R`, `notebook/sc_analysis.qmd` | R pipeline + notebook | active | Scores cluster marker modules, removes configured confident microglia/photoreceptor and Cdkn1b-high source clusters, regenerates PFlog MG-selected clustering/figures, plots MG-selected UMAPs and descriptive pooled cell-level cluster abundance Fisher/CLR summaries, runs Mouse × Condition pseudobulk DE/DD/enrichment, and runs descriptive `FindAllMarkers()` cluster marker ranking. Current notebook includes no-cell-cycle and cell-cycle-filtered MG-selected UMAP/marker/feature figures, a six-bin marker dot plot, an abundance enrichment bar plot, and a DE-vs-DD effect scatter. |
| MG-selected manuscript write-up plan | `analysis/MG_SELECTED_WRITEUP_PLAN.md` | Markdown planning note | draft | Summarizes current MG-selected results, pipeline weaknesses, and planned edits for manuscript-facing prose in `notebook/sc_analysis.qmd`. |
| Tripwire checks | `tools/run-tripwires.R`, `analysis_labels.yml` | R QA runner | active | Checks cluster wrapper execution, branch artifact separation, report freshness, missing-input failure, metadata contract, label firewall, and future contrast direction. |

## For Future Mycelium Analyses

Create a subdirectory under `analysis/` only for new standalone analysis reports. Each subdirectory should include an UPPER_SNAKE_CASE `.md` note describing inputs, commands, outputs, and validation.
