# Analysis Manifest

ESPI's executable analysis pipeline currently lives in `scripts/` and `notebook/`, not under this Mycelium `analysis/` directory. Keep existing pipeline files in place unless the user explicitly asks to reorganize.

| Entry | Location | Type | Status | Notes |
|-------|----------|------|--------|-------|
| Preprocess Seurat object branches | `scripts/03-preprocess.R`, `scripts/03-preprocess-all.R` | R pipeline | active | Produces normalization and cell-cycle-filter branches plus QC/HVG/PCA diagnostics from either the legacy source object or the counts-derived, QC-filtered object selected at preprocessing. |
| Candidate clustering | `scripts/04-cluster.R`, `scripts/04-cluster-all.R`, `scripts/05-summarize-clusters.R` | R pipeline | active | Produces UMAP/clustree candidate clustering outputs, a 36-row supplemental grid summary table, a 12-panel clustree grid, and a representative UMAP resolution sweep from preprocessed objects. |
| Clustering criteria ideation | `analysis/ideas/2026-07-03-clustering-criteria-brainstorm/` | Mycelium ideation session | active | Persona-generated criteria ideas for label-blind selection of normalization, PC count, and Leiden resolution. |
| Cluster proportion testing ideation | `analysis/ideas/2026-07-05-cluster-proportion-testing/` | Mycelium ideation session | active | Targeted statistics, causal inference, compositional, single-cell methods, and experimental design ideas for testing E-Stim-associated MG-selected cluster proportion differences using Mouse × Condition samples rather than cell-pooled Fisher tests. |
| Single-cell analysis notebook | `notebook/sc_analysis.qmd` | Quarto notebook | active | Uses notebook-relative figure paths; rerender after source figure updates. |
| Cell type marker references | `data-raw/cell-type-marker-genes.R`, `data/cell_type_marker_genes.rda`, `data/cell_type_marker_labels.rda` | R package data | active | Curated broad retinal cell type marker lists and display labels for downstream single-cell annotation. |
| Cell type marker heatmap | `scripts/06-plot-marker-heatmap.R` | R figure script | active | Generates per-cell marker heatmap PNG/PDF outputs from any clustered Seurat object (branch tag derived from `sobj@misc$clustering$branch_tag`) and symlinks the PNG into `notebook/figures/`. |
| Cluster marker module p27 heatmaps | `scripts/10-plot-cluster-marker-heatmaps.R`, `R/cluster-marker-heatmap.R` | R figure script + helper | active | Generates per-cluster cell-type module score heatmaps with sample-aware within-Mouse × Condition p27 enrichment z-score strips, plus PNG/PDF, module-score TSV, p27-enrichment TSV, and notebook symlinks. |
| MG-selected analysis branch | `scripts/07-select-mg-subset.R`, `scripts/08-summarize-mg-clusters.R`, `scripts/09-plot-mg-figures.R`, `scripts/10-plot-cluster-marker-heatmaps.R`, `scripts/12-run-mg-de.R`, `scripts/11-find-mg-markers.R`, `R/cluster-abundance.R`, `R/cluster-marker-heatmap.R`, `notebook/sc_analysis.qmd` | R pipeline + notebook | active | Scores cluster marker modules, removes configured confident microglia/photoreceptor and p27-high source clusters, regenerates PFlog MG-selected clustering/figures, plots MG-selected UMAPs and descriptive pooled cell-level cluster abundance Fisher/CLR summaries, plots per-cluster module/p27 heatmaps, runs Mouse × Condition pseudobulk DE/DD/enrichment, and runs descriptive `FindAllMarkers()` cluster-marker ranking. |
| MG-selected manuscript write-up plan | `analysis/MG_SELECTED_WRITEUP_PLAN.md` | Markdown planning note | draft | Summarizes current MG-selected results, pipeline weaknesses, and planned edits for manuscript-facing prose in `notebook/sc_analysis.qmd`. |
| Tripwire checks | `tools/run-tripwires.R`, `analysis_labels.yml` | R QA runner | active | Checks cluster wrapper execution, branch artifact separation, report freshness, missing-input failure, metadata contract, label firewall, and future contrast direction. |

## For Future Mycelium Analyses

Create a subdirectory under `analysis/` only for new standalone analysis reports. Each subdirectory should include an UPPER_SNAKE_CASE `.md` note describing inputs, commands, outputs, and validation.
