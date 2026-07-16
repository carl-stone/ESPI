# Analysis Manifest

ESPI's active executable analysis consists of four phase scripts and the
Quarto notebook. Shared code lives in four focused package modules.

| Entry | Location | Type | Status | Notes |
|-------|----------|------|--------|-------|
| Routine publication interface | `justfile` | Just recipes | active/current | Five analysis commands are `just run [overwrite]`, `just figures [overwrite]`, `just markers [overwrite]`, `just de [overwrite]`, and deliberate `just regenerate-frozen`. Maintenance recipes are `load`, `document`, `readme`, `format`, and `lint`. |
| Frozen regeneration | `scripts/01-regenerate-frozen.R` | R maintenance script | active/current | Rebuilds counts, QC, four preprocessing branches, source and MG grids, summaries, and frozen-stage artifacts. It has no arguments, requires writable frozen-object directories, and refuses read-only directories. |
| Publication figures | `scripts/02-publication-figures.R` | R analysis script | active/current | Loads the final source, final MG-selected, and CC-filtered MG sensitivity objects once each; writes source/MG descriptive figures, heatmaps, UMAPs, abundance summaries, and required supplemental artifacts. |
| Marker analysis | `scripts/03-marker-analysis.R` | R analysis script | active/current | Runs fixed no-merge `FindAllMarkers()` analysis on the final MG object and writes the four marker tables plus dotplot. Its outputs do not feed phase 04. |
| DE and enrichment | `scripts/04-de-enrichment.R` | R analysis script | active/current | Independently loads the final MG object, builds six Mouse × Condition pseudobulk samples, runs primary and paired DE, rebuilds curated marker overlap, and writes GO/Bayesian enrichment artifacts. |
| Configuration and contracts | `R/config.R` | R package module | active/current | Owns paths, labels, seed, palettes, chosen object contracts, and fixed output/overwrite invariants. |
| Seurat methods | `R/seurat-methods.R` | R package module | active/current | Owns PFlog/log1p PCA and nonstandard cluster-grid summary/stability calculations. |
| Publication analysis | `R/publication-analysis.R` | R package module | active/current | Owns cluster abundance, sample proportions, exact randomization, module scores, and p27 enrichment computations. |
| Publication plots | `R/publication-plots.R` | R package module | active/current | Owns the publication theme, safe figure writer, curated marker heatmap, and module/p27 heatmap writers. |
| Single-cell analysis notebook | `notebook/sc_analysis.qmd` | Quarto notebook | active/current | Consumer only. Visible prose, captions, order, and values remain locked; rendering uses notebook-relative regular-file figures. |
| Cell type marker references | `data-raw/cell-type-marker-genes.R`, `data/cell_type_marker_genes.rda`, `data/cell_type_marker_labels.rda` | R package data | active | Curated broad retinal cell-type marker lists and display labels for annotation and curated marker overlap. |
| Clustering criteria ideation | `analysis/ideas/2026-07-03-clustering-criteria-brainstorm/` | Mycelium ideation session | active | Persona-generated criteria ideas for label-blind selection of normalization, PC count, and Leiden resolution. |
| Cluster proportion testing ideation | `analysis/ideas/2026-07-05-cluster-proportion-testing/` | Mycelium ideation session | active | Methods and design ideas for Mouse × Condition cluster proportion comparisons rather than cell-pooled inference. |
| MG-selected manuscript write-up plan | `analysis/MG_SELECTED_WRITEUP_PLAN.md` | Markdown planning note | implemented/current | Records current MG-selected results, curated DE effects, primary volcano specification, enrichment themes, interpretation limits, and notebook endpoint. |

Notebook mirrors reject symlink destinations, copy through a temporary regular
sibling, verify hashes and dimensions, and atomically replace regular files.

## For Future Mycelium Analyses

Create a subdirectory under `analysis/` only for new standalone analysis reports. Each subdirectory should include an UPPER_SNAKE_CASE `.md` note describing inputs, commands, outputs, and validation.
