# Analysis Manifest

The publication pipeline consists of four phases across five scripts and the
Quarto notebook. Shared code lives in four focused package modules. This is a
reference inventory, not a required bookkeeping step; code defines the current
workflow.

| Entry | Location | Type | Status | Notes |
|-------|----------|------|--------|-------|
| Routine publication interface | `justfile` | Just recipes | active/current | Five analysis commands are `just run [overwrite]`, `just figures [overwrite]`, `just markers [overwrite]`, `just de [overwrite]`, and deliberate `just regenerate-frozen [start]`, where `start` is `all` or `mg-selection`; regeneration continues through phases 02–04 and notebook rendering with overwrite enabled. Maintenance recipes are `load`, `document`, `readme`, `format`, and `lint`. |
| Frozen preprocessing and MG selection | `scripts/01-regenerate-frozen.R` | R maintenance script | active/current | Full mode rebuilds counts, QC, four source preprocessing branches, the source grid, MG selection, and both saved MG preprocessing branches. MG-selection mode loads the selected source RDS and resumes at MG selection. |
| MG clustering sensitivity | `scripts/01b-cluster-mg-sensitivity.R` | R maintenance script | active/current | Loads both saved MG preprocessing branches, rebuilds the existing Leiden/UMAP grid and summaries, and writes the selected frozen MG objects. The chosen no-CC 20-PC/resolution-0.3 result uses seed 2847. |
| Publication figures | `scripts/02-publication-figures.R` | R analysis script | active/current | Loads the final source, final MG-selected, and CC-filtered MG sensitivity objects once each; writes source/MG descriptive figures, heatmaps, UMAPs, abundance summaries, and required supplemental artifacts. |
| Marker analysis | `scripts/03-marker-analysis.R` | R analysis script | active/current | Runs fixed no-merge `FindAllMarkers()` analysis on the final MG object and writes the four marker tables plus dotplot. Its outputs do not feed phase 04. |
| DE and enrichment | `scripts/04-de-enrichment.R` | R analysis script | active/current | Independently loads the final MG object, builds six Mouse × Condition pseudobulk samples, runs primary and paired DE, rebuilds curated marker overlap, and writes GO/Bayesian enrichment artifacts. |
| Configuration and contracts | `R/config.R` | R package module | active/current | Owns paths, labels, seed, palettes, chosen object contracts, and fixed output/overwrite invariants. Frozen inputs are guarded by expected cell counts, cluster columns, and cluster counts rather than byte-level hashes. |
| Seurat methods | `R/seurat-methods.R` | R package module | active/current | Owns PFlog/log1p PCA and nonstandard cluster-grid summary/stability calculations. |
| Publication analysis | `R/publication-analysis.R` | R package module | active/current | Owns cluster abundance, sample proportions, exact randomization, module scores, and p27 enrichment computations. |
| Publication plots | `R/publication-plots.R` | R package module | active/current | Owns the publication theme, safe figure writer, curated marker heatmap, and module/p27 heatmap writers. |
| Single-cell analysis notebook | `notebook/sc_analysis.qmd` | Quarto notebook | active/current | Consumes saved analysis outputs and notebook-relative regular-file figures. Update prose and captions with the relevant evidence, then render the HTML deliverable. |
| Cell type marker references | `data-raw/cell-type-marker-genes.R`, `data/cell_type_marker_genes.rda`, `data/cell_type_marker_labels.rda` | R package data | active | Curated broad retinal cell-type marker lists and display labels for annotation and curated marker overlap. |
| MG-selected manuscript write-up plan | `analysis/MG_SELECTED_WRITEUP_PLAN.md` | Markdown planning note | implemented/current | Records current MG-selected results, curated DE effects, primary volcano specification, enrichment themes, interpretation limits, and notebook endpoint. |
| Sequential analysis write-up | `analysis/SEQUENTIAL_ANALYSIS_WRITEUP.md` | Markdown scientific prose | complete/current | Provides 17 standalone Methods and Results snippets in execution order, with hidden evidence comments and no added interpretation or conclusion. |
| Targeted neurogenic proportion | `analysis/targeted-neurogenic-proportion/` | Standalone sensitivity analysis | complete/exploratory | Tests six prespecified progenitor-high, proliferation-low gates with sample-level beta-binomial models, null bootstraps, optimizer diagnostics, and fixed tabular outputs. |
| Exploratory plotting sandbox | `analysis/exploratory/` | Interactive R workspace | exploratory/noncanonical | Uses the final MG object and current DE tables for last-mile gene-expression and detection-model exploration. It has no stable CLI or publication output contract. |

Notebook mirrors reject symlink destinations, copy through a temporary regular
sibling, verify hashes and dimensions, and atomically replace regular files.
