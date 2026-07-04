# Decision Log

Append-only log of non-obvious decisions and their rationale.

### [2026-07-03] Enable Mycelium without restructuring ESPI

**Tags**: mycelium, repo-structure, r-package

**Context**: ESPI is already organized as a minimal R package with executable analysis scripts and a Quarto notebook.

**Decision**: Initialize Mycelium's `.living/` memory layer, manifests, convention packs, and local Claude hooks, but keep existing `R/`, `scripts/`, and `notebook/` layout intact.

**Alternatives considered**: Moving existing scripts into Mycelium `analysis/` folders would satisfy generic Mycelium structure but would obscure the current package workflow.

**Rationale**: The repo's scientific and executable source of truth already lives in stable R package locations. Mycelium should add memory and provenance, not force a new layout.

**Consequences**: Future agents should use `analysis/` for new standalone Mycelium reports only. Existing pipeline work remains in `scripts/` and `R/`.

### [2026-07-03] Install bioinformatics conventions by default

**Tags**: mycelium, conventions, bioinformatics, single-cell

**Context**: Mycelium core conventions were installed manually because the initializer did not auto-detect the plugin network path from the OMP plugin cache. ESPI is a single-cell bioinformatics analysis.

**Decision**: Install `robust-analysis`, `report-generator`, `idea-generator`, and the domain-specific `bioinformatics` convention pack.

**Alternatives considered**: Installing only core packs would be more conservative but would omit relevant single-cell and genomics guidance.

**Rationale**: The project context clearly fits the bioinformatics convention pack.

**Consequences**: Generic conventions apply only when compatible with ESPI-specific rules in `AGENTS.md` and `.living/conventions.md`.

### [2026-07-03] Treat Mycelium restructure audit as advisory only

**Tags**: mycelium, repo-structure, audit, r-package

**Context**: The Mycelium init protocol says existing repos should run restructure audit before moving files. The audit reported 13 reference-material candidates and 64 unclassified package/notebook/generated files.

**Decision**: Do not apply the restructure moves automatically.

**Alternatives considered**: Moving root docs and agent files into `reference_material/` would align with generic Mycelium output but would break or obscure existing ESPI/agent conventions.

**Rationale**: ESPI is a working R package; package files, agent files, Quarto notebook files, generated Rd files, and existing docs should stay in their current locations unless the user explicitly requests a layout migration.

**Consequences**: `reference_material/` remains a manifest/index area. Any future restructuring requires explicit user approval.

### [2026-07-03] Enable skill-bridge after cloning available skillpacks

**Tags**: mycelium, skillpacks, skill-bridge, conventions

**Context**: The user asked to set up the external skillpacks referenced by Mycelium. `scientific-agent-skills` and `bioSkills` cloned successfully; the documented `Autonomous-Science` URL was not found.

**Decision**: Install the `skill-bridge` convention pack now and document the missing `Autonomous-Science` repo instead of blocking use of the available skillpacks.

**Alternatives considered**: Waiting for all three external repositories would keep setup symmetric but would leave the two available references unused.

**Rationale**: `skill-bridge` can route to the installed skillpack libraries immediately, and the missing repo is a source URL problem outside this repo.

**Consequences**: `skillpacks/README.md` records exact installed commits and the missing URL. Add `Autonomous-Science` later only after a corrected public or accessible URL is available.

### [2026-07-03] Use Rscript orchestration and Seurat-safe cluster branch tags

**Tags**: clustering, seurat, scripts, reproducibility

**Context**: `scripts/cluster-all.R` used a shell loop around `CURRENT_OBJECT_DIR`, and Seurat rewrote hyphenated UMAP reduction names such as `umap_pflog_no-filter-cc_dims20` to dotted names during storage. `DimPlot()` then looked up the unsanitized name and failed.

**Decision**: Keep `cluster-all.R` as an Rscript orchestrator that loads ESPI constants in R, supports `--dry-run`, and calls `cluster-sobj.R` with `system2()`. Use underscore branch tags (`filter_cc`, `no_filter_cc`) for clustered Seurat reductions, cluster metadata columns, clustree tags, and clustered RDS names.

**Alternatives considered**: Exporting `CURRENT_OBJECT_DIR` to the shell would preserve the old wrapper shape but keep path state split across R and shell. Keeping hyphen tags for reductions would continue relying on names Seurat does not store literally.

**Rationale**: A single R orchestration layer removes shell environment coupling. Underscore branch tags preserve the normalization + cell-cycle distinction while matching Seurat's object-name constraints.

**Consequences**: Preprocessed input filenames can keep the existing `filter-cc` / `no-filter-cc` names. New clustered artifacts use names such as `cluster_pflog_no_filter_cc_elbow20.rds`.

### [2026-07-03] Summarize clustering sensitivity with tables plus representative figures

**Tags**: clustering, sensitivity, supplemental-figures, reporting

**Context**: The clustering grid spans normalization method, cell-cycle-HVG policy, PC count, and Leiden resolution, which creates 36 clustering configurations before marker overlays or downstream analyses.

**Decision**: Generate a supplemental 36-row grid summary table with ARI and best-overlap Jaccard against `cluster_pflog_filter_cc_dims50_res0.3`, show the full clustree grid as a 12-panel figure, and show a representative PFlog filtered 50-PC UMAP resolution sweep instead of displaying every possible UMAP panel in the notebook.

**Alternatives considered**: Showing every UMAP and marker plot would be exhaustive but unreadable. Reporting only a qualitative statement would be too hand-wavy for a parameter-sensitive clustering choice.

**Rationale**: The table makes the full parameter grid auditable, while the clustree and representative UMAP figures show the relevant structure without overwhelming the report.

**Consequences**: `scripts/summarize-cluster-grid.R` must be rerun after regenerating clustered objects, and `notebook/sc_analysis.qmd` must be rerendered after the supplemental figures change.

### [2026-07-03] Store cell type markers as package data

**Tags**: marker-genes, annotation, package-data, single-cell

**Context**: Broad retinal cell type annotation needs reusable marker lists and display labels that can feed Seurat and ComplexHeatmap-style plots.

**Decision**: Keep the editable marker source in `data-raw/cell-type-marker-genes.R` and expose `cell_type_marker_genes` plus `cell_type_marker_labels` as small R package data objects in `data/`.

**Alternatives considered**: Keeping markers only in notebook code would hide reusable annotation inputs; using display labels such as `Müller glia` as list keys would make later column names and file names more fragile.

**Rationale**: ASCII marker-list keys are stable for code, while labels preserve human-readable names for plots. Package data loads naturally through `devtools::load_all()` and mirrors the existing `mouse_cell_cycle_genes` pattern.

**Consequences**: Marker provenance is documented as Ed and Megan's domain/literature curation, but marker-by-marker rationale is not recorded in the repo.

### [2026-07-03] Keep manuscript marker heatmap as a script

**Tags**: marker-genes, heatmap, plotting, r-package

**Context**: The manuscript needs one per-cell marker heatmap using the curated marker gene lists and the selected PFlog filtered clustering.

**Decision**: Implement the per-cell heatmap as `scripts/big-heatmap-plot.R` instead of adding a new `R/` helper file, and default its expression layer to PFlog to match the selected clustering branch.

**Alternatives considered**: A package helper would make sense if marker heatmaps become reusable across analyses, but this figure is currently one standard manuscript plot.

**Rationale**: Keeping the code in a script preserves locality for figure-specific choices while avoiding a shallow package module. Using the PFlog layer keeps the marker visualization on the same normalization scale used for the selected PCA/UMAP/clustering branch.

**Consequences**: If the same heatmap workflow gets reused in multiple scripts or reports, promote the repeated logic into `R/`.

### [2026-07-03] Use top-to-bottom executable script style

**Tags**: scripts, r, interactivity, reproducibility

**Context**: The analysis scripts need to be callable from Bash while staying easy to step through interactively in RStudio.

**Decision**: Structure scripts as purpose/usage docs, package loading, parameter defaults and CLI overrides near the top, validation/work sections, and output side effects at the end. Avoid wrapping simple scripts in `main()`.

**Alternatives considered**: `main()` wrappers are cleaner for software packages, but they make interactive line-by-line execution worse for this analysis workflow.

**Rationale**: Top-to-bottom scripts preserve shell reproducibility without hiding inspectable intermediate objects from RStudio.

**Consequences**: Promote repeated logic to `R/` only when it hides a real conceptual operation or removes substantial repetition.

### [2026-07-04] Use all Mouse × Condition samples as primary DE unit

**Tags**: differential-expression, pseudobulk, paired-design, mg-selected

**Context**: The `mg-selected` dataset has six Mouse × Condition pseudobulk samples: paired mice 10 and 3, mouse 30 as E-Stim only, and mouse 33 as control only.

**Decision**: Run primary DESeq2 and differential-detection models with `~ condition` across all six samples, and save `~ mouse + condition` results as a paired sensitivity analysis on mice 10 and 3 only.

**Alternatives considered**: A paired-only primary model would model mouse blocking, but would discard one control-only and one E-Stim-only mouse.

**Rationale**: The requested condition-level unit is Mouse × Condition pseudobulk. The unpaired primary design preserves all biological samples and documents the limitation that pairing is not modeled.

**Consequences**: Report primary results with the design limitation. Use the paired sensitivity files to check whether headline marker directions depend on the unpaired mice.

### [2026-07-04] Keep MG-selected clustering at 30 PCs and resolution 0.3

**Tags**: clustering, mg-selected, pFlog, parameters

**Context**: After removing selected-source clusters 4 (microglia score plus Cdkn1b-high), 7 (microglia score), and 9/10 (Cdkn1b-high), the MG-selected PFlog object was reclustered across 20, 30, and 50 PCs; Leiden resolutions 0.3, 0.5, and 0.8; and with or without cell-cycle HVG filtering. The marker exclusion rule required top class microglia/photoreceptor, top score ≥ 0.5, and score margin ≥ 0.25. The Cdkn1b rule required expression and detection BH-adjusted Wilcoxon q < 0.05 plus detection fraction ≥ 0.20.

**Decision**: Use PFlog, cell-cycle HVGs retained, 30 PCs, and Leiden resolution 0.3 for downstream MG-selected figures and DE.

**Alternatives considered**: 20 PCs is close to the elbow; 50 PCs carries more tail variance; 0.5 and 0.8 split major structures more aggressively.

**Rationale**: Thirty PCs is safely past the elbow, and resolution 0.3 preserves broad MG-selected structure for manuscript-scale figures.

**Consequences**: Higher resolutions remain available for later subclustering, but not as the current manuscript-level branch.
