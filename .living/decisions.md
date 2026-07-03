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
