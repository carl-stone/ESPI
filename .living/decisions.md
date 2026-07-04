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

**Reporting addendum 2026-07-04**: For both DE and DD, treat concordant primary and paired-sensitivity effects as the strongest design-robust signals. Interpret paired-only hits as within-paired-mice sensitivity results, not as globally more accurate than the all-sample primary condition model.

### [2026-07-04] Keep MG-selected clustering at 30 PCs and resolution 0.3

**Tags**: clustering, mg-selected, pFlog, parameters

**Context**: After removing selected-source clusters 4 (microglia score plus Cdkn1b-high), 7 (microglia score), and 9/10 (Cdkn1b-high), the MG-selected PFlog object was reclustered across 20, 30, and 50 PCs; Leiden resolutions 0.3, 0.5, and 0.8; and with or without cell-cycle HVG filtering. The marker exclusion rule required top class microglia/photoreceptor, top score ≥ 0.5, and score margin ≥ 0.25. The Cdkn1b rule required expression and detection BH-adjusted Wilcoxon q < 0.05 plus detection fraction ≥ 0.20.

**Decision**: Use PFlog, cell-cycle HVGs retained, 30 PCs, and Leiden resolution 0.3 for downstream MG-selected figures and DE.

**Alternatives considered**: 20 PCs is close to the elbow; 50 PCs carries more tail variance; 0.5 and 0.8 split major structures more aggressively.

**Rationale**: Thirty PCs is safely past the elbow, and resolution 0.3 preserves broad MG-selected structure for manuscript-scale figures.

**Consequences**: Higher resolutions remain available for later subclustering, but not as the current manuscript-level branch.

### [2026-07-04] Keep repo-local Mycelium skills synced from OMP

**Tags**: mycelium, hooks, skills, reproducibility

**Context**: The Mycelium command protocols refer to repo-relative paths such as `skills/core/scripts/validate_structure.py`, but the installed OMP plugin stores those files under a machine-local plugin cache. A symlink would work locally but would break for other clones.

**Decision**: Keep a real repo-local copy of `skills/core/`, add `tools/sync-mycelium-skills-core.py`, and run it quietly at session start before local Mycelium hooks run.

**Alternatives considered**: Using only the OMP cache avoids duplication but breaks repo-relative protocol paths. A symlink is simpler but points at a machine-local cache path. Running the sync before every hook would stay fresher but add unnecessary hook overhead and context noise.

**Rationale**: A session-start sync keeps `skills/core/` current when the OMP plugin changes while preserving portable repo-relative paths for scripts and hooks.

**Consequences**: Local hooks now call `skills/core/hooks/...`; if the OMP plugin updates, the next session start refreshes `skills/core/` excluding nuisance files such as `.DS_Store`, `__pycache__`, `.pytest_cache`, and Python bytecode.

### [2026-07-04] Bridge Mycelium hooks through OMP extensions

**Tags**: mycelium, hooks, omp, session-state

**Context**: The repo's `.claude/settings.local.json` configures Claude Code shell hooks, but OMP does not run those hooks directly. That left Mycelium read-access logging and `.claude/last-session.md` updates inactive in OMP sessions.

**Decision**: Add a project-local OMP extension that mirrors the relevant Mycelium hook events: session start runs the repo-local skill sync and health hook, read tool results log `.living/` access, write/edit tool results record session activity, bash tool results invoke the Mycelium post-action hook, and session stop invokes the stop hook plus a deterministic five-section `last-session.md` fallback.

**Alternatives considered**: Rewriting `.claude/settings.local.json` is insufficient because OMP ignores Claude Code hook config. Editing the synced `skills/core/hooks/` files would be overwritten by the next skill sync. A repo-local extension preserves the OMP-specific adapter outside synced Mycelium source.

**Rationale**: OMP extensions are the native event surface, and the existing generated-file guard already uses the same adapter pattern.

**Consequences**: OMP sessions now get Mycelium read tracking and session-resume state without depending on Claude Code hook execution.

### [2026-07-04] Filter Mycelium maintenance commands before post-action hooks

**Tags**: mycelium, hooks, omp, session-state

**Context**: The OMP Mycelium adapter and Claude hook settings could forward Bash commands to `skills/core/hooks/mycelium-post-action.sh`. That synced hook treats Python script execution as significant work, so maintenance commands such as `python3 skills/core/scripts/generate_index.py` and `python3 skills/core/scripts/validate_structure.py` could recreate `.claude/mycelium-reminded.tmp` after triage was already complete.

**Decision**: Add a repo-owned post-action wrapper at `.agents/hooks/adapters/mycelium-post-action-wrapper.sh` that skips only Mycelium maintenance commands under `skills/core/scripts/` and `tools/sync-mycelium-skills-core.py`, then delegates all other Bash payloads to the synced hook. Point both `.claude/settings.local.json` and the OMP adapter at the wrapper.

**Alternatives considered**: Editing `skills/core/hooks/mycelium-post-action.sh` would fix the source hook but would be overwritten by the repo-local skill sync. Clearing `.claude/mycelium-reminded.tmp` manually fixes one session but does not prevent recurrence. Suppressing all Python or Bash post-action forwarding would hide real analysis work. An adapter-only guard did not cover all live hook paths in this session.

**Rationale**: The wrapper and adapter are repo-owned and survive skill sync. A narrow path-based guard prevents Mycelium bookkeeping from starting a new post-action cycle while preserving reminders for legitimate R and analysis scripts.

**Consequences**: Future OMP or Claude-hook sessions should no longer get false-positive stop blocks after running Mycelium maintenance commands. The wrapper was directly verified to skip `validate_structure.py` and pass through an `Rscript` analysis command; the running OMP session may still need a reload to pick up adapter changes.

### [2026-07-04] Mirror complete Mycelium hook behavior in OMP

**Tags**: mycelium, hooks, omp, data-lineage, r

**Context**: A behavioral check showed the OMP adapter ran the post-action wrapper but discarded its `additionalContext`, did not call the data-lineage hook bundle, and the synced lineage extractor detected `Rscript` commands without recognizing ESPI's R I/O calls.

**Decision**: Return post-action wrapper context from the OMP `tool_result` handler, run the data tracker and data-lineage stop hook from the OMP adapter, and route data tracking through a repo-owned wrapper that adds R I/O regexes while leaving synced `skills/core/` files untouched.

**Alternatives considered**: Reporting the gap only would leave Mycelium reminders and lineage incomplete under OMP. Editing `skills/core/` directly would be overwritten by the next sync from the plugin cache.

**Rationale**: Repo-owned adapters survive skill sync and are already the project-local seam between OMP and Mycelium.

**Consequences**: OMP and Claude-hook sessions now share the same wrapper path for data tracking. R lineage records include source-level I/O expressions such as `input` and `out_path`; they do not resolve every runtime path computed from CLI arguments.

### [2026-07-04] Reset prior-session Mycelium sentinels at SessionStart

**Tags**: mycelium, hooks, omp, session-state

**Context**: OMP can run the Mycelium stop hook at turn boundaries. A prior session's `.claude/mycelium-reminded.tmp` and `.claude/mycelium-session-activity.tmp` could survive when a new session started less than one hour later, so the first stop event of the new session inherited an already-old reminder and falsely blocked with `.living/ not updated`.

**Decision**: Treat work sentinels older than the current `.claude/session-start-ts.tmp` as prior-session state. Clear those sentinels during `mycelium-health.sh` SessionStart and make `mycelium-stop-check.sh` self-heal by ignoring and removing older reminder/activity files if SessionStart cleanup did not run.

**Alternatives considered**: Increasing the five-minute stop-hook debounce would mask the symptom but still allow cross-session state bleed. Clearing all sentinels unconditionally would risk erasing live primary-session work when a subagent or nested hook runs.

**Rationale**: The session-start timestamp is the boundary that distinguishes current work from prior-session leftovers. Comparing sentinel timestamps to that boundary fixes the false positive without weakening enforcement for real current-session work.

**Consequences**: Back-to-back OMP sessions should no longer block immediately from stale reminder/activity files. The stop hook still blocks current-session work older than the debounce window when no `.living/` triage file was updated.

### [2026-07-04] Match differential-detection and DE gene universes

**Tags**: differential-expression, differential-detection, mg-selected, filtering

**Context**: The MG-selected DESeq2 analysis filters genes to pseudobulk row sums >= 10 before DE testing, but the limma differential-detection analysis originally tested every gene in the Seurat counts matrix.

**Decision**: Filter primary differential detection to the primary DESeq2 tested genes, and filter paired-sensitivity differential detection to the paired-sensitivity DESeq2 tested genes.

**Alternatives considered**: Testing all detected/undetected gene rows would maximize the DD search space but use a different multiple-testing background from DE. Adding a separate detection-specific prevalence filter would be defensible but would make DE and DD less directly comparable.

**Rationale**: The notebook compares DE and DD as paired summaries of the same MG-selected branch. Matching the tested gene universe keeps the multiple-testing background aligned and removes genes that failed the DE count floor.

**Consequences**: Primary DD now tests 24,514 genes, matching primary DE. Paired-sensitivity DD now tests 22,663 genes, matching paired DE. `numbers.json` records both DD tested-gene counts and DD hit counts. The primary DD result remains negative, while paired-sensitivity DD becomes exploratory because it has only four samples and one residual degree of freedom.

**Superseded 2026-07-04 by "Replace limma logit DD with muscat edgeR_NB_optim":** DD gene universe is now defined by muscat's internal 90% detection filter, not the DE tested-gene set. See the following decision entry.

### [2026-07-04] Replace limma logit DD with muscat edgeR_NB_optim

**Tags**: differential-detection, mg-selected, muscat, mycelium

**Context**: The prior DD implementation used a limma empirical-logit model with a 0.5 pseudocount. It did not adjust for per-cell library depth or per-sample cellular detection rate, which is a known confound in low-depth PipSeq data.

**Decision**: Replace primary and paired-sensitivity DD in `scripts/run-mg-selected-de.R` with `muscat::pbDS(method = "DD")` (equivalently `muscat::pbDD()`), the `edgeR_NB_optim` workflow from Gilis et al., BMC Genomics 26:886 (2025).

**Alternatives considered**: MAST hurdle with `(1 | mouse)`; a hand-rolled edgeR QL on detected-cell counts with CDR offset. Both are cell-level or bespoke; muscat implements the published, benchmarked workflow directly at Mouse × Condition sample level.

**Rationale**: The muscat DD workflow uses an internal 90%-detection filter and CDR normalization offset (`log(nc * of)`) with edgeR QL robust dispersion, which is the published answer to the per-cell depth confound and preserves the Mouse × Condition statistical unit.

**Consequences**: DD outputs at `DEG_DIR/mg_selected/detection_*` use edgeR result columns (`logFC`, `logCPM`, `F`, `p_val`, `p_adj.loc` mapped to `pvalue`, `padj`). The DD tested-gene universe is set by muscat's internal 90% detection filter, not the DE gene set (supersedes "Match differential-detection and DE gene universes"). Verification showed this muscat-native universe is larger than the DE universe for this sparse PipSeq dataset: primary DD tested 36,468 genes with 0 FDR-significant hits, and paired-sensitivity DD tested 34,880 genes with 40 FDR-significant hits. The prior limma paired-sensitivity DD hit count was 108, so the hit set changed materially. `numbers.json` records `dd_method` and per-analysis DD tested-gene counts.
