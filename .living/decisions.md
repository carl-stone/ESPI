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
**Superseded 2026-07-12 by "Use 20 PCs and Leiden resolution 0.3 for the counts-derived PFlog run":** The 30-PC choice and associated cluster IDs describe the earlier legacy-source analysis. The current counts-derived analysis uses the no-cell-cycle-HVG PFlog branch at 20 PCs and resolution 0.3; see the 2026-07-13 rebuild decision for current results.

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

### [2026-07-04] Run MG-selected cluster markers on data layer with positive detection filter

**Tags**: marker-analysis, findallmarkers, mg-selected, seurat

**Context**: The MG-selected chosen clustering is PFlog-derived, but `FindAllMarkers()` on the PFlog layer produced sign-incoherent positive markers because Seurat's default fold-change math assumes log-normalized expression.

**Decision**: Use the retained PFlog-derived Leiden labels as the marker identities, run `Seurat::FindAllMarkers()` on the RNA `data` layer, and keep only rows with `pct.1 > pct.2` for descriptive positive marker ranking.

**Alternatives considered**: Running `FindAllMarkers()` directly on PFlog would keep the visualization and clustering scale aligned but was empirically invalid. Supplying a PFlog-specific mean function fixed some ranking behavior but still left clusters without clean positive detection markers. Manually merging cluster 2 before marker ranking lacked a specific supported merge partner.

**Rationale**: The clustering choice and marker expression scale answer different questions. The PFlog branch defines the graph structure; the log-normalized `data` layer gives Seurat-compatible marker fold changes. The positive detection filter prevents anti-markers from entering manuscript-facing marker tables.

**Consequences**: Cluster 2 has no retained positive detection-enriched marker genes and should not be treated as a marker-defined interpreted identity. Cluster 8 has only 15 cells, so its many retained markers are descriptive and potentially unstable. Marker p-values remain cell-level descriptive ranks, not Mouse × Condition condition-effect evidence.

### [2026-07-04] Show MG-selected cell-cycle-filtered figures as complementary analysis

**Tags**: mg-selected, plotting, notebook, differential-expression

**Context**: The manuscript notebook used the PFlog, no-cell-cycle-filtered MG-selected branch for downstream clustering and DE, but the revised manuscript needs matching visual evidence from the cell-cycle-filtered MG-selected branch.

**Decision**: Keep the no-cell-cycle-filtered branch as the primary downstream DE branch, show the cell-cycle-filtered cluster-identity UMAP and feature UMAP grid as complementary clustering figures, drop the redundant cell-cycle-filtered marker heatmap from the notebook, and summarize DE/DD jointly with an effect-size scatter plot that inner-joins genes tested in both workflows.

**Alternatives considered**: Replacing the primary branch with the cell-cycle-filtered branch would change the DE input after downstream interpretation was already grounded. Plotting DE and DD in separate panels would hide whether expression-magnitude and detection-fraction effects agree for the same genes.

**Rationale**: Complementary filtered UMAP figures show sensitivity to cell-cycle HVG handling without repeating a marker heatmap that did not add interpretive signal. The explicit inner join makes the DE/DD comparison honest because the muscat DD and DESeq2 gene universes differ.

**Consequences**: `notebook/sc_analysis.qmd` now includes both MG-selected visual branches, while `scripts/run-mg-selected-de.R` writes `mg_selected_de_dd_effect_scatter.(png|pdf)` and a notebook symlink. Genes tested by only one workflow are not plotted in the DE/DD scatter.

**Revision 2026-07-04**: Retain and regenerate the cell-cycle-filtered marker heatmap in the notebook rather than dropping it; the missing symlink was an output gap, not a scientific decision to remove the panel.

### [2026-07-04] Scale MG-selected feature UMAP panels per gene

**Tags**: mg-selected, plotting, notebook, umap

**Context**: The MG-selected feature UMAP grids compare spatial expression patterns across marker genes, but raw PFlog ranges differ by gene and gave each panel its own color scale.

**Decision**: Min-max scale each feature UMAP panel to 0–1 within gene and branch, use identical color limits across panels, collect the patchwork guides into one shared legend labeled `Scaled expression`, and use square coordinate limits for every feature panel.

**Alternatives considered**: Keeping raw PFlog color scales preserves absolute expression magnitude but makes the 3 × 3 grid harder to scan. Using a global range across all genes would be dominated by high-range genes and hide lower-range spatial patterns.

**Rationale**: The grid is a pattern-visualization figure, not a quantitative expression-magnitude comparison. Per-gene scaling makes within-gene spatial localization comparable across panels, the shared 0–1 legend keeps the figure visually clean, and square coordinate limits prevent the panels from reading as tall rectangles.

**Consequences**: The color intensity in the feature UMAP grids now represents relative expression within each gene, not absolute PFlog expression across genes. Constant-expression genes fall back to scaled value 0; the current nine plotted features are not constant in either MG-selected branch.


### [2026-07-04] Use binned diverging guide for MG-selected marker dot plot

**Tags**: mg-selected, plotting, notebook, marker-analysis

**Context**: The MG-selected marker dot plot used continuous row z-score colors, which made it hard to read which expression range a dot occupied.

**Decision**: Keep dot size as detected-cell percentage, map color to clipped mean expression row z-scores, and display the color scale with a binned diverging guide whose end labels are `<= -2` and `>= 2`.

**Alternatives considered**: A discrete unordered palette would make bins easy to count but would lose low-to-high direction. A continuous gradient preserves direction but makes range lookup harder.

**Rationale**: The dot plot is a descriptive marker-summary figure. Binned diverging colors preserve the low/high expression direction while making the z-score range visually explicit.

**Consequences**: Figure 11 now communicates approximate z-score bins rather than continuous color values. Values outside ±2 remain clipped into the endpoint bins.

**Superseded 2026-07-04 by "Use six endpoint-aware bins for MG-selected marker dot plot":** row z-scores are no longer pre-clipped before plotting; endpoint squishing now happens in the color scale.

### [2026-07-04] Use six endpoint-aware bins for MG-selected marker dot plot

**Tags**: mg-selected, plotting, notebook, marker-analysis

**Context**: The previous Figure 11 binned guide collapsed values at or beyond `+/-2` into the same endpoint value before plotting, so the darkest endpoint colors could not distinguish `1` to `2` from `>= 2` or `-2` to `-1` from `<= -2`.

**Decision**: Keep raw row z-scores after zero-SD/NA handling, use a six-color diverging binned scale, and squish out-of-range values into endpoint bins at plot time.

**Alternatives considered**: Pre-clipping z-scores to `+/-2` keeps the scale bounded but hides meaningful extremes. A continuous gradient makes outliers visible but makes the requested range bins harder to read.

**Rationale**: The marker dot plot is a descriptive figure. Six diverging bins preserve direction, show `1` to `2` separately from `>= 2`, and show the three requested blue-side ranges.

**Consequences**: The guide now shows endpoint labels `<= -2` and `>= 2`, and values outside the display limits render with the endpoint colors instead of disappearing as missing values.

### [2026-07-04] Treat cluster abundance Fisher CLR as descriptive

**Tags**: mg-selected, abundance, fisher-exact, clr, plotting

**Context**: The requested cluster abundance summary mirrors `compute_cluster_abundance()` from the separate `megan_sc` repo, which pools raw cell counts across mice within each condition before running Fisher exact tests.

**Decision**: Port the Fisher/CLR logic into ESPI-owned `R/` code, plot vertical CLR log2 enrichment bars with E-Stim-enriched clusters in `#e31a8c` and E-Stim-depleted clusters in `#2166ac`, and label the output as a pooled cell-level descriptive summary.

**Alternatives considered**: Replacing Fisher with a per-mouse model would better match the Mouse × Condition statistical unit used for DE/DD, but it would not implement the method the user asked to calculate.

**Rationale**: The plot answers a descriptive cluster-composition question and should be visually available in the notebook, while the caveat prevents it from being presented as equivalent to pseudobulk condition-level inference.

**Consequences**: The abundance figure is useful for cluster-level pattern review, but DE/DD remain the condition-level analyses with Mouse × Condition pseudobulk samples as the statistical unit.

### [2026-07-05] Use Mouse × Condition cluster counts for future abundance inference

**Tags**: mg-selected, abundance, paired-design, causal-inference, compositional-data

**Context**: A targeted ideation panel reviewed how to test whether E-Stim changes MG-selected cluster proportions without using the existing pooled cell-level Fisher/CLR plot as condition-level evidence.

**Decision**: Treat future inferential cluster-abundance analyses as Mouse × Condition sample-level analyses. Start from a sample × cluster count table, use paired mice 10 and 3 as the anchor, use mice 30 and 33 as unpaired context or sensitivity, and keep pooled Fisher/CLR as descriptive only.

**Alternatives considered**: Reusing pooled Fisher tests would be simple but violates the biological replicate unit. Fully model-based multinomial or neighborhood differential-abundance approaches are better robustness layers but are heavier than the immediate paired/sample-level summaries.

**Rationale**: The ESPI statistical unit is Mouse × Condition. Sample-level proportions or log-ratios align the abundance question with the DE/DD design, while explicit caveats about frozen/data-derived clusters address circularity from testing a clustering learned on the same dataset.

**Consequences**: Any future manuscript-facing cluster-proportion claim should cite a mouse-level analysis, not the descriptive Fisher/CLR plot. Stronger causal or cell-state claims need frozen labels, cross-fit/alternative clustering sensitivity, or validation data.

### [2026-07-05] Make cluster-proportion inference design-restricted and exact

**Tags**: mg-selected, abundance, paired-design, randomization, notebook

**Context**: The pooled Fisher/CLR cluster-abundance plot was already labeled descriptive, but the notebook still needed a primary sample-level screen for MG-selected cluster-proportion shifts.

**Decision**: Add exported helpers for Mouse × Condition sample cluster proportions, exact paired sign-flip randomization, and by-mouse proportion plotting. Use paired mice 10 and 3 as the primary statistic and report the mouse 30 E-Stim-only plus mouse 33 control-only block as a clearly labeled paired-plus-singleton sensitivity.

**Alternatives considered**: Keeping only pooled Fisher/CLR would be simpler but would keep cell counts as the inferential unit. A fuller model-based differential-abundance method is heavier than the current small design supports as a first screen.

**Rationale**: Exact enumeration matches the tiny design and makes the p-value floor transparent: 0.5 for paired-only and 0.25 for paired-plus-singleton. Reporting effect size and directional consistency is more honest than treating coarse p-values as conventional significance evidence.

**Consequences**: The notebook now treats Fisher/CLR bars as descriptive only and points per-cluster inference to `TABLE_DIR/mg_selected/mg_selected_cluster_proportion_randomization_pflog_mg_selected_no_filter_cc_dims30_res0.3.tsv`.

### [2026-07-05] Batch open TODOs by shared code surfaces

**Tags**: todo, planning, mg-selected, plotting, notebook

**Context**: The open TODO registry contained five items: shared plotting palette, condition-label standardization, stable cross-references, cluster-abundance plot-helper split, and MG-selected marker/p27 heatmaps. Palette and label changes both touched plotting code, while the heatmap task added new scientific/report outputs.

**Decision**: Batch the palette, contrast-label, plot-helper split, and stable-reference cleanup together as one presentation-infrastructure pass, then implement the MG-selected marker/p27 heatmaps as a separate analysis/report pass.

**Alternatives considered**: Batching strictly by registry category would separate infrastructure, refactor, and writing tasks even though they touch the same plot labels and notebook prose. Combining all five TODOs would bury a new biological heatmap analysis inside refactor and style changes.

**Rationale**: The first batch can normalize shared plotting semantics and captions once before regenerating figures. The second batch stays reviewable as a scientific addition and can reuse the standardized palette and contrast labels.

**Consequences**: Future implementation should produce at least two focused commits: one for presentation infrastructure/reference cleanup and one for MG-selected heatmap analysis.

### [2026-07-05] Require stable cross-references going forward

**Tags**: conventions, cross-references, notebook, reporting, session-logs

**Context**: Batch 1 presentation cleanup found fragile auto-numbered figure references in records after notebook panels had shifted.

**Decision**: Use stable Quarto figure IDs (`#fig-...`), file basenames, or decision headings in new notebook prose, session logs, review reports, and Mycelium records. Preserve historical `.living/log/` entries, review evidence, and old decision entries unless their content is otherwise being edited; do not rewrite old records solely to remove `Figure N`.

**Alternatives considered**: Rewriting every historical mention would make records look cleaner but would erase audit and review context. Continuing to use auto-numbered labels would keep new records brittle.

**Rationale**: Stable handles survive panel insertion, removal, and reorderings, while historical records should remain faithful to what the session or review observed at the time.

**Consequences**: New prose and Mycelium records should cite stable handles. Existing `Figure N` mentions may remain in older logs, reviews, and decisions as historical evidence.

### [2026-07-05] Centralize presentation palette and contrast display labels

**Tags**: plotting, mg-selected, notebook, contrast-labels, palette

**Context**: Batch 1 presentation cleanup needed a shared three-color plotting palette and consistent rendered contrast text across the MG-selected abundance and DE/DD figures.

**Decision**: Keep metadata condition labels unchanged (`p27CKO +EStim`, `p27CKO`), move them to `R/conditions.R`, add display-only labels for human-facing text, and render plot contrast labels from `CONTRAST_DISPLAY_LABEL`. Use `palette_analysis_three` as the shared blue/gray/pink low/mid/high palette and derive `palette_dotplot_pair` from it.

**Alternatives considered**: Renaming metadata values would make display text prettier but would break existing Seurat metadata and tripwire contrast declarations. Leaving each script with local colors and labels would preserve current output but keep drift-prone hard-coded strings.

**Rationale**: Separating metadata values from display labels keeps computation stable while giving plots a single source of truth for manuscript-facing labels and colors.

**Consequences**: Future plotting code should import display labels and palette values from the package namespace instead of hard-coding condition text or endpoint colors.

### [2026-07-05] Protect Mycelium hook provenance with project-owned wrappers

**Tags**: mycelium, hooks, provenance, session-logs, omp

**Context**: Synced Mycelium hooks can regenerate `LOG_REGISTRY.md`, `.living/INDEX.md`, and `.claude/last-session.md` with generic placeholders after manual semantic records have been repaired.

**Decision**: Keep synced `skills/core/*` untouched and route OMP hook execution through project-owned wrappers in `.agents/hooks/adapters/`. The wrappers snapshot provenance files, run the synced hook, then invoke `tools/mycelium-provenance-guard.py` to restore semantic registry rows, INDEX knowledge-summary ordering, and last-session content only when the hook output regresses them.

**Alternatives considered**: Patching `skills/core/*` would be overwritten by `tools/sync-mycelium-skills-core.py`. Relying on manual post-hook inspection had already failed repeatedly. Committing `.claude/settings.local.json` would conflict with the repo rule that `.claude/` is machine-local.

**Rationale**: A project-owned wrapper survives core syncs, keeps the synced plugin cache as the source for upstream hook behavior, and adds a narrow ESPI guard around provenance files without changing analysis code.

**Consequences**: Future OMP hook changes should preserve wrapper routing. Claude Code `settings.local.json` routing remains local-only unless hook setup regenerates it on another clone.

### [2026-07-05] Use within-sample p27 permutation z-scores in cluster heatmaps

**Tags**: mg-selected, plotting, heatmap, marker-analysis, randomization

**Context**: The MG-selected notebook needed compact per-cluster heatmaps that show cell-type marker-module context and p27 status for both the full source clustering and the retained MG-selected clustering.

**Decision**: Plot cluster-level `cell_type_marker_genes` module scores as row z-scores, and add a top strip for p27 enrichment computed by permuting cluster labels within each Mouse × Condition sample. Keep the permutation statistic in exported R helpers and keep ComplexHeatmap drawing and artifact writing in `scripts/plot-cluster-marker-heatmaps.R`.

**Alternatives considered**: A pooled p27 mean or pooled between-cluster permutation would be simpler but would ignore sample composition. Putting all drawing and statistics in a single script would match the older marker heatmap, but the p27 permutation is reusable, testable analysis logic.

**Rationale**: Within-sample label shuffling preserves each Mouse × Condition sample's p27 expression distribution and cluster-size composition while asking whether a cluster has more or less p27 expression than expected under sample-aware relabeling. The helper/script split keeps statistical code unit-testable without turning one-off figure layout into package API.

**Consequences**: The notebook now has two per-cluster module/p27 heatmaps with p27 strip scales fitted per figure. Future shared-scale comparisons would need a coupled two-object script or an explicit shared z-limit option.

### [2026-07-06] Label DE/DD scatter only with significant curated markers

**Tags**: differential-expression, differential-detection, mg-selected, plotting, marker-genes

**Context**: The MG-selected DE/DD effect-size scatter had been changed to use `ggrepel::geom_text_repel()` for every gene significant in one or both DE/DD tests. That rendered hundreds of eligible labels and the visible labels were not from `cell_type_marker_genes`, while the notebook prose now frames the labels as curated marker genes.

**Decision**: Restrict DE/DD scatter labels to the intersection of genes directly listed in `cell_type_marker_genes` and genes with FDR < 0.05 in DE, DD, or both. Keep all tested genes as points and keep FDR category as point color.

**Alternatives considered**: Labeling every significant gene preserves a broad hit overview but hides the curated marker narrative. Labeling every curated marker would add nonsignificant genes back to the plot. Reusing `make_marker_table()` would also include standalone `Cdkn1b`, which is useful for overlap reports but is not part of the curated cell-type marker list.

**Rationale**: The figure is a marker-context plot, not a top-hit labeling plot. Intersecting curated markers with the existing significance filter preserves the user's requested marker focus without labeling nonsignificant genes.

**Consequences**: `scripts/run-mg-selected-de.R` now derives `curated_marker` directly from `cell_type_marker_genes`, and the rendered scatter labels only `Ccn1`, `Glul`, `Grm6`, `Hes6`, `Mcm2`, `Mcm6`, `Pcna`, and `Serpina3n` for the current outputs.

### [2026-07-06] Show condition-colored MG-selected UMAPs beside cluster UMAPs

**Tags**: mg-selected, plotting, notebook, umap, condition

**Context**: The notebook needed condition-colored UMAPs next to the existing MG-selected cluster UMAPs for both the cell-cycle-HVG-retained and cell-cycle-HVG-filtered PFlog branches.

**Decision**: Add `mg_selected_condition_umap_*` PNG/PDF outputs to `scripts/plot-mg-selected-figures.R`, color cells by `Condition`, use the shared `palette_dotplot_pair` colors with display labels, and insert each condition UMAP immediately after its matching cluster UMAP in `notebook/sc_analysis.qmd`.

**Alternatives considered**: Reusing cluster UMAP legend placement would put a wide condition legend outside the panel and waste space. Making a combined two-panel figure would be compact but would move away from the user's requested exact notebook placement after each existing figure.

**Rationale**: Adjacent cluster and condition UMAPs make it easy to see whether condition structure aligns with the chosen clustering without changing the existing figure order. In-panel legends use empty space in each branch-specific UMAP: top right for the cell-cycle-retained branch and bottom right for the cell-cycle-filtered branch.

**Consequences**: Running `scripts/plot-mg-selected-figures.R` now regenerates and links an additional condition UMAP for whichever MG-selected branch is requested, and `notebook/sc_analysis.html` must be rerendered after these figure outputs change.

### [2026-07-06] Track Ed DEG and enrichment handoff as a TODO

**Tags**: todo, mg-selected, differential-expression, enrichment, collaboration

**Context**: The user asked to capture a future work item for sending Ed a compact set of current MG-selected DEG and GSEA/GO result tables.

**Decision**: Add an open medium-priority analysis TODO with explicit acceptance criteria and source output locations instead of generating the bundle in the same turn.

**Alternatives considered**: Creating the CSV bundle immediately would have produced new analysis artifacts without an explicit build request. Leaving the request only in chat would make the handoff easy to lose.

**Rationale**: The TODO registry is the project-visible place for actionable future work. Recording source DEG and enrichment directories plus acceptance criteria makes the later bundle task reproducible and easy to scope.

**Consequences**: Future work should convert the selected MG-selected TSV outputs to clearly named CSVs in a small sendable directory or archive, then update the TODO status when complete.

### [2026-07-06] Bundle Ed handoff as significant DE/DD subsets plus explicit GO/GSEA term sets

**Tags**: todo, mg-selected, differential-expression, differential-detection, enrichment, collaboration

**Context**: The Ed handoff TODO asked for a small sendable CSV bundle of current MG-selected DE/DD and GSEA/GO outputs.

**Decision**: Put the bundle under the Box `exports/` directory with a matching zip archive. For DE/DD, include small design/sample/run summaries, significant DESeq2 subsets, marker-overlap tables, and FDR-filtered differential-detection subsets rather than the multi-MB all-gene result dumps. For GO/GSEA, include both FDR-filtered term CSVs and full all-term conversions, with filenames and `manifest.csv` distinguishing filtered from full tables.

**Alternatives considered**: Including every all-gene DE/DD result table would be complete but unnecessarily heavy for an external handoff. Shipping only filtered GO/GSEA hits would be smaller but would lose full pathway rankings and the GSEA Entrez mapping needed to interpret core-enrichment IDs.

**Rationale**: This keeps the archive small while preserving the pieces Ed is most likely to review: significant gene-level results, marker-context tables, pathway hits, and full pathway context where rankings/mappings matter.

**Consequences**: The bundle explicitly contains a zero-row primary differential-detection significant-gene CSV because no primary DD genes pass `padj < 0.05` in the current outputs; the manifest records that this is a genuine null result, not a conversion failure.

### [2026-07-09] Use just as the ESPI command interface

**Tags**: tooling, reproducibility, scripts, r-package

**Context**: ESPI had many raw `Rscript`, `devtools`, Quarto, Air, and tripwire commands spread across docs and script headers. The user asked to install `just` and route scripts and pipeline build steps through `just` commands.

**Decision**: Add a root lowercase `justfile` as the preferred command interface, with recipes for package loading/documentation, README rebuilds, formatting, preprocessing, clustering dry-runs/runs, marker heatmaps, MG-selected marker/figure/DE outputs, notebook rendering, and tripwires. Keep the recipes as direct wrappers over existing scripts rather than moving pipeline logic into `just`.

**Alternatives considered**: Leaving raw command snippets in docs would preserve the status quo but keep the workflow harder to discover. Moving orchestration into new R helpers would add package surface area for a one-off analysis pipeline.

**Rationale**: `just --list` gives a small discoverable interface while existing R scripts remain the implementation and source of scientific behavior.

**Consequences**: Future routine command examples should prefer `just` recipes. New pipeline scripts should get a matching recipe when they become a repeated build step.

### [2026-07-09] Derive raw-count inputs from `DATA_ROOT_DIR` and stop at an in-memory object

**Tags**: scripts, single-cell, reproducibility, data-lineage

**Context**: The new first pipeline stage loads six 10X count-matrix samples and requires stable sample metadata without creating a downstream artifact prematurely.

**Decision**: Build the raw-input directory from `DATA_ROOT_DIR` plus `data/input/Raw Matrices`, use the metadata table as the sample manifest, and validate the combined Seurat object in memory without saving it.

**Alternatives considered**: A hardcoded machine path, filesystem discovery without metadata, or an immediate RDS write would make the stage less portable, less explicit, or commit to an artifact before the next stage is defined.

**Rationale**: The configured data root is the portable source of location, and the metadata table is the authoritative mapping of sample identity and fields.

**Consequences**: `scripts/process-counts.R` currently has no artifact output; saving or preprocessing begins only with a defined next pipeline stage.

### [2026-07-09] Persist the raw Seurat object at the user-chosen input path

**Tags**: scripts, data-lineage, reproducibility, provenance

**Context**: The raw-count loader now produces a verified Seurat object that preprocessing must consume.

**Decision**: Persist `sobj_raw.rds` at `DATA_ROOT_DIR/data/input/sobj_raw.rds`; invoke preprocessing with that path through explicit `--input`.

**Alternatives considered**: Changing preprocessing's default would silently redefine its existing input contract.

**Rationale**: The user selected this durable raw-object location, while preprocessing's default points elsewhere.

**Consequences**: This is an explicit handoff, not a new preprocessing default.

### [2026-07-09] Limit QC filtering to low-complexity count thresholds
**Superseded 2026-07-09 by "Use complete mixed-label mitochondrial metric with a data-specific >20% cutoff":** The prior `^mt-` metric counted only `mt-Rnr1` and `mt-Rnr2`; it was not the complete mitochondrial set.


**Tags**: qc-filtering, scripts, data-lineage, validation, provenance

**Context**: The raw Seurat object must pass through a documented cell-QC stage before preprocessing, but the custom reference represents mitochondrial genes only as `mt-Rnr1` and `mt-Rnr2`.

**Decision**: Retain cells only when `nFeature_RNA >= 50` and `nCount_RNA >= 100`. Keep `percent.mt` and `percent.ribo` as diagnostic metadata; apply no mitochondrial cutoff. Apply no sample-, empty-droplet-, ambient-RNA-, or droplet-doublet gate. Treat PIPseeker upstream matrix selection as separate from this cell-QC stage.

**Alternatives considered**: A mitochondrial threshold or additional sample/droplet/doublet gates would add filtering criteria that the available reference or project inputs do not support.

**Rationale**: The two count thresholds remove low-complexity observations while avoiding an unsupported mitochondrial criterion from an incomplete mitochondrial reference.

**Consequences**: The verified output is `INPUT_OBJECT_DIR/sobj_qc_filtered.rds`, retaining 22,751 of 983,903 cells across all six samples. Preprocessing must receive this path explicitly with `--input`; this verification creates no scientific result.

### [2026-07-09] Use complete mixed-label mitochondrial metric with a data-specific >20% cutoff
**Superseded 2026-07-13 by "Use liberal floors before sample-specific MAD and doublet filtering":** The global thresholds remain the first gate, but the saved QC object now also requires the sample-specific MAD criteria and a scDblFinder singlet call.

**Tags**: qc-filtering, mitochondrial, provenance, validation, data-lineage

**Context**: The prior `^mt-` metric incorrectly treated its two matches, `mt-Rnr1` and `mt-Rnr2`, as the complete mitochondrial set. A feature audit found 37 mitochondrial rows in the raw object: those two rRNAs, 13 unprefixed uppercase protein-coding identifiers (`ND1`, `ND2`, `COX1`, `COX2`, `ATP8`, `ATP6`, `COX3`, `ND3`, `ND4L`, `ND4`, `ND5`, `ND6`, `CYTB`), and 22 `Trn*` identifiers. The earlier lossless `Read10X()` audit establishes the R import but not PIPseeker upstream cell calling; this feature audit shows no mitochondrial features were lost by the reference or import.

**Decision**: Calculate `percent.mt` from all 37 observed mitochondrial features and retain cells with `nFeature_RNA >= 50`, `nCount_RNA >= 100`, and `percent.mt <= 20`. Keep `percent.ribo` diagnostic only. Apply no ribosomal, high-complexity, sample, empty-droplet, ambient-RNA, or doublet filter.

**Alternatives considered**: Keeping the two-rRNA metric would preserve the earlier implementation but would understate mitochondrial content. Adding a ribosomal cutoff or further gates would add unsupported filtering criteria.

**Rationale**: In the 22,751 complexity-passing cells, complete-mito P95/P97.5/P99 are 16.038/19.313/27.666%; the >20% tail removes 503 cells (2.211%). The cutoff is data-specific, not a scientific finding.

**Consequences**: `sobj_qc_filtered.rds` retains 22,248 of 983,903 cells across S2/S3/S4/S5/S7/S8. An independent test reconstructed the exact saved cell IDs. The previous two-rRNA provenance statement remains historical but is erroneous and superseded.

### [2026-07-12] Select preprocessing source at the pipeline entrypoint

**Tags**: scripts, data-lineage, reproducibility, provenance

**Context**: The pipeline needs to run either the established legacy Seurat object or the Seurat object reconstructed from counts and filtered by QC without changing later clustering, figure, or differential-analysis commands.

**Decision**: Select `legacy` or `counts-qc` with `--input-source` in preprocessing and its wrappers. Keep `--input` for an explicit custom object. Record the selected source in preprocessing metadata; downstream stages consume the regenerated current artifacts.

**Rationale**: Preprocessing is the only shared upstream seam. Selecting there keeps every downstream artifact contract unchanged and prevents parallel source-selection logic from drifting across later scripts.

**Consequences**: A run replaces current branch artifacts. Regenerate preprocessing before clustering and downstream analyses whenever the source changes.

### [2026-07-12] Number executable pipeline stages

**Tags**: scripts, workflow, reproducibility

**Decision**: Prefix each executable stage with its pipeline order. Keep all-branch preprocessing and clustering as numbered companion scripts because they preserve an explicit all-branch interface and clustering's dry-run contract.

**Consequences**: Documentation, recipes, and tripwires use numbered script paths. Retain `R/` helpers that carry substantial, reusable, testable scientific computation; inline the one-call input-source path selector in preprocessing.

### [2026-07-12] Use 20 PCs and Leiden resolution 0.3 for the counts-derived PFlog run
**Run-specific results superseded 2026-07-13 by "Preserve the chosen PFlog parameters for the emptyDrops/log-MAD rebuild":** The 20-PC, resolution-0.3 parameter choice remains current, but the source/MG cluster counts and excluded cluster IDs below describe the prior QC object.

**Tags**: clustering, pflog, mg-selected, reproducibility

**Context**: The counts-derived QC-filtered object was rebuilt through the PFlog preprocessing and clustering grids. The full source candidate grid contained 13, 29, and 37 clusters at 20, 30, and 50 PCs respectively at resolution 0.3; the MG-selected grid contained 19, 36, and 32 clusters.

**Decision**: Use the no-cell-cycle-HVG PFlog candidate with 20 PCs and Leiden resolution 0.3 for source-cluster filtering and downstream MG-selected figures, marker ranking, pseudobulk DE/DD, and enrichment. Retain 30- and 50-PC candidates as sensitivity outputs.

**Rationale**: Twenty PCs are past the elbow and preserve broad separated groups without the additional fragmentation seen at higher PC counts or higher resolutions.

**Consequences**: Source clusters 10 (microglia) and 13 (p27-high) are excluded; no source cluster met the configured photoreceptor exclusion criterion. Notebook paths and cluster references use `dims20_res0.3`.

### [2026-07-12] Canonicalize the E-Stim condition label at count ingestion

**Tags**: metadata, contrast, data-lineage, reproducibility

**Context**: Raw sample metadata uses `p27CKO + EStim`, while the analysis contract declares `p27CKO +EStim`. The mismatch blocked condition-aware figures and pseudobulk DE.

**Decision**: In `scripts/01-process-counts.R`, normalize only the optional whitespace between `+` and `EStim` before creating the Seurat object.

**Consequences**: Counts-derived objects use exactly `p27CKO` and `p27CKO +EStim`, matching `analysis_labels.yml` and the DE contrast.

### [2026-07-13] Use liberal floors before sample-specific MAD and doublet filtering
**Superseded 2026-07-13 by "Use emptyDrops calls and log-MAD QC flags in the rewritten QC stage":** Carl replaced this implementation intentionally. The current script does not run scDblFinder or apply the fixed liberal thresholds to `pass_qc`.

**Tags**: qc-filtering, doublets, mitochondrial, data-lineage, reproducibility

**Context**: The counts-derived QC stage needs a permissive global gate before robust sample-specific quality thresholds and doublet detection. Each `Sample` represents an independently processed input and therefore the physical scope in which a doublet can arise.

**Decision**: First retain cells with `nFeature_RNA >= 50`, `nCount_RNA >= 100`, and `percent.mt <= 20`. Calculate lower three-MAD thresholds for counts and detected features and an upper three-MAD threshold for mitochondrial percentage within each sample, bounded by those global floors and ceiling. Run `scDblFinder` independently by `Sample` after the liberal gate. Save cells that pass all three effective MAD criteria and are classified as singlets. Keep `percent.ribo` diagnostic only; apply no sample, empty-drop, or ambient-RNA filter.

**Alternatives considered**: Saving only the liberal-floor object would preserve every candidate cell but would not apply the requested stricter thresholds. Two-sided MAD limits for counts and features would remove high-complexity cells despite separate doublet classification and the strongly right-skewed count distributions.

**Rationale**: The liberal gate removes the extreme low-information population before robust summaries and doublet modeling. One-sided quality limits match the expected failure direction for low counts/features and high mitochondrial percentage, while sample-specific scDblFinder calls avoid impossible cross-capture doublets.

**Consequences**: The liberal gate retains 22,248 of 983,903 cells; 22,010 pass all effective MAD criteria; 1,566 liberal-floor cells are called doublets; and `sobj_qc_filtered.rds` retains 20,459 MAD-passing singlets across all six samples. The effective count and feature lower limits equal the liberal floors in every sample, while mitochondrial ceilings tighten in S2, S5, S7, and S8.

### [2026-07-13] Use emptyDrops calls and log-MAD QC flags in the rewritten QC stage

**Tags**: qc-filtering, emptydrops, metadata, data-lineage, reproducibility

**Context**: Carl rewrote `scripts/02-qc-filtering.R` to simplify the Seurat metadata and replace the prior scDblFinder-plus-fixed-floor implementation. He confirmed that the omitted scDblFinder step and unused fixed liberal thresholds are intentional.

**Decision**: Run `barcodeRanks()` and `emptyDrops()` on the raw RNA count matrix, store cell-call probability/FDR and `is_cell`, calculate sample-specific count and feature lower limits on the log10 scale and a mitochondrial upper limit among called cells, and define `pass_qc` from those three MAD criteria. Save all cells with the compact QC metadata as `sobj_raw_with_qc.rds` and save the `pass_qc` subset as `sobj_qc_filtered.rds`.

**Alternatives considered**: Restoring scDblFinder, adding the prior fixed-floor gate to `pass_qc`, or preserving the larger earlier metadata set would override intentional analysis choices in the rewrite.

**Rationale**: The rewritten script makes cell calling, threshold estimation, and saved QC flags explicit while retaining a compact handoff. The raw annotated object preserves enough information to reconstruct alternate subsets without carrying every prior diagnostic field.

**Consequences**: Downstream `counts-qc` preprocessing remains path-compatible and receives RNA counts plus `Mouse` and `Condition`. The filtered object is defined by `pass_qc`, not by scDblFinder, the fixed 50-feature/100-count/20%-mitochondrial gate, or `is_cell`; those omissions are intentional.

### [2026-07-13] Preserve the chosen PFlog parameters for the emptyDrops/log-MAD rebuild

**Tags**: clustering, pflog, mg-selected, reproducibility, emptydrops

**Context**: Carl regenerated preprocessing and full clustering from the rewritten emptyDrops/log-MAD QC object, then requested the remaining downstream pipeline. The prior chosen analysis identifiers remain the no-cell-cycle-HVG PFlog branch, 20 PCs, and Leiden resolution 0.3.

**Decision**: Continue the downstream rebuild with `cluster_pflog_no_filter_cc_dims20_res0.3` as the source clustering and `cluster_pflog_mg_selected_no_filter_cc_dims20_res0.3` as the chosen MG clustering. Regenerate 30- and 50-PC MG candidates as sensitivity outputs without changing the chosen identifiers.

**Alternatives considered**: Re-selecting normalization, PC count, or resolution from the new candidate grids would expand the requested pipeline rerun into a new parameter-selection exercise.

**Rationale**: The user asked to continue the established downstream pipeline after rebuilding preprocessing and clustering, not to revise the chosen clustering contract.

**Consequences**: The revised QC object contains 4,145 source cells in 9 chosen source clusters. Source clusters 2, 7, and 8 meet configured exclusion criteria, leaving 3,460 MG-selected cells; the chosen MG clustering contains 7 clusters. Notebook prose and marker notes now use the regenerated cluster IDs.

### [2026-07-13] Plan one canonical run interface while preserving scientific stages

**Tags**: workflow, scripts, planning, reproducibility

**Context**: Reproducing the current counts-derived manuscript analysis requires more than a dozen mixed `Rscript` and `just` commands, manual object-path handoffs, and repeated branch, PC-count, and resolution identifiers. Historical 30- and 50-PC defaults remain exposed beside the chosen 20-PC result.

**Decision**: Plan one deep orchestration module behind `just run`. The module will default to the counts-derived source, retain legacy and explicit-object sources, own stage ordering and generated handoffs, expose a side-effect-free dry run, and require explicit overwrite intent for protected statistical outputs. Existing low-level scripts and recipes remain available for focused reruns.

**Alternatives considered**: Three coarse human-run stage commands would reduce but not eliminate manual sequencing. A general workflow engine or persistent resume system would add machinery that this one-off analysis does not need.

**Rationale**: One small human interface removes repeated operational knowledge without changing the scientific modules or hiding the complete stage plan from dry-run review.

**Consequences**: GitHub issue #1 records the incremental implementation and testing plan. The existing long command sequence remains authoritative until the canonical command completes one verified counts-derived run. QC, clustering choices, MG selection, Mouse × Condition pseudobulk analysis, sensitivity grids, and low-level expert commands remain unchanged.

### [2026-07-13] Keep cell calling and doublet status separate from the saved MAD selector

**Tags**: qc-filtering, emptydrops, scDblFinder, metadata, data-lineage

**Context**: The current QC script combines emptyDrops cell calling, per-sample scDblFinder calls, and sample-specific log-MAD thresholds, but the saved `pass_qc` subset intentionally has a narrower selector than the threshold-estimation population.

**Decision**: Run scDblFinder per sample on emptyDrops-called barcodes above the count and feature floors, estimate sample-specific count, feature, and mitochondrial thresholds among called singlets, and define the saved `pass_qc` subset from those three MAD criteria only. Preserve `is_cell`, `doublet_call`, and `is_singlet` as separate metadata rather than adding them to `pass_qc`.

**Alternatives considered**: Requiring `is_cell` or `is_singlet` in `pass_qc` would change the user-authored selector. Omitting scDblFinder from the documented threshold-estimation path would misdescribe the current script.

**Rationale**: Separate flags preserve the current analysis choice and allow later sensitivity subsets without hiding how the MAD thresholds were estimated.

**Consequences**: The counts-derived source contains 4,146 `pass_qc` cells. Documentation now states that called singlets define the MAD thresholds while cell-call and singlet flags remain separate from `pass_qc`. This supersedes the earlier description of a QC stage without scDblFinder.

### [2026-07-13] Select source 30 PCs and MG-selected 20 PCs at resolution 0.5

**Tags**: clustering, pflog, mg-selected, stability, reproducibility

**Context**: Regenerating the post-QC source and MG-selected 20/30/50-PC × 0.3/0.5/0.8 grids changed the candidate structure relative to the earlier 20-PC/resolution-0.3 choices.

**Decision**: Use the PFlog no-cell-cycle-HVG source candidate at 30 PCs and resolution 0.3, then use the PFlog no-cell-cycle-HVG MG-selected candidate at 20 PCs and resolution 0.5. Continue computing 50 MG PCs and retaining all sensitivity candidates.

**Alternatives considered**: Keeping the prior 20-PC/resolution-0.3 identifiers would ignore the regenerated stability summaries. Choosing the 50-PC MG candidate would retain more tail variance without a compensating stability or sample-support advantage.

**Rationale**: The source candidate has the strongest local stability in the regenerated source grid. The MG candidate balances local stability, cluster count, and sample support; all eight clusters span both conditions and seven span every mouse.

**Consequences**: The source has 4,146 cells in 9 clusters. Excluding source clusters 2, 7, and 8 retains 3,456 cells; the chosen MG clustering has 8 clusters. Source figures use dims30/resolution 0.3 and MG figures, markers, and DE use dims20/resolution 0.5. This supersedes the preserved 20-PC/resolution-0.3 choices recorded earlier on 2026-07-13.

### [2026-07-13] Use one canonical just interface for complete runs

**Tags**: workflow, just, orchestration, reproducibility, validation

**Context**: The manuscript pipeline required a long sequence of low-level commands with repeated paths and clustering identifiers.

**Decision**: Make `just run [source] [overwrite]` the routine interface and `just run-dry-run [source] [overwrite]` its side-effect-free plan. Default to `counts-qc` and `overwrite=false`; support `legacy` and explicit RDS paths. Keep low-level recipes for checkpoint recovery.

**Alternatives considered**: Three coarse stage commands would still expose handoff state. A workflow engine, cache, resume database, or parallel scheduler would add machinery that this one-off analysis does not need.

**Rationale**: The orchestration module owns stage order, paths, chosen identifiers, output validation, report rendering, and QA behind two human inputs.

**Consequences**: A complete counts-derived run executes 24 validated stages through notebook rendering and tripwires. Legacy and explicit-object plans contain 22 stages because they skip count ingestion and QC. Protected statistical outputs require explicit overwrite intent.

### [2026-07-14] Remove differential detection and report primary DE with a volcano plot

**Tags**: differential-expression, mg-selected, plotting, reproducibility

**Context**: The MG-selected condition analysis combined DESeq2 differential expression with a muscat differential-detection workflow and a joint effect-size scatter. The differential-detection branch added a separate gene universe and output contract without contributing to the retained condition-response interpretation.

**Decision**: Remove differential detection completely. Keep the six-sample primary DESeq2 model, the paired-mouse DE sensitivity, marker overlap, and enrichment. Report the primary model with a volcano plot using shrunken log2 fold change and adjusted P value, and label the ten most significant genes deterministically.

**Alternatives considered**: Keeping DD as a sensitivity analysis would preserve extra outputs but retain the dependency and mixed gene-universe interpretation. Keeping the joint scatter without DD would not provide a standard significance view.

**Rationale**: A DE-only contract is simpler and matches the retained scientific claim. The primary-model volcano communicates effect direction and FDR evidence directly, while paired DE remains available as a sensitivity table.

**Consequences**: `muscat`, `pbDS`, detection tables, DD report values, and the DE/DD scatter are removed from active code, contracts, documentation, and artifacts. `mg_selected_de_volcano.(png|pdf)` and its notebook link replace the scatter. This decision supersedes earlier active DD implementation and plotting decisions while preserving their historical record.

### [2026-07-14] Permute condition labels only after preserving sample identity

**Tags**: tripwires, label-permutation, preprocessing, reproducibility

**Context**: `scripts/03-preprocess.R` validates `Mouse` and `Condition`, then derives `sample_id` from both fields before blind HVG and PCA work. Permuting `Condition` in the input object would also change the derived identity and could test metadata inconsistency rather than label leakage.

**Decision**: Restrict condition permutation to explicit tripwire mode after `sample_id` derivation. Permute globally across intact sample units while preserving `sample_id`, `Mouse`, cells, counts, and all other metadata. Require a checkpoint log and `STOP_AFTER_CHECKPOINT=blind_qc_complete`; skip pre-checkpoint figure writes and compare exact retained-HVG and PCA-standard-deviation fingerprints.

**Alternatives considered**: Permuting the input metadata before preprocessing would conflate label leakage with sample-identity validation. Cell-level permutation would destroy sample structure. Allowing the seed environment variable without a test-mode guard could silently alter a production run.

**Rationale**: The metamorphic test should change only the interpretation label while leaving every blind-analysis input invariant. The strict guard and stop boundary prevent permuted labels or scratch diagnostics from reaching production outputs.

**Consequences**: `label-permutation` now executes two guarded scratch runs and passes only when the blind fingerprints match and production paths remain unchanged. Missing metadata reconciliation and provenance semantics also run as executable tripwires.

### [2026-07-14] Balance volcano labels by effect direction

**Tags**: differential-expression, plotting, mg-selected

**Context**: Ranking one pooled set of volcano labels by adjusted P value favored the stronger negative-expression signal and did not summarize both response directions evenly.

**Decision**: Label up to 20 genes with positive shrunken log2 fold change and 20 genes with negative shrunken log2 fold change. Rank each direction independently by raw DESeq2 P value, then by absolute effect size and gene name for deterministic ties. Keep the y-axis and significance colors based on adjusted P values.

**Rationale**: Separate ranking guarantees representation of both directions without changing the FDR-based visual evidence.

**Consequences**: The primary-model volcano now carries at most 40 labels, balanced by fold-change direction; notebook prose distinguishes raw-P-value label ranking from adjusted-P-value plotting and significance.

### [2026-07-14] Use apeglm shrinkage and simplified GO summaries

**Tags**: differential-expression, enrichment, plotting, mg-selected, reproducibility

**Context**: The MG-selected DE stage defaulted to normal-prior LFC shrinkage and wrote only raw GO ORA/GSEA tables. The notebook lacked a dedicated enrichment section and a direct Ascl1/Hes6 coexpression view.

**Decision**: Default DESeq2 LFC shrinkage to `apeglm`. Filter GO results to FDR-significant terms before `clusterProfiler::simplify()`, and compare ORA results with deterministic `enrichit::bayes_enrich()` selection before simplifying. Keep these helpers inline in the executable scripts because they serve one analysis stage.

**Rationale**: `apeglm` is the requested shrinkage model; simplifying only enriched terms preserves the inferential scope and avoids unnecessary semantic-similarity work. Inline helpers keep the one-off analysis flow visible.

**Consequences**: The pipeline now protects five simplified tables, five PNG/PDF dotplots, their notebook links, and both branch-specific Ascl1/Hes6 scatter outputs. Notebook section 4 reports the enrichment results.
### [2026-07-14] Lint first-party R scopes, not inert skillpacks

**Tags**: scilintr, lint, tooling, skillpacks

**Context**: The scilintr integration checks ESPI's first-party R analysis scopes. `skillpacks/` contains inert external reference material and is not part of the executable project workflow.

**Decision**: Make `just lint` cover `R/`, `scripts/`, `data-raw/`, `tools/`, `notebook/sc_analysis.qmd`, and `config.local.example.R`, while excluding `skillpacks/`.

**Alternatives considered**: Linting every R-like file would include copied external skill content and produce findings unrelated to ESPI's analysis code.

**Rationale**: The lint boundary should match code owned and executed by this project, not inert vendored references.

**Consequences**: Findings in covered scopes must be fixed or recorded with a structured `ANALYSIS_OK` waiver; excluded skillpack content does not affect the project lint gate.

### [2026-07-15] Freeze the pipeline through MG clustering

**Decision**: Routine `just run` starts from the existing clustered MG-selected RDS objects. Scripts `01` through `07`, every `04-cluster.R` call, source summaries, and script `06` marker heatmaps run only through the explicit `just regenerate-frozen` interface or deliberate low-level recovery.

**Rationale**: The upstream cell selection and clustering choices now define the fixed analysis cohort. Routine downstream work must not silently regenerate those artifacts.

**Consequences**: `just run [overwrite]` begins at the MG cluster summary. Full regeneration requires an explicit recipe name and source argument.
### [2026-07-15] Cluster marker heatmap only at the cluster-mean level

**Tags**: plotting, heatmap, clustering

**Decision**: Preserve input cell order within each cluster. Cluster the scaled-expression means with `hclust(dist(t(cluster_means)))`, use the resulting dendrogram leaf order directly for the column-split factor, and draw only that dendrogram above the heatmap.

**Rationale**: ComplexHeatmap couples its parent slice dendrogram to within-slice clustering; disabling cell clustering also removes the parent dendrogram. Computing and drawing the parent dendrogram separately keeps the intended cluster relationships without implying a meaningful ordering among individual cells.

**Consequences**: The marker heatmap no longer shows per-cell dendrogram branches or the dotted parent/child boundary. Cluster labels remain in the colored block annotation.


### [2026-07-16] Freeze all existing Seurat objects and designate two final objects

**Tags**: seurat, reproducibility, frozen, clustering

**Context**: The publication analysis now has fixed source and MG-selected clustering choices. Existing alternative and intermediate Seurat objects remain useful as provenance and sensitivity artifacts but must not be regenerated or silently replaced.

**Decision**: Freeze all 18 existing RDS files under the external `seurat_objects/` root without regeneration. Designate `current/cluster_pflog_no_filter_cc_elbow20.rds` with `cluster_pflog_no_filter_cc_dims30_res0.3` as the final source object and `current/cluster_pflog_mg_selected_no_filter_cc_elbow20.rds` with `cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5` as the final MG-selected object. Here, `no_filter_cc` means cell-cycle genes remain eligible as HVGs.

**Rationale**: Fixed object bytes and cluster columns make every downstream publication analysis conditional on one explicit cell cohort and clustering definition.

**Consequences**: `FROZEN_OBJECTS.tsv` records SHA-256 checksums, byte sizes, and roles. All listed RDS files, the manifest, and the `seurat_objects/`, `input/`, and `current/` directories are read-only. Deliberate regeneration now requires explicitly restoring write permissions before using the existing regeneration interface.

### [2026-07-16] Define boundaries for the publication-pipeline restructuring proposal

**Tags**: pipeline, planning, scripts, reproducibility

**Context**: The current analysis framework supports more inputs, stage-level options, validation, and orchestration paths than the final publication analysis needs.

**Draft proposal pending approval**: Use a fixed entire-pipeline cutover with four scientific phase scripts and a minimal R package. Keep one explicit frozen-stage regeneration script; consolidate downstream figures, marker analysis, and DE/enrichment into three fixed scripts. Preserve exact tabular results, visually equivalent figures, and the visible notebook content. Keep only checks that protect scientific meaning. Leave existing artifacts in place.

**Rationale**: A computational biologist should be able to read the scripts in order without learning a custom execution framework, while the frozen objects and publication results remain unchanged.

**Status and consequences**: This is not an adopted architecture. If approved, implement incrementally beside the existing code, prove equivalence phase by phase, then delete replaced scripts, shallow helpers, broad tripwires, compatibility paths, and stage-level CLI options. The proposal is documented in `architecture-restructure-proposal.html`; no restructuring has been implemented.

### [2026-07-16] Adopt the clean four-phase publication cutover

**Tags**: pipeline, architecture, reproducibility

**Context**: The approved publication workflow needs a readable fixed interface
without dynamic stage selection or framework-only validation.

**Decision**: Use exactly four phase scripts (`01-regenerate-frozen.R`,
`02-publication-figures.R`, `03-marker-analysis.R`, and `04-de-enrichment.R`)
and four focused R modules. The five analysis commands are `just run`,
`just figures`, `just markers`, `just de`, and deliberate
`just regenerate-frozen`; `load`, `document`, `readme`, `format`, and `lint`
remain maintenance recipes. Phase 01 refuses writable frozen-object
directories, and downstream phases use fixed contracts plus one overwrite
flag.

**Rationale**: A computational biologist can follow the scientific work from
frozen inputs to publication artifacts without learning an orchestration
framework, while the existing analysis remains the equivalence oracle.

**Consequences**: Legacy execution paths and broad instrumentation are not
active interfaces. Notebook figure mirrors are regular files updated through
temporary copies with hash and dimension checks.

### [2026-07-16] Resolve the source-versus-visible MG heatmap oracle conflict

**Tags**: plotting, heatmap, sensitivity, notebook

**Context**: The visible notebook contains both the final MG heatmap and the
CC-filtered MG sensitivity heatmap, while the designated final MG object is
the source for marker and DE analyses.

**Decision**: Phase 02 loads the final source, final MG-selected, and
CC-filtered MG sensitivity objects once each and preserves both MG heatmap
branches. The notebook's ordered visible figures, dimensions, captions, and
artifact names are the acceptance oracle. Phase 03 marker results remain
descriptive and phase 04 independently rebuilds curated marker overlap.

**Rationale**: Treating the source object as the only heatmap oracle would
silently drop a visible sensitivity result; treating marker output as a DE
input would change the curated scientific method.

**Consequences**: Equivalence review checks both MG heatmap branches and the
rendered notebook sequence. Safe regular-file mirrors may change encoding
bytes but must preserve content, dimensions, and names.
