# Learnings

Append-only log of gotchas, surprises, and reusable workflow lessons.

### [2026-07-03] Mycelium hook paths are local plugin-cache paths

**Tags**: mycelium, hooks, setup, portability

**Category**: Tooling setup

**What happened**: Mycelium generated `.claude/settings.local.json` with hook commands pointing at `/Users/carlstone/.omp/plugins/cache/plugins/mycelium___mycelium___0.0.0/`.

**Why it matters**: These paths can break after plugin upgrades, cache cleanup, or a fresh clone on another machine.

**Resolution**: Keep `.claude/` gitignored and document the need to rerun Mycelium hook setup after plugin changes.

**mitigation_type**: documentation

**structural_mitigation_candidate**: Prefer stable project-local hook wrappers if Mycelium supports them later.

### [2026-07-03] Quarto embedded HTML must be rerendered after figure regeneration

**Tags**: quarto, figures, box, notebook

**Category**: Reproducible reporting

**What happened**: `notebook/sc_analysis.qmd` uses `embed-resources: true`, so regenerated Box figures do not update `notebook/sc_analysis.html` until Quarto rerenders the document.

**Why it matters**: Symlink targets can be current while the HTML still shows stale embedded image bytes.

**Resolution**: Run `quarto render notebook/sc_analysis.qmd` after regenerating source figures.

**mitigation_type**: workflow

**structural_mitigation_candidate**: Add a notebook render check to any future figure-regeneration script.

### [2026-07-03] Documented Autonomous-Science skillpack URL is unavailable

**Tags**: mycelium, skillpacks, github, setup

**Category**: Tooling setup

**What happened**: `git clone https://github.com/arjunrajlaboratory/Autonomous-Science.git` failed with `Repository not found`, and GitHub/web searches found no matching public repository.

**Why it matters**: The Mycelium skillpacks README references a third external skill library that cannot currently be installed from the documented URL.

**Resolution**: Clone the two available libraries, install `skill-bridge`, and document the missing repository in `skillpacks/README.md`.

**mitigation_type**: documentation

**structural_mitigation_candidate**: Update the Mycelium plugin's skillpack template if a corrected `Autonomous-Science` URL becomes available.

### [2026-07-03] Seurat rewrites hyphens in reduction names

**Tags**: seurat, clustering, reductions, naming

**Category**: Pipeline robustness

**What happened**: `RunUMAP(reduction.name = "umap_pflog_no-filter-cc_dims20")` stored the reduction as `umap_pflog_no.filter.cc_dims20`, so a later lookup for the original hyphenated name failed.

**Why it matters**: Cluster artifact names that look valid as file names may not be valid as stable Seurat object keys.

**Resolution**: Use underscore-only branch tags for Seurat reductions and metadata columns, and keep a tripwire that requires branch tags to match `[A-Za-z0-9_]+`.

**mitigation_type**: code-and-tripwire

**structural_mitigation_candidate**: Centralize branch-tag construction if more scripts start creating Seurat object keys.

### [2026-07-03] Pass clustree only cluster columns for large Seurat metadata

**Tags**: clustree, seurat, plotting, warnings

**Category**: Plotting robustness

**What happened**: Passing the full Seurat object into `clustree()` emitted repeated warnings about unrelated metadata column name sanitization.

**Why it matters**: Those warnings obscure real plotting problems and repeat once per panel in multi-panel clustree figures.

**Resolution**: Build clustree plots from a metadata data frame containing only the cluster columns matching the requested prefix.

**mitigation_type**: code

**structural_mitigation_candidate**: Keep plotting helpers narrow: pass only the columns a plotting package needs when the source object has large heterogeneous metadata.

### [2026-07-03] Use ASCII keys plus labels for marker sets

**Tags**: marker-genes, unicode, plotting, r-package

**Category**: Data object design

**What happened**: The marker set includes Müller glia, whose display name contains a non-ASCII character and a space.

**Why it matters**: R can store quoted Unicode list names, but those names become awkward when reused as column names, file tags, factor keys, or programmatic accessors.

**Resolution**: Store marker lists under ASCII keys such as `muller_glia`, and store plot/report names separately in `cell_type_marker_labels`.

**mitigation_type**: data-design

**structural_mitigation_candidate**: Validate that marker-list keys and label-vector names match whenever marker data is regenerated.

### [2026-07-03] Use unlink for notebook figure symlink replacement

**Tags**: r, symlink, notebook, figures

**Category**: Figure output workflow

**What happened**: Replacing an existing notebook figure symlink with `file.remove()` produced a warning when the path was a broken symlink.

**Why it matters**: Figure scripts should be rerunnable without noisy warnings or stale notebook figure links.

**Resolution**: Detect either an existing file or symlink with `file.exists(path) || nzchar(Sys.readlink(path))`, then remove the path with `unlink(path)` before recreating the symlink.

**mitigation_type**: structural

**structural_mitigation_candidate**: Use this symlink replacement idiom in future scripts that link Box figures into `notebook/figures/`.

### [2026-07-04] Read TSV row counts with a count tool when output looks truncated

**Tags**: read-tool, tsv, validation

**Category**: Tooling validation

**What happened**: A line-selected read of `detection_full_results.tsv` appeared to show only 44 lines even though `wc -l` and R both showed 50,765 file lines.

**Why it matters**: Large TSV previews can mislead downstream interpretation if the apparent displayed range is treated as a row count.

**Resolution**: Use `wc -l` or R `nrow(readr::read_tsv(...))` for row-count verification when a preview seems inconsistent with the expected analysis scale.

**mitigation_type**: validation

**structural_mitigation_candidate**: Keep reportable output counts in `numbers.json` and re-read that file directly before summarizing final DE results.

### [2026-07-04] Make DE output overwrites explicit

**Tags**: differential-expression, reproducibility, output-provenance

**Category**: Analysis output safety

**What happened**: The MG-selected DE script wrote canonical output filenames while also accepting alternate inputs and cluster columns.

**Why it matters**: Exploratory reruns can silently replace manuscript-facing DE, detection, enrichment, and reportable-value files if overwrite behavior is implicit.

**Resolution**: Require `--overwrite` before replacing existing DE/enrichment outputs, record the input object, cluster column, counts layer, shrinkage method, and GSEA seed in `design_summary.tsv` and `numbers.json`, and keep a GSEA symbol-to-Entrez mapping ledger.

**mitigation_type**: structural

**structural_mitigation_candidate**: For future analysis scripts with configurable inputs and canonical output paths, add an explicit overwrite gate and write run-provenance fields beside the main results.

### [2026-07-04] Completed LOG_REGISTRY rows still need semantic fields

**Tags**: mycelium, session-logs, hooks, reproducibility

**Category**: Tooling setup

**What happened**: `.living/log/LOG_REGISTRY.md` had a completed `2026-07-04-001` row with empty Summary and Key Outputs after the deterministic stop-hook path finalized the row without log-scribe enhancement.

**Why it matters**: The registry is the quick cross-session map, so filename stubs or empty semantic fields make it harder to resume work and audit what a session accomplished.

**Resolution**: Treat `LOG_REGISTRY.md` as a per-session registry, manually enhance completed rows when automatic log-scribe output is missing, and add a matching `## Session Summary` section to the linked log.

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: Add a validation check that flags completed registry rows with empty Summary or Key Outputs.

### [2026-07-04] Clear stale Mycelium reminder files after false-positive stop blocks

**Tags**: mycelium, hooks, session-state, reproducibility

**Category**: Tooling setup

**What happened**: The stop hook repeatedly reported `.living/ not updated` even after the review session was triaged and logged. The remaining `.claude/mycelium-reminded.tmp` timestamp was newer than `.living/learnings.md`, while later updates touched only `.living/log/` and `.claude/last-session.md`, which the stop hook does not count as satisfying triage.

**Why it matters**: A stale reminder sentinel can create a false-positive stop block loop and push agents toward padding `.living/` with unnecessary entries.

**Resolution**: Confirm the requested triage is complete, then remove the stale `.claude/mycelium-reminded.tmp` sentinel instead of adding duplicate decisions or learnings. Prevent recurrence by keeping Mycelium maintenance commands from being forwarded from the OMP adapter to the synced post-action hook.

**mitigation_type**: structural

**structural_mitigation_candidate**: Keep the OMP adapter guard that skips post-action forwarding for `skills/core/scripts/{generate_index,validate_structure,recall_lessons,...}.py` and `tools/sync-mycelium-skills-core.py`, or move the same exclusion upstream into Mycelium.

### [2026-07-04] OMP hook adapters must return modified tool content

**Tags**: mycelium, hooks, omp, data-lineage, r

**Category**: Tooling setup

**What happened**: The Mycelium OMP adapter invoked the post-action wrapper but dropped the wrapper's JSON `additionalContext`, so real analysis commands could create sentinels without surfacing the mandatory Mycelium post-action protocol. The synced data-lineage extractor also detected `Rscript` commands but only scanned Python-style I/O calls.

**Why it matters**: Silent hooks make the framework look wired while key user-visible behavior is missing. In ESPI, Python-only lineage regexes make R pipeline lineage effectively empty.

**Resolution**: Return appended tool content from the OMP `tool_result` handler, call the data tracker and lineage stop hook from OMP, and use a repo-owned data-tracker wrapper with R I/O expression detection instead of editing synced `skills/core/` files.

**mitigation_type**: structural

**structural_mitigation_candidate**: Keep Mycelium customizations in repo-owned adapters or wrappers, then verify them with synthetic hook JSON after any OMP or Mycelium plugin update.

### [2026-07-04] Stop hooks need session-boundary sentinel checks

**Tags**: mycelium, hooks, session-state, omp, tests

**Category**: Tooling setup

**What happened**: False-positive stop blocks occurred at the end of early turns because `mycelium-health.sh` only deleted stale reminder/activity files after one hour. A new OMP session that started within one hour of the previous session inherited the previous `.claude/mycelium-reminded.tmp`, and `mycelium-stop-check.sh` treated the old timestamp as current-session work.

**Why it matters**: OMP can call session-stop behavior at turn boundaries, so cross-session sentinel bleed creates blocks before the agent has had a chance to do meaningful session-end triage.

**Resolution**: Reset reminder/activity sentinels when they predate the current session-start timestamp, and add stop-hook regression tests for prior-session reminder and activity files that are less than one hour old.

**mitigation_type**: structural

**structural_mitigation_candidate**: Keep `test_stop_hook.sh` cases 14b/14c and `test_hooks_stress.sh` case 27 as regression coverage for cross-session sentinel bleed.

### [2026-07-04] Differential-detection gene filters can change empirical-Bayes behavior

**Tags**: differential-detection, limma, filtering, mg-selected

**Category**: Analysis robustness

**What happened**: Filtering MG-selected differential detection to the same gene universe as DE left primary DD with zero FDR-significant genes, but changed paired-sensitivity DD from zero to 108 FDR-significant genes.

**Why it matters**: limma's empirical-Bayes moderation depends on the tested gene set. In the paired sensitivity model, only four samples and one residual degree of freedom remain, so changes in the gene universe can strongly affect moderated detection p-values.

**Resolution**: Keep the matched DE/DD gene universe for comparability, but treat paired-sensitivity DD hits as exploratory unless they survive an additional sample-level DD method or a larger paired design. Marker-list overlap remained zero after the filter change.

**mitigation_type**: structural

**structural_mitigation_candidate**: Keep DD tested-gene counts and primary/paired DD hit counts in `numbers.json`, and add a targeted DD sensitivity using a sample-level count model before making strong detection-fraction claims.

### [2026-07-04] Muscat DD changes the paired-sensitivity hit set

**Tags**: differential-detection, muscat, filtering, mg-selected

**Category**: Analysis robustness

**What happened**: Replacing limma empirical-logit DD with muscat `edgeR_NB_optim` kept primary DD at zero FDR-significant genes but changed paired-sensitivity DD from the prior limma 108 hits to 40 hits. The muscat-native tested-gene universe also changed from the prior matched-DE 22,663 paired genes to 34,880 paired DD genes.

**Why it matters**: The approved no-extra-prefilter muscat workflow lets the internal 90%-detection filter define the DD universe. In sparse PipSeq data, that filter can keep more genes than the DE count floor, so plan expectations that DD tested genes must be fewer than DE tested genes can be false.

**Resolution**: Keep the muscat-native DD universe unless the analysis plan explicitly chooses an additional prefilter. Report primary and paired DD tested-gene counts and hit counts from `numbers.json`; treat paired-sensitivity DD as exploratory because it still uses only four samples.

**mitigation_type**: analysis-provenance

**structural_mitigation_candidate**: Record the DD method, tested-gene counts, and hit counts in `numbers.json` and `design_summary.tsv` whenever the DD method or gene universe changes.

### [2026-07-04] Do not report planned marker analyses as completed

**Tags**: reporting, marker-analysis, notebook, review

**Category**: Reproducible reporting

**What happened**: Static review found that `notebook/sc_analysis.qmd` introduced a `FindAllMarkers()` cluster-marker section saying the analysis was used, while `scripts/find-markers-mg-selected.R` was still only an untracked comment scaffold with no marker table or dot-plot outputs.

**Why it matters**: In ESPI, notebook prose is treated as the scientific narrative. Claiming an ungenerated marker analysis can make planned exploratory work look like completed evidence, and `FindAllMarkers()` p-values would also need explicit descriptive/non-confirmatory framing because the same cells define the clusters and the marker tests.

**Resolution**: Remove or mark planned marker prose until the executable script and output artifacts exist. When adding cluster markers, describe `FindAllMarkers()` as descriptive marker ranking unless a sample-aware or held-out validation is added.

**mitigation_type**: report-freshness

**structural_mitigation_candidate**: Extend `tools/run-tripwires.R` so report-number/artifact freshness checks compare notebook claims directly against generated outputs and `numbers.json`, without hard-coding expected biological counts in the runner.

### [2026-07-04] Seurat FoldChange defaults are invalid for PFlog marker ranking

**Tags**: seurat, findallmarkers, pflog, marker-analysis, mg-selected

**Category**: Analysis robustness

**What happened**: Running `Seurat::FindAllMarkers()` directly on the PFlog layer produced sign-incoherent positive markers: several top-ranked genes had `pct.1 < pct.2`, and the dot plot showed those genes lower in their assigned cluster than in other clusters.

**Why it matters**: Seurat's default fold-change calculation assumes log-normalized expression and can silently turn PFlog-scale anti-markers into positive `avg_log2FC` markers. That would mislabel clusters and overstate cluster identity evidence.

**Resolution**: Run the marker workflow on Seurat's standard `data` layer for fold-change ranking while keeping the PFlog-derived clustering, and enforce `pct.1 > pct.2` for every retained positive marker row. The marker script now drops identities with no positive detection-enriched markers instead of forcing anti-markers into the top table.

**mitigation_type**: structural

**structural_mitigation_candidate**: Keep the `pct_diff > --min-diff-pct` assertion in `scripts/find-markers-mg-selected.R` and add the same invariant to any future `FindAllMarkers()` workflow that uses non-standard assay layers.

### [2026-07-04] DE and DD effect plots need explicit gene joins

**Tags**: differential-expression, differential-detection, plotting, mg-selected

**Category**: Reproducible reporting

**What happened**: The MG-selected DESeq2 and muscat DD workflows test different gene universes: primary DE tested 24,514 genes, while primary DD tested 36,468 genes. A DE-vs-DD effect-size scatter cannot assume row alignment or full overlap.

**Why it matters**: Joining by row order would silently pair effects from different genes. Outer-joining without a clear policy would place genes missing one effect on an arbitrary or NA axis.

**Resolution**: Build DE/DD scatter plot data with an explicit inner join on `gene` for each design, and state in the caption/prose that only genes tested in both workflows are plotted.

**mitigation_type**: structural

**structural_mitigation_candidate**: Keep the join policy in the plotting code and include a count of plotted shared genes in future `numbers.json` report fragments if this scatter becomes a manuscript figure.

### [2026-07-05] Use Rscript files for heavy ESPI smoke tests

**Tags**: r, validation, mg-selected, tooling

**Category**: Workflow validation

**What happened**: The MCP R runner timed out while loading the package and reading the MG-selected Seurat RDS object, and then also timed out on a later lightweight parse request in the same session.

**Why it matters**: ESPI smoke tests that call `devtools::load_all()` and touch large Seurat objects can exceed short interactive tool timeouts. Retrying the same MCP call can waste time or leave uncertainty about whether the session is still busy.

**Resolution**: For heavy validation, write a temporary R script and run it with `Rscript <file>` through the Bash tool with an explicit timeout. Keep `Rscript -e` out of Bash unless a plan explicitly requires it and no safer evaluator works.

**mitigation_type**: workflow

**structural_mitigation_candidate**: Keep reusable smoke-test scripts in `tools/` only if the same validation becomes routine across multiple analysis changes.

### [2026-07-05] Mycelium hook summaries can overwrite manual semantic log rows

**Tags**: mycelium, hooks, session-logs, provenance

**Category**: Workflow provenance

**What happened**: During Batch 1 review cleanup, manual semantic edits to `.living/log/LOG_REGISTRY.md`, `.living/log/2026-07-05-004-espi.md`, and `.claude/last-session.md` were replaced by hook-generated file-list placeholders, leaving completed records with blank Key Outputs and Tags.

**Why it matters**: The session registry is the audit trail for analysis and review work. If hook-generated placeholders overwrite semantic summaries after manual cleanup, future agents can falsely treat incomplete provenance as complete and reviewers must rediscover the same issue.

**Resolution**: Re-check the current session registry row and last-session summary after any Mycelium hook activity, then repair Summary, Key Outputs, Tags, duration, and file counts before committing or yielding.

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: Add a stop-hook invariant that rejects completed `LOG_REGISTRY.md` rows with blank Summary, blank Key Outputs, blank Tags, or file-list-only summaries.

### [2026-07-05] Persist hook-provenance guard tests with the guard

**Tags**: mycelium, hooks, provenance, tests, session-logs

**Category**: Workflow provenance

**What happened**: A throwaway eval fixture caught that the initial provenance guard misclassified sentence summaries ending in a period as file lists, and a later review caught missing persisted coverage for INDEX restoration.

**Why it matters**: Hook-provenance fixes protect the audit trail itself. If their tests live only in scrollback, a future edit can silently reintroduce registry or INDEX clobbering.

**Resolution**: Keep `tools/test_mycelium_provenance_guard.py` next to `tools/mycelium-provenance-guard.py` and run `python3 -m unittest tools/test_mycelium_provenance_guard.py` after changing the guard or wrappers.

**mitigation_type**: structural

**structural_mitigation_candidate**: Keep provenance guard behavior covered by persisted `tools/` tests whenever the guard learns a new clobbering pattern.

### [2026-07-05] Review sessions can still expose Mycelium provenance clobbering

**Tags**: mycelium, hooks, provenance, review, session-logs, last-session

**Category**: Workflow provenance

**What happened**: The static review of the MG module/p27 heatmap working tree found that `.living/log/LOG_REGISTRY.md` and the linked session log had again been overwritten with file-list-only hook summaries after a semantic log row had been written; the later stop hook also rewrote the current session log and `.claude/last-session.md` after validation passed.

**Why it matters**: The provenance guard and prior awareness reduce this failure mode but do not remove the need to inspect completed log rows and last-session summaries before committing. A final valid verification pass can still be followed by provenance clobbering that loses decisions, outputs, and validation evidence.

**Resolution**: Treat completed `LOG_REGISTRY.md` rows with blank Key Outputs or Tags as a blocking review finding, then re-run session-end triage after stop-hook activity: repair the current log and registry row, regenerate `.living/INDEX.md`, and rewrite `.claude/last-session.md` until the hook path preserves semantic summaries reliably.

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: Extend `tools/mycelium-provenance-guard.py` or its tests to cover this exact file-list-only overwrite after review/session activity, including the current session log and `.claude/last-session.md`.

### [2026-07-05] Treat git status as authority for todo-only stop-hook triage

**Tags**: mycelium, hooks, todo, session-logs, provenance

**Category**: Workflow provenance

**What happened**: After commit `0a4e88e`, the stop hook reported stale prior heatmap files even though `git status --short` showed only the TODO registry and two new TODO files before session-end repair. Later in the same implementation thread, the stop hook appended a generic file-list tail to `2026-07-05-011-espi.md` and rewrote `.claude/last-session.md` after the semantic session-end record had already been written.

**Why it matters**: Todo-only turns and follow-on implementation turns should not inherit scientific or code-change summaries from prior sessions. The existing hook-summary clobbering lessons still apply, but these recurrences add post-commit, todo-only, and post-verification no-new-code cases.

**Resolution**: Use actual `git status --short` as the authoritative changed-file set for session-end provenance; after stop-hook activity, re-read the current log and last-session file, then remove stale file-list appendages and restore semantic summaries.

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: Extend the provenance guard to reject stop-hook summaries whose changed-file list is not a subset of the current Git-visible working tree plus expected local Mycelium files.

### [2026-07-06] Use direct marker vectors for curated-only plot filters

**Tags**: marker-genes, plotting, differential-expression, differential-detection, mg-selected

**Category**: Plotting robustness

**What happened**: The DE/DD scatter needed labels limited to curated cell-type marker genes that were significant in DE or DD. `make_marker_table()` looked tempting because it wraps `cell_type_marker_genes`, but it also adds standalone `Cdkn1b` for marker-overlap reports.

**Why it matters**: Plot filters that claim to use the curated cell-type marker list should not silently inherit report-specific standalone genes or other marker-table decorations.

**Resolution**: For curated-only label filters, derive the set with `unique(unlist(cell_type_marker_genes, use.names = FALSE))`, then intersect with the plot-specific significance rule.

**mitigation_type**: workflow

**structural_mitigation_candidate**: If more plots need marker-list membership checks, add a small helper that explicitly distinguishes direct `cell_type_marker_genes` membership from marker-overlap tables that include standalone genes.

### [2026-07-06] Wide condition legends fit better inside UMAP panels

**Tags**: plotting, notebook, umap, mg-selected, condition

**Category**: Plotting presentation

**What happened**: The first condition-colored MG-selected UMAP used the default outside legend, and visual inspection showed the rendered title/canvas was clipped on the left. The condition legend also occupied a large amount of horizontal space because the condition labels are much wider than cluster-number legend entries.

**Why it matters**: Condition UMAPs should show the same cell geometry clearly while keeping the legend close to the color encoding. Wide external legends can waste notebook space, and rendered figures need visual checks because plot layout issues may not appear as R errors or warnings.

**Resolution**: For condition UMAPs, keep `Seurat::DimPlot()` but move the legend inside empty plot space; use the top right for the cell-cycle-retained branch and the bottom right for the cell-cycle-filtered branch.

**mitigation_type**: workflow

**structural_mitigation_candidate**: If more UMAP overlays with long labels are added, consider a small branch-aware legend-placement helper instead of repeating `theme(legend.position = ...)` logic.

### [2026-07-06] DimPlot point stroke changes apparent UMAP point size

**Tags**: plotting, umap, seurat, mg-selected

**Category**: Plotting presentation

**What happened**: MG-selected condition UMAPs drawn with `Seurat::DimPlot()` looked heavier than explicit ggplot UMAPs even when the nominal point size matched.

**Why it matters**: `DimPlot()` defaults to a circle shape with a stroke, while explicit ggplot UMAP layers commonly use shape 19, a filled circle without stroke. The stroke increases the apparent point size and can make dense UMAPs look overplotted.

**Resolution**: When matching apparent UMAP point size between `DimPlot()` and explicit ggplot layers, set the `DimPlot()` point stroke to `0` in addition to matching `pt.size`.

**mitigation_type**: workflow

**structural_mitigation_candidate**: Keep stroke handling explicit in future Seurat UMAP plotting helpers that need to match custom ggplot UMAP layers.

### [2026-07-09] APFS case-insensitive paths can hide justfile casing mistakes

**Tags**: tooling, just, macos, reproducibility

**Category**: Tooling setup

**What happened**: During `just` setup, `read Justfile` and `read justfile` resolved to the same underlying file on the default case-insensitive APFS filesystem, so treating them as duplicates risked deleting the only command file.

**Why it matters**: Casing-only filename checks on macOS can report misleading paths, and cleanup decisions can remove the real file rather than a duplicate.

**Resolution**: Keep a single lowercase `justfile`, confirm it with a path glob after changes, and avoid casing-only cleanup unless the filesystem behavior is known.

**mitigation_type**: workflow

**structural_mitigation_candidate**: If this recurs, add a repository check that flags case-colliding paths before release or commit.

### [2026-07-09] Treat load_all as normal package loading

**Tags**: r, documentation, just, workflow

**Category**: Documentation terminology

**What happened**: README and `justfile` wording described `devtools::load_all()` / `just load` as a package-load smoke test.

**Why it matters**: Loading the package is a legitimate development workflow step in ESPI, not only a validation check, so "smoke test" underdescribes the command and can mislead future docs.

**Resolution**: Describe `just load` as loading the package after setup or dependency changes, and reserve smoke-test language for actual validation scenarios.

**mitigation_type**: documentation

**structural_mitigation_candidate**: Keep command descriptions action-oriented in the `justfile` comments and generated docs.

### [2026-07-09] Raw sample files may sit below the user-described input root

**Tags**: scripts, data-lineage, portability, validation

**Category**: Pipeline portability

**What happened**: The actual 10X sample directories were one level below the user-described raw-input root.

**Why it matters**: Hardcoded paths can silently target the wrong level of a configured data tree and make a script machine-specific.

**Resolution**: Build configuration-root paths from `DATA_ROOT_DIR` and confirm the on-disk nesting before loading sample files.

**mitigation_type**: workflow

**structural_mitigation_candidate**: Keep raw-input path construction at the script parameter boundary rather than embedding machine paths in loading code.

### [2026-07-09] Non-default raw-object handoffs require explicit preprocessing input

**Tags**: scripts, workflow, data-lineage, provenance

**Category**: Pipeline workflow

**What happened**: The persisted raw Seurat object location differs from preprocessing's default input path.

**Why it matters**: Omitting `--input` can make documentation or scripts invoke preprocessing on the wrong path.

**Resolution**: Pass the persisted raw-object path explicitly with `--input` wherever this handoff is documented or invoked.

**mitigation_type**: workflow

**structural_mitigation_candidate**: Keep non-default handoff paths visible at the command boundary.

### [2026-07-09] Read10X integrity does not establish upstream PIPseeker cell calling

**Tags**: qc-filtering, data-lineage, validation, pipseeker, provenance

**Category**: Pipeline validation

**What happened**: A named-vector `Read10X()` import reconciled losslessly, but the raw source matrices still contained a large low-count PIP/barcode population. No PIPseeker selected matrix or selection report was present under the mounted project data root.

**Why it matters**: Lossless import verifies the R handoff, not whether PIPseeker already selected cells upstream.

**Resolution**: Reconcile named-vector imports and validate PIPseeker's selected matrix or report before treating raw source matrices as upstream cell-called input.

**mitigation_type**: workflow

**structural_mitigation_candidate**: Preserve the selected-matrix path and PIPseeker report with future raw-input provenance.

### [2026-07-09] A `^mt-` match can be a partial mitochondrial metric

**Tags**: qc-filtering, mitochondrial, provenance, validation, data-lineage

**Category**: QC metric validation

**What happened**: `^mt-` matched only `mt-Rnr1` and `mt-Rnr2`, while the raw object contained a complete 37-feature mitochondrial set with unprefixed uppercase protein-coding and `Trn*` identifiers. The prior two-rRNA metric was therefore erroneous, not evidence of missing mitochondrial features.

**Why it matters**: Naming patterns cannot establish a QC feature set when a reference uses mixed labels. A lossless `Read10X()` audit verifies the R import but does not establish PIPseeker upstream cell calling.

**Resolution**: Audit feature identifiers directly, calculate `percent.mt` from all 37 observed mitochondrial features, apply the data-specific <=20% cutoff, leave `percent.ribo` diagnostic only, and independently reconstruct saved cell IDs from the selector.

**mitigation_type**: validation

**structural_mitigation_candidate**: Record the explicit mitochondrial feature count and label classes beside every mitochondrial QC threshold.


### [2026-07-09] Keep barcode-rank plots diagnostic

**Tags**: qc-filtering, plotting, pipseeker, workflow, provenance

**Category**: QC visualization workflow

**What happened**: Per-sample PipSeq barcode-rank plots used UMI count and barcode rank to show complexity tails.

**Why it matters**: These plots visualize QC structure but do not distinguish cells from empty droplets.

**Resolution**: Keep barcode-rank plots diagnostic; do not use them for cell calling or filtering in place of PIPseeker selection or a droplet empty-drop method.

**mitigation_type**: workflow

**structural_mitigation_candidate**: Preserve the diagnostic-versus-calling boundary whenever QC visualizations are added.

### [2026-07-12] Source selection belongs at preprocessing

**Tags**: scripts, workflow, data-lineage, reproducibility

**Category**: Pipeline interface

**What happened**: The counts-derived QC object and the established object first converge at preprocessing; later stages already consume branch artifacts from the current object directory.

**Why it matters**: Adding source choices downstream would duplicate state and allow mixed-source artifacts.

**Resolution**: Select the source once with `--input-source legacy|counts-qc`, retain `--input` for custom objects, and regenerate current branches before continuing.

**mitigation_type**: workflow

**structural_mitigation_candidate**: Keep source selection and provenance metadata at the preprocessing seam.

### [2026-07-12] PFlog PCA diagnostics must use the PCA source layer

**Tags**: pflog, plotting, preprocessing, seurat

**Category**: Assay-layer contract

**What happened**: `Seurat::DimHeatmap()` defaults to `scale.data`, which is absent from PFlog PCA objects and caused count-derived PFlog preprocessing to fail.

**Why it matters**: The PFlog branch has a valid PCA reduction but stores its diagnostic expression values in the `pflog` layer.

**Resolution**: Read `misc$preprocessing$pca_source_layer`, validate it, and pass it to `DimHeatmap(slot = ...)`. PFlog preprocessing also materializes the standard log-normalized `data` layer used intentionally by descriptive `FindAllMarkers()` output.

**mitigation_type**: code

**structural_mitigation_candidate**: Treat preprocessing metadata as the authoritative assay-layer contract for downstream diagnostic plots.

### [2026-07-13] MAD-derived thresholds are not necessarily stricter than fixed floors
**Historical implementation note:** The quantitative thresholds below came from the superseded D-45 implementation. The general lesson remains valid, but the current QC script estimates log-scale MAD thresholds among emptyDrops-called scDblFinder singlets above count and feature floors and does not clamp the estimates to fixed floors; the final `pass_qc` selector remains independent of cell-call and singlet flags.

**Tags**: qc-filtering, mitochondrial, doublets, validation, data-lineage

**Category**: QC threshold interpretation

**What happened**: After applying the liberal QC gate, the per-sample median-minus-three-MAD thresholds for `nFeature_RNA` and `nCount_RNA` fell at or below the fixed floors in every sample. Bounding them by the liberal floors made the effective lower thresholds 50 features and 100 counts. Median-plus-three-MAD mitochondrial ceilings tightened below 20% in S2, S5, S7, and S8.

**Why it matters**: A MAD calculation does not guarantee a stricter cutoff. The observed distribution and the direction of the failure mode determine whether it adds filtering beyond the fixed gate.

**Resolution**: Export both the calculated and effective sample-specific thresholds, report retention at the liberal, MAD, doublet, and combined stages, and encode MAD-pass status in the all-cell QC plot.

**mitigation_type**: validation

**structural_mitigation_candidate**: Keep effective threshold tables and staged retention summaries beside future QC-filtered objects.

### [2026-07-13] MG candidate PCA must cover the largest requested clustering dimension

**Tags**: clustering, pflog, mg-selected, workflow, validation

**Category**: Pipeline parameter contract

**What happened**: Rebuilding the MG-selected preprocessing objects with `--dims 20` produced only 20 principal components, while `04-cluster.R` retained the default 20-, 30-, and 50-PC candidate grid. Both MG clustering commands failed when `FindNeighbors()` reached dimensions that had not been computed.

**Why it matters**: The clustering elbow parameter selects the primary candidate but does not limit the default sensitivity grid. The preprocessing PCA depth must cover every requested candidate dimension.

**Resolution**: Re-run `07-select-mg-subset.R` with `--dims 50`, then cluster with elbow 20 and the default 30/50 sensitivity dimensions.

**mitigation_type**: workflow

**structural_mitigation_candidate**: Validate requested clustering dimensions against available PCA dimensions before starting candidate clustering.

### [2026-07-13] MG validators must derive branch-specific chosen columns

**Tags**: clustering, mg-selected, orchestration, validation

**Category**: Branch contract

**What happened**: The first full MG orchestration validated both clustered branches against the no-filter-CC chosen column. The valid filter-CC object failed because its metadata contained the filter-CC branch column instead.

**Why it matters**: A shared validator can silently encode one branch's identity even when every branch follows the same structural contract.

**Resolution**: Derive the expected chosen column from each object's expected branch tag, chosen dimensions, and resolution, then validate that branch-specific column.

**mitigation_type**: test

**structural_mitigation_candidate**: Test each branch handoff independently with exact expected branch tags and chosen columns.

### [2026-07-13] Recursive just calls can suppress canonical dry-run output

**Tags**: just, orchestration, dry-run, testing

**Category**: Command interface

**What happened**: Public recipes initially invoked a hidden helper with a nested `just` process. `just --quiet run-dry-run` propagated quiet behavior to the child process and returned no pipeline plan, so the public interface differed from the direct CLI.

**Why it matters**: Dry-run output is the observable contract for human review and automated testing; a wrapper that hides it defeats the interface.

**Resolution**: Express `run` and `run-dry-run` as parameterized dependencies on the hidden helper in the same Just process. Compare their complete output with the direct runner for counts-qc, legacy, explicit-object, and overwrite cases.

**mitigation_type**: test

**structural_mitigation_candidate**: Keep one public-interface tripwire that executes recipes rather than inspecting the justfile.

### [2026-07-14] Label-permutation tripwires must preserve derived sample identity

**Tags**: tripwires, label-permutation, preprocessing, sample-identity

**Category**: Metamorphic testing

**What happened**: ESPI derives `sample_id` from `Mouse` and `Condition` before blind preprocessing. Permuting `Condition` in the input would therefore change both the interpretation label and sample identity; pre-checkpoint plotting would also write production figures unless explicitly suppressed.

**Why it matters**: A metamorphic test is interpretable only when the tested label is the sole changed input. Otherwise a failure can reflect invalid metadata or side effects rather than label leakage into HVG/PCA.

**Resolution**: Derive and preserve `sample_id` first, then permute `Condition` across intact sample units only in a strictly guarded tripwire mode. Require a checkpoint log and stop boundary, skip figure writes, fix the PCA seed, and compare exact blind-output fingerprints.

**mitigation_type**: test

**structural_mitigation_candidate**: Keep scientific label perturbations behind explicit test-mode guards and assert that identity, counts, non-label metadata, and production paths remain unchanged.

### [2026-07-14] Anchored Markdown wrappers must be removed as a pair

**Tags**: provenance, markdown, parsing, tripwires

**Category**: Text parsing

**What happened**: The provenance validator used `sub("^`|`$", "", target)` on backtick-wrapped paths. Because `sub()` replaces only the first match, it removed the opening backtick and left the closing one, making every valid path appear missing.

**Why it matters**: Alternation between start and end anchors does not express paired-wrapper removal when the replacement function stops after one match.

**Resolution**: Match the complete wrapped token with `^`(.*)`$` and replace it with the captured interior. The finalized provenance tripwire then validated all current strict rows.

**mitigation_type**: test

**structural_mitigation_candidate**: Include wrapped and unwrapped path fixtures in Markdown provenance parsers.

### [2026-07-14] Volcano label ranking and significance encoding can use different statistics

**Tags**: differential-expression, plotting, p-value, fdr

**Category**: Statistical communication

**What happened**: The volcano uses adjusted P values for its y-axis and significance colors, while the requested top-gene labels use raw P values ranked separately within positive and negative fold-change directions.

**Why it matters**: Describing labels only as “most significant” hides which statistic and grouping rule selected them.

**Resolution**: Retain both `pvalue` and `padj` in the plot data, state each role explicitly in the caption, and use deterministic effect-size and gene-name tie-breakers.

**mitigation_type**: documentation

**structural_mitigation_candidate**: Keep selection-statistic wording beside plot captions whenever ranking differs from the displayed significance scale.

### [2026-07-14] Simplify only significant enrichment terms

**Tags**: enrichment, clusterprofiler, enrichit, reproducibility

**Category**: Enrichment workflow

**What happened**: Raw GO ORA/GSEA objects contained thousands of tested terms, while `simplify()` should collapse only the terms retained as enriched. Bayesian term selection is provided by `enrichit::bayes_enrich()` rather than the installed clusterProfiler release.

**Why it matters**: Simplifying all tested terms is both slow and scientifically broader than the reported FDR-significant enrichment set. Package ownership also determines the runtime dependency and namespace call.

**Resolution**: Subset each enrichment object's `@result` to `p.adjust < 0.05`, call `clusterProfiler::simplify()` on that subset, and call seeded `enrichit::bayes_enrich()` before simplifying the exploratory ORA comparison. Do not assign `@pvalueCutoff` because `gseaResult` has no such slot.

**mitigation_type**: code

**structural_mitigation_candidate**: Keep significance filtering, Bayesian selection, simplification, table writing, and plotting in one explicit enrichment stage.
