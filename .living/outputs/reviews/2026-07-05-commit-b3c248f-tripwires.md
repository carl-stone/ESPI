# Tripwire audit & plan — commit b3c248f — 2026-07-05

**Static review**: `.living/outputs/reviews/2026-07-05-commit-b3c248f.md`  
**Mode**: audit — this document describes what would be tested; no code was run and nothing was modified  
**Hooks the project has**: full/near-full for current ESPI tripwires: checkpoint emission via `emit_tripwire_checkpoint()`, `STOP_AFTER_CHECKPOINT`, `analysis_labels.yml`, drop-ledger helpers via `write_tripwire_drop_ledger()`, and `tools/run-tripwires.R`.

## What this audit does

This audit names the behavioral checks worth applying to the MG-selected figure revision. Static review found no Major issues, so these tests are mostly guardrails against stale reports, missing inputs, label leakage, and flipped contrast direction.

## Tests that would apply here

### tripwire1. The report-numbers-still-match check

**Watches for**: notebook prose or rendered HTML values that no longer match the source tables or registered numbers.  
**How it'd work**: rerun the report freshness check against `notebook/sc_analysis.qmd`, `notebook/sc_analysis.html`, referenced figures, and `degs/mg_selected/numbers.json`; fail if the notebook is older than sources or quotes values that disagree with registered outputs.  
**Slug**: `report-values-freshness` (freshness category)  
**Related static finding**: none directly; this is the standard documentation-change check for this diff.  
**Today**: already supported by `tools/run-tripwires.R`; the commit's verification run passed it with 21 referenced figure targets.  
**Starter check**: rerun `Rscript tools/run-tripwires.R` after any notebook prose or figure-source change.

### tripwire2. The missing-input-file blow-up check

**Watches for**: an analysis script silently continuing when a required Seurat/counts input is missing.  
**How it'd work**: point a script such as `scripts/preprocess-sobj.R`, `scripts/plot-mg-selected-figures.R`, or `scripts/run-mg-selected-de.R` at a deliberately missing input path; the pipeline should stop before writing downstream artifacts.  
**Slug**: `missing-counts-file` (fault-injection category)  
**Related static finding**: none; static review found the touched code uses explicit file/path validation.  
**Today**: partially covered by the existing `missing-counts-file` tripwire in `tools/run-tripwires.R`; extending it to the new abundance-producing script would be the next specific addition.  
**Starter check**: run `Rscript scripts/plot-mg-selected-figures.R --input /tmp/does-not-exist.rds` in a scratch environment and confirm it exits non-zero without writing figures.

### tripwire3. The missing-metadata-row check

**Watches for**: sample/cell metadata and expression matrices falling out of alignment without a hard failure.  
**How it'd work**: remove or corrupt a metadata row in a scratch Seurat object, rerun the relevant plotting or DE step, and confirm the code stops at the explicit metadata-alignment checks rather than silently dropping or recycling cells.  
**Slug**: `missing-metadata-sample` (fault-injection category)  
**Related static finding**: none; static review found explicit checks in the touched abundance and UMAP paths.  
**Today**: the project has checkpoint and drop-ledger hooks, but this exact perturbation is not currently run for `scripts/plot-mg-selected-figures.R` or `R/cluster-abundance.R`.  
**Starter check**: add a scratch-object test that deletes one `sobj@meta.data` row before `compute_cluster_abundance()` and expects an alignment or metadata error.

### tripwire4. The labels-don't-leak-into-blind-analysis check

**Watches for**: condition labels influencing blind preprocessing, HVG selection, PCA, or UMAP/clustering before the intended label-aware stages.  
**How it'd work**: permute declared label columns from `analysis_labels.yml`, rerun blind preprocessing/clustering checkpoints, and confirm blind artifacts stay unchanged.  
**Slug**: `label-permutation` (metamorphic category)  
**Related static finding**: none; this remains a high-value guardrail for single-cell analyses.  
**Today**: `tools/run-tripwires.R` performs a static firewall and reports the full permutation as skipped because the project lacks a scratch-output hook for safe blind HVG/PCA/UMAP reruns.  
**Starter check**: keep the current static firewall until a scratch-output mode exists.

### tripwire5. The comparison-direction-not-flipped check

**Watches for**: DE/DD effect signs flipping because the reference and target condition levels were reversed.  
**How it'd work**: run a tiny known-answer contrast where the E-Stim group has a known positive marker effect and assert the reported DE/DD direction is positive for E-Stim vs control.  
**Slug**: `toy-contrast-direction` (known-answer category)  
**Related static finding**: none; Figure 12 axis labels and DE/DD contrast wording were consistent in static review.  
**Today**: already supported by `tools/run-tripwires.R`; the commit's verification run passed `toy-contrast-direction`.  
**Starter check**: rerun the tripwire runner after any change to condition labels, design formulas, or Figure 12 axis text.

### tripwire6. The figure-reference-still-points-to-the-right-panel check

**Watches for**: audit logs or manuscript prose using stale auto-numbered figure labels after inserting or restoring panels.  
**How it'd work**: parse rendered Quarto figure IDs/captions and scan `.living/log/*.md` plus notebook prose for `Figure N` references; flag any numeric reference that does not also include a stable Quarto figure ID or caption.  
**Slug**: no existing standard slug; this is a documentation-fidelity audit suggested by F1.  
**Related static finding**: F1 in the static review.  
**Today**: not implemented in `tools/run-tripwires.R`; the cheap fix is to avoid numeric figure labels in session logs.  
**Starter check**: grep session logs for `Figure [0-9]+` after notebook figure insertions and replace with stable IDs/captions.

## Hooks this project would need to run the rest

The repo already has the main hooks needed for many tripwires:

- **Checkpoint log**: `emit_tripwire_checkpoint()` writes named checkpoints when `CHECKPOINT_LOG` is set.
- **Stop-after option**: `STOP_AFTER_CHECKPOINT` stops after a named checkpoint.
- **Label declaration**: `analysis_labels.yml` declares `Mouse`, `Condition`, and derived `sample_id` boundaries.
- **Drop ledger**: `write_tripwire_drop_ledger()` writes drop records when `DROP_LEDGER` is set.
- **Runner**: `tools/run-tripwires.R` executes the current project tripwires.

The main missing piece for the skipped label-permutation test is a scratch-output mode that can rerun blind HVG/PCA/UMAP stages safely without touching canonical Box outputs.

## What this audit does NOT cover

- Visual taste judgments beyond whether figure references and encodings stay synchronized.
- Runtime performance of regenerating the large notebook and figures.
- Numerical snapshot regression against every prior figure pixel.
