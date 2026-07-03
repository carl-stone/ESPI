# Tripwire audit & plan — working tree — 2026-07-03

**Static review**: `.living/outputs/reviews/2026-07-03-working-tree.md`
**Mode**: audit — this document describes what would be tested; no code was run and nothing was modified.
**Hooks the project has**: none executable yet. I found convention text mentioning checkpoints, but no first-party checkpoint log, `--stop-after` / `STOP_AFTER_CHECKPOINT`, `analysis_labels.yml`, drop ledger, or tripwire runner.

## What this audit does

Walks the ESPI preprocessing / clustering / notebook path, names the failure modes worth testing, says what each test would actually do, and points at the cheapest check available today without adding instrumentation.

## Tests that would apply here

### tripwire1. The branch-outputs-stay-separated check

**Watches for**: filtered and unfiltered cell-cycle branches overwriting each other or becoming indistinguishable in object names, reduction names, cluster columns, and plot filenames.

**How it'd work**: run clustering on two Seurat objects that differ only in `sobj@misc$preprocessing$filtered_cell_cycle`, with the same normalization/elbow/resolution settings. The tripwire would assert that every persisted artifact contains both normalization and cell-cycle branch tags, and that no output path is reused by both branches.

**Slug**: `branch-artifact-collision` (freshness / artifact-provenance category)

**Related static finding**: F1 — cluster outputs can silently collide across cell-cycle branches.

**Today**: cannot run end-to-end without instrumentation, but the static review already found the likely failure: `scripts/cluster-sobj.R` names cluster columns, UMAP reductions, clustree tags, and clustered RDS files without a `filter-cc` / `no-filter-cc` tag.

**Starter check**: a small static check can parse `scripts/cluster-sobj.R` and assert each persisted name built from `norm` also includes `sobj@misc$preprocessing$filtered_cell_cycle` or a derived `cc_tag`. This would catch the current collision class without running Seurat.

### tripwire2. The report-is-current check

**Watches for**: rendered HTML or figure prose that no longer matches the QMD source, figure targets, or plotting parameters.

**How it'd work**: fingerprint `notebook/sc_analysis.qmd`, all referenced `notebook/figures/*.png` symlink targets, and the rendered `notebook/sc_analysis.html`; fail if the rendered HTML is older than any source or if registered figure-parameter claims disagree with the code that generated the figures.

**Slug**: `report-values-freshness` (freshness category)

**Related static findings**: F2, F3, F4, F9 — stale embedded HTML, HVG top-10/top-20 mismatch, DimHeatmap `cells = 500` prose mismatch, and validation claims without attached logs.

**Today**: partially checkable without instrumentation. File mtimes already show `notebook/sc_analysis.html` is older than `notebook/sc_analysis.qmd` and at least one referenced figure target.

**Starter check**: a 20-40 line R or Python script can read the QMD, resolve every `figures/*.png` symlink target, compare mtimes/fingerprints against the HTML, and grep for known prose/code contracts such as `n_top = 20`, `top 20`, `cells = 500`, and “500 cells”.

### tripwire3. The missing-input-file blow-up check

**Watches for**: a required Seurat object path being missing or wrong, but the pipeline silently falling back to another input.

**How it'd work**: call `scripts/preprocess-sobj.R --input <missing-file.rds>` and `scripts/cluster-sobj.R --input <missing-file.rds> --elbow-n <N>`. The pipeline should fail before a `raw_data_available` checkpoint and should not emit new objects or figures.

**Slug**: `missing-counts-file` (fault-injection category)

**Related static finding**: none direct; this is one of the starter checks because missing-input fallbacks are a common silent scientific-data failure mode.

**Today**: likely to fail via `readRDS()` for explicit missing inputs, but there is no checkpoint log to prove where it failed or to assert no partial artifacts were written.

**Starter check**: run a dry static audit of both scripts and path constants: verify every source path flows to `readRDS()` and no `tryCatch`, alternate cache, or fallback write path is present. A future executable version should use a temporary output directory and assert it remains empty on failure.

### tripwire4. The missing-metadata-sample check

**Watches for**: missing or malformed sample metadata passing through preprocessing and creating misleading `sample_id` values or unrecorded sample/cell drops.

**How it'd work**: clone a small input Seurat object, blank or remove one `Mouse` / `Condition` metadata value, then run preprocessing until a `samples_reconciled` checkpoint. The checkpoint should fail before `sample_id` construction or record an explicit drop in a ledger.

**Slug**: `missing-metadata-sample` (fault-injection category)

**Related static finding**: none direct; this is a starter check because `sample_id` and Mouse × Condition are load-bearing for pseudobulk interpretation.

**Today**: cannot run cleanly because there is no sample-reconciliation checkpoint, no declared required metadata schema, and no drop ledger.

**Starter check**: add a static assertion review around `scripts/preprocess-sobj.R` lines 58-63: required columns are `Mouse` and `Condition`; values should be non-missing before `sample_id <- paste0(...)`; generated `sample_id` values should be unique or intentionally many-to-one by cell. This can be checked against the current input object once, outside the full pipeline.

### tripwire5. The labels-don't-leak-into-blind-analysis check

**Watches for**: condition labels influencing HVG selection, PCA, UMAP, or clustering before labels are supposed to be used for interpretation.

**How it'd work**: declare `Condition` as a label column, permute it row-wise in the Seurat metadata, rerun blind preprocessing through HVG/PCA, and compare branch artifacts that should not depend on labels. HVG lists, PCA embeddings/loadings, and blind QC summaries should be invariant apart from metadata-derived labels/plots.

**Slug**: `label-permutation` (metamorphic category)

**Related static finding**: none direct; this is a starter check because the analysis uses condition labels and blind preprocessing/clustering diagnostics.

**Today**: cannot run as a tripwire because there is no `analysis_labels.yml`, no blind-stage policy, and no stop-after checkpoint for HVG/PCA.

**Starter check**: create a one-page label declaration by hand for review purposes: `Condition` is a label; `Mouse` is a biological replicate/blocking variable; `sample_id` is derived metadata. Then statically confirm `Condition` is not used by `FindVariableFeatures`, `run_log1p_pca()`, `run_pflog_pca()`, or `RunUMAP()` before interpretation plots.

### tripwire6. The comparison-direction-not-flipped check

**Watches for**: future differential-expression contrasts reversing target/reference or reporting the opposite sign from the intended condition comparison.

**How it'd work**: run the DE module on a tiny synthetic pseudobulk dataset where one marker gene is known to be higher in `p27CKO +EStim` than `p27CKO`. The tripwire would assert the reported log fold-change sign and contrast label match the intended direction.

**Slug**: `toy-contrast-direction` (known-answer category)

**Related static finding**: none in current code; this is a starter check because condition-level DE is a planned/load-bearing analysis and contrast direction errors are common.

**Today**: cannot run because the current working-tree changes do not add a DE script/contract to execute.

**Starter check**: before writing DE code, record the contrast convention in a small spec: numerator / target = `p27CKO +EStim`; denominator / reference = `p27CKO`; statistical unit = Mouse × Condition pseudobulk sample. The toy-data tripwire can then be added at the same time as the DE function.

## Hooks this project would need to execute tripwires

These tests can be executed only when the analysis pipeline emits enough signal to observe. The minimum ESPI-specific hooks would be:

- A **checkpoint log** under the Box output root or a temp run directory, with one structured record per boundary: `raw_data_available`, `metadata_complete`, `variable_features_selected`, `pca_ready`, `cluster_artifacts_written`, `report_values_ready`.
- A **stop-after** option, probably via `STOP_AFTER_CHECKPOINT`, so `preprocess-sobj.R` and `cluster-sobj.R` can halt after HVG/PCA or clustering without running the full downstream path.
- An **analysis label declaration**, e.g. `analysis_labels.yml`, naming `Condition` as the treatment label, `Mouse` as the replicate/blocking variable, and `sample_id` as derived metadata.
- A **drop ledger** for any dropped sample/cell/object during metadata checks or filtering, with stage, ID, reason, and whether the drop is allowed by policy.
- A **run scratch/output override** so fault-injection tripwires write into a temp directory rather than the real Box artifact tree.

## What this audit does NOT cover

- Full numerical regression against a prior accepted ESPI run.
- Performance or runtime regressions.
- A biological correctness review of the raw Trailmaker object.
- Running Seurat/PFlog/Quarto; this audit describes tests only.
