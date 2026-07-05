# Tripwire audit & plan — Working tree MG module/p27 heatmaps — 2026-07-05

**Static review**: `.living/outputs/reviews/2026-07-05-working-tree-mg-module-p27-heatmaps.md`  
**Mode**: audit — this document describes what would be tested; no code was run and nothing was modified  
**Hooks the project has**: full — checkpoint emission (`emit_tripwire_checkpoint()`), stop-after support (`STOP_AFTER_CHECKPOINT`), label declarations (`analysis_labels.yml`), drop-ledger hook (`DROP_LEDGER`), and `tools/run-tripwires.R`

## What this audit does

Walks the heatmap-related working-tree changes, names the failure modes worth testing, says what each test would do, and points at the cheapest check available today.

## Tests that would apply here

### tripwire1. The report-numbers-still-match check

**Watches for**: generated notebook/report text or Mycelium provenance that no longer matches source artifacts.

**How it'd work**: fingerprint `notebook/sc_analysis.qmd`, `notebook/sc_analysis.html`, the two heatmap PNG symlinks, the four heatmap TSVs, `.living/log/2026-07-05-007-espi.md`, and `.living/log/LOG_REGISTRY.md`; fail if the notebook or registry claims outputs/validation that are stale or missing.

**Slug**: `report-values-freshness` (freshness category)

**Related static finding**: F5 in the static review.

**Today**: partly runnable through `tools/run-tripwires.R`, which already checks report freshness for notebook figures. It does not yet validate Mycelium registry row semantics.

**Starter check**: compare the registry row and linked log against the expected artifact basenames: both heatmap stems, `_module_scores.tsv`, `_p27_enrichment.tsv`, and the rendered notebook path.

### tripwire2. The missing-input-file blow-up check

**Watches for**: a heatmap script silently falling back to a different Seurat object when the requested input path is wrong.

**How it'd work**: run `scripts/plot-cluster-marker-heatmaps.R` with a deliberately missing `--input`; assert no PNG/PDF/TSV/symlink is written and the process fails before the first analysis boundary.

**Slug**: `missing-counts-file` (fault-injection category)

**Related static finding**: none; starter check.

**Today**: the existing tripwire runner exercises missing-input failure for preprocessing, not specifically this heatmap script.

**Starter check**: run the heatmap script with `--input /tmp/does-not-exist.rds --out-dir <scratch>` and confirm it exits non-zero with `Input Seurat object does not exist` and no files under the scratch output directory.

### tripwire3. The labels-don't-leak-into-blind-analysis check

**Watches for**: condition labels influencing blind preprocessing or clustering steps that should remain label-independent.

**How it'd work**: permute the declared `Condition` labels before blind preprocessing/clustering, rerun to the blind checkpoint, and assert HVGs/PCA/cluster-independent artifacts are unchanged.

**Slug**: `label-permutation` (metamorphic category)

**Related static finding**: none; starter check.

**Today**: `tools/run-tripwires.R` includes a static firewall and marks full label permutation as skipped because there is no scratch-output hook for safe blind reruns.

**Starter check**: keep the existing static firewall result visible; scaffold scratch-output isolation before executing the full metamorphic run.

### tripwire4. The comparison-direction-not-flipped check

**Watches for**: target/reference signs reversing in condition contrasts.

**How it'd work**: run a tiny known-answer dataset through the DE/DD contrast path and assert the marker effect is positive for `p27CKO + E-Stim vs. p27CKO`.

**Slug**: `toy-contrast-direction` (known-answer category)

**Related static finding**: none; starter check.

**Today**: already runnable through `tools/run-tripwires.R`; this review did not run it as part of the audit mode.

**Starter check**: run `Rscript tools/run-tripwires.R` and inspect the `toy-contrast-direction` row.

### tripwire5. The helper-does-not-disturb-randomness check

**Watches for**: exported helper code changing the caller's random-number stream.

**How it'd work**: set a seed, draw one random value, call `compute_cluster_p27_enrichment()` on a tiny Seurat object, draw another random value, and compare it with the value expected when the helper preserves `.Random.seed`.

**Slug**: `rng-state-preservation` (known-answer/metamorphic category; project-specific)

**Related static finding**: F1 in the static review.

**Today**: not covered by the existing tripwire runner.

**Starter check**: add a small helper-level R check alongside the synthetic helper test used during implementation; it does not need the full pipeline.

## Hooks this project would need to run the rest

The repo already has the four core hooks. Two targeted additions would make these heatmap-specific checks executable rather than manual:

- A scratch-output option for heatmap scripts so fault-injection checks cannot touch Box outputs or notebook symlinks.
- A Mycelium provenance semantic check that validates completed `LOG_REGISTRY.md` rows have non-empty Summary, Key Outputs, Status, Tags, and a matching linked log.

## What this audit does NOT cover

- Visual aesthetics of the heatmap panels.
- Numerical regression against a prior biological result.
- Performance of 2,000 p27 permutations on much larger future objects.

## Implementation status

- Added executable `tools/run-tripwires.R` checks for heatmap missing-input failure, p27 RNG-state preservation, and review-scoped Mycelium provenance semantics.
- The Mycelium provenance check intentionally covers only heatmap session rows `2026-07-05-007` and `2026-07-05-008`; broader historical registry backfill is a separate task.
- `Rscript tools/run-tripwires.R` passes these added checks.
