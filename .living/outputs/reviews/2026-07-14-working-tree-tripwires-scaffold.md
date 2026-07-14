# Tripwire scaffold proposal — working tree — 2026-07-14

ESPI already has the four core hooks needed for behavioral checks: structured checkpoints, stop-after support, `analysis_labels.yml`, and a drop ledger, plus `tools/run-tripwires.R`. The missing piece is safe scratch-input/output routing for the two checks that must perturb real analysis inputs. This proposal adds that routing and extends the existing runner; it does not apply patches or run the analysis.

**Static review**: [2026-07-14-working-tree.md](2026-07-14-working-tree.md)

## Existing hooks to reuse

- `R/tripwire-hooks.R::emit_tripwire_checkpoint()` already writes append-only TSV checkpoints through `CHECKPOINT_LOG` and honors `STOP_AFTER_CHECKPOINT`.
- `R/tripwire-hooks.R::write_tripwire_drop_ledger()` already records rejected metadata rows through `DROP_LEDGER`.
- `analysis_labels.yml` already declares `Mouse`, `Condition`, derived `sample_id`, blind-stage policy, and contrast direction.
- `tools/run-tripwires.R` already implements the missing-input-file, report-freshness, metadata-contract, label-firewall, and toy-contrast checks.
- Current limitation: `label-permutation` stops after a static firewall and reports `SKIP` because preprocessing has no scratch-safe output contract.

## Proposed patches

### Patch 1 — Add scratch-safe raw-input overrides

**Files**: `scripts/01-process-counts.R`, `scripts/run-pipeline.R`

Add optional expert/test-only arguments to `scripts/01-process-counts.R`:

```r
--raw-counts-dir <path>
--metadata <path>
--output <path>
```

Keep the current `DATA_ROOT_DIR` paths as defaults. Validate that explicit paths exist before reading. Write only to the explicit `--output` when supplied. Emit these boundaries through the existing helper:

```r
emit_tripwire_checkpoint(
  "raw_data_available",
  raw_counts_dir = raw_counts_dir,
  metadata_path = metadata_path,
  n_metadata_rows = nrow(metadata)
)

validate_required_metadata(
  metadata,
  required_metadata_columns,
  stage = "samples_reconciled"
)

emit_tripwire_checkpoint(
  "samples_reconciled",
  n_samples = nrow(metadata),
  sample_ids = paste(sort(metadata$Sample), collapse = ",")
)
```

Do not expose these overrides through routine `just run`; the tripwire runner invokes the script directly against a temporary directory. `scripts/run-pipeline.R` only needs to preserve `CHECKPOINT_LOG`, `STOP_AFTER_CHECKPOINT`, and `DROP_LEDGER` in child processes, which `system2()` already inherits unless explicitly scrubbed.

**Enables**: the missing-metadata-sample check (`missing-metadata-sample`, fault-injection category) without touching Box inputs.

### Patch 2 — Add one scratch output root for preprocessing

**Files**: `scripts/03-preprocess.R`, optionally `R/paths.R`

Add an optional `--output-root <path>` argument. When present, derive the object and figure destinations beneath that root rather than `CURRENT_OBJECT_DIR` and `FIGURE_DIR`; keep production defaults unchanged. Every figure helper used by this script must receive the scratch destination or run under a temporary `DATA_ROOT_DIR` resolved before `devtools::load_all()`.

Preferred project-specific shape:

```r
output_root <- arg_value("--output-root", default = NULL)
object_dir <- if (is.null(output_root)) CURRENT_OBJECT_DIR else file.path(output_root, "objects")
figure_dir <- if (is.null(output_root)) file.path(FIGURE_DIR, "preprocess") else file.path(output_root, "figures")
```

Do not add a second checkpoint implementation. Continue using `emit_tripwire_checkpoint()` and stop after `pca_ready`, before persistent object output, for the label permutation comparison.

At `pca_ready`, record a deterministic blind-output fingerprint computed from label-independent outputs, for example the sorted HVG vector plus rounded PCA singular values. Store the algorithm name with the digest:

```r
emit_tripwire_checkpoint(
  "blind_qc_complete",
  fingerprint_algorithm = "sha256",
  hvg_hash = digest::digest(sort(VariableFeatures(sobj)), algo = "sha256"),
  pca_sdev_hash = digest::digest(round(Stdev(sobj[[reduction]]), 10), algo = "sha256")
)
```

Declare `digest` in `DESCRIPTION` if it is not already available through an imported dependency.

**Enables**: the labels-don't-leak-into-blind-analysis check (`label-permutation`, metamorphic category) with no production writes.

### Patch 3 — Extend the existing runner instead of creating another one

**File**: `tools/run-tripwires.R`

Add two executable functions to the current `checks` list.

#### The missing-metadata-sample check

1. Copy the six-row metadata file into `tempdir()`.
2. Remove one sample row while leaving all six count directories linked read-only.
3. Set per-test `CHECKPOINT_LOG` and `DROP_LEDGER` paths.
4. Invoke `scripts/01-process-counts.R` with explicit scratch paths and `STOP_AFTER_CHECKPOINT=samples_reconciled`.
5. Pass only if the process exits non-zero, `samples_reconciled` is absent, and the ledger records the missing sample or count/metadata mismatch as `allowed_by_policy = FALSE`.

Expected scientific boundary: the count directories and metadata must reconcile by `Sample`; positional or partial matching must never continue.

#### The labels-don't-leak-into-blind-analysis check

1. Copy the QC-filtered Seurat object into a temporary workspace.
2. Run `scripts/03-preprocess.R` once with original `Condition` values and once with `Condition` permuted across Mouse × Condition samples.
3. Set separate checkpoint logs, one scratch output root per run, and `STOP_AFTER_CHECKPOINT=blind_qc_complete`.
4. Pass only if the declared HVG and PCA fingerprints match exactly.
5. Fail if either run writes beneath production `CURRENT_OBJECT_DIR` or `FIGURE_DIR`.

The existing static label firewall remains useful as a fast preflight, but should no longer be the final assertion.

### Patch 4 — Strengthen provenance semantics around F1

**File**: `tools/run-tripwires.R`

Extend `tripwire_mycelium_provenance_semantics()` from its two historical review-scoped rows to every newly added `complete` registry row. For each complete row, require:

- non-empty Summary, Key Outputs, Status, and Tags fields;
- a linked log that exists;
- frontmatter duration/files matching the registry;
- `## Session Summary`, `## Key Outputs`, and `## Status` sections;
- every repository path in `Files Modified` to exist or be explicitly marked deleted;
- no `.log-scribe-*` authentication-error artifact included as session provenance.

Use a tiny known-answer fixture inside `tempdir()` containing one valid and one malformed registry/log pair. The check should prove that the malformed pair fails before applying the same validator to `.living/log/LOG_REGISTRY.md`.

**Linked static finding**: F1 in the static report.

## Starter behavioral checks after these patches

| Plain-English check | Stable slug | Category | Perturbation | Expected result |
|---|---|---|---|---|
| Missing-input-file blow-up check | `missing-counts-file` | fault-injection | Point preprocessing at a nonexistent RDS | Non-zero exit; no `raw_data_available` checkpoint |
| Missing-metadata-sample check | `missing-metadata-sample` | fault-injection | Remove one metadata row in a scratch copy | Reconciliation fails; rejected sample is in the ledger |
| Labels-don't-leak-into-blind-analysis check | `label-permutation` | metamorphic | Permute `Condition` in a scratch object | HVG and PCA fingerprints remain identical |
| Comparison-direction-not-flipped check | `toy-contrast-direction` | known-answer | Run the existing synthetic two-condition contrast | Positive effect means E-Stim > control |
| Report-numbers-still-match check | `report-values-freshness` | freshness | Compare QMD/HTML/figure timestamps and active values | Rendered report is current and claims match sources |
| Session-records-are-auditable check | `mycelium-provenance-semantics` | known-answer | Validate a malformed fixture, then current records | Fixture fails and every complete current row passes |

## Acceptance criteria for a later implementation

- No tripwire reads from or writes to production Box outputs except through read-only source links.
- Every run gets unique checkpoint, ledger, and output paths under `tempdir()`.
- Environment variables are restored with `on.exit()` even when a child process fails.
- The two new behavioral checks fail against deliberately broken fixtures and pass against current inputs.
- `label-permutation` changes from `SKIP` to `PASS` only after comparing real blind-output fingerprints.
- Existing `just tripwires` behavior and concise PASS/FAIL/SKIP table remain unchanged for callers.
- `devtools::document()` runs after any `R/tripwire-hooks.R` change; this proposal does not require changing that helper.

## Not proposed

- No second runner, checkpoint helper, label file, or drop-ledger format.
- No template-generated Python harness in this R analysis.
- No automatic patch application.
- No full canonical pipeline run for each perturbation; stop at the named boundary.
