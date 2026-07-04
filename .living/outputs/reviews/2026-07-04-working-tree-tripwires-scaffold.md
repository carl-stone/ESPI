# Tripwire scaffold — working tree — 2026-07-04

**Mode**: Propose project-specific patches only; no code applied.
**Static review linked**: `.living/outputs/reviews/2026-07-04-working-tree.md`
**Prompted by findings**: F1, F3, F4.

## Plain-English summary

Two checks would help here: one that proves report numbers still match generated outputs, and one that prevents DD cells from being assigned to the wrong Mouse × Condition sample if metadata order drifts.

## Existing hooks found

- **Checkpoint emission**: `R/tripwire-hooks.R::emit_tripwire_checkpoint()` with `CHECKPOINT_LOG` support.
- **Stop-after support**: `STOP_AFTER_CHECKPOINT` in `emit_tripwire_checkpoint()`.
- **Label declaration**: `analysis_labels.yml` exists and declares `Condition`, `Mouse`, `sample_id`, and contrast direction.
- **Drop ledger**: `R/tripwire-hooks.R::write_tripwire_drop_ledger()` exists.
- **Runner**: `tools/run-tripwires.R` exists.

The repo has enough tripwire infrastructure for small runner additions. The patches below are proposals, not applied changes.

## Proposed patches

### Patch 1 — Register the primary DD total-hit count

**Why**: F4 notes that `notebook/sc_analysis.qmd` says primary muscat DD found no FDR-significant genes, but `numbers.json` currently lacks `n_detection_hits`.

**Target**: `scripts/run-mg-selected-de.R`, reportable-values block.

**Proposed shape**:

```r
n_detection_hits <- sum(
  !is.na(full_detection$padj) & full_detection$padj < 0.05
)

reportable_values <- list(
  ...,
  n_detection_hits = n_detection_hits,
  ...
)
```

**Expected effect**: `numbers.json` can become the single source for all quoted primary DD hit counts.

### Patch 2 — Extend the report-numbers-still-match check

**Why**: The existing `tripwire_report_values_freshness()` checks HTML/QMD/figure freshness and two prose contracts, but not MG-selected DD/DE quoted counts.

**Target**: `tools/run-tripwires.R`.

**Add helper**:

```r
read_json_number <- function(path, key) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  pattern <- sprintf('"%s"[[:space:]]*:[[:space:]]*([0-9]+)', key)
  m <- regexec(pattern, text, perl = TRUE)
  hit <- regmatches(text, m)[[1]]
  if (length(hit) < 2) return(NA_integer_)
  as.integer(hit[[2]])
}
```

**Add check inside or beside `tripwire_report_values_freshness()`**:

```r
numbers_path <- file.path(
  Sys.getenv(
    "ESPI_BOX_ROOT",
    file.path(path.expand("~"), "Library/CloudStorage/Box-Box/megan_sc_data")
  ),
  "degs", "mg_selected", "numbers.json"
)

reported_counts <- list(
  n_tested_genes = extract_qmd_number(
    qmd_text,
    "Primary DESeq2 tested ([0-9,]+) genes"
  ),
  n_degs = extract_qmd_number(
    qmd_text,
    "found ([0-9,]+) FDR-significant genes for"
  ),
  n_detection_tested_genes = extract_qmd_number(
    qmd_text,
    "Primary muscat DD tested ([0-9,]+) genes"
  ),
  n_detection_hits = extract_zero_or_number(
    qmd_text,
    "Primary muscat DD tested [0-9,]+ genes and found (no|[0-9,]+) FDR-significant genes"
  ),
  paired_sensitivity_detection_tested_genes = extract_qmd_number(
    qmd_text,
    "paired sensitivity DD tested ([0-9,]+) genes"
  ),
  paired_sensitivity_detection_hits = extract_qmd_number(
    qmd_text,
    "paired sensitivity DD tested [0-9,]+ genes and found ([0-9,]+) FDR-significant genes"
  )
)

for (key in names(reported_counts)) {
  json_value <- read_json_number(numbers_path, key)
  qmd_value <- reported_counts[[key]]
  if (is.na(json_value) || is.na(qmd_value) || !identical(json_value, qmd_value)) {
    problems <- c(problems, sprintf(
      "%s mismatch between numbers.json (%s) and notebook prose (%s)",
      key, json_value, qmd_value
    ))
  }
}
```

Do not hard-code today's biological result values in the runner. `numbers.json`
is the source of truth; the tripwire should fail only when the notebook's quoted
numbers diverge from `numbers.json`, or when a quoted number is missing from
`numbers.json`. If both change together after a legitimate rerun, the tripwire
should pass.

**Expected outcome**: The check fails if someone edits the notebook number, reruns the pipeline without refreshing notebook prose, or forgets to register a quoted number.

### Patch 3 — Add an explicit DD cell-alignment guard

**Why**: F3 notes that `run_detection_muscat_dd()` currently relies on positional alignment between `counts` columns and `cell_metadata` rows after a caller-level check.

**Target**: `scripts/run-mg-selected-de.R`, inside `run_detection_muscat_dd()`.

**Proposed shape**:

```r
assert_cell_metadata_alignment <- function(counts, cell_metadata, stage) {
  if (!identical(colnames(counts), rownames(cell_metadata))) {
    stop(
      stage,
      ": count matrix columns must exactly match cell_metadata rownames.",
      call. = FALSE
    )
  }
}

assert_cell_metadata_alignment(counts, cell_metadata, "muscat DD input")

sample_ids <- rownames(sample_table)
keep_cells <- cell_metadata$pseudobulk_sample_id %in% sample_ids
cell_metadata <- cell_metadata[keep_cells, , drop = FALSE]
counts <- counts[, rownames(cell_metadata), drop = FALSE]

assert_cell_metadata_alignment(counts, cell_metadata, "muscat DD filtered input")
```

**Expected outcome**: Reordered or partially missing cell metadata fails before muscat sees the data.

### Patch 4 — Add the sample-order-does-not-matter check

**Why**: The scientific boundary is “metadata order must not change DD sample labels.” This is a metamorphic check tied to F3.

**Target**: `tools/run-tripwires.R`.

**Practical first version**: static runner check that fails until Patch 3 is present. This avoids running the full DD pipeline and still guards the risky boundary.

```r
tripwire_muscat_dd_alignment_guard <- function(root) {
  slug <- "shuffled-sample-order"
  path <- file.path(root, "scripts", "run-mg-selected-de.R")
  if (!file.exists(path)) {
    return(fail(slug, "scripts/run-mg-selected-de.R is missing."))
  }

  lines <- read_text(path)
  text <- paste(lines, collapse = "\n")
  has_input_guard <- grepl(
    "identical\\(colnames\\(counts\\), rownames\\(cell_metadata\\)\\)",
    text
  )
  has_name_subset <- grepl(
    "counts <- counts\\[, rownames\\(cell_metadata\\), drop = FALSE\\]",
    text
  )

  if (!has_input_guard || !has_name_subset) {
    return(fail(
      slug,
      "run_detection_muscat_dd() must assert count/metadata alignment and subset counts by metadata rownames before assigning muscat sample labels."
    ))
  }

  pass(
    slug,
    "run_detection_muscat_dd() guards count/metadata alignment and avoids positional metadata subsetting."
  )
}
```

Add it to the runner's `results <- c(...)` list.

**Future stronger version**: if the DD helper is later moved into `R/` or made sourceable without running the pipeline, replace the static check with a tiny toy `dgCMatrix` whose metadata rows are deliberately shuffled. The tripwire should pass only if identical DD setup occurs after name-based reordering, or fail before model fitting when names disagree.

## Expected runner output after patches

- `report-values-freshness`: FAIL until `n_detection_hits` is registered and the notebook prose matches all registered values.
- `shuffled-sample-order`: FAIL until `run_detection_muscat_dd()` asserts alignment and subsets counts by metadata rownames.
- Existing checks should keep their current PASS/SKIP behavior unless these patches expose a real freshness or alignment issue.

## Related static findings

- F1 — ungenerated `FindAllMarkers()` prose: `report-values-freshness` should fail if future report prose claims outputs before source artifacts exist.
- F3 — DD count/metadata positional alignment: `shuffled-sample-order` should guard the muscat sample-label boundary.
- F4 — missing primary DD hit registration: `report-values-freshness` should require `n_detection_hits` before accepting the quoted zero.
