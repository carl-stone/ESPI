# Tripwire scaffold — unstaged working tree — 2026-07-04

**Mode**: Propose patches only; no patches applied.
**Static review linked**: `.living/outputs/reviews/2026-07-04-unstaged-working-tree.md`

## Plain-English summary

The review found one behavioral boundary worth instrumenting: value-bearing command-line flags must fail when present without a value. Today that boundary is inconsistent across refactored scripts.

## Instrumentation already present

- **Checkpoint emission**: `emit_tripwire_checkpoint()` in `R/tripwire-hooks.R` and calls in preprocessing/clustering scripts.
- **Stop-after hook**: `STOP_AFTER_CHECKPOINT` support in `emit_tripwire_checkpoint()`.
- **Label contract**: `analysis_labels.yml` exists.
- **Drop ledger**: `write_tripwire_drop_ledger()` exists.
- **Tripwire runner**: `tools/run-tripwires.R` exists.

The repo has enough infrastructure to add this as a normal runner check.

## Proposed patch set

### Patch 1 — Add consistent value-argument helpers in scripts that parse CLI args

Target files:

- `scripts/preprocess-sobj.R`
- `scripts/cluster-all.R`
- `scripts/cluster-sobj.R`

Add a small local helper near each existing `arg()` function:

```r
arg_value <- function(name, default = NULL, required = FALSE) {
  value <- arg(name)
  if (identical(value, TRUE)) {
    stop("Missing value for ", name, call. = FALSE)
  }
  if (is.null(value)) {
    if (required) {
      stop("Missing required argument ", name, call. = FALSE)
    }
    return(default)
  }
  value
}
```

Keep boolean flags separate:

```r
arg_flag <- function(name) {
  identical(arg(name), TRUE)
}
```

Then replace risky parsing patterns:

```r
# preprocess-sobj.R
input <- arg_value(
  "--input",
  default = file.path(INPUT_OBJECT_DIR, "pipseq_processed_matrix_with_egfp.rds")
)
normalization <- arg_value("--normalization", default = "log1p")
filter_cc <- arg_flag("--filter-cell-cycle") ||
  identical(tolower(arg_value("--filter-cell-cycle", default = "false")), "true")
```

```r
# cluster-all.R
elbow_n <- arg_value("--elbow-n", default = "20")
input_dir <- arg_value("--input-dir", default = CURRENT_OBJECT_DIR)
extra_dims <- arg_value("--extra-dims", default = NULL)
resolutions <- arg_value("--resolutions", default = NULL)
dry_run <- arg_flag("--dry-run")
```

```r
# cluster-sobj.R
input <- arg_value("--input", required = TRUE)
elbow_n <- as.integer(arg_value("--elbow-n", required = TRUE))
extra_dims <- parse_csv_int(arg_value("--extra-dims", default = NULL), default = c(30, 50))
resolutions <- parse_csv_num(arg_value("--resolutions", default = NULL), default = c(0.3, 0.5, 0.8))
```

Notes:

- For `--filter-cell-cycle`, preserve current behavior where bare presence means `TRUE` and explicit `true` also means `TRUE`.
- For all value-bearing flags, absent means default; present-without-value means hard failure.

### Patch 2 — Add a CLI missing-value tripwire to `tools/run-tripwires.R`

Add a new static/behavioral check named `cli-value-boundaries`.

Suggested test cases:

```r
tripwire_cli_value_boundaries <- function(root) {
  slug <- "cli-value-boundaries"
  rscript <- file.path(R.home("bin"), "Rscript")

  current_object_dir <- file.path(
    path.expand("~/Library/CloudStorage/Box-Box"),
    "megan_sc_data",
    "seurat_objects",
    "current"
  )
  preprocess_input <- file.path(
    current_object_dir,
    "preprocess_pflog_filter-cc.rds"
  )

  commands <- list(
    list(
      script = file.path(root, "scripts", "cluster-all.R"),
      args = c("--dry-run", "--elbow-n"),
      expected = "Missing value for --elbow-n"
    ),
    list(
      script = file.path(root, "scripts", "preprocess-sobj.R"),
      args = c("--input", "--normalization", "pflog"),
      expected = "Missing value for --input"
    )
  )

  if (file.exists(preprocess_input)) {
    commands <- c(commands, list(list(
      script = file.path(root, "scripts", "cluster-sobj.R"),
      args = c("--input", preprocess_input, "--elbow-n"),
      expected = "Missing value for --elbow-n"
    )))
  }

  failures <- character()
  for (command in commands) {
    output <- system2(
      rscript,
      c(command$script, command$args),
      stdout = TRUE,
      stderr = TRUE
    )
    status <- attr(output, "status")
    if (is.null(status) || identical(as.integer(status), 0L)) {
      failures <- c(failures, paste("unexpected success:", command$script))
    } else if (!any(grepl(command$expected, output, fixed = TRUE))) {
      failures <- c(
        failures,
        paste("missing expected error for", command$script, ":", command$expected)
      )
    }
  }

  if (length(failures) > 0) {
    return(fail(slug, paste(failures, collapse = "; ")))
  }
  pass(slug, "Value-bearing CLI flags fail when present without values.")
}
```

Register it in the runner's result list near the existing script-contract checks:

```r
results <- c(
  results,
  list(tripwire_cli_value_boundaries(root))
)
```

Implementation notes:

- The `cluster-sobj.R` case should skip if the preprocess object is absent, to avoid requiring large data on machines without Box mounted.
- `cluster-all.R --dry-run --elbow-n` is safe and should not run analysis.
- `preprocess-sobj.R --input --normalization pflog` should fail before reading data after Patch 1.

### Patch 3 — Add the heatmap symlink side effect to the script header

Target file:

- `scripts/big-heatmap-plot.R`

Add an `Outputs:` block to the top docs:

```r
# Outputs:
#   <out-dir>/cell_type_marker_heatmap_<layer>_cells_dims<dims>_res<resolution>.png
#   <out-dir>/cell_type_marker_heatmap_<layer>_cells_dims<dims>_res<resolution>.pdf
#   notebook/figures/<same PNG basename> as a symlink to the PNG output.
#   Existing notebook figure paths with the same name are unlinked first.
```

## Expected outcomes after applying patches

The following malformed invocations should all exit nonzero before reading/writing analysis artifacts:

```sh
Rscript scripts/cluster-all.R --dry-run --elbow-n
Rscript scripts/preprocess-sobj.R --input --normalization pflog
Rscript scripts/cluster-sobj.R --input <valid-preprocess-object.rds> --elbow-n
```

The normal smoke tests should still pass:

```sh
Rscript scripts/cluster-all.R --dry-run --elbow-n 20
Rscript scripts/big-heatmap-plot.R
Rscript scripts/summarize-cluster-grid.R
Rscript tools/run-tripwires.R
```

## Related static findings

- F1: bare `--elbow-n` becomes one PC in `scripts/cluster-sobj.R`.
- F2: valueless preprocessing flags fall back to defaults in `scripts/preprocess-sobj.R`.
- F3: batch clustering wrapper ignores missing values for value-bearing flags in `scripts/cluster-all.R`.
- F4: heatmap docs omit notebook symlink side effect.
