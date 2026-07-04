# Review — unstaged working tree — 2026-07-04

**Scope**: Unstaged working-tree diff against `HEAD`
**Files reviewed**: 10
**Sub-agents run**: 6 (`stats-causal`, `data-pipeline-leakage`, `bioinformatics`, `llm-failure-modes`, `doc-schema-fidelity`, `code-quality`)

## Key decisions in this analysis

- **Top-to-bottom script style** — The diff adopts RStudio-step-friendly executable scripts with package loading, parameters, validation/work sections, and output side effects near the end.
- **`here::i_am()` for script root discovery** — Scripts now rely on `here::i_am()` plus `devtools::load_all(here::here(), export_all = FALSE, quiet = TRUE)` instead of ad hoc working-directory handling.
- **`system2()` orchestration** — Wrapper scripts now invoke child R scripts with argument vectors and explicit nonzero-exit checks instead of shell command strings.
- **Inline argument parsing remains duplicated** — Each script still owns a small local parser; this is acceptable for simple scripts but created the missing-value handling findings below.
- **Analytical defaults unchanged for valid invocations** — Preprocessing branches, clustering grid defaults, PFlog marker heatmap defaults, and summary-grid helper calls remain unchanged for well-formed commands.

## Questions for the analyst

- Should every value-bearing CLI flag fail when present without a value, even when that flag has a documented default if omitted?
- Should script defaults such as `--elbow-n 20`, `--extra-dims 30,50`, and heatmap `--layer pflog` be treated as reportable analysis parameters in the notebook or only implementation defaults?
- Is the main interactive workflow to edit the parameter block in RStudio, pass CLI flags from Bash, or both equally?
- Do you want a tiny shared CLI helper in `scripts/`/`R/`, or should each script keep local parsing as long as the missing-value behavior is consistent?

## Findings

### Statistics & causal inference
#### Major
No separate statistics-only findings. The consequential statistical risk is the same CLI-parameter validation issue listed under Data pipeline & leakage: malformed commands can silently change the preprocessing branch, input object, or PC count.

#### Minor
None.

### Data pipeline & leakage
#### Major
##### F1. Bare `--elbow-n` becomes one PC instead of failing
`scripts/cluster-sobj.R:67-86`
```r
input <- arg("--input")
elbow_n <- as.integer(arg("--elbow-n"))
...
if (
  length(elbow_n) != 1 ||
  is.na(elbow_n) ||
  !is.finite(elbow_n) ||
  elbow_n <= 0
) {
```
**Why it matters here**: `--elbow-n` is the primary PC count selected from elbow diagnostics. If a command accidentally omits the value, `arg("--elbow-n")` returns `TRUE`; `as.integer(TRUE)` becomes `1`; the script can generate plausible `cluster_*_elbow1.rds`, UMAP, and clustree artifacts from one PC instead of failing.
**Fix**: Parse the raw value first, reject `NULL` and `TRUE`, then coerce to integer and keep the existing positive-integer checks.

##### F2. Valueless preprocessing flags silently fall back to defaults
`scripts/preprocess-sobj.R:52-67`
```r
input <- file.path(INPUT_OBJECT_DIR, "pipseq_processed_matrix_with_egfp.rds")
normalization <- "log1p"
...
if (!is.null(arg("--input")) && !identical(arg("--input"), TRUE)) {
  input <- arg("--input")
}
```
**Why it matters here**: `--input` and `--normalization` define the raw Seurat object and normalization branch. A malformed command such as `--input --normalization pflog` or `--normalization --filter-cell-cycle` can run the default raw object or `log1p` branch and still write valid-looking `preprocess_*` artifacts.
**Fix**: Use defaults only when the flag is absent. If a value-bearing flag is present and `arg()` returns `TRUE`, stop with a missing-value error.

##### F3. Batch clustering wrapper ignores missing values for value-bearing flags
`scripts/cluster-all.R:56-75`
```r
elbow_n <- "20"
input_dir <- CURRENT_OBJECT_DIR
...
if (!is.null(arg("--elbow-n")) && !identical(arg("--elbow-n"), TRUE)) {
  elbow_n <- arg("--elbow-n")
}
```
**Why it matters here**: `cluster-all.R` fans out clustering across every preprocessed object. A malformed command such as `--input-dir --elbow-n 30` can ignore the missing `--input-dir` value and run against `CURRENT_OBJECT_DIR`; a bare `--elbow-n` keeps default `20`. This can regenerate a complete clustering grid for unintended inputs or parameters.
**Fix**: Reject present-without-value for `--elbow-n`, `--input-dir`, `--extra-dims`, and `--resolutions` before listing inputs or launching child scripts.

#### Minor
None.

### Bioinformatics
#### Major
No additional bioinformatics-specific findings beyond F1/F2. Those findings matter biologically because PC count, input object, and normalization branch define the scRNA-seq clustering and marker-annotation context.

#### Minor
None.

### LLM coding antipatterns
#### Major
No separate finding beyond F1-F3. These are the LLM-relevant pattern too: hidden default-parameter drift introduced during a refactor that otherwise reads as stylistic.

#### Minor
None.

### Documentation & schema fidelity
#### Major
Covered by F1 and F2: script usage docs require `<positive integer>`, `<raw-seurat-object.rds>`, and `<log1p|pflog>`, but present-without-value invocations are not rejected consistently.

#### Minor
##### F4. Heatmap header omits notebook symlink side effect
`scripts/big-heatmap-plot.R:24-28,262-268`
```r
#   --out-dir
#     Directory for PNG/PDF outputs. Defaults to FIGURE_DIR/annotation.
...
notebook_figure_dir <- here::here("notebook", "figures")
notebook_png_path <- file.path(notebook_figure_dir, basename(png_path))
link_created <- file.symlink(png_path, notebook_png_path)
```
**Why it matters here**: Running the heatmap script updates both Box figure outputs and the notebook figure symlink. A future analyst reading only the script-local docs could miss that rerendering/staleness decisions may depend on this side effect.
**Fix**: Add an `Outputs:` block to the header listing the PNG/PDF under `--out-dir` and the symlink written into `notebook/figures/`.

### Code quality
#### Major
F1-F3 are also the code-quality issue: duplicated local CLI parsing has inconsistent behavior for missing values.

#### Minor
None. The move to `system2()` argument vectors and top-to-bottom script layout is otherwise a maintainability improvement.

## What was checked but is fine

- **Statistics & causal inference**: No changed code alters tests, model choices, pseudobulk statistical unit, or causal claims.
- **Data pipeline & leakage**: `system2()` command construction uses argument vectors; no shell quoting or injection issue was found.
- **Bioinformatics**: The refactor preserves preprocessing branches, optional cell-cycle HVG filtering, Seurat-safe cluster artifact tags, PFlog heatmap layer default, and marker-list handling for valid commands.
- **LLM coding antipatterns**: No hallucinated package APIs or fabricated outputs were found; `here`, `system2()`, and Seurat/ComplexHeatmap calls are real and previously exercised.
- **Documentation & schema fidelity**: `summarize-cluster-grid.R` documents no arguments and now explicitly rejects extra CLI args, so that schema matches behavior.
- **Code quality**: Repeated small `arg()` helpers are acceptable for these scripts if missing-value behavior is made consistent.

## Notes

- The review found one root remediation path: make value-bearing CLI flags distinguish absent from present-without-value everywhere. Fixing that should close F1-F3.
- The diff is intended as a refactor; for well-formed commands, sub-agents did not find evidence that analytical defaults changed.
