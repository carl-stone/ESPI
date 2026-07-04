# Review — working tree — 2026-07-04

**Scope**: Working tree (`git diff HEAD` plus untracked files)
**Files reviewed**: 14
**Sub-agents run**: 6 (`stats-causal`, `data-pipeline-leakage`, `bioinformatics`, `llm-failure-modes`, `doc-schema-fidelity`, `code-quality`)

## Key decisions in this analysis

- **DD statistical kernel** — Differential detection now uses `muscat::pbDS(method = "DD")` / `edgeR_NB_optim` on detected-cell pseudobulk counts with a cellular-detection-rate offset and edgeR quasi-likelihood robust dispersion. See F2, F5.
- **DD gene universe** — The DD tested-gene universe is defined by muscat's internal 90%-detection filter, not by the DESeq2 tested-gene set.
- **Primary condition design** — Primary DE and DD use `~ condition` across all six Mouse × Condition pseudobulk samples.
- **Paired sensitivity design** — Paired DE and DD use `~ mouse + condition` only on mice 10 and 3 and are interpreted as within-paired-mice sensitivity results.
- **Report-number provenance** — `numbers.json`, `design_summary.tsv`, and the notebook prose carry the DD tested-gene and hit counts. See F4.
- **Cluster-marker narrative state** — The notebook introduces a `FindAllMarkers()` cluster-marker section before a working marker script or marker output exists. See F1.

## Questions for the analyst

- Should the cluster-marker section be committed now as planned work, or should it be removed until `FindAllMarkers()` outputs exist?
- If cluster markers are added, should their p-values be shown at all, or should the report restrict them to descriptive marker rankings for annotation?
- Is the primary DD claim meant to support a paper result, or only to rule out obvious detection-fraction artifacts in this branch?
- Do you want report freshness to require every quoted headline count, including “zero primary DD hits,” to live in `numbers.json`?
- Should the DD method label be a stable provenance constant now, or is this script unlikely to change methods again before submission?

## Findings

### Statistics & causal inference

#### Major

No major findings.

#### Minor

##### F2. DD prose overstates what the CDR-offset model establishes
`notebook/sc_analysis.qmd:194-197`
```r
offset, and fit edgeR quasi-likelihood models with robust dispersion. This
adjusts for per-cell library depth and per-sample detection intensity so DD hits
reflect changes in the fraction of expressing cells rather than PipSeq depth
variation. Primary DD uses `~ condition` across all six pseudobulk samples;
```
**Why it matters here**: The muscat offset reduces depth-driven detection artifacts, but six Mouse × Condition samples with partial pairing cannot prove every DD signal is free of residual technical or sample-level confounding.
**Fix**: Soften to “intended to reduce PipSeq depth-driven detection artifacts and estimate condition-associated changes in detected-cell fraction.”

### Data pipeline & leakage

#### Major

No major findings.

#### Minor

##### F3. Muscat DD helper lacks a local counts/metadata alignment assertion
`scripts/run-mg-selected-de.R:535-552`
```r
sample_ids <- rownames(sample_table)
keep_cells <- cell_metadata$pseudobulk_sample_id %in% sample_ids
cell_metadata <- cell_metadata[keep_cells, , drop = FALSE]
counts <- counts[, keep_cells, drop = FALSE]
...
rownames(col_data) <- colnames(counts)
```
**Why it matters here**: The current script checks global Seurat metadata/count alignment before calling this helper, but the helper itself assigns muscat sample labels by position. A future caller with reordered metadata could silently assign cells to the wrong Mouse × Condition sample.
**Fix**: Add `identical(colnames(counts), rownames(cell_metadata))` checks before and after filtering, or subset/reorder `counts` by metadata rownames.

##### F4. Primary DD total-hit count is quoted but not registered
`scripts/run-mg-selected-de.R:1301-1325`
```r
reportable_values <- list(
  n_samples = nrow(sample_table),
  n_cells = ncol(counts),
  n_tested_genes = length(primary_de$tested_genes),
  n_detection_tested_genes = nrow(full_detection),
```
**Why it matters here**: The notebook says primary muscat DD found no FDR-significant genes. `numbers.json` registers primary DD tested genes and marker hits, but not total primary DD hits, so a report-values freshness check cannot guard that quoted zero.
**Fix**: Add `n_detection_hits <- sum(!is.na(full_detection$padj) & full_detection$padj < 0.05)` and register it in `reportable_values`.

### Bioinformatics

#### Major

##### F1. Notebook presents an ungenerated `FindAllMarkers()` analysis as completed evidence
`notebook/sc_analysis.qmd:166-170`; `scripts/find-markers-mg-selected.R:1-5`
```r
### Cluster marker genes

TODO before this: hand-combine clusters if they look the same before FindAllMarkers().

Next, moving forward with the PFlog, 30-PC, resolution-0.3 dataset, I used FindAllMarkers() to identify uniquely differentially expressed genes in each cluster.
```
```r
# The purpose of this script is to run FindAllMarkers on the selected clusters of the mg-selected datasets (FindAllMarkers doesn't use HVGs so doesn't matter if it's cc filtered or not) with selected number of PCs and clustering resolution.

# The script should save a .csv file that's the FindAllMarkers results object, and it should save dot plots of the top ~some number~ markers for each cluster. (like 5 each)
```
**Why it matters here**: The notebook is the scientific narrative, but the only matching uncommitted script is a comment stub. Also, `FindAllMarkers()` on cells reused for clustering is exploratory and cell-level; its p-values should not be treated as confirmatory Mouse × Condition evidence.
**Fix**: Remove the section until marker outputs exist, or reframe it as planned/descriptive marker ranking and add the actual sample-aware validation or explicit non-confirmatory warning before reporting p-values.

#### Minor

No minor findings beyond F1.

### LLM coding antipatterns

#### Major

See F1. This is also an LLM-failure-mode issue: ungrounded report prose states “I used FindAllMarkers()” before the corresponding analysis artifact exists.

#### Minor

No minor findings.

### Documentation & schema fidelity

#### Major

See F1.

#### Minor

No additional documentation findings. The DD method prose and counts match the regenerated `numbers.json` values reviewed by the agents.

### Code quality

#### Major

No major findings.

#### Minor

##### F5. DD method provenance label is hard-coded in three places
`scripts/run-mg-selected-de.R:616-618,1226-1230,1317-1318`
```r
contrast = CONTRAST_DIRECTION,
design = design_label,
method = "muscat_edgeR_NB_optim",
```
```r
method = c(
  "deseq2_wald",
  "muscat_edgeR_NB_optim",
```
```r
dd_method = "muscat_edgeR_NB_optim",
```
**Why it matters here**: The method string is provenance in DD result rows, `design_summary.tsv`, and `numbers.json`. If the label changes later, one missed occurrence would make downstream reporting inconsistent.
**Fix**: Define one script-level `DD_METHOD <- "muscat_edgeR_NB_optim"` near the other constants and reuse it everywhere.

## What was checked but is fine

- **Statistics & causal inference**: The primary `~ condition` and paired `~ mouse + condition` split matches the documented Mouse × Condition design decision; paired-only hits are framed as sensitivity results.
- **Data pipeline & leakage**: Current primary and paired DD call paths preserve count/metadata alignment and validate Mouse × Condition sample construction before model fitting.
- **Bioinformatics**: The muscat DD call, DD tested-gene universe, marker-overlap universe, and paired-sensitivity framing match the recorded decisions.
- **LLM coding antipatterns**: The muscat `min_cells = 0L`, `filter = "none"`, and `verbose = FALSE` arguments match the approved muscat DD workflow rather than smuggled defaults.
- **Documentation & schema fidelity**: `detection_full_results.tsv` and `design_summary.tsv` schemas match the changed script's documented DD output structure.
- **Code quality**: The new `run_detection_muscat_dd()` helper has two real callers and mirrors the existing `run_deseq2()` helper pattern; it is not a speculative abstraction.

## Notes

- F1 is cross-cutting: one remediation should fix the bioinformatics, LLM-failure, documentation-fidelity, and code-quality versions of the same issue.
- Suggested tripwires emitted by reviewers: `report-value-freshness` for F1/F4 and `shuffled-sample-order` for F3.
