# Review — Working tree MG module/p27 heatmaps — 2026-07-05

**Scope**: Working tree against `HEAD`, including untracked files  
**Files reviewed**: 19  
**Sub-agents run**: 6 — stats-causal, data-pipeline-leakage, bioinformatics, llm-failure-modes, doc-schema-fidelity, code-quality

## Key decisions in this analysis

The consequential analytical choices in this work. Some have associated findings below; others are informational so you can decide whether to revisit them.

- **p27 enrichment estimand** — The p27 strip is a cluster-level descriptive enrichment score, not a condition-level Mouse × Condition inference. See F2.
- **Permutation null** — The p27 null permutes cluster labels within each Mouse × Condition sample, preserving sample-local p27 distributions and cluster-size composition.
- **Expression layers** — The p27 strip defaults to RNA `pflog`, while module scoring defaults to RNA `data`. See F3.
- **Reusable-statistics boundary** — Module aggregation and p27 permutation live in exported R helpers; ComplexHeatmap layout and artifact writing live in `scripts/plot-cluster-marker-heatmaps.R`. See F1.
- **Report integration** — The notebook embeds generated PNG symlinks and must be rerendered after figure regeneration.
- **Artifact provenance** — Mycelium session logs and registry rows are intended to record generated heatmap artifacts and validation evidence. See F5.

## Questions for the analyst

Things the diff alone cannot settle, whose answers change which cleanup matters most.

- Should the module-score body intentionally use Seurat `data` while the p27 strip uses `pflog`, or should both use PFlog for this combined figure?
- Is the p27 strip meant only as a visual annotation of selected clusters, or should it support a manuscript-facing statistical statement?
- Should exported analysis helpers preserve the caller RNG state as a package convention, even when the immediate script has no later randomized step?
- Do you want environment docs to use machine-local absolute Box paths, or should every reproducibility command resolve paths through ESPI constants?
- Should hook-generated in-progress Mycelium logs be cleaned before commit, or kept as evidence of review/session activity?

## Findings

### Statistics & causal inference

#### Major

No major findings.

#### Minor

##### F1. Exported p27 helper resets the caller RNG stream

`R/cluster-marker-heatmap.R:270-276`
```r
  set.seed(seed)
  for (perm_idx in seq_len(n_perm)) {
    permuted_cluster <- cluster
    for (idx in sample_indices) {
      permuted_cluster[idx] <- sample(
```
**Why it matters here**: `compute_cluster_p27_enrichment()` is exported package code. The current plotting script does not run a later random step, but a future analysis that calls this helper before another stochastic step would make downstream randomness depend on call order.
**Fix**: Save and restore `.Random.seed`, or use a local seed context so the helper is deterministic without mutating the caller's RNG stream.

##### F2. Notebook prose does not state the p27 strip is descriptive

`notebook/sc_analysis.qmd:144`
```markdown
This heatmap summarizes the same selection signals per source cluster. The body z-scores each cell type's `cell_type_marker_genes` module score across clusters; the top strip is the p27 (`Cdkn1b`) enrichment z-score whose null permutes cluster labels within each Mouse × Condition sample, so the null preserves per-sample composition. It visualizes the microglia, cell-cycle, and p27 signals behind the cluster exclusions described above.
```
**Why it matters here**: The p27 strip is computed from cell-level cluster labels within samples; it is useful annotation, but it is not the Mouse × Condition condition-level inference used elsewhere in the notebook. The same line also drifts from `CONTEXT.md`, which says to use `p27` in prose and reserve gene naming for `CDKN1B`.
**Fix**: Reword the paragraph to call the strip a descriptive cluster-level p27 annotation and use glossary-approved p27/CDKN1B terminology.

### Data pipeline & leakage

#### Major

No major findings.

#### Minor

No minor findings.

### Bioinformatics

#### Major

No major findings.

#### Minor

##### F3. Combined heatmap mixes module and p27 expression layers without rationale

`scripts/plot-cluster-marker-heatmaps.R:124-128`
```r
expression_layer <- get_arg(cli_args, "--layer", "pflog")
if (!nzchar(expression_layer)) {
  stop("--layer must not be empty.", call. = FALSE)
}
score_slot <- get_arg(cli_args, "--slot", "data")
```
**Why it matters here**: The figure places marker-module rows and the p27 strip in one visual summary, but their defaults use different RNA layers. That may be intentional because `AddModuleScore()` follows the existing selection script, but the notebook/caption does not tell readers that the body and strip use different expression scales.
**Fix**: Either run the module scores on `pflog` for this figure or document that module scores use Seurat `data` while p27 uses `pflog`, with a one-sentence rationale.

### LLM coding antipatterns

#### Major

No major findings.

#### Minor

No minor findings.

### Documentation & schema fidelity

#### Major

No major findings.

#### Minor

##### F4. Environment docs hard-code a workstation-specific Box path

`ENVIRONMENTS_INSTALLATIONS.md:67-70`
```sh
Rscript scripts/plot-cluster-marker-heatmaps.R --dims 50 --resolution 0.3
Rscript scripts/plot-cluster-marker-heatmaps.R \
  --input /Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data/seurat_objects/current/cluster_pflog_mg_selected_no_filter_cc_elbow20.rds \
  --dims 30 --resolution 0.3
```
**Why it matters here**: ESPI intentionally uses a local Box root, but this command only works on this workstation and bypasses the package constants the scripts otherwise use. A future clone following the environment guide would hit a path error before testing the figure script.
**Fix**: Show an `Rscript -e 'cat(file.path(ESPI:::CURRENT_OBJECT_DIR, ...))'` form or a clear `<CURRENT_OBJECT_DIR>/...` placeholder instead of the personal absolute path.

##### F5. Mycelium provenance row was clobbered back to placeholder content

`.living/log/LOG_REGISTRY.md:29`
```markdown
| 2026-07-05 | 2026-07-05-007 | espi | main | 11m | 10 | espi_cluster_marker_heatmap_check.R,ANALYSIS_MANIFEST.md,ENVIRONMENTS_INSTALLATIONS.md (+6 more) | | complete | | [log](2026-07-05-007-espi.md) |
```
**Why it matters here**: The registry no longer names the heatmap outputs, validation, or tags, and the linked log also has an appended auto-summary listing `/tmp/espi_cluster_marker_heatmap_check.R` and `local://mg-module-p27-heatmaps-plan.md` as modified files. That makes the audit trail less reliable for the generated analysis artifacts.
**Fix**: Restore the semantic registry row/log summary and remove or quarantine hook-generated in-progress review logs before commit.

### Code quality

#### Major

No major findings.

#### Minor

No additional minor findings beyond F1 and F4.

## What was checked but is fine

- **Statistics & causal inference**: The permutation construction itself preserves Mouse × Condition groups, and the diff does not add new p-value or multiple-comparison claims.
- **Data pipeline & leakage**: Cell/expression alignment uses `colnames(sobj)`, p27 rows are matched to module-score columns, and generated notebook symlinks were verified by the implementation step.
- **Bioinformatics**: The review found no gene-list corruption, missing p27 feature, or sample-vs-cell condition-effect claim in the new helper/script code.
- **LLM coding antipatterns**: No hallucinated APIs, broad silent fallbacks, fabricated verification output, or plan-drift implementation were identified.
- **Documentation & schema fidelity**: Generated Rd exports and the TODO completion match the implemented helper/script work; README.Rmd does not enumerate these scripts.
- **Code quality**: The helper/script split follows the recorded decision; one-off ComplexHeatmap layout remains in the script rather than becoming package API.

## Notes

- The pipeline-leakage reviewer flagged `module_z[is.na(module_z)] <- 0`; synthesis did not include it because this behavior was in the approved plan, matches the existing marker-heatmap convention, and the generated module-score TSVs currently have no zero-variance cell-type rows.
- `notebook/sc_analysis.html` is large because `embed-resources: true` embeds figure bytes; this is expected for this repo.

## Resolution

- F1 fixed: `compute_cluster_p27_enrichment()` now saves/restores caller RNG state while keeping seeded permutations deterministic.
- F2 fixed: notebook prose now frames the p27 strip as a descriptive cluster-level annotation.
- F3 fixed: script and notebook text document the intentional layer split: module scores use RNA `data`; p27 uses RNA `pflog`.
- F4 fixed: environment docs resolve the MG-selected input through ESPI package constants instead of a personal absolute path.
- F5 fixed: review-scoped `2026-07-05-007` and `2026-07-05-008` registry rows/logs carry semantic provenance; older historical blank rows remain out of scope.
- Guarded by `Rscript tools/run-tripwires.R`: heatmap missing input, RNG preservation, and review-scoped Mycelium provenance all pass.
