# Tripwire run — working tree — 2026-07-03

**Static review**: `.living/outputs/reviews/2026-07-03-working-tree.md`
**Audit/plan**: `.living/outputs/reviews/2026-07-03-working-tree-tripwires.md`
**Runner**: `tools/run-tripwires.R`
**Command**: `Rscript tools/run-tripwires.R`
**Result**: 5 PASS, 0 FAIL, 2 SKIP

## Results

| Status | Slug | Message |
|---|---|---|
| PASS | `cluster-wrapper-contract` | `cluster-all.R` loads package constants in R, discovers preprocess inputs without a shell loop, and exposes a non-executing command preview. |
| PASS | `branch-artifact-collision` | Cluster columns, UMAP reductions, clustree tags, and clustered RDS names use a validated normalization + cell-cycle branch_tag with Seurat-safe characters. |
| PASS | `report-values-freshness` | HTML is newer than QMD and 12 referenced figure target(s); prose states top 20 HVGs and 500 cells. |
| PASS | `missing-counts-file` | `preprocess-sobj.R` returns non-zero for a deliberately missing `--input` path and surfaces the missing-file boundary. |
| PASS | `metadata-contract` | `analysis_labels.yml` declares Mouse, Condition, and derived sample_id; `preprocess-sobj.R` validates and derives them before downstream use. |
| SKIP | `label-permutation` | Static firewall passed: Condition is confined to metadata/sample_id boundaries in early scripts. Full label permutation is skipped because the project has no scratch-output hook for safe blind HVG/PCA/UMAP reruns. |
| SKIP | `toy-contrast-direction` | Contrast direction is encoded in `analysis_labels.yml`; toy known-answer DE run is skipped because no differential-expression entry point exists yet. |

## Notes

- SKIP rows are not failures; they mark scientific boundaries that need future execution hooks or a future DE entry point.
- The freshness tripwire passed only after rerendering `notebook/sc_analysis.qmd`.
- The missing-input tripwire runs the preprocessing script with a deliberately absent input path and treats the expected non-zero child process as a PASS.
- The cluster wrapper tripwire catches the previous shell-loop / exported-`CURRENT_OBJECT_DIR` failure mode.
- The branch artifact tripwire catches hyphenated Seurat reduction tags before `DimPlot()` can look up a name Seurat rewrote.
