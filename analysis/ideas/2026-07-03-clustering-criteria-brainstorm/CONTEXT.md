# ESPI Clustering Criteria Ideation Context

**Date**: 2026-07-03
**Session**: clustering-criteria-brainstorm
**Focus**: Choose defensible, label-blind criteria for selecting normalization method, number of PCs, and Leiden clustering resolution for downstream scRNA-seq analysis.

## Project

ESPI is a minimal R package plus executable scientific analysis pipeline for a PipSeq V T2 single-cell RNA-seq experiment using a custom reference including eGFP. The biological contrast is `p27CKO` vs `p27CKO +EStim`.

Experimental constraints from the repo context:

- Conditions: `p27CKO` vs `p27CKO +EStim`.
- Replicates: six Mouse × Condition pseudobulk samples; mice 10 and 3 are paired, mouse 30 is E-Stim only, mouse 33 is control only.
- BrdU: added at E-Stim + 24 h, washed after 48 h, fixed at 5 d post E-Stim.
- Primary condition-level differential expression unit: Mouse × Condition pseudobulk sample, not cell.
- Early preprocessing/clustering parameter selection must be label-blind: do not use `Condition`, E-Stim separation, BrdU effects, or desired biological effects to choose normalization, PC count, or clustering resolution.

## Data and outputs

Large analysis inputs and outputs live outside the repo under:

```text
/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data
```

Current relevant outputs:

```text
/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data/tables/cluster/cluster_grid_summary.tsv
/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data/figures/cluster/cluster_grid_clustree_12_panel.png
/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data/figures/cluster/umap_resolution_sweep_pflog_filter_cc_dims30.png
```

Notebook symlinks:

```text
notebook/figures/cluster_grid_clustree_12_panel.png
notebook/figures/umap_resolution_sweep_pflog_filter_cc_dims30.png
```

## Existing pipeline

Manifests list these active analyses and methods:

- Preprocessing: `scripts/preprocess-sobj.R`, `scripts/preprocess-all.R`
  - Produces normalization and cell-cycle-HVG-filtering branches plus QC/HVG/PCA diagnostics.
- Candidate clustering: `scripts/cluster-sobj.R`, `scripts/cluster-all.R`, `scripts/summarize-cluster-grid.R`
  - Produces UMAP/clustree candidate clustering outputs, a 36-row clustering grid summary, a 12-panel clustree grid, and a representative UMAP resolution sweep.
- Single-cell notebook: `notebook/sc_analysis.qmd`
  - Currently documents the clustering grid and supplemental figures.
- Tripwires: `tools/run-tripwires.R`, `analysis_labels.yml`
  - Checks cluster wrapper execution, branch artifact separation, report freshness, missing-input failure, metadata contract, label firewall, and future contrast direction.

Reusable R helpers:

- `R/dim-reduction.R`: PCA branch helpers.
- `R/preprocess-plots.R`: QC/HVG/PCA diagnostics.
- `R/cluster-plots.R`: UMAP overlays and clustree plots.
- `R/cluster-grid-summary.R`: Writes grid summaries with ARI/Jaccard to the current reference and saves multi-panel diagnostics.

## Parameter grid

Current grid:

```text
2 normalizations × 2 CC-HVG policies × 3 PC counts × 3 resolutions = 36 configurations
```

- Normalization: `log1p`, `pflog`.
- Cell-cycle-HVG policy: retained / filtered.
- PCs: 20, 30, 50.
- Leiden resolutions: 0.3, 0.5, 0.8.

Current reference clustering column:

```text
cluster_pflog_filter_cc_dims30_res0.3
```

Reference details:

```text
normalization = pflog
filtered_cell_cycle = TRUE
dims = 30
resolution = 0.3
branch_tag = pflog_filter_cc
n_cells = 5538
n_clusters = 11
```

## Current grid summary signals

The summary table has 36 rows and includes:

- `n_clusters`
- `min_cluster_n`
- `n_small_clusters`
- `fraction_cells_in_small_clusters`
- `ari_vs_reference`
- `mean_best_jaccard_to_reference`
- `min_best_jaccard_to_reference`

Observed numeric patterns from `cluster_grid_summary.tsv`:

- PFlog branches are closer to the current reference than log1p branches.
  - PFlog mean ARI across grid rows: ~0.776; minimum ARI: ~0.608.
  - log1p mean ARI: ~0.560; minimum ARI: ~0.372.
- Resolution 0.3 gives the highest agreement with the current reference and least fragmentation.
  - At 20 PCs / res 0.3: all four branches have 10 clusters; mean ARI ~0.823.
  - At 30 PCs / res 0.3: all four branches have 11 clusters; mean ARI ~0.799.
  - At 50 PCs / res 0.3: all four branches have 11 clusters; mean ARI ~0.818.
- Higher resolution increases fragmentation.
  - Res 0.5 gives 12–14 clusters depending on PCs/branch.
  - Res 0.8 gives 13–17 clusters and a larger small-cluster burden in several branches.
- Small-cluster burden is low overall but rises at high resolution.
  - Maximum fraction of cells in clusters smaller than 50 cells is ~1.45%.
- Top non-reference ARI rows are PFlog-based.
  - `cluster_pflog_no_filter_cc_dims50_res0.3`: 11 clusters, ARI ~0.948, minimum best Jaccard ~0.787.
  - `cluster_pflog_filter_cc_dims50_res0.3`: 11 clusters, ARI ~0.911, minimum best Jaccard ~0.669.
  - `cluster_pflog_filter_cc_dims20_res0.3`: 10 clusters, ARI ~0.897, minimum best Jaccard ~0.595.

## Visual diagnostics already reviewed

The 12-panel clustree grid spans:

- Rows: 20, 30, 50 PCs.
- Columns: log1p retained, log1p filtered, PFlog retained, PFlog filtered.
- Each panel shows the 0.3 → 0.5 → 0.8 resolution sweep.
- Node position is only a graph-layout aid; edge structure, node size, and split pattern are the interpretable parts.

The representative UMAP sweep uses PFlog, CC-HVG filtered, 30 PCs:

- Res 0.3: 11 clusters; coarse structure with a few small isolated clusters.
- Res 0.5: 13 clusters; subdivides main manifolds.
- Res 0.8: 15 clusters; further subdivisions and small clusters.

## Current criterion problem

The next decision is not “which clustering looks best by condition separation.” The goal is to define criteria that make a defensible, label-blind choice among normalization, PC count, and resolution before marker annotation and pseudobulk DE.

Strong candidate criteria should include:

1. Stability across adjacent PC counts and preprocessing branches.
2. Low fragmentation: avoid many tiny clusters unless marker coherence later justifies them.
3. Interpretable cluster-size distribution.
4. Agreement metrics: ARI and best-overlap Jaccard to a stated reference, while acknowledging that any reference choice can bias interpretation.
5. Consistent broad structure across CC-HVG retained vs filtered branches.
6. Biological marker coherence later, but only after a label-blind clustering choice or as a pre-declared tie-breaker that does not use condition labels.
7. Reproducibility and transparency: chosen criteria must be documented in the notebook and `.living/decisions.md`.

## Requested ideation output

Generate 2 concrete ideas from the assigned persona for how ESPI should choose or justify criteria for selecting:

- normalization (`log1p` vs `pflog`),
- number of PCs (20 vs 30 vs 50), and
- Leiden resolution (0.3 vs 0.5 vs 0.8).

Ideas should be grounded in the existing 36-row table, clustree grid, UMAP sweep, and future marker/coherence checks. Do not use condition labels or desired E-Stim biology as selection criteria.
