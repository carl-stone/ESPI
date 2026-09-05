# Data and outputs

This repo stores small R package data in `data/`. Large analysis inputs and outputs live outside the repo in Box Drive.

Paths in this document are relative to the repository root unless absolute.

| Dataset | Location | Type | Status | Notes |
|---------|----------|------|--------|-------|
| Mouse cell-cycle genes | `data/mouse_cell_cycle_genes.rda` | R package data | active | Mouse mapping used to remove known Seurat cell-cycle HVGs in preprocessing branches. |
| Cell type marker genes | `data/cell_type_marker_genes.rda`, `data/cell_type_marker_labels.rda` | R package data | active | Curated marker-gene lists and display labels for broad retinal cell type annotation. Generated from `data-raw/cell-type-marker-genes.R`; marker choices were provided by Ed and Megan from domain knowledge/literature review, with no marker-by-marker rationale recorded in the repo. |
| UMAP feature gene list | `data/umap_feature_list.rda` | R package data | active | Selected manuscript feature list for MG-selected UMAP expression plots. Generated from `data-raw/umap-feature-list.R`. |
| Volcano gene list | `data/volcano_genes.rda` | R package data | available | Gene vector saved by `data-raw/volcano_genes.R`; consult the plotting script for whether it uses this list. |
| ESPI Seurat objects and figures | Resolved by `R/config.R`; default: `/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data` | external Box data root | active external | Contains Seurat objects, preprocessing figures, clustering figures, and derived artifacts. An untracked `config.local.R` may define `MEGAN_SC_DATA_DIR` directly or `BOX_PATH` as its parent. |

Do not commit large Box data or generated figure outputs here.

## Repo-local figures and results

- `notebook/figures/`: report figures copied by the scripts, with primary pipeline
  outputs retained in Box. Grid stages copy selected views, not all candidates;
  there is no automatic selection based on notebook text.
- `notebook/figures/custom-markers/`: custom-marker manuscript plots and retained
  candidate exports. `scripts/custom-marker-plots.R` writes the per-cell heatmap
  and dot plot here; the cluster-average pair is a retained candidate, not an
  output of the current script.
- `exploratory/neurogenic-proportion/outputs/`: preserved threshold-grid results
  and bootstrap diagnostics, with the original run settings unchanged.
- `exploratory/plotting/outputs/`: prior exports moved from the repo root; see
  the adjacent notes for their provenance limits.
