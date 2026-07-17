# Data Manifest

This repo stores small R package data in `data/`. Large analysis inputs and outputs live outside the repo in Box Drive.

| Dataset | Location | Type | Status | Notes |
|---------|----------|------|--------|-------|
| Mouse cell-cycle genes | `data/mouse_cell_cycle_genes.rda` | R package data | active | Mouse mapping used to remove known Seurat cell-cycle HVGs in preprocessing branches. |
| Cell type marker genes | `data/cell_type_marker_genes.rda`, `data/cell_type_marker_labels.rda` | R package data | active | Curated marker-gene lists and display labels for broad retinal cell type annotation. Generated from `data-raw/cell-type-marker-genes.R`; marker choices were provided by Ed and Megan from domain knowledge/literature review, with no marker-by-marker rationale recorded in the repo. |
| UMAP feature gene list | `data/umap_feature_list.rda` | R package data | active | Selected 3 × 3 manuscript feature list for MG-selected UMAP expression plots. Generated from `data-raw/umap-feature-list.R`. |
| ESPI Seurat objects and figures | Resolved by `R/config.R`; default: `/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data` | external Box data root | active external | Contains current Seurat objects, preprocessing figures, clustering figures, and derived artifacts. An untracked `config.local.R` may define `MEGAN_SC_DATA_DIR` directly or `BOX_PATH` as its parent; loading fails if the resolved root does not exist. |

## Mycelium Data Directories

- `data/raw/`: placeholder for future tracked metadata or small immutable inputs.
- `data/processed/`: placeholder for future small derived data.
- `data/metadata/`: placeholder for schemas/provenance notes.

Do not commit large Box data or generated figure outputs here.
