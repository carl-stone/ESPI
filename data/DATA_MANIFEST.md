# Data Manifest

This repo stores small R package data in `data/`. Large analysis inputs and outputs live outside the repo in Box Drive.

| Dataset | Location | Type | Status | Notes |
|---------|----------|------|--------|-------|
| Mouse cell-cycle genes | `data/mouse_cell_cycle_genes.rda` | R package data | active | Mouse mapping used to remove known Seurat cell-cycle HVGs in preprocessing branches. |
| Cell type marker genes | `data/cell_type_marker_genes.rda`, `data/cell_type_marker_labels.rda` | R package data | active | Curated marker-gene lists and display labels for broad retinal cell type annotation. Generated from `data-raw/cell-type-marker-genes.R`; marker choices were provided by Ed and Megan from domain knowledge/literature review, with no marker-by-marker rationale recorded in the repo. |
| ESPI Seurat objects and figures | `/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data` | external Box data root | active external | Contains current Seurat objects, preprocessing figures, clustering figures, and derived artifacts. Do not add fallback paths. |

## Mycelium Data Directories

- `data/raw/`: placeholder for future tracked metadata or small immutable inputs.
- `data/processed/`: placeholder for future small derived data.
- `data/metadata/`: placeholder for schemas/provenance notes.

Do not commit large Box data or generated figure outputs here.
