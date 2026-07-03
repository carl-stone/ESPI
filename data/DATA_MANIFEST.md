# Data Manifest

This repo stores small R package data in `data/`. Large analysis inputs and outputs live outside the repo in Box Drive.

| Dataset | Location | Type | Status | Notes |
|---------|----------|------|--------|-------|
| Mouse cell-cycle genes | `data/mouse_cell_cycle_genes.rda` | R package data | active | Mouse mapping used to remove known Seurat cell-cycle HVGs in preprocessing branches. |
| ESPI Seurat objects and figures | `/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data` | external Box data root | active external | Contains current Seurat objects, preprocessing figures, clustering figures, and derived artifacts. Do not add fallback paths. |

## Mycelium Data Directories

- `data/raw/`: placeholder for future tracked metadata or small immutable inputs.
- `data/processed/`: placeholder for future small derived data.
- `data/metadata/`: placeholder for schemas/provenance notes.

Do not commit large Box data or generated figure outputs here.
