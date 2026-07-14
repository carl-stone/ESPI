# Environments & Installations

## Primary Environment

- **Project type**: Minimal R package for ESPI single-cell analysis.
- **Command runner**: `just` from the repo root; run `just --list` to
  discover recipes.
- **Package load**: `just load` (`devtools::load_all()`).
- **Formatter**: `just format <paths>` for R files.
- **Notebook renderer**: `just notebook`.
- **External data root**: `/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data`.

## Setup From Scratch

Install the R package dependencies, Quarto, Air, and `just` (for example,
`brew install just`). Then load the package from the repo root:

```sh
just load
```

Common package and document tasks:

```sh
just document
just readme
just format R/<file>.R scripts/<script>.R
just notebook
```

## Routine Pipeline

Use the canonical wrapper for a complete analysis run:

```sh
just run [source] [overwrite]
just run-dry-run [source] [overwrite]
```

`source` defaults to `counts-qc` and accepts `counts-qc`, `legacy`, or a
quoted explicit RDS path. `overwrite` defaults to `false`; set it to `true`
only when intentionally replacing protected marker or DE outputs. The
dry-run prints the deterministic plan without changing files. A full run
validates each stage, renders the notebook, and then runs tripwires.

`counts-qc` runs count processing and QC before preprocessing the source.
`legacy` uses the existing
`INPUT_OBJECT_DIR/pipseq_processed_matrix_with_egfp.rds`. An explicit RDS
path bypasses named-source selection and uses that object as the source.

The current selected branches are:

- **Source**: PFlog, no-cell-cycle-HVG (no-filter-CC), 30 PCs, resolution
  0.3. The source has 4,146 cells and 9 clusters; exclude source clusters
  2, 7, and 8, retaining 3,456 cells.
- **MG-selected**: reselect HVGs and compute a 50-PC PCA so 20-, 30-, and
  50-PC candidates remain available. The selected branch is PFlog,
  no-filter-CC, 20 PCs, resolution 0.5, with 8 clusters.

Current cluster columns are
`cluster_pflog_no_filter_cc_dims30_res0.3` for the source and
`cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5` for MG-selected.

## Expert Recovery and Maintenance

Use the commands in this section only to recover a failed stage, inspect
intermediate outputs, or replace a single artifact. They are not the routine
workflow.

### Raw counts and QC

Create the combined raw Seurat object from the six 10X directories and sample
metadata:

```sh
Rscript scripts/01-process-counts.R
```

This writes `data/input/sobj_raw.rds` beneath the Box data root
(`DATA_ROOT_DIR/data/input/sobj_raw.rds`).

Run QC filtering:

```sh
Rscript scripts/02-qc-filtering.R
```

The script writes QC figures to `FIGURE_DIR/qc/*.png`, QC tables to
`TABLE_DIR/qc/*.tsv`, the annotated raw object to
`INPUT_OBJECT_DIR/sobj_raw_with_qc.rds`, and the filtered object to
`INPUT_OBJECT_DIR/sobj_qc_filtered.rds`. `percent.mt` uses all 37 observed
mitochondrial features, whose labels are mixed (including `mt-Rnr1`,
`mt-Rnr2`, and non-`mt-` labels). `DropletUtils::emptyDrops()` supplies the
cell-call FDR and `is_cell` flag. `scDblFinder` runs per sample on called
barcodes above the count and feature floors and supplies doublet scores and
singlet/doublet calls. Called singlets define sample-specific lower three-MAD
thresholds for log10 counts and detected features and an upper three-MAD
threshold for mitochondrial percentage. The filtered object contains the
4,146 cells that pass those three MAD criteria; `is_cell` and `is_singlet`
remain separate from `pass_qc`. `percent.ribo` remains diagnostic only.

### Low-level preprocessing

The low-level `preprocess` recipe intentionally defaults to the historical
`legacy` source, unlike canonical `just run`, whose default is `counts-qc`:

```sh
just preprocess
just preprocess counts-qc
just preprocess-one pflog false counts-qc
```

Pass an explicit object when recovering another input:

```sh
Rscript scripts/03-preprocess-all.R --input-source legacy
Rscript scripts/03-preprocess.R --input /path/to/object.rds --normalization pflog
```

Preprocessing replaces the current branch artifacts, so downstream recovery
commands consume the source selected for that run.

### Clustering and summaries

```sh
just cluster-dry-run
just cluster
just summarize-clusters
just summarize-mg-selected 20
```

Clustered outputs use underscore branch tags such as `pflog_no_filter_cc`
because Seurat rewrites hyphens in reduction names. The wrapper loads ESPI
path constants in R; do not rely on a shell-exported `CURRENT_OBJECT_DIR`.

For direct figure or marker recovery, pass the current dimensions, resolution,
and object explicitly:

```sh
just marker-heatmap 30 0.3 /path/to/source-clustered-object.rds
just marker-heatmap 20 0.5 /path/to/mg-selected-object.rds
just cluster-marker-heatmaps 30 0.3 /path/to/source-clustered-object.rds
just cluster-marker-heatmaps 20 0.5 /path/to/mg-selected-object.rds
just mg-markers
just mg-figures
just mg-de
```

### Notebook and tripwires

Render the notebook or run the lightweight scientific-boundary checks
separately when recovering those final stages:

```sh
just notebook
just tripwires
```

The tripwire runner uses `analysis_labels.yml` for label and contrast
declarations and exits non-zero for `FAIL` rows. It exercises missing input
and metadata failures, report freshness, provenance semantics, contrast
direction, and a guarded scratch-only label permutation. The permutation runs
after `sample_id` derivation, preserves sample identity, skips figure writes,
and stops at `blind_qc_complete`.

## Required R Packages

Declared in `DESCRIPTION`:

- `clustree`
- `DropletUtils`
- `gghalves` from `erocoar/gghalves`
- `ggplot2`
- `ggraph`
- `ggrepel`
- `here`
- `igraph`
- `mclust`
- `patchwork`
- `scclrR` from `cleartools/scclrR`
- `Seurat`
- `SeuratObject`
- `SingleCellExperiment`

Suggested:

- `biomaRt`
- `devtools`

## Data and Output Paths

The analysis expects Box Drive data at:

```text
/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data
```

Do not add fallback data paths. Missing Box paths should fail explicitly.

Important subdirectories:

```text
seurat_objects/current/
figures/preprocess/
figures/cluster/
```

Notebook figures should be symlinked into `notebook/figures/` and referenced
with notebook-relative paths.

## Mycelium Local Hooks

Mycelium created local Claude Code hooks in `.claude/settings.local.json`. The
file points at the installed plugin cache:

```text
/Users/carlstone/.omp/plugins/cache/plugins/mycelium___mycelium___0.0.0/
```

`.claude/` is gitignored because those paths are local and can break after
plugin upgrades or cache cleanup. On a fresh clone or after upgrading the
plugin, rerun Mycelium initialization or hook setup before expecting hooks to
fire.
