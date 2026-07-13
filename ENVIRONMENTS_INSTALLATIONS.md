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

## Preprocessing Pipeline

First create the combined raw Seurat object from the six 10X directories and
sample metadata:

```sh
Rscript scripts/process-counts.R
```

This writes `data/input/sobj_raw.rds` beneath the Box data root
(`DATA_ROOT_DIR/data/input/sobj_raw.rds`).

Run QC filtering before preprocessing:

```sh
Rscript scripts/qc-filtering.R
```

The script writes QC figures to `FIGURE_DIR/qc/*.png`, QC tables to
`TABLE_DIR/qc/*.tsv`, and `INPUT_OBJECT_DIR/sobj_qc_filtered.rds` under the
Box data root. `percent.mt` uses all 37 observed mitochondrial features,
whose labels are mixed (including `mt-Rnr1`, `mt-Rnr2`, and non-`mt-` labels);
none were lost upstream or from the custom reference. It retains 22,248 of
983,903 cells across S2, S3, S4, S5, S7, and S8 only when `nFeature_RNA >= 50`,
`nCount_RNA >= 100`, and `percent.mt <= 20`. This is a data-specific cutoff,
not a universal PipSeq rule: among the 22,751 cells meeting the complexity
thresholds, complete-mitochondrial P95/P97.5/P99 were
16.038/19.313/27.666%, and the >20% sparse extreme tail removed 503 cells
(2.211%). `percent.ribo` is diagnostic only. This step applies no sample,
droplet/empty-drop, ambient-RNA, doublet, or high-complexity filter;
PIPseeker source-matrix selection and cell calling remain upstream.

Choose one preprocessing input for an analysis run. The default `legacy`
source selects the original
`INPUT_OBJECT_DIR/pipseq_processed_matrix_with_egfp.rds`; `counts-qc` selects
the counts-derived and QC-filtered `INPUT_OBJECT_DIR/sobj_qc_filtered.rds`:

```sh
just preprocess counts-qc
```

For one branch, use:

```sh
just preprocess-one pflog false counts-qc
```

An explicit alternative object remains available through
`Rscript scripts/preprocess-sobj.R --input /path/to/object.rds --normalization pflog`.
Preprocessing replaces the current branch artifacts, so clustering and all
downstream commands consume the source selected for the current run.

## Clustering Pipeline

Preview all current clustering commands without running clustering:

```sh
just cluster-dry-run
```

Run the all-branch clustering wrapper from the repo root:

```sh
just cluster
```

Generate the supplemental cluster grid table, 12-panel clustree grid, and
representative UMAP resolution sweep after clustered objects exist:

```sh
just summarize-clusters
```

The wrapper loads ESPI path constants in R; do not rely on a shell-exported
`CURRENT_OBJECT_DIR`. Clustered outputs use underscore branch tags such as
`pflog_no_filter_cc` because Seurat rewrites hyphens in reduction names.

## Marker Annotation Figures

Generate the per-cell cell type marker heatmap from the repo root:

```sh
just marker-heatmap
```

The script defaults to `cluster_pflog_filter_cc_dims50_res0.3` and the PFlog
expression layer, writes per-cell PNG/PDF outputs under `figures/annotation/`
in the Box data root, and symlinks the PNG into `notebook/figures/`.

Generate the per-cluster cell-type module and p27 enrichment heatmaps:

```sh
just cluster-marker-heatmaps 50 0.3
just summarize-mg-selected
just cluster-marker-heatmaps 30 0.3 /path/to/mg-selected-object.rds
```

The script writes PNG/PDF heatmaps under `figures/annotation/`, writes module
score and p27 enrichment TSVs under `tables/annotation/`, and symlinks PNGs
into `notebook/figures/`.


## Tripwire Checks

Run the lightweight scientific-boundary checks from the repo root:

```sh
just tripwires
```

The runner uses `analysis_labels.yml` for label/contrast declarations and exits
non-zero only for `FAIL` rows. `SKIP` rows mark checks that need future pipeline
stages or scratch-output instrumentation.


## Required R Packages

Declared in `DESCRIPTION`:

- `clustree`
- `ggplot2`
- `ggraph`
- `here`
- `igraph`
- `mclust`
- `patchwork`
- `scclrR` from `cleartools/scclrR`
- `Seurat`
- `SeuratObject`

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

Notebook figures should be symlinked into `notebook/figures/` and referenced with notebook-relative paths.

## Mycelium Local Hooks

Mycelium created local Claude Code hooks in `.claude/settings.local.json`. The file points at the installed plugin cache:

```text
/Users/carlstone/.omp/plugins/cache/plugins/mycelium___mycelium___0.0.0/
```

`.claude/` is gitignored because those paths are local and can break after plugin upgrades or cache cleanup. On a fresh clone or after upgrading the plugin, rerun Mycelium initialization or hook setup before expecting hooks to fire.
