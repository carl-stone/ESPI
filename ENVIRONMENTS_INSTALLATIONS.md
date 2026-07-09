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
`brew install just`). Then smoke-test the package load from the repo root:

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

Run all preprocessing branches:

```sh
just preprocess
```

Run one branch directly, with defaults matching the R script:

```sh
just preprocess-one pflog false
```

Pass an explicit input path to override the default Trailmaker input object:
`just preprocess-one pflog false /path/to/object.rds`.

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
just cluster-marker-heatmaps 30 0.3
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
