# Environments & Installations

## Primary Environment

- **Project type**: Minimal R package for ESPI single-cell analysis.
- **Package load**: `devtools::load_all()` from the repo root.
- **Formatter**: Air (`air format <file>`).
- **Notebook renderer**: Quarto (`quarto render notebook/sc_analysis.qmd`).
- **External data root**: `/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data`.

## Setup From Scratch

```r
# From an R session with working directory at the repo root
devtools::load_all()
devtools::document()
```

```sh
# Format changed R files
air format R/<file>.R scripts/<script>.R

# Render the notebook after figure source changes
quarto render notebook/sc_analysis.qmd
```

## Clustering Pipeline

Preview all current clustering commands without running clustering:

```sh
Rscript scripts/cluster-all.R --dry-run --elbow-n 20
```

Run the all-branch clustering wrapper from the repo root:

```sh
Rscript scripts/cluster-all.R --elbow-n 20
```

Generate the supplemental cluster grid table, 12-panel clustree grid, and
representative UMAP resolution sweep after clustered objects exist:

```sh
Rscript scripts/summarize-cluster-grid.R
```

The wrapper loads ESPI path constants in R; do not rely on a shell-exported
`CURRENT_OBJECT_DIR`. Clustered outputs use underscore branch tags such as
`pflog_no_filter_cc` because Seurat rewrites hyphens in reduction names.

## Marker Annotation Figures

Generate the per-cell cell type marker heatmap from the repo root:

```sh
Rscript scripts/big-heatmap-plot.R
```

The script defaults to `cluster_pflog_filter_cc_dims50_res0.3` and the PFlog
expression layer, writes per-cell PNG/PDF outputs under `figures/annotation/`
in the Box data root, and symlinks the PNG into `notebook/figures/`.


## Tripwire Checks

Run the lightweight scientific-boundary checks from the repo root:

```sh
Rscript tools/run-tripwires.R
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
