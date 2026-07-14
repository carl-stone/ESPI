
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ESPI

<!-- badges: start -->

<!-- badges: end -->

This is the companion code to Stone<sub>\textit{et</sub>al.}~2026—.

## Installation

You can install the development version of ESPI from the [project
repository](https://github.com/carl-stone/ESPI) with:

``` r
# install.packages("pak")
pak::pak("carl-stone/ESPI")
```

## End-to-end workflow

Use `just` from the repo root for package tasks, pipeline steps,
notebook rendering, and tripwire checks:

``` sh
just --list
```

Load the package after setup or dependency changes:

``` sh
just load
```

Create the combined raw Seurat object from the six 10X directories and
sample metadata:

``` sh
Rscript scripts/01-process-counts.R
```

This saves `sobj_raw.rds` at `data/input/sobj_raw.rds` beneath the Box
data root (`DATA_ROOT_DIR/data/input/sobj_raw.rds`).

Run QC filtering before preprocessing:

``` sh
Rscript scripts/02-qc-filtering.R
```

The script writes QC figures to `FIGURE_DIR/qc/*.png`, QC tables to
`TABLE_DIR/qc/*.tsv`, the annotated raw object to
`INPUT_OBJECT_DIR/sobj_raw_with_qc.rds`, and the filtered object to
`INPUT_OBJECT_DIR/sobj_qc_filtered.rds`. `percent.mt` uses all 37
observed mitochondrial features, whose labels are mixed (including
`mt-Rnr1`, `mt-Rnr2`, and non-`mt-` labels).
`DropletUtils::emptyDrops()` supplies the cell-call FDR and `is_cell`
flag. Called cells define sample-specific lower three-MAD thresholds for
log10 counts and detected features and an upper three-MAD threshold for
mitochondrial percentage. The saved filtered object contains the 4,145
cells that pass those three MAD criteria; `is_cell` remains separate
from `pass_qc`. `percent.ribo` remains diagnostic only.

Choose one preprocessing input for an analysis run. The default `legacy`
source is the original
`INPUT_OBJECT_DIR/pipseq_processed_matrix_with_egfp.rds`. The
`counts-qc` source is the counts-derived and QC-filtered
`INPUT_OBJECT_DIR/sobj_qc_filtered.rds`:

``` sh
just preprocess counts-qc
```

Run a single branch with the same source choice:

``` sh
just preprocess-one pflog false counts-qc
```

Use `--input` only for another explicit Seurat object. Preprocessing
replaces the current branch artifacts, so cluster and downstream
commands consume the source selected for the current run.

Preview clustering commands before running them, then cluster all
preprocessed objects and summarize the clustering grid:

``` sh
just cluster-dry-run
just cluster
just summarize-clusters
```

The current chosen counts-derived analysis uses the no-cell-cycle-HVG
PFlog branch at 20 PCs and Leiden resolution 0.3. The downstream
commands below reproduce the source and MG-selected figures, marker
ranking, and Mouse × Condition pseudobulk analysis. MG preprocessing
computes 50 PCs so the 30- and 50-PC sensitivity candidates remain
available.

``` sh
DATA_ROOT="/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data"
OBJECT_DIR="$DATA_ROOT/seurat_objects/current"

FULL_OBJECT="$OBJECT_DIR/cluster_pflog_no_filter_cc_elbow20.rds"
MG_PREPROCESS="$OBJECT_DIR/preprocess_pflog_mg_selected_no_filter_cc.rds"
MG_FILTER_PREPROCESS="$OBJECT_DIR/preprocess_pflog_mg_selected_filter_cc.rds"
MG_OBJECT="$OBJECT_DIR/cluster_pflog_mg_selected_no_filter_cc_elbow20.rds"
MG_FILTER_OBJECT="$OBJECT_DIR/cluster_pflog_mg_selected_filter_cc_elbow20.rds"

Rscript scripts/07-select-mg-subset.R \
  --input "$FULL_OBJECT" \
  --cluster-column cluster_pflog_no_filter_cc_dims20_res0.3 \
  --dims 50

just cluster-one "$MG_PREPROCESS" 20
just cluster-one "$MG_FILTER_PREPROCESS" 20
just summarize-mg-selected 20

just marker-heatmap 20 0.3 "$FULL_OBJECT"
just marker-heatmap 20 0.3 "$MG_OBJECT"
just marker-heatmap 20 0.3 "$MG_FILTER_OBJECT"
just cluster-marker-heatmaps 20 0.3 "$FULL_OBJECT"
just cluster-marker-heatmaps 20 0.3 "$MG_OBJECT"
just cluster-marker-heatmaps 20 0.3 "$MG_FILTER_OBJECT"

just mg-figures \
  "$MG_OBJECT" pflog_mg_selected_no_filter_cc 20 20 0.3
just mg-figures \
  "$MG_FILTER_OBJECT" pflog_mg_selected_filter_cc 20 20 0.3

Rscript scripts/11-find-mg-markers.R \
  --input "$MG_OBJECT" \
  --branch-tag pflog_mg_selected_no_filter_cc \
  --elbow-n 20 \
  --dims 20 \
  --resolution 0.3 \
  --layer data \
  --counts-layer counts \
  --confirm-no-merge

just mg-de \
  "$MG_OBJECT" \
  cluster_pflog_mg_selected_no_filter_cc_dims20_res0.3
```

The shorter MG recipes still default to historical 30- or 50-PC
candidate columns. Pass the explicit 20-PC identifiers above for the
current chosen analysis. Use `just --list` to inspect every
parameterized figure, marker, and MG-selected recipe.

Rebuild generated documentation and render the embedded-resource
notebook after source or figure changes:

``` sh
just document
just readme
just notebook
```

Run tripwire checks before handing off analysis changes:

``` sh
just tripwires
```

The underlying scripts still accept raw flags when a recipe does not
fit.
