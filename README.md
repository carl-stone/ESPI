
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ESPI

<!-- badges: start -->

<!-- badges: end -->

This is the companion code to Stone<sub>\textit{et</sub>al.}~2026—.

## Installation

You can install the development version of ESPI from
[GitHub](https://github.com/) with:

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
Rscript scripts/process-counts.R
```

This saves `sobj_raw.rds` at `data/input/sobj_raw.rds` beneath the Box
data root (`DATA_ROOT_DIR/data/input/sobj_raw.rds`).

Run QC filtering before preprocessing:

``` sh
Rscript scripts/qc-filtering.R
```

The script writes QC figures to `FIGURE_DIR/qc/*.png`, QC tables to
`TABLE_DIR/qc/*.tsv`, and the filtered Seurat object to the Box path
`INPUT_OBJECT_DIR/sobj_qc_filtered.rds`. `percent.mt` uses all 37
observed mitochondrial features, whose labels are mixed (including
`mt-Rnr1`, `mt-Rnr2`, and non-`mt-` labels); none were lost upstream or
from the custom reference. It retains 22,248 of 983,903 cells across S2,
S3, S4, S5, S7, and S8 only when `nFeature_RNA >= 50`,
`nCount_RNA >= 100`, and `percent.mt <= 20`. This is a data-specific
cutoff, not a universal PipSeq rule: among the 22,751 cells meeting the
complexity thresholds, complete-mitochondrial P95/P97.5/P99 were
16.038/19.313/27.666%, and the \>20% sparse extreme tail removed 503
cells (2.211%). `percent.ribo` is diagnostic only. This step applies no
sample, droplet/empty-drop, ambient-RNA, doublet, or high-complexity
filter, and does not replace upstream PIPseeker source-matrix selection
or cell calling.

Preprocess the QC-filtered object with an explicit input path:

``` sh
Rscript scripts/preprocess-sobj.R \
  --input INPUT_OBJECT_DIR/sobj_qc_filtered.rds \
  --normalization pflog
```

The equivalent `just` command is
`just preprocess-one pflog false INPUT_OBJECT_DIR/sobj_qc_filtered.rds`.
Repeat with the desired normalization and cell-cycle-HVG settings.

Preview clustering commands before running them, then cluster all
preprocessed objects and summarize the clustering grid:

``` sh
just cluster-dry-run
just cluster
just summarize-clusters
```

Use `just --list` to discover figure, marker, and MG-selected follow-up
recipes.

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
