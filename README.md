
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

Run all preprocessing branches, or one branch with an explicit
normalization and cell-cycle-HVG filter setting:

``` sh
just preprocess
just preprocess-one log1p false
just preprocess-one log1p true
just preprocess-one pflog false
just preprocess-one pflog true
```

Pass an explicit input path to override the default Trailmaker input
object: `just preprocess-one pflog false /path/to/object.rds`.

Preview clustering commands before running them, then cluster all
preprocessed objects and summarize the clustering grid:

``` sh
just cluster-dry-run
just cluster
just summarize-clusters
```

Generate marker-annotation outputs and MG-selected follow-up outputs:

``` sh
just marker-heatmap
just cluster-marker-heatmaps
just summarize-mg-selected
just mg-markers-no-merge
just mg-figures
just mg-de-overwrite
```

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
