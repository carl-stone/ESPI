
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

## Preprocessing

Run all preprocessing branches:

``` sh
Rscript scripts/preprocess-all.R
```

Or run one branch directly:

``` sh
Rscript scripts/preprocess-sobj.R --normalization log1p
Rscript scripts/preprocess-sobj.R --normalization log1p --filter-cell-cycle
Rscript scripts/preprocess-sobj.R --normalization pflog
Rscript scripts/preprocess-sobj.R --normalization pflog --filter-cell-cycle
```

Use `--input /path/to/object.rds` to override the default Trailmaker
input object.
