
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

## Command interface

Use `just` from the repo root for common package, pipeline, notebook,
and tripwire commands:

``` sh
just --list
```

## Preprocessing

Run all preprocessing branches:

``` sh
just preprocess
```

Or run one branch directly:

``` sh
just preprocess-one log1p false
just preprocess-one log1p true
just preprocess-one pflog false
just preprocess-one pflog true
```

Pass an explicit input path to override the default Trailmaker input
object: `just preprocess-one pflog false /path/to/object.rds`. The
underlying script still accepts its raw flags when needed.
