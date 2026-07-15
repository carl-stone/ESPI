
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

## Frozen analysis workflow

The count processing, QC, source preprocessing and clustering, source
summaries and marker heatmaps, MG selection, and MG clustering are
frozen. Routine commands consume the existing clustered MG-selected RDS
objects; they do not run scripts `01` through `07` or `04-cluster.R`.

Inspect the downstream plan, then run it:

``` sh
just run-dry-run
just run
```

These commands start at `scripts/08-summarize-mg-clusters.R`, regenerate
the downstream MG figures, marker tables, DE and enrichment outputs,
render `notebook/sc_analysis.qmd`, and run `tools/run-tripwires.R`. The
optional `overwrite` argument is `false` or `true`; existing protected
marker and DE outputs stop the run unless replacement is explicitly
allowed:

``` sh
just run-dry-run true
just run true
```

Only regenerate frozen data deliberately. Inspect the complete
regeneration plan before executing it:

``` sh
just regenerate-frozen-dry-run counts-qc false
just regenerate-frozen counts-qc false
```

For explicit regeneration, `source` may be `counts-qc`, `legacy`, or a
quoted Seurat RDS path. The regeneration recipe runs count processing
when applicable, source preprocessing and clustering, source summaries
and marker heatmaps, MG selection and clustering, then the downstream
analysis.

Current selected identifiers and counts:

| Role | Selected branch and cluster column | Result |
|----|----|----|
| Source | `pflog_no_filter_cc`; `cluster_pflog_no_filter_cc_dims30_res0.3` | 4,146 cells; 9 clusters; exclude source clusters 2, 7, and 8; retain 3,456 |
| MG-selected | `pflog_mg_selected_no_filter_cc`; `cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5` | MG PCA/candidate depth 50; 8 chosen clusters |

Use the low-level entry points below only for intentional checkpoint
recovery.

| Checkpoint | Existing recipe or script |
|----|----|
| Combined raw object | `Rscript scripts/01-process-counts.R` |
| QC annotation and filtered object | `Rscript scripts/02-qc-filtering.R` |
| Preprocessing branches | `just preprocess counts-qc` (`scripts/03-preprocess-all.R`) or `just preprocess-one` |
| Source clustering and summary | `just cluster` (`scripts/04-cluster-all.R`), then `just summarize-clusters` |
| MG subset selection | `Rscript scripts/07-select-mg-subset.R` with source column `cluster_pflog_no_filter_cc_dims30_res0.3` and `--dims 50` |
| MG clustering and summary | `just cluster-one` / `just summarize-mg-selected` |
| Figures and marker tables | `just marker-heatmap`, `just cluster-marker-heatmaps`, `just mg-figures`, `just mg-markers-no-merge` |
| MG DE and enrichment | `just mg-de` (`scripts/12-run-mg-de.R`); defaults to apeglm shrinkage and writes raw plus simplified GO BP ORA/GSEA tables, seeded Bayesian ORA comparisons, and notebook §4 dotplots |
| Notebook and final checks | `just notebook`, then `just tripwires` |

Use `just --list` to inspect expert and maintenance recipes. Load
package code after setup or dependency changes with `just load`.
