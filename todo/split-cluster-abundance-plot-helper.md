# Split cluster abundance plot helper into plot file

| Field | Value |
|-------|-------|
| **Date** | 2026-07-05 |
| **Author** | OMP agent |
| **Priority** | low |
| **Status** | complete |
| **Category** | refactor |
| **Related analyses** | `R/cluster-abundance.R`; `scripts/plot-mg-selected-figures.R` |
| **Related data** | — |

## Description

Move `plot_clr_fisher_enrichment()` into a new `R/cluster-abundance-plots.R` file while keeping computation code in `R/cluster-abundance.R`.

## Motivation

Separating plot helpers by analysis step makes it easier to find figure-producing code associated with each workflow stage.

## Proposed Approach

Create `R/cluster-abundance-plots.R`, move the plot helper and plot-only support code there, leave `compute_cluster_abundance()` in the compute file, and rerun roxygen documentation.

## Acceptance Criteria

- [x] `plot_clr_fisher_enrichment()` lives in `R/cluster-abundance-plots.R`.
- [x] `compute_cluster_abundance()` remains in `R/cluster-abundance.R`.
- [x] `NAMESPACE` and Rd files are regenerated with `devtools::document()`.
- [x] Changed R files are formatted and linted.

## Completion

Completed in Batch 1 presentation cleanup by separating cluster-abundance plotting helpers from compute/test code.

## Notes

Created from the review follow-up after commit `b3c248f`.
