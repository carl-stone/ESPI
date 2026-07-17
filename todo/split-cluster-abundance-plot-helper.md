# Split cluster abundance plot helper into plot file

| Field | Value |
|-------|-------|
| **Date** | 2026-07-05 |
| **Author** | OMP agent |
| **Priority** | low |
| **Status** | complete |
| **Category** | refactor |
| **Related analyses** | `R/publication-analysis.R`; `scripts/02-publication-figures.R` |
| **Related data** | — |

## Description

This completed task originally separated cluster-abundance computation from plotting. The later publication-pipeline consolidation superseded that file split: computation now lives in `R/publication-analysis.R`, while the one-use abundance plot remains visible in `scripts/02-publication-figures.R`.

## Motivation

Separating plot helpers by analysis step makes it easier to find figure-producing code associated with each workflow stage.

## Implemented Approach

The original split was completed, then simplified during consolidation. The
current design keeps `compute_cluster_abundance()` in
`R/publication-analysis.R` and the single-use ggplot construction in
`scripts/02-publication-figures.R`; no dedicated cluster-abundance plot module
remains.

## Acceptance Criteria

- [x] Cluster-abundance computation is separate from plot construction.
- [x] `compute_cluster_abundance()` lives in `R/publication-analysis.R`.
- [x] The one-use abundance plot stays visible in phase 02.
- [x] Changed R files were documented, formatted, and linted.

## Completion

Completed in Batch 1 presentation cleanup. The later four-phase consolidation preserved the computation/plot separation while folding the code into the current files listed above.

## Notes

Created from the review follow-up after commit `b3c248f`.
