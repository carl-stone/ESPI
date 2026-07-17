# Standardize condition labels and contrast text

| Field | Value |
|-------|-------|
| **Date** | 2026-07-05 |
| **Author** | OMP agent |
| **Priority** | medium |
| **Status** | complete |
| **Category** | refactor |
| **Related analyses** | `R/config.R`; `scripts/02-publication-figures.R`; `scripts/04-de-enrichment.R` |
| **Related data** | — |

## Description

Centralize project-wide metadata and display condition labels in `R/config.R` and standardize displayed contrasts as `p27CKO + E-Stim vs. p27CKO`.

## Motivation

Megan and Ed prefer explicit comparison labels. Programmatic contrast labels reduce drift across plots, captions, and tables.

## Implemented Approach

`R/config.R` defines the metadata labels, display labels, and constructed contrast text consumed by publication figures and DE/enrichment. Hard-coded variants such as `(E-Stim vs control)` were replaced where appropriate.

## Acceptance Criteria

- [x] Project-wide condition labels and display constants live in `R/config.R`.
- [x] Plot labels programmatically render `(p27CKO + E-Stim vs. p27CKO)`.
- [x] Hard-coded `(E-Stim vs control)` or similar variants are replaced where appropriate.
- [x] Changed R files are formatted, documented, linted, and affected figures/notebook are regenerated if labels change rendered outputs.

## Completion

Completed in Batch 1 presentation cleanup with metadata labels preserved and display contrast text standardized as `(p27CKO + E-Stim vs. p27CKO)`.

## Notes

Created from the review follow-up after commit `b3c248f`.
