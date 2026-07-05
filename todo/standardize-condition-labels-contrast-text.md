# Standardize condition labels and contrast text

| Field | Value |
|-------|-------|
| **Date** | 2026-07-05 |
| **Author** | OMP agent |
| **Priority** | medium |
| **Status** | complete |
| **Category** | refactor |
| **Related analyses** | `R/paths.R`; MG-selected figure and DE scripts |
| **Related data** | — |

## Description

Move project-wide condition labels out of the paths-focused constants file and standardize displayed contrasts as `p27CKO + E-Stim vs. p27CKO`.

## Motivation

Megan and Ed prefer explicit comparison labels. Programmatic contrast labels reduce drift across plots, captions, and tables.

## Proposed Approach

Create a project-values/constants file for condition labels and display labels. Update `CTRL_LABEL`, `ESTIM_LABEL`, or companion display constants so plot text can construct labels exactly as `(p27CKO + E-Stim vs. p27CKO)`. Replace hard-coded variants such as `(E-Stim vs control)` with constructed labels.

## Acceptance Criteria

- [x] Project-wide condition labels/display constants live in a clearly named R file instead of only `R/paths.R`.
- [x] Plot labels programmatically render `(p27CKO + E-Stim vs. p27CKO)`.
- [x] Hard-coded `(E-Stim vs control)` or similar variants are replaced where appropriate.
- [x] Changed R files are formatted, documented, linted, and affected figures/notebook are regenerated if labels change rendered outputs.

## Completion

Completed in Batch 1 presentation cleanup with metadata labels preserved and display contrast text standardized as `(p27CKO + E-Stim vs. p27CKO)`.

## Notes

Created from the review follow-up after commit `b3c248f`.
