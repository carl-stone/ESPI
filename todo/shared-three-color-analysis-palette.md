# Create shared three-color analysis palette

| Field | Value |
|-------|-------|
| **Date** | 2026-07-05 |
| **Author** | OMP agent |
| **Priority** | medium |
| **Status** | complete |
| **Category** | infrastructure |
| **Related analyses** | `R/config.R`; `scripts/02-publication-figures.R`; `scripts/03-marker-analysis.R`; `scripts/04-de-enrichment.R` |
| **Related data** | — |

## Description

Create a project-wide palette for three-color low/medium/high or depleted/not-significant/enriched encodings using `#2166ac`, `grey75`, and `#e31a8c`.

## Motivation

ESPI already exports `palette_dotplot_pair`. Standardizing the three-color palette prevents duplicated hard-coded color values and keeps manuscript figures consistent.

## Implemented Approach

Define the named blue/gray/pink palette in `R/config.R` and expose it through `publication_config()` for the publication phases. Use the gray entry wherever plots encode a not-significant midpoint and use the appropriate palette subset for low-to-high scales.

## Acceptance Criteria

- [x] `R/config.R` defines the named three-color analysis palette.
- [x] Existing hard-coded not-significant gray colors are replaced with the palette gray.
- [x] Existing hard-coded `#2166ac` / `#e31a8c` direction colors use the shared palette where appropriate.
- [x] Changed R files are formatted, documented, and linted.

## Completion

Completed in Batch 1 presentation cleanup with a shared low/mid/high analysis palette used by the affected plotting code.

## Notes

Created from the review follow-up after commit `b3c248f`.
