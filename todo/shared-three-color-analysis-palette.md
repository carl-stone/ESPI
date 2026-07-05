# Create shared three-color analysis palette

| Field | Value |
|-------|-------|
| **Date** | 2026-07-05 |
| **Author** | OMP agent |
| **Priority** | medium |
| **Status** | open |
| **Category** | infrastructure |
| **Related analyses** | `R/themes.R`; MG-selected figure scripts |
| **Related data** | — |

## Description

Create a project-wide palette for three-color low/medium/high or depleted/not-significant/enriched encodings using `#2166ac`, `grey75`, and `#e31a8c`.

## Motivation

ESPI already exports `palette_dotplot_pair`. Standardizing the three-color palette prevents duplicated hard-coded color values and keeps manuscript figures consistent.

## Proposed Approach

Add an exported palette in `R/themes.R`, e.g. blue/gray/pink values for low/mid/high. Use the gray entry anywhere plots currently hard-code a not-significant gray. For low-to-high scales, use the gray-to-pink subset where appropriate, such as feature UMAP-style plots.

## Acceptance Criteria

- [ ] `R/themes.R` exports a named three-color analysis palette.
- [ ] Existing hard-coded not-significant gray colors are replaced with the palette gray.
- [ ] Existing hard-coded `#2166ac` / `#e31a8c` direction colors use the shared palette where appropriate.
- [ ] Changed R files are formatted, documented, and linted.

## Notes

Created from the review follow-up after commit `b3c248f`.
