# Filter DE/DD effect size labels to significant genes

| Field | Value |
|-------|-------|
| **Date** | 2026-07-05 |
| **Author** | OMP agent |
| **Priority** | medium |
| **Status** | complete |
| **Category** | analysis |
| **Related analyses** | — |
| **Related data** | — |

## Description

For the DE and DD effect size comparison plot, label genes with `geom_text_repel()` only when they are significant in one or both DE/DD tests. Do not label nonsignificant genes.

## Motivation

Filtering labels to significant genes keeps the plot readable and focuses attention on results supported by the tests.

## Proposed Approach

Find the DE/DD effect size comparison plot source and restrict the label data passed to `ggrepel::geom_text_repel()` to genes that meet the existing significance criteria in either the DE test, the DD test, or both.

## Acceptance Criteria

- [x] Significant DE-only, DD-only, and shared significant genes can be labeled.
- [x] Nonsignificant genes are not labeled.
- [x] Label placement uses `ggrepel::geom_text_repel()`.
- [x] The plot output is regenerated from the updated labeling rule.

## Notes

Requested as a future plotting cleanup item.
