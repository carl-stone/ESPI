# Order CLR log2 enrichment columns by contrast enrichment

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

For the CLR log2 enrichment figure, order the columns by enrichment for the `p27CKO + EStim` versus `p27CKO` contrast.

## Motivation

Ordering columns by the contrast effect makes the strongest enrichment visually prominent and makes the left-to-right pattern easier to compare.

## Proposed Approach

Find the CLR log2 enrichment plot source and sort the plotted columns by the existing enrichment value for `p27CKO + EStim` versus `p27CKO`, descending from highest enrichment on the left to lowest enrichment on the right.

## Acceptance Criteria

- [x] The highest `p27CKO + EStim` versus `p27CKO` enrichment appears in the leftmost column.
- [x] The lowest enrichment appears in the rightmost column.
- [x] The figure output is regenerated from the updated ordering.

## Notes

Requested as a future plotting cleanup item.
