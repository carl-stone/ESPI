# Standardize stable references for cross-referenced artifacts

| Field | Value |
|-------|-------|
| **Date** | 2026-07-05 |
| **Author** | OMP agent |
| **Priority** | medium |
| **Status** | complete |
| **Category** | writing |
| **Related analyses** | `notebook/sc_analysis.qmd`; `.living/log/` |
| **Related data** | — |

## Description

Use stable identifiers for figures, tables, decisions, and other cross-referenced artifacts instead of fragile auto-numbered labels where possible.

## Motivation

Auto-numbered figure labels can shift when manuscript panels are inserted, removed, or restored. Stable references make communication and audit logs clearer.

## Proposed Approach

Prefer Quarto figure IDs, explicit captions, file basenames, decision headings, or other stable handles in notebook prose, session logs, review reports, and Mycelium records.

## Acceptance Criteria

- [x] Document the project convention for stable cross-references.
- [x] Supersede the original criterion, `Replace current fragile Figure N references in review/session records where stable IDs are available`, via the `Require stable cross-references going forward` decision and historical-fidelity exception; do not rewrite existing historical review/session evidence solely to replace `Figure N` text.
- [x] Future notebook/log prose uses stable IDs or captions for cross-referenced artifacts.

## Completion

Completed in Batch 1 presentation cleanup for forward-going records. Added the Cross-References convention and recorded the historical-record exception; the original criterion, `Replace current fragile Figure N references in review/session records where stable IDs are available`, was superseded rather than performed, so review/log evidence remains intact by design.

## Notes

Created from the review follow-up after commit `b3c248f`.
