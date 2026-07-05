# Standardize stable references for cross-referenced artifacts

| Field | Value |
|-------|-------|
| **Date** | 2026-07-05 |
| **Author** | OMP agent |
| **Priority** | medium |
| **Status** | open |
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

- [ ] Document the project convention for stable cross-references.
- [ ] Replace current fragile `Figure N` references in review/session records where stable IDs are available.
- [ ] Future notebook/log prose uses stable IDs or captions for cross-referenced artifacts.

## Notes

Created from the review follow-up after commit `b3c248f`.
