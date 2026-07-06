# Bundle DEG and enrichment CSVs for Ed

| Field | Value |
|-------|-------|
| **Date** | 2026-07-06 |
| **Author** | OMP agent |
| **Priority** | medium |
| **Status** | open |
| **Category** | analysis |
| **Related analyses** | MG-selected DE/DD and enrichment outputs |
| **Related data** | `DEG_DIR/mg_selected/`; `ENRICHMENT_DIR/mg_selected/` |

## Description

Create a small sendable bundle of differential-expression and enrichment CSVs for Ed.

## Motivation

Ed needs the relevant DEG and GSEA/GO result tables in a portable format without needing to navigate the full analysis output tree.

## Proposed Approach

- Select the current MG-selected DEG tables that are useful for external review.
- Select the current GO/GSEA enrichment tables that summarize the same comparison.
- Export or copy them as CSVs into a small clearly named bundle directory or archive.
- Include only lightweight tabular outputs, not large objects or figures unless explicitly requested.

## Acceptance Criteria

- [ ] Bundle includes the relevant DEG CSVs.
- [ ] Bundle includes the relevant GSEA/GO CSVs.
- [ ] File names are self-explanatory for someone outside the repo.
- [ ] Bundle location is easy to send to Ed.

## Notes

The current source outputs are TSVs under the MG-selected DEG and enrichment output directories. Convert to CSVs for the bundle unless Ed prefers TSVs.
