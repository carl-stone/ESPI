# Bundle DEG and enrichment CSVs for Ed

| Field | Value |
|-------|-------|
| **Date** | 2026-07-06 |
| **Author** | OMP agent |
| **Priority** | medium |
| **Status** | complete |
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

- [x] Bundle includes the relevant DEG CSVs.
- [x] Bundle includes the relevant GSEA/GO CSVs.
- [x] File names are self-explanatory for someone outside the repo.
- [x] Bundle location is easy to send to Ed.

## Notes

The current source outputs are TSVs under the MG-selected DEG and enrichment output directories. Convert to CSVs for the bundle unless Ed prefers TSVs.

## Completion

- Bundle directory: `/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data/exports/ed_mg_selected_de_dd_gsea_go_2026-07-06/`
- Sendable archive: `/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data/exports/ed_mg_selected_de_dd_gsea_go_2026-07-06.zip`
- Contents: 18 result CSVs plus `manifest.csv`.
- Note: `mg_selected_detection_primary_significant_genes.csv` has zero rows because no primary differential-detection genes pass `padj < 0.05` in the current output.
- Note: GO/GSEA includes both FDR-filtered term CSVs and full all-term conversion CSVs; filenames and `manifest.csv` label the difference.
