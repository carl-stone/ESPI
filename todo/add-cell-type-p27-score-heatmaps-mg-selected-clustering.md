# Add cell-type and p27 score heatmaps to MG-selected clustering

| Field | Value |
|-------|-------|
| **Date** | 2026-07-05 |
| **Author** | OMP agent |
| **Priority** | medium |
| **Status** | complete |
| **Category** | analysis |
| **Related analyses** | `notebook/sc_analysis.qmd` section 2.2, `### MG-selected clustering` |
| **Related data** | MG-selected clustering outputs |

## Description

Add compact heatmaps to the MG-selected clustering section that summarize cell-type module scores and p27 enrichment across clusters.

## Motivation

The MG-selected clustering section currently describes how clusters were selected and shows the MG-selected UMAP. Small heatmaps would make the cluster annotation signals easier to compare in place, including both cell-type marker programs and p27-related enrichment.

## Proposed Approach

- Add a small full-dataset heatmap in `notebook/sc_analysis.qmd` section 2.2, `### MG-selected clustering`, which starts at line 128. Insert it immediately after the first paragraph (currently lines 130-140) and before the `I reselected HVGs...` paragraph.
- Add the matching MG-selected heatmap immediately after the MG-selected UMAP figure `#fig-mg-selected-cluster-umap`, currently around line 150.
- Each heatmap should show module scores for `cell_type_marker_genes` and a sample-aware cluster-permuted p27 enrichment score.
- Keep the heatmaps compact and focused on cluster-level interpretation; do not replace the existing UMAP or abundance figures.

## Acceptance Criteria

- [x] A full-dataset heatmap appears immediately after the first paragraph of section 2.2, `### MG-selected clustering` (currently lines 130-140), and before the `I reselected HVGs...` paragraph.
- [x] A matching MG-selected heatmap appears immediately after the MG-selected UMAP figure `#fig-mg-selected-cluster-umap`, currently around line 150.
- [x] Both heatmaps include module scores for `cell_type_marker_genes`.
- [x] Both heatmaps include a sample-aware cluster-permuted p27 enrichment score.
- [x] Notebook prose briefly explains what the heatmaps show without treating them as a replacement for the existing clustering or abundance analyses.

## Notes

Requested as a future notebook/figure task for the MG-selected clustering section. The requested placement is immediately after the first paragraph of `### MG-selected clustering` (section starts at `notebook/sc_analysis.qmd` line 128; paragraph currently lines 130-140) for the full-dataset heatmap and immediately after Figure 6 / `#fig-mg-selected-cluster-umap` (currently around line 150) for the MG-selected heatmap.
