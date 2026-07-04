---
topic: population-structure
description: Empirical findings about marker-defined cellular structure within focused analysis branches. Entries treat cluster marker rankings as descriptive cell-level evidence unless paired with replicate-level tests.
created: 2026-07-04
last_updated: 2026-07-04
status: active
---

# Population Structure

## F-003: Cluster 2 lacks a positive marker-defined identity
**Status:** preliminary
**Claim:** In the MG-selected chosen clustering, cluster 2 contains 961 cells but has zero retained positive marker genes under the `FindAllMarkers` positive detection-enrichment filter (`pct.1 > pct.2`).
**Implications:** Cluster 2 should not yet be treated as a defensible interpreted identity; the result argues for caution or later merge review, but does not identify a specific merge partner.
**Tags:** MG-selected, clustering, markers, interpreted identity, negative result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-04 | notebook/sc_analysis.qmd cluster marker section; tables/mg_selected/find_all_markers_summary_data_pflog_mg_selected_no_filter_cc_dims30_res0.3.csv | MG-selected chosen clustering, PFlog branch, 30 PCs, resolution 0.3 | ESPI | Marker summary reports cluster 2 with 961 cells, `n_retained_markers = 0`, and `n_top_markers = 0`; notebook notes this argues against treating cluster 2 as a marker-defined identity. | supports |

### Open Questions
- Does cluster 2 merge with a neighboring marker profile under a curated merge map, or does it reflect a non-marker technical/state axis?
- Are there negative markers or continuous gradients that explain cluster 2 despite the absence of positive markers?

## F-004: Cluster 3 carries a neurogenic progenitor marker program
**Status:** preliminary
**Claim:** In the MG-selected chosen clustering, cluster 3 is enriched for neurogenic-progenitor marker genes, including `Ascl1`, `Hes6`, and `Dll1`, in descriptive cell-level marker rankings.
**Implications:** Cluster 3 is the strongest marker-ranked neurogenic-progenitor-like subpopulation in the selected branch, but this does not by itself show a Mouse × Condition condition effect.
**Tags:** MG-selected, clustering, markers, neurogenic progenitor, cell-level ranking

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-04 | notebook/sc_analysis.qmd cluster marker section; tables/mg_selected/find_all_markers_summary_data_pflog_mg_selected_no_filter_cc_dims30_res0.3.csv; tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims30_res0.3.csv | MG-selected chosen clustering, PFlog branch, 30 PCs, resolution 0.3 | ESPI | Cluster 3 has 446 cells and 2,614 retained positive markers; top-five marker rows include `Ascl1` (avg_log2FC 2.80, pct.1 0.881 vs pct.2 0.137), `Hes6` (2.57, 0.975 vs 0.436), and `Dll1` (3.11, 0.751 vs 0.162). | supports |

### Open Questions
- Does the cluster 3 marker program correspond to a stable interpreted identity across normalization branches or cluster resolutions?
- How much of the cluster 3 program overlaps the primary DE `Hes6` signal versus cell-composition or within-cell abundance changes?

## F-005: Cluster 8 is a very small outlier with strong marker separation
**Status:** preliminary
**Claim:** In the MG-selected chosen clustering, cluster 8 is a 15-cell outlier with strong descriptive marker separation, including `Cnmd`, `Cldn10`, `Scrg1`, `Htr1b`, and `Flt4`.
**Implications:** Cluster 8 may represent a rare state or outlier population, but its small size makes it unsuitable for broad identity or condition-effect claims without additional validation.
**Tags:** MG-selected, clustering, markers, outlier, rare population

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-04 | notebook/sc_analysis.qmd cluster marker section; tables/mg_selected/find_all_markers_summary_data_pflog_mg_selected_no_filter_cc_dims30_res0.3.csv; tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims30_res0.3.csv | MG-selected chosen clustering, PFlog branch, 30 PCs, resolution 0.3 | ESPI | Cluster 8 has 15 cells and 6,825 retained positive markers; top-five rows show high marker separation for `Htr1b` (avg_log2FC 9.92, pct.1 0.467 vs pct.2 0), `Flt4` (9.52, 0.600 vs 0.001), `Cnmd` (8.99, 1.000 vs 0.006), `Cldn10` (8.91, 0.867 vs 0.004), and `Scrg1` (8.67, 0.867 vs 0.003). | supports |

### Open Questions
- Is cluster 8 reproducible across preprocessing/clustering choices, or is it driven by a few high-signal outlier cells?
- Do cluster 8 cells have QC, doublet, or sample-origin patterns that explain the strong marker separation?
