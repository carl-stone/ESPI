---
topic: population-structure
description: Empirical findings about marker-defined cellular structure within focused analysis branches. Entries treat cluster marker rankings as descriptive cell-level evidence unless paired with replicate-level tests.
created: 2026-07-04
last_updated: 2026-07-13
status: active
---

# Population Structure

## F-003: Cluster 2 lacks a positive marker-defined identity
**Status:** superseded
**Claim:** The prior 30-PC MG-selected clustering contained a cluster 2 without retained positive markers; this does not hold in the current 20-PC emptyDrops/log-MAD rebuild.
**Implications:** Treat the earlier marker-negative cluster 2 as a run-specific historical result, not a current interpreted population.
**Tags:** MG-selected, clustering, markers, interpreted identity, negative result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-04 | notebook/sc_analysis.qmd cluster marker section; tables/mg_selected/find_all_markers_summary_data_pflog_mg_selected_no_filter_cc_dims30_res0.3.csv | MG-selected chosen clustering, PFlog branch, 30 PCs, resolution 0.3 | ESPI | Marker summary reports cluster 2 with 961 cells, `n_retained_markers = 0`, and `n_top_markers = 0`; notebook notes this argues against treating cluster 2 as a marker-defined identity. | supports |
| 2026-07-13 | session 2026-07-13-003; `tables/mg_selected/find_all_markers_summary_data_pflog_mg_selected_no_filter_cc_dims20_res0.3.csv` | MG-selected emptyDrops/log-MAD rebuild, 20 PCs, resolution 0.3 | ESPI | Current cluster 2 contains 937 cells and 894 retained markers, including `Aldh1a1`, `Nxph1`, `Nr2f2`, `Efnb2`, and `Tes`; this invalidates the prior marker-negative cluster-number claim for the current run. | supersedes |

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
| 2026-07-13 | session 2026-07-13-003; `tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims20_res0.3.csv` | MG-selected emptyDrops/log-MAD rebuild, 20 PCs, resolution 0.3 | ESPI | Current cluster 3 contains 453 cells and retains the neurogenic program: `Ascl1` (avg_log2FC 3.32), `Hes6` (2.67), and `Dll1` (2.98) are enriched. | supports |

### Open Questions
- Does the cluster 3 marker program correspond to a stable interpreted identity across normalization branches or cluster resolutions?
- How much of the cluster 3 program overlaps the primary DE `Hes6` signal versus cell-composition or within-cell abundance changes?

## F-005: Cluster 7 is a very small outlier with strong marker separation
**Status:** preliminary
**Claim:** In the current MG-selected chosen clustering, cluster 7 is a 17-cell outlier with strong descriptive marker separation, including `Cnmd`, `Scrg1`, `Gja1`, `Ccdc190`, and `Cldn19`.
**Implications:** Cluster 7 may represent a rare state or outlier population, but its small size makes it unsuitable for broad identity or condition-effect claims without additional validation.
**Tags:** MG-selected, clustering, markers, outlier, rare population

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-04 | notebook/sc_analysis.qmd cluster marker section; tables/mg_selected/find_all_markers_summary_data_pflog_mg_selected_no_filter_cc_dims30_res0.3.csv; tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims30_res0.3.csv | MG-selected chosen clustering, PFlog branch, 30 PCs, resolution 0.3 | ESPI | Cluster 8 has 15 cells and 6,825 retained positive markers; top-five rows show high marker separation for `Htr1b` (avg_log2FC 9.92, pct.1 0.467 vs pct.2 0), `Flt4` (9.52, 0.600 vs 0.001), `Cnmd` (8.99, 1.000 vs 0.006), `Cldn10` (8.91, 0.867 vs 0.004), and `Scrg1` (8.67, 0.867 vs 0.003). | supports |
| 2026-07-13 | session 2026-07-13-003; `tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims20_res0.3.csv` | MG-selected emptyDrops/log-MAD rebuild, 20 PCs, resolution 0.3 | ESPI | Current cluster 7 contains 17 cells; top markers are `Cnmd`, `Scrg1`, `Gja1`, `Ccdc190`, and `Cldn19`, preserving the small strongly separated population under a new cluster ID. | refines |

### Open Questions
- Is cluster 7 reproducible across preprocessing/clustering choices, or is it driven by a few high-signal outlier cells?
- Do cluster 7 cells have QC or sample-origin patterns that explain the strong marker separation?
