---
topic: population-structure
description: Empirical findings about marker-defined cellular structure within focused analysis branches. Entries treat cluster marker rankings as descriptive cell-level evidence unless paired with replicate-level tests.
created: 2026-07-04
last_updated: 2026-07-17
status: active
---

# Population Structure

## F-003: The prior marker-negative cluster 2 was run-specific
**Status:** superseded
**Claim:** The prior 30-PC MG-selected clustering contained a cluster 2 without retained positive markers. That result does not hold in the current 20-PC/resolution-0.5 clustering, where cluster 2 has positive markers.
**Implications:** Treat the earlier marker-negative cluster 2 as a run-specific historical result, not a current interpreted population.
**Tags:** MG-selected, clustering, markers, interpreted identity, negative result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-04 | notebook/sc_analysis.qmd cluster marker section; tables/mg_selected/find_all_markers_summary_data_pflog_mg_selected_no_filter_cc_dims30_res0.3.csv | MG-selected chosen clustering, PFlog branch, 30 PCs, resolution 0.3 | ESPI | Marker summary reports cluster 2 with 961 cells, `n_retained_markers = 0`, and `n_top_markers = 0`; notebook notes this argues against treating cluster 2 as a marker-defined identity. | supports |
| 2026-07-13 | session 2026-07-13-003; `tables/mg_selected/find_all_markers_summary_data_pflog_mg_selected_no_filter_cc_dims20_res0.3.csv` | MG-selected emptyDrops/log-MAD rebuild, 20 PCs, resolution 0.3 | ESPI | Current cluster 2 contains 937 cells and 894 retained markers, including `Aldh1a1`, `Nxph1`, `Nr2f2`, `Efnb2`, and `Tes`; this invalidates the prior marker-negative cluster-number claim for the current run. | supersedes |
| 2026-07-13 | session 2026-07-13-007; `tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims20_res0.5.csv` | Current MG-selected emptyDrops/scDblFinder/log-MAD rebuild, 20 PCs, resolution 0.5 | ESPI | Cluster 2 contains 796 cells; its top markers are `Bcat1`, `Htra1`, `Car14`, `Gbe1`, and `Etnppl`, confirming that the earlier marker-negative cluster-number result remains superseded. | supersedes |
| 2026-07-17 | session 2026-07-17-001; `tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims20_res0.5.csv` | QC-corrected MG-selected branch, 20 PCs, resolution 0.5 | ESPI | Cluster 2 contains 727 cells and has positive top markers `Car14`, `Bcat1`, `Nudt4`, `Htra1`, and `Espn`, so the prior marker-negative result remains superseded. | supersedes |

### Open Questions
- What biological state does current cluster 2 represent, and is it stable across normalization branches and resolutions?
- How do its marker profile and sample composition compare with neighboring clusters?

## F-004: Cluster 4 carries a neurogenic progenitor marker program
**Status:** preliminary
**Claim:** In the current MG-selected chosen clustering, cluster 4 is enriched for neurogenic-progenitor marker genes, including `Ascl1` and `Hes6`, in descriptive cell-level marker rankings.
**Implications:** Cluster 4 is the strongest marker-ranked neurogenic-progenitor-like subpopulation in the selected branch, but this does not by itself show a Mouse × Condition condition effect.
**Tags:** MG-selected, clustering, markers, neurogenic progenitor, cell-level ranking

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-04 | notebook/sc_analysis.qmd cluster marker section; tables/mg_selected/find_all_markers_summary_data_pflog_mg_selected_no_filter_cc_dims30_res0.3.csv; tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims30_res0.3.csv | MG-selected chosen clustering, PFlog branch, 30 PCs, resolution 0.3 | ESPI | Cluster 3 has 446 cells and 2,614 retained positive markers; top-five marker rows include `Ascl1` (avg_log2FC 2.80, pct.1 0.881 vs pct.2 0.137), `Hes6` (2.57, 0.975 vs 0.436), and `Dll1` (3.11, 0.751 vs 0.162). | supports |
| 2026-07-13 | session 2026-07-13-003; `tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims20_res0.3.csv` | MG-selected emptyDrops/log-MAD rebuild, 20 PCs, resolution 0.3 | ESPI | Current cluster 3 contains 453 cells and retains the neurogenic program: `Ascl1` (avg_log2FC 3.32), `Hes6` (2.67), and `Dll1` (2.98) are enriched. | supports |
| 2026-07-13 | session 2026-07-13-007; `tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims20_res0.5.csv` | Current MG-selected emptyDrops/scDblFinder/log-MAD rebuild, 20 PCs, resolution 0.5 | ESPI | Cluster 4 contains 470 cells; top markers include `Chrna4`, `Ascl1`, `Miat`, `Hes6`, and `Ncald`, preserving the descriptive neurogenic-progenitor program under a new cluster ID. | refines |
| 2026-07-17 | session 2026-07-17-001; `tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims20_res0.5.csv` | QC-corrected MG-selected branch, 20 PCs, resolution 0.5 | ESPI | Cluster 4 contains 443 cells; top markers include `Ascl1`, `Chrna4`, `Miat`, `Hes6`, and `Qsox1`, preserving the descriptive neurogenic-progenitor program. | refines |

### Open Questions
- Does the cluster 4 marker program correspond to a stable interpreted identity across normalization branches or cluster resolutions?
- How much of the cluster 4 program overlaps the primary DE `Hes6` signal versus cell-composition or within-cell abundance changes?

## F-005: Cluster 8 is a very small outlier with strong marker separation
**Status:** preliminary
**Claim:** In the current MG-selected chosen clustering, cluster 8 is a 15-cell outlier with strong descriptive marker separation, including `Cnmd`, `Scrg1`, `Gja1`, `Ccdc190`, and `Cldn19`.
**Implications:** Cluster 8 may represent a rare state or outlier population, but its small size makes it unsuitable for broad identity or condition-effect claims without additional validation.
**Tags:** MG-selected, clustering, markers, outlier, rare population

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-04 | notebook/sc_analysis.qmd cluster marker section; tables/mg_selected/find_all_markers_summary_data_pflog_mg_selected_no_filter_cc_dims30_res0.3.csv; tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims30_res0.3.csv | MG-selected chosen clustering, PFlog branch, 30 PCs, resolution 0.3 | ESPI | Cluster 8 has 15 cells and 6,825 retained positive markers; top-five rows show high marker separation for `Htr1b` (avg_log2FC 9.92, pct.1 0.467 vs pct.2 0), `Flt4` (9.52, 0.600 vs 0.001), `Cnmd` (8.99, 1.000 vs 0.006), `Cldn10` (8.91, 0.867 vs 0.004), and `Scrg1` (8.67, 0.867 vs 0.003). | supports |
| 2026-07-13 | session 2026-07-13-003; `tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims20_res0.3.csv` | MG-selected emptyDrops/log-MAD rebuild, 20 PCs, resolution 0.3 | ESPI | Current cluster 7 contains 17 cells; top markers are `Cnmd`, `Scrg1`, `Gja1`, `Ccdc190`, and `Cldn19`, preserving the small strongly separated population under a new cluster ID. | refines |
| 2026-07-13 | session 2026-07-13-007; `tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims20_res0.5.csv` | Current MG-selected emptyDrops/scDblFinder/log-MAD rebuild, 20 PCs, resolution 0.5 | ESPI | Cluster 8 contains 17 cells; top markers are `Cnmd`, `Scrg1`, `Gja1`, `Ccdc190`, and `Cldn19`, preserving the small strongly separated population under its current cluster ID. | refines |
| 2026-07-17 | session 2026-07-17-001; `tables/mg_selected/find_all_markers_top5_data_pflog_mg_selected_no_filter_cc_dims20_res0.5.csv` | QC-corrected MG-selected branch, 20 PCs, resolution 0.5 | ESPI | Cluster 8 contains 15 cells; top markers are `Gja1`, `Cldn19`, `Ccdc190`, `Scrg1`, and `Cnmd`, preserving the small strongly separated population. | refines |

### Open Questions
- Is cluster 8 reproducible across preprocessing/clustering choices, or is it driven by a few high-signal outlier cells?
- Do cluster 8 cells have QC or sample-origin patterns that explain the strong marker separation?
