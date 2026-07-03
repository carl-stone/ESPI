# Persistence Plateau for Resolution and PC Selection

## Persona
**Topologist / Geometer** — Treats each clustering as a coarse cover of an underlying expression manifold and asks which parameter choices preserve connected components without turning continuous folds into artificial islands.

## Motivation
From a topological lens, a defensible clustering should sit on a persistence plateau: broad components remain recognizable as PCs and resolution vary, while short-lived splits appear only at finer scales. The ESPI grid already contains a small filtration-like sweep over resolutions 0.3 → 0.5 → 0.8 and PC counts 20 → 30 → 50, so the clustree can be read as a persistence diagram in practical form: durable trunks are more trustworthy than branches that appear only at one parameter setting.

## Connection to Existing Data
This is feasible from the current 36-row `cluster_grid_summary.tsv`, the 12-panel `cluster_grid_clustree_12_panel.png`, and the representative `umap_resolution_sweep_pflog_filter_cc_dims30.png`. The relevant existing metrics are `n_clusters`, `min_cluster_n`, `n_small_clusters`, `fraction_cells_in_small_clusters`, `ari_vs_reference`, `mean_best_jaccard_to_reference`, and `min_best_jaccard_to_reference`. Current signals already suggest that resolution 0.3 has the least fragmentation, with 10 clusters at 20 PCs, 11 clusters at 30 and 50 PCs, and high agreement across adjacent PC counts; higher resolutions produce 12–17 clusters and more small-cluster burden.

## Approach
1. Treat each normalization × cell-cycle-HVG policy panel in the clustree as a resolution filtration and score every 0.3 cluster by how cleanly its cells persist into 0.5 and 0.8 descendants, using edge concentration or descendant entropy rather than visual impression alone.
2. For each candidate PC count, compute a topology-stability score from the summary table: reward adjacent-PC agreement in `n_clusters`, high `mean_best_jaccard_to_reference`, and high `min_best_jaccard_to_reference`; penalize jumps in `n_clusters`, `n_small_clusters`, and `fraction_cells_in_small_clusters`.
3. Select resolution by looking for the coarsest level before widespread short-lived branching: in the current grid, explicitly compare whether 0.3 preserves the same broad trunks seen at 0.5/0.8 while avoiding the 12–17-cluster fragmentation seen at higher resolution.
4. Select PC count by plateau behavior rather than maximum dimensionality: require that 20, 30, and 50 PCs give a stable component count and overlap pattern, then prefer the middle of the stable interval if 30 PCs gives the same broad topology as 20/50 without extra fragmentation.
5. Use future marker/coherence checks only after the label-blind choice to audit whether persistent trunks correspond to coherent marker programs; if a short-lived split lacks marker coherence, treat it as geometric over-segmentation rather than a selected resolution.

## Expected Insights
This would turn the current visual clustree review into an executable criterion: choose the normalization, PC count, and resolution whose clusters behave like persistent connected components across the ESPI parameter filtration. It would likely clarify whether `pflog_filter_cc_dims30_res0.3` is a defensible central point on a stable plateau, or whether a neighboring PFlog/resolution-0.3 configuration such as 50 PCs has better minimum Jaccard support without changing the broad manifold partition.

## Feasibility
- **Effort**: Medium
- **Data ready**: Mostly
- **Methods available**: Standard tools
- **Key risk**: The current clustree outputs may not expose enough edge-level counts for a formal entropy score without reading or regenerating the underlying clustering assignments; if only the PNG and summary table are used, part of the criterion remains semi-quantitative.

# Mapper-Style Cover Consistency Across Normalization Branches

## Persona
**Topologist / Geometer** — Views each normalization branch as a different coordinate chart on the same cell-state manifold, then asks whether clustering decisions are invariant under that change of chart.

## Motivation
A good normalization should not dramatically change the intrinsic shape of the data. In differential geometry terms, `log1p` and `pflog` are alternative coordinate systems; the chosen clustering should reflect structures that survive the coordinate change, not artifacts introduced by one chart. Mapper-style thinking asks whether local neighborhoods and broad connected regions are covered consistently across preprocessing choices before committing to a particular partition.

## Connection to Existing Data
The ESPI grid directly supports this because it crosses 2 normalizations, 2 cell-cycle-HVG policies, 3 PC counts, and 3 Leiden resolutions. The summary table reports that PFlog rows are closer to the current reference than log1p rows, with PFlog mean ARI around 0.776 versus log1p around 0.560, and that top non-reference agreement rows are PFlog-based at resolution 0.3. The clustree grid gives branch-wise split structure for log1p retained, log1p filtered, PFlog retained, and PFlog filtered, while the representative PFlog/filtered/30-PC UMAP sweep shows how 0.3, 0.5, and 0.8 progressively subdivide the same embedded geometry.

## Approach
1. Build a branch-invariance table from `cluster_grid_summary.tsv` for each PC × resolution pair: compare `n_clusters`, `ari_vs_reference`, `mean_best_jaccard_to_reference`, `min_best_jaccard_to_reference`, and small-cluster metrics across the four preprocessing branches.
2. Define a coordinate-chart consistency criterion for normalization: prefer the normalization whose retained-vs-filtered cell-cycle-HVG branches have smaller changes in cluster count, lower small-cluster burden, and stronger best-overlap Jaccard across 20/30/50 PCs at the same resolution.
3. Use the 12-panel clustree grid as a Mapper cover diagnostic: mark splits as trustworthy only when the same trunk-to-branch pattern appears under both retained and filtered branches of a normalization, and mark splits as chart artifacts when they appear in only one preprocessing branch or only one PC count.
4. Apply the criterion to resolution by requiring that selected clusters form a stable cover of the manifold: favor resolution 0.3 if it preserves shared trunks across charts while 0.5/0.8 mainly add local subdivisions and increase small-cluster burden.
5. After selecting parameters without labels, use future marker/coherence checks as a geometry audit: persistent clusters should have internally coherent marker structure, while chart-specific or high-resolution fragments should need unusually strong marker evidence to be reconsidered.

## Expected Insights
This would produce a concrete rule for choosing between `log1p` and `pflog`: select the normalization that behaves most like a stable coordinate chart of the same manifold across cell-cycle-HVG policy and PC count. It would also explain PC and resolution selection in geometric terms: a chosen setting is not simply close to the current reference, but preserves an invariant cover of the data across nearby parameterizations while avoiding small artificial islands.

## Feasibility
- **Effort**: Low to Medium
- **Data ready**: Yes
- **Methods available**: Standard tools
- **Key risk**: Because the current agreement metrics are all measured against a stated reference, the chart-consistency score could inherit reference bias unless paired with direct branch-to-branch overlap calculations or clearly framed as a first-pass summary-table criterion.
