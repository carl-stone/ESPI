# Phase-Diagram Plateau Criterion

## Persona
**Statistical Physicist** — sees the clustering grid as a finite phase diagram: normalization, PC count, and Leiden resolution are control parameters; cluster count, overlap stability, and small-cluster burden are order parameters; abrupt jumps mark phase boundaries rather than trustworthy biological discoveries.

## Motivation
From this lens, a defensible parameter choice is not the point with the most attractive story; it is a point inside a stable phase. The best clustering should sit on a plateau where small perturbations to preprocessing or dimensionality do not change the macroscopic state. Resolution then behaves like a coupling constant: too low may under-resolve, but too high crosses into a fragmented phase with high susceptibility.

## Connection to Existing Data
The existing `cluster_grid_summary.tsv` is already a 36-row phase diagram over `log1p` vs `pflog`, cell-cycle HVG retained vs filtered, 20/30/50 PCs, and resolutions 0.3/0.5/0.8. Concrete order-parameter signals are available: `n_clusters`, `min_cluster_n`, `n_small_clusters`, `fraction_cells_in_small_clusters`, `ari_vs_reference`, `mean_best_jaccard_to_reference`, and `min_best_jaccard_to_reference`. The observed grid already shows a low-fragmentation plateau at resolution 0.3: all four branches have 10 clusters at 20 PCs and 11 clusters at 30/50 PCs, while resolution 0.5 rises to 12–14 clusters and resolution 0.8 rises to 13–17 clusters. PFlog rows are closer to the current reference than log1p rows, with mean ARI around 0.776 for PFlog versus around 0.560 for log1p, and the top non-reference rows are PFlog-based, including `cluster_pflog_no_filter_cc_dims50_res0.3` with ARI ~0.948 and minimum best Jaccard ~0.787. The clustree 12-panel grid and representative PFlog/filtered/30-PC UMAP sweep can show whether the numeric plateau corresponds to coherent broad structure rather than graph-layout artifacts.

## Approach
1. Treat each grid row as a point in control-parameter space and compute local susceptibilities: finite differences in `n_clusters`, `fraction_cells_in_small_clusters`, `ari_vs_reference`, `mean_best_jaccard_to_reference`, and `min_best_jaccard_to_reference` across adjacent PC counts and adjacent resolutions within each normalization/CC-HVG branch.
2. Add a branch-invariance score for each candidate normalization and PC count by comparing retained vs filtered CC-HVG branches at the same normalization, PCs, and resolution; favor choices whose macroscopic order parameters change least under this preprocessing perturbation.
3. Penalize fragmentation as an instability term: increasing `n_clusters` without a compensating improvement in local overlap stability, plus any rise in `fraction_cells_in_small_clusters`, counts as crossing into a high-susceptibility phase.
4. Select the normalization, PC count, and resolution from the interior of the most stable low-fragmentation plateau, not from a single best reference-overlap row. Operationally, this likely favors PFlog, resolution 0.3, and the lowest PC count that has entered the 11-cluster plateau if the 20-to-30 transition is judged a real stabilization rather than a one-cluster boundary effect.
5. Use the clustree grid and representative UMAP sweep as visual checks of the same label-blind criterion: stable plateaus should show mostly persistent trunks from 0.3 to 0.5, while candidate choices beyond the plateau should show avalanche-like splitting of large manifolds or proliferation of small nodes. Future marker/coherence checks can be pre-declared as post-selection validation or a tie-breaker among plateau-equivalent candidates, without using condition labels.

## Expected Insights
This would convert the current qualitative preference for low resolution and PFlog into an explicit rule: choose the point with low response to perturbations in PCs and CC-HVG policy, low fragmentation, and broad agreement with nearby grid points. It would also distinguish two cases that look similar in a static table: a parameter setting that is stable because it lies in a broad basin versus one that has high ARI only because it is close to the chosen reference row.

## Feasibility
- **Effort**: Medium
- **Data ready**: Yes
- **Methods available**: Standard tools
- **Key risk**: The current grid is coarse, so finite-difference susceptibility may overinterpret a three-level PC axis and three-level resolution axis; if the plateau boundary lies between tested values, the criterion can rank candidates but not precisely locate the transition.

# Metastable Attractor and Split-Entropy Criterion

## Persona
**Statistical Physicist** — treats clusters as metastable attractors in an energy landscape and the Leiden resolution sweep as changing the effective temperature or interaction scale. A useful clustering preserves deep basins while avoiding shallow, noise-driven substates.

## Motivation
In statistical mechanics, not every split of a macrostate is a new phase; many are thermal fluctuations or finite-size artifacts. The analogous clustering criterion is to accept subdivisions only when they are robust attractor basins: they should carry enough mass, persist across nearby PC counts and normalization branches, and later show marker coherence. Resolution selection becomes the search for the last coarse-grained scale before unstable splitting dominates.

## Connection to Existing Data
The 12-panel clustree grid already traces cluster lineages from resolution 0.3 to 0.5 to 0.8 across 20, 30, and 50 PCs and across log1p/PFlog with retained/filtered CC-HVG branches. The representative PFlog, filtered, 30-PC UMAP sweep shows a concrete candidate trajectory: 11 clusters at resolution 0.3, 13 clusters at 0.5, and 15 clusters at 0.8, with higher resolutions subdividing main manifolds and adding small clusters. The table quantifies the same process: high resolution raises the cluster count to 13–17 and the maximum fraction of cells in clusters smaller than 50 cells to ~1.45%, while PFlog resolution-0.3 candidates have low small-cluster burden around 0.67% and strong overlap with the current reference. For example, PFlog filtered 30 PCs at resolution 0.3 has 11 clusters and is the current reference; PFlog filtered 50 PCs at resolution 0.3 also has 11 clusters with ARI ~0.911 and minimum best Jaccard ~0.669, while PFlog no-filter 50 PCs at resolution 0.3 reaches ARI ~0.948 and minimum best Jaccard ~0.787.

## Approach
1. From the cluster assignments underlying the clustree panels, compute a split-entropy table for each parent cluster as resolution increases: parent mass, child mass fractions, largest-child retention, number of children above a minimum size, and entropy of the child-size distribution.
2. Define stable attractor basins as parents whose largest-child retention remains high and whose secondary children are either reproducible across adjacent PC counts/CC-HVG policies or large enough to avoid being treated as finite-size droplets. Define unstable splitting as high split entropy, many small children, or children that appear only in one branch.
3. Score each normalization and PC count by the reproducibility of its attractor tree: PFlog vs log1p and 20/30/50 PCs should be compared on whether the same broad basins persist, not on whether they separate any experimental condition.
4. Choose the Leiden resolution just before the split-entropy curve accelerates. With the current observations, resolution 0.3 is the default candidate because it preserves 10–11 broad basins and low small-cluster burden, whereas 0.5 and 0.8 add subdivisions; a higher resolution would need a pre-declared, label-blind justification that its split children are robust across branches and later marker-coherent.
5. Use future marker/coherence checks as an attractor-depth assay: accepted split children should have internally coherent marker structure and interpretable size, while children lacking coherence should be folded back into their parent for downstream analysis planning.

## Expected Insights
This criterion would say not merely “avoid too many small clusters,” but which splits are physically analogous to real phases versus shallow substates. It could justify choosing a coarse resolution for primary pseudobulk annotation while preserving a ranked list of candidate subclusters whose stability and marker coherence make them worth later inspection.

## Feasibility
- **Effort**: Medium
- **Data ready**: Mostly
- **Methods available**: Needs custom implementation
- **Key risk**: Clustree-derived split entropy depends on having or reconstructing the full cluster-assignment cross-tabs, not just the summary table; if only plotted images are available, the analysis must be regenerated from the existing clustering outputs before it can be scored quantitatively.
