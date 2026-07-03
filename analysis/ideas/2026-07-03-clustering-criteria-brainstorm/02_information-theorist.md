# Minimum Description Length Score for Grid Partitions

## Persona
**Information Theorist** — treats each candidate clustering as a lossy code for 5,538 cell-level expression states, asking which normalization / PC / resolution setting gives the shortest honest description without spending bits on unstable fragmentation.

## Motivation
From an information-theory view, a clustering is useful when it compresses the expression manifold: many cells can share a compact cluster label, but the labels still preserve reproducible structure. Too few clusters underfit and throw away information; too many clusters spend extra code length on splits that are not stable across nearby preprocessing choices. This suggests a label-blind minimum-description-length criterion: prefer the partition that is cheap to encode, stable under small perturbations of the analysis pipeline, and not dominated by singleton-like exception codes.

## Connection to Existing Data
The existing `cluster_grid_summary.tsv` already gives the core ingredients for a first MDL-style score across the 36-row grid: `n_clusters`, `min_cluster_n`, quartiles of cluster size, `n_small_clusters`, `fraction_cells_in_small_clusters`, `ari_vs_reference`, `mean_best_jaccard_to_reference`, and `min_best_jaccard_to_reference`. It shows that resolution 0.3 is much less fragmented than 0.5 or 0.8: at 20 PCs / res 0.3 all four branches have 10 clusters, at 30 and 50 PCs / res 0.3 all four branches have 11 clusters, while res 0.8 ranges from 13 to 17 clusters. The small-cluster burden remains low but rises at high resolution, with a maximum fraction of cells in clusters smaller than 50 cells of about 1.45%. The PFlog rows are closer to the current reference than log1p rows, with PFlog mean ARI across grid rows around 0.776 versus log1p around 0.560, and the strongest non-reference rows are PFlog-based, including `cluster_pflog_no_filter_cc_dims50_res0.3` with 11 clusters, ARI about 0.948, and minimum best Jaccard about 0.787.

## Approach
1. Start with `cluster_grid_summary.tsv` and compute, for each grid row, the empirical cluster-label entropy `H(C) = -sum_k p_k log2(p_k)` from the candidate cluster sizes; use this as the per-cell code length for storing the partition.
2. Add an exception-code penalty for fragmentation: `lambda_small * fraction_cells_in_small_clusters` plus a monotone penalty for `n_small_clusters`, using the existing `<50 cells` summary so that res 0.8 rows with 70–80 cells in small clusters are penalized more than res 0.3 rows with 37–47 cells in small clusters.
3. Add a redundancy/stability reward from existing overlap metrics: high `mean_best_jaccard_to_reference` and high `min_best_jaccard_to_reference` reduce the score, while low minimum best Jaccard marks clusters that are hard to transmit across reference/candidate partitions. Treat the current reference as one alignment anchor, not biological truth.
4. For each normalization × PC count × resolution, compute an MDL-like score such as `H(C) + fragmentation_penalty + instability_penalty`; then summarize score ranks by normalization, by PC count, and by resolution. The expected first-pass criterion is to retain settings near the minimum score plateau rather than the absolute minimum if the plateau favors simpler structure.
5. Use the 12-panel clustree grid as a qualitative audit of the same code-length idea: choose the lowest resolution before large unstable splits appear, and check whether the representative PFlog / filtered / 30-PC UMAP sweep shows res 0.5 or 0.8 buying visibly coherent information rather than just extra labels.

## Expected Insights
This would turn the current informal preference for low-fragmentation, high-stability configurations into a pre-declared numeric rule. It should tell whether the information-efficient region is broad around PFlog with 20–50 PCs at resolution 0.3, or whether a specific 50-PC PFlog row earns its extra complexity by preserving enough overlap information. It would also make explicit the tradeoff between cluster entropy and robustness: a row with more clusters must justify every extra bit by increasing stable overlap, not merely by subdividing the manifold.

## Feasibility
- **Effort**: Medium
- **Data ready**: Mostly
- **Methods available**: Standard tools
- **Key risk**: The current table gives summary size statistics but not the full cluster-size vector or all pairwise partition overlaps, so the first implementation may approximate the code length unless it reloads the candidate cluster assignments.

# Variation-of-Information Stability Map Across Neighboring Grid Settings

## Persona
**Information Theorist** — sees normalization, PC count, and Leiden resolution as noisy channel settings, where the chosen parameter set should transmit the same latent structure through adjacent analysis perturbations with minimal information loss.

## Motivation
Adjusted Rand index and best-overlap Jaccard are useful, but an information theorist would ask a slightly different question: how many bits are needed to translate one clustering into another? Variation of Information, `VI(A,B) = H(A|B) + H(B|A)`, directly measures partition instability as residual uncertainty after one clustering is known. A defensible clustering criterion is therefore to select a parameter region with low local VI to neighboring PC counts, CC-HVG policies, and normalizations, while rejecting resolutions that increase cluster entropy without reducing translation uncertainty.

## Connection to Existing Data
The 36-row grid is explicitly structured for this analysis: two normalizations (`log1p`, `pflog`), two CC-HVG policies, three PC counts (20, 30, 50), and three resolutions (0.3, 0.5, 0.8). The existing summary table already shows the local-stability signal indirectly: res 0.3 gives consistent cluster counts across branches at each PC count (10 clusters at 20 PCs; 11 clusters at 30 and 50 PCs), while res 0.5 and 0.8 introduce 12–17 clusters and more small-cluster burden. The top non-reference agreement rows are PFlog-based, especially `cluster_pflog_no_filter_cc_dims50_res0.3` and `cluster_pflog_filter_cc_dims50_res0.3`, and the clustree grid shows each branch's 0.3 → 0.5 → 0.8 split structure. The representative UMAP sweep for PFlog / filtered / 30 PCs gives a visual check that res 0.5 and res 0.8 subdivide main manifolds beyond the 11-cluster res 0.3 baseline.

## Approach
1. Reload the candidate cluster columns for all 36 configurations and compute pairwise Variation of Information and Normalized Mutual Information for every pair of clusterings; keep the analysis label-blind by using only cluster labels and cell IDs.
2. Define local neighborhoods in parameter space: same normalization / CC-HVG policy / resolution with adjacent PC counts; same normalization / PC count / resolution across CC-HVG retained vs filtered; same PC count / CC-HVG policy / resolution across log1p vs PFlog; and same branch / PC count across adjacent resolutions 0.3→0.5 and 0.5→0.8.
3. For each grid row, compute a local channel-stability score: low mean VI to adjacent PC and CC-HVG neighbors, low worst-neighbor VI, and high NMI to its closest preprocessing neighbor. Require that a preferred normalization and PC count lie inside a stable basin, not as a single isolated high-ARI row.
4. For resolution, compute the information gain from each split step as `H(C_high_res) - H(C_low_res)` and compare it with the VI cost between the two resolutions. A split is accepted only if it adds stable mutual information across neighboring branches; otherwise, it is treated as overcoding.
5. Cross-check the resulting basin against the clustree grid and the representative UMAP sweep, then reserve future marker/coherence checks as a label-blind validation of whether accepted splits have internally coherent marker programs rather than arbitrary extra bits.

## Expected Insights
This would produce a stability map rather than a single-reference ranking. It could justify PFlog if PFlog rows form a lower-VI basin than log1p, choose 30 or 50 PCs depending on which has lower local translation cost, and choose resolution 0.3 if the 0.5 and 0.8 splits have high VI cost relative to their entropy gain. It would also reveal whether the strong ARI of `cluster_pflog_no_filter_cc_dims50_res0.3` is part of a stable channel-capacity plateau or merely close to the current reference by construction.

## Feasibility
- **Effort**: Medium
- **Data ready**: Mostly
- **Methods available**: Standard tools
- **Key risk**: VI/NMI require access to the full per-cell cluster assignments for every grid row; if only summary tables are available, the method cannot distinguish stable large-cluster refinements from unstable small-cluster swaps.
