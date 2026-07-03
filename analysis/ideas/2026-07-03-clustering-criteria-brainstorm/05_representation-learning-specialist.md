# Invariance-First Latent Partition Selection

## Persona
**Representation Learning Specialist** — evaluates clustering parameters by whether they recover stable latent factors of variation rather than artifacts of a particular preprocessing view.

## Motivation
A useful representation should preserve the same coarse latent partition when nuisance choices change. Here, normalization method, cell-cycle HVG policy, PC truncation, and Leiden resolution are all different views or bottlenecks on the same cell manifold. From a representation-learning lens, the best clustering is not the one that maximizes any sample-label separation; it is the one whose assignments are invariant under label-blind perturbations of the representation while still retaining enough capacity to separate real modes.

## Connection to Existing Data
The current 36-row `cluster_grid_summary.tsv` already contains `n_clusters`, `min_cluster_n`, `n_small_clusters`, `fraction_cells_in_small_clusters`, `ari_vs_reference`, `mean_best_jaccard_to_reference`, and `min_best_jaccard_to_reference` for the full grid: `log1p` vs `pflog`, CC-HVG retained vs filtered, 20/30/50 PCs, and resolutions 0.3/0.5/0.8. The observed table patterns are directly relevant: PFlog rows are closer to the current reference than log1p rows, resolution 0.3 gives the least fragmentation, 20-PC/resolution-0.3 rows consistently produce 10 clusters, 30- and 50-PC/resolution-0.3 rows consistently produce 11 clusters, and resolution 0.8 increases the cluster count to 13–17 with a higher small-cluster burden. The 12-panel clustree grid adds visual evidence about whether splits are stable hierarchical refinements or brittle fragmentations, and the representative PFlog/filter-CC/30-PC UMAP sweep shows 11, 13, and 15 clusters as resolution increases.

## Approach
1. Reframe the existing grid summary as a label-blind invariance table: for each candidate normalization × PC count × resolution, compute how stable its clustering is against adjacent representation choices, not just against the current reference column. Use pairwise ARI and best-overlap Jaccard across matched candidates that differ by one axis: 20↔30 PCs, 30↔50 PCs, retained↔filtered CC-HVG policy, and `log1p`↔`pflog` where the PC count and resolution are fixed.
2. Define an explicit representation-invariance score, for example: median adjacent-pair ARI + median `mean_best_jaccard` − penalty for low `min_best_jaccard` − penalty for `fraction_cells_in_small_clusters` and `n_small_clusters`. Keep all inputs label-blind and report the component scores rather than hiding them in a black-box rank.
3. Use the existing clustree grid to qualify the score: prefer candidates whose 0.3→0.5→0.8 transitions look like nested refinement of large nodes, and demote candidates where the same parent repeatedly sheds tiny unstable leaves. This uses edge structure and node size from the clustree panels, not graph-layout position.
4. Apply a decision rule: choose the lowest resolution that reaches a high invariance plateau across adjacent PC counts and CC-HVG policies. With the currently observed signals, this rule would likely favor a PFlog-based, resolution-0.3 solution, then use 20 vs 30 vs 50 PCs as the remaining question: pick the PC count whose adjacent-pair stability is high without losing the extra broad state implied by the 10-cluster vs 11-cluster difference.
5. After the label-blind choice is fixed, run marker/coherence checks as an audit rather than selection evidence: each chosen cluster should have internally coherent marker programs and should not exist only as a tiny split with poor best-overlap Jaccard across adjacent representations.

## Expected Insights
This would turn “PFlog/resolution 0.3 looks stable” into a defensible criterion: the selected parameters should be those whose clusters are invariant to reasonable perturbations of the latent representation. It would also separate two different questions that are currently entangled: whether `pflog` gives a more robust representation than `log1p`, and whether 30 or 50 PCs add meaningful latent structure beyond 20 PCs instead of just increasing sensitivity to noise.

## Feasibility
- **Effort**: Medium
- **Data ready**: Mostly
- **Methods available**: Standard tools
- **Key risk**: The current summary table reports agreement to one reference clustering, so true all-pairs invariance may require regenerating or extending the summary to compute pairwise ARI/Jaccard across all 36 candidate cluster columns; relying only on `ari_vs_reference` could bake in the current reference choice.

# Factor-Basis Cluster Coherence Audit

## Persona
**Representation Learning Specialist** — asks whether clusters correspond to compact, interpretable factors in latent space, with separable axes and coherent local neighborhoods rather than arbitrary graph cuts.

## Motivation
Clusters should act like discrete regions carved out of a meaningful continuous representation. If changing PC count or resolution creates clusters that are not compact in PC/UMAP space, lack coherent marker directions, or appear only as thin slices of a broader manifold, then those clusters are likely over-parameterized graph artifacts. A representation-learning criterion can therefore combine latent-space compactness, split coherence, and post hoc marker-program coherence while remaining blind to sample labels.

## Connection to Existing Data
The ESPI grid already spans 20, 30, and 50 PCs, so it can test whether added latent dimensions produce stable factors or noisy axes. The representative UMAP sweep for `pflog_filter_cc_dims30` shows a concrete coarse-to-fine path: 11 clusters at resolution 0.3, 13 at 0.5, and 15 at 0.8, with higher resolutions subdividing main manifolds. The grid summary shows the same fragmentation pressure numerically: resolution 0.5 produces 12–14 clusters and resolution 0.8 produces 13–17 clusters, while small-cluster burden remains low but peaks around 1.45% of cells. The clustree grid provides the split graph needed to ask whether new clusters are coherent latent factors or unstable fragments.

## Approach
1. For each candidate clustering in the 36-row grid, compute label-blind latent compactness within the same representation used for clustering: within-cluster dispersion in the selected PC space, nearest-neighbor purity, and separation from adjacent clusters. Normalize these scores within each PC count so that 50-PC candidates are not rewarded or punished merely for having more dimensions.
2. Build a split-coherence table from the clustree grid data: for each parent cluster at resolution 0.3, record whether higher-resolution children have substantial size, high parent-to-child Jaccard, and clear separation in PC/UMAP coordinates. Penalize splits that create small children without improving compactness or separation.
3. Use the current grid summary metrics as hard guardrails: reject candidates with elevated `n_small_clusters` or `fraction_cells_in_small_clusters`, reject candidates with poor `min_best_jaccard_to_reference` unless they are strong across the all-pairs invariance audit, and prefer rows near the stable cluster-count plateau seen at resolution 0.3 over rows where resolution alone inflates cluster count.
4. Define a predeclared marker/coherence audit that runs only after choosing a short list label-blind: for each candidate cluster, test whether top marker genes form coherent expression programs within the cluster and whether those marker programs are stable across adjacent PC counts or CC-HVG retained/filtered branches. Keep the audit restricted to intrinsic cluster coherence and explicitly exclude sample-label or treatment-derived contrasts from the scoring rule.
5. Select parameters by the simplest factor basis that passes all audits: enough PCs to represent the stable broad factors, a normalization whose clusters remain compact and invariant across CC-HVG policy, and the lowest Leiden resolution before split-coherence starts declining. Given the existing outputs, this would likely nominate PFlog with resolution 0.3 and then adjudicate 20 vs 30 vs 50 PCs by whether the 11-cluster solutions at 30/50 PCs show an additional coherent factor absent from the 10-cluster 20-PC solution.

## Expected Insights
This would reveal whether the extra clusters introduced by higher resolution or higher PC count behave like meaningful latent factors or like fragmented graph partitions. It would also provide a transparent way to justify keeping a slightly richer representation, such as 30 or 50 PCs, only if the extra cluster has compact latent geometry and stable marker coherence rather than merely improving agreement with the current reference.

## Feasibility
- **Effort**: Medium
- **Data ready**: Mostly
- **Methods available**: Standard tools
- **Key risk**: Compactness and marker coherence can favor broad, easy-to-separate states and underweight rare but real populations; the criterion needs explicit small-cluster guardrails rather than an automatic bias against every rare cluster.
