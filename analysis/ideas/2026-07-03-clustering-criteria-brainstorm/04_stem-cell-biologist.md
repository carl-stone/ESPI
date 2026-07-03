# Fate-State Persistence Score for Choosing the Coarsest Stable Manifold

## Persona
**Stem Cell Biologist** — Applies a stem-cell-state lens: pluripotency networks, lineage priming, cellular memory, commitment thresholds, and bistable transitions should appear as reproducible expression-state basins rather than parameter-induced shards.

## Motivation
Stem-cell systems often contain metastable attractor-like states: broad self-renewal, primed, transitional, and lineage-leaning states can coexist, but true commitment points should be robust to modest analysis perturbations. From this lens, the best clustering parameters are not the ones that create the most subtypes; they are the parameters that preserve broad fate basins across preprocessing choices while avoiding artificial fragmentation of continuous priming trajectories. A defensible clustering should capture persistent cell-state memory without over-interpreting every local density wrinkle as a fate decision.

## Connection to Existing Data
ESPI already has the needed label-blind grid outputs: the 36-row `cluster_grid_summary.tsv`, the 12-panel clustree grid spanning 20/30/50 PCs, log1p/PFlog, retained/filtered cell-cycle HVGs, and resolution 0.3/0.5/0.8, plus the representative PFlog filtered 30-PC UMAP sweep. The current summary shows that resolution 0.3 has the least fragmentation and high agreement across branches: 10 clusters at 20 PCs, 11 clusters at 30 and 50 PCs, with mean ARI around 0.799–0.823 depending on PC count. Higher resolutions create 12–17 clusters and raise the small-cluster burden; the UMAP sweep similarly shows res 0.5 and 0.8 subdividing main manifolds beyond the 11-cluster coarse structure.

## Approach
1. Start from `cluster_grid_summary.tsv` and assign each candidate a fate-state persistence score combining: high ARI to neighboring PC counts within the same preprocessing branch, high mean best Jaccard to the stated reference, low `n_small_clusters`, low `fraction_cells_in_small_clusters`, and stable `n_clusters` across CC-HVG retained vs filtered branches.
2. Treat clustree splits from 0.3 to 0.5 to 0.8 as putative fate-branching events; score a split as credible only if the parent-to-child structure is preserved across adjacent PC counts and both CC-HVG policies, not if it appears only in one normalization/PC/resolution corner.
3. Prefer the lowest Leiden resolution whose clusters are stable across adjacent PC counts and preprocessing branches; operationally, require the chosen resolution to avoid a step increase in `n_clusters` without a matching improvement in minimum best Jaccard or clustree-consistent split reproducibility.
4. Use the PFlog/log1p contrast as a normalization robustness check: if PFlog rows consistently show higher ARI and Jaccard while retaining the same broad cluster count and lower fragmentation, choose PFlog; if log1p preserves a split that PFlog loses, require later label-blind marker coherence before accepting that split.
5. After parameter choice is frozen, run marker/coherence checks for each chosen cluster: stem-state, priming, differentiation, cell-cycle, stress, and mitochondrial/ribosomal signatures can annotate or challenge clusters, but not retroactively select the parameters.

## Expected Insights
This would produce a concrete selection rule: choose the normalization, PC count, and resolution that maximize fate-basin persistence while minimizing fragile subclustering. It should clarify whether the apparent 11-cluster structure at PFlog, filtered cell-cycle HVGs, 30 or 50 PCs, res 0.3 behaves like a stable set of broad cellular-memory states, and whether the extra clusters at res 0.5/0.8 look like reproducible commitment branches or analysis-induced fragmentation.

## Feasibility
- **Effort**: Medium
- **Data ready**: Mostly
- **Methods available**: Standard tools plus small custom scoring from existing grid metrics and clustree edges
- **Key risk**: A score anchored to the current reference clustering can inherit reference bias; the criterion should therefore emphasize adjacent-PC and preprocessing-branch stability, not only ARI to the reference.

# Lineage-Priming Coherence Gate for Accepting Finer Resolution

## Persona
**Stem Cell Biologist** — Uses lineage priming, epigenetic memory, reprogramming barriers, and commitment points as a model for deciding when a finer partition is biologically coherent enough to keep.

## Motivation
In stem-cell biology, a meaningful split should resemble a commitment or priming boundary: daughter states should have internally coherent regulatory programs and interpretable transitions from a parent state. If a higher Leiden resolution subdivides a manifold, those new clusters should show coherent lineage-priming or cell-state-memory signatures later, not just small differences in local graph density. This idea treats higher resolution as guilty until proven coherent by label-blind biology.

## Connection to Existing Data
The current grid already shows a natural test case. In `cluster_grid_summary.tsv`, res 0.3 has 10–11 clusters across PC counts and branches, while res 0.5 gives 12–14 clusters and res 0.8 gives 13–17 clusters with higher small-cluster burden, up to about 1.45% of cells in clusters under 50 cells. The representative PFlog filtered 30-PC UMAP sweep shows res 0.3 as a coarse 11-cluster structure, res 0.5 subdividing main manifolds into 13 clusters, and res 0.8 further subdividing into 15 clusters. The clustree grid provides the parent-child split map needed to ask which finer clusters are stable daughters rather than ephemeral shards.

## Approach
1. Define a pre-declared split-acceptance gate for moving from res 0.3 to 0.5 or 0.8: a finer split is accepted only if it has stable parent-child mapping in the clustree grid, adequate child size using `min_cluster_n` and `fraction_cells_in_small_clusters`, and high best-overlap Jaccard to a corresponding cluster in neighboring PC/preprocessing branches.
2. For each candidate normalization and PC count, compute a split ledger: parent cluster, child clusters, child sizes, Jaccard consistency across adjacent PC counts, whether the split appears under both CC-HVG policies, and whether it is unique to one normalization.
3. Choose normalization and PC count by the branch whose split ledger is most parsimonious: enough PCs to preserve stable coarse states, but not so many that res 0.5/0.8 creates one-off daughter states absent from 20/30/50-PC neighbors. Current grid patterns make PFlog with 30 or 50 PCs and res 0.3 the likely default candidate to challenge.
4. Freeze a primary coarse clustering at the lowest passing resolution, then reserve finer clusters as provisional only if later marker/coherence analysis shows internally coherent stem-state, priming, differentiation, cell-cycle, stress, or metabolic programs that align with the clustree split rather than merely recapitulating technical or size artifacts.
5. Document the rule as a decision table: “accept finer resolution only when clustree reproducibility + Jaccard overlap + size floor + post hoc label-blind marker coherence all pass”; otherwise collapse to the parent resolution for downstream analysis.

## Expected Insights
This would separate two questions that are often conflated: which coarse clustering is stable enough for primary analysis, and which higher-resolution splits deserve interpretation as candidate lineage-priming or commitment substates. It would make res 0.5/0.8 useful as a discovery audit while preventing small or unstable daughter clusters from steering the main clustering choice.

## Feasibility
- **Effort**: Medium
- **Data ready**: Mostly
- **Methods available**: Standard tools, using existing grid metrics, clustree parent-child structure, UMAP sweep review, and later marker/coherence summaries
- **Key risk**: True rare commitment states may be small and could fail a size or reproducibility gate; the rule should allow them to be flagged for follow-up marker review without letting them determine the primary parameter selection.
