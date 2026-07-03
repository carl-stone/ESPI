# Sensitivity Frontier for Analysis Interventions

## Persona
**Causal Inference Researcher** — treats normalization, PC count, and Leiden resolution as analysis interventions whose effects on cluster assignments should be identified using only pre-declared, label-blind diagnostics.

## Motivation
From a causal perspective, clustering parameters are not passive settings; they are interventions on the analysis pipeline. A defensible choice should sit on a stability frontier: changing nuisance analysis choices nearby should not strongly change the estimand-defining object, namely the cluster partition. This reframes parameter selection as a sensitivity analysis rather than an aesthetic choice. If a resolution or PC count creates clusters only under one nearby specification, it resembles a non-identifiable effect: too dependent on arbitrary analytic intervention to support downstream inference.

## Connection to Existing Data
The current 36-row `cluster_grid_summary.tsv` already contains label-blind outcomes for each analysis intervention: `n_clusters`, `min_cluster_n`, `n_small_clusters`, `fraction_cells_in_small_clusters`, `ari_vs_reference`, `mean_best_jaccard_to_reference`, and `min_best_jaccard_to_reference`. The reviewed patterns are directly useful: PFlog rows are closer to the current reference than log1p rows, with mean ARI about 0.776 versus about 0.560; resolution 0.3 has the highest agreement and least fragmentation; 20 PCs/resolution 0.3 yields 10 clusters across all four preprocessing branches with mean ARI about 0.823; 30 and 50 PCs/resolution 0.3 yield 11 clusters across all four branches with mean ARI about 0.799 and 0.818. The 12-panel clustree grid provides the corresponding split structure from 0.3 to 0.5 to 0.8, and the representative PFlog/filtered/30-PC UMAP sweep shows the move from 11 clusters at 0.3 to 13 and 15 clusters at higher resolutions.

## Approach
1. Treat each grid row as an intervention node in a DAG: normalization, cell-cycle-HVG policy, PC count, and resolution point to the cluster partition; latent technical structure points to both PCA geometry and observed fragmentation; downstream biological analyses are descendants and are not allowed into selection.
2. Build a label-blind stability score for each candidate using the existing table: high `ari_vs_reference`, high `mean_best_jaccard_to_reference`, high `min_best_jaccard_to_reference`, low `n_small_clusters`, low `fraction_cells_in_small_clusters`, and stable `n_clusters` across adjacent PC counts and the retained/filtered branches.
3. Add an executable pairwise extension: compute ARI and best-overlap Jaccard for every pair of candidate cluster columns, then summarize each candidate's average agreement with its local neighborhood, defined as same normalization/resolution with adjacent PC counts and same normalization/PC count with retained versus filtered cell-cycle-HVG policy.
4. Use the clustree grid as a graphical sensitivity diagnostic: penalize candidate resolutions whose 0.3 to 0.5 to 0.8 trajectory shows many thin one-off branches or repeated subdivision of large parent nodes without compensating marker coherence later.
5. Select the lowest-resolution candidate on the broad stability plateau, with a predeclared tie-breaker favoring the branch that maximizes local neighborhood agreement while keeping `fraction_cells_in_small_clusters` near the observed low range.

## Expected Insights
This would turn the current qualitative observation that resolution 0.3 is less fragmented into an explicit causal-sensitivity criterion. It should distinguish whether PFlog/30 PCs/resolution 0.3 is chosen because it is the current reference or because it lies on a wider equivalence plateau that also includes PFlog/50 PCs/resolution 0.3 and possibly PFlog/20 PCs/resolution 0.3. The output would be a transparent statement like: "choose the coarsest partition whose clusters are invariant to adjacent analysis interventions and whose worst matched cluster remains acceptably overlapping."

## Feasibility
- **Effort**: Medium
- **Data ready**: Mostly
- **Methods available**: Standard tools
- **Key risk**: The current table only reports agreement to one reference, so the criterion could inherit reference-choice bias unless the pairwise-neighborhood ARI/Jaccard extension is implemented.

# Negative-Control Nuisance Dependence Screen

## Persona
**Causal Inference Researcher** — separates pre-analysis nuisance variables from downstream biological descendants and uses negative-control logic to reject clusterings that appear driven by technical or preprocessing artifacts.

## Motivation
A causal graph for this decision has nuisance variables upstream of the expression representation: library size, detected genes, mitochondrial fraction, cell-cycle scores, batch-like processing variation if present, and preprocessing branch. These variables can confound the relationship between the observed expression manifold and cluster labels. A good clustering criterion should not ask whether clusters separate the later scientific contrast; it should ask whether clusters are unnecessarily predictable from nuisance causes after the chosen normalization and PC representation. This is a negative-control design: if a high-resolution setting creates clusters that mainly encode technical or cell-cycle nuisance axes, that is evidence against the parameter setting even before marker annotation.

## Connection to Existing Data
The existing grid summary gives the first nuisance-dependence warning signs: higher resolution increases fragmentation, with resolution 0.5 producing 12-14 clusters and resolution 0.8 producing 13-17 clusters, while the maximum fraction of cells in clusters smaller than 50 cells reaches about 1.45%. The 12-panel clustree grid can identify whether those extra clusters arise as repeated splits of broad parents within specific preprocessing branches. The representative PFlog/filtered/30-PC UMAP sweep shows that resolution 0.5 and 0.8 subdivide main manifolds beyond the 11-cluster resolution 0.3 structure. Future marker and coherence checks can be added as label-blind descendants of the selected partition, using cluster marker consistency and not sample labels as evidence.

## Approach
1. For each candidate clustering in the 36-row grid, fit label-blind nuisance-dependence models: predict cluster assignment from QC covariates such as total counts, number of detected genes, mitochondrial percentage, cell-cycle scores, and any available processing covariates that are not scientific labels.
2. Summarize each candidate with a nuisance predictability index, such as cross-validated multinomial deviance improvement or pseudo-$R^2$, and report this alongside `n_clusters`, `fraction_cells_in_small_clusters`, `ari_vs_reference`, and `min_best_jaccard_to_reference` from `cluster_grid_summary.tsv`.
3. Combine this with the clustree grid by flagging splits that appear only at higher resolution and whose child clusters differ mainly on nuisance covariates rather than on broad expression-marker coherence.
4. Use the representative UMAP sweep as a sanity check: for PFlog/filtered/30 PCs, compare whether the additional clusters at resolutions 0.5 and 0.8 have nuisance-heavy signatures or coherent marker programs; do not use scientific condition labels in this comparison.
5. Define a rejection rule before final selection: exclude parameter settings with high nuisance predictability, elevated small-cluster burden, or clustree splits lacking marker coherence; among the remaining candidates, choose the most stable low-fragmentation setting by ARI/Jaccard and local-neighborhood stability.

## Expected Insights
This would identify whether extra resolution is buying real expression structure or merely opening collider-like paths from technical variation into cluster labels. It would also make the cell-cycle-HVG retained versus filtered comparison more principled: if retained branches create clusters whose labels are strongly predicted by cell-cycle scores, filtering has a causal-nuisance justification independent of desired biology. The likely product is a shortlist of clusterings that are stable by ARI/Jaccard and weakly dependent on nuisance variables, with marker coherence reserved as a label-blind confirmation rather than a post hoc rationale.

## Feasibility
- **Effort**: Medium
- **Data ready**: Mostly
- **Methods available**: Standard tools
- **Key risk**: Some genuine cell states may legitimately correlate with QC or cell-cycle covariates, so an overly aggressive nuisance screen could discard biologically meaningful structure; the screen should flag and rank candidates rather than act as an automatic veto.
