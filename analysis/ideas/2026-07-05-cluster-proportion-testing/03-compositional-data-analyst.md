# Compositional Data Analyst Ideas

## Idea 1: Sample-level additive log-ratio models for targeted cluster contrasts

### Title
Sample-level ALR/CLR linear models for E-Stim cluster-composition shifts

### Motivation
Cluster proportions are parts of a whole: if one cluster expands, at least one other cluster must contract even when total cell yield changes. Modeling raw proportions cluster-by-cluster ignores this simplex constraint and can turn closure artifacts into apparent biological effects. A practical first-pass analysis should therefore use Mouse × Condition cluster count vectors and test log-ratio changes, not pooled cells.

### Connection to Existing Data
Use the existing MG-selected cluster assignments and sample metadata to form one count vector per Mouse × Condition sample: control mice 10, 3, and 33; E-Stim mice 10, 3, and 30. The existing pooled Fisher/CLR plot can guide visualization, but not primary inference.

### Approach
- Build a sample × cluster count matrix and keep sample metadata for mouse, condition, pairing status, total cells, and clustering branch.
- Handle zeros with a documented small-count replacement before log-ratio transforms, preferably a multiplicative/simplex-preserving replacement or, for a near-term implementation, a fixed pseudocount such as 0.5 with sensitivity checks over 0.25, 0.5, and 1.
- Compute log-ratio outcomes at the sample level:
  - practical primary: additive log-ratio (ALR), `log2((cluster + replacement) / (reference composition + replacement))`, using a stable reference such as all other MG-selected cells or a biologically stable cluster set;
  - companion display: sample-level CLR coordinates for heatmaps and effect-size plots, interpreted only relative to the sample geometric mean.
- Fit one linear model per tested cluster/log-ratio coordinate with condition as the treatment term and mouse blocking for paired mice 10 and 3 where possible. For the partially paired design, use a model that includes paired within-mouse contrasts for mice 10 and 3 plus unpaired sample information from mouse 30 and mouse 33, rather than discarding either singleton.
- Treat clusters as multiple related tests. Report log-ratio effect sizes, confidence intervals, raw p-values, and FDR-adjusted q-values across clusters. Emphasize that a positive ALR effect means the cluster increased relative to the chosen denominator, not necessarily in absolute cell count.

### Expected Output
A table with one row per cluster containing E-Stim log2-ratio effect, standard error, confidence interval, p-value, q-value, zero count flag, and sensitivity across pseudocount choices. A companion plot would show sample-level log-ratio values with paired mice connected and singleton mice shown separately.

### Feasibility
High. This can be implemented with base R/tidyverse model formulas and the current Mouse × Condition metadata. It avoids cell-level pseudoreplication and produces directly interpretable effect sizes for each cluster.

### Risks/Limitations
- With only six Mouse × Condition samples, inference is fragile and confidence intervals will be wide.
- ALR results depend on the denominator; different reference choices can change which biological contrast looks strongest.
- Zero handling can dominate rare clusters, so sparse clusters need sensitivity labels or filtering.
- Testing clusters defined from the same dataset can be circular if condition-associated structure influenced the clustering. Results should be framed as post hoc composition summaries unless validated with held-out labels, an external reference, or a frozen clustering definition.

### Next Steps
1. Define the sample × cluster count table and pre-specify the cluster denominator/reference.
2. Choose a zero-replacement rule and rare-cluster filter before looking at treatment effects.
3. Fit sample-level log-ratio models and run pseudocount/reference sensitivity checks.
4. Compare the sample-level effect directions to the existing pooled descriptive plot without using the pooled p-values as evidence.

## Idea 2: Dirichlet-multinomial or logistic-normal compositional model with partial pooling

### Title
Hierarchical compositional count model for treatment effects across all clusters

### Motivation
A more robust analysis should model the whole cluster count composition jointly. This respects that all cluster counts in a sample are coupled, accounts for different total cell numbers per sample, and borrows information across clusters while avoiding the false precision of cell-level tests.

### Connection to Existing Data
The same Mouse × Condition count vectors can feed a joint model: each sample contributes total cells and counts over all MG-selected clusters. Pairing is encoded through mouse-level effects for mice observed in both conditions, while mice 30 and 33 contribute as condition-specific singleton samples.

### Approach
- Model the observed cluster count vector for each Mouse × Condition sample as a multinomial-like draw with overdispersion:
  - practical robust option: Dirichlet-multinomial regression with condition effects on cluster composition;
  - more flexible option: logistic-normal multinomial model with cluster-specific treatment effects and mouse-level random effects.
- Use one cluster or the composition center as the reference and estimate treatment effects as log-ratio shifts for every cluster. Interpret each coefficient as the E-Stim change in that cluster relative to the reference balance of the other clusters.
- Encode partial pairing directly: include mouse random intercepts or paired offsets for mice 10 and 3, while allowing unpaired mice 30 and 33 to inform the condition contrast through the population-level condition effect.
- Use shrinkage/partial pooling across cluster treatment effects so rare clusters do not produce unstable extreme estimates solely because of small counts or zeros.
- Handle zeros through the count likelihood rather than ad hoc transformed zero replacement where possible; if the chosen package requires positive transformed values, use a zero-replacement sensitivity analysis and flag clusters whose posterior/effect estimates depend on the replacement.
- Summarize multiple clusters with posterior intervals, local false-sign rates, or FDR-adjusted model-based contrasts. Avoid selecting only clusters that look enriched in the pooled descriptive plot before fitting the model.

### Expected Output
A joint model summary with cluster-specific treatment log-ratio effects, uncertainty intervals, shrinkage-aware rankings, and posterior/sign probabilities or adjusted q-values. A composition-level visualization could show estimated E-Stim and control compositions with uncertainty, plus cluster log-ratio effect intervals.

### Feasibility
Moderate. The design is small, but the model matches the data-generating structure better than separate proportion tests. It may require packages such as `DirichletMultinomial`, `brms`, `rstanarm`, or a custom logistic-normal workflow, so it is better as a robust follow-up after the near-term ALR analysis.

### Risks/Limitations
- Six samples limit the amount of hierarchy the data can identify; priors or shrinkage choices will matter.
- Dirichlet-multinomial assumptions may be too restrictive if only a few clusters drive overdispersion.
- Bayesian/logistic-normal models are harder to explain than simple log-ratio linear models.
- The model still cannot remove circularity from data-derived clusters. It reduces pseudoreplication, but cluster definitions should be frozen before testing, sensitivity-tested against alternative resolutions, or validated on an independent/held-out labeling scheme.

### Next Steps
1. Start from the same pre-specified sample × cluster count matrix used for the practical ALR analysis.
2. Fit a simple Dirichlet-multinomial regression and compare effect directions with the ALR results.
3. If results are promising and diagnostics are acceptable, fit a logistic-normal hierarchical model with weakly informative shrinkage.
4. Report only sample-level composition effects and clearly label any post hoc cluster-selection caveats.

## Recommendation
Prefer Idea 1 as the immediate ESPI analysis: sample-level ALR models are transparent, quick to implement, and align with the Mouse × Condition inferential unit. Use Idea 2 as a robustness analysis if the practical log-ratio results suggest biologically important composition shifts and the team wants a joint model that better respects the full simplex and count uncertainty.
