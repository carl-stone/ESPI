# Causal Inference Researcher

## Idea 1: Design-Based Mouse-Level Proportion Contrasts With Pre-Specified Clusters

### Title
Design-based cluster proportion contrasts using Mouse × Condition samples.

### Motivation
The causal question is whether E-Stim changes the distribution of MG-selected cells across clusters. The estimand should live at the mouse-sample level, not the cell level: for cluster $k$, the average treatment effect on the mouse-level cluster proportion, $E[p_{ik}(1) - p_{ik}(0)]$, where $p_{ik}(a)$ is the fraction of MG-selected cells in mouse $i$ that would fall in cluster $k$ under condition $a$.

### Connection to Existing Data
Use one count vector per Mouse × Condition sample. The design has two paired mice, 10 and 3, plus one E-Stim-only mouse 30 and one control-only mouse 33. The existing pooled Fisher/CLR output can nominate patterns to inspect, but it should not define the inferential claim.

### Approach
- Freeze the cluster definition before testing: use the current MG-selected cluster labels as an explicit, pre-specified measurement rule, and state that the test targets these observed clusters rather than undiscovered cell states.
- For each Mouse × Condition sample, compute $n_{ik}$ cells in cluster $k$, total MG-selected cells $N_i$, and proportion $p_{ik} = n_{ik}/N_i$.
- Primary contrast: estimate the condition effect for each cluster using mouse-level proportions.
  - Paired component: within-mouse differences for mice 10 and 3.
  - Unpaired component: condition means for mouse 30 and mouse 33, reported as supportive because each condition has only one unpaired sample.
- Fit a small-sample design-based model per cluster, such as a linear model on transformed proportions with condition plus a paired-mouse blocking term for mice 10 and 3 and explicit weights or precision notes from $N_i$.
- Report effect sizes first: absolute percentage-point change and log-ratio or centered log-ratio change at the mouse-sample level. Use confidence intervals from permutation/sign-flip or bootstrap over mouse samples only, not over cells.
- Treat multiplicity descriptively: rank clusters by effect size and uncertainty rather than making strong discovery claims from six samples.

### Expected Output
A table per cluster with control mean, E-Stim mean, paired differences for mice 10 and 3, overall estimated E-Stim minus control contrast, uncertainty interval, and a clear flag that inference uses Mouse × Condition samples.

### Feasibility
High. This only needs the existing cluster labels, mouse IDs, condition labels, and cluster counts. It can be implemented later in R/tidyverse without changing the clustering pipeline.

### Risks/Limitations
- Threats to validity include low biological replication, partial pairing, compositional dependence, and post-selection cluster definitions.
- Very small sample size limits formal testing and makes intervals wide.
- Partial pairing means the estimand blends strong within-mouse evidence from mice 10 and 3 with weak between-mouse information from mice 30 and 33.
- Cluster proportions are compositional, so an increase in one cluster can force decreases elsewhere.
- If clusters were learned in a way that used condition-associated structure, the estimand is partly post-selection: it asks whether E-Stim changes membership in these data-derived clusters, not whether E-Stim changes an externally validated cell type.

### Next Steps
- Define the exact frozen cluster-label input and sample-level count table.
- Choose one transformation and one primary effect-size scale before seeing results.
- Write the result text to emphasize design-based effect estimation, not cell-level hypothesis testing.
- Keep the pooled Fisher/CLR plot as descriptive context only.

## Idea 2: Cross-Fit Cluster Definition and Sensitivity Analysis for Post-Selection Validity

### Title
Cross-fit and sensitivity-check cluster proportion effects under alternative cluster definitions.

### Motivation
The strongest threat is circularity: clusters learned from the same cells may encode condition differences, so testing E-Stim effects on those clusters can overstate evidence. The target estimand is a validation estimand: for a cluster-assignment rule $g^{(-s)}$ trained without sample $s$ or without the held-out condition contrast, estimate $E[p_{ik}^{g}(1) - p_{ik}^{g}(0)]$, where $p_{ik}^{g}(a)$ is the mouse-level proportion assigned to cluster $k$ by a rule not trained on that mouse-sample's cells.

### Connection to Existing Data
The current MG-selected Seurat clusters provide the reference labeling system. The same six Mouse × Condition samples can support leave-one-sample-out or leave-one-mouse-out assignment checks, with special care for paired mice 10 and 3 and single-condition mice 30 and 33.

### Approach
- Freeze a reference feature set and clustering recipe before testing.
- Build cluster-assignment rules in a cross-fit way:
  - Leave-one-Mouse × Condition sample out, learn or map clusters on the other five samples, then assign held-out cells to clusters without using their condition labels for cluster discovery.
  - For paired mice 10 and 3, also try leave-one-mouse-out assignment so neither condition from the same mouse defines the held-out mouse's clusters.
- Compute held-out mouse-level cluster proportions from assigned labels, then estimate E-Stim contrasts using the same Mouse × Condition design as Idea 1.
- Add sensitivity analyses that vary the cluster-definition rule:
  - Current clusters as the reference rule.
  - Condition-blind clusters learned after balancing cells per mouse-condition.
  - Cell-cycle-filtered branch as a complementary rule, not the primary rule.
- Separate two estimands in reporting:
  - Measurement-rule estimand: E-Stim effect on membership in the current frozen clusters.
  - Robust cell-state estimand: effects that remain directionally consistent under cross-fit and alternative condition-blind cluster definitions.

### Expected Output
A robustness matrix with clusters as rows and cluster-definition rules as columns, showing mouse-level E-Stim effect estimates, uncertainty, and whether the direction persists under cross-fit assignment and alternative definitions.

### Feasibility
Moderate. It requires more pipeline work than the near-term design, especially if reclustering or label transfer must be scripted carefully, but it remains feasible with current Seurat outputs and R tooling.

### Risks/Limitations
- Threats to validity include fold-specific cluster drift, weak replication after cross-fitting, sensitivity-rule multiplicity, and possible loss of true E-Stim-specific states.
- With only six Mouse × Condition samples, cross-fitting protects against circular cluster selection but does not create more biological replication.
- Label transfer can blur rare or borderline clusters and may change the meaning of cluster IDs across folds.
- If a cluster exists only because of E-Stim, condition-blind cross-fit rules may split or miss it; that result is scientifically important but hard to summarize as one stable cluster proportion.
- Multiple sensitivity rules can invite selective interpretation unless the primary rule and success criteria are pre-specified.

### Next Steps
- Decide whether leave-one-sample-out or leave-one-mouse-out is the primary circularity safeguard.
- Define stable cluster matching criteria across folds.
- Pre-specify that only effects consistent across the frozen-rule and cross-fit-rule analyses support robust causal language.
- Use disagreements between rules to label findings as cluster-definition-sensitive rather than failed.

## Recommendation
Prefer Idea 1 for the immediate ESPI revision because it aligns with the Mouse × Condition inferential unit, respects partial pairing, and can clearly replace cell-level pseudoreplication with sample-level effect estimates. Use Idea 2 as the stronger validation plan when making any claim that E-Stim changes a reproducible cell-state composition rather than membership in the current data-derived clusters.
