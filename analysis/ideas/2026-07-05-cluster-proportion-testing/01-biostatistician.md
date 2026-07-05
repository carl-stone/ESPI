# Biostatistician ideas: cluster-proportion testing

## Idea 1: Design-restricted sample-level randomization tests

### Title
Partially paired randomization tests on Mouse × Condition cluster proportions

### Motivation
Use the experiment's real replication unit and avoid treating cells as independent animals. With only six Mouse × Condition samples, asymptotic model p-values will be fragile, so a design-restricted randomization test gives a transparent near-term inferential check.

### Connection to Existing Data
For each MG-selected cluster, use one count per Mouse × Condition sample: cells in the cluster and total MG-selected cells in that sample. The statistical unit is the Mouse × Condition sample, not individual cells. The two complete mice, 10 and 3, supply within-mouse condition contrasts; mouse 30 contributes E-Stim only and mouse 33 contributes control only.

### Approach
- Build a sample-level table with `mouse`, `condition`, `cluster`, `cluster_n`, `sample_total`, and `proportion`.
- For each cluster, compute a treatment effect on a stabilized sample-level scale, such as logit proportion with a 0.5 count correction: `logit((cluster_n + 0.5) / (sample_total + 1))`.
- Use a partially paired test statistic that prioritizes within-mouse contrasts: mean E-Stim minus control among paired mice, optionally combined with the singleton E-Stim-vs-control contrast from mouse 30 and mouse 33 as a sensitivity statistic.
- Generate the null distribution by preserving the design: flip condition labels within mice 10 and 3, and treat the unmatched E-Stim/control singleton assignment as one exchangeable unpaired block only if that assumption is accepted. Report the paired-only result and the paired-plus-singleton result separately.
- Correct across tested clusters with Benjamini-Hochberg FDR, while emphasizing effect sizes and the coarse p-value resolution from the small randomization space.

### Expected Output
A cluster-level table with effect estimate, paired-only statistic, optional paired-plus-singleton statistic, exact or permutation p-value, BH-adjusted q-value, and sample-level proportions. A companion plot can show each mouse's proportion by condition with paired mice connected and singleton mice shown separately.

### Feasibility
High. This is implementable with tidyverse data summaries and base R randomization enumeration. It needs no new model-fitting dependency and can sit beside the current pooled Fisher/CLR plot as the primary sample-level inferential screen.

### Risks/Limitations
The null distribution is very small, so p-values are coarse and power is low. The singleton comparison relies on cross-mouse exchangeability and should not override the paired-mouse result. Testing clusters defined from the same data can still be circular if condition structure helped create or select the clusters; frame results as conditional on the chosen clustering and avoid claiming discovery of treatment-created clusters without validation.

### Next Steps
1. Freeze the cluster set and state that the tests condition on those fixed labels.
2. Create the Mouse × Condition cluster-count table.
3. Implement paired-only and paired-plus-singleton randomization tests per cluster.
4. Apply BH FDR across clusters and review effect sizes before interpreting q-values.

## Idea 2: Hierarchical multinomial model for partially paired cluster composition

### Title
Bayesian logistic-normal multinomial model with mouse-level partial pooling

### Motivation
Cluster proportions are compositional: increasing one cluster necessarily decreases others. A hierarchical multinomial model can analyze the full cluster-count vector per Mouse × Condition sample, estimate treatment effects with uncertainty, borrow strength across clusters, and naturally handle partial pairing.

### Connection to Existing Data
Each Mouse × Condition sample contributes one vector of cluster counts and one sample total. The model treats mice 10 and 3 as paired through mouse-specific effects, while mice 30 and 33 still contribute information through the condition effect and the population-level mouse variance. Cells only contribute to sample-level counts; they are not independent replicates.

### Approach
- Model each sample's cluster-count vector as multinomial or Dirichlet-multinomial to allow overdispersion beyond simple multinomial sampling.
- Use a logistic-normal linear predictor with cluster-specific intercepts, cluster-specific E-Stim effects, and mouse-level random intercepts or deviations. The paired mice estimate within-mouse treatment shifts; the unpaired mice inform between-mouse variability and condition means without pretending to be pairs.
- Put weakly regularizing priors on treatment effects and hierarchical shrinkage across clusters to stabilize estimates under $n = 6$ samples.
- Summarize each cluster by posterior treatment log-odds ratio, posterior interval, posterior probability of increase/decrease, and a multiplicity-aware quantity such as posterior expected false discovery rate or local false sign rate.
- Add a design sensitivity analysis: fit the model with all six samples, paired mice only, and with wider mouse-variance priors to test whether conclusions depend on singleton assumptions.

### Expected Output
A posterior summary table for all clusters plus compositional treatment-effect plots. The strongest outputs would be calibrated uncertainty intervals and probability statements, for example `Pr(E-Stim increases cluster k proportion)`, rather than binary significant/non-significant calls.

### Feasibility
Moderate. The model is statistically appropriate but heavier than the near-term randomization analysis. It likely needs a Bayesian engine such as brms, rstanarm, Stan, or nimble, and careful prior/predictive checks. It is best used after the cluster labels and sample-level count table are finalized.

### Risks/Limitations
With six samples, priors and overdispersion assumptions will matter. Multinomial models can be sensitive to rare clusters and zero counts. The model still conditions on the observed clustering, so it does not remove circularity if clusters were chosen because they looked treatment-associated. A stronger version would define clusters on a treatment-blinded reference, external reference, or held-out data before testing proportions.

### Next Steps
1. Decide whether the main estimand is cluster-specific log-odds shift or absolute proportion shift.
2. Build the sample-by-cluster count matrix and metadata with mouse and condition.
3. Fit a simple hierarchical multinomial/Dirichlet-multinomial prototype and run posterior predictive checks.
4. Compare all-sample, paired-only, and prior-sensitivity results before making claims.

## Recommendation
Prefer Idea 1 for the immediate ESPI revision because it is transparent, sample-level, and honest about the tiny partially paired design. Use Idea 2 as the robust follow-up if the project needs model-based uncertainty for the full cluster composition or stronger multiplicity-aware probability statements.
