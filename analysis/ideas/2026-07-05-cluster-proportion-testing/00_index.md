# Cluster Proportion Testing Ideas

## Context

The panel evaluated how to test whether p27CKO + E-Stim changes MG-selected cluster proportions without treating cells as independent replicates. All ideas use Mouse × Condition samples as the inferential unit and treat the current pooled Fisher/CLR result as descriptive only.

## Recommended synthesis

Start with a mouse-level cluster count table and report sample-level effect sizes before p-values. The strongest near-term path is a transparent, partially paired analysis:

1. Build cluster counts/proportions for each Mouse × Condition sample.
2. Use paired mice 10 and 3 as the primary treatment contrast.
3. Use mouse 30 E-Stim-only and mouse 33 control-only as explicit unpaired context or sensitivity, not as substitutes for missing pairs.
4. Test or model transformed mouse-level proportions/log-ratios, with effect-size tables and paired-line plots.
5. Label cluster-inference claims as conditional on the frozen MG-selected clustering; use cross-fit or validation designs before making stronger claims about reproducible treatment-induced cell states.

## Idea table

| Persona | Practical idea | Ambitious idea | Feasibility | One-line summary |
|---|---|---|---|---|
| [Biostatistician](01-biostatistician.md) | Partially paired randomization tests on Mouse × Condition cluster proportions | Bayesian logistic-normal multinomial model with mouse-level partial pooling | High / Moderate | Use design-restricted sample-level randomization now; use hierarchical composition modeling if stronger uncertainty estimates are needed. |
| [Causal inference researcher](02-causal-inference-researcher.md) | Design-based cluster proportion contrasts using Mouse × Condition samples | Cross-fit cluster definition and sensitivity analysis | High / Moderate | Define the mouse-level treatment estimand, freeze clusters before testing, and use cross-fit labels to probe post-selection circularity. |
| [Compositional data analyst](03-compositional-data-analyst.md) | Sample-level ALR/CLR linear models | Dirichlet-multinomial or logistic-normal compositional model | High / Moderate | Analyze cluster count vectors as compositions, not independent raw proportions; run zero and reference sensitivity checks. |
| [Single-cell methods statistician](04-single-cell-methods-statistician.md) | Mouse-level pseudobulk cluster differential abundance with `speckle::propeller()` and edgeR sensitivity | Graph-neighborhood DA with `miloR` and cross-cluster sensitivity | High / Medium | Use pseudobulk cluster counts as the primary cluster test; use neighborhood DA to reduce dependence on one clustering resolution. |
| [Experimental design statistician](05-experimental-design-statistician.md) | Mouse-blocked log-ratio contrasts using paired mice as the anchor | Balanced paired validation with prespecified cluster definitions | High / Moderate | Make the current claim exploratory and mouse-level; use a balanced paired validation study for definitive cluster-proportion inference. |

## By feasibility

### Low effort / near-term

- Mouse-level paired contrasts with unpaired context.
- Design-restricted randomization over paired mice, with explicit p-value coarseness.
- Sample-level ALR/CLR models with pseudocount/reference sensitivity.
- `speckle::propeller()` or edgeR count sensitivity on cluster × sample counts.

### Medium effort / robustness

- Cross-fit or leave-one-mouse-out cluster assignment before abundance testing.
- `miloR` neighborhood differential abundance as a resolution-sensitivity layer.
- Dirichlet-multinomial or logistic-normal joint composition models.

### High effort / validation

- New balanced paired validation experiment with frozen reference labels or prespecified marker gates.
- Simulation-based power planning using current mouse-level variance and cluster baseline frequencies.

## Main risks to carry forward

- Six Mouse × Condition samples give low power; paired-only inference has two within-mouse contrasts.
- Mouse 30 and mouse 33 cannot replace missing paired samples.
- Cluster proportions are compositional; an apparent decrease can reflect expansion elsewhere.
- Rare clusters and zeros can dominate transformed-scale effects.
- Testing data-derived clusters is conditional on the chosen clustering unless cluster definitions are frozen, cross-fit, externally mapped, or validated.

## Files

- [01-biostatistician.md](01-biostatistician.md)
- [02-causal-inference-researcher.md](02-causal-inference-researcher.md)
- [03-compositional-data-analyst.md](03-compositional-data-analyst.md)
- [04-single-cell-methods-statistician.md](04-single-cell-methods-statistician.md)
- [05-experimental-design-statistician.md](05-experimental-design-statistician.md)
