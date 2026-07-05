# Experimental Design Statistician

## Idea 1: Mouse-blocked cluster proportion contrasts for the current partially paired design

### Title

Mouse-blocked log-ratio contrasts using paired mice as the anchor.

### Motivation

The current design has six Mouse × Condition samples, not many biological replicates. The safest near-term analysis should estimate treatment-associated cluster proportion shifts at the mouse level, preserve the two within-mouse comparisons, and avoid claiming that thousands of cells create thousands of replicates.

### Connection to Existing Data

Use one cluster count/proportion vector per Mouse × Condition sample: paired mice 10 and 3, E-Stim-only mouse 30, and control-only mouse 33. Treat the pooled Fisher/CLR plot as a descriptive screen only, not as evidence for a condition-level effect.

### Approach

- Build a sample-by-cluster count table for the MG-selected branch.
- Convert each sample to cluster proportions or log-ratios, using a small pseudocount to avoid undefined values.
- Primary estimate: for each cluster, compute the within-mouse E-Stim minus control log-ratio change for mice 10 and 3.
- Secondary context: compare the paired estimate against the unpaired contrast between mouse 30 E-Stim and mouse 33 control, and against leave-one-mouse sensitivity summaries.
- Report effect sizes and uncertainty with design-aware intervals from mouse-level bootstrap/permutation where possible, not cell-level resampling.
- Use blocking language: paired mice provide the cleanest treatment contrast; unpaired mice widen the external biological context but cannot resolve mouse-specific effects.
- Adjust or rank clusters cautiously because cluster definitions came from the same dataset. Treat this as confirmatory only for predeclared clusters and exploratory for clusters selected after seeing treatment structure.

### Expected Output

A per-cluster table and plot showing paired log-ratio treatment shifts, unpaired context points, direction agreement, and sensitivity to each mouse. Claims should read as mouse-level evidence of large, consistent shifts in this small design, not definitive population-level proof.

### Feasibility

High. The design uses existing cluster labels and Mouse × Condition metadata. It needs only sample-level aggregation and simple model/sensitivity summaries in R.

### Risks/Limitations

- Power is low: the paired effect estimate has only two within-mouse contrasts.
- Mouse 30 and mouse 33 are not interchangeable substitutes for missing paired samples.
- Cluster-level multiple testing will be weak; emphasize estimation and ranking over binary discovery.
- If clustering was influenced by condition-associated abundance or expression structure, abundance testing of those same clusters is partly circular.

### Next Steps

1. Define the exact Mouse × Condition cluster count table.
2. Predeclare which cluster labels are biologically interpretable enough for primary reporting.
3. Implement paired contrasts first, then add unpaired sensitivity displays.
4. Phrase current-data conclusions as exploratory or hypothesis-generating unless effects are large and robust across the paired mice.

## Idea 2: Future blocked validation study with prespecified cluster gates and balanced pairing

### Title

Balanced paired validation with prespecified cluster definitions.

### Motivation

A stronger cluster-proportion claim needs more independent mice, balanced treatment allocation, and cluster definitions that are not tuned to the same treatment effect being tested. The validation design should separate cluster discovery from inference.

### Connection to Existing Data

Use the current ESPI data to choose candidate clusters, estimate plausible effect sizes, and plan sample size. Do not use the same data as both the source of cluster hypotheses and the final confirmatory test.

### Approach

- Design a validation experiment with matched control and E-Stim samples for every mouse or litter/block where feasible.
- Treat mouse, litter, processing batch, and sequencing run as planned blocking factors.
- Prespecify cluster identities before testing: either map validation cells onto a frozen reference from current data or define marker-based cluster rules independent of validation treatment labels.
- Analyze each validation sample as a Mouse × Condition count vector with a paired or blocked model: cluster counts out of total cells per sample, with mouse/block effects and treatment effect.
- Power the study on mouse-level variance from the current data, using simulations over plausible intramouse correlation, total cell yield, and cluster baseline frequency.
- Include sensitivity targets: enough mice to detect large shifts in common clusters, and explicit acknowledgment that rare clusters may require many more mice or targeted enrichment.

### Expected Output

A design/power table giving proposed numbers of paired mice, expected detectable changes by baseline cluster frequency, and a confirmatory analysis plan. The final validation claim could support treatment-associated changes in prespecified cluster proportions across mice, with batch/block adjustment.

### Feasibility

Moderate. The statistical plan is straightforward, but it requires new biological samples and disciplined separation between discovery and validation. It becomes much stronger if all mice are paired and processed in balanced batches.

### Risks/Limitations

- A frozen reference can miss new treatment-induced states; marker gates can oversimplify continuous biology.
- Power depends more on mouse-to-mouse variability than cell count once each cluster has adequate cells.
- Batch imbalance could mimic treatment effects if control and E-Stim samples are not processed together.
- Validation can confirm prespecified cluster proportion changes, but it cannot retroactively make current-data cluster-selection claims non-circular.

### Next Steps

1. Use current mouse-level proportions to simulate power under several mouse counts and effect sizes.
2. Choose a small set of prespecified clusters for validation rather than testing every discovered cluster equally.
3. Plan paired control/E-Stim collection and balanced processing blocks before sequencing.
4. Register the sample-level model, sensitivity analyses, and claim language before looking at validation outcomes.

## Recommendation

Prefer Idea 1 for the current ESPI analysis: report mouse-blocked paired contrasts with unpaired samples as sensitivity/context, and limit claims to exploratory mouse-level evidence. Use Idea 2 as the validation path if the goal is a defensible treatment-effect claim about cluster proportions beyond this small partially paired dataset.
