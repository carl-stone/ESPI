# Single-cell Methods Statistician Ideas

## Idea 1: Mouse-level pseudobulk cluster differential abundance

### Motivation

Test whether E-Stim changes MG-selected cluster proportions using the same biological unit used for DE: one observation per Mouse × Condition sample. This converts the current pooled Fisher/CLR pattern into an inferential workflow that does not treat cells as independent replicates.

### Connection to Existing Data

Use the current MG-selected Seurat outputs, interpreted cluster labels, Mouse ID, condition, and per-cell cluster assignment. Build a cluster × sample count table for the six Mouse × Condition samples: mice 10 and 3 paired, mouse 30 E-Stim only, and mouse 33 control only.

### Approach

1. Count cells per cluster per Mouse × Condition sample and record each sample's total MG-selected cell count.
2. Primary practical model: run `speckle::propeller()` on cluster proportions with a design matrix such as `~ paired_block + condition`, where `paired_block` has separate levels for mice 10 and 3 plus an `unpaired` level for mice 30 and 33. This lets the paired mice contribute within-mouse comparisons while retaining the two unpaired samples.
3. Add a paired-only sensitivity using mice 10 and 3 with `~ mouse + condition`; treat agreement in sign and approximate magnitude as stronger evidence than either fit alone.
4. Add a count-based sensitivity using `edgeR::DGEList()` on cluster counts, sample library sizes equal to total MG-selected cells, robust dispersion estimation, and `edgeR::glmQLFit()` / `edgeR::glmQLFTest()` with the same partial-pairing design. Avoid interpreting this as gene-style expression; use it as a cluster-count uncertainty check.
5. Flag clusters with very small counts, absent counts in multiple samples, or unstable estimates. Report raw per-sample proportions beside model estimates.
6. Treat the existing pooled Fisher/CLR plot as descriptive only; do not cite it as primary evidence.

### Expected Output

A table with cluster, per-sample counts/proportions, E-Stim effect estimate, standard error or confidence interval, p-value, FDR, minimum sample count, and sensitivity labels: full partial-pair fit, paired-only fit, and edgeR count check. A compact plot would show each mouse-level proportion with paired lines for mice 10 and 3 and model-estimated direction.

### Feasibility

High. This is a straightforward R/tidyverse implementation from current cell metadata. `speckle` is purpose-built for cell-type or cluster proportion testing, and the edgeR sensitivity uses familiar Bioconductor GLM machinery.

### Risks/Limitations

Only six Mouse × Condition samples limit power, and the paired-only sensitivity has only two paired mice. The `unpaired` block is a pragmatic encoding, not a substitute for more balanced replication. Compositional coupling means an increase in one cluster can force apparent decreases elsewhere. Small clusters can have high binomial/count uncertainty. Cluster labels were defined from the same dataset, so this tests abundance of the chosen cluster partition rather than proving a condition-independent biological subtype changed.

### Next Steps

Prototype the cluster × sample count table, run `propeller` and edgeR sensitivity fits, predefine low-count filters, and compare full partial-pair results against paired-only signs before writing any biological claim.

## Idea 2: Graph-neighborhood DA with cross-cluster sensitivity and circularity checks

### Motivation

Reduce dependence on one selected clustering resolution and separate abundance testing from post hoc cluster interpretation. Neighborhood differential abundance tests can detect local cell-state shifts while still aggregating cells to Mouse × Condition counts, avoiding cell-level pseudoreplication.

### Connection to Existing Data

Use the current MG-selected single-cell object, UMAP/PCA or nearest-neighbor graph, Mouse × Condition metadata, and existing cluster annotations only as labels for interpreting significant neighborhoods after testing. Repeat on the primary MG-selected PFlog/no-cell-cycle branch and, if available, the cell-cycle-filtered branch as a sensitivity analysis.

### Approach

1. Build a `SingleCellExperiment` from current Seurat outputs and run `miloR` on the MG-selected graph or reduced dimensions.
2. Define graph neighborhoods without using cluster-level enrichment results as a selection rule. Count cells per neighborhood per Mouse × Condition sample.
3. Test neighborhood DA with `miloR::testNhoods()` using an edgeR-backed design such as `~ paired_block + condition`; run a paired-only sensitivity with `~ mouse + condition` for mice 10 and 3.
4. Map significant neighborhoods back to current clusters only after testing, reporting cluster overlap as annotation rather than as the tested hypothesis. This directly addresses circularity from clusters that may have been selected or named because they separate by condition.
5. Quantify cluster-definition uncertainty by repeating the workflow across nearby clustering resolutions, primary versus cell-cycle-filtered branches, and leave-one-mouse-out or downsampled-cell resamples. Summarize whether the same biological region/cluster family remains DA.
6. Quantify count uncertainty by tracking neighborhood size, number of contributing samples, and bootstrap stability of per-sample neighborhood proportions. Suppress claims for neighborhoods driven by one mouse or one tiny local region.

### Expected Output

A neighborhood-level DA table with spatial FDR, log-fold abundance effect, contributing sample counts, overlap with current clusters, branch/resolution stability, and leave-one-mouse-out sensitivity. A summary heatmap could show which current clusters contain stable DA neighborhoods without making the original clusters the inferential units.

### Feasibility

Medium. `miloR` is an established R package for single-cell differential abundance using sample-level neighborhood counts and edgeR models. The workflow needs more setup than pseudobulk cluster proportions but can reuse current Seurat embeddings, graph structure, and metadata.

### Risks/Limitations

With four mice and six Mouse × Condition samples, neighborhood DA may still be underpowered and sensitive to graph construction. Many neighborhoods create a multiple-testing burden. Cross-resolution and branch sensitivity can expose instability rather than solve it. Partial pairing remains imperfect: full models borrow information from unpaired samples, while paired-only models depend on two mice. Results should be framed as robust localization of abundance shifts, not proof that a preselected cluster is a condition-caused cell type.

### Next Steps

Create a small `miloR` prototype on the current MG-selected object, confirm that neighborhood sample counts match Mouse × Condition totals, run full and paired-only designs, then compare significant neighborhoods with current cluster labels and the practical pseudobulk results.

## Recommendation

Prefer Idea 1 as the near-term primary analysis because it directly answers the cluster-proportion question with Mouse × Condition replication and clear partial-pair sensitivities. Use Idea 2 as the robustness layer before making strong claims about any cluster whose definition may be condition-influenced or unstable across clustering choices.
