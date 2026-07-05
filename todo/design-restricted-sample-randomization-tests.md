# Design-restricted sample-level randomization tests

| Field | Value |
|-------|-------|
| **Date** | 2026-07-05 |
| **Author** | OMP agent |
| **Priority** | medium |
| **Status** | done |
| **Category** | analysis |
| **Related analyses** | [Cluster proportion testing ideas](../analysis/ideas/2026-07-05-cluster-proportion-testing/00_index.md) |
| **Related data** | ESPI Mouse × Condition MG-selected cluster counts |

## Description

Implement the biostatistician practical idea: design-restricted randomization tests on Mouse × Condition MG-selected cluster proportions.

## Motivation

The current pooled Fisher/CLR abundance plot is descriptive because it treats cells as pooled condition-level counts. A sample-level randomization test would use the experiment's biological replicate unit, preserve the paired mice, and provide a transparent primary inferential screen for cluster-proportion shifts.

## Proposed Approach

- Build a sample-level table with `mouse`, `condition`, `cluster`, `cluster_n`, `sample_total`, and `proportion`.
- For each cluster, compute a stabilized sample-level effect scale, such as `logit((cluster_n + 0.5) / (sample_total + 1))`.
- Use paired mice 10 and 3 as the primary contrast: mean E-Stim minus control within mouse.
- Optionally report a paired-plus-singleton sensitivity that includes mouse 30 E-Stim-only and mouse 33 control-only as an exchangeable unpaired block only if that assumption is accepted.
- Generate the null distribution by preserving the design: flip labels within paired mice and keep paired-only and paired-plus-singleton results separate.
- Apply Benjamini-Hochberg FDR across clusters, while emphasizing effect sizes and the coarse p-value resolution from the small design.

## Acceptance Criteria

- [x] Sample × cluster count table is generated from the MG-selected branch with one row per Mouse × Condition × cluster.
- [x] Paired-only randomization results are written with effect estimates, raw p-values, BH-adjusted q-values, and sample-level proportions.
- [x] Paired-plus-singleton sensitivity is either implemented and clearly labeled or explicitly deferred with rationale.
- [x] Output plots show mouse-level proportions, connect paired mice, and show singleton mice separately.
- [x] Notebook/report text states that inference conditions on the chosen cluster labels and does not use cell-pooled Fisher p-values as primary evidence.

## Notes

Source idea: `analysis/ideas/2026-07-05-cluster-proportion-testing/01-biostatistician.md`, Idea 1.

Completed 2026-07-05 in `R/cluster-abundance.R`, `scripts/plot-mg-selected-figures.R`,
and `notebook/sc_analysis.qmd`. Outputs are written under
`TABLE_DIR/mg_selected/` and `FIGURE_DIR/mg_selected/`.
