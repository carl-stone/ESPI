# Idea Generation Session: Clustering Criteria Brainstorm

**Date**: 2026-07-03
**Project**: ESPI
**Personas used**: 6
**Ideas per persona**: 2
**Total ideas generated**: 12
**Focus**: Label-blind criteria for selecting normalization, PC count, and Leiden clustering resolution.

## Criteria synthesis

The personas converge on the same decision shape:

1. **Choose from a stability plateau, not a single best-looking row.** Treat normalization, PC count, CC-HVG policy, and resolution as perturbations of the analysis pipeline.
2. **Use resolution 0.3 as the default candidate unless finer splits pass a predeclared gate.** Current grid evidence shows 0.3 has 10-11 clusters, low small-cluster burden, and higher agreement; 0.5/0.8 add 12-17 clusters and more fragmentation.
3. **Prefer PFlog if branch-invariance holds after pairwise checks.** Existing reference-relative summaries favor PFlog over log1p, but all-pairs or local-neighborhood agreement should reduce reference bias.
4. **Adjudicate 20 vs 30 vs 50 PCs by local neighborhood stability.** The key question is whether the 11-cluster 30/50-PC plateau adds a reproducible broad state beyond the 10-cluster 20-PC result, not whether more PCs look richer.
5. **Use marker/coherence only as a label-blind audit or predeclared tie-breaker.** Do not use condition separation, E-Stim effects, BrdU response, or desired biology to pick parameters.

## Ideas at a Glance

| # | Persona | Idea Title | Effort | Data Ready | One-Line Summary |
|---|---------|------------|--------|------------|------------------|
| 1 | Statistical Physicist | [Phase-Diagram Plateau Criterion](01_statistical-physicist.md#phase-diagram-plateau-criterion) | Medium | Yes | Choose the point with low response to perturbations in PCs and CC-HVG policy, low fragmentation, and broad agreement with nearby grid points. |
| 2 | Statistical Physicist | [Metastable Attractor and Split-Entropy Criterion](01_statistical-physicist.md#metastable-attractor-and-split-entropy-criterion) | Medium | Mostly | Classify higher-resolution splits as deep basins versus shallow substates using split entropy, size, cross-branch persistence, and later marker coherence. |
| 3 | Information Theorist | [Minimum Description Length Score for Grid Partitions](02_information-theorist.md#minimum-description-length-score-for-grid-partitions) | Medium | Mostly | Penalize unstable fragmentation and extra label bits unless added clusters preserve reproducible structure. |
| 4 | Information Theorist | [Variation-of-Information Stability Map Across Neighboring Grid Settings](02_information-theorist.md#variation-of-information-stability-map-across-neighboring-grid-settings) | Medium | Mostly | Build a stability map using VI/NMI across neighboring grid settings instead of ranking only by one reference clustering. |
| 5 | Topologist / Geometer | [Persistence Plateau for Resolution and PC Selection](03_topologist-geometer.md#persistence-plateau-for-resolution-and-pc-selection) | Medium | Mostly | Choose parameters whose clusters behave like persistent connected components across the ESPI parameter filtration. |
| 6 | Topologist / Geometer | [Mapper-Style Cover Consistency Across Normalization Branches](03_topologist-geometer.md#mapper-style-cover-consistency-across-normalization-branches) | Low to Medium | Yes | Select the normalization that behaves most like a stable coordinate chart across CC-HVG policy and PC count. |
| 7 | Stem Cell Biologist | [Fate-State Persistence Score for Choosing the Coarsest Stable Manifold](04_stem-cell-biologist.md#fate-state-persistence-score-for-choosing-the-coarsest-stable-manifold) | Medium | Mostly | Maximize broad fate-basin persistence while minimizing fragile subclustering. |
| 8 | Stem Cell Biologist | [Lineage-Priming Coherence Gate for Accepting Finer Resolution](04_stem-cell-biologist.md#lineage-priming-coherence-gate-for-accepting-finer-resolution) | Medium | Mostly | Separate primary coarse clustering from provisional higher-resolution substates accepted only after stable split and marker-coherence checks. |
| 9 | Representation Learning Specialist | [Invariance-First Latent Partition Selection](05_representation-learning-specialist.md#invariance-first-latent-partition-selection) | Medium | Mostly | Select parameters whose clusters are invariant to reasonable label-blind perturbations of the latent representation. |
| 10 | Representation Learning Specialist | [Factor-Basis Cluster Coherence Audit](05_representation-learning-specialist.md#factor-basis-cluster-coherence-audit) | Medium | Mostly | Test whether extra clusters from higher resolution or higher PC count behave like meaningful latent factors or fragmented graph partitions. |
| 11 | Causal Inference Researcher | [Sensitivity Frontier for Analysis Interventions](06_causal-inference-researcher.md#sensitivity-frontier-for-analysis-interventions) | Medium | Mostly | Treat parameters as interventions and pick the coarsest partition stable to adjacent analysis interventions. |
| 12 | Causal Inference Researcher | [Negative-Control Nuisance Dependence Screen](06_causal-inference-researcher.md#negative-control-nuisance-dependence-screen) | Medium | Mostly | Reject or flag clusterings whose extra clusters are strongly predicted by QC, cell-cycle, or other nuisance covariates. |

## By Feasibility

### Low to Medium Effort

- [Mapper-Style Cover Consistency Across Normalization Branches](03_topologist-geometer.md#mapper-style-cover-consistency-across-normalization-branches) — Uses existing grid summaries and clustree panels to compare normalization branch consistency.

### Medium Effort, Data Ready

- [Phase-Diagram Plateau Criterion](01_statistical-physicist.md#phase-diagram-plateau-criterion) — Uses existing 36-row grid metrics to score plateau stability.

### Medium Effort, Mostly Data Ready

- [Metastable Attractor and Split-Entropy Criterion](01_statistical-physicist.md#metastable-attractor-and-split-entropy-criterion)
- [Minimum Description Length Score for Grid Partitions](02_information-theorist.md#minimum-description-length-score-for-grid-partitions)
- [Variation-of-Information Stability Map Across Neighboring Grid Settings](02_information-theorist.md#variation-of-information-stability-map-across-neighboring-grid-settings)
- [Persistence Plateau for Resolution and PC Selection](03_topologist-geometer.md#persistence-plateau-for-resolution-and-pc-selection)
- [Fate-State Persistence Score for Choosing the Coarsest Stable Manifold](04_stem-cell-biologist.md#fate-state-persistence-score-for-choosing-the-coarsest-stable-manifold)
- [Lineage-Priming Coherence Gate for Accepting Finer Resolution](04_stem-cell-biologist.md#lineage-priming-coherence-gate-for-accepting-finer-resolution)
- [Invariance-First Latent Partition Selection](05_representation-learning-specialist.md#invariance-first-latent-partition-selection)
- [Factor-Basis Cluster Coherence Audit](05_representation-learning-specialist.md#factor-basis-cluster-coherence-audit)
- [Sensitivity Frontier for Analysis Interventions](06_causal-inference-researcher.md#sensitivity-frontier-for-analysis-interventions)
- [Negative-Control Nuisance Dependence Screen](06_causal-inference-researcher.md#negative-control-nuisance-dependence-screen)

## Recommended next implementation path

1. **Immediate, table-first criterion**: Implement a branch-invariance / plateau score from existing `cluster_grid_summary.tsv`.
2. **Reduce reference bias**: Extend the summary to compute all-pairs or local-neighborhood ARI/Jaccard/VI across all 36 clusterings.
3. **Resolution gate**: Use clustree-derived parent-child cross-tabs to quantify split entropy and child-size stability.
4. **Nuisance screen**: Check whether finer splits are predicted by QC/cell-cycle covariates before accepting them.
5. **Marker coherence audit**: After choosing the label-blind primary clustering or shortlist, test whether clusters and accepted splits have coherent marker programs.

## Files

- [CONTEXT.md](CONTEXT.md)
- [01_statistical-physicist.md](01_statistical-physicist.md)
- [02_information-theorist.md](02_information-theorist.md)
- [03_topologist-geometer.md](03_topologist-geometer.md)
- [04_stem-cell-biologist.md](04_stem-cell-biologist.md)
- [05_representation-learning-specialist.md](05_representation-learning-specialist.md)
- [06_causal-inference-researcher.md](06_causal-inference-researcher.md)
