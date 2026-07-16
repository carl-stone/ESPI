# Cluster Proportion Testing Ideas

## Purpose

Generate expert analysis ideas for testing whether p27CKO + E-Stim changes MG-selected cluster proportions without treating cells as independent biological replicates.

## Status

active

## Datasets used

- `data/DATA_MANIFEST.md`: ESPI Seurat objects and generated MG-selected
  outputs under Box Drive.
- `scripts/01-regenerate-frozen.R`, `scripts/02-publication-figures.R`,
  `scripts/03-marker-analysis.R`, and `scripts/04-de-enrichment.R` are the
  four active phase scripts.
- `R/config.R`, `R/seurat-methods.R`, `R/publication-analysis.R`, and
  `R/publication-plots.R` are the four focused active modules. `R/config.R`
  owns the fixed paths and contracts used here.

## Algorithms used

- `R/publication-analysis.R::compute_cluster_abundance()` computes pooled
  cluster × condition Fisher/CLR summaries for descriptive plotting.
- `R/publication-analysis.R::compute_sample_cluster_proportions()` builds the
  Mouse × Condition sample-level proportions.
- `R/publication-analysis.R::test_cluster_proportion_randomization()` performs
  the deterministic exact sample-level randomization analysis.
- The four phase scripts and four focused modules are the active
  implementation; this page's persona proposals remain ideation, not runtime
  code.

## Key constraints

- Primary condition-level inference must use Mouse × Condition samples, not cells.
- Paired mice 10 and 3 should inform paired comparisons where methods support partial pairing.
- Cell-pooled Fisher tests can remain descriptive but should not support inferential treatment-effect claims.
- Cluster-definition circularity must be addressed explicitly.

## Outputs

- `00_index.md`: compiled idea index.
- Persona files: two ideas per selected expert persona.
