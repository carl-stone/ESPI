# Targeted Neurogenic Proportion

## Purpose

Estimate whether E-Stim samples contain a larger fraction of cells with high
neurogenic-progenitor scores and low proliferation scores, and test whether the
estimated direction is stable across prespecified score thresholds.

## Status

**Status**: complete exploratory analysis

The threshold grid was selected before running the grid bootstrap. This analysis
supports estimation and sensitivity assessment; it is not part of the canonical
four-phase publication pipeline.

## Datasets

- The frozen MG-selected Seurat object configured by `publication_config()`.
- `module_scores.tsv` from `scripts/module-score-milo-da.R` with primary settings
  `k = 60` and `prop = 0.04`.

## Methods

The six Mouse × Condition samples are the analysis units. Control-derived score
quantiles give each control sample equal total weight. The grid crosses
progenitor percentiles 80, 90, and 95 with proliferation percentiles 25 and 50.

Each gate uses beta-binomial regression with a condition-only logit model. A
likelihood-ratio statistic compares the condition model with an intercept-only
model. Each null distribution uses 1,000 parametric-bootstrap simulations.
Failed Hessian checks are retried with BFGS; remaining failures are retained in
the output so lower and upper p-value bounds can be audited. Bootstrap p-values
are adjusted across the six gates with the Benjamini-Hochberg method.

## Key Findings

- Every gate estimated a higher target-cell proportion under E-Stim.
- E-Stim/control mean-proportion ratios ranged from 1.8 to 5.2.
- Bootstrap calibration did not yield an adjusted value below 0.10.
- The sparsest gate had the most optimizer failures and should not be interpreted
  from its raw bootstrap p-value alone.

## Reproducibility

First generate the primary module-score table:

```bash
ESPI_MODULE_SCORE_MILO_K=60 \
ESPI_MODULE_SCORE_MILO_PROP=0.04 \
Rscript scripts/module-score-milo-da.R
```

Then run this analysis:

```bash
Rscript analysis/targeted-neurogenic-proportion/analysis.R
```

Set `ESPI_MODULE_SCORE_TABLE` to use a module-score table at another path. Set
`ESPI_OVERWRITE=true` only when deliberately replacing the saved outputs.

## Outputs

| File | Description |
|------|-------------|
| `outputs/threshold_grid_results.tsv` | Gate definitions, effects, bootstrap diagnostics, p-values, and BH-adjusted values. |
| `outputs/threshold_grid_sample_proportions.tsv` | Responder counts and proportions for every gate and biological sample. |
| `outputs/bootstrap_null_statistics.tsv` | All simulated likelihood-ratio statistics and optimizer status values. |
| `outputs/bootstrap_settings.tsv` | Input hash, thresholds, seeds, model specification, and package versions. |
| `outputs/threshold_grid_sample_proportions.png` | Faceted sample-level gate proportions with paired mice connected. |
| `outputs/threshold_grid_sample_proportions.pdf` | Vector version of the faceted sample-level plot. |
