# Exploratory plotting sandbox

`plot-sandbox.R` is an interactive, noncanonical workspace for last-mile manuscript plot exploration. It loads the configured final MG-selected Seurat object plus the current full and significant DESeq2 tables.

The script explores condition-stratified gene-expression violins, gene-pair scatterplots, cluster-focused views, and sample-level quasibinomial detection models for selected cluster 4 and 5 marker combinations.

Run sections interactively after `devtools::load_all()`. The file has no stable command-line interface or output contract and is not part of the four-phase publication pipeline. Promote any accepted analysis into `R/` or `scripts/`, add fixed output paths and validation, and regenerate the notebook before using it for publication results.
