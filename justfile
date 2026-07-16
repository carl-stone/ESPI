set positional-arguments

# List available recipes
_default:
    just --list

# Run publication figures, marker analysis, DE/enrichment, and render the notebook
[group: "Publication analysis"]
run overwrite="false":
    ESPI_OVERWRITE={{overwrite}} Rscript scripts/02-publication-figures.R
    ESPI_OVERWRITE={{overwrite}} Rscript scripts/03-marker-analysis.R
    ESPI_OVERWRITE={{overwrite}} Rscript scripts/04-de-enrichment.R
    quarto render notebook/sc_analysis.qmd

# Generate publication figures and descriptive tables
[group: "Publication analysis"]
figures overwrite="false":
    ESPI_OVERWRITE={{overwrite}} Rscript scripts/02-publication-figures.R

# Run fixed MG-selected marker analysis
[group: "Publication analysis"]
markers overwrite="false":
    ESPI_OVERWRITE={{overwrite}} Rscript scripts/03-marker-analysis.R

# Run fixed MG-selected pseudobulk DE and enrichment
[group: "Publication analysis"]
de overwrite="false":
    ESPI_OVERWRITE={{overwrite}} Rscript scripts/04-de-enrichment.R

# Intentionally rebuild frozen objects after restoring directory write access
[group: "Expert and maintenance"]
regenerate-frozen:
    ESPI_OVERWRITE=true Rscript scripts/01-regenerate-frozen.R

# Load package code
[group: "Expert and maintenance"]
load:
    Rscript -e 'devtools::load_all(".", export_all = FALSE, quiet = TRUE)'

# Update package documentation
[group: "Expert and maintenance"]
document:
    Rscript -e 'devtools::document()'

# Rebuild README.md from README.Rmd
[group: "Expert and maintenance"]
readme:
    Rscript -e 'devtools::build_readme()'

# Format first-party R code with Air
[group: "Expert and maintenance"]
format:
    air format R scripts

# Run scilintr over first-party analysis code
[group: "Expert and maintenance"]
lint:
    #!/usr/bin/env bash
    set -euo pipefail
    for scope in R scripts data-raw notebook/sc_analysis.qmd config.local.example.R; do
        Rscript -e 'scope <- commandArgs(trailingOnly = TRUE)[1]; quit(status = scilintr::main(scope))' "$scope"
    done
