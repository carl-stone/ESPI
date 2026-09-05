
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ESPI

This is the companion code to Stone et al. (2026), *Electrical
stimulation combined with p27Kip1 inactivation drives proliferative
neurogenic reprogramming of Muller glia in the adult mouse retina*.

## Getting started

This repository is a minimal R package plus manuscript and exploratory
analyses. Install R and the dependencies in `DESCRIPTION`, Quarto, and
`just`; configure the external data directory as described in
[setup](dev/setup.md), then load local package code:

``` sh
just load
```

To install the package from GitHub instead:

``` r
# install.packages("pak")
pak::pak("carl-stone/ESPI")
```

## Where things live

| Directory | Purpose |
|----|----|
| `R/` | Shared configuration, Seurat methods, statistics, and plotting functions. |
| `scripts/` | Manuscript analyses: the numbered pipeline plus reprogramming, custom markers, and QC summaries. |
| `notebook/` | Analysis report (`sc_analysis.qmd`), rendered HTML, and report/manuscript figures. |
| `exploratory/` | Milo, plotting experiments, and neurogenic-proportion sensitivity analysis. |
| `data-raw/` | Scripts that build the small package datasets. |
| `data/` | Package datasets; large inputs and primary pipeline outputs remain in Box. |
| `man/` | Generated R documentation. |
| `dev/` | Setup, study design, data provenance, manuscript drafts, and reference material. |

## Publication pipeline

Four phases span five scripts. Routine publication commands consume
saved Seurat objects; they do not recompute clustering:

``` sh
just run [overwrite]       # phases 02 → 03 → 04, then render the notebook
just figures [overwrite]   # publication figures and descriptive tables
just markers [overwrite]   # FindAllMarkers tables and dotplot
just de [overwrite]        # pseudobulk DE and enrichment
just regenerate-frozen [start] # deliberate regeneration, then downstream phases
```

`overwrite` is `false` by default. Set it to `true` only when replacing
existing publication outputs. `start` is `all` (default) or
`mg-selection`. Full regeneration requires writable input and
current-object directories; `mg-selection` resumes from the selected
source object and requires only the current-object directory. Both
regeneration modes run the MG clustering grid, then phases 02–04 and
notebook rendering with overwrite enabled.

| Step | Script |
|----|----|
| Counts, QC, source preprocessing, and MG selection | `scripts/01-regenerate-frozen.R` |
| MG clustering and sensitivity grid | `scripts/01b-cluster-mg-sensitivity.R` |
| Publication figures and descriptive tables | `scripts/02-publication-figures.R` |
| Cluster marker analysis | `scripts/03-marker-analysis.R` |
| Pseudobulk DE and enrichment | `scripts/04-de-enrichment.R` |

Phase 02 loads the selected source, MG, and cell-cycle-filtered
sensitivity objects. Phase 03 marker outputs are descriptive and do not
feed phase 04. Phase 04 independently loads the MG object and rebuilds
curated marker overlap from package marker data plus `Cdkn1b`.

**Change clustering settings in `publication_config()` in
`R/config.R`.** The selected MG clustering uses 20 PCs, resolution 0.5,
and seed 2847. Column names and plot labels derive from those settings.
The source stays at resolution 0.3; the cell-cycle-filtered sensitivity
uses resolution 0.5 and seed 1312.

Changing settings selects a column in a saved object; it does not
recompute that column. A seed change requires rerunning clustering. With
both MG preprocessing objects already available, rerun just that stage
with:

``` sh
ESPI_OVERWRITE=true Rscript scripts/01b-cluster-mg-sensitivity.R
```

Then update downstream outputs and the notebook deliberately. Existing
reports and drafting notes are snapshots, not evidence that a new
configuration has run.

The other package modules have distinct jobs: `R/seurat-methods.R` owns
PCA and cluster-grid summaries; `R/publication-analysis.R` owns
abundance, sample proportions, randomization, and module/p27 statistics;
`R/publication-plots.R` owns plot writers and notebook mirroring.

### Additional manuscript work

These scripts support the manuscript but are not automatically run by
the numbered pipeline. Keep interactive analyses section-by-section;
their presence here does not imply a shared execution order.

| Work | Script |
|----|----|
| Reprogramming scores and score decomposition | `scripts/reprog_scoring.R` |
| Marker-score violins and UMAPs | `scripts/plot-marker-scores.R` |
| Custom curated marker heatmaps and dot plot | `scripts/custom-marker-plots.R` |
| Interactive marker analysis with a cluster merge | `scripts/findmarkers.R` |
| Sample-level QC summary | `scripts/qc_table.R` |

Custom-marker figures live in `notebook/figures/custom-markers/`,
including retained manuscript candidates not yet embedded in the report.
The custom plotting script writes there. File writers use
`output_path()` to create parent directories and require overwrite
opt-in for existing files.

### Exploratory work

- [Milo](exploratory/milo/notes.md): neighborhood differential abundance
  in PCA and module-score space.
- [Plotting](exploratory/plotting/notes.md): interactive expression and
  detection exploration, plus preserved prior exports.
- [Neurogenic proportion](exploratory/neurogenic-proportion/notes.md):
  the completed threshold-grid sensitivity analysis, with saved results
  alongside its code. It consumes the module-score Milo table, not the
  PCA Milo results.

## Documentation

- [Setup and external paths](dev/setup.md)
- [Study design and terminology](dev/study.md)
- [Data provenance and output locations](dev/data.md)
- [Methods and Results drafting material](dev/methods-results.md)
- [Manuscript plan and interpretation notes](dev/manuscript-plan.md)
- [Reference figure PDF](dev/references/pnas_sc_figures.pdf)

Drafting notes retain their evidence annotations; they are not
automatically refreshed result reports.

## Editing and rendering

Use `just document` after changing `R/`, and `just readme` after
changing `README.Rmd`. `just format` formats first-party R code;
`just lint` is an optional on-demand diagnostic. No full analysis run is
required for a documentation edit.

Render the report with `quarto render notebook/sc_analysis.qmd` after
changing its prose or figure inputs when updating the HTML deliverable.
It embeds image bytes. Pipeline writers retain their primary PNG/PDF
outputs in Box and copy report figures into `notebook/figures/`. The
scripts choose which figures to copy; they do not parse the notebook,
hash images, or maintain rollback files. Grid stages copy the selected
clustering views rather than every candidate.

`output_path()` checks overwrite permission at each write and rejects
symlink destinations. There is no separate output inventory; if a later
write fails, earlier outputs from that run remain. Use
`ESPI_OVERWRITE=true` for deliberate replacement when running scripts
directly, or `output_path(..., overwrite = TRUE)` for an explicitly
approved write in an interactive session.
