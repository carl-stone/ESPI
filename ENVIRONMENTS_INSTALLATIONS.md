# Environments & Installations

## Primary environment

- **Project type:** Minimal R package for ESPI single-cell analysis.
- **Command runner:** `just` from the repository root.
- **Package load:** `just load` (`devtools::load_all()`).
- **Notebook renderer:** `quarto render notebook/sc_analysis.qmd`.
- **External data root:** resolved by `R/config.R`; the default is `/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data`.

## Setup from scratch

Install R and the packages declared in `DESCRIPTION`, Quarto, Air, and
`just` (for example, `brew install just`). From the repository root:

```sh
just load
```

Use `just document` after changing package code, `just readme` after changing
`README.Rmd`, `just format` for first-party R code, and `just lint` for
scilintr checks.

## Five analysis commands

The fixed publication interface has one optional run-level overwrite value
(`false` by default):

```sh
just run [overwrite]       # phases 02 → 03 → 04, then Quarto
just figures [overwrite]   # phase 02
just markers [overwrite]   # phase 03
just de [overwrite]        # phase 04
just regenerate-frozen     # phase 01; deliberate maintenance only
```

`just regenerate-frozen` requires
`seurat_objects/{input,current}` to exist and be writable; it refuses to run
when either directory is read-only. It has no source selector or other
scientific arguments. The routine publication phases load the fixed frozen
objects and protect existing outputs unless `overwrite` is `true`.

The three downstream inputs are loaded and validated by phase 02 as follows:

- final source:
  `current/cluster_pflog_no_filter_cc_elbow20.rds`,
  `cluster_pflog_no_filter_cc_dims30_res0.3`, 4,146 cells;
- final MG-selected:
  `current/cluster_pflog_mg_selected_no_filter_cc_elbow20.rds`,
  `cluster_pflog_mg_selected_no_filter_cc_dims20_res0.5`, 3,456 cells;
- CC-filtered MG sensitivity:
  `current/cluster_pflog_mg_selected_filter_cc_elbow20.rds`,
  `cluster_pflog_mg_selected_filter_cc_dims20_res0.5`, 3,456 cells.

Phase 03 marker tables are descriptive outputs only. Phase 04 independently
loads the final MG object and rebuilds curated marker overlap from package
marker data plus `Cdkn1b`; phase 03 output is not a phase-04 input.

## Expert maintenance

Retained maintenance recipes are:

```sh
just load
just document
just readme
just format
just lint
```

Run a phase directly only when diagnosing a fixed-output failure:

```sh
Rscript scripts/02-publication-figures.R
Rscript scripts/03-marker-analysis.R
Rscript scripts/04-de-enrichment.R
```

The active executable surface is exactly
`scripts/01-regenerate-frozen.R`, `scripts/02-publication-figures.R`,
`scripts/03-marker-analysis.R`, and `scripts/04-de-enrichment.R`.
The package surface is exactly `R/config.R`, `R/seurat-methods.R`,
`R/publication-analysis.R`, and `R/publication-plots.R`.

## Notebook and output safety

Notebook figure inputs are regular files, not symlinks. Each phase mirrors a
figure by refusing a symlink destination, copying to a temporary regular
sibling, checking SHA-256 and image dimensions, then atomically replacing the
existing regular destination. Never write through a symlink or copy directly
over a destination.

## Required R packages

Runtime dependencies are declared in `DESCRIPTION`, including Seurat,
SeuratObject, DropletUtils, scDblFinder, BiocParallel, DESeq2, apeglm,
clusterProfiler, enrichplot, enrichit, ComplexHeatmap, circlize, clustree,
scclrR, magick, digest, dplyr, ggplot2, gghalves, ggrepel, here, Matrix,
mclust, org.Mm.eg.db, patchwork, readr, scales, S4Vectors, tibble, and tidyr.
`biomaRt`, `devtools`, and `scilintr` are suggested packages.

## Data and output paths

`R/config.R` resolves the external data root in this order:

1. `MEGAN_SC_DATA_DIR` in an untracked `config.local.R` file;
2. `<BOX_PATH>/megan_sc_data` from `config.local.R`; or
3. the default `~/Library/CloudStorage/Box-Box/megan_sc_data`.

Copy `config.local.example.R` to `config.local.R` and set one of those values
when the default Box Drive location is unavailable. The resolved directory must
already exist; package loading fails explicitly when it does not.

Important subdirectories include `seurat_objects/current/`,
`figures/preprocess/`, `figures/cluster/`, `figures/mg_selected/`,
`tables/mg_selected/`, `degs/mg_selected/`, and
`enrichment/mg_selected/`.

## Mycelium local hooks

Mycelium created local Claude Code hooks in
`.claude/settings.local.json`. The file points at the installed plugin cache;
`.claude/` is gitignored because that path is machine-local. On a fresh clone
or after a plugin upgrade, rerun Mycelium initialization or hook setup before
expecting hooks to fire.
