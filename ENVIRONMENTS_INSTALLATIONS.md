# Environments & Installations

## Primary Environment

- **Project type**: Minimal R package for ESPI single-cell analysis.
- **Package load**: `devtools::load_all()` from the repo root.
- **Formatter**: Air (`air format <file>`).
- **Notebook renderer**: Quarto (`quarto render notebook/sc_analysis.qmd`).
- **External data root**: `/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data`.

## Setup From Scratch

```r
# From an R session with working directory at the repo root
devtools::load_all()
devtools::document()
```

```sh
# Format changed R files
air format R/<file>.R scripts/<script>.R

# Render the notebook after figure source changes
quarto render notebook/sc_analysis.qmd
```

## Tripwire Checks

Run the lightweight scientific-boundary checks from the repo root:

```sh
Rscript tools/run-tripwires.R
```

The runner uses `analysis_labels.yml` for label/contrast declarations and exits
non-zero only for `FAIL` rows. `SKIP` rows mark checks that need future pipeline
stages or scratch-output instrumentation.


## Required R Packages

Declared in `DESCRIPTION`:

- `clustree`
- `ggplot2`
- `ggraph`
- `here`
- `igraph`
- `scclrR` from `cleartools/scclrR`
- `Seurat`
- `SeuratObject`

Suggested:

- `biomaRt`
- `devtools`

## Data and Output Paths

The analysis expects Box Drive data at:

```text
/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data
```

Do not add fallback data paths. Missing Box paths should fail explicitly.

Important subdirectories:

```text
seurat_objects/current/
figures/preprocess/
figures/cluster/
```

Notebook figures should be symlinked into `notebook/figures/` and referenced with notebook-relative paths.

## Mycelium Local Hooks

Mycelium created local Claude Code hooks in `.claude/settings.local.json`. The file points at the installed plugin cache:

```text
/Users/carlstone/.omp/plugins/cache/plugins/mycelium___mycelium___0.0.0/
```

`.claude/` is gitignored because those paths are local and can break after plugin upgrades or cache cleanup. On a fresh clone or after upgrading the plugin, rerun Mycelium initialization or hook setup before expecting hooks to fire.
