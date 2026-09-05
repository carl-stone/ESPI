# Setup and external paths

Run commands from the repository root. This is an R package and analysis repo;
no agent-specific tooling is needed to run it.

## Software

Install R and the packages declared in `DESCRIPTION`, Quarto, and `just`
(for example, `brew install just`). Air provides optional R formatting;
`scilintr` provides optional lint diagnostics.

Additional manuscript scripts use packages such as `ggview` and `ggplotify`.
The Bayesian reprogramming-model section also needs `cmdstanr` and a working
CmdStan installation; those are not needed for the numbered pipeline.
Exploratory Milo scripts use `miloR` and `SingleCellExperiment`; the
neurogenic-proportion analysis uses `glmmTMB`, and the plotting sandbox uses
`lme4`. Consult each script's package-loading section when running that work.
These scripts are not automatically executed by `just run`.

## External data directory

`R/config.R` resolves the data root in this order:

1. `MEGAN_SC_DATA_DIR` in an untracked `config.local.R` file;
2. `<BOX_PATH>/megan_sc_data` from `config.local.R`; or
3. `~/Library/CloudStorage/Box-Box/megan_sc_data`.

Copy `config.local.example.R` to `config.local.R` and set one of those values
when the default location is unavailable. Package loading does not require the
data directory; an analysis needs its actual input files. Do not commit
`config.local.R` or large external data.

Important subdirectories include `seurat_objects/{input,current}/`,
`figures/{preprocess,cluster,mg_selected}/`, `tables/mg_selected/`,
`degs/mg_selected/`, and `enrichment/mg_selected/`. The restructure does not move
or rewrite these external data or outputs. See [data and outputs](data.md).

## Load and run

```sh
just load                 # devtools::load_all()
just --list               # available recipes
just figures              # refuses to replace existing publication outputs
just figures true         # deliberately replaces existing publication outputs
```

The [README](../README.md) maps the numbered pipeline, manuscript scripts, and
exploratory analyses. `just run` executes phases 02–04 and renders the notebook;
it does not regenerate the input Seurat objects.

## Deliberate frozen regeneration

- `just regenerate-frozen` starts from raw counts and requires
  `seurat_objects/input/` and `seurat_objects/current/` to exist and be writable.
- `just regenerate-frozen mg-selection` resumes from the configured source RDS.
  It requires only `seurat_objects/current/` to be writable and skips count
  construction, QC, and source preprocessing/clustering.

Both modes run `01b-cluster-mg-sensitivity.R` and all downstream phases, then
render the notebook with overwrite enabled. These start modes select execution
scope, not a different scientific analysis.

## Changing clustering choices

Edit the `selected` list in `publication_config()` in `R/config.R`. Selected MG
uses 20 PCs, resolution 0.5, and seed 2847. Source clustering retains resolution
0.3; the cell-cycle-filtered MG sensitivity uses resolution 0.5 and seed 1312.
Column names derive from these settings. No fixed cell or cluster count blocks
a different clustering choice.

For a seed change, rerun `01b-cluster-mg-sensitivity.R` before downstream
analyses. That stage only loads the two MG preprocessing objects; it no longer
loads unrelated source branches to produce `frozen_object_numbers.tsv`.
The existing summary table, figures, and manuscript notes remain earlier-run
artifacts until deliberately replaced; they are not automatically relabeled.

## Maintenance

```sh
just document             # regenerate R documentation
just readme               # rebuild README.md from README.Rmd
just format               # format R/, scripts/, data-raw/, and exploratory/
just lint                 # optional scilintr diagnostics
quarto render notebook/sc_analysis.qmd
```

The notebook embeds figure bytes. Scripts copy chosen figures with ordinary
file copying, without notebook parsing, hashing, or backup/rollback machinery.
`output_path()` creates output directories and checks overwrite permission at
each write. A failed run may leave earlier completed outputs in place.
The custom-marker script also retains manuscript figures under
`notebook/figures/custom-markers/` independently of the numbered pipeline.
