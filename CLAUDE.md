# CLAUDE.md — ESPI Mycelium Setup

This repository is Mycelium-enabled. `AGENTS.md` remains the source of truth for ESPI scientific context, R style, and package workflow. Use this file for the Mycelium memory layer and Claude Code hook workflow.

## Start Here

1. Read `AGENTS.md` for project rules and experimental design.
2. Read `.living/INDEX.md` for accumulated Mycelium decisions, learnings, and conventions.
3. Read `CONTEXT.md` for the ESPI glossary before writing scientific prose or analysis labels.
4. Read `ENVIRONMENTS_INSTALLATIONS.md` before running the R package or Quarto notebook.
5. Read the relevant manifest when touching Mycelium-managed areas:
   - `analysis/ANALYSIS_MANIFEST.md`
   - `algorithms/ALGORITHM_MANIFEST.md`
   - `data/DATA_MANIFEST.md`
   - `reference_material/REFERENCE_MANIFEST.md`

## Repository Shape

ESPI is a minimal R package plus executable analysis scripts. Do not reorganize existing package code into generic Mycelium folders.

- `R/`: package helpers and plotting functions.
- `scripts/`: executable preprocessing and clustering pipeline.
- `notebook/`: Quarto analysis notebook and symlinked figure inputs.
- `data/`: R package data plus Mycelium metadata/raw/processed placeholders.
- `analysis/`: reserved for future standalone Mycelium analysis reports, not a replacement for `scripts/`.
- `algorithms/`: reserved for reusable method writeups; current code remains in `R/`.
- `.living/`: decisions, learnings, conventions, findings, and session logs.

## R Workflow

Use R-specific workflow, not Python/marimo defaults from generic Mycelium docs.

- Load package code with `devtools::load_all()`.
- Run executable pipeline steps with `Rscript scripts/<script>.R`.
- Run `devtools::document()` after modifying files in `R/`.
- Edit `README.Rmd`, then run `devtools::build_readme()`; never edit `README.md` directly.
- Format R code with `air format <file>`.
- Render the notebook with `quarto render notebook/sc_analysis.qmd` after figure source changes because `embed-resources: true` embeds image bytes.
- Run tripwire checks with `Rscript tools/run-tripwires.R` after changing analysis paths, report prose, or scientific-boundary code.
- Prefer tidyverse-style R where useful, but keep helpers narrow and avoid one-off helper functions for a few commands.
- Routine workflow rule: use `just run [source] [overwrite]` or `just run-dry-run [source] [overwrite]`; defaults are `counts-qc` and `false`, source also accepts `legacy` or a quoted explicit RDS path. Treat raw scripts and `just preprocess*` as expert recovery only; full runs render the notebook before tripwires.

## Installed Mycelium Convention Packs

Check `.living/conventions/ACTIVE_CONVENTIONS.yaml` for the registry.

- `robust-analysis`: validation, spot checks, sensitivity analysis, null controls.
- `report-generator`: structured reports and report QC.
- `idea-generator`: persona-based ideation.
- `bioinformatics`: single-cell and genomics analysis conventions.
- `skill-bridge`: routes analyses to inert external skill repositories under `skillpacks/`.

Apply these only where they fit the ESPI R package. Project-specific rules in `AGENTS.md` and `.living/conventions.md` override generic convention text.

## After Significant Work

For analysis, data, or method changes:

1. Update the relevant manifest if a tracked artifact changes.
2. Append `.living/decisions.md` for non-obvious choices.
3. Append `.living/learnings.md` for gotchas, bugs, or reusable workflow lessons.
4. Append `.living/findings/<topic>.md` only for scientific findings with evidence.
5. Update `.living/INDEX.md` with Mycelium's index generator.
6. Update `.claude/last-session.md` locally even if not using Claude Code.

Do not let Mycelium bookkeeping override the direct ESPI workflow. Keep entries concise and evidence-based.

## Local Hook Note

Mycelium generated `.claude/settings.local.json` with hooks pointing at the local plugin cache under `/Users/carlstone/.omp/plugins/cache/plugins/mycelium___mycelium___0.0.0/`. `.claude/` is gitignored because those paths are machine- and plugin-version-specific. Re-run Mycelium initialization or hook setup after a plugin upgrade, cache cleanup, or fresh clone.

