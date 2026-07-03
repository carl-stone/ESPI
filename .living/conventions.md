# Repo-Specific Conventions

Project-specific rules here override generic Mycelium convention packs.

## ESPI R Package Shape

- Keep ESPI as a minimal R package plus executable `scripts/` pipeline.
- Do not move existing `R/`, `scripts/`, or `notebook/` code into generic Mycelium `analysis/` or `algorithms/` directories unless the user explicitly asks.
- Keep helper functions only when they remove real repetition or clarify a deep operation. Do not create one-off helpers for a few commands.

## R and Documentation Workflow

- Use `<-` assignment and the native base pipe `|>`.
- Prefer explicit namespaces over conflict-resolution side effects.
- Run `air format <file>` after R edits.
- Run `devtools::document()` after editing files in `R/`.
- Never edit generated `.Rd` files directly.
- Edit `README.Rmd`; regenerate `README.md` with `devtools::build_readme()`.

## Data and Figures

- Box Drive is the data/output source of truth: `/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data`.
- Do not add fallback paths.
- Use notebook-relative figure paths in Quarto.
- Symlink Box figures into `notebook/figures/`; do not copy them unless explicitly requested.
- Re-render `notebook/sc_analysis.qmd` after regenerating figures because HTML embeds image bytes.

## Statistical Unit

The primary condition-level differential-expression unit is Mouse × Condition pseudobulk sample, not cell.
