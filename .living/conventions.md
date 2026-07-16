# Repo-Specific Conventions

Project-specific rules here override generic Mycelium convention packs.

## ESPI R Package Shape

- Keep ESPI as a minimal R package plus executable `scripts/` pipeline.
- Do not move existing `R/`, `scripts/`, or `notebook/` code into generic Mycelium `analysis/` or `algorithms/` directories unless the user explicitly asks.
- Keep helper functions only when they remove real repetition, enforce an invariant, or clarify a deep operation. Do not create one-off helpers for a few commands, and do not wrap one- or two-line R expressions in helper functions; inline those at the call site.

## R and Documentation Workflow

- Use `<-` assignment and the native base pipe `|>`.
- Prefer explicit namespaces over conflict-resolution side effects.
- Use `just --list` to discover project commands. Prefer `just` recipes for routine package, pipeline, notebook, formatting, and tripwire steps; use raw `Rscript` only when a recipe does not fit.
- Describe `devtools::load_all()` / `just load` as loading the package, not as a smoke test; loading is a normal development workflow step.
- Run `devtools::document()` after editing files in `R/`.
- Never edit generated `.Rd` files directly.
- Edit `README.Rmd`; regenerate `README.md` with `devtools::build_readme()`.
- Keep executable scripts in `scripts/` top-to-bottom and RStudio-step-friendly: purpose/usage docs, package loading, `# ---- parameters ----`, validation/work sections, and side effects at the end. Avoid wrapping simple scripts in `main()`.

## Data and Figures

- Box Drive is the data/output source of truth: `/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data`.
- Do not add fallback paths.
- Use notebook-relative figure paths in Quarto.
- Symlink Box figures into `notebook/figures/`; do not copy them unless explicitly requested.
- Re-render `notebook/sc_analysis.qmd` after regenerating figures because HTML embeds image bytes.

## Statistical Unit

The primary condition-level differential-expression unit is Mouse × Condition pseudobulk sample, not cell.

## Cross-References

Refer to figures, tables, decisions, and other artifacts by stable handles — Quarto figure IDs (`#fig-...`), file basenames, or decision headings — not auto-numbered labels such as `Figure N` in notebook prose, session logs, review reports, and Mycelium records. Auto-numbers shift when panels are inserted or removed.

## Mycelium Session-End Provenance

Before yielding or committing after Mycelium-tracked work, verify semantic session-end records against the current `git status --short` output:

- `LOG_REGISTRY.md` current row has a sentence Summary, non-empty Key Outputs, useful Tags, and no filename-stub summary.
- The linked session log has `## Session Summary`, `## Key Outputs`, and `## Status` sections that describe the work and verification.
- `.living/INDEX.md` is regenerated after `.living/` decision, learning, convention, finding, or log-registry changes.
- `.claude/last-session.md` covers the full session and reports the current commit/uncommitted state.
- False-positive hook artifacts, such as `.log-scribe-*` authentication-failure logs or file-list-only session logs from read-only Q&A, are removed rather than backfilled into prior closed sessions.
- Generated activity inventories may include disposable-worktree paths and `xd:/` tool-device names; record their provenance, but never treat them as live-repository changes without confirmation from the main working tree.

Source: `.living/generated-conventions/session-end-provenance-integrity/` from L-10, L-20, L-21, L-22, and L-23.
