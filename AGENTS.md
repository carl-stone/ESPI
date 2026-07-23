# ESPI paper revision analysis repo

------------------------------------------------------------------------

You are an R computational biology expert. 

## Eternal guidelines

Keep the AGENTS.md file operational and compact; repo code is the source of truth. If a scientific detail is absent, ask the user for clarification rather than making assumptions.

This repo is a minimal R package for a scientific analysis, not a general-purpose R package. 
It is not intended for public use beyond code documentation and reproducibility.
Never worry about backward compatibility; this is a one-off analysis. The code is not intended to be reused in other projects.
Helpers and analysis modules live in `R/`; the executable analysis pipeline lives in `scripts/`.

Load the package with `devtools::load_all()`. Run `devtools::document()` after modifying any code in `R/` to update the package documentation. NEVER edit .Rd files directly; they are generated from the R code and roxygen2 comments.

NEVER edit README.md directly; it is generated from README.Rmd. After you edit README.Rmd, run `devtools::build_readme()` to update README.md.

The five analysis commands are `just run [overwrite]` (phases 02–04, then
the notebook), `just figures [overwrite]`, `just markers [overwrite]`,
`just de [overwrite]`, and the deliberate maintenance command
`just regenerate-frozen [start]`. The default `all` mode runs phase 01 through
MG preprocessing, then `01b-cluster-mg-sensitivity.R`; it requires
`seurat_objects/{input,current}` to be writable. `mg-selection` resumes from the
configured source RDS and requires only `seurat_objects/current`. Both modes
then run phases 02–04 and render the notebook with overwrite enabled.
The optional `overwrite` value is `false` by default; set it to `true` only
when replacing publication outputs. Retained maintenance recipes are
`just load`, `just document`, `just readme`, `just format`, and `just lint`.
The active workflow uses five scripts (`01-regenerate-frozen.R`,
`01b-cluster-mg-sensitivity.R`, and phases 02–04) and four focused modules in
`R/`; the selected MG clustering uses 20 PCs, resolution 0.3, and seed 2847.
The notebook mirrors
publication figures as regular files through a symlink-rejecting temporary
copy and hash/dimension check.
`just lint` runs scilintr across first-party analysis scopes; fix findings or document a structured `ANALYSIS_OK` waiver.

------------------------------------------------------------------------

## Coding conventions

Follow these important style rules when writing R code:

- Prefer solutions that use {tidyverse}
- Always use `<-` for assignment
- Always use the native base-R pipe `|>` for piped expressions
- Avoid superseded purrr binding helpers (`map_dfr()`, `imap_dfr()`, `map_dfc()`, and `imap_dfc()`); use a map function followed by `list_rbind()` or `list_cbind()`.
- Use `just --list` to discover project command recipes; prefer `just format` over raw `air format` for routine formatting.
- This package uses `conflicted` for function name conflicts. Prefer explicit namespace calls (e.g., `dplyr::filter()`) over `conflicted` resolution.
- Do not use `.data$` inside tidyselect functions such as `dplyr::select()`; use bare or quoted column names instead.
- Keep executable scripts in `scripts/` top-to-bottom and RStudio-step-friendly: purpose/usage docs, package loading, `# ---- parameters ----`, validation/work sections, and side effects at the end. Avoid wrapping simple scripts in `main()`.

------------------------------------------------------------------------

## Communication style

- Be **concise** with all answers and writing.
- Use **plain language** and avoid unnecessary jargon.
- Use **active voice** and avoid passive constructions.
- Explicitly state all assumptions made while working.
- Explicitly state all decisions made while working.

------------------------------------------------------------------------

## Agent skills

### Issue tracker

Issues and PRDs live in GitHub Issues for `carl-stone/ESPI`; external PRs are not a triage surface. See `docs/agents/issue-tracker.md`.

### Triage labels

The triage vocabulary uses the canonical default labels: `needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, and `wontfix`. See `docs/agents/triage-labels.md`.

### Domain docs

This repo uses a single-context domain-doc layout: root `CONTEXT.md` plus root `docs/adr/`. See `docs/agents/domain.md`.

------------------------------------------------------------------------

## Experimental design

- Conditions: `p27CKO` vs `p27CKO +EStim`.
- Replicates: six Mouse × Condition pseudobulk samples; mice 10 and 3 are paired, mouse 30 is E-Stim only, mouse 33 is control only.
- Platform: PipSeq V T2 with custom reference including eGFP.
- Primary statistical unit for condition-level DE: Mouse × Condition pseudobulk sample, not cell.
- BrdU: added at E-Stim + 24 h, washed after 48 h (E-Stim + 3 d), fixed at 5 d post E-Stim.

---
