# ESPI paper revision analysis repo

## Scope and workflow

This is a one-off scientific analysis in a minimal R package, not a reusable
software framework. Helpers live in `R/`; manuscript analyses live in `scripts/`;
Milo and other exploratory work live in `exploratory/`. Do not add
backward-compatibility layers or unnecessary abstractions.

- Edit the requested source files directly, run checks relevant to the change,
  and report what changed and what you actually verified.
- Code is the source of truth. Ask Carl when a necessary scientific detail is
  missing rather than inventing an assumption.
- Preserve scientific input checks, the Mouse × Condition statistical unit,
  and explicit permission to replace publication outputs.
- Do not add or restore agent hooks, edit guards, or automatic background agents.
  No session logs, registries, manifest updates, or convention generation are
  required to edit files or finish a task. Record scientific methods and
  non-obvious choices in the relevant code or analysis prose, not duplicate logs.
- Use `dev/study.md` when relevant to scientific terminology and `dev/setup.md`
  for environment setup. The README maps the workflow and manuscript scripts.

## Commands and generated files

Use `just --list` to discover recipes. The publication commands are:

- `just run [overwrite]`: phases 02–04, then render the notebook.
- `just figures [overwrite]`, `just markers [overwrite]`, `just de [overwrite]`:
  run one publication phase.
- `just regenerate-frozen [start]`: deliberate regeneration from `all` or
  `mg-selection`, followed by downstream phases and notebook rendering.
  See `dev/setup.md` for required writable directories.

`overwrite` defaults to `false`; use `true` only when replacing publication
outputs. The selected MG clustering uses 20 PCs, resolution 0.3, and seed 2847.
Edit `publication_config()` in `R/config.R` to change the selected settings;
column names derive from those settings. Do not add hard-coded cell/cluster counts.

- Load the package with `devtools::load_all()` (`just load`).
- After editing `R/`, run `devtools::document()` (`just document`). Edit roxygen
  comments, not generated `man/*.Rd` or `NAMESPACE`.
- Edit `README.Rmd`, not `README.md`; rebuild with `devtools::build_readme()`
  (`just readme`).
- Use `just format` for routine R formatting. `just lint` is an optional,
  on-demand scilintr diagnostic, not a prerequisite for every edit or completion.
- After changing notebook prose or figure inputs, render with
  `quarto render notebook/sc_analysis.qmd` when updating the HTML deliverable;
  it embeds image bytes. Scripts choose which figures to copy with
  `copy_notebook_figure()`; do not add notebook parsing, image hashing, or rollback
  machinery. Use `output_path()` at file writes for overwrite opt-in and directory
  creation. Do not maintain duplicate output inventories or write through symlinks.

## R conventions

- Prefer tidyverse solutions, `<-` assignment, and the native `|>` pipe.
- Use `theme_stone()` by default for every plot in this repo, including exploratory
  plots. Set it with `ggplot2::theme_set(theme_stone())` for plotting sessions.
- Keep plot text minimal: short axis labels with units, concise facet labels, and
  legends only when needed to decode mappings. Avoid subtitles, captions, verbose
  legend text, and redundant titles unless requested or essential to interpretation.
- Do not put routine statistical caveats or methodological disclaimers in plots
  (e.g. "descriptive, not mouse-level inference"). Omit reminders of familiar
  limitations; put genuinely necessary qualifications in analysis prose or chat,
  not on the image.
- Use a map function followed by `list_rbind()` or `list_cbind()`, not superseded
  purrr helpers such as `map_dfr()` or `map_dfc()` and their indexed variants.
- Prefer explicit namespaces such as `dplyr::filter()` over conflict-resolution
  side effects. Use bare or quoted columns, not `.data$`, in tidyselect calls.
- Keep scripts top-to-bottom and RStudio-step-friendly: purpose/usage, package
  loading, `# ---- parameters ----`, validation/work, then side effects. Avoid
  wrapping simple scripts in `main()` or trivial expressions in one-off helpers.

## Experimental design

- Conditions: `p27CKO` vs `p27CKO +EStim`.
- Replicates: six Mouse × Condition pseudobulk samples; mice 10 and 3 are paired,
  mouse 30 is E-Stim only, mouse 33 is control only.
- Platform: PipSeq V T2 with custom reference including eGFP.
- Primary statistical unit for condition-level DE: Mouse × Condition pseudobulk
  sample, not cell.
- BrdU: added at E-Stim + 24 h, washed after 48 h (E-Stim + 3 d), fixed at 5 d
  post E-Stim.

## Communication

Be concise, use plain language and active voice, and state assumptions,
decisions, inferences, and verification limits explicitly. Keep this file operational and
compact; do not turn it into a session journal.
