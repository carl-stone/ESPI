# AGENTS.md - ESPI paper revision analysis

------------------------------------------------------------------------

## Eternal guidelines

Keep the AGENTS.md file operational and compact; repo code is the source of truth. If a scientific detail is absent, ask the user for clarification rather than making assumptions.

This repo is a minimal R package for a scientific analysis, not a general-purpose R package. 
It is not intended for public use beyond code documentation and reproducibility.
Never worry about backward compatibility; this is a one-off analysis. The code is not intended to be reused in other projects.
Helpers and analysis modules live in `R/`; the executable analysis pipeline lives in `scripts/`.

Load the package with `devtools::load_all()`. Run `devtools::document()` after modifying any code in `R/` to update the package documentation. NEVER edit .Rd files directly; they are generated from the R code and roxygen2 comments.

NEVER edit README.md directly; it is generated from README.Rmd. After you edit README.Rmd, run `devtools::build_readme()` to update README.md.

Update `README.Rmd`, `AGENTS.md`, and other package documentation as needed to reflect changes in the analysis pipeline.

------------------------------------------------------------------------

## Coding conventions

Follow these important style rules when writing R code:

- Prefer solutions that use {tidyverse}
- Always use `<-` for assignment
- Always use the native base-R pipe `|>` for piped expressions
- This package uses Air for code formatting; use `air format` in bash to format code.
- This package uses `conflicted` for function name conflicts. Prefer explicit namespace calls (e.g., `dplyr::filter()`) over `conflicted` resolution.

------------------------------------------------------------------------

## Communication style

- Be **concise** with all answers and writing.
- Use **plain language** and avoid unnecessary jargon.
- Use **active voice** and avoid passive constructions.

------------------------------------------------------------------------

## Experimental design

- Conditions: `p27CKO` vs `p27CKO +EStim`.
- Replicates: six Mouse × Condition pseudobulk samples; mice 10 and 3 are paired, mouse 30 is E-Stim only, mouse 33 is control only.
- Platform: PipSeq V T2 with custom reference including eGFP.
- Primary statistical unit for condition-level DE: Mouse × Condition pseudobulk sample, not cell.
- BrdU: added at E-Stim + 24 h, washed after 48 h (E-Stim + 3 d), fixed at 5 d post E-Stim.

---
