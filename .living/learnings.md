# Learnings

Append-only log of gotchas, surprises, and reusable workflow lessons.

### [2026-07-03] Mycelium hook paths are local plugin-cache paths

**Tags**: mycelium, hooks, setup, portability

**Category**: Tooling setup

**What happened**: Mycelium generated `.claude/settings.local.json` with hook commands pointing at `/Users/carlstone/.omp/plugins/cache/plugins/mycelium___mycelium___0.0.0/`.

**Why it matters**: These paths can break after plugin upgrades, cache cleanup, or a fresh clone on another machine.

**Resolution**: Keep `.claude/` gitignored and document the need to rerun Mycelium hook setup after plugin changes.

**mitigation_type**: documentation

**structural_mitigation_candidate**: Prefer stable project-local hook wrappers if Mycelium supports them later.

### [2026-07-03] Quarto embedded HTML must be rerendered after figure regeneration

**Tags**: quarto, figures, box, notebook

**Category**: Reproducible reporting

**What happened**: `notebook/sc_analysis.qmd` uses `embed-resources: true`, so regenerated Box figures do not update `notebook/sc_analysis.html` until Quarto rerenders the document.

**Why it matters**: Symlink targets can be current while the HTML still shows stale embedded image bytes.

**Resolution**: Run `quarto render notebook/sc_analysis.qmd` after regenerating source figures.

**mitigation_type**: workflow

**structural_mitigation_candidate**: Add a notebook render check to any future figure-regeneration script.

### [2026-07-03] Documented Autonomous-Science skillpack URL is unavailable

**Tags**: mycelium, skillpacks, github, setup

**Category**: Tooling setup

**What happened**: `git clone https://github.com/arjunrajlaboratory/Autonomous-Science.git` failed with `Repository not found`, and GitHub/web searches found no matching public repository.

**Why it matters**: The Mycelium skillpacks README references a third external skill library that cannot currently be installed from the documented URL.

**Resolution**: Clone the two available libraries, install `skill-bridge`, and document the missing repository in `skillpacks/README.md`.

**mitigation_type**: documentation

**structural_mitigation_candidate**: Update the Mycelium plugin's skillpack template if a corrected `Autonomous-Science` URL becomes available.
