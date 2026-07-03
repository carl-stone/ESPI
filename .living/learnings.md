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

### [2026-07-03] Seurat rewrites hyphens in reduction names

**Tags**: seurat, clustering, reductions, naming

**Category**: Pipeline robustness

**What happened**: `RunUMAP(reduction.name = "umap_pflog_no-filter-cc_dims20")` stored the reduction as `umap_pflog_no.filter.cc_dims20`, so a later lookup for the original hyphenated name failed.

**Why it matters**: Cluster artifact names that look valid as file names may not be valid as stable Seurat object keys.

**Resolution**: Use underscore-only branch tags for Seurat reductions and metadata columns, and keep a tripwire that requires branch tags to match `[A-Za-z0-9_]+`.

**mitigation_type**: code-and-tripwire

**structural_mitigation_candidate**: Centralize branch-tag construction if more scripts start creating Seurat object keys.

### [2026-07-03] Pass clustree only cluster columns for large Seurat metadata

**Tags**: clustree, seurat, plotting, warnings

**Category**: Plotting robustness

**What happened**: Passing the full Seurat object into `clustree()` emitted repeated warnings about unrelated metadata column name sanitization.

**Why it matters**: Those warnings obscure real plotting problems and repeat once per panel in multi-panel clustree figures.

**Resolution**: Build clustree plots from a metadata data frame containing only the cluster columns matching the requested prefix.

**mitigation_type**: code

**structural_mitigation_candidate**: Keep plotting helpers narrow: pass only the columns a plotting package needs when the source object has large heterogeneous metadata.
