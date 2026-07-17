<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-07-17

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 6 sections | 2026-07-17 | ESPI R Package Shape, R and Documentation Workflow, Data and Figures, Statistical Unit, Cross-References |
| decisions.md | 62 entries (large — read selectively) | 2026-07-16 | Enable Mycelium without restructuring ESPI, Install bioinformatics conventions by default, Treat Mycelium restructure audit as advisory only, Enable skill-bridge after cloning available skillpacks, Use Rscript orchestration and Seurat-safe cluster branch tags |
| learnings.md | 50 entries (large — read selectively) | 2026-07-16 | Mycelium hook paths are local plugin-cache paths, Quarto embedded HTML must be rerendered after figure regeneration, Documented Autonomous-Science skillpack URL is unavailable, Seurat rewrites hyphens in reduction names, Pass clustree only cluster columns for large Seurat metadata |
| log/ | 73 sessions | 2026-07-16 | espi (73) |
| findings/ | 4 findings across 2 topics | 2026-07-14 | condition-response, population-structure |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-07-17 (heuristic)

## Tag clusters

- **mg-selected** (35 entries) — D-47, D-50, D-52, D-54, D-55
- **reproducibility** (29 entries) — D-53, D-55, D-59, D-60, D-61
- **plotting** (26 entries) — D-52, D-54, D-55, D-58, D-62
- **mycelium** (21 entries) — D-14, D-15, D-16, D-18, D-30
- **data-lineage** (17 entries) — D-41, D-44, D-45, D-46, D-49
- **notebook** (17 entries) — D-27, D-28, D-29, D-33, D-62

## Most recent (10)

- [2026-07-16] L-47: Separate source contracts from visible sensitivity oracles
- [2026-07-16] L-48: Publication cutovers need layered equivalence checks
- [2026-07-16] L-49: Review maintenance paths separately from routine runs
- [2026-07-16] L-50: Generated session inventories can include non-repository paths
- [2026-07-16] D-59: Freeze all existing Seurat objects and designate two final objects
- [2026-07-16] D-60: Define boundaries for the publication-pipeline restructuring proposal
- [2026-07-16] D-61: Adopt the clean four-phase publication cutover
- [2026-07-16] D-62: Resolve the source-versus-visible MG heatmap oracle conflict
- [2026-07-15] L-45: Pipeline interface changes require both contract tripwires
- [2026-07-15] L-46: ComplexHeatmap cannot independently retain the split-mean dendrogram

## By tag

- `mg-selected`: L-14, L-15, L-17, L-18, L-19, L-24, L-25, L-26, L-37, L-38, D-10, D-11, D-17, D-18, D-19, D-20, D-21, D-22, D-23, D-24, D-25, D-26, D-27, D-29, D-31, D-32, D-33, D-34, D-35, D-43, D-47, D-50, D-52, D-54, D-55
- `reproducibility`: L-9, L-10, L-11, L-27, L-34, L-43, L-48, D-5, D-9, D-12, D-36, D-37, D-38, D-41, D-42, D-43, D-44, D-45, D-46, D-47, D-48, D-50, D-51, D-52, D-53, D-55, D-59, D-60, D-61
- `plotting`: L-5, L-6, L-18, L-24, L-25, L-26, L-33, L-35, L-42, L-46, D-8, D-20, D-21, D-22, D-23, D-24, D-27, D-29, D-31, D-32, D-33, D-52, D-54, D-55, D-58, D-62
- `mycelium`: L-1, L-3, L-10, L-11, L-12, L-13, L-20, L-21, L-22, L-23, D-1, D-2, D-3, D-4, D-12, D-13, D-14, D-15, D-16, D-18, D-30
- `data-lineage`: L-12, L-29, L-30, L-31, L-32, L-34, L-36, D-15, D-37, D-38, D-39, D-40, D-41, D-44, D-45, D-46, D-49
- `notebook`: L-2, L-7, L-16, L-25, L-47, L-48, L-49, D-20, D-21, D-22, D-23, D-26, D-27, D-28, D-29, D-33, D-62
- `hooks`: L-1, L-10, L-11, L-12, L-13, L-20, L-21, L-22, L-23, L-50, D-12, D-13, D-14, D-15, D-16, D-30
- `provenance`: L-20, L-21, L-22, L-23, L-30, L-31, L-32, L-33, L-41, D-30, D-38, D-39, D-40, D-41
- `differential-expression`: L-9, L-18, L-24, L-42, D-10, D-17, D-20, D-32, D-34, D-35, D-52, D-54, D-55
- `scripts`: L-29, L-30, L-34, D-5, D-9, D-36, D-37, D-38, D-39, D-41, D-42, D-48, D-60
- `clustering`: L-4, L-37, L-38, L-46, D-5, D-6, D-11, D-43, D-47, D-50, D-58, D-59
- `validation`: L-8, L-19, L-29, L-31, L-32, L-36, L-37, L-38, L-48, D-39, D-40, D-51
- `workflow`: L-28, L-30, L-33, L-34, L-37, L-44, L-45, D-42, D-48, D-51
- `qc-filtering`: L-31, L-32, L-33, L-36, D-39, D-40, D-45, D-46, D-49
- `differential-detection`: L-14, L-15, L-18, L-24, D-17, D-18, D-32, D-35
- `seurat`: L-4, L-5, L-17, L-26, L-35, D-5, D-19, D-59
- `omp`: L-12, L-13, D-13, D-14, D-15, D-16, D-30
- `session-logs`: L-10, L-20, L-21, L-22, L-23, D-28, D-30
- `heatmap`: L-46, L-47, D-8, D-31, D-58, D-62
- `marker-analysis`: L-16, L-17, D-19, D-22, D-23, D-31
- `pflog`: L-17, L-35, L-37, D-43, D-47, D-50
- `pipeline`: L-45, L-47, L-48, L-49, D-60, D-61
- `r`: L-7, L-12, L-19, L-28, D-9, D-15
- `marker-genes`: L-6, L-24, D-7, D-8, D-32
- `r-package`: L-6, D-1, D-3, D-8, D-36
- `session-state`: L-11, L-13, D-13, D-14, D-16
- `enrichment`: L-43, D-34, D-35, D-55
- `just`: L-27, L-28, L-39, D-51
- `mitochondrial`: L-32, L-36, D-40, D-45
- `todo`: L-23, D-27, D-34, D-35
- `tooling`: L-19, L-27, D-36, D-56
- `tripwires`: L-40, L-41, L-45, D-53
- `umap`: L-25, L-26, D-21, D-33
- `abundance`: D-24, D-25, D-26
- `conventions`: D-2, D-4, D-28
- `emptydrops`: D-46, D-47, D-49
- `filtering`: L-14, L-15, D-17
- `metadata`: D-44, D-46, D-49
- `orchestration`: L-38, L-39, D-51
- `paired-design`: D-10, D-25, D-26
- `planning`: D-27, D-48, D-60
- `preprocessing`: L-35, L-40, D-53
- `reporting`: L-16, D-6, D-28
- `sensitivity`: L-47, D-6, D-62
- `single-cell`: D-2, D-7, D-37
- `skillpacks`: L-3, D-4, D-56
- `cli`: L-44, L-45
- `collaboration`: D-34, D-35
- `condition`: L-25, D-33
- `doublets`: L-36, D-45
- `figures`: L-2, L-7
- `findallmarkers`: L-17, D-19
- `label-permutation`: L-40, D-53
- `lint`: L-44, D-56
- `muscat`: L-15, D-18
- `pipseeker`: L-31, L-33
- `portability`: L-1, L-29
- `randomization`: D-26, D-31
- `repo-structure`: D-1, D-3
- `review`: L-16, L-22
- `scilintr`: L-44, D-56
- `setup`: L-1, L-3
- `tests`: L-13, L-21
- `annotation`: D-7
- `architecture`: D-61
- `audit`: D-3
- `bioinformatics`: D-2
- `box`: L-2
- `causal-inference`: D-25
- `clr`: D-24
- `clusterprofiler`: L-43
- `clustree`: L-5
- `compositional-data`: D-25
- `contrast`: D-44
- `contrast-labels`: D-29
- `cross-references`: D-28
- `dependencies`: L-49
- `documentation`: L-28
- `dry-run`: L-39
- `enrichit`: L-43
- `fdr`: L-42
- `fisher-exact`: D-24
- `frozen`: D-59
- `git`: L-50
- `github`: L-3
- `interactivity`: D-9
- `last-session`: L-22
- `limma`: L-14
- `macos`: L-27
- `markdown`: L-41
- `naming`: L-4
- `output-provenance`: L-9
- `overwrite`: L-49
- `p-value`: L-42
- `pFlog`: D-11
- `package-data`: D-7
- `palette`: D-29
- `parameters`: D-11
- `parsing`: L-41
- `pseudobulk`: D-10
- `quarto`: L-2
- `read-tool`: L-8
- `reductions`: L-4
- `sample-identity`: L-40
- `scDblFinder`: D-49
- `session-log`: L-50
- `skill-bridge`: D-4
- `skills`: D-12
- `stability`: D-50
- `supplemental-figures`: D-6
- `symlink`: L-7
- `testing`: L-39
- `tool-devices`: L-50
- `tsv`: L-8
- `unicode`: L-6
- `warnings`: L-5

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
