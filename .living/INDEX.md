<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-07-12

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 6 sections | 2026-07-09 | ESPI R Package Shape, R and Documentation Workflow, Data and Figures, Statistical Unit, Cross-References |
| decisions.md | 44 entries (large — read selectively) | 2026-07-12 | Enable Mycelium without restructuring ESPI, Install bioinformatics conventions by default, Treat Mycelium restructure audit as advisory only, Enable skill-bridge after cloning available skillpacks, Use Rscript orchestration and Seurat-safe cluster branch tags |
| learnings.md | 35 entries (large — read selectively) | 2026-07-12 | Mycelium hook paths are local plugin-cache paths, Quarto embedded HTML must be rerendered after figure regeneration, Documented Autonomous-Science skillpack URL is unavailable, Seurat rewrites hyphens in reduction names, Pass clustree only cluster columns for large Seurat metadata |
| log/ | 51 sessions | 2026-07-12 | espi (51) |
| findings/ | 5 findings across 2 topics | 2026-07-04 | population-structure, condition-response |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-07-12 (heuristic)

## Tag clusters

- **mg-selected** (28 entries) — D-32, D-33, D-34, D-35, D-43
- **mycelium** (21 entries) — D-14, D-15, D-16, D-18, D-30
- **plotting** (19 entries) — D-27, D-29, D-31, D-32, D-33
- **hooks** (15 entries) — D-13, D-14, D-15, D-16, D-30
- **reproducibility** (15 entries) — D-38, D-41, D-42, D-43, D-44
- **data-lineage** (13 entries) — D-38, D-39, D-40, D-41, D-44

## Most recent (10)

- [2026-07-12] L-34: Source selection belongs at preprocessing
- [2026-07-12] L-35: PFlog PCA diagnostics must use the PCA source layer
- [2026-07-12] D-41: Select preprocessing source at the pipeline entrypoint
- [2026-07-12] D-42: Number executable pipeline stages
- [2026-07-12] D-43: Use 20 PCs and Leiden resolution 0.3 for the counts-derived PFlog run
- [2026-07-12] D-44: Canonicalize the E-Stim condition label at count ingestion
- [2026-07-09] L-27: APFS case-insensitive paths can hide justfile casing mistakes
- [2026-07-09] L-28: Treat load_all as normal package loading
- [2026-07-09] L-29: Raw sample files may sit below the user-described input root
- [2026-07-09] L-30: Non-default raw-object handoffs require explicit preprocessing input

## By tag

- `mg-selected`: L-14, L-15, L-17, L-18, L-19, L-24, L-25, L-26, D-10, D-11, D-17, D-18, D-19, D-20, D-21, D-22, D-23, D-24, D-25, D-26, D-27, D-29, D-31, D-32, D-33, D-34, D-35, D-43
- `mycelium`: L-1, L-3, L-10, L-11, L-12, L-13, L-20, L-21, L-22, L-23, D-1, D-2, D-3, D-4, D-12, D-13, D-14, D-15, D-16, D-18, D-30
- `plotting`: L-5, L-6, L-18, L-24, L-25, L-26, L-33, L-35, D-8, D-20, D-21, D-22, D-23, D-24, D-27, D-29, D-31, D-32, D-33
- `hooks`: L-1, L-10, L-11, L-12, L-13, L-20, L-21, L-22, L-23, D-12, D-13, D-14, D-15, D-16, D-30
- `reproducibility`: L-9, L-10, L-11, L-27, L-34, D-5, D-9, D-12, D-36, D-37, D-38, D-41, D-42, D-43, D-44
- `data-lineage`: L-12, L-29, L-30, L-31, L-32, L-34, D-15, D-37, D-38, D-39, D-40, D-41, D-44
- `notebook`: L-2, L-7, L-16, L-25, D-20, D-21, D-22, D-23, D-26, D-27, D-28, D-29, D-33
- `provenance`: L-20, L-21, L-22, L-23, L-30, L-31, L-32, L-33, D-30, D-38, D-39, D-40, D-41
- `scripts`: L-29, L-30, L-34, D-5, D-9, D-36, D-37, D-38, D-39, D-41, D-42
- `differential-expression`: L-9, L-18, L-24, D-10, D-17, D-20, D-32, D-34, D-35
- `differential-detection`: L-14, L-15, L-18, L-24, D-17, D-18, D-32, D-35
- `omp`: L-12, L-13, D-13, D-14, D-15, D-16, D-30
- `session-logs`: L-10, L-20, L-21, L-22, L-23, D-28, D-30
- `seurat`: L-4, L-5, L-17, L-26, L-35, D-5, D-19
- `validation`: L-8, L-19, L-29, L-31, L-32, D-39, D-40
- `marker-analysis`: L-16, L-17, D-19, D-22, D-23, D-31
- `r`: L-7, L-12, L-19, L-28, D-9, D-15
- `clustering`: L-4, D-5, D-6, D-11, D-43
- `marker-genes`: L-6, L-24, D-7, D-8, D-32
- `qc-filtering`: L-31, L-32, L-33, D-39, D-40
- `r-package`: L-6, D-1, D-3, D-8, D-36
- `session-state`: L-11, L-13, D-13, D-14, D-16
- `workflow`: L-28, L-30, L-33, L-34, D-42
- `todo`: L-23, D-27, D-34, D-35
- `umap`: L-25, L-26, D-21, D-33
- `abundance`: D-24, D-25, D-26
- `conventions`: D-2, D-4, D-28
- `filtering`: L-14, L-15, D-17
- `paired-design`: D-10, D-25, D-26
- `pflog`: L-17, L-35, D-43
- `reporting`: L-16, D-6, D-28
- `single-cell`: D-2, D-7, D-37
- `tooling`: L-19, L-27, D-36
- `collaboration`: D-34, D-35
- `condition`: L-25, D-33
- `enrichment`: D-34, D-35
- `figures`: L-2, L-7
- `findallmarkers`: L-17, D-19
- `heatmap`: D-8, D-31
- `just`: L-27, L-28
- `mitochondrial`: L-32, D-40
- `muscat`: L-15, D-18
- `pipseeker`: L-31, L-33
- `portability`: L-1, L-29
- `randomization`: D-26, D-31
- `repo-structure`: D-1, D-3
- `review`: L-16, L-22
- `setup`: L-1, L-3
- `skillpacks`: L-3, D-4
- `tests`: L-13, L-21
- `annotation`: D-7
- `audit`: D-3
- `bioinformatics`: D-2
- `box`: L-2
- `causal-inference`: D-25
- `clr`: D-24
- `clustree`: L-5
- `compositional-data`: D-25
- `contrast`: D-44
- `contrast-labels`: D-29
- `cross-references`: D-28
- `documentation`: L-28
- `fisher-exact`: D-24
- `github`: L-3
- `interactivity`: D-9
- `last-session`: L-22
- `limma`: L-14
- `macos`: L-27
- `metadata`: D-44
- `naming`: L-4
- `output-provenance`: L-9
- `pFlog`: D-11
- `package-data`: D-7
- `palette`: D-29
- `parameters`: D-11
- `planning`: D-27
- `preprocessing`: L-35
- `pseudobulk`: D-10
- `quarto`: L-2
- `read-tool`: L-8
- `reductions`: L-4
- `sensitivity`: D-6
- `skill-bridge`: D-4
- `skills`: D-12
- `supplemental-figures`: D-6
- `symlink`: L-7
- `tsv`: L-8
- `unicode`: L-6
- `warnings`: L-5

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
