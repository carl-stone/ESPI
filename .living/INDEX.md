<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-07-09

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 6 sections | 2026-07-09 | ESPI R Package Shape, R and Documentation Workflow, Data and Figures, Statistical Unit, Cross-References |
| decisions.md | 40 entries (large — read selectively) | 2026-07-09 | Enable Mycelium without restructuring ESPI, Install bioinformatics conventions by default, Treat Mycelium restructure audit as advisory only, Enable skill-bridge after cloning available skillpacks, Use Rscript orchestration and Seurat-safe cluster branch tags |
| learnings.md | 33 entries (large — read selectively) | 2026-07-09 | Mycelium hook paths are local plugin-cache paths, Quarto embedded HTML must be rerendered after figure regeneration, Documented Autonomous-Science skillpack URL is unavailable, Seurat rewrites hyphens in reduction names, Pass clustree only cluster columns for large Seurat metadata |
| log/ | 46 sessions | 2026-07-09 | espi (46) |
| findings/ | 5 findings across 2 topics | 2026-07-04 | population-structure, condition-response |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-07-09 (heuristic)

## Tag clusters

- **mg-selected** (27 entries) — D-31, D-32, D-33, D-34, D-35
- **mycelium** (21 entries) — D-14, D-15, D-16, D-18, D-30
- **plotting** (18 entries) — D-27, D-29, D-31, D-32, D-33
- **hooks** (15 entries) — D-13, D-14, D-15, D-16, D-30
- **notebook** (13 entries) — D-26, D-27, D-28, D-29, D-33
- **provenance** (12 entries) — L-33, D-30, D-38, D-39, D-40

## Most recent (10)

- [2026-07-09] L-27: APFS case-insensitive paths can hide justfile casing mistakes
- [2026-07-09] L-28: Treat load_all as normal package loading
- [2026-07-09] L-29: Raw sample files may sit below the user-described input root
- [2026-07-09] L-30: Non-default raw-object handoffs require explicit preprocessing input
- [2026-07-09] L-31: Read10X integrity does not establish upstream PIPseeker cell calling
- [2026-07-09] L-32: A `^mt-` match can be a partial mitochondrial metric
- [2026-07-09] L-33: Keep barcode-rank plots diagnostic
- [2026-07-09] D-36: Use just as the ESPI command interface
- [2026-07-09] D-37: Derive raw-count inputs from `DATA_ROOT_DIR` and stop at an in-memory object
- [2026-07-09] D-38: Persist the raw Seurat object at the user-chosen input path

## By tag

- `mg-selected`: L-14, L-15, L-17, L-18, L-19, L-24, L-25, L-26, D-10, D-11, D-17, D-18, D-19, D-20, D-21, D-22, D-23, D-24, D-25, D-26, D-27, D-29, D-31, D-32, D-33, D-34, D-35
- `mycelium`: L-1, L-3, L-10, L-11, L-12, L-13, L-20, L-21, L-22, L-23, D-1, D-2, D-3, D-4, D-12, D-13, D-14, D-15, D-16, D-18, D-30
- `plotting`: L-5, L-6, L-18, L-24, L-25, L-26, L-33, D-8, D-20, D-21, D-22, D-23, D-24, D-27, D-29, D-31, D-32, D-33
- `hooks`: L-1, L-10, L-11, L-12, L-13, L-20, L-21, L-22, L-23, D-12, D-13, D-14, D-15, D-16, D-30
- `notebook`: L-2, L-7, L-16, L-25, D-20, D-21, D-22, D-23, D-26, D-27, D-28, D-29, D-33
- `provenance`: L-20, L-21, L-22, L-23, L-30, L-31, L-32, L-33, D-30, D-38, D-39, D-40
- `data-lineage`: L-12, L-29, L-30, L-31, L-32, D-15, D-37, D-38, D-39, D-40
- `reproducibility`: L-9, L-10, L-11, L-27, D-5, D-9, D-12, D-36, D-37, D-38
- `differential-expression`: L-9, L-18, L-24, D-10, D-17, D-20, D-32, D-34, D-35
- `differential-detection`: L-14, L-15, L-18, L-24, D-17, D-18, D-32, D-35
- `scripts`: L-29, L-30, D-5, D-9, D-36, D-37, D-38, D-39
- `omp`: L-12, L-13, D-13, D-14, D-15, D-16, D-30
- `session-logs`: L-10, L-20, L-21, L-22, L-23, D-28, D-30
- `validation`: L-8, L-19, L-29, L-31, L-32, D-39, D-40
- `marker-analysis`: L-16, L-17, D-19, D-22, D-23, D-31
- `r`: L-7, L-12, L-19, L-28, D-9, D-15
- `seurat`: L-4, L-5, L-17, L-26, D-5, D-19
- `marker-genes`: L-6, L-24, D-7, D-8, D-32
- `qc-filtering`: L-31, L-32, L-33, D-39, D-40
- `r-package`: L-6, D-1, D-3, D-8, D-36
- `session-state`: L-11, L-13, D-13, D-14, D-16
- `clustering`: L-4, D-5, D-6, D-11
- `todo`: L-23, D-27, D-34, D-35
- `umap`: L-25, L-26, D-21, D-33
- `abundance`: D-24, D-25, D-26
- `conventions`: D-2, D-4, D-28
- `filtering`: L-14, L-15, D-17
- `paired-design`: D-10, D-25, D-26
- `reporting`: L-16, D-6, D-28
- `single-cell`: D-2, D-7, D-37
- `tooling`: L-19, L-27, D-36
- `workflow`: L-28, L-30, L-33
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
- `contrast-labels`: D-29
- `cross-references`: D-28
- `documentation`: L-28
- `fisher-exact`: D-24
- `github`: L-3
- `interactivity`: D-9
- `last-session`: L-22
- `limma`: L-14
- `macos`: L-27
- `naming`: L-4
- `output-provenance`: L-9
- `pFlog`: D-11
- `package-data`: D-7
- `palette`: D-29
- `parameters`: D-11
- `pflog`: L-17
- `planning`: D-27
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
