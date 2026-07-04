<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-07-04

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 4 sections | 2026-07-04 | ESPI R Package Shape, R and Documentation Workflow, Data and Figures, Statistical Unit |
| decisions.md | 19 entries | 2026-07-04 | Enable Mycelium without restructuring ESPI, Install bioinformatics conventions by default, Treat Mycelium restructure audit as advisory only, Enable skill-bridge after cloning available skillpacks, Use Rscript orchestration and Seurat-safe cluster branch tags |
| learnings.md | 17 entries | 2026-07-04 | Mycelium hook paths are local plugin-cache paths, Quarto embedded HTML must be rerendered after figure regeneration, Documented Autonomous-Science skillpack URL is unavailable, Seurat rewrites hyphens in reduction names, Pass clustree only cluster columns for large Seurat metadata |
| log/ | 13 sessions | 2026-07-04 | espi (13) |
| findings/ | 5 findings across 2 topics | 2026-07-04 | population-structure, condition-response |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-07-04 (heuristic)

## Tag clusters

- **mycelium** (16 entries) — D-13, D-14, D-15, D-16, D-18
- **hooks** (10 entries) — D-12, D-13, D-14, D-15, D-16
- **mg-selected** (8 entries) — D-10, D-11, D-17, D-18, D-19
- **omp** (6 entries) — L-13, D-13, D-14, D-15, D-16
- **reproducibility** (6 entries) — L-10, L-11, D-5, D-9, D-12
- **session-state** (5 entries) — L-11, L-13, D-13, D-14, D-16

## Most recent (10)

- [2026-07-04] L-8: Read TSV row counts with a count tool when output looks truncated
- [2026-07-04] L-9: Make DE output overwrites explicit
- [2026-07-04] L-10: Completed LOG_REGISTRY rows still need semantic fields
- [2026-07-04] L-11: Clear stale Mycelium reminder files after false-positive stop blocks
- [2026-07-04] L-12: OMP hook adapters must return modified tool content
- [2026-07-04] L-13: Stop hooks need session-boundary sentinel checks
- [2026-07-04] L-14: Differential-detection gene filters can change empirical-Bayes behavior
- [2026-07-04] L-15: Muscat DD changes the paired-sensitivity hit set
- [2026-07-04] L-16: Do not report planned marker analyses as completed
- [2026-07-04] L-17: Seurat FoldChange defaults are invalid for PFlog marker ranking

## By tag

- `mycelium`: L-1, L-3, L-10, L-11, L-12, L-13, D-1, D-2, D-3, D-4, D-12, D-13, D-14, D-15, D-16, D-18
- `hooks`: L-1, L-10, L-11, L-12, L-13, D-12, D-13, D-14, D-15, D-16
- `mg-selected`: L-14, L-15, L-17, D-10, D-11, D-17, D-18, D-19
- `omp`: L-12, L-13, D-13, D-14, D-15, D-16
- `reproducibility`: L-9, L-10, L-11, D-5, D-9, D-12
- `session-state`: L-11, L-13, D-13, D-14, D-16
- `seurat`: L-4, L-5, L-17, D-5, D-19
- `clustering`: L-4, D-5, D-6, D-11
- `differential-detection`: L-14, L-15, D-17, D-18
- `r`: L-7, L-12, D-9, D-15
- `r-package`: L-6, D-1, D-3, D-8
- `differential-expression`: L-9, D-10, D-17
- `filtering`: L-14, L-15, D-17
- `marker-analysis`: L-16, L-17, D-19
- `marker-genes`: L-6, D-7, D-8
- `notebook`: L-2, L-7, L-16
- `plotting`: L-5, L-6, D-8
- `conventions`: D-2, D-4
- `data-lineage`: L-12, D-15
- `figures`: L-2, L-7
- `findallmarkers`: L-17, D-19
- `muscat`: L-15, D-18
- `repo-structure`: D-1, D-3
- `reporting`: L-16, D-6
- `scripts`: D-5, D-9
- `setup`: L-1, L-3
- `single-cell`: D-2, D-7
- `skillpacks`: L-3, D-4
- `annotation`: D-7
- `audit`: D-3
- `bioinformatics`: D-2
- `box`: L-2
- `clustree`: L-5
- `github`: L-3
- `heatmap`: D-8
- `interactivity`: D-9
- `limma`: L-14
- `naming`: L-4
- `output-provenance`: L-9
- `pFlog`: D-11
- `package-data`: D-7
- `paired-design`: D-10
- `parameters`: D-11
- `pflog`: L-17
- `portability`: L-1
- `pseudobulk`: D-10
- `quarto`: L-2
- `read-tool`: L-8
- `reductions`: L-4
- `review`: L-16
- `sensitivity`: D-6
- `session-logs`: L-10
- `skill-bridge`: D-4
- `skills`: D-12
- `supplemental-figures`: D-6
- `symlink`: L-7
- `tests`: L-13
- `tsv`: L-8
- `unicode`: L-6
- `validation`: L-8
- `warnings`: L-5

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
