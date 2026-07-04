<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-07-04

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 4 sections | 2026-07-03 | ESPI R Package Shape, R and Documentation Workflow, Data and Figures, Statistical Unit |
| decisions.md | 16 entries | 2026-07-04 | Enable Mycelium without restructuring ESPI, Install bioinformatics conventions by default, Treat Mycelium restructure audit as advisory only, Enable skill-bridge after cloning available skillpacks, Use Rscript orchestration and Seurat-safe cluster branch tags |
| learnings.md | 13 entries | 2026-07-04 | Mycelium hook paths are local plugin-cache paths, Quarto embedded HTML must be rerendered after figure regeneration, Documented Autonomous-Science skillpack URL is unavailable, Seurat rewrites hyphens in reduction names, Pass clustree only cluster columns for large Seurat metadata |
| log/ | 5 sessions | 2026-07-04 | espi (5) |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-07-04 (heuristic)

## Tag clusters

- **mycelium** (15 entries) — D-12, D-13, D-14, D-15, D-16
- **hooks** (10 entries) — D-12, D-13, D-14, D-15, D-16
- **omp** (6 entries) — L-13, D-13, D-14, D-15, D-16
- **reproducibility** (6 entries) — L-10, L-11, D-5, D-9, D-12
- **session-state** (5 entries) — L-11, L-13, D-13, D-14, D-16
- **clustering** (4 entries) — L-4, D-5, D-6, D-11

## Most recent (10)

- [2026-07-04] L-8: Read TSV row counts with a count tool when output looks truncated
- [2026-07-04] L-9: Make DE output overwrites explicit
- [2026-07-04] L-10: Completed LOG_REGISTRY rows still need semantic fields
- [2026-07-04] L-11: Clear stale Mycelium reminder files after false-positive stop blocks
- [2026-07-04] L-12: OMP hook adapters must return modified tool content
- [2026-07-04] L-13: Stop hooks need session-boundary sentinel checks
- [2026-07-04] D-10: Use all Mouse × Condition samples as primary DE unit
- [2026-07-04] D-11: Keep MG-selected clustering at 30 PCs and resolution 0.3
- [2026-07-04] D-12: Keep repo-local Mycelium skills synced from OMP
- [2026-07-04] D-13: Bridge Mycelium hooks through OMP extensions

## By tag

- `mycelium`: L-1, L-3, L-10, L-11, L-12, L-13, D-1, D-2, D-3, D-4, D-12, D-13, D-14, D-15, D-16
- `hooks`: L-1, L-10, L-11, L-12, L-13, D-12, D-13, D-14, D-15, D-16
- `omp`: L-12, L-13, D-13, D-14, D-15, D-16
- `reproducibility`: L-9, L-10, L-11, D-5, D-9, D-12
- `session-state`: L-11, L-13, D-13, D-14, D-16
- `clustering`: L-4, D-5, D-6, D-11
- `r`: L-7, L-12, D-9, D-15
- `r-package`: L-6, D-1, D-3, D-8
- `marker-genes`: L-6, D-7, D-8
- `plotting`: L-5, L-6, D-8
- `seurat`: L-4, L-5, D-5
- `conventions`: D-2, D-4
- `data-lineage`: L-12, D-15
- `differential-expression`: L-9, D-10
- `figures`: L-2, L-7
- `mg-selected`: D-10, D-11
- `notebook`: L-2, L-7
- `repo-structure`: D-1, D-3
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
- `naming`: L-4
- `output-provenance`: L-9
- `pFlog`: D-11
- `package-data`: D-7
- `paired-design`: D-10
- `parameters`: D-11
- `portability`: L-1
- `pseudobulk`: D-10
- `quarto`: L-2
- `read-tool`: L-8
- `reductions`: L-4
- `reporting`: D-6
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
