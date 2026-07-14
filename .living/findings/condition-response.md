---
topic: condition-response
description: Empirical findings about condition-associated expression responses within focused analysis branches. Entries distinguish Mouse × Condition pseudobulk evidence from descriptive cell-level marker rankings.
created: 2026-07-04
last_updated: 2026-07-14
status: active
---

# Condition Response

## F-001: E-Stim shifts curated marker abundance toward glial activation and away from proliferation
**Status:** supported
**Claim:** In the MG-selected branch, Mouse × Condition pseudobulk DE shows E-Stim-associated increases in `Glul` and `Ccn1` and decreases in `Mcm2` and `Mcm6`; these four curated marker effects are concordant and FDR-significant in both the all-sample primary design and the paired-mice sensitivity design.
**Implications:** The most design-consistent marker-level expression response in this selected branch is glial/activated-marker upregulation with reduced proliferative-marker abundance, not a cell-level marker-ranking result.
**Tags:** MG-selected, E-Stim, pseudobulk, curated markers, glial activation, proliferation

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-04 | notebook/sc_analysis.qmd MG-selected DE section; degs/mg_selected/deseq2_marker_overlap.tsv | MG-selected, 4,713 cells, 6 Mouse × Condition pseudobulk samples | ESPI | Primary `~ condition` DE tested 24,514 genes and found 453 FDR-significant genes; curated marker hits included increased `Glul` (log2FC 0.776, padj 1.57e-05) and `Ccn1` (0.606, padj 0.00441), and decreased `Mcm2` (-1.015, padj 4.59e-07) and `Mcm6` (-0.407, padj 0.0185). | supports |
| 2026-07-04 | degs/mg_selected/deseq2_paired_sensitivity_marker_overlap.tsv | MG-selected paired-mice sensitivity, mice 10 and 3 | ESPI | Paired `~ mouse + condition` sensitivity retained the same directions and FDR significance for `Glul` (log2FC 0.767, padj 0.00280), `Ccn1` (0.652, padj 0.0219), `Mcm2` (-0.854, padj 0.00117), and `Mcm6` (-0.975, padj 0.000246). | supports |
| 2026-07-13 | session 2026-07-13-003; `degs/mg_selected/deseq2_marker_overlap.tsv`; `deseq2_paired_sensitivity_marker_overlap.tsv` | MG-selected emptyDrops/log-MAD rebuild, 3,460 cells, 6 Mouse × Condition pseudobulk samples | ESPI | Primary `~ condition` DE found 444 FDR-significant genes and retained significant, concordant effects for `Glul` (log2FC 0.769), `Ccn1` (0.603), `Mcm2` (-0.983), and `Mcm6` (-0.385). Paired sensitivity retained the same directions and significance: 0.768, 0.652, -0.834, and -0.977 respectively. | supports |
| 2026-07-13 | session 2026-07-13-007; canonical `just run counts-qc true`; `degs/mg_selected/deseq2_marker_overlap.tsv` | MG-selected emptyDrops/scDblFinder/log-MAD rebuild, 3,456 cells, 6 Mouse × Condition pseudobulk samples | ESPI | Primary `~ condition` DE found 442 FDR-significant genes; paired sensitivity found 141. `Glul` (log2FC 0.770), `Ccn1` (0.603), `Mcm2` (-0.985), and `Mcm6` (-0.384) remained directionally concordant and FDR-significant in both designs. | supports |

### Open Questions
- Do the glial-activation and lower proliferative-abundance signals hold in a branch that is not conditioned on low/negative `Cdkn1b` and microglia-removal rules?
- Are the primary-only curated marker hits (`Hes6`, `Pcna`, `Grm6`) stable enough for interpretation after broader sensitivity checks?

