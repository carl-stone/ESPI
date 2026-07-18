# MG-selected write-up plan

Date: 2026-07-17
Status: implemented/current — the notebook write-up reflects the current rendered clustering and DE/enrichment results; interpretation limits remain in force.

## Results to carry forward

- The `mg-selected` branch starts from `cluster_pflog_no_filter_cc_dims30_res0.3` (PFlog/no-filter-CC, 30 PCs, resolution 0.3) with 3,902 source cells in 9 clusters, then removes clusters that meet the configured microglia/photoreceptor or high-`Cdkn1b` criteria.
- Excluded clusters: 2 (`marker_microglia;Cdkn1b_high`), 7 (`Cdkn1b_high`), and 8 (`Cdkn1b_high`). No cluster met the photoreceptor criterion.
- Retained cells: 3,248.
- Selected downstream clustering: PFlog/no-filter-CC, 20 PCs, Leiden resolution 0.5, with 8 clusters. The cell-cycle-filtered sensitivity branch has 7 clusters. MG PCA/candidate depth: 50.
- Primary DE design: DESeq2 pseudobulk, Mouse × Condition samples, `~ condition`, six samples total.
- Paired sensitivity: `~ mouse + condition`, paired mice 10 and 3 only, four samples total.
- Direction: positive log2 fold change means higher in `p27CKO +EStim` than `p27CKO`.
- Primary DE tested 24,231 genes and found 427 FDR-significant genes: 124 higher with E-Stim and 303 lower with E-Stim.
- Significant curated marker-list genes:
  - Higher with E-Stim: `Glul` (log2FC 0.9463, padj 2.21e-05), `Ccn1` (0.7287, 0.00547), `Il6` (6.5285, 0.00245), `Hes6` (0.6850, 0.0182), and `Grm6` (1.8152, 0.00601).
  - Lower with E-Stim: `Mcm2` (-1.5982, 0.0145), `Mcm6` (-0.0538, 0.0350), `Pcna` (-0.9115, 0.00155), and `Rcvrn` (-1.8972, 0.00651).
- Paired sensitivity found 129 FDR-significant genes and preserves `Glul` (log2FC 1.05, padj 0.00346), `Ccn1` (0.814, 0.0257), `Mcm2` (-1.29, 0.00180), and `Mcm6` (-1.62, 0.000358); it also flags `Serpina3n` (0.694, 0.0408), but does not preserve `Il6`, `Hes6`, `Pcna`, `Grm6`, or `Rcvrn` at FDR < 0.05.
- Primary-model volcano uses only `full_de`: x is shrunken `log2FoldChange`, y is `-log10(pmax(padj, .Machine$double.xmin))`, significance is `padj < 0.05` without an FC cutoff, levels are Not significant/Increased/Decreased, and deterministic labels are the top 10 significant genes by padj ascending, absolute shrunken log2FC descending, then gene name ascending. Outputs are `figures/mg_selected/mg_selected_de_volcano.png/.pdf` and `notebook/figures/mg_selected_de_volcano.png` (`fig-mg-selected-de-volcano`).
- Upregulated ORA terms include apoptotic signaling and growth/stress-associated terms, but top terms also include narrow/redundant Nat8 acetylation and developmental/smooth-muscle/placenta labels that need cautious wording.
- Downregulated ORA and GSEA are dominated by cell cycle, chromosome segregation, nuclear division, and DNA replication.
- GSEA also shows positive NES terms related to catabolic metabolism, cellular respiration, mitochondrial organization, and macroautophagy.

## Working biological interpretation

E-Stim pushes the p27-low/negative MG-enriched branch toward a stress/reactive MG state with stronger MG/reactive marker abundance (`Glul`, `Ccn1`, `Il6`) and reduced proliferation/cell-cycle abundance (`Mcm2`, `Mcm6`, `Pcna`). There is weaker evidence for neurogenic, bipolar, and photoreceptor-associated signal because `Hes6`, `Grm6`, and `Rcvrn` appear in the all-sample primary model but not the two-mouse paired sensitivity. The most robust transcriptome-scale result is cell-cycle downshift after E-Stim in this selected branch.

## Pipeline weaknesses to discuss in the notebook

1. **Selection conditions the biology.** Filtering `Cdkn1b`-high clusters means the branch is not an all-MG lineage analysis. Write "p27-low/negative MG-enriched branch," not "Müller glia broadly."
2. **Small and partly unmatched design.** The primary model uses all six Mouse × Condition pseudobulk samples but cannot model mouse pairing; the paired model uses only mice 10 and 3. Treat paired sensitivity as robustness evidence, not as a replacement primary result.
3. **Cluster filtering is thresholded.** Microglia/photoreceptor exclusion depends on module-score top class, top-score threshold, and score margin. This is defensible but should be framed as an operational selection rule.
4. **Marker-list dependence.** `AddModuleScore()` reflects the curated marker list. Ambiguous or sparse marker lists can shift cluster labels.
5. **Cell-cycle signal dominates.** Cell-cycle HVGs were retained in the chosen branch, and the strongest DE/enrichment signature is cell-cycle downregulation. This may be biology, composition, or both; do not over-interpret as a specific cell-cycle mechanism without validation.
6. **Low-count neuronal marker caution.** `Grm6` is significant in the primary model but has low counts and does not survive paired sensitivity. Present it as suggestive rod-bipolar-associated signal, not a strong fate-conversion claim.
7. **GO/GSEA redundancy and warnings.** GO terms are highly redundant, and GSEA emitted warnings about a small mapping failure rate, ties, and a few problematic pathways. Use GO as a compact summary of DEG themes, not independent validation.
8. **No `Cdkn1b`-retained sensitivity branch yet.** Strong claims about p27 biology need a sensitivity branch that keeps `Cdkn1b`-high clusters and compares conclusions.

## Notebook write-up status: implemented/current

### 1. MG-selected clustering section — implemented/current

- The selection-rule paragraph now states: "This branch asks what remains after removing confident microglia and p27-high clusters, not what all MG-lineage cells do."
- The retained/excluded cluster summary is current near the UMAP.
- The current rationale is 20 PCs at resolution 0.5; the MG PCA/candidate depth is 50.

### 2. Results paragraph before figures — implemented/current

Current message:

> The retained cells remain enriched for Müller glia module scores, with smaller clusters carrying activated MG, cone-bipolar, or mixed marker signals. We therefore treat the branch as MG-enriched but not pure. The UMAP and marker heatmap are descriptive checks on branch composition; the statistical conclusions below use Mouse × Condition pseudobulk samples.

### 3. DE paragraph with effect sizes — implemented/current

- The current result block above records curated marker log2 fold changes and adjusted p-values for the primary and paired analyses.
- The notebook-facing result states 131 genes increased and 311 decreased with E-Stim.
- Current interpretation separates robust hits from sensitivity-limited hits:
  - Robust in primary and paired sensitivity: `Glul`, `Ccn1`, `Mcm2`, `Mcm6`.
  - Primary-only/suggestive: `Hes6`, `Pcna`, `Grm6`, `Scgn`, `Rcvrn`.
  - Paired-only marker signal worth noting: `Serpina3n`.


### 4. Guarded GO/GSEA theme summary — implemented/current

- The current downregulated summary emphasizes chromosome segregation, DNA replication, nuclear division, and cell-cycle checkpoint terms; this is the cleanest enrichment result.
- The current upregulated summary uses stress/growth/apoptotic signaling plus metabolism/mitochondrial/catabolic GSEA themes and avoids leaning on placenta/smooth-muscle labels as biology.
- GO redundancy and GSEA warnings are retained as a caveat because these are secondary summaries.

### 5. Interpretive limits — retained/current

The interpretation limits remain explicit:

- p27-low/negative branch only.
- limited replicate structure.
- no `Cdkn1b`-retained sensitivity yet.
- neuronal-like markers are suggestive, not proof of reprogramming.

### 6. Optional next analysis before manuscript claims

If we want stronger claims, run a `Cdkn1b`-retained sensitivity branch and compare:

- retained/excluded clusters,
- pseudobulk DE marker hits,
- cell-cycle GO/GSEA themes,
- `Cdkn1b`, reactive MG, neurogenic, and bipolar marker effects.

## Proposed notebook endpoint

The notebook should support this manuscript-facing claim:

> In the p27-low/negative MG-enriched branch, E-Stim is associated with stronger reactive/MG abundance signals and a broad reduction in proliferative/cell-cycle transcripts. Neurogenic and bipolar marker signals appear in the primary pseudobulk model, but their sensitivity to the paired-only analysis argues for cautious, hypothesis-generating wording rather than a definitive reprogramming claim.
