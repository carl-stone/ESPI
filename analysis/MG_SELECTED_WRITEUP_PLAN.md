# MG-selected write-up plan

Date: 2026-07-13

## Results to carry forward

- The `mg-selected` branch starts from `cluster_pflog_no_filter_cc_dims20_res0.3` and removes clusters that meet the configured microglia/photoreceptor or high-`Cdkn1b` criteria.
- Excluded clusters: 2 (`marker_microglia;Cdkn1b_high`), 7 (`Cdkn1b_high`), and 8 (`Cdkn1b_high`). No cluster met the photoreceptor criterion.
- Retained cells: 3,460.
- Selected downstream clustering: PFlog, cell-cycle HVGs retained, 20 PCs, Leiden resolution 0.3.
- Primary DE design: DESeq2 pseudobulk, Mouse × Condition samples, `~ condition`, six samples total.
- Paired sensitivity: `~ mouse + condition`, paired mice 10 and 3 only, four samples total.
- Direction: positive log2 fold change means higher in `p27CKO +EStim` than `p27CKO`.
- Primary DE tested 24,612 genes and found 444 FDR-significant genes: 133 higher with E-Stim and 311 lower with E-Stim.
- Significant curated marker-list genes:
  - Higher with E-Stim: `Glul` (log2FC 0.769, padj 1.77e-05), `Ccn1` (0.603, 0.00455), `Hes6` (0.547, 0.0216), `Grm6` (0.563, 0.00730), and `Scgn` (0.531, 0.0207).
  - Lower with E-Stim: `Mcm2` (-0.983, 1.12e-06), `Mcm6` (-0.385, 0.0240), `Pcna` (-0.677, 0.00127), and `Rcvrn` (-0.626, 0.00184).
- Paired sensitivity found 140 FDR-significant genes and preserves `Glul`, `Ccn1`, `Mcm2`, and `Mcm6`; it also flags `Serpina3n`, but does not preserve `Hes6`, `Pcna`, `Grm6`, `Scgn`, or `Rcvrn` at FDR < 0.05.
- Primary differential detection tested 35,386 genes and found no FDR-significant genes or marker-list hits. Paired sensitivity tested 33,808 genes and found 44 FDR-significant genes but no marker-list hits. Interpret curated marker changes as pseudobulk abundance shifts, not binary on/off detection shifts.
- Upregulated ORA terms include apoptotic signaling and growth/stress-associated terms, but top terms also include narrow/redundant Nat8 acetylation and developmental/smooth-muscle/placenta labels that need cautious wording.
- Downregulated ORA and GSEA are dominated by cell cycle, chromosome segregation, nuclear division, and DNA replication.
- GSEA also shows positive NES terms related to catabolic metabolism, cellular respiration, mitochondrial organization, and macroautophagy.

## Working biological interpretation

E-Stim pushes the p27-low/negative MG-enriched branch toward a stress/reactive MG state with stronger MG/reactive marker abundance (`Glul`, `Ccn1`) and reduced proliferation/cell-cycle abundance (`Mcm2`, `Mcm6`, `Pcna`). There is weaker evidence for neurogenic, bipolar, and photoreceptor-associated signal because `Hes6`, `Grm6`, `Scgn`, and `Rcvrn` appear in the all-sample primary model but not the two-mouse paired sensitivity. The most robust transcriptome-scale result is cell-cycle downshift after E-Stim in this selected branch.

## Pipeline weaknesses to discuss in the notebook

1. **Selection conditions the biology.** Filtering `Cdkn1b`-high clusters means the branch is not an all-MG lineage analysis. Write "p27-low/negative MG-enriched branch," not "Müller glia broadly."
2. **Small and partly unmatched design.** The primary model uses all six Mouse × Condition pseudobulk samples but cannot model mouse pairing; the paired model uses only mice 10 and 3. Treat paired sensitivity as robustness evidence, not as a replacement primary result.
3. **Cluster filtering is thresholded.** Microglia/photoreceptor exclusion depends on module-score top class, top-score threshold, and score margin. This is defensible but should be framed as an operational selection rule.
4. **Marker-list dependence.** `AddModuleScore()` reflects the curated marker list. Ambiguous or sparse marker lists can shift cluster labels.
5. **Cell-cycle signal dominates.** Cell-cycle HVGs were retained in the chosen branch, and the strongest DE/enrichment signature is cell-cycle downregulation. This may be biology, composition, or both; do not over-interpret as a specific cell-cycle mechanism without validation.
6. **Pseudobulk abundance vs detection.** Primary detection testing found no FDR-significant genes. The paired sensitivity found 44 genes overall but no curated marker hits. Curated marker changes are quantitative abundance shifts within detected expression, not clear changes in the fraction of expressing cells.
7. **Low-count neuronal marker caution.** `Grm6` is significant in the primary model but has low counts and does not survive paired sensitivity. Present it as suggestive rod-bipolar-associated signal, not a strong fate-conversion claim.
8. **GO/GSEA redundancy and warnings.** GO terms are highly redundant, and GSEA emitted warnings about a small mapping failure rate, ties, and a few problematic pathways. Use GO as a compact summary of DEG themes, not independent validation.
9. **No `Cdkn1b`-retained sensitivity branch yet.** Strong claims about p27 biology need a sensitivity branch that keeps `Cdkn1b`-high clusters and compares conclusions.

## Plan for `notebook/sc_analysis.qmd`

### 1. Tighten the MG-selected clustering section

- Keep the selection-rule paragraph, but add one sentence up front: "This branch asks what remains after removing confident microglia and p27-high clusters, not what all MG-lineage cells do."
- Add a compact retained/excluded cluster summary near the UMAP.
- Keep the 20-PC/resolution-0.3 rationale, but make it shorter.

### 2. Add a results paragraph before the figures

Draft message:

> The retained cells remain enriched for Müller glia module scores, with smaller clusters carrying activated MG, cone-bipolar, or mixed marker signals. We therefore treat the branch as MG-enriched but not pure. The UMAP and marker heatmap are descriptive checks on branch composition; the statistical conclusions below use Mouse × Condition pseudobulk samples.

### 3. Expand the DE paragraph with effect sizes

- Add exact marker-hit effect sizes and adjusted p-values in one sentence or small inline table.
- Explicitly state 133 genes increased and 311 decreased with E-Stim.
- Separate robust hits from sensitivity-limited hits:
  - Robust in primary and paired sensitivity: `Glul`, `Ccn1`, `Mcm2`, `Mcm6`.
  - Primary-only/suggestive: `Hes6`, `Pcna`, `Grm6`, `Scgn`, `Rcvrn`.
  - Paired-only marker signal worth noting: `Serpina3n`.

### 4. Add a detection-analysis sentence immediately after marker DE

Draft message:

> Primary differential detection did not identify FDR-significant genes, including among curated markers. The paired sensitivity identified 44 genes overall but no curated marker hits. The curated marker shifts should therefore be read as pseudobulk abundance changes rather than clear changes in the fraction of cells expressing each marker.

### 5. Rewrite GO/GSEA as a guarded theme summary

- Downregulated side: emphasize chromosome segregation, DNA replication, nuclear division, and cell-cycle checkpoint terms; this is the cleanest enrichment result.
- Upregulated side: summarize as stress/growth/apoptotic signaling plus metabolism/mitochondrial/catabolic GSEA themes; avoid leaning on placenta/smooth-muscle labels as biology.
- Add one caveat sentence that GO redundancy and GSEA warnings make these secondary summaries.

### 6. Add a short "Interpretive limits" paragraph

Include:

- p27-low/negative branch only.
- limited replicate structure.
- no detection-fraction hits.
- no `Cdkn1b`-retained sensitivity yet.
- neuronal-like markers are suggestive, not proof of reprogramming.

### 7. Optional next analysis before manuscript claims

If we want stronger claims, run a `Cdkn1b`-retained sensitivity branch and compare:

- retained/excluded clusters,
- pseudobulk DE marker hits,
- cell-cycle GO/GSEA themes,
- `Cdkn1b`, reactive MG, neurogenic, and bipolar marker effects.

## Proposed notebook endpoint

The notebook should support this manuscript-facing claim:

> In the p27-low/negative MG-enriched branch, E-Stim is associated with stronger reactive/MG abundance signals and a broad reduction in proliferative/cell-cycle transcripts. Neurogenic and bipolar marker signals appear in the primary pseudobulk model, but their sensitivity to the paired-only analysis argues for cautious, hypothesis-generating wording rather than a definitive reprogramming claim.
