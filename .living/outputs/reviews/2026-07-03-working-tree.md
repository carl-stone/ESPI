# Review — working tree — 2026-07-03

**Scope**: Working tree first-party changes. External cloned skillpacks were excluded; copied Mycelium convention-pack contents were treated as context.
**Files reviewed**: 40 first-party files listed; 28 text files diffed; generated/binary artifacts listed but not body-reviewed.
**Sub-agents run**: 6 — stats-causal, data-pipeline-leakage, bioinformatics, llm-failure-modes, doc-schema-fidelity, code-quality.

## Key decisions in this analysis

The consequential analytical choices in this work. Some have associated findings below; others are informational so you can decide whether to revisit them.

- **Cell-cycle branch comparison** — The analysis now treats HVG-filtered and non-filtered cell-cycle branches as complementary preprocessing branches that must remain distinguishable through clustering and reporting. See F1.
- **Notebook HTML embeds resources** — `notebook/sc_analysis.qmd` uses `embed-resources: true`, so the HTML must be regenerated after source or figure changes. See F2.
- **Normalization branch comparison** — The workflow compares Seurat log1p normalization with PFlog, while the notebook currently argues strongly for PFlog. See F5 and F6.
- **Cell-cycle HVG filtering point** — `mouse_cell_cycle_genes` are removed from `VariableFeatures(sobj)` after HVG selection and before HVG scatter plotting/PCA.
- **HVG label count** — The current preprocessing script labels the top 20 retained HVGs in scatter plots, while the helper default remains `n_top = 10`. See F3.
- **PCA diagnostics and PC grid** — The notebook uses elbow plots and PCA heatmaps to motivate candidate PC dimensions, with downstream comparisons planned at 20, 30, and 50 PCs. See F4 and F6.
- **Candidate clustering grid** — Clustering uses Leiden over candidate PC dimensions and resolutions 0.3, 0.5, and 0.8.
- **Condition-level statistical unit** — Project conventions define Mouse × Condition pseudobulk sample, not individual cell, as the primary condition-level DE unit.
- **R-package-first Mycelium integration** — Mycelium adds `.living/`, manifests, and convention metadata without moving executable ESPI code out of `R/`, `scripts/`, or `notebook/`.
- **Box Drive as artifact root** — Large inputs and generated figures live under `/Users/carlstone/Library/CloudStorage/Box-Box/megan_sc_data`; notebook figures are symlinks.

## Questions for the analyst

Things the diff alone cannot tell us, whose answers would change which findings matter most.

- Is PFlog intended to be the primary normalization now, or should log1p and PFlog remain equal sensitivity branches until downstream conclusions agree?
- Should every clustering artifact be comparable across both cell-cycle branches, or is clustering only planned for one branch at a time?
- For the PCA heatmaps, do you want Seurat's default feature display, or should the figure explicitly show a fixed number of top-loading genes per PC?
- Is condition-level DEG strictly pseudobulk-by-Mouse × Condition, independent of the HVG-derived clustering branches, or will cluster definitions from these branches feed any DE summaries?
- Should rendered `notebook/sc_analysis.html` be committed/reviewed now, or should it stay out of scope until prose and figure contracts are final?

## Findings

### Statistics & causal inference

#### Major

No major findings.

#### Minor

##### F5. PFlog normalization prose overstates the evidence and includes an unsupported gene-length claim

`notebook/sc_analysis.qmd:29,58-61`
```markdown
29:1. Normalize counts: this normalizes the counts of genes by cell to account for differences in sequencing depth and technical factors. May also include normalization by gene length, so that longer genes do not appear more highly expressed simply because they have more reads.
58:  * This is also known as a shifted centered log-ratio transformation, which is 
59:  mathematically the best way to deal with compositional data on this scale.
61:**Choice 1**: log1p normalization is the default in Seurat, but it's not the best option for compositional data, and our relatively small library size exacerbates the problem. PFlog normalization is mathematically the best option for compositional data, and there is a benchmarking paper showing it outperforms log1p normalization in virtually every situation. It's also made by a group that is well-known for making read mapping and analysis pipelines bowtie, cufflinks, and kallisto. The only downside is that this exact implementation is relatively new, but the underlying math is very well-established. I think we should go with this because it stabilizes the variance better than the log1p normalization. We just may need to cite a few papers to explain why we're using it, but the single-cell GOAT says this is the best so that's good enough for me. There's kind of a single-cell methods war that's been going on for literally 2 decades, so every single option we pick in this analysis will have someone with a strong opinion on it, so I think we should just use this because it's the best.
```
**Why it matters here**: Normalization feeds HVG selection, PCA, UMAP, and clustering. The current text makes PFlog look settled by authority and broad superiority rather than justified for this ESPI dataset; the gene-length sentence also does not match the shown UMI scRNA-seq preprocessing path.
**Fix**: Remove the gene-length sentence, cite the exact PFlog benchmark/scope, and phrase PFlog as a chosen or compared branch with stated assumptions rather than as universally best.

##### F6. The PFlog elbow explanation makes a causal claim the plot does not establish

`notebook/sc_analysis.qmd:77-80`
```markdown
77:PFlog normalization has a steeper elbow, because it is better at removing noise 
78:from low-depth cells. I'm doing downstream analyses with both at
79:$\mathrm{PC} = \{20, 30, 50\}$ to see how it affects clustering and marker 
80:segregation. Those will be supplemental figures.
```
**Why it matters here**: The elbow plot can show a variance profile, but it does not identify why the profile changed. That explanation can bias the normalization decision before the planned branch comparison is interpreted.
**Fix**: Rephrase as an observed pattern, or add a targeted diagnostic linking variance/noise behavior to low-depth cells.

### Data pipeline & leakage

#### Major

##### F1. Cluster outputs can silently collide across cell-cycle branches

`scripts/cluster-sobj.R:69,84,109,124`
```r
69:    name <- sprintf("cluster_%s_dims%d_res%s", norm, d, res_tag(r))
84:  reduction_name <- sprintf("umap_%s_dims%d", norm, d)
109:    out_tag = sprintf("%s_dims%d", norm, d)
124:  sprintf("cluster_%s_elbow%d.rds", norm, cli_args$elbow_n)
```
**Why it matters here**: The notebook now frames `filter-cc` and `no-filter-cc` as complementary branches. Clustering object names, reduction names, UMAP/clustree filenames, and cluster columns carry `norm`/dims/resolution but not `filtered_cell_cycle`, so running the same normalization for both branches can overwrite or misattribute the last branch written.
**Fix**: Derive `cc_tag <- if (isTRUE(sobj@misc$preprocessing$filtered_cell_cycle)) "filter-cc" else "no-filter-cc"` and include it in clustered RDS names, UMAP reduction names, clustree/UMAP output tags, and persisted candidate cluster names.
**Behavioral check**: `report-values-freshness` / artifact-freshness tripwire should verify that report references point to the branch-specific objects and plots actually generated.

#### Minor

No minor findings beyond F7, which is the helper-level version of the same output-naming risk.

### Bioinformatics

#### Major

No major findings.

#### Minor

##### F3. HVG scatter prose says top 10 labels, but the script now labels top 20 retained HVGs

`notebook/sc_analysis.qmd:98-99; scripts/preprocess-sobj.R:86-87`
```markdown
98:We'll use the default 2000 HVGs, with and without cell cycle genes removed. 
99:These plots label the top 10 HVGs in each condition.
86:# Plot gene mean-vs-variance scatter with top 20 retained HVGs labeled.
87:splot_hvg_scatter(sobj, n_top = 20)
```
**Why it matters here**: The labeled genes are the biology readers will inspect first. Reporting the wrong label count makes the figure/method text inconsistent with the regenerated diagnostic and can confuse which retained HVGs were prioritized after cell-cycle filtering.
**Fix**: Update the notebook to say top 20 retained HVGs, or change the script back to `n_top = 10` and regenerate figures/HTML.
**Behavioral check**: `report-values-freshness` should compare report text against the plotting parameters used to generate the embedded figures.

##### F4. PCA heatmap prose misdescribes `cells = 500` as 500 genes

`notebook/sc_analysis.qmd:92; R/preprocess-plots.R:140-143`
```markdown
92:For each PC, the heatmap shows the top 10 genes that contribute to that PC. The 500 genes with the strongest expression on either side of the PC are shown. Plots where 
140:  plot <- Seurat::DimHeatmap(
142:    dims = 1:6,
143:    cells = 500,
# no `nfeatures = 10` argument is set in this call
```
**Why it matters here**: In `Seurat::DimHeatmap`, `cells` controls displayed/sampled cells, not the number of genes. The current text reverses what the diagnostic is showing and states a feature-count contract the code does not request.
**Fix**: Rewrite the text to describe PCs 1-6 and 500 cells, or set and document `nfeatures = 10` explicitly if that is the intended figure contract.
**Behavioral check**: `report-values-freshness` should catch figure prose that no longer matches plotting arguments.

### LLM coding antipatterns

#### Major

No major findings.

#### Minor

##### F9. The session log says validation passed without preserving a log or command artifact

`.living/log/LOG_REGISTRY.md:5`
```markdown
5:| 2026-07-03 | mycelium-init-audit | ESPI | not checked | current session | `.living/`; `CLAUDE.md`; `ENVIRONMENTS_INSTALLATIONS.md`; manifests; `todo/`; `.claude/settings.local.json` | Initialized Mycelium and audited every init checklist step. | Core plus bioinformatics conventions installed; `.living/INDEX.md` generated; validation passed | complete | mycelium, setup, init | — |
```
**Why it matters here**: Future agents may treat the Mycelium setup as verified, but the registry row has no durable log pointer for the validation command/output.
**Fix**: Attach the validation command/output artifact in the `Log` column, or soften the wording to a setup summary without “validation passed.”

### Documentation & schema fidelity

#### Major

##### F2. Rendered notebook HTML is stale relative to the QMD source and regenerated figures

`notebook/sc_analysis.qmd:3-10; notebook/sc_analysis.html`
```yaml
3:format: 
4:  html:
5:    toc: true
10:    embed-resources: true
# verified mtimes: sc_analysis.html 2026-07-03T14:18:31 < sc_analysis.qmd 2026-07-03T14:26:54; dim_heatmap_pflog_filter-cc target 2026-07-03T14:20:00
```
**Why it matters here**: The HTML embeds image bytes, so opening `notebook/sc_analysis.html` can show older figures/prose even though the QMD and linked Box figures changed. This is exactly the stale-report failure documented in `.living/learnings.md`.
**Fix**: Rerender `notebook/sc_analysis.qmd` after final source/figure edits, then verify the HTML is newer than the QMD and all referenced figure targets.
**Behavioral check**: `report-values-freshness` should compare rendered-report fingerprints/mtimes against source QMD and figure targets.

#### Minor

No additional minor findings after dedupe; F3, F4, F5, F6, and F9 are documentation-contract problems placed under the categories where their consequences matter most.

### Code quality

#### Major

No major findings.

#### Minor

##### F7. `splot_umap_by()` takes a reduction but omits that reduction from filenames

`R/cluster-plots.R:15-28`
```r
15:  plot <- Seurat::DimPlot(
17:    reduction = umap,
18:    group.by = color_by,
22:    file.path(out_dir, sprintf("umap_%s.png", color_by)),
28:    file.path(out_dir, sprintf("umap_%s.pdf", color_by)),
```
**Why it matters here**: The helper is exported and its contract is “plot a specified UMAP reduction.” A caller comparing the same metadata column across multiple reductions can silently overwrite the first plot; F1 is the current cluster-script manifestation of the same naming risk.
**Fix**: Include a sanitized `umap` value in the basename again, or add an explicit `out_tag` parameter that must encode the reduction and branch.

##### F8. `patchwork` is loaded by the notebook but not declared or used

`notebook/sc_analysis.qmd:13-17`
```r
13:```{r}
14:#| label: package-load
15:#| include: false
16:library(Seurat)
17:library(patchwork)
```
**Why it matters here**: A clean notebook render can depend on a transitive install rather than the repo’s declared environment. The current QMD contains no other `patchwork` references.
**Fix**: Remove `library(patchwork)` if unused, or add `patchwork` to `DESCRIPTION` and `ENVIRONMENTS_INSTALLATIONS.md` before relying on it.

## What was checked but is fine

- **Statistics & causal inference**: No changed R code introduced statistical tests, p-values, multiple-comparison procedures, or a new condition-level DE unit; the Mouse × Condition pseudobulk convention remains intact.
- **Data pipeline & leakage**: Moving `splot_hvg_scatter()` after cell-cycle HVG removal reduces diagnostic mismatch; no new joins, missing-value coercions, sample-alignment paths, or fallback data paths were introduced.
- **Bioinformatics**: The current cell-cycle HVG filter is applied to retained `VariableFeatures` before PCA, and no new gene-reference, doublet, ambient-RNA, or cell-type annotation claim was introduced in code.
- **LLM coding antipatterns**: No silent try/fallbacks, warning suppression, hallucinated load-bearing APIs, or fake portable-path claims were found; local Box/plugin-cache paths are explicitly documented as local.
- **Documentation & schema fidelity**: Roxygen and generated Rd for `splot_umap_by()` match the R source; Mycelium active convention paths, manifests, skillpack status, and figure symlinks are consistent.
- **Code quality**: The Mycelium setup respects the R-package-first layout; `.claude/` is ignored for the right reason; generated notebook figures are symlinks rather than copied binaries.

## Notes

- The first fixes I would make are F1 and F2, because they can make downstream artifacts point at the wrong branch or stale report bytes.
- F3-F6 are prose/figure-contract issues; they are fast to fix but should be done before the next Quarto render.
- F7 can be fixed together with F1 by deciding on one branch/reduction-aware output tag convention and applying it consistently.
- No new `.living/learnings.md` entry was added: the stale-HTML pattern is already recorded as L-2, and this review confirms it is still easy to trip.
