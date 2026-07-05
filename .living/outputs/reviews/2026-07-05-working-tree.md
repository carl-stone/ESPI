# Review — Working tree — 2026-07-05

**Scope**: Working tree against `HEAD`, including untracked files  
**Files reviewed**: 32  
**Sub-agents run**: 6 — stats-causal, data-pipeline-leakage, bioinformatics, llm-failure-modes, doc-schema-fidelity, code-quality

## Key decisions in this analysis

The consequential analytical choices in this work. Some have associated findings below; others are informational so you can decide whether to revisit them.

- **Metadata labels vs display labels** — Computational condition labels stay `p27CKO +EStim`/`p27CKO`; human-facing plots use `p27CKO + E-Stim` through display constants.
- **Contrast direction** — Plot labels now say `p27CKO + E-Stim vs. p27CKO`, while DESeq2 and muscat still use `estim_vs_control` with positive effects higher in E-Stim.
- **Palette source of truth** — Direction/neutral colors now come from `palette_analysis_three`; `palette_dotplot_pair` derives from that palette.
- **Cluster-abundance file split** — `R/cluster-abundance.R` retains compute/randomization code; `R/cluster-abundance-plots.R` owns plotting helpers.
- **Abundance interpretation boundary** — Pooled Fisher/CLR abundance remains descriptive; Mouse × Condition sample-level randomization remains the inferential abundance screen.
- **Stable cross-reference policy** — New records use stable IDs, but historical `Figure N` references are preserved by current decision. See F2.
- **Living-record audit trail** — The current tree records regenerated figures, DE counts, notebook render, and tripwires as completion evidence. See F1 and F3.

## Questions for the analyst

Things the diff alone cannot settle, whose answers change which cleanup matters most.

- Should historical `Figure N` references be preserved as audit history, or should the original TODO remain open until a targeted stable-ID cleanup is done?
- Do you want figure legends to spell out `p27CKO + E-Stim` everywhere, or is shorthand `E-Stim` acceptable once the axis names the full contrast?
- Should Mycelium hook-generated session logs be commit artifacts, or should in-progress/failed hook logs be cleaned before committing?
- Is the regenerated DE/enrichment output intended to be part of this commit, or only used as verification that labels and report values remain fresh?

## Findings

### Statistics & causal inference

#### Major

No major findings.

#### Minor

##### F1. Abundance legend still uses shorthand condition labels

`R/cluster-abundance-plots.R:165-168`
```r
  plot_data$direction <- factor(
    plot_data$direction,
    levels = c("Enriched in E-Stim", "Depleted in E-Stim", "Not significant")
  )
```
**Why it matters here**: The same abundance figure now has a full contrast axis, `p27CKO + E-Stim vs. p27CKO`, but a shorthand `E-Stim` legend. The calculation is unaffected, but the display-label cleanup is incomplete and could confuse readers about whether `E-Stim` means the full `p27CKO + E-Stim` condition.  
**Fix**: Derive direction labels from `ESTIM_DISPLAY_LABEL`/`CTRL_DISPLAY_LABEL`, or centralize direction labels beside the condition display constants.

**Resolution status**: Resolved by adding `direction_labels` derived from `ESTIM_DISPLAY_LABEL` in `plot_clr_fisher_enrichment()`, regenerating the MG-selected abundance figure, rerendering the notebook, and rerunning tripwires.

### Data pipeline & leakage

#### Major

No major findings.

#### Minor

No minor findings. The reviewers checked condition-label integrity, output freshness, gene joins, Mouse × Condition unit handling, and palette source-of-truth drift.

### Bioinformatics

#### Major

No major findings.

#### Minor

No minor findings. The reviewers confirmed the metadata labels still match `analysis_labels.yml`, DE/DD still use Mouse × Condition pseudobulk samples, and the notebook still frames pooled Fisher/CLR abundance as descriptive.

### LLM coding antipatterns

#### Major

##### F2. Session log claims verification after the recorded session ended

`.living/log/2026-07-05-003-espi.md:6-31`
```markdown
ended: 2026-07-05T14:31:40-0500
### 14:36 — Regenerated outputs and verified
- Ran `Rscript scripts/plot-mg-selected-figures.R`; regenerated MG-selected UMAP, feature UMAP, abundance, and cluster-proportion outputs.
- Checked DE row counts after rerun: 453 significant DEGs and 24,514 primary DE genes, unchanged from notebook prose.
- Ran `Rscript tools/run-tripwires.R`; all non-skipped checks passed (`label-permutation` skipped by design).
```
**Why it matters here**: This log is the audit trail for the regenerated figures, DE row-count check, notebook render, and tripwires. The underlying commands did run in this session, but the committed record is internally impossible, so a reader cannot trust the chronology without re-checking external transcripts.  
**Fix**: Correct the session end time/duration and remove the stale `14:31 — Session ended` block, or regenerate the session log from the verified command outputs.

**Resolution status**: Resolved by correcting `.living/log/2026-07-05-003-espi.md` to end after the 14:36 verification block and removing the stale generated end/file-list section.

#### Minor

No minor LLM-specific findings after deduplication; the TODO-status issue is reported under Documentation & schema fidelity.

### Documentation & schema fidelity

#### Major

No major findings beyond F2.

#### Minor

##### F3. Stable-reference TODO was closed by changing the criterion

`todo/stable-cross-references.md:27-33`
```markdown
- [x] Document the project convention for stable cross-references.
- [x] Record the historical-fidelity exception: current `Figure N` references in historical review/session records are intentionally left intact unless the record is otherwise being edited.
- [x] Future notebook/log prose uses stable IDs or captions for cross-referenced artifacts.
Completed in Batch 1 presentation cleanup. Added the forward-going Cross-References convention and recorded the historical-record exception; the original review/log evidence was left intact by design.
```
**Why it matters here**: The original TODO asked to replace current fragile `Figure N` references where stable IDs are available. The current file instead marks a narrower historical-exception policy complete, which hides that the original cleanup was consciously superseded rather than executed.  
**Fix**: Restore the original criterion as open/deferred or explicitly mark it superseded by the stable-reference decision; close only the convention/historical-exception work that actually happened.

**Resolution status**: Resolved by restoring the original historical-rewrite criterion as explicitly superseded by the `Require stable cross-references going forward` decision instead of silently marking it performed.

##### F4. Completed LOG_REGISTRY row still has placeholder fields

`.living/log/LOG_REGISTRY.md:25`
```markdown
| 2026-07-05 | 2026-07-05-003 | espi | main | 7m | 22 | conventions.md,decisions.md,INDEX.md (+16 more) | | complete | | [log](2026-07-05-003-espi.md) |
```
**Why it matters here**: The registry is the searchable schema for Mycelium sessions. A row marked complete with a file-list summary, blank Key Outputs, and blank Tags makes this presentation-cleanup session hard to find and contradicts the actual work described in the log.  
**Fix**: Fill Summary, Key Outputs, and Tags with semantic Batch 1 content, and make duration/files changed match the corrected session log.

**Resolution status**: Resolved for row `2026-07-05-003` by replacing placeholder fields with semantic Summary, Key Outputs, Tags, and matching duration/files-changed values.

### Code quality

#### Major

No major findings.

#### Minor

##### F5. Hidden log-scribe authentication failure is untracked

`.living/log/.log-scribe-2026-07-05-003.log:1`
```text
Failed to authenticate. API Error: 401 Invalid authentication credentials
```
**Why it matters here**: This is transient tool output, not durable provenance. Committing it would add noisy external-service failure detail under `.living/log/` without improving reproducibility.  
**Fix**: Delete the scratch log before commit, or add the log-scribe scratch pattern to `.gitignore` if these files are generated routinely.

**Resolution status**: Resolved by removing `.living/log/.log-scribe-2026-07-05-003.log` from the working tree.

##### F6. In-progress session log is untracked with blank metadata

`.living/log/2026-07-05-004-espi.md:5-8`
```yaml
started: 2026-07-05T14:31:41-0500
ended:
duration_minutes:
files_changed:
```
**Why it matters here**: If committed as-is, this creates an unfinished living-record artifact that confuses audit and resume tooling. It is probably the current review session stub, but it should not be committed incomplete.  
**Fix**: Finalize this log before commit, or remove it if the review report and `last-session.md` are the intended records.

**Resolution status**: Resolved after the review by finalizing `.living/log/2026-07-05-004-espi.md`; retained here as a point-in-time finding.

## What was checked but is fine

- **Statistics & causal inference**: DESeq2 contrast direction remains target-over-reference, multiple-comparison handling is unchanged, and pooled Fisher/CLR is still explicitly descriptive.
- **Data pipeline & leakage**: The display-label split does not mutate metadata labels; regenerated DE counts still match notebook claims; DE/DD scatter still discloses its inner gene join.
- **Bioinformatics**: Mouse × Condition remains the condition-level unit, `analysis_labels.yml` still matches package metadata labels, and no gene-symbol or scRNA-seq preprocessing drift was introduced.
- **LLM coding antipatterns**: Removed palette fallbacks now fail fast; new symbols are exported; no hallucinated ESPI constants were found.
- **Documentation & schema fidelity**: `NAMESPACE`/Rd files match the new exports; notebook abundance prose and plot source agree on the full contrast label.
- **Code quality**: The plot-helper split left no duplicate definitions, generated Rd/HTML files are expected tracked artifacts, and the remaining `#4daf4a` is an intentional fourth DE/DD category outside the low/mid/high palette.

## Notes

The static review found one fix-before-commit issue: the Mycelium session log/registry audit trail is currently stale relative to the actual verified commands. The scientific computations and regenerated report values checked out; the problems are presentation consistency and provenance hygiene.
