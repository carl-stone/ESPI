# Review — working tree — 2026-07-14

**Scope**: All tracked and untracked working-tree changes since `HEAD`
**Files reviewed**: 24
**Sub-agents run**: 6

## Key decisions in this analysis

- **Remove differential detection** — Retire the muscat detection-fraction estimand, dependency, outputs, notebook text, and combined DE/DD scatter while preserving historical decision records.
- **Keep DESeq2 as the condition-level analysis** — Use Mouse × Condition pseudobulk samples with the six-sample `~ condition` model as primary and the four-sample `~ mouse + condition` model as paired sensitivity.
- **Use a primary-model-only volcano** — Plot shrunken `log2FoldChange` against `-log10(padj)` from `full_de`; retain paired results in tables and prose rather than a second panel.
- **Define significance only by FDR** — Color genes at `padj < 0.05` without an absolute fold-change cutoff; clamp zero adjusted P values only for plotting.
- **Label the top ten significant genes** — Rank by adjusted P value, then absolute shrunken effect, then gene name; force all ten selected labels to render.
- **Drop unplottable rows** — Exclude missing/empty genes and non-finite shrunken effects or adjusted P values from the volcano without modifying stored DE tables.

## Questions for the analyst

- Was differential detection ever preregistered, promised to collaborators, or included in a submitted manuscript? If so, should its observed null result remain in a supplement even though the active pipeline no longer computes it?
- Is the volcano intended for a manuscript figure, exploratory review, or target selection? That determines whether the caption should report the plotted-gene count and omitted non-estimable rows.
- Should the ten labels remain purely rank-based, or will downstream readers expect biologically curated labels or explicit validation priorities?
- Is the paired analysis only a direction-of-effect sensitivity check, or should agreement with the primary model become a formal reporting criterion?

## Findings

### Statistics & causal inference

No findings.

### Data pipeline & leakage

No findings.

### Bioinformatics

No findings.

### LLM coding antipatterns

No findings.

### Documentation & schema fidelity

#### Minor

##### F1. Repair malformed Mycelium session provenance

`.living/log/LOG_REGISTRY.md:64-66`
```markdown
| 2026-07-14 | 2026-07-14-001 | espi | main | 51m | 1 |  | | complete | | [log](2026-07-14-001-espi.md) |
| 2026-07-14 | 2026-07-14-002 | espi | main | 5m | 2 | remove-differential-detection-plan.md,remove-differential-detection-plan.md | | complete | | [log](2026-07-14-002-espi.md) |
| 2026-07-14 | 2026-07-14-003 | espi | main | 35m | 16 | Removed differential detection and replaced the combined scatter with a primary DE volcano. | Canonical counts-QC pipeline, focused DE regeneration, notebook render, and tripwires passed; ten volcano labels verified | complete | differential-expression, mg-selected, plotting, reproducibility | [log](2026-07-14-003-espi.md) |
```

**Why it matters here**: Sessions 001 and 002 are marked complete without semantic summaries, outputs, or tags. Session 003 conflicts with its linked log (`37m`, `18` files), and that log records several incorrect root-level paths; the untracked `.log-scribe-2026-07-14-003.log` is only a 401 authentication failure. These records weaken the audit trail for a scientifically consequential output cutover.

**Fix**: Remove false-positive sessions 001/002 and the authentication artifact, then repair session 003 metadata, semantic sections, and file paths before regenerating `.living/INDEX.md`.

### Code quality

No additional findings.

## What was checked but is fine

- **Statistics & causal inference**: The pseudobulk unit, all-sample primary model, paired-mouse sensitivity, BH threshold, and no-fold-change-cutoff policy are explicit; removing DD is surfaced above as a consequential scope decision rather than treated as a defect without evidence of a prespecified obligation.
- **Data pipeline & leakage**: Volcano filtering, zero-P plotting clamp, deterministic labels, sample identity, output protection, and current external artifact cleanup are internally consistent.
- **Bioinformatics**: DESeq2 replicate handling, shrinkage semantics, gene ranking, FDR interpretation, and notebook claims match the implementation.
- **LLM coding antipatterns**: No fabricated analysis output, hidden model fallback, hallucinated API, or load-bearing undocumented default was established. Excluding non-estimable rows is an explicit plotting contract.
- **Documentation & schema fidelity**: Active notebook source/render, manifests, dependency declarations, volcano names, thresholds, and tripwire expectations agree with code and generated outputs.
- **Code quality**: The volcano helpers are cohesive, old DD paths are absent from active contracts, and historical DD decisions remain clearly superseded rather than acting as duplicate configuration.

## Notes

- The completed null DD endpoint raises a selective-reporting question only if it was prespecified or previously committed as a reporting endpoint; the working-tree diff alone does not establish that condition.
- Reviewers considered stale external DD artifacts, but the canonical Box directories and obsolete export bundle were already cleaned during this cutover, so no current output-lifecycle defect was retained.
