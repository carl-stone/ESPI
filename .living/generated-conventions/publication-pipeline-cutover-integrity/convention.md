---
id: publication-pipeline-cutover-integrity
title: Prove publication-pipeline cutovers with independent artifact oracles and safe replacement.
status: proposed
created: 2026-07-25
source_learnings:
  - L-47
  - L-48
  - L-49
  - L-51
source_decisions:
  - D-59
  - D-60
  - D-61
  - D-62
description: Preserve scientific meaning and report-visible artifacts when replacing a publication pipeline.
---

## Statement

For a publication-pipeline cutover or maintenance change:

1. Define separate oracles for source-contract outputs and report-visible artifacts. Do not use a source object as the only oracle for visible sensitivity outputs.
2. Apply layered scientific equivalence checks: verify fixed inputs; compare canonicalized tables; compare decoded figure pixels and dimensions; normalize only documented volatile report identifiers; and compare the rendered report at a fixed viewport.
3. Execute deliberate maintenance paths separately from routine runs. Exercise frozen-input regeneration and declare each direct optional dependency, output-directory requirement, and replacement authorization.
4. Replace artifacts non-destructively. Create a temporary regular-file replacement, verify its hash, dimensions, and content, then atomically replace the prior artifact. Preserve the prior regular file until verification succeeds.
5. Record the allowed normalizations and the evidence for each oracle before declaring the cutover complete.

## Rationale

Byte equality can report harmless path, ordering, identifier, or encoding changes as scientific drift, while visual review can miss table and contract changes. Routine downstream success can also hide a broken frozen-regeneration path. Independent source and visible-artifact ownership, layered checks, separate maintenance execution, and verified replacement protect both scientific meaning and the published report.

## Correct Application

In ESPI, retain the frozen source, MG-selected, and sensitivity objects; use the source contracts for analysis outputs and the ordered notebook figures for visible artifacts. Compare canonical tables, decoded figures, dimensions, normalized report identifiers, and a fixed-viewport render. Run `just regenerate-frozen` as its own maintenance review, then mirror each notebook figure through a temporary regular file and replace it only after hash and dimension checks pass.

The same pattern applies to another scientific pipeline when a new script set replaces an established analysis: define input/output contracts and report oracles independently, test the rarely used rebuild path, and keep the old publication artifact until the verified replacement is ready.

## Incorrect Application

```text
- Declare equivalence from byte-for-byte output comparison alone.
- Use the primary analysis object as the oracle for every visible sensitivity figure.
- Treat a successful routine run as proof that frozen-input regeneration works.
- Delete the existing report artifact before verifying its replacement.
- Allow undocumented path, ordering, DOM, or encoding normalization.
```

These shortcuts can hide scientific drift, omit visible results, or turn a transient replacement failure into data loss.

## Exceptions

- **Documented serialization differences**: Allow path, ordering, identifier, or encoding differences only when the canonicalization is defined before comparison and does not change scientific content.
- **New or intentionally changed science**: Update the source contract and report oracle, record the intended change, and do not call it an equivalence-preserving cutover.
- **Read-only or protected frozen inputs**: Do not make them writable for routine runs; use the explicit, authorized maintenance path only.

## Review

This convention is proposed for ESPI. Reassess it after the publication pipeline, frozen regeneration path, and report-mirroring checks are maintained by structural safeguards; then shorten the manual checklist only if those safeguards preserve all four oracle and replacement requirements.
