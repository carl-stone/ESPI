# Origin

## Source Learnings

- **2026-07-16 — L-47, "Separate source contracts from visible sensitivity oracles"**: showed that scientific input ownership and report-visible output ownership can differ, requiring independent source and visible-artifact oracles.
- **2026-07-16 — L-48, "Publication cutovers need layered equivalence checks"**: showed that byte equality and visual inspection alone cannot distinguish serialization drift from scientific drift.
- **2026-07-16 — L-49, "Review maintenance paths separately from routine runs"**: showed that routine publication parity can hide broken frozen regeneration, undeclared dependencies, and destructive artifact replacement.
- **2026-07-17 — L-51, "Exercise frozen regeneration before relying on routine phases"**: showed that deliberate phase-01 execution exposes failures that downstream routine phases do not exercise.

## Source Decisions

- **2026-07-16 — D-59, "Freeze all existing Seurat objects and designate two final objects"**: fixed the source and MG-selected objects, their roles, checksums, and protected directories as explicit contracts.
- **2026-07-16 — D-60, "Define boundaries for the publication-pipeline restructuring proposal"**: required phase-by-phase equivalence proof, preserved artifacts, and retention of only checks that protect scientific meaning.
- **2026-07-16 — D-61, "Adopt the clean four-phase publication cutover"**: established fixed phase interfaces, a deliberate frozen-stage regeneration path, and safe regular-file notebook mirrors.
- **2026-07-16 — D-62, "Resolve the source-versus-visible MG heatmap oracle conflict"**: required both MG heatmap branches and the ordered visible notebook artifacts as acceptance oracles.

## Pattern

Four recent learnings and four supporting decisions share `pipeline`, `reproducibility`, `notebook`, `validation`, or `maintenance` concerns. Together they define a recurring publication-pipeline integrity pattern: source contracts and report-visible artifacts need separate ownership; equivalence needs scientific layers; maintenance needs independent execution; and replacement must be verified without destroying the prior artifact.

The pattern meets the threshold because it recurs across cutover design, figure and notebook parity, frozen-object protection, and deliberate regeneration. It is transferable to scientific publication pipelines that replace scripts or rebuild report artifacts, while the ESPI object names and commands remain project-specific examples.

## Contributing Projects

- ESPI: four source learnings and four supporting decisions.
