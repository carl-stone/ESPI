# Origin

## Source Learnings

- **2026-07-04 — L-10**: "Completed LOG_REGISTRY rows still need semantic fields" — established that completed registry rows can be structurally present but semantically empty.
- **2026-07-05 — L-20**: "Mycelium hook summaries can overwrite manual semantic log rows" — showed that hook-generated file-list placeholders can replace manually repaired summaries.
- **2026-07-05 — L-21**: "Persist hook-provenance guard tests with the guard" — identified persisted tests as the structural mitigation for provenance guards.
- **2026-07-05 — L-22**: "Review sessions can still expose Mycelium provenance clobbering" — showed that review and stop-hook activity can still clobber registry/log/last-session records after validation.
- **2026-07-05 — L-23**: "Treat git status as authority for todo-only stop-hook triage" — added the rule that current `git status --short` should anchor session-end provenance rather than stale hook state.

## Source Decisions

- **2026-07-04**: "Filter Mycelium maintenance commands before post-action hooks" — created a wrapper to avoid maintenance commands starting new post-action cycles.
- **2026-07-04**: "Mirror complete Mycelium hook behavior in OMP" — made the OMP adapter return hook context and call data-lineage hooks.
- **2026-07-04**: "Reset prior-session Mycelium sentinels at SessionStart" — addressed stale sentinel bleed across sessions.

## Pattern

Five recent learnings share `mycelium`, `hooks`, `session-logs`, or `provenance` tags and describe the same operational risk: automated session hooks can create, overwrite, or preserve syntactically valid but semantically weak provenance records. The consequences recur at session boundaries, after validation, after commits, and in todo-only work.

This warrants a convention rather than another ad hoc reminder because several entries remain `ambient-awareness` and depend on the agent remembering to inspect the records. The convention turns that memory into a checklist and keeps the structural guard tests linked to the process.

## Contributing Projects

- ESPI: 5 source learnings.
