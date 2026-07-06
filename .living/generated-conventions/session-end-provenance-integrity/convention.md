---
id: session-end-provenance-integrity
title: Verify semantic Mycelium session-end records before yielding or committing.
status: active
created: 2026-07-06
source_learnings:
  - L-10
  - L-20
  - L-21
  - L-22
  - L-23
description: Prevent hook-generated file-list placeholders and stale stop-hook state from replacing the semantic audit trail.
---

## Statement

Before yielding or committing after Mycelium-tracked work, verify that the current session records are semantic and grounded in the current Git-visible state:

1. `git status --short` is the authoritative changed-file set.
2. `.living/log/LOG_REGISTRY.md` has a completed row for the current session with a sentence Summary, non-empty Key Outputs, useful Tags, and no filename-stub summary.
3. The linked `.living/log/<session>.md` has `## Session Summary`, `## Key Outputs`, and `## Status` sections that describe the work, not just files touched.
4. `.claude/last-session.md` covers the full session and reports the current commit/uncommitted state.
5. `.living/INDEX.md` has been regenerated after any `.living/` decision, learning, convention, finding, or log-registry change.
6. Any hook-created `.log-scribe-*` authentication-failure logs or false-positive session logs are removed before yielding.

## Rationale

This project repeatedly saw Mycelium hook output overwrite or dilute the audit trail: completed registry rows with blank semantic fields, file-list-only summaries appended after validation, stale stop-hook reminders from prior work, and local log-scribe authentication failures. These failures make future resume, review, and scientific audit harder because the repo appears complete while losing the actual decisions, outputs, and validation evidence.

**Source learnings:**
- L-10: completed `LOG_REGISTRY` rows still need semantic fields.
- L-20: Mycelium hook summaries can overwrite manual semantic log rows.
- L-21: hook-provenance guard tests must persist with the guard.
- L-22: review sessions can still expose provenance clobbering.
- L-23: `git status --short` should anchor todo-only stop-hook triage.

## Correct Application

```text
1. Finish implementation and verification.
2. Run git status --short.
3. Repair current LOG_REGISTRY row and linked session log if they contain file-list stubs or blank semantic fields.
4. Regenerate .living/INDEX.md.
5. Update .claude/last-session.md.
6. Re-run git status --short and yield/commit only after the records match the actual changed-file set.
```

This is correct because it treats hooks as useful assistants, not as the source of truth. The changed-file set and verified outputs drive the provenance record.

## Incorrect Application

```text
- Accept a completed LOG_REGISTRY row with blank Key Outputs.
- Leave a current session log that only says "Modified: file1,file2 (+N more)".
- Report a clean or complete state without re-checking after hooks update .living/ or .claude/ files.
- Backdate a new session's work into a closed prior session log.
```

These patterns make the audit trail look complete while hiding what actually changed and why.

## Exceptions

- **Read-only Q&A with no Git-visible changes**: Do not manufacture `.living/` entries just to appease a false-positive hook. Remove false-positive hook artifacts if they appear.
- **Trivial local-only `.claude/` updates**: Keep them out of Git unless the project explicitly tracks `.claude/` files.

## Review

This convention is active for ESPI because the same failure mode recurred across five learnings. If the hook/provenance guard later enforces every checklist item structurally, downgrade the manual checklist to a brief verification note.
