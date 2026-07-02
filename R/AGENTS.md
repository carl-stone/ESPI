# AGENTS.md - R/

Keep helper functions narrow, but do not wrap a few commands in a new function unless the helper is called repeatedly or names a real analysis concept.

Prefer inline code for one-off operations. A helper that only replaces one or two obvious statements makes the analysis harder to read by hiding the work behind another name.
