# Domain Docs

How the engineering skills should consume this repo's domain documentation when exploring the codebase.

## Layout

This repo uses a **single-context** domain-doc layout.

Expected locations:

- `CONTEXT.md` at the repo root for the project glossary and domain language.
- `docs/adr/` for architectural decision records that apply to this repo.

## Before exploring, read these

- Read root `CONTEXT.md` when it exists.
- Read relevant ADRs under `docs/adr/` when they touch the area you're about to change.

If these files don't exist, proceed silently. Don't flag their absence or suggest creating them upfront. The domain-modeling workflow creates them lazily when terms or decisions get resolved.

## File structure

```text
/
├── CONTEXT.md
└── docs/adr/
    ├── 0001-example-decision.md
    └── 0002-example-decision.md
```

## Use the glossary's vocabulary

When your output names a domain concept, use the term as defined in `CONTEXT.md`. Don't drift to synonyms the glossary explicitly avoids.

If the concept you need isn't in the glossary yet, either reconsider the term or note the gap for domain modeling.

## Flag ADR conflicts

If your output contradicts an existing ADR, surface the conflict explicitly rather than silently overriding it.
