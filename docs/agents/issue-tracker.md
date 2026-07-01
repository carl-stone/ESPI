# Issue tracker: GitHub

Issues and PRDs for this repo live in GitHub Issues for `carl-stone/ESPI`. Use the `gh` CLI for all issue operations.

## Conventions

- **Create an issue**: `gh issue create --title "..." --body "..."`. Use a heredoc for multi-line bodies.
- **Read an issue**: `gh issue view <number> --comments`, including labels and relevant comments.
- **List issues**: `gh issue list --state open --json number,title,body,labels,comments` with appropriate `--label` and `--state` filters.
- **Comment on an issue**: `gh issue comment <number> --body "..."`.
- **Apply / remove labels**: `gh issue edit <number> --add-label "..."` / `--remove-label "..."`.
- **Close an issue**: `gh issue close <number> --comment "..."`.

Infer the repo from `git remote -v`; `gh` does this automatically when run inside this clone.

## Pull requests as a triage surface

External PRs are **not** a request surface for this repo. `/triage` should process GitHub Issues, not pull requests.

GitHub shares one number space across issues and PRs. If a referenced number is ambiguous, check whether it is an issue before treating it as a triage ticket.

## When a skill says "publish to the issue tracker"

Create a GitHub issue in `carl-stone/ESPI`.

## When a skill says "fetch the relevant ticket"

Run `gh issue view <number> --comments` in this repo.
