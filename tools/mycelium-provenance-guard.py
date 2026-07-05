#!/usr/bin/env python3
"""Repair Mycelium hook placeholder provenance without editing synced core.

Synced Mycelium hooks may upsert a placeholder LOG_REGISTRY row, rewrite
.claude/last-session.md with a generic resume stub, or regenerate the heuristic
INDEX.md summary in a way that drops recent semantic ordering. This guard
compares pre-hook snapshots with post-hook files and restores prior semantic
content that was clobbered by those placeholders.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable


REGISTRY_REL = Path(".living/log/LOG_REGISTRY.md")
INDEX_REL = Path(".living/INDEX.md")
LAST_SESSION_REL = Path(".claude/last-session.md")

REGISTRY_COLUMNS = {
    "session": 1,
    "summary": 6,
    "key_outputs": 7,
    "tags": 9,
}

GENERIC_DECISIONS = "- See `.living/decisions.md` for decisions logged during the session."
GENERIC_BLOCKERS = "- No blocker was recorded by the OMP Mycelium hook."
GENERIC_NEXT_STEPS = "- Continue from the current user request and rerun targeted validation after further changes."


def read_text_if_present(path: Path) -> str | None:
    if not path.is_file():
        return None
    return path.read_text(encoding="utf-8")


def split_registry_row(line: str) -> list[str] | None:
    stripped = line.strip()
    if not stripped.startswith("|") or not stripped.endswith("|"):
        return None
    cells = [cell.strip() for cell in stripped.strip("|").split("|")]
    if len(cells) <= max(REGISTRY_COLUMNS.values()):
        return None
    if all(set(cell) <= {"-", ":", " "} for cell in cells):
        return None
    return cells


def registry_session_id(line: str) -> str | None:
    cells = split_registry_row(line)
    if cells is None:
        return None
    session_id = cells[REGISTRY_COLUMNS["session"]].strip()
    return session_id or None


def nonempty_cell(text: str) -> bool:
    stripped = text.strip()
    return stripped not in {"", "—", "-", "_"}


def looks_like_file_list(text: str) -> bool:
    stripped = text.strip().strip("`")
    if not stripped:
        return True
    lowered = stripped.lower()
    if lowered.startswith("modified ") or lowered.startswith("session: "):
        return True
    tokens = [part.strip().strip("`") for part in stripped.replace(";", ",").split(",")]
    tokens = [token for token in tokens if token]
    if not tokens:
        return True
    fileish = 0
    for token in tokens:
        if "(+" in token and "more" in token:
            fileish += 1
        elif "/" in token or token.endswith(
            (".R", ".r", ".md", ".qmd", ".py", ".sh", ".ts", ".json", ".yml", ".yaml")
        ):
            fileish += 1
    return fileish == len(tokens)


def registry_row_is_semantic(line: str) -> bool:
    cells = split_registry_row(line)
    if cells is None:
        return False
    summary = cells[REGISTRY_COLUMNS["summary"]]
    key_outputs = cells[REGISTRY_COLUMNS["key_outputs"]]
    tags = cells[REGISTRY_COLUMNS["tags"]]
    return nonempty_cell(summary) and not looks_like_file_list(summary) and (
        nonempty_cell(key_outputs) or nonempty_cell(tags)
    )


def registry_row_lost_semantics(line: str) -> bool:
    cells = split_registry_row(line)
    if cells is None:
        return False
    summary = cells[REGISTRY_COLUMNS["summary"]]
    key_outputs = cells[REGISTRY_COLUMNS["key_outputs"]]
    tags = cells[REGISTRY_COLUMNS["tags"]]
    return (
        not nonempty_cell(summary)
        or looks_like_file_list(summary)
        or not nonempty_cell(key_outputs)
        or not nonempty_cell(tags)
    )


def registry_rows_by_session(lines: Iterable[str]) -> dict[str, str]:
    rows: dict[str, str] = {}
    for line in lines:
        session_id = registry_session_id(line)
        if session_id:
            rows[session_id] = line
    return rows


def repair_registry(before_text: str | None, after_path: Path) -> bool:
    after_text = read_text_if_present(after_path)
    if before_text is None or after_text is None:
        return False

    before_rows = registry_rows_by_session(before_text.splitlines())
    after_lines = after_text.splitlines()
    changed = False
    repaired: list[str] = []

    for line in after_lines:
        session_id = registry_session_id(line)
        before_row = before_rows.get(session_id or "")
        if (
            before_row is not None
            and registry_row_is_semantic(before_row)
            and registry_row_lost_semantics(line)
        ):
            repaired.append(before_row)
            changed = True
        else:
            repaired.append(line)

    if not changed:
        return False

    trailing_newline = "\n" if after_text.endswith("\n") else ""
    after_path.write_text("\n".join(repaired) + trailing_newline, encoding="utf-8")
    return True

def is_generated_last_session_stub(text: str) -> bool:
    if text.lstrip().startswith("# Session resume"):
        return GENERIC_DECISIONS in text and GENERIC_BLOCKERS in text

    if not text.lstrip().startswith("## What was worked on"):
        return False
    if "## Current state" not in text:
        return False
    if "## Key decisions made" in text or "## Blockers & surprises" in text:
        return False
    return all(is_generic_or_file_only_line(line) for line in text.splitlines())


def is_generic_or_file_only_line(line: str) -> bool:
    stripped = line.strip()
    if stripped in {"", GENERIC_DECISIONS, GENERIC_BLOCKERS, GENERIC_NEXT_STEPS}:
        return True
    if stripped.startswith("#"):
        return True
    if stripped.startswith("- Modified `") and stripped.endswith("`"):
        return True
    if stripped.startswith("- `") and stripped.endswith("`"):
        return True
    if stripped.startswith("- Session: "):
        return True
    if stripped == "- No Edit/Write activity was tracked.":
        return True
    if stripped.startswith("- Branch: "):
        return True
    return False


def last_session_is_semantic(text: str | None) -> bool:
    if text is None or not text.strip():
        return False
    return any(not is_generic_or_file_only_line(line) for line in text.splitlines())


def repair_last_session(before_text: str | None, after_path: Path) -> bool:
    after_text = read_text_if_present(after_path)
    if before_text is None or after_text is None:
        return False
    if not last_session_is_semantic(before_text):
        return False
    if not is_generated_last_session_stub(after_text):
        return False
    if before_text == after_text:
        return False
    after_path.write_text(before_text, encoding="utf-8")
    return True


def markdown_block(text: str, begin: str, end: str) -> str | None:
    start = text.find(begin)
    if start < 0:
        return None
    stop = text.find(end, start + len(begin))
    if stop < 0:
        return None
    return text[start : stop + len(end)]


def replace_markdown_block(text: str, begin: str, end: str, block: str) -> str:
    start = text.find(begin)
    stop = text.find(end, start + len(begin))
    if start < 0 or stop < 0:
        return text
    return text[:start] + block + text[stop + len(end) :]


def most_recent_ids(block: str) -> list[str]:
    ids: list[str] = []
    in_recent = False
    for line in block.splitlines():
        if line.strip() == "## Most recent (10)":
            in_recent = True
            continue
        if in_recent and line.startswith("## "):
            break
        if in_recent and line.startswith("- ["):
            parts = line.split("] ", 1)
            if len(parts) == 2:
                ids.append(parts[1].split(":", 1)[0].strip())
    return ids


def index_summary_regressed(before_block: str, after_block: str) -> bool:
    before_ids = most_recent_ids(before_block)
    after_ids = most_recent_ids(after_block)
    if not before_ids or not after_ids:
        return False
    common_before = [entry_id for entry_id in before_ids if entry_id in after_ids]
    common_after = [entry_id for entry_id in after_ids if entry_id in before_ids]
    return common_before != common_after


def repair_index(before_text: str | None, after_path: Path) -> bool:
    after_text = read_text_if_present(after_path)
    if before_text is None or after_text is None:
        return False
    begin = "<!-- BEGIN KNOWLEDGE SUMMARY -->"
    end = "<!-- END KNOWLEDGE SUMMARY -->"
    before_block = markdown_block(before_text, begin, end)
    after_block = markdown_block(after_text, begin, end)
    if before_block is None or after_block is None or before_block == after_block:
        return False
    if not index_summary_regressed(before_block, after_block):
        return False
    after_path.write_text(
        replace_markdown_block(after_text, begin, end, before_block),
        encoding="utf-8",
    )
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--registry-before", type=Path, required=True)
    parser.add_argument("--index-before", type=Path, required=True)
    parser.add_argument("--last-session-before", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    before_registry = read_text_if_present(args.registry_before)
    before_index = read_text_if_present(args.index_before)
    before_last_session = read_text_if_present(args.last_session_before)

    repair_registry(before_registry, repo_root / REGISTRY_REL)
    repair_index(before_index, repo_root / INDEX_REL)
    repair_last_session(before_last_session, repo_root / LAST_SESSION_REL)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
