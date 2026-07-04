#!/usr/bin/env python3
"""Sync repo-local skills/core from the OMP-installed Mycelium plugin.

Default behavior is intentionally quiet: no output when the copy is already
current, one line when files are updated, and a non-zero exit on errors.
"""

from __future__ import annotations

import argparse
import fnmatch
import hashlib
import os
import shutil
import sys
from pathlib import Path

DEFAULT_SOURCE_CANDIDATES = (
    "~/.omp/plugins/cache/marketplaces/mycelium/skills/core",
    "~/.omp/plugins/cache/plugins/mycelium___mycelium___0.0.0/skills/core",
)
EXCLUDE_NAMES = {".DS_Store", "__pycache__", ".pytest_cache"}
EXCLUDE_PATTERNS = {"*.pyc", "*.pyo"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sync skills/core from the OMP Mycelium plugin cache."
    )
    parser.add_argument(
        "--source",
        default=os.environ.get("MYCELIUM_SKILLS_CORE_SOURCE"),
        help="Source skills/core path. Defaults to known OMP Mycelium plugin cache paths.",
    )
    parser.add_argument(
        "--dest",
        default="skills/core",
        help="Destination skills/core path, relative to the current repo by default.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Only check whether sync is needed; do not update files.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress the normal updated/already-current message.",
    )
    return parser.parse_args()


def is_excluded(path: Path) -> bool:
    return any(part in EXCLUDE_NAMES for part in path.parts) or any(
        fnmatch.fnmatch(path.name, pattern) for pattern in EXCLUDE_PATTERNS
    )


def resolve_source(source_arg: str | None) -> Path:
    candidates = [source_arg] if source_arg else list(DEFAULT_SOURCE_CANDIDATES)
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate).expanduser().resolve()
        if path.is_dir():
            return path
    searched = ", ".join(str(Path(p).expanduser()) for p in candidates if p)
    raise FileNotFoundError(f"Could not find OMP Mycelium skills/core. Searched: {searched}")


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_files(root: Path) -> dict[Path, Path]:
    files: dict[Path, Path] = {}
    for path in root.rglob("*"):
        rel = path.relative_to(root)
        if is_excluded(rel):
            continue
        if path.is_file():
            files[rel] = path
    return files


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def sync(source: Path, dest: Path, check: bool) -> tuple[bool, list[str]]:
    source_files = iter_files(source)
    dest_files = iter_files(dest) if dest.exists() else {}
    changed = False
    actions: list[str] = []

    for rel, source_path in sorted(source_files.items()):
        dest_path = dest / rel
        needs_copy = rel not in dest_files
        if not needs_copy and file_hash(source_path) != file_hash(dest_path):
            needs_copy = True
        if needs_copy:
            changed = True
            actions.append(f"copy {rel}")
            if not check:
                ensure_parent(dest_path)
                shutil.copy2(source_path, dest_path)

    for rel, dest_path in sorted(dest_files.items()):
        if rel not in source_files:
            changed = True
            actions.append(f"remove {rel}")
            if not check:
                dest_path.unlink()

    if not check and dest.exists():
        for directory in sorted(
            [p for p in dest.rglob("*") if p.is_dir()],
            key=lambda p: len(p.parts),
            reverse=True,
        ):
            try:
                directory.rmdir()
            except OSError:
                pass

    return changed, actions


def main() -> int:
    args = parse_args()
    source = resolve_source(args.source)
    dest = Path(args.dest).resolve()

    changed, actions = sync(source, dest, args.check)

    if args.check and changed:
        if not args.quiet:
            print(f"skills/core differs from {source} ({len(actions)} action(s) needed).")
        return 1

    if not args.quiet:
        if changed:
            print(f"Synced skills/core from {source} ({len(actions)} action(s)).")
        else:
            print("skills/core already matches OMP Mycelium plugin copy.")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:  # noqa: BLE001 - hook should surface actionable failures.
        print(f"sync-mycelium-skills-core: {error}", file=sys.stderr)
        raise SystemExit(2)
