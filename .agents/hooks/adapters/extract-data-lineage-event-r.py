#!/usr/bin/env python3
"""ESPI wrapper for Mycelium data-lineage event extraction.

The synced Mycelium extractor detects Rscript invocations but only has Python
I/O regexes. This wrapper keeps the synced extractor as the implementation
source and adds lightweight R I/O expression detection for ESPI's R scripts.
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
CORE_EXTRACTOR = REPO_ROOT / "skills/core/scripts/extract_data_lineage_event.py"

spec = importlib.util.spec_from_file_location("mycelium_extract_data_lineage_event", CORE_EXTRACTOR)
if spec is None or spec.loader is None:
    print(f"Could not load {CORE_EXTRACTOR}", file=sys.stderr)
    raise SystemExit(2)

core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)

BASE_SCAN_SOURCE = core.scan_source

R_INPUT_REGEXES = [
    re.compile(r"""\breadRDS\s*\(\s*([^,\)\n]+)"""),
    re.compile(r"""\b(?:readr::)?read_(?:csv|tsv|rds|delim|table)\s*\(\s*(?:file\s*=\s*)?([^,\)\n]+)"""),
    re.compile(r"""\b(?:utils::)?read\.(?:csv|table|delim)\s*\(\s*(?:file\s*=\s*)?([^,\)\n]+)"""),
]

R_OUTPUT_REGEXES = [
    re.compile(r"""\bsaveRDS\s*\(\s*[^,\n]+,\s*(?:file\s*=\s*)?([^,\)\n]+)"""),
    re.compile(r"""\b(?:ggplot2::)?ggsave\s*\(\s*(?:filename\s*=\s*)?([^,\)\n]+)"""),
    re.compile(r"""\b(?:utils::)?write\.(?:table|csv)\s*\(\s*[^,\n]+,\s*(?:file\s*=\s*)?([^,\)\n]+)"""),
    re.compile(r"""\bwrite_tsv\s*\(\s*[^,\n]+,\s*([^,\)\n]+)"""),
]


def clean_r_expr(expr: str) -> str:
    expr = expr.strip()
    if (expr.startswith('"') and expr.endswith('"')) or (expr.startswith("'") and expr.endswith("'")):
        return expr[1:-1]
    return expr


def scan_source_with_r(source: str) -> tuple[list[str], list[str], list[str], list[int]]:
    inputs, outputs, filters, seeds = BASE_SCAN_SOURCE(source)

    for pattern in R_INPUT_REGEXES:
        inputs.extend(clean_r_expr(match.group(1)) for match in pattern.finditer(source))

    for pattern in R_OUTPUT_REGEXES:
        outputs.extend(clean_r_expr(match.group(1)) for match in pattern.finditer(source))

    return core._dedupe(inputs), core._dedupe(outputs), filters, seeds


core.scan_source = scan_source_with_r

if __name__ == "__main__":
    raise SystemExit(core.main())
