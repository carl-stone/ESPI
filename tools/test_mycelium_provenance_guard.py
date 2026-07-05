#!/usr/bin/env python3
"""Regression tests for tools/mycelium-provenance-guard.py."""

from __future__ import annotations

import importlib.util
import tempfile
import textwrap
import unittest
from pathlib import Path


GUARD_PATH = Path(__file__).with_name("mycelium-provenance-guard.py")
SPEC = importlib.util.spec_from_file_location("mycelium_provenance_guard", GUARD_PATH)
assert SPEC is not None and SPEC.loader is not None
GUARD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GUARD)


class MyceliumProvenanceGuardTests(unittest.TestCase):
    def test_sentence_summary_is_semantic(self) -> None:
        row = (
            "| 2026-07-05 | 2026-07-05-005 | espi | main | 22m | 0 | "
            "Checked MG-selected UMAP point-size settings and regenerated existing outputs. | "
            "Confirmed point-size settings and rendered notebook | complete | "
            "mg-selected, plotting, umap, notebook | [log](2026-07-05-005-espi.md) |"
        )
        self.assertTrue(GUARD.registry_row_is_semantic(row))

    def test_index_restores_regressed_recent_order(self) -> None:
        before = textwrap.dedent(
            """\
            # Index
            <!-- BEGIN KNOWLEDGE SUMMARY -->
            ## Most recent (10)
            - [2026-07-05] L-20: Preserve semantic hook records
            - [2026-07-05] L-19: Use Rscript files for smoke tests
            ## By tag
            <!-- END KNOWLEDGE SUMMARY -->
            """
        )
        after = textwrap.dedent(
            """\
            # Index
            <!-- BEGIN KNOWLEDGE SUMMARY -->
            ## Most recent (10)
            - [2026-07-05] L-19: Use Rscript files for smoke tests
            - [2026-07-05] L-20: Preserve semantic hook records
            ## By tag
            <!-- END KNOWLEDGE SUMMARY -->
            """
        )
        with tempfile.TemporaryDirectory() as tmp:
            index_path = Path(tmp) / "INDEX.md"
            index_path.write_text(after, encoding="utf-8")
            changed = GUARD.repair_index(before, index_path)
            self.assertTrue(changed)
            self.assertEqual(index_path.read_text(encoding="utf-8"), before)

    def test_index_noop_when_recent_order_unchanged(self) -> None:
        before = textwrap.dedent(
            """\
            # Index
            <!-- BEGIN KNOWLEDGE SUMMARY -->
            ## Most recent (10)
            - [2026-07-05] L-20: Preserve semantic hook records
            - [2026-07-05] L-19: Use Rscript files for smoke tests
            ## By tag
            - `old`: L-20
            <!-- END KNOWLEDGE SUMMARY -->
            """
        )
        after = textwrap.dedent(
            """\
            # Index
            <!-- BEGIN KNOWLEDGE SUMMARY -->
            ## Most recent (10)
            - [2026-07-05] L-20: Preserve semantic hook records
            - [2026-07-05] L-19: Use Rscript files for smoke tests
            ## By tag
            - `new`: L-20, L-19
            <!-- END KNOWLEDGE SUMMARY -->
            """
        )
        with tempfile.TemporaryDirectory() as tmp:
            index_path = Path(tmp) / "INDEX.md"
            index_path.write_text(after, encoding="utf-8")
            changed = GUARD.repair_index(before, index_path)
            self.assertFalse(changed)
            self.assertEqual(index_path.read_text(encoding="utf-8"), after)


if __name__ == "__main__":
    unittest.main()
