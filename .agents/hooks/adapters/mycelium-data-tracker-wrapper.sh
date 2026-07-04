#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel 2>/dev/null || pwd)
export MYCELIUM_DATA_HELPER="$REPO_ROOT/.agents/hooks/adapters/extract-data-lineage-event-r.py"
exec "$REPO_ROOT/skills/core/hooks/mycelium-data-tracker.sh"
