#!/usr/bin/env bash
set -euo pipefail

INPUT=$(cat)
COMMAND=$(printf '%s' "$INPUT" | python3 -c 'import json, sys
try:
    data = json.load(sys.stdin)
except Exception:
    print("")
    raise SystemExit
for key in ("tool_input", "input"):
    value = data.get(key)
    if isinstance(value, dict) and isinstance(value.get("command"), str):
        print(value["command"])
        raise SystemExit
print("")
' 2>/dev/null || true)

SHOULD_SKIP=$(COMMAND="$COMMAND" python3 -c 'import os, re
cmd = os.environ.get("COMMAND", "").replace("\\\\", "/").strip()
if re.search(r"[;&|]", cmd):
    print("false")
    raise SystemExit
core = "|".join([
    "generate_index",
    "validate_structure",
    "recall_lessons",
    "detect_recurrence",
    "crystallize_findings",
    "init_knowledge",
    "migrate_existing_repos",
    "install_convention",
    "init_repo",
])
pattern = re.compile(rf"^(?:python3?|uv\s+run\s+python3?)\s+(?:(?:[^\s]*/)?skills/core/scripts/(?:{core})\.py|(?:[^\s]*/)?tools/sync-mycelium-skills-core\.py)(?:\s+.*)?$")
print("true" if pattern.search(cmd) else "false")
')

if [ "$SHOULD_SKIP" = "true" ]; then
  exit 0
fi

REPO_ROOT=$(git rev-parse --show-toplevel 2>/dev/null || pwd)
HOOK="$REPO_ROOT/skills/core/hooks/mycelium-post-action.sh"
if [ ! -x "$HOOK" ]; then
  exit 0
fi

printf '%s' "$INPUT" | "$HOOK"
