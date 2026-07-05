#!/usr/bin/env bash
# Project-owned wrapper for the synced Mycelium SessionStart health hook.
# It preserves the heuristic INDEX.md knowledge summary when the synced health
# hook regenerates it in a way that regresses recent semantic ordering.
set -uo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel 2>/dev/null || pwd)
CORE_HOOK="$REPO_ROOT/skills/core/hooks/mycelium-health.sh"
GUARD="$REPO_ROOT/tools/mycelium-provenance-guard.py"

TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/mycelium-health-wrapper.XXXXXX") || exit 1
cleanup() {
  rm -rf "$TMP_DIR"
}
trap cleanup EXIT INT TERM

STDIN_FILE="$TMP_DIR/stdin.json"
CORE_STDOUT="$TMP_DIR/core.stdout"
REGISTRY_BEFORE="$TMP_DIR/LOG_REGISTRY.before.md"
INDEX_BEFORE="$TMP_DIR/INDEX.before.md"
LAST_SESSION_BEFORE="$TMP_DIR/last-session.before.md"

cat > "$STDIN_FILE"
: > "$CORE_STDOUT"

if [ -f "$REPO_ROOT/.living/log/LOG_REGISTRY.md" ]; then
  cp "$REPO_ROOT/.living/log/LOG_REGISTRY.md" "$REGISTRY_BEFORE"
fi
if [ -f "$REPO_ROOT/.living/INDEX.md" ]; then
  cp "$REPO_ROOT/.living/INDEX.md" "$INDEX_BEFORE"
fi
if [ -f "$REPO_ROOT/.claude/last-session.md" ]; then
  cp "$REPO_ROOT/.claude/last-session.md" "$LAST_SESSION_BEFORE"
fi

CORE_STATUS=0
if [ -x "$CORE_HOOK" ]; then
  "$CORE_HOOK" < "$STDIN_FILE" > "$CORE_STDOUT"
  CORE_STATUS=$?
elif [ -f "$CORE_HOOK" ]; then
  bash "$CORE_HOOK" < "$STDIN_FILE" > "$CORE_STDOUT"
  CORE_STATUS=$?
fi

if [ -f "$GUARD" ]; then
  python3 "$GUARD" \
    --repo-root "$REPO_ROOT" \
    --registry-before "$REGISTRY_BEFORE" \
    --index-before "$INDEX_BEFORE" \
    --last-session-before "$LAST_SESSION_BEFORE" \
    >/dev/null 2>"$TMP_DIR/guard.stderr" || true
fi

cat "$CORE_STDOUT"
exit "$CORE_STATUS"
