#!/usr/bin/env bash
# Codex PreToolUse adapter for generated-file protection.
#
# Reads Codex hook JSON on stdin, extracts write-risk paths from direct file
# fields or apply_patch command text, then delegates to the shared policy.

set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
policy="$script_dir/../policies/generated-files.sh"

input=$(cat)

if ! command -v jq >/dev/null 2>&1; then
  echo "WARNING: jq is required to inspect Codex hook input JSON. Blocking as a precaution." >&2
  exit 2
fi

tool_name=$(jq -r '.tool_name // empty' <<<"$input") || {
  echo "WARNING: Failed to parse Codex hook input JSON. Blocking as a precaution." >&2
  exit 2
}

paths=()

add_path() {
  local path="$1"
  if [ -n "$path" ]; then
    paths+=("$path")
  fi
}

while IFS= read -r path; do
  add_path "$path"
done < <(
  jq -r '
    .tool_input.file_path
    // .tool_input.path
    // .tool_input.filename
    // .tool_input.filePath
    // .file_path
    // .path
    // .filename
    // empty
  ' <<<"$input"
)

command=$(jq -r '
  .tool_input.command
  // .tool_input.patch
  // .tool_input.input
  // .tool_input.content
  // .tool_input.text
  // empty
' <<<"$input")

# Codex file edits are commonly apply_patch calls whose paths live inside the
# patch command text. Support the common patch headers, OMP hashline headers
# when routed through Codex-compatible tools, plus a few diff-like forms for
# resilience.
if [ -n "$command" ]; then
  while IFS= read -r path; do
    add_path "$path"
  done < <(
    printf '%s\n' "$command" | awk '
      /^\*\*\* (Add|Update|Delete) File: / {
        sub(/^\*\*\* (Add|Update|Delete) File: /, "")
        print
        next
      }
      /^\*\*\* Rename from: / {
        sub(/^\*\*\* Rename from: /, "")
        print
        next
      }
      /^\*\*\* Rename to: / {
        sub(/^\*\*\* Rename to: /, "")
        print
        next
      }
      /^\[[^#\]\r\n]+#[0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f]\]$/ {
        sub(/^\[/, "")
        sub(/#[0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f]\]$/, "")
        print
        next
      }
      /^--- a\// {
        sub(/^--- a\//, "")
        if ($0 != "/dev/null") print
        next
      }
      /^\+\+\+ b\// {
        sub(/^\+\+\+ b\//, "")
        if ($0 != "/dev/null") print
        next
      }
    '
  )
fi

if [ "${#paths[@]}" -eq 0 ]; then
  exit 0
fi

# De-duplicate paths while preserving order.
unique_paths=()
for path in "${paths[@]}"; do
  seen=false
  for existing in "${unique_paths[@]}"; do
    if [ "$existing" = "$path" ]; then
      seen=true
      break
    fi
  done
  if [ "$seen" = false ]; then
    unique_paths+=("$path")
  fi
done

args=(--operation write)
for path in "${unique_paths[@]}"; do
  args+=(--path "$path")
done

set +e
output=$("$policy" "${args[@]}" 2>&1)
status=$?
set -e

case "$status" in
  0)
    exit 0
    ;;
  2)
    echo "$output" >&2
    exit 2
    ;;
  *)
    echo "WARNING: generated-file policy failed for Codex tool '$tool_name'. Blocking as a precaution." >&2
    if [ -n "$output" ]; then
      echo "$output" >&2
    fi
    exit 2
    ;;
esac
