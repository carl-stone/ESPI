#!/usr/bin/env bash
# Shared generated-file edit policy.
#
# Harness adapters call this script with normalized write-risk paths. The script
# knows repository policy only; it does not parse Codex, OMP, Claude, or any
# other harness payload.
#
# Usage:
#   generated-files.sh --operation write --path man/foo.Rd [--path README.md]
#
# Exit codes:
#   0: allow
#   2: block

set -euo pipefail

operation=""
paths=()

while [ "$#" -gt 0 ]; do
  case "$1" in
    --operation)
      if [ "$#" -lt 2 ]; then
        echo "generated-files policy: --operation requires a value" >&2
        exit 2
      fi
      operation="$2"
      shift 2
      ;;
    --path)
      if [ "$#" -lt 2 ]; then
        echo "generated-files policy: --path requires a value" >&2
        exit 2
      fi
      paths+=("$2")
      shift 2
      ;;
    --)
      shift
      break
      ;;
    *)
      echo "generated-files policy: unknown argument '$1'" >&2
      exit 2
      ;;
  esac
done

case "$operation" in
  write|edit|delete|move)
    ;;
  "")
    echo "generated-files policy: --operation is required" >&2
    exit 2
    ;;
  *)
    exit 0
    ;;
esac

if [ "${#paths[@]}" -eq 0 ]; then
  exit 0
fi

repo_root=$(git rev-parse --show-toplevel 2>/dev/null || pwd -P)
current_working_dir=$(pwd -P)

normalize_path() {
  local path="$1"
  local absolute_path
  local parent_dir
  local base_name

  # Remove URI-style local prefix if an adapter passed one through.
  path="${path#local://}"

  # Interpret relative paths from the harness session cwd, then reduce paths
  # with ../ segments when the parent directory exists. Generated-file targets
  # in this repo all have existing parents (repo root or man/).
  case "$path" in
    /*)
      absolute_path="$path"
      ;;
    *)
      absolute_path="$current_working_dir/$path"
      ;;
  esac

  parent_dir=$(dirname -- "$absolute_path")
  base_name=$(basename -- "$absolute_path")
  if cd "$parent_dir" 2>/dev/null; then
    absolute_path="$(pwd -P)/$base_name"
    cd "$current_working_dir" >/dev/null
  fi

  case "$absolute_path" in
    "$repo_root"/*)
      path="${absolute_path#"$repo_root"/}"
      ;;
    *)
      path="$absolute_path"
      ;;
  esac

  while [[ "$path" == ./* ]]; do
    path="${path#./}"
  done

  printf '%s\n' "$path"
}

remediation_message() {
  local path="$1"

  if [[ "$path" =~ ^man/.*\.Rd$ ]]; then
    echo "Edit roxygen comments in R/ instead, then run: devtools::document()"
  elif [ "$path" = "NAMESPACE" ]; then
    echo "Edit roxygen imports/exports in R/ instead, then run: devtools::document()"
  elif [ "$path" = "README.md" ]; then
    echo "Edit README.Rmd instead, then run: devtools::build_readme()"
  else
    echo "This file is generated and should not be manually edited."
  fi
}

is_generated_file() {
  local path="$1"

  if [[ "$path" =~ ^man/.*\.Rd$ ]]; then
    return 0
  fi

  case "$path" in
    NAMESPACE|README.md)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

blocked=()

for raw_path in "${paths[@]}"; do
  path=$(normalize_path "$raw_path")
  if is_generated_file "$path"; then
    blocked+=("BLOCKED: '$path' is generated. $(remediation_message "$path")")
  fi
done

if [ "${#blocked[@]}" -gt 0 ]; then
  printf '%s\n' "${blocked[@]}" >&2
  exit 2
fi

exit 0
