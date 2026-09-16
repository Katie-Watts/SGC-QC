#!/usr/bin/env bash
# Remove the trailing "<variant> not found in dataset" lines that plink appends
# to *.verbose.clumped files. Edits files in place across a directory.
#
# These lines list variants from the clump input that weren't in the reference
# panel; they carry no clump information and sit at the end of the file.
# Any line containing "not found in dataset" is removed.
#
# Usage:   ./clean_verbose.sh [dir]
# Example: ./clean_verbose.sh clumped_1mb/results
#          ./clean_verbose.sh clumped_1mb

set -uo pipefail

DIR="${1:-.}"
[[ -d "$DIR" ]] || { echo "ERROR: directory not found: $DIR" >&2; exit 1; }

shopt -s nullglob
files=("$DIR"/*.verbose.clumped)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "No *.verbose.clumped files in $DIR" >&2; exit 0
fi

total_removed=0
for f in "${files[@]}"; do
  before=$(grep -c "not found in dataset" "$f" || true)
  if [[ "$before" -gt 0 ]]; then
    tmp="${f}.tmp.$$"
    grep -v "not found in dataset" "$f" > "$tmp" && mv "$tmp" "$f"
    echo "[clean] $(basename "$f"): removed $before line(s)"
    total_removed=$((total_removed + before))
  fi
done

echo "Done. Removed $total_removed 'not found in dataset' line(s) across ${#files[@]} file(s)."
