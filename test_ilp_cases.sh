#!/bin/bash
# Small-instance ILP sanity suite: opt must be <= every heuristic's |C|.
set -u
cases=(
  "ACGTCTTAAACCCTCGTCTTAAACCCAACGTCTTAAACCC|#.#|note"
  "GCCTTTAAAGGCCTTTAAAGGCCTTTAAAG|#.#|GCCTTTAAAGx3"
  "ACACACACACAC|#.#|ACx6"
  "ACGTACGTACGTACGTACGT|#.#|ACGTx5"
  "GGGCGGCGGC|##|GGGCGGCGGC"
  "ACGTCTTAAACCCTCGTCTTAAACCCAACGTCTTAAACCC|##.##|note-##.##"
  "ACGTCTTAAACCCTCGTCTTAAACCCAACGTCTTAAACCC|##|note-##"
  "AAAAAAAAAAAAAAAA|#.#|A16"
  "ACGTACGTACGTACGTACGTACGTACGT|##.#|ACGTx7"
  "GCCTTTAAAGGCCTTTAAAGGCCTTTAAAGGCCTTTAAAG|#..#|GCCTTTAAAGx4"
)
fail=0
for c in "${cases[@]}"; do
  IFS='|' read -r text shape tag <<< "$c"
  echo "### $tag  (n=${#text} shape=$shape)"
  ./ilp_baseline "$text" "$shape" 2>&1 \
    | grep -E "universe=|^optimal|^heuristic|^self-check|^decode-check|^bounds"
  rc=${PIPESTATUS[0]}
  [ "$rc" != "0" ] && { echo "  -> ilp_baseline exit=$rc"; fail=1; }
  echo
done
exit $fail
