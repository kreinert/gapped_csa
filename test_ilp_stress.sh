#!/bin/bash
# Randomized stress: for many small instances the ILP optimum must decode
# correctly and be <= every heuristic's |C|.
set -u
N=${1:-120}
SEED=${2:-1}
shapes=("#.#" "##" "##.#" "#.#.#" "#..#" "###")
fail=0
for ((i=0; i<N; i++)); do
  len=$(( (RANDOM % 30) + 10 ))
  sh=${shapes[$(( RANDOM % ${#shapes[@]} ))]}
  # alphabet size 2 or 4 -> more repetition on the small alphabet
  if (( RANDOM % 2 )); then alpha="AC"; else alpha="ACGT"; fi
  text=""
  for ((j=0; j<len; j++)); do
    k=$(( RANDOM % ${#alpha} ))
    text+="${alpha:$k:1}"
  done
  out=$(./ilp_baseline "$text" "$sh" 2>&1)
  rc=$?
  if [ "$rc" != "0" ]; then
    echo "FAIL rc=$rc text=$text shape=$sh"
    echo "$out" | grep -E "^optimal|^heuristic|^self-check|^decode-check|^bounds"
    fail=$((fail+1))
  fi
done
echo "stress: $N instances, failures=$fail"
exit $(( fail > 0 ))
