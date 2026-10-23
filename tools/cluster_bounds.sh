#!/bin/bash
# One-shot sweep_bounds.py call over run_suite.py's SA cache (reads the
# sidecars in SA_CACHE/meta, needs no FASTA):
#   tools/cluster_bounds.sh SA_CACHE OUT_CSV [DATASET ...]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
CACHE=$1; OUT=$2; shift 2
python3 "$HERE/sweep_bounds.py" --sa-cache "$CACHE" --out "$OUT" --adds 8,16,32,64,128 --iters 30 \
  ${1:+--only "$@"}
