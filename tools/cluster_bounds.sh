#!/bin/bash
# Lower bounds for the three datasets too large to bound off-cluster
# (cod_freq_kmers, human_freq_kmers, pangenome_ecoli_real_n100), reusing the
# 2026-09-23 sweep's prepared SA blobs, so no SA is rebuilt.
#
#   sbatch/srun ... tools/cluster_bounds.sh RUNS_CSV DATA_DIR SA_CACHE OUT_CSV
#
# DATA_DIR must hold <dataset>.fasta (only read to derive the content hash of
# the blob name -- the same rule as gcsa --sa-cache). Memory: ~60 B per SA
# entry on top of the blob, i.e. ~35 GB for human_freq_kmers / n100.
# These three ran with Phase II OFF, so no run can hold an add > --max-add and
# the bound at each run's own max-add is valid for it.
set -euo pipefail
RUNS=$1; DATA=$2; CACHE=$3; OUT=$4
HERE=$(cd "$(dirname "$0")" && pwd)
python3 "$HERE/sweep_bounds.py" --runs "$RUNS" --data-dir "$DATA" --sa-cache "$CACHE" \
    --out "$OUT" --adds 8,16,32 --iters 30 \
    --only cod_freq_kmers human_freq_kmers pangenome_ecoli_real_n100
