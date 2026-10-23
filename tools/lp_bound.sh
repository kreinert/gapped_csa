#!/bin/bash
# LP-relaxation lower bound of ilp_baseline's full-universe ILP:  lp_bound.sh FASTA SHAPE MAXADD [TIME]
# (writes the .lp without solving -- no cbc/glpsol on PATH -- then solves the relaxation with HiGHS)
D=$(cd "$(dirname "$0")/.." && pwd); T=$(grep -v ">" "$1" | tr -d '\n'); LP=/tmp/lpb_$$.lp
env PATH=/usr/bin:/bin "$D/ilp_baseline" "$T" "$2" --max-add "$3" --lp-out $LP > /tmp/lpb_$$.out 2>&1
C=$(grep -o "candidates=[0-9]*" /tmp/lpb_$$.out)
R=$(python3 "$D/tools/solve_lp_highs.py" $LP "${4:-600}" --relax)
echo "fasta=$(basename $1) shape=$2 max_add=$3 $C $R"; rm -f $LP /tmp/lpb_$$.out
