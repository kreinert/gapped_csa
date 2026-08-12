#!/usr/bin/env python3
"""Compare tree-dp and tree-dp4 (|C| and wall time) over inputs x shapes x max_add.

Reports the minimum of several runs: wall time here is noisy enough that a
single sample was showing +47% for a pass that is provably O(#cycles).
"""
import argparse
import re
import subprocess
import sys
import time

ALGOS = ["tree-dp", "tree-dp4"]


def run_once(fasta, shape, algo, max_add):
    cmd = ["./gcsa", "-g", fasta, "-s", shape, "--algo", algo,
           "--max-add", str(max_add)]
    t0 = time.time()
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    dt = time.time() - t0
    m = re.search(r"\(\|C\|\): (\d+)", out)
    return (int(m.group(1)) if m else None), dt


def best_of(fasta, shape, algo, max_add, reps, budget):
    c, best = run_once(fasta, shape, algo, max_add)
    spent = best
    while reps > 1 and spent < budget:
        _, dt = run_once(fasta, shape, algo, max_add)
        best = min(best, dt)
        spent += dt
        reps -= 1
    return c, best


ap = argparse.ArgumentParser()
ap.add_argument("--files", nargs="+", default=["r300-3000.fasta", "repeats-long.fasta"])
ap.add_argument("--shapes", nargs="+", default=["#.#", "#####", "####.####"])
ap.add_argument("--max-adds", nargs="+", type=int, default=[8, 16])
ap.add_argument("--reps", type=int, default=3)
ap.add_argument("--budget", type=float, default=45.0, help="seconds per cell")
args = ap.parse_args()

cols = "".join(f"{a:>12}{'sec':>7}" for a in ALGOS)
hdr = f"{'file':<20}{'shape':<11}{'ma':<4}{cols}{'d|C|':>8}{'ovh%':>7}"
print(hdr)
print("-" * len(hdr))
for f in args.files:
    for s in args.shapes:
        for ma in args.max_adds:
            res = [best_of(f, s, a, ma, args.reps, args.budget) for a in ALGOS]
            row = f"{f:<20}{s:<11}{ma:<4}" + "".join(f"{c:>12}{t:>7.2f}" for c, t in res)
            row += f"{res[1][0] - res[0][0]:>8}"
            row += f"{100 * (res[1][1] - res[0][1]) / res[0][1]:>7.1f}"
            print(row)
            sys.stdout.flush()
