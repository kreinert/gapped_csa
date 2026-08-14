#!/usr/bin/env python3
"""Trends for candidate benchmark families around cycle-edge repair (tree-dp4).

Cycle repair pays when the preference graph has a cycle whose last-processed
member (largest packed k-mer name) is not the one holding the cheapest edge.
It can also *hurt*: it improves the forest DP's local objective but can hand
Phase II a worse-rooted chain.
"""
import argparse
import re
import subprocess
import sys

ALGOS = ["greedy", "dep-order", "tree-dp", "tree-dp3", "tree-dp4"]


def heur(text, shape, want_opt):
    out = subprocess.run(["./ilp_baseline", text, shape],
                         capture_output=True, text=True).stdout
    vals = dict(re.findall(r"(greedy|dep-order|tree-dp\d?)=(\d+)",
                           re.search(r"^heuristic \|C\|:.*$", out, re.M).group(0)))
    m = re.search(r"^optimal \|C\| = (\d+)", out, re.M)
    opt = int(m.group(1)) if (m and want_opt) else None
    bounds = "bounds all algos: YES" in out
    return {k: int(v) for k, v in vals.items()}, opt, bounds


ap = argparse.ArgumentParser()
ap.add_argument("--motifs", nargs="+", default=["AC", "TAAA", "GACA", "GTAC", "ACGT"])
ap.add_argument("--shape", default="#.#")
ap.add_argument("--lengths", nargs="+", type=int, default=[12, 20, 28, 36, 48])
ap.add_argument("--opt", action="store_true", help="also report the ILP optimum")
args = ap.parse_args()

hdr = f"{'motif':<8}{'n':<5}{'opt':>5}" + "".join(f"{a:>10}" for a in ALGOS) + f"{'d(t4-t)':>9}{'bnd':>5}"
print(hdr)
print("-" * len(hdr))
for motif in args.motifs:
    for n in args.lengths:
        text = (motif * (n // len(motif) + 1))[:n]
        vals, opt, bounds = heur(text, args.shape, args.opt)
        row = f"{motif:<8}{n:<5}{(opt if opt is not None else '-'):>5}"
        row += "".join(f"{vals.get(a, 0):>10}" for a in ALGOS)
        row += f"{vals['tree-dp4'] - vals['tree-dp']:>9}"
        row += f"{'YES' if bounds else 'no':>5}"
        print(row)
        sys.stdout.flush()
