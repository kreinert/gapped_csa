#!/usr/bin/env python3
"""Compare tree-dp and tree-dp3 (|C| and wall time) over inputs x shapes x max_add."""
import re
import subprocess
import sys
import time

FILES = ["r50-600.fasta", "r300-300.fasta", "r600-60.fasta", "repeats-long.fasta"]
SHAPES = ["#.#", "####.####", "#####"]
MAX_ADDS = [8, 16]


def run(fasta, shape, algo, max_add):
    cmd = ["./gcsa", "-g", fasta, "-s", shape, "--algo", algo, "--max-add", str(max_add)]
    t0 = time.time()
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    dt = time.time() - t0
    m = re.search(r"\(\|C\|\): (\d+)", out)
    return (int(m.group(1)) if m else None), dt


hdr = f"{'file':<20}{'shape':<11}{'ma':<4}{'tree-dp':>10}{'sec':>7}{'tree-dp3':>10}{'sec':>7}{'gain%':>8}"
print(hdr)
print("-" * len(hdr))
for f in FILES:
    for s in SHAPES:
        for ma in MAX_ADDS:
            c1, t1 = run(f, s, "tree-dp", ma)
            c3, t3 = run(f, s, "tree-dp3", ma)
            if c1 is None or c3 is None:
                print(f"{f:<20}{s:<11}{ma:<4}  ERROR")
                continue
            gain = 100.0 * (c1 - c3) / c1
            print(f"{f:<20}{s:<11}{ma:<4}{c1:>10}{t1:>7.2f}{c3:>10}{t3:>7.2f}{gain:>8.2f}")
            sys.stdout.flush()
