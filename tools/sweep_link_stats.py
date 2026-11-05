#!/usr/bin/env python3
"""Re-run sweep cells through ./link_stats to record the largest add each run
actually ends with (Phase II's retarget pass composes adds past --max-add, so
a lower bound is only valid for a run at an A >= this value).

  sweep_link_stats.py --sa-cache $S/sa --out linkstats.csv \
      [--only DATASET ...] [--algos dep-order pseudoforest-dp ...] \
      [--max-adds 8 16 32] [--shard I/N]

Uses the SAME environment as the sweep job -- set the GCSA_* variables in the
job script exactly as for run_suite.py (they are inherited).  Set
GCSA_CROSS_LCP=1 in the environment to test the wider candidate universe; the
column `cross_lcp` records which it was.  Resumable on --out.
"""
import argparse, csv, glob, json, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
SHORT = {"greedy-size": "gs", "greedy-degree": "gd", "dep-order": "dep",
         "pseudoforest-dp": "pfdp", "pseudoforest-dp-iterate": "pfdpi", "run-dp": "rundp"}
ap = argparse.ArgumentParser()
ap.add_argument("--sa-cache", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--only", nargs="*")
ap.add_argument("--algos", nargs="*", default=list(SHORT))
ap.add_argument("--max-adds", nargs="*", type=int, default=[8, 16, 32])
ap.add_argument("--shard", default="0/1")
ap.add_argument("--link-stats", default=os.path.join(ROOT, "link_stats"))
ap.add_argument("--plain-k", nargs="*",
                help="k values to run on each dataset's shared plain-SA blob "
                     "(default: the integer shapes in bench/datasets.py SHAPES)")
a = ap.parse_args()

from sa_cache_cells import sidecar_cells, default_plain_ks
cells = sidecar_cells(a.sa_cache, a.plain_k if a.plain_k is not None else default_plain_ks(), a.only)
jobs = [(ds, sh, b, al, A) for (ds, sh, b) in cells for al in a.algos for A in a.max_adds]
si, sn = (int(x) for x in a.shard.split("/"))
jobs = jobs[si::sn]
xl = os.environ.get("GCSA_CROSS_LCP", "0") not in ("", "0", "false", "no", "off")

done = set()
if os.path.exists(a.out):
    for r in csv.DictReader(open(a.out)):
        done.add((r["dataset"], r["shape"], r["algo"], r["max_add"], r["cross_lcp"]))
new = not os.path.exists(a.out) or os.path.getsize(a.out) == 0
f = open(a.out, "a", newline="")
cols = ["dataset", "shape", "algo", "max_add", "cross_lcp", "C", "m", "links",
        "max_add_used", "links_over_cap", "positions_covered_by_over_cap", "secs"]
w = csv.DictWriter(f, fieldnames=cols)
if new:
    w.writeheader()
print(f"[sweep_link_stats] shard {si}/{sn}: {len(jobs)} run(s), cross_lcp={xl}", flush=True)
for ds, sh, blob, al, A in jobs:
    if (ds, sh, al, str(A), str(xl)) in done:
        continue
    t0 = time.time()
    p = subprocess.run([a.link_stats, blob, sh, str(A), SHORT.get(al, al)],
                       capture_output=True, text=True)
    k = dict(re.findall(r"(\w+)=(\S+)", p.stdout))
    if p.returncode != 0 or "C" not in k:
        print(f"[FAIL] {ds} {sh} {al} {A}: rc={p.returncode} {p.stderr.strip()[-300:]}", flush=True)
        continue
    w.writerow(dict(dataset=ds, shape=sh, algo=al, max_add=A, cross_lcp=xl,
                    secs=round(time.time() - t0, 1), **{c: k.get(c, "") for c in cols[5:11]}))
    f.flush()
    print(f"{ds} {sh} {al} A={A}: C={k['C']} max_add_used={k['max_add_used']}", flush=True)
