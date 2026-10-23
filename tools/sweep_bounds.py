#!/usr/bin/env python3
"""Lower bounds for every (dataset, shape) cell of a sweep CSV.

  sweep_bounds.py --runs 20260923_runs.csv --data-dir DIR --out bounds.csv \
      [--sa-cache DIR] [--adds 8,16,32,64,128] [--only DATASET ...]

Looks for DIR/<dataset>.fasta (or uses --sa-cache blobs, which is what the
cluster sweep has).  One SA build per cell; one output row per (cell, max_add).
"""
import argparse, csv, os, re, subprocess, sys, time
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
ap = argparse.ArgumentParser()
ap.add_argument("--runs", required=True); ap.add_argument("--data-dir", required=True)
ap.add_argument("--out", required=True); ap.add_argument("--sa-cache")
ap.add_argument("--adds", default="8,16,32,64,128"); ap.add_argument("--iters", default="300")
ap.add_argument("--only", nargs="*")
a = ap.parse_args()
d = pd.read_csv(a.runs); d = d[d.complete == True]
cells = d.groupby(["dataset", "shape"]).c_final.min().reset_index()
done = set()
if os.path.exists(a.out):
    for r in csv.DictReader(open(a.out)): done.add((r["dataset"], r["shape"]))
f = open(a.out, "a", newline=""); w = None
for _, c in cells.iterrows():
    if a.only and c.dataset not in a.only: continue
    if (c.dataset, c["shape"]) in done: continue
    fasta = os.path.join(a.data_dir, c.dataset + ".fasta")
    if not os.path.exists(fasta):   # "provided" datasets keep their own file name (bench/datasets.py)
        sys.path.insert(0, os.path.join(ROOT, "bench"))
        try:
            import datasets
            ent = next((x for x in datasets.DATASETS if x["name"] == c.dataset), {})
            if ent.get("path"):
                fasta = os.path.join(a.data_dir, ent["path"])
        except Exception:
            pass
    cmd = [os.path.join(ROOT, "lower_bound"), "-g", fasta, "-s", c["shape"], "--max-add", a.adds,
           "--ub", str(int(c.c_final)), "--iters", a.iters]
    if a.sa_cache: cmd += ["--sa-cache", a.sa_cache]
    t0 = time.time(); out = subprocess.run(cmd, capture_output=True, text=True).stdout
    for line in out.strip().splitlines():
        row = dict(dataset=c.dataset, shape=c["shape"], **dict(re.findall(r"(\w+)=(\S+)", line)))
        if w is None:
            w = csv.DictWriter(f, fieldnames=list(row.keys()))
            if os.path.getsize(a.out) == 0: w.writeheader()
        w.writerow(row); f.flush()
    print(c.dataset, c["shape"], f"{time.time()-t0:.0f}s", flush=True)
