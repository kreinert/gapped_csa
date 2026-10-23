#!/usr/bin/env python3
"""Lower bounds on |C| for every (dataset, shape) cell of a sweep.

Cluster use -- straight from run_suite.py's prepared SA cache, no FASTA needed:

  sweep_bounds.py --sa-cache $S/sa --out bounds.csv \
      [--only DATASET ...] [--adds 8,16,32,64,128] [--iters 30] \
      [--runs runs.csv] [--shard I/N]

  Cells are read from the cache's sidecars ($S/sa/meta/<dataset>__<shape>.json,
  written by run_suite.py --prepare-run); each one names its blob, which is
  passed to ./lower_bound --sa-load.  --runs (any CSV with dataset, shape,
  c_final columns, e.g. data/20260923_runs.csv) is optional: the best c_final
  per cell only steers the Lagrangian step size and limits the cells to ones
  that were actually run.

Local use -- from FASTAs (the SA is built per cell):

  sweep_bounds.py --runs runs.csv --data-dir DIR --out bounds.csv

One output row per (cell, max_add).  The add scan is cumulative, so asking for
8,16,32,64,128 costs one scan to 128 plus one Lagrangian pass per value.
Re-running with the same --out skips cells already present (resumable).
"""
import argparse, csv, glob, json, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--sa-cache", help="run_suite.py's --sa-cache directory (holds *.gsa and meta/)")
ap.add_argument("--data-dir", help="FASTA directory (only without --sa-cache)")
ap.add_argument("--runs", help="CSV with dataset,shape,c_final (optional with --sa-cache)")
ap.add_argument("--adds", default="8,16,32,64,128")
ap.add_argument("--iters", default="30")
ap.add_argument("--only", nargs="*", help="dataset names to include")
ap.add_argument("--shard", default="0/1", help="I/N: this task's slice of the cells")
ap.add_argument("--lower-bound", default=os.path.join(ROOT, "lower_bound"))
a = ap.parse_args()
if not a.sa_cache and not a.data_dir:
    sys.exit("need --sa-cache or --data-dir")

# best known |C| per cell (optional)
best = {}
if a.runs:
    for r in csv.DictReader(open(a.runs)):
        if r.get("complete", "True") not in ("True", "true", "1"):
            continue
        try:
            c = int(float(r["c_final"]))
        except (KeyError, ValueError):
            continue
        k = (r["dataset"], r["shape"])
        best[k] = min(best.get(k, c), c)

# cell list: (dataset, shape, blob-or-None)
cells = []
if a.sa_cache:
    for p in sorted(glob.glob(os.path.join(a.sa_cache, "meta", "*.json"))):
        try:
            m = json.load(open(p))
        except (OSError, ValueError):
            continue
        if not all(k in m for k in ("blob", "dataset", "shape")):
            continue
        blob = os.path.join(a.sa_cache, m["blob"])
        if os.path.isfile(blob):
            cells.append((m["dataset"], m["shape"], blob))
    if a.runs:
        cells = [c for c in cells if (c[0], c[1]) in best]
else:
    if not a.runs:
        sys.exit("without --sa-cache, --runs is needed to list the cells")
    cells = [(ds, sh, None) for (ds, sh) in sorted(best)]
if a.only:
    cells = [c for c in cells if c[0] in a.only]

si, sn = (int(x) for x in a.shard.split("/"))
cells = cells[si::sn]

done = set()
if os.path.exists(a.out):
    for r in csv.DictReader(open(a.out)):
        done.add((r["dataset"], r["shape"]))


def fasta_for(ds):
    p = os.path.join(a.data_dir, ds + ".fasta")
    if os.path.exists(p):
        return p
    # "provided" datasets keep their own file name (bench/datasets.py)
    sys.path.insert(0, os.path.join(ROOT, "bench"))
    try:
        import datasets
        ent = next((x for x in datasets.DATASETS if x["name"] == ds), {})
        if ent.get("path"):
            return os.path.join(a.data_dir, ent["path"])
    except Exception:
        pass
    return p


print(f"[sweep_bounds] shard {si}/{sn}: {len(cells)} cell(s), {len(done)} already in {a.out}", flush=True)
f = open(a.out, "a", newline="")
w = None
fields = None
if os.path.getsize(a.out) > 0:
    fields = next(csv.reader(open(a.out)))
for ds, shape, blob in cells:
    if (ds, shape) in done:
        continue
    cmd = [a.lower_bound, "-s", shape, "--max-add", a.adds, "--iters", a.iters]
    cmd += ["--sa-load", blob] if blob else ["-g", fasta_for(ds)]
    if (ds, shape) in best:
        cmd += ["--ub", str(best[(ds, shape)])]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    lines = [l for l in p.stdout.splitlines() if l.startswith("m=")]
    if p.returncode != 0 or not lines:
        print(f"[FAIL] {ds} {shape} rc={p.returncode}: {p.stderr.strip()[-300:]}", flush=True)
        continue
    for line in lines:
        row = dict(dataset=ds, shape=shape, **dict(re.findall(r"(\w+)=(\S+)", line)))
        if w is None:
            w = csv.DictWriter(f, fieldnames=fields or list(row.keys()))
            if fields is None:
                w.writeheader()
        w.writerow(row)
        f.flush()
    print(f"{ds} {shape} {time.time() - t0:.0f}s", flush=True)
