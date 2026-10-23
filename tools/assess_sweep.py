#!/usr/bin/env python3
"""Grade every run of a sweep against the lower bounds.

  assess_sweep.py --runs runs.csv --bounds 'bounds_*.csv' \
      [--link-stats 'linkstats_*.csv'] [--out assessed.csv]

--runs     one or more CSVs with dataset, shape, algo, max_add, m, c_final
           (globs allowed; e.g. data/20260923_runs.csv)
--bounds   sweep_bounds.py output(s)
--link-stats  sweep_link_stats.py output(s) WITHOUT GCSA_CROSS_LCP: the
           largest add each run used.

For each run the bound is taken at the smallest computed A that is >= the
largest add the run can contain:
  measured     link-stats says so;
  no Phase II  cfg == nophase2 (or rt_disabled): adds cannot exceed max_add;
  assumed      otherwise the largest A in the bounds file (flagged).
and LB*(A) = max over computed A' >= A of LB(A') (OPT never increases in A).

deletion efficiency = (m - c_final) / (m - LB*)  -- a guaranteed lower bound on
the share of the optimum's deletions the run achieves; 100 % = proven optimal.
"""
import argparse, glob
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--runs", nargs="+", required=True)
ap.add_argument("--bounds", nargs="+", required=True)
ap.add_argument("--link-stats", nargs="*", default=[])
ap.add_argument("--out", default="assessed_runs.csv")
a = ap.parse_args()
rd = lambda pats: pd.concat([pd.read_csv(p) for pat in pats for p in sorted(glob.glob(pat))], ignore_index=True)

d = rd(a.runs)
if "complete" in d:
    d = d[d.complete.astype(str).isin(["True", "true", "1"])]
d = d.dropna(subset=["c_final"]).copy()
b = rd(a.bounds).drop_duplicates(["dataset", "shape", "max_add"])
b = b.sort_values(["dataset", "shape", "max_add"])
b["LBs"] = b.groupby(["dataset", "shape"]).LB.transform(lambda s: s[::-1].cummax()[::-1])

if a.link_stats:
    ls = rd(a.link_stats)
    ls = ls[~ls.cross_lcp.astype(str).isin(["True", "true", "1"])]
    d = d.merge(ls[["dataset", "shape", "algo", "max_add", "max_add_used"]],
                on=["dataset", "shape", "algo", "max_add"], how="left")
else:
    d["max_add_used"] = np.nan
grid = {k: sorted(v) for k, v in b.groupby(["dataset", "shape"]).max_add.unique().items()}


def pick(r):
    g = grid.get((r.dataset, r["shape"]))
    if not g:
        return pd.Series([np.nan, "no bound"])
    nophase2 = str(r.get("cfg", "")) == "nophase2" or str(r.get("rt_disabled", "")) == "True"
    if not np.isnan(r.max_add_used):
        need, src = max(r.max_add, r.max_add_used), "measured"
    elif nophase2:
        need, src = r.max_add, "no Phase II"
    else:
        return pd.Series([max(g), f"assumed <= {max(g)}"])
    ok = [x for x in g if x >= need]
    return pd.Series([ok[0], src] if ok else [max(g), f"INVALID: used {int(need)} > {max(g)}"])


d[["A_bound", "A_src"]] = d.apply(pick, axis=1)
d = d.merge(b[["dataset", "shape", "max_add", "LB0", "LB1", "DOM", "LAG", "LBs"]]
            .rename(columns={"max_add": "A_bound"}), on=["dataset", "shape", "A_bound"], how="left")
d["eff"] = np.where(d.m > d.LBs, (d.m - d.c_final) / (d.m - d.LBs) * 100, 100.0)
d.loc[d.LBs.isna(), "eff"] = np.nan
d["proven_optimal"] = d.c_final == d.LBs
d.to_csv(a.out, index=False)

pd.set_option("display.width", 200)
print(d.groupby("A_src").size().to_string(), "\n")
best = d.loc[d.dropna(subset=["eff"]).groupby(["dataset", "shape"]).c_final.idxmin()]
t = d.dropna(subset=["eff"]).groupby("dataset").agg(
    runs=("eff", "size"), proven_opt=("proven_optimal", "sum"),
    eff_median=("eff", "median"), eff_min=("eff", "min"))
t = t.join(best.groupby("dataset").eff.agg(best_median="median", best_min="min", best_max="max"))
print(t.round(1).to_string())
print(f"\nper-run table: {a.out}")
