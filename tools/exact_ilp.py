#!/usr/bin/env python3
"""Exact |C| for one (text, shape, max-add) instance via a *compact* ILP.

Reads the run dump written by `./lower_bound --dump-runs FILE` and solves the
same problem as ilp_baseline (the full link universe of compress.hpp: every
contiguous, clash-free sub-run of length >= kMinCoverage of every maximal
source run, adds 1..A, one link per k-mer), but without enumerating the
O(L^2) sub-runs.  Per maximal run R (fixed add, fixed target word) of length L:

    u_R          R is the run this word's single link reads from
    s_{R,t}      rank lo_R + t is a source of that link   (covers tgt_R[t])
    b_{R,t}      the source block starts at t             (contiguity)

    sum_{R -> c} u_R <= 1                    one H_offset per word      (L5)
    s_{R,t} <= u_R ;  sum_t b_{R,t} <= u_R ;  b_{R,t} >= s_{R,t} - s_{R,t-1}
    sum_t s_{R,t} >= kMinCoverage * u_R                                  (L4)
    s_{R,t} <= x_{lo_R+t}                    sources are kept            (L1)
    x_r + sum_{(R,t): tgt_R[t] = r} s_{R,t} <= 1   covered => dropped    (L3)

    maximise sum s  (|C| = m - sum s)

A clashing sub-run (one that covers one of its own sources) is infeasible
here automatically: its source would need x=1 and x=0 at once.  Variables
O(sum L) instead of O(sum L^2) candidates and O(sum L^3) pin constraints.

Usage:  exact_ilp.py RUNS.bin [--solver cpsat|highs] [--time 600] [--threads N]
Prints one key=value line: opt (or best/bound if not proven), status, time.
"""
import argparse
import struct
import sys
import time

import numpy as np


def read_dump(path):
    with open(path, "rb") as f:
        m, K, A = struct.unpack("<3q", f.read(24))
        runs = []
        while True:
            a, w, s, L = struct.unpack("<4q", f.read(32))
            if a == 0:
                break
            tg = np.frombuffer(f.read(4 * L), dtype=np.int32).copy()
            runs.append((a, w, s, tg))
        word_of_rank = np.frombuffer(f.read(4 * m), dtype=np.int32)
    return m, K, A, runs, word_of_rank


def solve_cpsat(m, runs, min_cov, time_limit, threads, hint=None):
    from ortools.sat.python import cp_model
    md = cp_model.CpModel()
    by_word = {}
    cover_terms = {}   # rank -> list of s vars covering it
    src_terms = {}     # rank -> list of s vars using it as a source
    all_s = []
    for (a, w, s, tg) in runs:
        L = len(tg)
        u = md.NewBoolVar("")
        by_word.setdefault(w, []).append(u)
        sv = [md.NewBoolVar("") for _ in range(L)]
        bv = [md.NewBoolVar("") for _ in range(L)]
        for t in range(L):
            md.AddImplication(sv[t], u)
            if t == 0:
                md.AddImplication(sv[0], bv[0])
            else:
                md.Add(bv[t] >= sv[t] - sv[t - 1])
        md.Add(sum(bv) <= u)
        md.Add(sum(sv) >= min_cov * u)
        for t in range(L):
            cover_terms.setdefault(int(tg[t]), []).append(sv[t])
            src_terms.setdefault(s + t, []).append(sv[t])
        all_s.extend(sv)
    for w, us in by_word.items():
        if len(us) > 1:
            md.AddAtMostOne(us)
    for r, cov in cover_terms.items():
        srcs = src_terms.get(r)
        if not srcs:
            if len(cov) > 1:
                md.AddAtMostOne(cov)
            continue
        x = md.NewBoolVar("")
        for sv in srcs:
            md.AddImplication(sv, x)
        md.Add(x + sum(cov) <= 1)
    md.Maximize(sum(all_s))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_workers = threads
    t0 = time.time()
    st = solver.Solve(md)
    dt = time.time() - t0
    name = solver.StatusName(st)
    best = m - int(round(solver.ObjectiveValue())) if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    bound = m - int(np.floor(solver.BestObjectiveBound() + 1e-6))
    return dict(status=name, best=best, bound=bound, secs=dt,
                nvars=len(all_s), proven=(st == cp_model.OPTIMAL))


def solve_highs(m, runs, min_cov, time_limit, threads):
    import highspy
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", float(time_limit))
    h.setOptionValue("threads", threads)
    inf = highspy.kHighsInf
    col_cost, col_lo, col_hi = [], [], []
    rows = []   # (lo, hi, [(col, coef)])

    def var(cost=0.0):
        col_cost.append(cost); col_lo.append(0.0); col_hi.append(1.0)
        return len(col_cost) - 1

    by_word, cover_terms, src_terms = {}, {}, {}
    for (a, w, s, tg) in runs:
        L = len(tg)
        u = var()
        by_word.setdefault(w, []).append(u)
        sv = [var(-1.0) for _ in range(L)]
        bv = [var() for _ in range(L)]
        for t in range(L):
            rows.append((-inf, 0.0, [(sv[t], 1), (u, -1)]))
            if t == 0:
                rows.append((-inf, 0.0, [(sv[0], 1), (bv[0], -1)]))
            else:
                rows.append((-inf, 0.0, [(sv[t], 1), (sv[t - 1], -1), (bv[t], -1)]))
        rows.append((-inf, 0.0, [(b, 1) for b in bv] + [(u, -1)]))
        rows.append((0.0, inf, [(x, 1) for x in sv] + [(u, -min_cov)]))
        for t in range(L):
            cover_terms.setdefault(int(tg[t]), []).append(sv[t])
            src_terms.setdefault(s + t, []).append(sv[t])
    for w, us in by_word.items():
        if len(us) > 1:
            rows.append((-inf, 1.0, [(u, 1) for u in us]))
    for r, cov in cover_terms.items():
        srcs = src_terms.get(r)
        if not srcs:
            if len(cov) > 1:
                rows.append((-inf, 1.0, [(c, 1) for c in cov]))
            continue
        x = var()
        for sv_ in srcs:
            rows.append((-inf, 0.0, [(sv_, 1), (x, -1)]))
        rows.append((-inf, 1.0, [(x, 1)] + [(c, 1) for c in cov]))
    n = len(col_cost)
    lp = highspy.HighsLp()
    lp.num_col_ = n
    lp.num_row_ = len(rows)
    lp.col_cost_ = np.array(col_cost)
    lp.col_lower_ = np.array(col_lo)
    lp.col_upper_ = np.array(col_hi)
    lp.row_lower_ = np.array([r[0] for r in rows])
    lp.row_upper_ = np.array([r[1] for r in rows])
    starts, idx, val = [], [], []
    for r in rows:
        starts.append(len(idx))
        for c, v in r[2]:
            idx.append(c); val.append(float(v))
    lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
    lp.a_matrix_.start_ = np.array(starts + [len(idx)], dtype=np.int32)
    lp.a_matrix_.index_ = np.array(idx, dtype=np.int32)
    lp.a_matrix_.value_ = np.array(val)
    lp.integrality_ = [highspy.HighsVarType.kInteger] * n
    h.passModel(lp)
    t0 = time.time()
    h.run()
    dt = time.time() - t0
    info = h.getInfo()
    st = h.modelStatusToString(h.getModelStatus())
    best = m + int(round(info.objective_function_value)) if info.primal_solution_status == 2 else None
    bound = m + int(np.ceil(info.mip_dual_bound - 1e-6))
    return dict(status=st, best=best, bound=bound, secs=dt, nvars=n,
                proven=(st == "Optimal"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--solver", default="cpsat", choices=["cpsat", "highs"])
    ap.add_argument("--time", type=float, default=600)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--min-coverage", type=int, default=3)
    a = ap.parse_args()
    m, K, A, runs, _ = read_dump(a.dump)
    fn = solve_cpsat if a.solver == "cpsat" else solve_highs
    if not runs:   # nothing is coverable: |C| = m
        r = dict(status="Optimal(empty)", best=m, bound=m, secs=0.0, nvars=0, proven=True)
    else:
        r = fn(m, runs, a.min_coverage, a.time, a.threads)
    r.update(m=m, max_add=A, runs=len(runs), solver=a.solver)
    print(" ".join(f"{k}={v}" for k, v in r.items()))


if __name__ == "__main__":
    main()
