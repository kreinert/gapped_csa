#!/usr/bin/env python3
"""Run the full benchmark matrix (datasets x shapes x algorithms x max_add)
through ./gcsa and write one combined CSV.

For each dataset in datasets.py this:
  1. resolves a local FASTA path -- generating it (synthetic), expecting it
     in --data-dir (fetched/provided), or skipping with a clear reason if
     it isn't available yet;
  2. computes that file's gzip ratio (raw bytes / gzip -9 bytes) as the
     input-compressibility baseline from the design doc;
  3. runs every (shape, algorithm, max_add) combination through ./gcsa and
     parses its stdout report into one CSV row.

Every row is appended and flushed immediately, so a run that's interrupted
or crashes partway through (a real risk once real-genome-scale inputs are
in the matrix) still leaves a usable partial CSV -- rerun and it just adds
the remaining rows (use --resume to skip combos already present in --out).

Usage:
  ./run_suite.py                                  # everything in datasets.py
  ./run_suite.py --only-category random repetitive pangenome
  ./run_suite.py --only-dataset repeat_50pct_x50_y200
  ./run_suite.py --dry-run                        # show the plan, run nothing
  ./run_suite.py --resume                         # skip rows already in --out
  ./run_suite.py --log-dir logs                    # keep every ./gcsa
                                                    # invocation's full
                                                    # stdout+stderr for
                                                    # inspection (see --out's
                                                    # log_path column)
  ./run_suite.py --disable-retarget                # skip the unpin/retarget
                                                    # sweep (the LNS still
                                                    # runs)
  ./run_suite.py --disable-lns                     # skip the cluster LNS (the
                                                    # retarget sweep still
                                                    # runs)
  ./run_suite.py --disable-retarget --disable-lns  # Phase I / leftover only.
                                                    # All four combinations
                                                    # land in the 'phase2'
                                                    # column as on /
                                                    # no-retarget / no-lns /
                                                    # off
  ./run_suite.py --env GCSA_PHASE2_MAX_ITERS=50 GCSA_QUIET=1
                                                    # arbitrary extra env vars
                                                    # for every ./gcsa call --
                                                    # see compress.hpp for the
                                                    # full GCSA_PHASE2_* knob
                                                    # list.

Two ways to run this concurrently on a cluster -- combine them freely:

  ./run_suite.py --jobs 8                          # up to 8 ./gcsa processes
                                                    # at once within this one
                                                    # run_suite.py process
                                                    # (match to your allocated
                                                    # CPUs -- each ./gcsa is
                                                    # its own OS process, so
                                                    # this multiplies real CPU
                                                    # usage, not just
                                                    # wall-clock).
  ./run_suite.py --shard 0/4 --out results/s0.csv  # this is shard 0 of 4 --
  ./run_suite.py --shard 1/4 --out results/s1.csv  # run the other 3 as
  ...                                               # separate cluster
                                                    # jobs/array tasks, then
                                                    # concatenate the CSVs
                                                    # (see README). Shards
                                                    # are balanced by
                                                    # estimated input size
                                                    # (fetch_data.py first
                                                    # for accurate weights),
                                                    # not split round-robin
                                                    # by index, so one
                                                    # shard doesn't get
                                                    # stuck with every big
                                                    # dataset. Splits by
                                                    # dataset (not by row) so
                                                    # each dataset's fetch/
                                                    # generate + gzip cost is
                                                    # paid once total, not
                                                    # once per shard.
  # a SLURM array job putting both together:
  #   --shard $SLURM_ARRAY_TASK_ID/$SLURM_ARRAY_TASK_COUNT --jobs $SLURM_CPUS_PER_TASK
"""
import argparse
import concurrent.futures
import csv
import gzip
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import config
from datasets import (ALGOS, DATASETS, INGREDIENT_ONLY_CATEGORIES, MAX_ADDS,
                      SHAPES, SKIP_SELFTEST_ABOVE_BYTES)

CSV_FIELDS = [
    "dataset", "category", "path", "raw_bytes", "gzip_bytes", "gzip_ratio",
    "shape", "span", "weight", "algo", "max_add", "phase2",
    "distinct_kmers", "m", "C", "keep_pct",
    "bytes_total", "full_sa_bytes", "size_pct",
    "self_test", "roundtrip", "wall_ms", "status", "log_path",
]

# Characters a shape string may legitimately contain (see main.cpp) --
# used only to sanity-check before dropping a shape straight into a log
# filename; SHAPES in datasets.py never produces anything outside this.
_SHAPE_FILENAME_SAFE = re.compile(r"^[#.]+$")


# The two Phase II passes -- the greedy unpin/retarget sweep and the exact
# cluster LNS -- switch on and off independently, so a run is one of four
# configurations rather than "phase2 on/off". One table maps that choice to
# the CSV label, the env ./gcsa actually sees, and the log-filename suffix, so
# the three can never drift apart.
#
# Note for anyone comparing against a CSV written before this existed: the old
# --disable-phase2 set GCSA_DISABLE_PHASE2, which only ever skipped the
# retarget pass -- the LNS still ran. So a historical phase2="off" row is what
# is now labelled "no-retarget", not "off".
PHASE2_LABEL = {(False, False): "on",     (True, False): "no-retarget",
                 (False, True):  "no-lns", (True, True):  "off"}
PHASE2_ENV = {
    "on":          {},
    "no-retarget": {"GCSA_DISABLE_RETARGET": "1"},
    "no-lns":      {"GCSA_DISABLE_LNS": "1"},
    "off":         {"GCSA_DISABLE_RETARGET": "1", "GCSA_DISABLE_LNS": "1"},
}
PHASE2_SUFFIX = {"on": "", "no-retarget": "__noretarget",
                 "no-lns": "__nolns", "off": "__nophase2"}


def log_filename(dataset: str, shape: str, algo: str, max_add: int,
                  phase2_label: str) -> str:
    shape_part = shape if _SHAPE_FILENAME_SAFE.match(shape) else re.sub(r"[^\w.-]", "_", shape)
    return f"{dataset}__{shape_part}__{algo}__ma{max_add}{PHASE2_SUFFIX[phase2_label]}.log"

GCSA_RE = {
    "span": re.compile(r"\bspan=(\d+)"),
    "weight": re.compile(r"\bweight=(\d+)"),
    "distinct_kmers": re.compile(r"distinct k-mers\s*:\s*(\d+)"),
    "m": re.compile(r"SA entries \(m\)\s*:\s*(\d+)"),
    "c_keep": re.compile(r"stored positions \(\|C\|\):\s*(\d+)\s*\(([\d.]+)% of full SA\)"),
    "self_test": re.compile(r"self-test \(vs brute force\):\s*(\S+)"),
    "roundtrip": re.compile(r"serialized round-trip\s*:\s*(\S+)"),
    "size": re.compile(r"total=(\d+)\s*\(full SA would be (\d+) bytes; ([\d.]+)%\)"),
}


def parse_gcsa_stdout(text: str) -> dict:
    out = {}
    if (m := GCSA_RE["span"].search(text)): out["span"] = int(m.group(1))
    if (m := GCSA_RE["weight"].search(text)): out["weight"] = int(m.group(1))
    if (m := GCSA_RE["distinct_kmers"].search(text)): out["distinct_kmers"] = int(m.group(1))
    if (m := GCSA_RE["m"].search(text)): out["m"] = int(m.group(1))
    if (m := GCSA_RE["c_keep"].search(text)):
        out["C"] = int(m.group(1)); out["keep_pct"] = float(m.group(2))
    if (m := GCSA_RE["self_test"].search(text)): out["self_test"] = m.group(1)
    if (m := GCSA_RE["roundtrip"].search(text)): out["roundtrip"] = m.group(1)
    if (m := GCSA_RE["size"].search(text)):
        out["bytes_total"] = int(m.group(1))
        out["full_sa_bytes"] = int(m.group(2))
        out["size_pct"] = float(m.group(3))
    return out


class _CountingWriter:
    """Discards bytes but counts them -- lets us stream gzip.GzipFile over a
    large FASTA without holding the compressed (or raw) output in memory."""
    def __init__(self):
        self.n = 0
    def write(self, b):
        self.n += len(b)
        return len(b)
    def flush(self):
        pass


def gzip_ratio(path: Path) -> "tuple[int, int, float]":
    raw = path.stat().st_size
    counter = _CountingWriter()
    with open(path, "rb") as fin, gzip.GzipFile(fileobj=counter, mode="wb",
                                                 compresslevel=9) as gz:
        while chunk := fin.read(1 << 20):
            gz.write(chunk)
    gz_bytes = counter.n
    ratio = raw / gz_bytes if gz_bytes else float("inf")
    return raw, gz_bytes, ratio


def resolve_path(name: str, by_name: dict, bin_dir: Path, data_dir: Path,
                  tmp_dir: Path, cache: dict, log) -> "Path | None":
    """Resolve one dataset name to a local FASTA path, generating/checking/
    concatenating as its `kind` requires. Memoized in `cache` so a dataset
    referenced via "@name" by others (e.g. "@ecoli_003", used by both the
    simulated pangenome sweep and the real-strain concat entry) is only
    resolved once per run."""
    if name in cache:
        return cache[name]
    ds = by_name.get(name)
    if ds is None:
        log(f"[skip] unknown dataset referenced: {name!r}")
        cache[name] = None
        return None

    path = None
    if ds["kind"] == "provided":
        candidate = data_dir / ds["path"]
        if candidate.exists():
            path = candidate
        else:
            log(f"[skip] {name}: not found at {candidate} "
                f"(expected to already be on disk -- see datasets.py)")

    elif ds["kind"] == "fetched":
        candidate = data_dir / f"{name}.fasta"
        if candidate.exists():
            path = candidate
        else:
            log(f"[skip] {name}: not fetched yet -- run "
                f"'./fetch_data.py --only {name}' first")

    elif ds["kind"] == "synthetic":
        gen_bin = bin_dir / ds["generator"]
        if not gen_bin.exists():
            log(f"[skip] {name}: generator binary {gen_bin} not built "
                f"(run `make {ds['generator']}` in gapped_csa/)")
        else:
            # Resolve any "@other_name" argument references first.
            resolved_args = []
            ok = True
            for a in ds["args"]:
                if isinstance(a, str) and a.startswith("@"):
                    ref = resolve_path(a[1:], by_name, bin_dir, data_dir, tmp_dir, cache, log)
                    if ref is None:
                        log(f"[skip] {name}: reference {a} could not be resolved")
                        ok = False
                        break
                    resolved_args.append(str(ref))
                else:
                    resolved_args.append(a)
            if ok:
                out_path = tmp_dir / f"{name}.fasta"
                cmd = [str(gen_bin), *resolved_args, "-o", str(out_path)]
                r = subprocess.run(cmd, capture_output=True, text=True)
                if r.returncode != 0:
                    log(f"[skip] {name}: generator failed: {r.stderr.strip()}")
                else:
                    path = out_path

    elif ds["kind"] == "concat":
        parts = []
        ok = True
        for ref in ds["refs"]:
            ref_name = ref[1:] if isinstance(ref, str) and ref.startswith("@") else ref
            p = resolve_path(ref_name, by_name, bin_dir, data_dir, tmp_dir, cache, log)
            if p is None:
                log(f"[skip] {name}: reference {ref} could not be resolved")
                ok = False
                break
            parts.append(p)
        if ok:
            out_path = tmp_dir / f"{name}.fasta"
            with open(out_path, "wb") as fout:
                for p in parts:
                    with open(p, "rb") as fin:
                        fout.write(fin.read())
                    fout.write(b"\n")  # guard against a missing trailing newline
            path = out_path

    else:
        log(f"[skip] {name}: unknown kind {ds['kind']!r}")

    cache[name] = path
    return path


def estimate_input_bytes(name: str, by_name: dict, data_dir: Path, cache: dict) -> int:
    """Estimate a dataset's resolved FASTA size in bytes, WITHOUT actually
    generating/fetching/concatenating it -- used to balance --shard
    partitions by workload rather than by raw dataset count (a
    pangenome_ecoli_real_n100-sized input takes vastly longer at every
    shape x algo x max_add combination than random_1e6 does, so splitting
    datasets round-robin by index can strand one cluster task with every
    big dataset while the others finish early and sit idle).

    "fetched"/"provided" sizes come from stat()-ing the cached/expected
    file -- 0 if it hasn't been fetched yet, so balancing is only as good
    as what's already on disk; run ./fetch_data.py first for accurate
    weights. "synthetic" sizes are computed directly from the generator's
    args (no need to actually run the generator). "concat" sizes are the
    sum of their resolved refs' sizes, recursing through `by_name` the same
    way resolve_path's own "@name" handling does. Memoized in `cache`,
    mirroring resolve_path's resolved_cache."""
    if name in cache:
        return cache[name]
    ds = by_name.get(name)
    if ds is None:
        cache[name] = 0
        return 0

    def arg_after(args, flag, default=None):
        return args[args.index(flag) + 1] if flag in args else default

    size = 0
    if ds["kind"] == "fetched":
        p = data_dir / f"{name}.fasta"
        size = p.stat().st_size if p.exists() else 0
    elif ds["kind"] == "provided":
        p = data_dir / ds["path"]
        size = p.stat().st_size if p.exists() else 0
    elif ds["kind"] == "synthetic":
        args, gen = ds["args"], ds["generator"]
        if gen == "simulate_random":
            size = int(arg_after(args, "-n", 0))
        elif gen == "simulate_repeats":
            x = int(arg_after(args, "-x", 0))
            y = int(arg_after(args, "-y", 0))
            frac = float(arg_after(args, "--repetitive-frac", 1.0))
            size = round(x * y / frac) if frac else x * y
        elif gen == "simulate_pangenome":
            n = int(arg_after(args, "-n", 1))
            ref = arg_after(args, "-r")
            if ref is not None:
                ref_name = ref[1:] if ref.startswith("@") else ref
                size = estimate_input_bytes(ref_name, by_name, data_dir, cache) * n
            else:
                size = int(arg_after(args, "-y", 0)) * n
        # else: unknown generator -- leave size=0 (under- rather than
        # over-balancing an unrecognized one is the safer failure mode)
    elif ds["kind"] == "concat":
        size = sum(
            estimate_input_bytes(r[1:] if isinstance(r, str) and r.startswith("@") else r,
                                  by_name, data_dir, cache)
            for r in ds["refs"]
        )
    cache[name] = size
    return size


def balance_shards(datasets: list, weights: dict, shard_n: int) -> "tuple[list, list]":
    """Greedy longest-processing-time-first bin packing: process datasets
    heaviest (by `weights`) first, always adding the next one to whichever
    of the `shard_n` bins currently has the smallest total weight.
    Deterministic given the same (datasets, weights, shard_n), so every
    independently-launched `--shard I/N` cluster task computes the
    identical partition on its own -- no coordination between tasks
    needed. Not optimal bin-packing (that's NP-hard) but a well-known
    good-enough heuristic, and simple enough to audit. Degrades to
    (unbalanced, but still deterministic) grouping if every weight is 0 --
    e.g. nothing has been fetched yet, see estimate_input_bytes.
    Returns (bins, bin_totals)."""
    order = sorted(range(len(datasets)), key=lambda i: (-weights[datasets[i]["name"]], i))
    totals = [0] * shard_n
    bins = [[] for _ in range(shard_n)]
    for i in order:
        b = min(range(shard_n), key=lambda b: (totals[b], b))
        bins[b].append(datasets[i])
        totals[b] += weights[datasets[i]["name"]]
    # Re-sort each shard's members back into their original relative order
    # -- purely cosmetic (readable --dry-run/log output), doesn't affect
    # which dataset ended up in which shard.
    order_index = {id(d): idx for idx, d in enumerate(datasets)}
    for b in bins:
        b.sort(key=lambda d: order_index[id(d)])
    return bins, totals


def already_done(out_csv: Path) -> set:
    if not out_csv.exists():
        return set()
    done = set()
    with open(out_csv, newline="") as f:
        for row in csv.DictReader(f):
            # row.get("phase2", "on"): a CSV written before the Phase II
            # switches existed has no "phase2" column at all -- treat every row
            # in it as a both-passes-on run, which is what it actually was. A
            # CSV written by the old --disable-phase2 says "off" where this
            # version would now say "no-retarget"; --resume will re-run those
            # rather than treat them as done, which is the safe direction.
            done.add((row["dataset"], row["shape"], row["algo"], row["max_add"],
                       row.get("phase2") or "on"))
    return done


def check_resumable_schema(out_csv: Path) -> None:
    """--resume appends to an existing CSV with a plain csv.DictWriter, which
    does not re-check the header -- if --out's on-disk header doesn't match
    CSV_FIELDS (e.g. an older suite.csv from before --log-dir/--disable-phase2
    added columns), appending would silently misalign columns. Fail loudly
    instead and say what to do."""
    if not out_csv.exists():
        return
    with open(out_csv, newline="") as f:
        header = next(csv.reader(f), None)
    if header is not None and header != CSV_FIELDS:
        sys.exit(
            f"[error] --resume target {out_csv} has an older column layout "
            f"than this version of run_suite.py expects.\n"
            f"  on disk : {header}\n"
            f"  expected: {CSV_FIELDS}\n"
            f"Point --out at a new file for this run, or rerun without "
            f"--resume to regenerate it from scratch."
        )


def run_experiment(bin_dir, path, ds, raw, gz, ratio, shape, algo, max_add,
                    phase2_label, skip_selftest, extra_env, log_dir):
    """Run one (shape, algo, max_add) experiment for an already-resolved
    dataset and return (row, log_path, log_text, fail_msg) -- pure aside
    from the ./gcsa subprocess call itself, so it's safe to run concurrently
    across a ThreadPoolExecutor. The caller does all file writing (the CSV
    row, the log file) on a single thread; this function never touches a
    file handle shared with anything else."""
    row: "dict[str, str | int | float]" = dict(
        dataset=ds["name"], category=ds["category"],
        path=str(path), raw_bytes=raw, gzip_bytes=gz,
        gzip_ratio=round(ratio, 4), shape=shape, algo=algo,
        max_add=max_add, phase2=phase2_label, status="ok",
        log_path="")
    # Precedence: inherited shell env, then the automatic size-based skip,
    # then the Phase II switches, then --env -- each step can override the one
    # before it, so --env is the final word if it names the same var.
    overrides = {}
    if skip_selftest:
        overrides["GCSA_SKIP_SELFTEST"] = "1"
    overrides.update(PHASE2_ENV[phase2_label])
    overrides.update(extra_env)
    env = dict(os.environ)
    env.update(overrides)

    cmd = [str(bin_dir / "gcsa"), "-g", str(path), "-s", shape,
           "--algo", algo, "--max-add", str(max_add)]
    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    row["wall_ms"] = round((time.perf_counter() - t0) * 1000, 1)

    fail_msg = None
    if r.returncode != 0:
        row["status"] = f"error: {r.stderr.strip()[:200]}"
        fail_msg = f"  [FAIL] {ds['name']} {shape} {algo} ma={max_add}: {row['status']}"
    else:
        row.update(parse_gcsa_stdout(r.stdout))
        bad = row.get("self_test") == "FAIL" or row.get("roundtrip") == "FAIL"
        if bad:
            row["status"] = "correctness_fail"
            fail_msg = (f"  [FAIL] {ds['name']} {shape} {algo} ma={max_add}: "
                        f"self_test={row.get('self_test')} "
                        f"roundtrip={row.get('roundtrip')}")

    log_path, log_text = None, None
    if log_dir is not None:
        log_path = log_dir / log_filename(ds["name"], shape, algo, max_add, phase2_label)
        env_note = " ".join(f"{k}={v}" for k, v in overrides.items()) or "(none)"
        log_text = (
            f"$ {' '.join(cmd)}\n"
            f"env overrides: {env_note}\n"
            f"exit_code: {r.returncode}\n"
            f"wall_ms: {row['wall_ms']}\n"
            f"\n=== stdout ===\n{r.stdout}"
            f"\n=== stderr ===\n{r.stderr}"
        )
        row["log_path"] = str(log_path)

    return row, log_path, log_text, fail_msg


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin-dir", default=None,
                     help="where gcsa/simulate_* live (default: see config.py)")
    ap.add_argument("--data-dir", default=None,
                     help="fetched/provided FASTAs (default: see config.py)")
    ap.add_argument("--tmp-dir", default=None,
                     help="synthetic-input scratch dir (default: see config.py)")
    ap.add_argument("--out", default=None,
                     help="default: results/suite.csv, or results/suite.shard<I>of<N>.csv "
                          "if --shard is given without an explicit --out (so parallel "
                          "shards don't silently share a file by accident).")
    ap.add_argument("--only-category", nargs="*")
    ap.add_argument("--only-dataset", nargs="*")
    ap.add_argument("--shapes", nargs="*", default=SHAPES)
    ap.add_argument("--algos", nargs="*", default=ALGOS)
    ap.add_argument("--max-adds", nargs="*", type=int, default=MAX_ADDS)
    ap.add_argument("--keep-synthetic", action="store_true",
                     help="don't delete generated synthetic FASTAs from --tmp-dir afterward")
    ap.add_argument("--resume", action="store_true",
                     help="skip (dataset, shape, algo, max_add, phase2) rows already in --out")
    ap.add_argument("--dry-run", action="store_true", help="print the plan, run nothing")
    ap.add_argument("--log-dir", default=None,
                     help="write each ./gcsa invocation's full stdout+stderr here, one "
                          "file per (dataset, shape, algo, max_add) -- see log_filename() "
                          "for the naming scheme. Off by default (nothing extra is kept "
                          "beyond the parsed CSV row and, on failure, a 200-char stderr "
                          "snippet in the 'status' column).")
    ap.add_argument("--disable-retarget", action="store_true",
                     help="set GCSA_DISABLE_RETARGET=1 for every ./gcsa invocation, "
                          "skipping the greedy unpin/retarget sweep. The cluster LNS "
                          "still runs -- combine with --disable-lns to skip both and see "
                          "Phase I's output on its own. Recorded in the 'phase2' CSV "
                          "column (on / no-retarget / no-lns / off) and in the log "
                          "filename, so all four configurations for the same (dataset, "
                          "shape, algo, max_add) coexist in one --out without colliding. "
                          "Use this flag rather than '--env GCSA_DISABLE_RETARGET=1': the "
                          "env route has the same effect on ./gcsa but is NOT reflected "
                          "in the column or the filename.")
    ap.add_argument("--disable-lns", action="store_true",
                     help="set GCSA_DISABLE_LNS=1 for every ./gcsa invocation, skipping "
                          "the exact bounded-cluster LNS. The retarget sweep still runs. "
                          "Same CSV/filename bookkeeping as --disable-retarget.")
    ap.add_argument("--disable-phase2", action="store_true",
                     help="DEPRECATED alias for --disable-retarget. Despite the name it "
                          "never disabled all of Phase II: it set GCSA_DISABLE_PHASE2, "
                          "which skips the retarget sweep and leaves the LNS running. "
                          "Kept so existing scripts keep working; new runs should say "
                          "--disable-retarget (and add --disable-lns for what the old "
                          "name suggests).")
    ap.add_argument("--env", nargs="*", default=[], metavar="KEY=VALUE",
                     help="extra environment variables for every ./gcsa invocation, e.g. "
                          "--env GCSA_PHASE2_MAX_ITERS=50 GCSA_QUIET=1 (see compress.hpp "
                          "for the full GCSA_PHASE2_* / GCSA_QUIET knob list). Applied "
                          "after, and so overrides, the automatic GCSA_SKIP_SELFTEST and "
                          "the Phase II switches if you name those same variables -- but "
                          "going that route leaves the 'phase2' CSV column and the log "
                          "filename saying otherwise, so prefer the dedicated flags. Not "
                          "a CSV column (arbitrary keys don't fit a fixed schema); use "
                          "--log-dir to keep a per-row record of exactly what was set.")
    ap.add_argument("--jobs", type=int, default=1, metavar="N",
                     help="run up to N ./gcsa invocations concurrently (across the shape "
                          "x algo x max_add matrix for one dataset at a time -- dataset "
                          "resolution/cleanup itself stays sequential). Match this to the "
                          "CPUs available to this process/job: each ./gcsa invocation is "
                          "its own OS process, so this multiplies real CPU usage rather "
                          "than just hiding I/O wait. Default 1 (sequential, unchanged "
                          "from before this flag existed).")
    ap.add_argument("--shard", default=None, metavar="I/N",
                     help="process only shard I of N (0-based I, e.g. --shard 0/4 .. "
                          "--shard 3/4 for four cluster array tasks). Datasets are "
                          "balanced across the N shards by estimated input size (greedy "
                          "longest-processing-time-first bin packing -- see "
                          "estimate_input_bytes/balance_shards), NOT split round-robin by "
                          "index, so one shard doesn't get stuck with every large dataset "
                          "while the others finish early and idle. Every --shard I/N "
                          "invocation recomputes the same partition independently and "
                          "deterministically -- no coordination between cluster tasks "
                          "needed. Estimates for fetched/provided datasets come from "
                          "stat()-ing the cached file, so run ./fetch_data.py first for "
                          "balancing to reflect real sizes (unfetched ones weigh 0 and "
                          "may land anywhere). Splits by dataset, not by row, so each "
                          "dataset's fetch/generate + gzip cost is paid once total, not "
                          "once per shard -- combine with --jobs for finer-grained "
                          "parallelism within one shard. Filters AFTER --only-category/"
                          "--only-dataset, over whatever's left. Each shard needs its own "
                          "--out (see --out's default); concatenate the resulting CSVs "
                          "afterward (same header on every shard, so `head -1 s0.csv > "
                          "merged.csv && tail -n +2 -q s*.csv >> merged.csv` works) -- "
                          "see README.")
    args = ap.parse_args()

    extra_env = {}
    for kv in args.env:
        if "=" not in kv:
            sys.exit(f"[error] --env entries must look like KEY=VALUE, got {kv!r}")
        k, _, v = kv.partition("=")
        if not k:
            sys.exit(f"[error] --env entry has an empty key: {kv!r}")
        extra_env[k] = v

    if args.jobs < 1:
        sys.exit(f"[error] --jobs must be >= 1, got {args.jobs}")

    shard_i, shard_n = None, None
    if args.shard is not None:
        try:
            i_str, n_str = args.shard.split("/")
            shard_i, shard_n = int(i_str), int(n_str)
        except ValueError:
            sys.exit(f"[error] --shard must look like I/N (e.g. 0/4), got {args.shard!r}")
        if shard_n < 1 or not (0 <= shard_i < shard_n):
            sys.exit(f"[error] --shard {args.shard!r}: need 0 <= I < N with N >= 1")

    out_path_str = args.out
    if out_path_str is None:
        out_path_str = (f"results/suite.shard{shard_i}of{shard_n}.csv"
                         if args.shard is not None else "results/suite.csv")

    bin_dir = config.resolve("BIN_DIR", args.bin_dir)
    data_dir = config.resolve("DATA_DIR", args.data_dir)
    tmp_dir = config.resolve("TMP_DIR", args.tmp_dir)
    out_csv = Path(out_path_str)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    log_dir = Path(args.log_dir) if args.log_dir else None
    if log_dir is not None:
        log_dir.mkdir(parents=True, exist_ok=True)

    by_name = {d["name"]: d for d in DATASETS}

    # INGREDIENT_ONLY_CATEGORIES (real_genome / real_genome_pool -- the
    # E. coli strain pool) are excluded from a default sweep: they're only
    # meant to be resolved via "@name" by the pangenome concat entries, not
    # each run through the full shape x algo x max_add grid on their own
    # (that's what used to add ~100 extra individual-strain experiments to
    # every run). by_name above is built from the unfiltered DATASETS, so
    # "@name" resolution for the concat entries still works regardless.
    # --only-category/--only-dataset explicitly naming one of these
    # categories or datasets overrides the exclusion.
    explicit_categories = set(args.only_category or [])
    explicit_names = set(args.only_dataset or [])
    datasets = [
        d for d in DATASETS
        if d["category"] not in INGREDIENT_ONLY_CATEGORIES
        or d["category"] in explicit_categories
        or d["name"] in explicit_names
    ]
    if args.only_category:
        cats = set(args.only_category)
        datasets = [d for d in datasets if d["category"] in cats]
    if args.only_dataset:
        names = set(args.only_dataset)
        datasets = [d for d in datasets if d["name"] in names]
    # Estimated input size per dataset (bytes), used to balance --shard
    # partitions by workload -- see estimate_input_bytes's docstring for
    # what "estimated" means per dataset kind and its one caveat (fetched/
    # provided datasets not yet on disk weigh 0). Computed unconditionally
    # (cheap -- filesystem stat()s and arithmetic, no subprocess calls) so
    # --dry-run can show it too, not just --shard.
    size_cache = {}
    weights = {d["name"]: estimate_input_bytes(d["name"], by_name, data_dir, size_cache)
               for d in datasets}

    shard_bytes = None
    if args.shard is not None:
        assert shard_i is not None and shard_n is not None  # set together, above
        shard_bins, shard_totals = balance_shards(datasets, weights, shard_n)
        datasets = shard_bins[shard_i]
        shard_bytes = shard_totals[shard_i]

    def log(msg):
        print(msg, file=sys.stderr)

    if args.disable_phase2 and not args.disable_retarget:
        print("note: --disable-phase2 is deprecated; it means --disable-retarget "
              "(the LNS still runs). Add --disable-lns to skip both passes.",
              file=sys.stderr)
    phase2_label = PHASE2_LABEL[(args.disable_retarget or args.disable_phase2,
                                  args.disable_lns)]
    log(f"bin-dir={bin_dir}  data-dir={data_dir}  tmp-dir={tmp_dir}  out={out_csv}  "
        f"(run 'python3 config.py' to see where each came from)")
    shard_log_suffix = ""
    if args.shard is not None:
        assert shard_bytes is not None  # set together with shard_i/shard_n, above
        shard_log_suffix = (f"  shard={shard_i}/{shard_n} ({len(datasets)} dataset(s), "
                             f"~{shard_bytes / 1e6:.1f} MB estimated input -- balanced by "
                             f"size, not by dataset count)")
    log(f"phase2={phase2_label}"
        + (f"  log-dir={log_dir}" if log_dir is not None else "  log-dir=(disabled)")
        + f"  jobs={args.jobs}"
        + shard_log_suffix)
    if extra_env:
        log("extra env: " + " ".join(f"{k}={v}" for k, v in extra_env.items())
            + "  (also GCSA_SKIP_SELFTEST=1 automatically on inputs over "
              f"{SKIP_SELFTEST_ABOVE_BYTES} bytes, unless overridden above)")

    if args.dry_run:
        print(f"{len(datasets)} datasets x {len(args.shapes)} shapes x "
              f"{len(args.algos)} algos x {len(args.max_adds)} max_adds = "
              f"{len(datasets) * len(args.shapes) * len(args.algos) * len(args.max_adds)} rows")
        for d in datasets:
            mb = weights[d["name"]] / 1e6
            print(f"  {d['name']:<28} category={d['category']:<12} kind={d['kind']:<10} "
                  f"~{mb:.2f} MB")
        return

    if args.resume:
        check_resumable_schema(out_csv)
    skip_done = already_done(out_csv) if args.resume else set()
    write_header = not (args.resume and out_csv.exists())
    resolved_cache = {}

    with open(out_csv, "a" if args.resume else "w", newline="") as fcsv:
        writer = csv.DictWriter(fcsv, fieldnames=CSV_FIELDS)
        if write_header:
            writer.writeheader()
            fcsv.flush()

        for ds in datasets:
            path = resolve_path(ds["name"], by_name, bin_dir, data_dir, tmp_dir,
                                 resolved_cache, log)
            if path is None:
                continue

            raw, gz, ratio = gzip_ratio(path)
            log(f"[ok] {ds['name']}: {path}  {raw} bytes, gzip_ratio={ratio:.3f}")
            skip_selftest = raw > SKIP_SELFTEST_ABOVE_BYTES

            todo_jobs = [
                (shape, algo, max_add)
                for shape in args.shapes
                for algo in args.algos
                for max_add in args.max_adds
                if (ds["name"], shape, algo, str(max_add), phase2_label) not in skip_done
            ]

            def emit(result):
                row, log_path, log_text, fail_msg = result
                if fail_msg:
                    log(fail_msg)
                if log_path is not None:
                    log_path.write_text(log_text)
                writer.writerow(row)
                fcsv.flush()

            if args.jobs <= 1:
                # Sequential -- identical behavior to before --jobs existed.
                for shape, algo, max_add in todo_jobs:
                    emit(run_experiment(bin_dir, path, ds, raw, gz, ratio, shape, algo,
                                         max_add, phase2_label, skip_selftest,
                                         extra_env, log_dir))
            else:
                # Concurrent subprocess dispatch; all file writes (CSV row,
                # log file) still happen only here on the main thread, in
                # completion order -- no locking needed since nothing else
                # ever touches fcsv/log files.
                with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
                    futures = [
                        pool.submit(run_experiment, bin_dir, path, ds, raw, gz, ratio,
                                    shape, algo, max_add, phase2_label, skip_selftest,
                                    extra_env, log_dir)
                        for shape, algo, max_add in todo_jobs
                    ]
                    for fut in concurrent.futures.as_completed(futures):
                        emit(fut.result())

            if by_name[ds["name"]]["kind"] in ("synthetic", "concat") and not args.keep_synthetic:
                path.unlink(missing_ok=True)

    print(f"wrote {out_csv}")


if __name__ == "__main__":
    main()
