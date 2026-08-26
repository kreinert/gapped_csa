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
  ./run_suite.py --disable-phase2                  # GCSA_DISABLE_PHASE2=1 for
                                                    # every run (Phase I /
                                                    # leftover only)
  ./run_suite.py --env GCSA_PHASE2_MAX_ITERS=50 GCSA_QUIET=1
                                                    # arbitrary extra env vars
                                                    # for every ./gcsa call --
                                                    # see compress.hpp for the
                                                    # full GCSA_PHASE2_* knob
                                                    # list. See --env's and
                                                    # --disable-phase2's --help
                                                    # for a real gotcha:
                                                    # GCSA_DISABLE_PHASE2 is
                                                    # presence-checked, not
                                                    # value-checked, so --env
                                                    # can't force it back off.

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
                                                    # (see README). Splits by
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
from datasets import ALGOS, DATASETS, MAX_ADDS, SHAPES, SKIP_SELFTEST_ABOVE_BYTES

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


def log_filename(dataset: str, shape: str, algo: str, max_add: int,
                  disable_phase2: bool) -> str:
    shape_part = shape if _SHAPE_FILENAME_SAFE.match(shape) else re.sub(r"[^\w.-]", "_", shape)
    suffix = "__nophase2" if disable_phase2 else ""
    return f"{dataset}__{shape_part}__{algo}__ma{max_add}{suffix}.log"

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
    referenced via "@name" by others (e.g. "@ecoli_k12", used by both the
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


def already_done(out_csv: Path) -> set:
    if not out_csv.exists():
        return set()
    done = set()
    with open(out_csv, newline="") as f:
        for row in csv.DictReader(f):
            # row.get("phase2", "on"): a CSV written before --disable-phase2
            # existed has no "phase2" column at all -- treat every row in it
            # as a phase2-on run, which is what it actually was.
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
                    phase2_label, skip_selftest, disable_phase2, extra_env,
                    log_dir):
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
    # then --disable-phase2, then --env -- each step can override the one
    # before it, so --env is the final word if it names the same var.
    overrides = {}
    if skip_selftest:
        overrides["GCSA_SKIP_SELFTEST"] = "1"
    if disable_phase2:
        overrides["GCSA_DISABLE_PHASE2"] = "1"
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
        log_path = log_dir / log_filename(ds["name"], shape, algo, max_add, disable_phase2)
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
    ap.add_argument("--disable-phase2", action="store_true",
                     help="set GCSA_DISABLE_PHASE2=1 for every ./gcsa invocation (Phase I "
                          "/ leftover-DP output only, no Phase II local-search pass). "
                          "Recorded in the 'phase2' CSV column so phase2-on and "
                          "phase2-off rows for the same (dataset, shape, algo, max_add) "
                          "can coexist in one --out without colliding. This is the flag "
                          "to use for that specific variable -- setting it instead via "
                          "'--env GCSA_DISABLE_PHASE2=1' has the same effect on ./gcsa but "
                          "will NOT be reflected in the 'phase2' column or the log filename. "
                          "Note ./gcsa checks this var for *presence*, not value (see "
                          "compress.hpp) -- GCSA_DISABLE_PHASE2=0 still disables Phase II, "
                          "so --env can't be used to force it back on; just omit this flag.")
    ap.add_argument("--env", nargs="*", default=[], metavar="KEY=VALUE",
                     help="extra environment variables for every ./gcsa invocation, e.g. "
                          "--env GCSA_PHASE2_MAX_ITERS=50 GCSA_QUIET=1 (see compress.hpp "
                          "for the full GCSA_PHASE2_* / GCSA_QUIET knob list). Applied "
                          "after, and so overrides, the automatic GCSA_SKIP_SELFTEST if "
                          "you name that same variable -- but NOT a reliable way to "
                          "override --disable-phase2 back off (see that flag's help). Not "
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
                     help="process only the datasets at index i where i %% N == I (0-based "
                          "I, e.g. --shard 0/4 .. --shard 3/4 for four cluster array "
                          "tasks). Splits by dataset, not by row, so each dataset's "
                          "fetch/generate + gzip cost is paid once total, not once per "
                          "shard -- combine with --jobs for finer-grained parallelism "
                          "within one shard. Filters AFTER --only-category/--only-dataset, "
                          "over whatever's left. Each shard needs its own --out (see "
                          "--out's default); concatenate the resulting CSVs afterward "
                          "(same header on every shard, so `head -1 s0.csv > merged.csv "
                          "&& tail -n +2 -q s*.csv >> merged.csv` works) -- see README.")
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
    datasets = DATASETS
    if args.only_category:
        cats = set(args.only_category)
        datasets = [d for d in datasets if d["category"] in cats]
    if args.only_dataset:
        names = set(args.only_dataset)
        datasets = [d for d in datasets if d["name"] in names]
    if args.shard is not None:
        assert shard_i is not None and shard_n is not None  # set together, above
        datasets = [d for idx, d in enumerate(datasets) if idx % shard_n == shard_i]

    def log(msg):
        print(msg, file=sys.stderr)

    phase2_label = "off" if args.disable_phase2 else "on"
    log(f"bin-dir={bin_dir}  data-dir={data_dir}  tmp-dir={tmp_dir}  out={out_csv}  "
        f"(run 'python3 config.py' to see where each came from)")
    log(f"phase2={phase2_label}"
        + (f"  log-dir={log_dir}" if log_dir is not None else "  log-dir=(disabled)")
        + f"  jobs={args.jobs}"
        + (f"  shard={shard_i}/{shard_n} ({len(datasets)} dataset(s))" if args.shard is not None else ""))
    if extra_env:
        log("extra env: " + " ".join(f"{k}={v}" for k, v in extra_env.items())
            + "  (also GCSA_SKIP_SELFTEST=1 automatically on inputs over "
              f"{SKIP_SELFTEST_ABOVE_BYTES} bytes, unless overridden above)")

    if args.dry_run:
        print(f"{len(datasets)} datasets x {len(args.shapes)} shapes x "
              f"{len(args.algos)} algos x {len(args.max_adds)} max_adds = "
              f"{len(datasets) * len(args.shapes) * len(args.algos) * len(args.max_adds)} rows")
        for d in datasets:
            print(f"  {d['name']:<28} category={d['category']:<12} kind={d['kind']}")
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
                                         args.disable_phase2, extra_env, log_dir))
            else:
                # Concurrent subprocess dispatch; all file writes (CSV row,
                # log file) still happen only here on the main thread, in
                # completion order -- no locking needed since nothing else
                # ever touches fcsv/log files.
                with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
                    futures = [
                        pool.submit(run_experiment, bin_dir, path, ds, raw, gz, ratio,
                                    shape, algo, max_add, phase2_label, skip_selftest,
                                    args.disable_phase2, extra_env, log_dir)
                        for shape, algo, max_add in todo_jobs
                    ]
                    for fut in concurrent.futures.as_completed(futures):
                        emit(fut.result())

            if by_name[ds["name"]]["kind"] in ("synthetic", "concat") and not args.keep_synthetic:
                path.unlink(missing_ok=True)

    print(f"wrote {out_csv}")


if __name__ == "__main__":
    main()
