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
"""
import argparse
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
    "shape", "span", "weight", "algo", "max_add",
    "distinct_kmers", "m", "C", "keep_pct",
    "bytes_total", "full_sa_bytes", "size_pct",
    "self_test", "roundtrip", "wall_ms", "status",
]

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
            done.add((row["dataset"], row["shape"], row["algo"], row["max_add"]))
    return done


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin-dir", default=None,
                     help="where gcsa/simulate_* live (default: see config.py)")
    ap.add_argument("--data-dir", default=None,
                     help="fetched/provided FASTAs (default: see config.py)")
    ap.add_argument("--tmp-dir", default=None,
                     help="synthetic-input scratch dir (default: see config.py)")
    ap.add_argument("--out", default="results/suite.csv")
    ap.add_argument("--only-category", nargs="*")
    ap.add_argument("--only-dataset", nargs="*")
    ap.add_argument("--shapes", nargs="*", default=SHAPES)
    ap.add_argument("--algos", nargs="*", default=ALGOS)
    ap.add_argument("--max-adds", nargs="*", type=int, default=MAX_ADDS)
    ap.add_argument("--keep-synthetic", action="store_true",
                     help="don't delete generated synthetic FASTAs from --tmp-dir afterward")
    ap.add_argument("--resume", action="store_true",
                     help="skip (dataset, shape, algo, max_add) rows already in --out")
    ap.add_argument("--dry-run", action="store_true", help="print the plan, run nothing")
    args = ap.parse_args()

    bin_dir = config.resolve("BIN_DIR", args.bin_dir)
    data_dir = config.resolve("DATA_DIR", args.data_dir)
    tmp_dir = config.resolve("TMP_DIR", args.tmp_dir)
    out_csv = Path(args.out)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    by_name = {d["name"]: d for d in DATASETS}
    datasets = DATASETS
    if args.only_category:
        cats = set(args.only_category)
        datasets = [d for d in datasets if d["category"] in cats]
    if args.only_dataset:
        names = set(args.only_dataset)
        datasets = [d for d in datasets if d["name"] in names]

    def log(msg):
        print(msg, file=sys.stderr)

    log(f"bin-dir={bin_dir}  data-dir={data_dir}  tmp-dir={tmp_dir}  "
        f"(run 'python3 config.py' to see where each came from)")

    if args.dry_run:
        print(f"{len(datasets)} datasets x {len(args.shapes)} shapes x "
              f"{len(args.algos)} algos x {len(args.max_adds)} max_adds = "
              f"{len(datasets) * len(args.shapes) * len(args.algos) * len(args.max_adds)} rows")
        for d in datasets:
            print(f"  {d['name']:<28} category={d['category']:<12} kind={d['kind']}")
        return

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

            for shape in args.shapes:
                for algo in args.algos:
                    for max_add in args.max_adds:
                        key = (ds["name"], shape, algo, str(max_add))
                        if key in skip_done:
                            continue

                        row = dict(dataset=ds["name"], category=ds["category"],
                                   path=str(path), raw_bytes=raw, gzip_bytes=gz,
                                   gzip_ratio=round(ratio, 4), shape=shape, algo=algo,
                                   max_add=max_add, status="ok")
                        env = dict(os.environ)
                        if skip_selftest:
                            env["GCSA_SKIP_SELFTEST"] = "1"

                        cmd = [str(bin_dir / "gcsa"), "-g", str(path), "-s", shape,
                               "--algo", algo, "--max-add", str(max_add)]
                        t0 = time.perf_counter()
                        r = subprocess.run(cmd, capture_output=True, text=True, env=env)
                        row["wall_ms"] = round((time.perf_counter() - t0) * 1000, 1)

                        if r.returncode != 0:
                            row["status"] = f"error: {r.stderr.strip()[:200]}"
                            log(f"  [FAIL] {ds['name']} {shape} {algo} ma={max_add}: "
                                f"{row['status']}")
                        else:
                            row.update(parse_gcsa_stdout(r.stdout))
                            bad = row.get("self_test") == "FAIL" or row.get("roundtrip") == "FAIL"
                            if bad:
                                row["status"] = "correctness_fail"
                                log(f"  [FAIL] {ds['name']} {shape} {algo} ma={max_add}: "
                                    f"self_test={row.get('self_test')} "
                                    f"roundtrip={row.get('roundtrip')}")

                        writer.writerow(row)
                        fcsv.flush()

            if by_name[ds["name"]]["kind"] in ("synthetic", "concat") and not args.keep_synthetic:
                path.unlink(missing_ok=True)

    print(f"wrote {out_csv}")


if __name__ == "__main__":
    main()
