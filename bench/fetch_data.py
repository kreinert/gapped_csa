#!/usr/bin/env python3
"""Fetch (and verify the cache of) every "kind": "fetched" dataset in
datasets.py.

Nothing this script writes is committed: everything lands under --data-dir
(see config.py for how that's resolved -- gapped_csa/.gitignore also
excludes *.fasta repo-wide as a second line of defense). Re-running is a
no-op once a dataset is cached; use --force to re-download.

Usage:
  ./fetch_data.py                        # fetch/verify everything
  ./fetch_data.py --only ecoli_k12       # just one
  ./fetch_data.py --print-hash ecoli_k12 # fetch, then print the sha256 to
                                          # pin in datasets.py
"""
import argparse
import gzip
import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

import config
from datasets import DATASETS


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_one(ds: dict, data_dir: Path, force: bool) -> Path:
    out = data_dir / f"{ds['name']}.fasta"
    if out.exists() and not force:
        print(f"[skip] {ds['name']}: already cached at {out}")
        return out

    data_dir.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".fasta.download")
    print(f"[fetch] {ds['name']} <- {ds['url']}")
    try:
        urllib.request.urlretrieve(ds["url"], tmp)
    except Exception as e:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"[error] {ds['name']}: download failed: {e}")

    # Ensembl-style URLs end in .fa.gz; NCBI eutils returns plain text.
    if ds["url"].endswith(".gz"):
        with gzip.open(tmp, "rb") as fin, open(out, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        tmp.unlink()
    else:
        tmp.rename(out)

    digest = sha256_of(out)
    expected = ds.get("sha256")
    if expected is None:
        print(f"  sha256={digest}")
        print(f"  (not pinned yet -- copy this into datasets.py as this "
              f"entry's sha256 so future fetches are verified)")
    elif digest != expected:
        out.unlink()
        raise SystemExit(
            f"[error] {ds['name']}: checksum mismatch (got {digest}, "
            f"expected {expected}). Deleted the download -- do not use it."
        )
    else:
        print(f"  sha256 OK ({digest[:12]}...)")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=None,
                     help="where fetched FASTAs are cached (default: see config.py)")
    ap.add_argument("--only", nargs="*", metavar="NAME",
                     help="only fetch these dataset names")
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    ap.add_argument("--print-hash", metavar="NAME",
                     help="fetch (or use the cache) and print NAME's sha256, then exit")
    args = ap.parse_args()

    data_dir = config.resolve("DATA_DIR", args.data_dir)
    todo = [d for d in DATASETS if d["kind"] == "fetched"]

    if args.print_hash:
        ds = next((d for d in todo if d["name"] == args.print_hash), None)
        if ds is None:
            sys.exit(f"no 'fetched' dataset named {args.print_hash!r}")
        out = fetch_one(ds, data_dir, force=False)
        print(sha256_of(out))
        return

    if args.only:
        missing = set(args.only) - {d["name"] for d in todo}
        if missing:
            sys.exit(f"not 'fetched' datasets: {sorted(missing)}")
        todo = [d for d in todo if d["name"] in args.only]

    if not todo:
        print("nothing to fetch (no 'fetched'-kind datasets matched)")
        return

    for ds in todo:
        fetch_one(ds, data_dir, args.force)


if __name__ == "__main__":
    main()
