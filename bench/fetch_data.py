#!/usr/bin/env python3
"""Fetch (and verify the cache of) every "kind": "fetched" dataset in
datasets.py.

Nothing this script writes is committed: everything lands under --data-dir
(see config.py for how that's resolved -- gapped_csa/.gitignore also
excludes *.fasta repo-wide as a second line of defense).

Re-running does NOT re-download once a dataset is cached, but it DOES
re-verify: the cached file's sha256 is recomputed and checked against
whatever's pinned in datasets.py, every time. That's the answer to "how do
I check the hashes I just pinned are right" -- just run this again with no
flags. Use --force to actually re-download instead of only re-hashing.

Usage:
  ./fetch_data.py                        # verify everything already cached
                                          # (fetch whatever isn't cached yet)
  ./fetch_data.py --only ecoli_k12       # just one
  ./fetch_data.py --print-hash ecoli_k12 # fetch/verify, then print the
                                          # sha256 to pin in datasets.py
  ./fetch_data.py --force                # re-download + re-verify everything
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


def _check_digest(name: str, out: Path, digest: str, expected, *, just_downloaded: bool) -> bool:
    """Print a verdict for one file's hash and return whether it's OK to
    trust (True = matches or nothing to check against; False = mismatch)."""
    if expected is None:
        print(f"  sha256={digest}")
        print(f"  (not pinned yet -- copy this into datasets.py as this "
              f"entry's sha256 so future runs verify it)")
        return True
    if digest == expected:
        print(f"  sha256 OK ({digest[:12]}...)")
        return True
    if just_downloaded:
        out.unlink()
        raise SystemExit(
            f"[error] {name}: checksum mismatch (got {digest}, "
            f"expected {expected}). Deleted the download -- do not use it."
        )
    print(f"[MISMATCH] {name}: cached file at {out} does not match the "
          f"pinned sha256\n  expected {expected}\n  got      {digest}")
    print(f"  Not deleting a pre-existing cached file automatically -- this "
          f"is as likely to be a copy/paste slip in datasets.py as a bad "
          f"download. Compare by hand before trusting either. If the file "
          f"really is bad, rerun with --force to re-download it.")
    return False


def fetch_one(ds: dict, data_dir: Path, force: bool) -> "tuple[Path, bool]":
    out = data_dir / f"{ds['name']}.fasta"
    if out.exists() and not force:
        print(f"[cached] {ds['name']}: {out} -- re-hashing to verify")
        digest = sha256_of(out)
        ok = _check_digest(ds["name"], out, digest, ds.get("sha256"), just_downloaded=False)
        return out, ok

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
    ok = _check_digest(ds["name"], out, digest, ds.get("sha256"), just_downloaded=True)
    return out, ok


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
        out, ok = fetch_one(ds, data_dir, force=False)
        print(sha256_of(out))
        if not ok:
            sys.exit(1)
        return

    if args.only:
        missing = set(args.only) - {d["name"] for d in todo}
        if missing:
            sys.exit(f"not 'fetched' datasets: {sorted(missing)}")
        todo = [d for d in todo if d["name"] in args.only]

    if not todo:
        print("nothing to fetch (no 'fetched'-kind datasets matched)")
        return

    results = [fetch_one(ds, data_dir, args.force) for ds in todo]
    failed = [ds["name"] for ds, (_, ok) in zip(todo, results) if not ok]
    if failed:
        print(f"\n[FAILED] {len(failed)} dataset(s) did not verify: {failed}")
        sys.exit(1)
    print(f"\nall {len(results)} dataset(s) verified OK")


if __name__ == "__main__":
    main()
