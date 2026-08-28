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

If --data-dir is a network/shared mount, a read can transiently fail with
"Device or resource busy" or similar (a sync client or another process
briefly holding the file, a mount hiccup) -- sha256_of() retries a few
times with backoff before giving up. If it's still unreadable after that,
the dataset is reported "unverified" (distinct from a real hash mismatch)
rather than crashing the whole run; just rerun later to retry it.

Usage:
  ./fetch_data.py                        # verify everything already cached
                                          # (fetch whatever isn't cached yet)
  ./fetch_data.py --only ecoli_003       # just one
  ./fetch_data.py --print-hash ecoli_003 # fetch/verify, then print the
                                          # sha256 to pin in datasets.py
  ./fetch_data.py --force                # re-download + re-verify everything
"""
import argparse
import gzip
import hashlib
import shutil
import sys
import time
import urllib.request
from pathlib import Path

import config
from datasets import DATASETS


def sha256_of(path: Path, retries: int = 4, base_delay: float = 0.5) -> str:
    """Hash `path`, retrying a few times with backoff on a transient OSError
    (EBUSY and friends -- typical of a network/shared mount that's mid-sync
    or briefly locked by another process right after a write). Raises the
    last OSError if it's still unreadable after `retries` attempts."""
    if retries < 1:
        raise ValueError("retries must be >= 1")
    last_err: "OSError | None" = None
    for attempt in range(retries):
        try:
            h = hashlib.sha256()
            with open(path, "rb") as f:
                for chunk in iter(lambda: f.read(1 << 20), b""):
                    h.update(chunk)
            return h.hexdigest()
        except OSError as e:
            last_err = e
            if attempt < retries - 1:
                delay = base_delay * (2 ** attempt)
                reason = e.strerror or str(e)
                print(f"  [retry {attempt + 1}/{retries - 1}] {path} not readable "
                      f"yet ({reason}, errno={e.errno}) -- retrying in {delay:.1f}s",
                      file=sys.stderr)
                time.sleep(delay)
    assert last_err is not None  # retries >= 1 guarantees the loop set this
    raise last_err


def _check_digest(name: str, out: Path, digest: str, expected, *, just_downloaded: bool) -> str:
    """Print a verdict for one file's hash. Returns "ok" (matches, or
    nothing to check against yet) or "mismatch"."""
    if expected is None:
        print(f"  sha256={digest}")
        print(f"  (not pinned yet -- copy this into datasets.py as this "
              f"entry's sha256 so future runs verify it)")
        return "ok"
    if digest == expected:
        print(f"  sha256 OK ({digest[:12]}...)")
        return "ok"
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
    return "mismatch"


def fetch_one(ds: dict, data_dir: Path, force: bool) -> "tuple[Path, str, 'str | None']":
    """Returns (path, status, digest). status is "ok", "mismatch", or
    "unverified" (couldn't even read the file to hash it, after retries --
    see sha256_of's docstring; not the same claim as a hash mismatch)."""
    out = data_dir / f"{ds['name']}.fasta"
    if out.exists() and not force:
        print(f"[cached] {ds['name']}: {out} -- re-hashing to verify")
        try:
            digest = sha256_of(out)
        except OSError as e:
            print(f"[UNVERIFIED] {ds['name']}: could not read {out} to verify "
                  f"({e.strerror or e}) after retries. If {data_dir} is a "
                  f"network/shared mount this is likely transient -- rerun "
                  f"later to retry. Not treated as a hash mismatch.")
            return out, "unverified", None
        status = _check_digest(ds["name"], out, digest, ds.get("sha256"), just_downloaded=False)
        return out, status, digest

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

    try:
        digest = sha256_of(out)
    except OSError as e:
        print(f"[UNVERIFIED] {ds['name']}: downloaded but could not read back "
              f"{out} to verify ({e.strerror or e}) after retries. The "
              f"download is left in place -- rerun later to verify it "
              f"(add --force only if you actually want to redo the download).")
        return out, "unverified", None
    status = _check_digest(ds["name"], out, digest, ds.get("sha256"), just_downloaded=True)
    return out, status, digest


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
        out, status, digest = fetch_one(ds, data_dir, force=False)
        if status == "unverified":
            sys.exit(1)  # fetch_one already printed why
        print(digest)
        if status == "mismatch":
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
    mismatched = [ds["name"] for ds, (_, status, _d) in zip(todo, results) if status == "mismatch"]
    unverified = [ds["name"] for ds, (_, status, _d) in zip(todo, results) if status == "unverified"]
    if mismatched or unverified:
        parts = []
        if mismatched:
            parts.append(f"{len(mismatched)} mismatched: {mismatched}")
        if unverified:
            parts.append(f"{len(unverified)} could not be read to verify "
                          f"(I/O error, rerun to retry): {unverified}")
        print(f"\n[FAILED] " + "; ".join(parts))
        sys.exit(1)
    print(f"\nall {len(results)} dataset(s) verified OK")


if __name__ == "__main__":
    main()
