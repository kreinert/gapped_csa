#!/usr/bin/env python3
"""Re-score an archived gcsa size line against the uncompressed baseline.

Current `gcsa` prints the baseline itself, on the line after the size report --
if a log has that line this script just reads it back. The estimator below is
for logs written *before* the size report grew that line, where re-running the
build just to get one number would cost hours (the human frequent-k-mer run took
32 minutes).

It reconstructs the baseline from the numbers already in the log:

    baseline = m*pb                            (the full suffix array)
             + sum_c [1 + pb + bytes(|I_c|)]   (one rest-only entry per k-mer)
             + EF(K entry starts in H_base)
             + names                           (identical name index; from the log)

Only the per-k-mer interval sizes |I_c| are missing, and they are pinned tightly
without them: every entry costs 1 marker + pb position bytes + 1..4 count bytes,
and sum_c |I_c| = m caps how many k-mers can possibly need a wide count field
(at most m/256 can have |I_c| >= 256, at most m/65536 can have |I_c| >= 65536,
...). That yields a lower and an upper bound on H_base a few parts in 10^4
apart -- far tighter than the difference the comparison is about.

EF COST MODEL. Old logs were written when serialize.hpp's ef_bits() used
n*(2 + floor(log2(u/n))), which assumed the bucket bitvector was exactly 2n bits
and so understated real Elias-Fano (by 14% for S_H on the human run, where
u/n = 6.2 = 2^2.64 is the worst case for that approximation; see the NOTE in
ef_bits()). This script keeps the legacy formula for legacy logs *on purpose*:
the log's own `total` was computed with it, so scoring the baseline the same way
keeps both sides of the ratio consistent. Both sides carried the same
understatement, so the ratio is stable across the change -- 88.45% under the
legacy model vs 88.51% under exact EF on the human run. Absolute S_H/names byte
counts, though, are not comparable across the two eras.

Usage:
    ./baseline_from_log.py run.log [more.log ...]     # parse logs
    ./baseline_from_log.py --numbers m=... K=... pb=4 H=... names=... total=...
"""
import argparse
import pathlib
import re
import sys


def legacy_ef_bytes(n: int, u: int) -> int:
    """serialize.hpp's ef_bits() as it stood before Sep 18 2026. See EF COST MODEL."""
    if n == 0:
        return 0
    lo = 0
    while (1 << (lo + 1)) * n <= u:
        lo += 1
    return (n * (2 + lo) + 7) // 8


def ef_bytes(n: int, u: int) -> int:
    """Exact Elias-Fano cost -- serialize.hpp's current ef_bits(), verbatim."""
    if n == 0:
        return 0
    u = max(u, n)
    q, base = u // n, 0
    while base < 62 and (1 << (base + 1)) <= q:
        base += 1
    bits = min(n * l + n + -(-u // (1 << l)) for l in range(max(0, base - 1), base + 2))
    return (bits + 7) // 8


def project_implicit(r: dict) -> None:
    """What dropping the stored interval positions would cost, on BOTH sides.

    SA intervals are contiguous and in name order, and so are the per-name kept
    runs in C. So on either side the i-th entry's position is derivable rather
    than stored: EF-code the K run boundaries and read
    `pos = select(i)`, `num = select(i+1) - select(i)`. That removes, per side,
    the marker bytes, the rest_pos bytes and the rest_num bytes -- and S_H with
    them, since what is left of H is a fixed-width array needing no entry index.

    This is a PROJECTION of a structure that is not implemented. It is kept out
    of `gcsa` itself deliberately, so the number cannot end up quoted as a
    measurement. Everything it uses except the two EF terms is measured.
    """
    K, m, C = r["distinct_kmers"], r["m"], r["C"]
    off_bytes = r["h_off_pos_bytes"] + r["h_off_add_bytes"] + r["h_off_num_bytes"]

    # index: C + names + EF(run boundaries in C) + the surviving offsets
    #        + EF(which names carry one)
    idx_ef   = ef_bytes(K, C)
    idx_which = ef_bytes(r["offsets"], K)
    idx_new  = r["bytes_C"] + r["bytes_names"] + idx_ef + off_bytes + idx_which
    # baseline: SA + names + EF(interval boundaries in the full SA). Nothing else.
    bas_ef   = ef_bytes(K, m)
    bas_new  = r["baseline_C"] + r["baseline_names"] + bas_ef

    print("  -- projection: implicit positions (EF-coded boundaries), NOT implemented --")
    print(f"     index    : C {r['bytes_C']:,} + names {r['bytes_names']:,}"
          f" + EF-bounds {idx_ef:,} + offsets {off_bytes:,} + EF-which {idx_which:,}"
          f" = {idx_new:,}  (now {r['bytes_total']:,})")
    print(f"     baseline : SA {r['baseline_C']:,} + names {r['baseline_names']:,}"
          f" + EF-bounds {bas_ef:,} = {bas_new:,}  (now {r['baseline_total']:,})")
    print(f"     => {100.0*idx_new/bas_new:.1f}% of the redesigned baseline"
          f"  (today: {100.0*r['bytes_total']/r['baseline_total']:.1f}%)"
          f"; both sides shrink {100*(1-idx_new/r['bytes_total']):.0f}%"
          f"/{100*(1-bas_new/r['baseline_total']):.0f}%")


def h_base_bounds(m: int, K: int, pb: int):
    """Bounds on the rest-only hash table's byte size.

    Each of the K entries costs 1 + pb + bytes_needed(|I_c|), and
    bytes_needed >= 1 always. The count field only widens for large intervals,
    and sum_c |I_c| = m bounds how many of those there can be.
    """
    fixed = K * (1 + pb + 1)
    extra = 0
    for width_threshold in (256, 1 << 16, 1 << 24):
        extra += min(K, m // width_threshold)   # each such k-mer costs 1 more byte
    return fixed, fixed + extra


def report(name, m, K, pb, total, names_bytes):
    lo_H, hi_H = h_base_bounds(m, K, pb)
    out = []
    for label, H in (("lower", lo_H), ("upper", hi_H)):
        sa = m * pb
        sh = legacy_ef_bytes(K, H)
        base = sa + H + sh + names_bytes
        out.append((label, H, sa, sh, base, 100.0 * total / base))
    print(f"\n{name}   [estimated, legacy EF cost model -- matches this log]")
    print(f"  m={m}  distinct k-mers K={K}  pb={pb}")
    print(f"  index total             : {total:,} bytes")
    print(f"  full SA (old comparison): {m * pb:,} bytes  -> {100.0 * total / (m * pb):.1f}%")
    for label, H, sa, sh, base, pct in out:
        print(f"  baseline ({label} bound) : {base:,} bytes"
              f"  [SA {sa:,} + H {H:,} + S_H {sh:,} + names {names_bytes:,}]"
              f"  -> {pct:.2f}%")
    lo_pct, hi_pct = out[1][5], out[0][5]
    print(f"  => index is {lo_pct:.1f}-{hi_pct:.1f}% of the like-for-like baseline "
          f"({100 - hi_pct:.1f}-{100 - lo_pct:.1f}% smaller)")


SIZE_RE = re.compile(
    r"pb\(bytes/pos\)=(\d+).*?H=(\d+)\s+S_H=(\d+)\s+names=(\d+)\s+total=(\d+)")
M_RE = re.compile(r"SA entries \(m\)\s*:\s*(\d+)")
K_RE = re.compile(r"distinct k-mers\s*:\s*(\d+)")
BASE_RE = re.compile(r"^baseline \(full SA \+ rest-only H\):.*$", re.M)
VS_RE = re.compile(r"^vs\. baseline:.*$", re.M)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("logs", nargs="*", help="gcsa log files to re-score")
    ap.add_argument("--numbers", nargs="+", metavar="KEY=VAL",
                    help="score one run given m=, K=, pb=, names=, total= directly")
    ap.add_argument("--project-implicit", action="store_true",
                    help="also project what EF-coding the interval boundaries "
                         "would cost on both sides (needs a log with the "
                         "composition lines; see project_implicit.__doc__)")
    args = ap.parse_args()

    if args.numbers:
        kv = dict(p.split("=", 1) for p in args.numbers)
        report("(command line)", int(kv["m"]), int(kv["K"]), int(kv["pb"]),
               int(kv["total"]), int(kv["names"]))
        return

    if not args.logs:
        ap.error("give a log file, or --numbers m=... K=... pb=... names=... total=...")

    for path in args.logs:
        text = open(path, encoding="utf-8", errors="replace").read()

        # Logs from a build that already prints the baseline need no estimating.
        already = BASE_RE.search(text)
        if already:
            print(f"\n{path}   [measured -- this log already carries the baseline]")
            print("  " + already.group(0))
            vs = VS_RE.search(text)
            if vs:
                print("  " + vs.group(0))
            if args.project_implicit:
                sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
                from run_suite import parse_gcsa_stdout
                r = parse_gcsa_stdout(text)
                need = ("distinct_kmers", "m", "C", "bytes_C", "bytes_names",
                        "bytes_total", "baseline_C", "baseline_names",
                        "baseline_total", "offsets", "h_off_pos_bytes",
                        "h_off_add_bytes", "h_off_num_bytes")
                if all(k in r for k in need):
                    project_implicit(r)
                else:
                    print("  (no composition lines in this log -- "
                          "rerun with a build from Sep 18 2026 or later)")
            continue

        size, mm, kk = SIZE_RE.search(text), M_RE.search(text), K_RE.search(text)
        if not (size and mm and kk):
            print(f"{path}: no complete gcsa size report found", file=sys.stderr)
            continue
        pb, _H, _SH, names_bytes, total = (int(g) for g in size.groups())
        report(path, int(mm.group(1)), int(kk.group(1)), pb, total, names_bytes)


if __name__ == "__main__":
    main()
