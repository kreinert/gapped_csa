#!/usr/bin/env python3
"""ILP scaling + bound tightness on growing instances.

For each (fasta, shape, max_add):
  * all five heuristics through ./link_stats (sweep's converged Phase II
    settings), recording |C| and the largest add actually used by a link;
  * ./lower_bound (LB0/LB1/DOM/LAG) at --max-add A and at the largest add any
    heuristic used (the bound that is valid for those runs);
  * the compact ILP (tools/exact_ilp.py, CP-SAT) with a time limit;
  * the original ./ilp_baseline (O(L^2) universe, CBC) with a time limit,
    skipped above --baseline-max-n.
Writes one CSV row per instance.
"""
import argparse, csv, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ENV = dict(os.environ, GCSA_RETARGET_ITERS="1000000", GCSA_LNS_ITERS="1000000",
           GCSA_DISABLE_RETARGET_TIMEOUT="1", GCSA_PFDP_REENUM="1", GCSA_QUIET="1")
ALGOS = ["gs", "gd", "dep", "pfdp", "pfdpi"]


def kv(s):
    return dict(re.findall(r"(\w+)=(\S+)", s))


def run(cmd, timeout=None, env=None):
    t0 = time.time()
    # own process group, so a timeout also kills cbc spawned by ilp_baseline
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
                         env=env, start_new_session=True)
    try:
        out, _ = p.communicate(timeout=timeout)
        return out, p.returncode, time.time() - t0
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, 9)
        p.communicate()
        return "", "timeout", time.time() - t0


def load_text(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--probe", default=os.path.join(ROOT, "link_stats"))
    ap.add_argument("--ilp-time", type=float, default=300)
    ap.add_argument("--baseline-time", type=float, default=300)
    ap.add_argument("--baseline-max-n", type=int, default=120000)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--no-highs", action="store_true")
    ap.add_argument("instances", nargs="+", help="fasta:shape:max_add")
    a = ap.parse_args()
    new = not os.path.exists(a.out)
    f = open(a.out, "a", newline="")
    w = None
    for spec in a.instances:
        fasta, shape, A = spec.rsplit(":", 2)
        A = int(A)
        row = dict(fasta=os.path.basename(fasta), shape=shape, max_add=A)
        # heuristics
        used = A
        for al in ALGOS:
            out, rc, dt = run([a.probe, fasta, shape, str(A), al], timeout=3600, env=ENV)
            d = kv(out)
            row[f"C_{al}"] = d.get("C", "")
            row[f"addmax_{al}"] = d.get("max_add_used", "")
            row[f"t_{al}"] = round(dt, 2)
            row["m"] = d.get("m", row.get("m", ""))
            if d.get("max_add_used"):
                used = max(used, int(d["max_add_used"]))
        Cs = [int(row[f"C_{al}"]) for al in ALGOS if row[f"C_{al}"] != ""]
        row["C_best"] = min(Cs) if Cs else ""
        row["add_used_max"] = used
        # bounds at A (and dump for the ILP)
        dump = f"/tmp/scal_{os.getpid()}.bin"
        out, rc, dt = run([os.path.join(ROOT, "lower_bound"), "-g", fasta, "-s", shape,
                           "--max-add", str(A), "--dump-runs", dump,
                           "--ub", str(row["C_best"] or 0)])
        d = kv(out)
        for k in ("LB0", "LB1", "DOM", "LAG", "LB", "K"):
            row[k] = d.get(k, "")
        row["t_lb"] = round(dt, 2)
        row["dump_bytes"] = os.path.getsize(dump) if os.path.exists(dump) else ""
        if used > A:
            out, rc, dt = run([os.path.join(ROOT, "lower_bound"), "-g", fasta, "-s", shape,
                               "--max-add", str(used)])
            row["LB_at_used"] = kv(out).get("LB", "")
        else:
            row["LB_at_used"] = row["LB"]
        # compact ILP
        out, rc, dt = run([sys.executable, os.path.join(HERE, "exact_ilp.py"), dump,
                           "--time", str(a.ilp_time), "--threads", str(a.threads)],
                          timeout=a.ilp_time * 3 + 60)
        d = kv(out)
        row["ilp_status"] = d.get("status", rc)
        row["ilp_best"] = d.get("best", "")
        row["ilp_bound"] = d.get("bound", "")
        row["ilp_nvars"] = d.get("nvars", "")
        row["ilp_secs"] = round(float(d["secs"]), 2) if "secs" in d else round(dt, 2)
        # original ILP
        n = len(load_text(fasta))
        row.update(base_opt="", base_status="", base_cands="", base_secs="",
                   hb_status="", hb_best="", hb_bound="", hb_secs="")
        if n <= a.baseline_max_n:
            lp = f"/tmp/scal_{os.getpid()}.lp"
            out, rc, dt = run([os.path.join(ROOT, "ilp_baseline"), load_text(fasta), shape,
                               "--max-add", str(A), "--lp-out", lp], timeout=a.baseline_time, env=ENV)
            mo = re.search(r"optimal \|C\| = (\d+)", out)
            mc = re.search(r"candidates=(\d+)", out)
            row["base_opt"] = mo.group(1) if mo else ""
            row["base_status"] = "ok" if mo else rc
            row["base_cands"] = mc.group(1) if mc else ""
            row["base_secs"] = round(dt, 2)
            row.update(hb_status="", hb_best="", hb_bound="", hb_secs="")
            if os.path.exists(lp) and not a.no_highs:
                out, rc, dt = run([sys.executable, os.path.join(HERE, "solve_lp_highs.py"), lp,
                                   str(a.baseline_time)], timeout=a.baseline_time * 2 + 60)
                d = kv(out)
                row.update(hb_status=d.get("status", rc), hb_best=d.get("best", ""),
                           hb_bound=d.get("bound", ""), hb_secs=d.get("secs", round(dt, 2)))
            for p in (lp, lp + ".cbc.sol"):
                if os.path.exists(p):
                    os.remove(p)
        else:
            row["base_status"] = "skipped"
        if os.path.exists(dump):
            os.remove(dump)
        if w is None:
            w = csv.DictWriter(f, fieldnames=list(row.keys()))
            if new:
                w.writeheader()
        w.writerow(row)
        f.flush()
        print(" ".join(f"{k}={v}" for k, v in row.items()), flush=True)


if __name__ == "__main__":
    main()
