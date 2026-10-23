#!/usr/bin/env python3
"""Split the heuristic-vs-optimum gap into a *universe* part and a *search* part.

For each fasta:shape:max_add, report
  opt_full    ilp_baseline --universe full    (the problem's optimum, adds <= A)
  opt_legacy  ilp_baseline --universe legacy  (optimum over the heuristics' own
              candidate view, enumerate_candidates_: maximal lcp>=add+1 runs)
  C_best      best of the five heuristics (defaults, converged Phase II)
  C_xlcp      best of the five with GCSA_CROSS_LCP=1
search gap = C_best - opt_legacy ; universe gap = opt_legacy - opt_full.

Usage: universe_gap.py ./link_stats TIME_LIMIT_S fasta:shape:max_add ... > out.csv
"""
import os, re, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
ENV = dict(os.environ, GCSA_RETARGET_ITERS="1000000", GCSA_LNS_ITERS="1000000",
           GCSA_DISABLE_RETARGET_TIMEOUT="1", GCSA_PFDP_REENUM="1", GCSA_QUIET="1")
probe, tlimit = sys.argv[1], float(sys.argv[2])
def run(cmd, env, timeout):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, env=env, start_new_session=True)
    try: return p.communicate(timeout=timeout)[0]
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, 9); p.communicate(); return ""
print("fasta,shape,max_add,m,opt_full,opt_legacy,C_best,C_xlcp,xlcp_add_used")
for spec in sys.argv[3:]:
    fasta, shape, A = spec.rsplit(":", 2)
    text = "".join(l.strip() for l in open(fasta) if not l.startswith(">")).upper()
    res = {}
    for u in ("full", "legacy"):
        out = run([os.path.join(ROOT, "ilp_baseline"), text, shape, "--max-add", A, "--universe", u,
                   "--lp-out", f"/tmp/ug_{os.getpid()}.lp"], ENV, tlimit)
        mo = re.search(r"optimal \|C\| = (\d+)", out); res[u] = mo.group(1) if mo else ""
        mm = re.search(r"m=(\d+)", out); res["m"] = mm.group(1) if mm else ""
    def best(env):
        cs, adds = [], []
        for al in ("gs", "gd", "dep", "pfdp", "pfdpi"):
            k = dict(re.findall(r"(\w+)=(\S+)", run([probe, fasta, shape, A, al], env, 3600)))
            if "C" in k: cs.append(int(k["C"])); adds.append(int(k["max_add_used"]))
        return min(cs), max(adds)
    cb, _ = best(ENV); cx, ax = best(dict(ENV, GCSA_CROSS_LCP="1"))
    print(f"{os.path.basename(fasta)},{shape},{A},{res['m']},{res['full']},{res['legacy']},{cb},{cx},{ax}", flush=True)
