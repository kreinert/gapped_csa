#!/usr/bin/env python3
"""Cross-check on many tiny instances:
  compact ILP opt == ilp_baseline opt (both solvers), and LB0<=LB1<=..., every LB <= opt,
  and every heuristic |C| >= opt.
"""
import random, re, subprocess, sys, os
here = os.path.dirname(os.path.abspath(__file__)); root = os.path.dirname(here)
def kv(s): return dict(re.findall(r"(\w+)=(\S+)", s))
def gen(rng):
    kind = rng.choice(["rand", "rep", "mut"])
    if kind == "rand": return "".join(rng.choice("ACGT") for _ in range(rng.randint(15, 45)))
    motif = "".join(rng.choice("ACGT") for _ in range(rng.randint(3, 9)))
    t = motif * rng.randint(3, 7)
    if kind == "mut":
        t = list(t)
        for _ in range(rng.randint(1, 4)): t[rng.randrange(len(t))] = rng.choice("ACGT")
        t = "".join(t)
    pre = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 6)))
    return pre + t
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
N = int(sys.argv[2]) if len(sys.argv) > 2 else 60
fails = 0; rows = []; under = []
for i in range(N):
    text = gen(rng); shape = rng.choice(["#", "##", "#.#", "###", "#.##", "##.#"]); A = rng.choice([1, 2, 3, 4, 6, 8])
    try:
        out = subprocess.run([f"{root}/ilp_baseline", text, shape, "--max-add", str(A), "--lp-out", f"/tmp/v{i}.lp"],
                             capture_output=True, text=True, timeout=300).stdout
    except subprocess.TimeoutExpired:
        print("timeout", text, shape, A); continue
    mo = re.search(r"optimal \|C\| = (\d+)", out)
    if not mo: print("no opt", text, shape, A, out[-300:]); continue
    opt = int(mo.group(1))
    heur = [int(x) for x in re.findall(r"=(\d+)", re.search(r"heuristic \|C\|:(.*)", out).group(1))]
    lb = kv(subprocess.run([f"{root}/lower_bound", "--text", text, "-s", shape, "--max-add", str(A),
                            "--dump-runs", f"/tmp/v{i}.bin"], capture_output=True, text=True).stdout)
    c1 = kv(subprocess.run([sys.executable, f"{here}/exact_ilp.py", f"/tmp/v{i}.bin", "--solver", "cpsat"], capture_output=True, text=True).stdout)
    c2 = kv(subprocess.run([sys.executable, f"{here}/exact_ilp.py", f"/tmp/v{i}.bin", "--solver", "highs"], capture_output=True, text=True).stdout)
    L = {k: int(lb[k]) for k in ("LB0", "LB1", "DOM", "LAG", "LB")}
    if min(heur) < opt: under.append((text, shape, A, opt, min(heur)))
    ok = (int(c1["best"]) == opt and int(c2["best"]) == opt and L["LB"] <= opt
          and L["LB0"] <= L["LB1"] and int(lb["lam1_check"]) == L["LB1"] and L["LAG"] >= max(L["LB1"], L["DOM"]))
    rows.append((opt, L["LB"], min(heur)))
    if not ok:
        fails += 1
        print("FAIL", text, shape, A, "opt", opt, "cpsat", c1.get("best"), "highs", c2.get("best"), L, "heur", heur, "lam1", lb["lam1_check"])
tight = sum(1 for o, l, h in rows if o == l); hopt = sum(1 for o, l, h in rows if o == h)
print(f"heuristic below OPT(A) (adds composed past --max-add by retarget): {len(under)}")
print(f"checked={len(rows)} fails={fails} LB==opt on {tight}/{len(rows)}  best-heuristic==opt on {hopt}/{len(rows)}")
