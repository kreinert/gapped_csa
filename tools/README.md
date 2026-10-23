# Assessment tools (2026-10-22)

See the project doc `claude/assessment-20261022.md` for the reasoning and results.

| tool | what it does |
|---|---|
| `../lower_bound` (`make lower_bound`) | LB0 / LB1 / DOM / Lagrangian lower bounds on \|C\|; `--max-add 8,16,32,...` reports each A from one SA build; `--sa-cache`/`--sa-load` reuse sweep blobs; `--dump-runs` feeds `exact_ilp.py` |
| `../link_stats` (`make link_stats`) | one gcsa build (from a FASTA or a `.gsa` blob); prints \|C\| and the **largest add in the final link set** (needed to pick a valid A) |
| `sweep_bounds.py` | bounds for every (dataset, shape) cell; with `--sa-cache` it reads run_suite's sidecars (`sa/meta/*.json`) and loads the blobs, no FASTA needed |
| `sweep_link_stats.py` | re-runs sweep cells through `link_stats` (from the blobs) to record each run's largest add; with `GCSA_CROSS_LCP=1` set, tests the wider candidate universe |
| `assess_sweep.py` | joins runs + bounds (+ link stats) into per-run deletion efficiency, picking a valid A per run |
| `cluster_bounds.sh` | one-shot `sweep_bounds.py` call for an interactive node |
| `exact_ilp.py` | compact exact ILP (CP-SAT or HiGHS) over the same universe as `ilp_baseline`; exact but weak LP |
| `lp_bound.sh` | LP-relaxation bound of `ilp_baseline`'s O(L^2) model via HiGHS (best bound up to m ~ 1e4–1e5) |
| `universe_gap.py` | splits heuristic-vs-OPT into universe gap (legacy vs full) and search gap |
| `scaling_study.py` | runs all of the above on growing instances (ILP scaling) |
| `validate_bounds.py` | random tiny instances: compact ILP == ilp_baseline, every LB <= OPT |

Python deps: `pandas numpy highspy ortools`; `ilp_baseline` needs `cbc` (e.g. `pip install cbcbox`).

## Cluster workflow

All three read run_suite.py's SA cache (`--sa-cache $S/sa`), so no FASTA is needed,
and all take `--only DATASET ...`, `--shard I/N` for array jobs, and resume on `--out`.

```
cd gapped_csa && make lower_bound link_stats && cd ..
# bounds (single-threaded; ~65 B per SA entry of memory)
python3 gapped_csa/tools/sweep_bounds.py --sa-cache $S/sa --out bounds.csv \
    --only cod_freq_kmers ... --adds 8,16,32,64,128 --iters 30 [--shard I/N]
# optional: largest add each run ends with (export the sweep's GCSA_* env first;
# add GCSA_CROSS_LCP=1 to test the wider candidate universe)
python3 gapped_csa/tools/sweep_link_stats.py --sa-cache $S/sa --out linkstats.csv \
    --only ... --max-adds 8 16 32 [--algos ...] [--shard I/N]
# grade every run
python3 gapped_csa/tools/assess_sweep.py --runs <runs.csv> --bounds 'bounds*.csv' \
    [--link-stats 'linkstats*.csv']
```

Runs with Phase II off cannot exceed their max-add, so bounds at the sweep's own
max-adds are enough for them; runs with Phase II on need the larger adds (64, 128).
