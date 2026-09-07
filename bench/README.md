# bench/

This directory holds only code -- no sequence data is ever committed here. See
`datasets.py`'s docstring for how each dataset's `kind` (synthetic / fetched
/ provided / concat) decides where its FASTA comes from.

## Setup

```bash
cd gapped_csa
make gcsa simulate_repeats simulate_random simulate_pangenome
```

## Configuring paths

Where the benchmark data actually lives is machine-specific and never
tracked by git. `config.py` resolves `DATA_DIR` / `BIN_DIR` / `TMP_DIR` in
this order (first match wins): a CLI flag (`--data-dir` etc., one-off), an
environment variable (`GCSA_BENCH_DATA_DIR` etc.), a gitignored
`config.local.py` (persistent, per-machine -- `cp config.local.py.example
config.local.py` and fill in what you need), or a generic built-in default
that assumes `bench/` sits next to `src/` with a `data` folder as a sibling
of the checkout. Run `python3 config.py` to see what each setting currently
resolves to and which of those four sources it came from.

## Running

```bash
cd bench
./run_suite.py --dry-run                     # see the full plan first
./fetch_data.py                               # download real-genome inputs
                                               # (URLs are pre-filled in
                                               # datasets.py; pin each sha256
                                               # after the first fetch with
                                               # --print-hash NAME)
./run_suite.py                                # everything -> results/suite.csv
```

Useful flags:

- `--only-category random repetitive` / `--only-dataset name1 name2` -- scope
  a run down while iterating. Note a plain, no-flags run already excludes
  `real_genome`/`real_genome_pool` (the individual E. coli strains that only
  exist to be concatenated into the `pangenome_ecoli_real_n{1,10,100}`
  entries -- see `INGREDIENT_ONLY_CATEGORIES` in `datasets.py`); pass
  `--only-category real_genome_pool` or `--only-dataset ecoli_042` etc. if
  you deliberately want to sweep one of those ~100 individual genomes.
- `--shapes ... --algos ... --max-adds ...` -- override the sweep in
  `datasets.py` for a quick pass (see `SHAPES`/`ALGOS`/`MAX_ADDS` there for
  the current full default matrix -- real-genome-scale inputs are slow
  enough that you'll usually want a narrower `--shapes`/`--algos` first).
- `--resume` -- skip `(dataset, shape, algo, max_add, phase2)` rows already
  present in `--out`, so an interrupted run (or one you're extending with a
  new dataset) doesn't redo everything. Appending to a `--out` written by an
  older run_suite.py (different columns) fails loudly rather than silently
  misaligning the CSV -- point `--out` at a new file in that case.
- `--keep-synthetic` -- don't delete generated FASTAs from `--tmp-dir`
  afterward, if you want to inspect one.
- `--log-dir DIR` -- keep every `./gcsa` invocation's full stdout+stderr as
  `DIR/<dataset>__<shape>__<algo>__ma<max_add>[__nophase2].log` (command,
  env overrides, exit code, wall time, then the raw output). Off by default
  -- without it you only get the parsed CSV row plus, on failure, a 200-char
  stderr snippet in `status`. Written on both success and failure, so this
  is the first place to look when a run errors or fails `self_test`/
  `roundtrip`; the row's own `log_path` column points straight at its file.
- `--disable-phase2` -- set `GCSA_DISABLE_PHASE2=1` for every invocation
  (Phase I / leftover-DP output only, no Phase II local-search pass).
  Recorded in the `phase2` CSV column (`on`/`off`), so a phase2-on and a
  phase2-off run of the same `(dataset, shape, algo, max_add)` can sit in
  the same `--out` without one being mistaken for a rerun of the other.
- `--env KEY=VALUE [KEY=VALUE ...]` -- any other environment variable for
  every `./gcsa` invocation, e.g. `--env GCSA_PHASE2_MAX_ITERS=50
  GCSA_QUIET=1` (see `compress.hpp` for the full `GCSA_PHASE2_*` /
  `GCSA_QUIET` knob list). Applied after, and so overrides, the automatic
  `GCSA_SKIP_SELFTEST`. Not its own CSV column (arbitrary keys don't fit a
  fixed schema) -- pair it with `--log-dir` if you want a per-row record of
  exactly what was set. One gotcha: `./gcsa` checks `GCSA_DISABLE_PHASE2`
  for *presence*, not value, so `--env GCSA_DISABLE_PHASE2=0` does NOT
  re-enable Phase II if `--disable-phase2` was also passed -- for that one
  variable, just omit `--disable-phase2` instead of trying to override it.
- `--jobs N` -- run up to N `./gcsa` invocations concurrently (within the
  shape x algo x max_add matrix for one dataset at a time; resolving/
  cleaning up each dataset's FASTA itself stays sequential). Match N to the
  CPUs actually available to this process -- on a shared or scheduled
  machine that's usually less than the box's total core count, e.g.
  `$SLURM_CPUS_PER_TASK` under SLURM, not `nproc`. Default 1 (sequential,
  same behavior as before this flag existed). Output is identical to a
  sequential run either way (same rows, same `--log-dir` files, same
  `--resume` behavior) -- only the wall-clock time and row order change.
- `--shard I/N` -- process only the datasets at index i where `i % N == I`
  (0-based I), e.g. `--shard 0/4` .. `--shard 3/4` to split one run across
  four independent cluster array tasks. Splits by dataset, not by row, so
  every `(shape, algo, max_add)` combination for a given dataset lands in
  the same shard. Filters *after* `--only-category`/`--only-dataset`, over
  whatever's left. Each shard needs its own `--out` so parallel jobs don't
  clobber each other's file -- if you don't pass `--out` explicitly while
  using `--shard`, one is picked for you (`results/suite.shard<I>of<N>.csv`).
  Combine `--shard` with `--jobs` for a two-level cluster setup (N array
  tasks, each running up to `--jobs` invocations at once) -- for example, a
  SLURM array job:
  ```bash
  ./run_suite.py --shard $SLURM_ARRAY_TASK_ID/$SLURM_ARRAY_TASK_COUNT \
                  --jobs $SLURM_CPUS_PER_TASK
  ```
  Afterward, merge the per-shard CSVs (they all share the same header):
  ```bash
  head -1 results/suite.shard0of4.csv > results/suite.csv
  tail -n +2 -q results/suite.shard*of4.csv >> results/suite.csv
  ```

Frequent k-mer datasets (and all other datasets with `kind = provided`) are
expected to already exist in `--data-dir`.

## Output

One row per `(dataset, shape, algorithm, max_add, phase2)` in
`results/suite.csv` (gitignored -- results are a run artifact, not source).
Columns: the dataset's `gzip_ratio` (input-compressibility baseline)
alongside `keep_pct`/`size_pct` (index-compressibility, from `./gcsa`'s own
report), `wall_ms`, `self_test`, `roundtrip`, whether Phase II ran
(`phase2`), and `log_path` (empty unless `--log-dir` was given). Load it
with pandas/whatever and plot `gzip_ratio` vs `size_pct` per the design
doc's main chart.

## Extending

Add a new dataset by adding one `dict(...)` to `DATASETS` in `datasets.py`.
New real-genome sources need a `url` (and, after the first successful fetch,
its printed `sha256` copied back in -- `./fetch_data.py --print-hash NAME`).
New synthetic families need a generator binary built the same way
`simulate_random`/`simulate_repeats`/`simulate_pangenome` are (a small
`.cpp` in `src/`, a Makefile rule, `-o <path>` as its output flag) --
`run_suite.py` just shells out to whatever `generator` name you give it. To
glue several already-resolved datasets into one FASTA (e.g. multiple real
genomes for a real, non-simulated pangenome) use `kind="concat"` with a
`refs` list of `"@name"` references instead -- see `pangenome_ecoli_real_n1`
/ `_n10` / `_n100` (prefixes of `REAL_ECOLI_STRAIN_POOL`, 100 real E. coli
strains deep) for an example.
