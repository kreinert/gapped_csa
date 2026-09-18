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
  afterward, if you want to inspect one. A kept file is reused rather than
  regenerated on the next run into the same scratch dir, so this doubles as a
  cache; that is safe because generated inputs are written to a temp name and
  `os.replace`d into position, making a file at the final path always
  complete.
- `--sa-cache DIR` / `--prepare-run` -- prepared gapped-SA blobs. The SA
  depends only on `(text, shape)`, so the ~15 runs of one cell otherwise
  rebuild it ~15 times; across the 2026-09-09/10 sweeps that was 47.9 h of
  130.2 h total, and on the two large datasets it is 40-49% of every run's wall
  time. `--prepare-run` materializes each dataset once, builds one blob per
  shape, records the dataset's raw/gzip sizes in a JSON sidecar beside it, and
  exits without running any experiment. The sweep then loads the blob and
  **never opens the genome at all** -- which also takes the
  generate/gzip/delete of a `concat` input out of the parallel path entirely,
  removing the failure mode that corrupted the 2026-09-10 run. Budget ~25 bytes
  per SA entry (~174 GB for the full 90-blob grid; 94% of that is
  `human_freq_kmers` and `pangenome_ecoli_real_n100`). A cell with no blob
  falls back to building from the genome, with a `[prepare-miss]` note; a blob
  that is corrupt or built for another shape is refused loudly rather than
  silently rebuilt.

  **Preparation parallelizes too**, and shards on its own axis: over the
  distinct `(dataset, shape)` pairs of the whole run, not over rows. Row
  sharding deliberately scatters one dataset's experiments across every task,
  which is right for the sweep and exactly wrong here -- each task would see
  most shapes of most datasets and rebuild most of the cache (8.0x the total
  work on the default 96-pair grid). Pair sharding builds each blob exactly
  once, weighted by input size so the two huge datasets spread across tasks
  rather than landing together (measured imbalance 1.22x heaviest:lightest over
  8 shards). Already-present entries are skipped and every write is atomic, so
  re-running or overlapping is safe.

  ```bash
  # prepare, in parallel -- same --array as the sweep
  ./run_suite.py --sa-cache /scratch/sa --prepare-run --shard $I/$N
  # then the sweep, also in parallel
  ./run_suite.py --sa-cache /scratch/sa --shard $I/$N ...
  ```

  The two must not overlap in time: a sweep task that starts before its blobs
  exist falls back to building from the genome. Sequence them with a SLURM
  dependency (`--dependency=afterok:<prepare job id>`), which also gives the
  prepare array its own, larger `--mem`.
- `--tmp-tag NAME` -- put this run's generated inputs in `--tmp-dir/NAME`
  instead of the automatic per-shard subdirectory. Needed whenever the shard
  index isn't unique across the processes running at once: several unsharded
  runs together, or two array jobs whose task indices coincide
  (`--tmp-tag $SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID` covers that).
- `--log-dir DIR` -- keep every `./gcsa` invocation's full stdout+stderr as
  `DIR/<dataset>__<shape>__<algo>__ma<max_add>[__noretarget|__nolns|__nophase2].log`
  (command,
  env overrides, exit code, wall time, then the raw output). Off by default
  -- without it you only get the parsed CSV row plus, on failure, a 200-char
  stderr snippet in `status`. Written on both success and failure, so this
  is the first place to look when a run errors or fails `self_test`/
  `roundtrip`; the row's own `log_path` column points straight at its file.
- `--disable-retarget` / `--disable-lns` -- set `GCSA_DISABLE_RETARGET=1` /
  `GCSA_DISABLE_LNS=1` for every invocation, independently, so all four Phase
  II configurations are reachable. The choice lands in the `phase2` CSV column
  (`on` / `no-retarget` / `no-lns` / `off`) and in the log filename, so the
  four can share one `--out`. `--disable-phase2` is a deprecated alias for
  `--disable-retarget`: despite the name it never disabled the LNS, so a
  `phase2=off` row in a CSV written before this change is what is now labelled
  `no-retarget`. `--resume` re-runs such a row rather than treating it as
  done, which is the safe direction.
- `--env KEY=VALUE [KEY=VALUE ...]` -- any other environment variable for
  every `./gcsa` invocation, e.g. `--env GCSA_PHASE2_MAX_ITERS=50
  GCSA_QUIET=1` (see `compress.hpp` for the full `GCSA_PHASE2_*` /
  `GCSA_QUIET` knob list). Applied after, and so overrides, the automatic
  `GCSA_SKIP_SELFTEST`. Not its own CSV column (arbitrary keys don't fit a
  fixed schema) -- pair it with `--log-dir` if you want a per-row record of
  exactly what was set. The current switches are value-checked, so
  `--env GCSA_DISABLE_LNS=0` really does leave the LNS running; the
  deprecated `GCSA_DISABLE_PHASE2` keeps its old *presence* test, so
  `--env GCSA_DISABLE_PHASE2=0` does NOT
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
- `--shard I/N` -- process only shard I of N (0-based I), e.g. `--shard 0/4`
  .. `--shard 3/4` to split one run across four independent cluster array
  tasks. Balances at the **row** level: every `(dataset, shape, algo,
  max_add)` experiment is weighted by its dataset's *estimated input size*
  and packed into the N shards with a greedy longest-processing-time-first
  bin packing (see `estimate_input_bytes`/`balance_rows` in
  `run_suite.py`), not split round-robin by index and not kept as whole
  per-dataset blocks -- so a dataset as large as `pangenome_ecoli_real_n100`
  or `human_freq_kmers` has its experiments spread across, and run
  concurrently on, several shards instead of stranding one shard with every
  one of that dataset's rows while the others finish early and idle (this
  is what a dataset-level split would still do, and what the older
  index-based split did before that). `--dry-run` prints each dataset's
  estimated size, and with `--shard` also prints how many of this shard's
  rows come from each dataset. A `--shard` run logs its shard's row count
  and total estimated row-weight. Every `--shard I/N` invocation recomputes
  the same partition independently and deterministically, so cluster array
  tasks don't need to coordinate. One caveat: size estimates for `fetched`/
  `provided` datasets come from stat-ing the cached file, so run
  `./fetch_data.py` first -- anything not yet on disk weighs 0 and can land
  in any shard. Trade-off of row-level splitting: a dataset whose rows land
  in more than one shard has its resolve/fetch/generate + gzip cost paid
  once *per shard it appears in*, not once total across the whole run --
  usually a good trade for the load balancing it buys, but worth knowing
  for an expensive `concat`/`fetched` dataset split across many shards.
  Each shard also gets its own scratch directory for generated inputs,
  `--tmp-dir/shard<I>of<N>/`, and removes it on the way out. This is not
  optional bookkeeping: the path a generated input is written to derives from
  the dataset name alone, so shards sharing one scratch dir truncate and
  delete each other's inputs mid-read. That shows up as `cannot open <fasta>`
  from `./gcsa`, as `OSError: [Errno 16] Device or resource busy` from
  `gzip_ratio` (on NFS an unlinked-but-open file is silly-renamed, and reads
  of it then fail), and -- worst, because it is silent -- as a run that exits
  0 having compressed a partly-written file. Budget the disk: a `concat`
  dataset split across N shards is materialized N times over. `--tmp-tag`
  overrides the subdirectory name when the shard index isn't unique enough.
  Filters *after* `--only-category`/`--only-dataset`, over whatever's left.
  Each shard needs its own `--out` so parallel jobs don't clobber each
  other's file -- if you don't pass `--out` explicitly while using `--shard`,
  one is picked for you (`results/suite.shard<I>of<N>.csv`). Combine
  `--shard` with `--jobs` for a two-level cluster setup (N array tasks, each
  running up to `--jobs` invocations at once) -- for example, a SLURM array
  job:
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
