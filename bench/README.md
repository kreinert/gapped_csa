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
  a run down while iterating.
- `--shapes ... --algos ... --max-adds ...` -- override the sweep in
  `datasets.py` for a quick pass (the full default matrix is 4 algorithms x
  2 max_adds x 4 shapes x every dataset -- real-genome-scale inputs are slow
  enough that you'll usually want a narrower `--shapes`/`--algos` first).
- `--resume` -- skip `(dataset, shape, algo, max_add)` rows already present
  in `--out`, so an interrupted run (or one you're extending with a new
  dataset) doesn't redo everything.
- `--keep-synthetic` -- don't delete generated FASTAs from `--tmp-dir`
  afterward, if you want to inspect one.

Frequent k-mer datasets (and all other datasets with `kind = provided`) are
expected to already exist in `--data-dir`.

## Output

One row per `(dataset, shape, algorithm, max_add)` in `results/suite.csv`
(gitignored -- results are a run artifact, not source). Columns: the
dataset's `gzip_ratio` (input-compressibility baseline) alongside `keep_pct`
/ `size_pct` (index-compressibility, from `./gcsa`'s own report) and
`wall_ms`, `self_test`, `roundtrip`. Load it with pandas/whatever and plot
`gzip_ratio` vs `size_pct` per the design doc's main chart.

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
`refs` list of `"@name"` references instead -- see `pangenome_ecoli_real_n3`
/ `_n8` / `_n16` (prefixes of `REAL_ECOLI_STRAIN_POOL`, 16 real E. coli
strains deep) for an example.
