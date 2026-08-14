# bench/

This directory holds only code -- no sequence data is ever committed here. See
`datasets.py`'s docstring for how each dataset's `kind` (synthetic /
fetched / provided) decides where its FASTA comes from. 

## Setup

```bash
cd gapped_csa
make gcsa simulate_repeats simulate_random simulate_pangenome
```

## Configuring paths

Where the benchmark data actually lives is machine-specific and never
tracked by git. `config.py` resolves `DATA_DIR` / `BIN_DIR` / `TMP_DIR` in
this order (first match wins): 
1. CLI flag (`--data-dir` etc., one-off), 
2. environment variable (`GCSA_BENCH_DATA_DIR` etc.)
3. gitignored config file `config.local.py`: The file `config.local.py.example` serves as a template. You can create a copy named `config.local.py` and modify the fields accordingly. 
4. generic built-in default
that assumes `bench/` sits next to `src/` and `data/` as siblings
under the checkout. 

Run `python3 config.py` to see what each setting currently
resolves to and which of those four sources it came from.

## Running

```bash
cd bench
./run_suite.py --dry-run                      # see the full plan first
./fetch_data.py                               # download real-genome inputs
                                              #   (see datasets.py -- URLs need
                                              #   picking/pinning before this
                                              #   does anything real)
./run_suite.py                                # everything -> results/suite.csv
```

Useful flags:

- `--only-category random repetitive` / `--only-dataset name1 name2`: scope
  a run down while iterating.
- `--shapes ... --algos ... --max-adds ...`: override the sweep in
  `datasets.py` for a quick pass (the full default matrix is 4 algorithms x
  2 max_adds x 4 shapes x every dataset. Real-genome-scale inputs are slow
  enough that you'll usually want a narrower `--shapes`/`--algos` first).
- `--resume`: if partial results exist at `--out`, skip `(dataset, shape, algo, max_add)` rows already present so an interrupted run (or one you're extending with a new dataset) doesn't redo everything.
- `--keep-synthetic` -- don't delete generated FASTAs from `--tmp-dir`
  afterward, if you want to inspect one.

Frequent k-mer datasets (and all other datasets with `kind = provided`) are expected to already exist in `DATA_DIR`: they are not generated or downloaded.

## Output

One row per `(dataset, shape, algorithm, max_add)` in (gitignored) `results/suite.csv`. The columns are: 
- the dataset's `gzip_ratio` (serves as an input-compressibility baseline)
- `keep_pct` / `size_pct` (index-compressibility, from `./gcsa`'s own report) 
- `wall_ms`
- `self_test`
- `roundtrip`


## Extending

Add a new dataset by adding one `dict(...)` to `DATASETS` in `datasets.py`.
- **Real-genome** sources need a `url` (and, after the first successful fetch,
its printed `sha256` copied back in -- `./fetch_data.py --print-hash NAME`).
- **Synthetic families** need a generator binary built the same way
`simulate_random`/`simulate_repeats`/`simulate_pangenome` are (a small
`.cpp` in `src/`, a Makefile rule, `-o <path>` as its output flag). `run_suite.py` just shells out to the provided `generator` name.
