# gapped_csa — Compressed suffix array for (gapped) k-mer → sorted positions

An implementation of the data structure from the note *"Compressed SA for gapped
shapes"* (K. Reinert): map every k-mer or **gapped** k-mer (arbitrary shape) to
the sorted list of its positions in a text, by

1. building a **gapped suffix array** for the shape via the DisLex transform
   (Horton 2008), using a linear-time integer-alphabet **SA-IS** for the suffix
   sorting, then
2. **differentially compressing** the suffix array using LCP intervals, storing
   the result as a byte-packed hash table with two entries per shape
   (`H_offset` + `H_rest`) over a compressed position array `C`.

## Layout

```
gapped_csa/
├── src/
│   ├── shape.hpp        # shape parsing + order-preserving gapped-k-mer naming
│   ├── sais.hpp         # linear-time SA-IS for integer alphabets
│   ├── gapped_sa.hpp    # DisLex lextext, gapped SA (SA-IS), Kasai LCP
│   ├── compress.hpp     # LCP-interval differential compression + query
│   ├── serialize.hpp    # byte-level H array + select bitvector (S_H) encoding
│   ├── main.cpp         # CLI: build / validate / query / size report
│   ├── compare_algos.cpp
│   ├── ilp_baseline.cpp # exact |C| ILP over candidate set (cbc/glpsol)
│   └── validate.cpp     # prints the full SA/LCP table (reproduces the note)
├── external/            # (legacy) libdivsufsort — no longer a build dependency
└── Makefile
```

## Build

Header-only; no external libraries are required:

```bash
make
```

## Usage

```bash
./gcsa                                  # runs the note's #.# example + all self-tests
./gcsa -g genome.fasta -s "##.##"       # build for a shape, report compression + verify
./gcsa -g genome.fasta -s "#.#" -q CTTAAC --table
./gcsa -g genome.fasta -s "#####" -r reads.fasta   # locate reads (like kmer/locate.cpp)
./validate "ACGT..." "#.#"              # print the gapped SA + LCP table

# Simulate FASTA: x concatenated copies of a random motif of length y
./simulate_repeats -x 50 -y 200 -o repeats.fasta
./simulate_repeats -x 100 -y 200 --seed 1 -o /tmp/rep.fasta

# Repetition benchmark: concat a 200-bp motif 10..100×; weight-30 shapes with ≥20 consecutive #
./bench_repetition
./bench_repetition --min-rep 10 --max-rep 100 --step 10 --seed 1

# Compare / choose compression heuristic (g / d / f / t / t2)
./compare_algos
./gcsa --algo greedy
./gcsa --algo dep-order
./gcsa --algo greedy-dfs          # DFS from best add=+1 hub, cumulative add
./gcsa --algo tree-dp             # preference-forest DP + Phase II
./gcsa --algo tree-dp2            # same as tree-dp without Phase II
GCSA_TRACE_DFS=1 ./gcsa -g /tmp/ex.fa -s "#.#" --algo greedy-dfs
GCSA_TRACE_DP=1  ./gcsa -g /tmp/ex.fa -s "#.#" --algo tree-dp
GCSA_TIMING=1    ./gcsa -g genome.fasta -s "#####" --algo tree-dp   # stage timings
GCSA_THREADS=8   ./gcsa -g genome.fasta -s "#####" --algo tree-dp   # parallel pref/DP (default=hw)
# Optional Phase II dirty-generation cap (default min(|I|+5, 32); may increase |C|):
GCSA_PHASE2_MAX_ITERS=2 ./gcsa -g genome.fasta -s "#####" --algo tree-dp

# Exact |C| ILP baseline on small texts (needs cbc or glpsol on PATH)
./ilp_baseline "ACGTCTTAAACCCTCGTCTTAAACCCAACGTCTTAAACCC" "#.#"
./ilp_baseline "GCCTTTAAAGGCCTTTAAAGGCCTTTAAAG" "#.#" --max-add 8
```

Shapes use `#` (care) and `.` (don't care), e.g. `#.#`, `##.##`, `#..#..#`.
Contiguous shapes (`#####`) give ordinary k-mers.

## How it works

### Naming (shape.hpp)
Each gapped k-mer of weight `w` is named by its lexicographic rank among all
valid k-mers over `{$,A,C,G,T}` with `$ < A < C < G < T` and `$` only allowed as
a suffix (it comes from padding past the text end). The closed form uses
`g(i)=1+σ·g(i+1)`; this is order-preserving and reproduces the note's arithmetic
codes exactly (`$$=0, AA=2, AC=3, C$=6, TT=20`, name space `21` for `#.#`), so no
separate rank-determination pass is needed.

### Gapped SA (gapped_sa.hpp + sais.hpp)
The DisLex "lextext" is the sequence of k-mer names laid out grouped by
`p mod span`. Its suffix array is the gapped suffix array. The names are a large
integer alphabet (up to ~2⁶²), so we sort them directly with a **linear-time
SA-IS** (induced sorting) for integer alphabets — no byte-encoding blow-up. The
distinct names are first **order-preservingly remapped to compact ranks
`0..k-1`** (sort-unique + `lower_bound`) so the SA-IS alphabet is dense and
small; the original arithmetic names are kept in `G.lex` for hashing/lookup.
SA-IS uses the end-of-string-is-smallest convention (a unique smallest sentinel,
handled internally), which is exactly libdivsufsort's ordering — a randomized
cross-check confirms `sais_int` produces byte-for-byte the same SA on 4000 random
integer strings. LCP is computed in symbol space with Kasai. This reproduces the
note's SA/orig/LCP table exactly.

### Differential compression (compress.hpp)
In the gapped SA, all positions with k-mer `c` form one interval `I_c` (its
"sorted positions"). Prepending `add` symbols to a suffix preserves relative
order and shifts the text position by `+add·span`. So an **LCP-interval** of
depth `≥ add+1` whose symbol at depth `add` is `c` maps (bijectively,
order-preserving) onto the subset of `I_c` that shares that left context. We drop
those positions and store a single pointer:

- `H_offset = (pos, add, num)` — recover `num` positions as `C[pos+t] + add·span`
- `H_rest   = (pos, num)`     — `num` positions stored directly in `C`

`C` is the compressed suffix array (kept positions in rank order). The v1 source
selection is an availability-aware greedy that prefers the largest coverage and,
on ties, the deeper LCP interval (so intervals fall back to deep, stable sources
— exactly the trick in the note). Links are accepted only when coverage is
**≥ 3** (`kMinCoverage` in `compress.hpp`; all heuristics + Phase II share
this floor). On the note's example it keeps **29/41** positions and returns
identical query results.

**Choosing the LCP intervals optimally is the open research core of the project**
(see below); the greedy here is a correct, reasonable baseline.

### Byte-level encoding (serialize.hpp)
The hash table is packed into a byte array `H` with a marker byte per entry:

```
bit7   : H_offset present     bit2   : H_rest present
bit6-5 : offset count bytes-1 bit1-0 : rest count bytes-1
bit4-3 : offset add bytes-1
```

Positions use `pb = ⌈log2(max)/8⌉` bytes; `add`/`num` use 1–4 bytes. A
select-supported bitvector `S_H` marks entry starts; `select(i)` locates entry
`i`, where `i` is the rank of the k-mer name (binary search over the sorted
names, i.e. the note's rank/EF-bitvector lookup). The marker byte reproduces the
note's worked example: `(2000000,16,300),(3000000000,24) → 10100100`, 13 bytes.

## Validation

`./gcsa` runs, on every build:
- brute-force check: `positions_of(c)` equals the true position set for every `c`;
- serialized round-trip: decoding straight from the packed bytes gives the same;
- the marker-byte check against the note's example.

All pass for contiguous and gapped shapes on the bundled genomes. Example
compression on a 180 kb repetitive input: `####.####` → **40% of the full SA**.

## Limitations & research extensions

- **Source selection**: five heuristics — `greedy`, `dep-order`, `greedy-dfs`,
  `tree-dp` (preference-forest DP KEEP vs COMPRESS, then Phase II
  unpin/retarget), and `tree-dp2` (same forest DP without Phase II).
  Use `./gcsa --algo tree-dp2` or `./compare_algos`.
  Tree-dp preference enumeration and per-root forest DP are parallelized via
  `std::thread` (`GCSA_THREADS=N`, default=`hardware_concurrency`). Set
  `GCSA_TIMING=1` for per-phase ms (pref / forest / dp / accept / leftover /
  Phase II). Leftover greedy reuses the static candidate cache (no re-enum).
  Forest cycle checks walk the parent chain; DP uses dense node ids.
- **ILP baseline (exact |C| on small instances)**: `make ilp_baseline` builds
  `./ilp_baseline`, which enumerates the same candidates as `compress.hpp`
  (`kMinCoverage=3`, `--max-add` default 8), writes a CPLEX `.lp`, and solves
  with `cbc` / `glpsol` on `PATH` (else brute-force when tiny, else LP-only).
  Reports optimal `|C|`, heuristic `|C|`, optimality gap, and a keep-set
  self-check:
  ```
  ./ilp_baseline "ACGTCTTAAACCCTCGTCTTAAACCCAACGTCTTAAACCC" "#.#"
  ./ilp_baseline "GCCTTTAAAGGCCTTTAAAGGCCTTTAAAG" "#.#"
  ```
- **Suffix array construction**: a linear-time integer-alphabet SA-IS
  (`src/sais.hpp`) is the default; libdivsufsort is no longer a build dependency.
  SA-IS and its indices are currently `int32_t`; for texts longer than 2³¹,
  switch the `Idx` alias in `sais.hpp` and `sais_int`'s return type to 64-bit (or
  drop in `divsufsort64`). The rest of the pipeline is unchanged.
- **Succinct structures**: `S_H` and the name index are counted with an
  Elias–Fano cost model but stored as plain arrays; wire in sdsl-lite for real
  EF/rank-select in a production build.
