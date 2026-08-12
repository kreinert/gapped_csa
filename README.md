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

# Compare / choose compression heuristic (g / d / t / t2)
./compare_algos
./gcsa --algo greedy
./gcsa --algo dep-order
./gcsa --algo tree-dp             # preference-forest DP + Phase II
./gcsa --algo tree-dp2            # same as tree-dp without Phase II
./gcsa --algo tree-dp3            # tree-dp + exact cluster LNS on top of Phase II
./gcsa --algo tree-dp4            # tree-dp + value-based cycle-edge repair (free)
GCSA_TRACE_DP=1  ./gcsa -g /tmp/ex.fa -s "#.#" --algo tree-dp
GCSA_TIMING=1    ./gcsa -g genome.fasta -s "#####" --algo tree-dp   # stage timings
GCSA_THREADS=8   ./gcsa -g genome.fasta -s "#####" --algo tree-dp   # parallel pref/DP (default=hw)
# Phase II runs a fixed 100 dirty generations (or fewer, if it reaches a fixed
# point first). Precedence: --phase2-iters > GCSA_PHASE2_MAX_ITERS > default.
# A budget below the fixed point may increase |C|:
./gcsa -g genome.fasta -s "#####" --algo tree-dp --phase2-iters 25
GCSA_PHASE2_MAX_ITERS=2 ./gcsa -g genome.fasta -s "#####" --algo tree-dp
# Optional adaptive early-stop, disabled by default (STALL=0). Stops when
# consecutive dirty generations each reduce |C| by < MIN_GAIN kept positions;
# the iteration budget still applies:
# GCSA_PHASE2_MIN_GAIN=2 GCSA_PHASE2_STALL=2 ./gcsa -g genome.fasta -s "#####" --algo tree-dp

# Exact |C| ILP baseline on small texts (needs cbc or glpsol on PATH)
./ilp_baseline "ACGTCTTAAACCCTCGTCTTAAACCCAACGTCTTAAACCC" "#.#"
./ilp_baseline "GCCTTTAAAGGCCTTTAAAGGCCTTTAAAG" "#.#" --max-add 8
# Smaller, faster, NOT a valid lower bound — for inspecting the heuristics' view:
./ilp_baseline "GCCTTTAAAGGCCTTTAAAGGCCTTTAAAG" "#.#" --universe legacy
# Drop intra-interval links from the heuristics' candidate sets (pre-unification):
GCSA_INTRA_LINKS=0 ./gcsa --algo greedy
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

- **Source selection**: six heuristics — `greedy`, `dep-order`, `tree-dp`
  (preference-forest DP KEEP vs COMPRESS, then Phase II unpin/retarget),
  `tree-dp2` (same forest DP without Phase II), `tree-dp3` (same forest DP,
  Phase II unpin/retarget followed by the exact cluster LNS below), and
  `tree-dp4` (`tree-dp` plus cycle-edge repair below). Use
  `./gcsa --algo <name>` or `./compare_algos`.
  Tree-dp preference enumeration and per-root forest DP are parallelized via
  `std::thread` (`GCSA_THREADS=N`, default=`hardware_concurrency`). Set
  `GCSA_TIMING=1` for per-phase ms (pref / forest / dp / accept / leftover /
  Phase II; Phase II also prints `stop=fixed-point|max-iters|adaptive-stall`).
  Leftover greedy reuses the static candidate cache (no re-enum).
  Forest cycle checks walk the parent chain; DP uses dense node ids.
  Phase II runs a fixed budget of 100 dirty generations
  (`kPhase2DefaultIters`), stopping earlier if it reaches a fixed point. The
  budget comes from `--phase2-iters N`, else `GCSA_PHASE2_MAX_ITERS`, else the
  default; `GCSA_TIMING=1` reports the effective value and its source. The
  adaptive early-stop
  (`GCSA_PHASE2_MIN_GAIN`, `GCSA_PHASE2_STALL`) ends generations when kept-drop
  plateaus, but is disabled by default (`GCSA_PHASE2_STALL=0`).
- **Phase II cluster LNS (`tree-dp3`)**: since `|C| = m -` total coverage and
  the only constraints are one link per name plus "a chosen link's sources are
  not covered by another chosen link" (both *pairwise* — there is no acyclicity
  requirement), a small set of names can be re-optimized *exactly* with the rest
  of the assignment frozen. Names are dependent when one's candidate sources
  from the other's interval; clusters are bounded BFS components of that graph,
  solved by branch-and-bound over each member's candidates plus "no link".
  Sharing a source is deliberately *not* a dependency — sources are read-only,
  so two links may share them. This recovers the ILP optimum on the two known
  forest-restriction counterexamples (`ACACACACACAC`, `ACGTACGTACGTACGTACGT`).
  It runs after the retarget loop rather than instead of it: retargeting reaches
  composite-`add` links the static candidate cache does not hold. Knobs:
  `GCSA_LNS_CLUSTER` (max names per cluster, default 8), `GCSA_LNS_OPTS`
  (candidates considered per name, 12), `GCSA_LNS_DEGREE` (dependency-graph
  degree cap, 16), `GCSA_LNS_NODES` (enumeration budget per cluster, 20000;
  clusters that exceed it retry at size 4), `GCSA_LNS_ONLY=1` (skip the retarget
  loop and run the LNS alone), `GCSA_TRACE_LNS=1`. The payoff grows with
  `--max-add`: on `./compare_algos --max-add 16` it beats `tree-dp` on 122 of
  139 configs and none worse (213078 vs 252329 total `|C|`).
  `GCSA_LNS_CLUSTER` and `GCSA_LNS_NODES` bind *jointly*, and raising them
  together is the single largest quality lever: a bigger `--max-add` densifies
  the dependency graph, so a fixed 8-name ball covers less of each name's real
  neighborhood, but a larger ball also overruns the enumeration budget and falls
  back to 4 names, which costs more than it gained. At cap 8 the node budget is
  irrelevant (0-1 aborts, 20x makes no difference); at cap 16 it is decisive
  (`--max-add 16`: 233826 at 20000 nodes vs 169388 at 400000). `GCSA_LNS_AUTO=1`
  sizes both from the measured mean degree. It is off by default because it
  costs 30-80x runtime — a dial, not a free win.
- **Phase I cycle-edge repair (`tree-dp4`)**: the preference forest gives each
  name one edge to its preferred source and drops any edge that would close a
  cycle — but *which* edge dies is decided by name order, not by value. On
  `ACACACACACAC` that keeps a coverage-3 edge and discards the coverage-4 one.
  This pass instead drops the cycle's least-valuable edge; adding the closing
  edge creates exactly one cycle, so removing any single edge of it keeps the
  forest acyclic. It runs inside forest construction (`forest=10.3ms` →
  `10.1ms` on a 22.7 s input, i.e. free) and reaches the optimum on
  `ACACACACACAC`. `GCSA_CYCLE_MIN_GAIN` (default 1) is the coverage improvement
  required before rewiring mid-chain.
- **Link universe (what a differential link may look like)**: an `H_offset`
  entry `(pos, add, num)` decodes to `{C[pos+t] + add*span}`, i.e. it reads
  `num` consecutive *kept* entries of `C` and shifts them. A link is therefore
  decodable iff its source ranks are all kept, its covered ranks are all
  dropped, every covered rank is the `add*span` successor of the matching
  source rank and carries the target k-mer, coverage is `>= kMinCoverage`, and
  the name has at most one offset entry. Keeping the source and dropping the
  covered set already forces `src ∩ covered = ∅`, and `C` holds literal
  positions, so decoding never recurses — there is no acyclicity condition.
  In particular **intra-interval** links (source overlapping `I_c` itself) and
  **non-maximal / low-lcp source runs** are all legal; the note's
  `lcp >= add+1` maximal-run rule is a pruning heuristic, not a correctness
  requirement. `GCSA_INTRA_LINKS=0` restores the old "source entirely outside
  `I_c`" rule for the heuristics.
- **ILP baseline (exact |C| on small instances)**: `make ilp_baseline` builds
  `./ilp_baseline`, which enumerates the **full link universe** above
  (`kMinCoverage=3`, `--max-add` default 8), writes a CPLEX `.lp`, and solves
  with `cbc` / `glpsol` on `PATH` (else brute-force when tiny, else LP-only).
  Since that universe is a superset of what any heuristic can build, the
  reported optimum is a valid lower bound for all of them. Reports optimal
  `|C|`, heuristic `|C|`, optimality gap, a keep-set self-check, an end-to-end
  decode check (materializes the ILP solution and compares `positions_of`
  against brute force), and whether the optimum bounds every algorithm:
  ```
  ./ilp_baseline "ACGTCTTAAACCCTCGTCTTAAACCCAACGTCTTAAACCC" "#.#"
  ./ilp_baseline "GCCTTTAAAGGCCTTTAAAGGCCTTTAAAG" "#.#"
  ```
  `--universe` / `GCSA_LINK_UNIVERSE` picks the candidate set: `full`
  (default, exact), `maximal` (maximal source runs only) or `legacy` (the
  heuristics' `enumerate_candidates_` view). `full` costs `O(L^2)` links per
  source run of length `L`, so on repetitive inputs the LP grows quickly —
  e.g. `GCCTTTAAAG`×12 with `#.#` goes from 77 candidates / 0.07 s (`legacy`)
  to 3625 candidates / 8.4 s (`full`). `maximal` and `legacy` are faster but
  are **not** lower bounds.
- **Suffix array construction**: a linear-time integer-alphabet SA-IS
  (`src/sais.hpp`) is the default; libdivsufsort is no longer a build dependency.
  SA-IS and its indices are currently `int32_t`; for texts longer than 2³¹,
  switch the `Idx` alias in `sais.hpp` and `sais_int`'s return type to 64-bit (or
  drop in `divsufsort64`). The rest of the pipeline is unchanged.
- **Succinct structures**: `S_H` and the name index are counted with an
  Elias–Fano cost model but stored as plain arrays; wire in sdsl-lite for real
  EF/rank-select in a production build.
