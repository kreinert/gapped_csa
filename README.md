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

# Compare / choose compression heuristic (g / gd / d / p / pi)
./compare_algos
./gcsa --algo greedy-size
./gcsa --algo greedy-degree        # richer graph (cand_cache), degree-based greedy MWIS
./gcsa --algo dep-order
./gcsa --algo pseudoforest-dp             # exact DP on the preference pseudoforest (single-fire) + Phase II
./gcsa --algo pseudoforest-dp-iterate     # same DP, iterated to a fixed point, + Phase II
# Phase II (dirty-set unpin/retarget, then the exact cluster LNS) runs after
# every algorithm's Phase I -- see run_phase2_and_lns_ in compress.hpp. Each
# pass has its own switch, so all four combinations are reachable and each
# pass's contribution can be isolated:
GCSA_DISABLE_RETARGET=1 ./gcsa --algo dep-order -g /tmp/ex.fa -s "#.#"   # LNS alone
GCSA_DISABLE_LNS=1      ./gcsa --algo dep-order -g /tmp/ex.fa -s "#.#"   # retarget alone
GCSA_DISABLE_RETARGET=1 GCSA_DISABLE_LNS=1 ./gcsa --algo dep-order -g /tmp/ex.fa -s "#.#"
                                                                          # Phase I output as-is
# Both are value-checked: =0, =false, =no and =off all mean "run this pass",
# so they can be forced back on from an inherited environment. GCSA_LNS_ONLY
# and GCSA_DISABLE_PHASE2 are deprecated spellings of GCSA_DISABLE_RETARGET
# (neither ever disabled the LNS); they still work, print a note, and keep
# their old presence-only semantics, and GCSA_DISABLE_RETARGET overrides both.
GCSA_TRACE_PFDP=1 ./gcsa -g /tmp/ex.fa -s "#.#" --algo pseudoforest-dp
GCSA_TIMING=1    ./gcsa -g genome.fasta -s "#####" --algo pseudoforest-dp   # stage timings
GCSA_THREADS=8   ./gcsa -g genome.fasta -s "#####" --algo pseudoforest-dp   # parallel pref/DP (default=hw)
# Phase II runs a fixed 100 dirty generations (or fewer, if it reaches a fixed
# point first). Precedence: --phase2-iters > GCSA_PHASE2_MAX_ITERS > default.
# A budget below the fixed point may increase |C|:
./gcsa -g genome.fasta -s "#####" --algo pseudoforest-dp --phase2-iters 25
GCSA_PHASE2_MAX_ITERS=2 ./gcsa -g genome.fasta -s "#####" --algo pseudoforest-dp
# Optional adaptive early-stop, disabled by default (STALL=0). Stops when
# consecutive dirty generations each reduce |C| by < MIN_GAIN kept positions;
# the iteration budget still applies:
# GCSA_PHASE2_MIN_GAIN=2 GCSA_PHASE2_STALL=2 ./gcsa -g genome.fasta -s "#####" --algo pseudoforest-dp

# Exact |C| ILP baseline on small texts (needs cbc or glpsol on PATH)
./ilp_baseline "ACGTCTTAAACCCTCGTCTTAAACCCAACGTCTTAAACCC" "#.#"
./ilp_baseline "GCCTTTAAAGGCCTTTAAAGGCCTTTAAAG" "#.#" --max-add 8
# Smaller, faster, NOT a valid lower bound — for inspecting the heuristics' view:
./ilp_baseline "GCCTTTAAAGGCCTTTAAAGGCCTTTAAAG" "#.#" --universe legacy
# Drop intra-interval links from the heuristics' candidate sets (pre-unification):
GCSA_INTRA_LINKS=0 ./gcsa --algo greedy-size
# Let the heuristics' own candidate search merge source runs across an LCP dip
# (not just numeric adjacency, which is what correctness actually needs; off
# by default -- see "Link universe" below):
GCSA_CROSS_LCP=1 ./gcsa -g genome.fasta -s "#####" --algo pseudoforest-dp
# --algo pseudoforest-dp-iterate above is the first-class way to iterate the
# DP; the old env-var toggle still works as a legacy override on plain
# pseudoforest-dp, for scripts that set it:
GCSA_PFDP_ITERATE=1 ./gcsa -g genome.fasta -s "#####" --algo pseudoforest-dp
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

- **Source selection**: five heuristics — `greedy-size` (availability-aware,
  largest-coverage-first; formerly named `greedy`), `greedy-degree` (degree-aware
  greedy MWIS on a richer graph, below), `dep-order`, `pseudoforest-dp`
  (exact DP on the full preference pseudoforest, below, single-fire), and
  `pseudoforest-dp-iterate` (the same DP, iterated to a fixed point, below).
  Every one of these builds its own initial accepted set (its "Phase I"),
  then hands it to the same shared Phase II: the dirty-set unpin/retarget
  pass followed by the exact cluster LNS (see `run_phase2_and_lns_` in
  `compress.hpp`). Both run by default and each can be switched off on its
  own -- `GCSA_DISABLE_RETARGET=1` runs the LNS alone, `GCSA_DISABLE_LNS=1`
  runs the retarget pass alone, and setting both exposes any algorithm's
  Phase I / DP output as-is. Disabling a pass skips its setup too, so
  `GCSA_DISABLE_LNS=1` also avoids building the LNS dependency graph, which
  is the bulk of that pass's cost on large inputs. (`GCSA_LNS_ONLY` and
  `GCSA_DISABLE_PHASE2` are deprecated spellings of `GCSA_DISABLE_RETARGET`
  -- despite its name the latter never disabled the LNS.) Use
  `./gcsa --algo <name>` or `./compare_algos`.
  Candidate/preference enumeration in `greedy-degree`, `dep-order`, and
  `pseudoforest-dp`/`pseudoforest-dp-iterate` is parallelized via
  `std::thread` (`GCSA_THREADS=N`, default=`hardware_concurrency`), as is
  Phase II's dirty-set re-enumeration for every algorithm, and as is
  `pseudoforest-dp`'s leftover greedy sweep -- that sweep re-enumerates each
  still-unresolved word's best *available* candidate, which dominates its
  cost (~43-54us/word, and ~30% of a `pseudoforest-dp` run's wall clock on
  the pangenome / frequent-k-mer inputs) and is independent per word, so it
  is enumerated in parallel and then accepted sequentially in the same
  size-descending order as before. This is exact rather than approximate:
  within the sweep `accept_` only ever raises `removed_`/`pin_count_`, so
  availability shrinks monotonically and a precomputed candidate that is
  still available is still that word's best; one that has gone is simply
  re-enumerated, exactly as the old sequential loop did for every word.
  `GCSA_LEFTOVER_PARALLEL=0` restores the original one-word-at-a-time sweep
  (for A/B timing), and `GCSA_LEFTOVER_CHUNK=N` (default 1048576) sets how
  many words are enumerated per batch, bounding the extra memory. With
  `GCSA_TIMING=1` the sweep reports `words=`, `accepted=` and `rechecked=`
  so the fast-path hit rate is visible. Set
  `GCSA_TIMING=1` for per-phase ms (Phase I / accept / Phase II; Phase II
  also prints the generation it reached out of the budget and why it
  stopped, as `gen X/Y (<budget source>, stop=<reason>)` with `reason` one
  of `fixed-point|max-iters|adaptive-stall`).
  Phase II runs a fixed budget of 100 dirty generations
  (`kPhase2DefaultIters`), stopping earlier if it reaches a fixed point. The
  budget comes from `--phase2-iters N`, else `GCSA_PHASE2_MAX_ITERS`, else the
  default; `GCSA_TIMING=1` reports the effective value and its source. The
  adaptive early-stop
  (`GCSA_PHASE2_MIN_GAIN`, `GCSA_PHASE2_STALL`) ends generations when kept-drop
  plateaus, but is disabled by default (`GCSA_PHASE2_STALL=0`).
- **Phase II cluster LNS**: since `|C| = m -` total coverage and the only
  constraints are one link per name plus "a chosen link's sources are
  not covered by another chosen link" (both *pairwise* — there is no acyclicity
  requirement), a small set of names can be re-optimized *exactly* with the rest
  of the assignment frozen. Names are dependent when one's candidate sources
  from the other's interval; clusters are bounded BFS components of that graph,
  solved by branch-and-bound over each member's candidates plus "no link".
  Sharing a source is deliberately *not* a dependency — sources are read-only,
  so two links may share them. This recovers the ILP optimum on the two known
  forest-restriction counterexamples (`ACACACACACAC`, `ACGTACGTACGTACGTACGT`).
  It runs after the retarget loop rather than instead of it, for every
  algorithm (see `run_phase2_and_lns_`): retargeting reaches composite-`add`
  links the static candidate cache does not hold. Knobs:
  `GCSA_LNS_CLUSTER` (max names per cluster, default 8), `GCSA_LNS_OPTS`
  (candidates considered per name, 12), `GCSA_LNS_DEGREE` (dependency-graph
  degree cap, 16), `GCSA_LNS_NODES` (enumeration budget per cluster, 20000;
  clusters that exceed it retry at size 4), `GCSA_DISABLE_RETARGET=1` (skip
  the retarget loop and run the LNS alone), `GCSA_DISABLE_LNS=1` (skip this
  pass entirely, graph build included), `GCSA_TRACE_LNS=1`. The payoff grows with
  `--max-add`; on a `--max-add 16` sweep the LNS beat the retarget loop alone
  on the large majority of configs and never lost.
  `GCSA_LNS_CLUSTER` and `GCSA_LNS_NODES` bind *jointly*, and raising them
  together is the single largest quality lever: a bigger `--max-add` densifies
  the dependency graph, so a fixed 8-name ball covers less of each name's real
  neighborhood, but a larger ball also overruns the enumeration budget and falls
  back to 4 names, which costs more than it gained. At cap 8 the node budget is
  irrelevant (0-1 aborts, 20x makes no difference); at cap 16 it is decisive
  (`--max-add 16`: 233826 at 20000 nodes vs 169388 at 400000). `GCSA_LNS_AUTO=1`
  sizes both from the measured mean degree. It is off by default because it
  costs 30-80x runtime — a dial, not a free win.
- **Exact pseudoforest DP (`pseudoforest-dp`)**: instead of a forest, builds
  every still-unresolved name's *preference* graph directly — its single
  best candidate, out-degree ≤ 1 (a name has at most one preferred source),
  with a real edge only when that candidate's source rank range actually
  overlaps the covered set of the source's own preferred candidate. This
  graph is a pseudoforest (every component is a tree, or a tree plus one
  extra edge closing a single cycle — see `doc/exact_pseudoforest_dp.pdf`,
  "Sparsifying the graph"), so no edge is ever dropped to force acyclicity
  (unlike a plain preference-forest build): tree components solve with a KEEP vs
  COMPRESS recurrence, and each unicyclic component's one cycle is solved
  exactly by fixing one cycle node's status and folding the ring into two
  forced scenarios (whichever is cheaper). A single greedy leftover sweep
  then mops up whatever the one-candidate-per-name sparsification couldn't
  see (every candidate for a name, not just its preferred one — including
  self-references), followed by Phase II. `GCSA_PFDP_RANK_AWARE=0` reverts
  to treating every candidate-with-a-preference as a real dependency edge
  (the pre-conflict-checked, more conservative graph). `GCSA_TRACE_PFDP=1`
  dumps each round's graph and cycle folds.
  `--algo pseudoforest-dp-iterate` (or the legacy `GCSA_PFDP_ITERATE=1`
  override on plain `pseudoforest-dp`) repeats the extract/solve step to a
  fixed point — each further round re-extracting the *currently available*
  best candidate for every still-unresolved name over the shrunk residual
  instance — before ever falling back to the leftover sweep (the doc's
  "Iterative DP" idea, `GCSA_PFDP_MAX_ROUNDS` caps the round count if ever
  needed, unbounded by default).
  Those later rounds pick from round 0's cached, availability-agnostic
  top-`kCandCacheDefaultCap` (64) list per word, so once accepts have consumed
  a word's whole cached list it leaves the DP for good even though the SA may
  still hold available lower-coverage candidates for it — which is exactly
  what the leftover sweep afterwards picks up, and why `pfdpi`'s leftover gain
  is small but never zero. `GCSA_PFDP_REENUM=1` (off by default) repairs
  exactly that case. Later rounds still read the cache first; a word whose
  cached list yields nothing available *and* was truncated by the cap then
  re-enumerates from the SA with the availability filter on, recovering the
  candidates the cap hid. Both conditions matter -- a word whose list was
  never capped already has its whole universe cached, so if none of it is
  available then none exists, and re-enumerating it would be pure waste.
  That targeting is what keeps this cheap: on `t_big`/`#####`/`--max-add
  1024` the repair fires for one word and adds +101ms / +134ms on rounds 1
  and 2 (1.08x the candidate-enumeration total), where blanket
  re-enumeration of every unresolved word cost +1494ms / +1450ms for the
  identical result (1.95x). Rounds enumerate over *intervals*, not over DP
  nodes, which is why the blanket form is so expensive: accepting a few
  hundred words barely shrinks the population walked next round, so round 1
  re-enumerates almost as many words as round 0 to discover that only a
  handful still have a usable candidate. The cache is never written with a
  filtered list: Phase II still reads round 0's unfiltered universe and
  filters it for itself. `GCSA_TIMING=1` reports how many words each round's
  repair recovered. Measured on 70 synthetic configurations it moved the
  leftover sweep's work into the exact DP (`n_leftover` dropping to 0
  wherever it had been non-zero) with final `|C|` unchanged in every one --
  i.e. it fixes the blindness, but on those inputs the greedy sweep was
  already reaching the same solution the DP does.
  The point of that is structural: if per-round re-enumeration means the
  leftover sweep can never accept anything, the pipeline collapses from four
  stages (DP / leftover / retarget / LNS) to three, which is materially easier
  to describe. That redundancy is *checked, not assumed*. It does not follow by
  construction -- a candidate whose source name is absent from the graph's
  `by_name` is dropped, and the round loop stops on "the DP chose nothing"
  rather than "no candidate remains" -- so an iterating, re-enumerating run
  that reaches its own fixed point and *then* sees the sweep accept a word
  prints a loud `[INVARIANT-VIOLATION]` line, and `GCSA_TIMING=1` records
  `invariant=held|VIOLATED|n/a` on every run (`n/a` when the premise does not
  apply: single-fire, re-enumeration off, or `GCSA_PFDP_MAX_ROUNDS` cutting
  the loop short). Across the 70 configurations the invariant held on all 70;
  the warning path itself was exercised by disabling the repair in a throwaway
  build, which correctly reported 37 surviving words. Note the synthetic set is
  a weak test -- 62 of those 70 already had `n_leftover == 0` without the
  change, so only 8 genuinely exercised it. The pangenome and frequent-k-mer
  inputs, where the sweep resolves thousands of words, are the real check. Worth re-measuring on the pangenome
  and frequent-k-mer inputs, where the leftover sweep resolves far more words
  (67k on `pangenome_ecoli_real_n10`, 10k on `kestrel_freq_kmers`) and so has
  more chance to be beaten by an exact solve. Kept as its own algorithm rather than a
  hidden env-var toggle so it can be swept and compared on equal footing
  with the rest. This is not a strict win despite every round being an exact
  solve: each round commits, in one batch snapshot, to one candidate per
  name and can't fall back to that name's second-best candidate within the
  round the way the live, one-at-a-time leftover sweep can. Measured on the
  `bench/` E. coli genomes across shapes and `--max-add`, the end-to-end
  effect is a coin flip — usually zero, otherwise a handful of `|C|`
  positions out of several million kept, in either direction — so it's off
  by default pending a wider sweep.
- **Degree-aware greedy MWIS (`greedy-degree`)**: `greedy-size` weighs a
  candidate only by its own coverage, ignoring how many *other* words'
  candidates it knocks out by consuming their shared source rows.
  `greedy-degree` builds a richer conflict graph -- every word's full
  `cand_cache` (up to 64 candidates each, the existing cap that keeps this
  from re-enumerating the full link universe and OOMing), not just each
  word's single preferred pick -- and scores each candidate as
  `coverage - sum(best still-conflicting candidate's coverage, one term per
  distinct conflicting word)`, popped in decreasing score order via a lazy
  priority queue (stale entries are recomputed on pop and re-pushed; a
  candidate's true score only rises as conflicting words resolve, so a
  popped entry matching its pushed score is provably the current max -- no
  separate leftover pass is needed). Same-word candidates are *not*
  materialized as graph edges (that would be `O(64^2)` per word and is
  redundant with `accept_()`'s own exclusivity bookkeeping); a word drops
  out of every remaining candidate's penalty the moment any of its own
  candidates is accepted, via a resolved-words set. The penalty term must be
  each conflicting word's *locally* best conflicting candidate, not that
  word's globally preferred one -- using the global preference would just
  reconstruct `pseudoforest-dp`'s own preference edges and add no
  information beyond what the exact DP already solves. Availability
  (`accept_()`/`cand_available_()`) is authoritative throughout, same as
  every other algorithm. `GCSA_TRACE_GREEDY_DEGREE=1` traces each accept.
  Empirically mixed so far (measured via `./compare_algos`, not yet swept as
  widely as the others): on shorter/lower-repetition synthetic inputs it can
  lose to plain `greedy-size` by a few percent `|C|`, since the richer graph
  sometimes talks itself out of a locally-large candidate over conflicts that
  never actually collide once resolution order is fixed; on more repetitive
  inputs (more reps of a fixed motif) it converges to within a fraction of a
  percent of `dep-order`/`pseudoforest-dp` and can beat
  `greedy-size` outright (e.g. ~3.4% fewer kept positions on a `motif-len 200`
  / `reps 40..100` sweep). Correctness (self-test + serialized round-trip) is
  clean in every configuration tried.
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
  `enumerate_candidates_` (the heuristics' search, as opposed to
  `enumerate_all_links_`/`collect_diff_problem`'s exact universe above) still
  applies that `lcp >= add+1` rule by default when grouping consecutive source
  ranks into one run, even though it's not required — every rank the grouping
  loop considers already has a verified-valid successor (`pred_lexpos_`/
  `succ_lexpos_` are exact inverses), so a low-LCP run between two otherwise
  fine, rank-adjacent sources can fragment one real link into pieces too small
  to individually clear `kMinCoverage`. `GCSA_CROSS_LCP=1` drops that clause
  (numeric adjacency alone decides grouping). Never changes correctness
  (self-test/round-trip pass either way — it only ever adds valid candidates
  to the search) and measured strictly better `|C|` on real genomic inputs
  across every algorithm, but is not a strict win in general: on a repeat-
  heavy synthetic input, Greedy/DepOrder got slightly *worse* (bigger merged
  runs make `emit_run_`'s all-or-nothing availability check more fragile —
  one already-pinned/removed rank anywhere in a larger window blocks the
  whole candidate, where a smaller LCP-bounded window would have left a
  neighboring piece usable) while the DP-based algorithms still improved.
  Off by default pending a sweep across `bench/`'s tunable-repetitiveness and
  pangenome categories, where this trade-off should actually get decided.
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
