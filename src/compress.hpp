// compress.hpp
//
// Differential compression of the gapped suffix array and the k-mer -> sorted
// positions index built on top of it, following the "Compressed SA for gapped
// shapes" note.
//
// ---------------------------------------------------------------------------
// Idea
// ---------------------------------------------------------------------------
// In the gapped SA the suffixes are grouped by their first symbol, which is the
// gapped k-mer name at that position.  So for a k-mer c the interval
// I_c = [lo, hi) of the SA lists exactly the text positions whose gapped k-mer
// is c (its "sorted positions").
//
// Reconstruction relation.  Let x be a lextext position and let its suffix be
// S = lex[x], lex[x+1], ...  The suffix at x-add is (add symbols) + S, and its
// original text position is orig(x) - add*span.  Hence, prepending `add` symbols
// preserves relative order and only shifts the original position by +add*span.
//
// Therefore, if an lcp-interval [i, j) of the SA has string depth >= add+1 and
// its symbol at depth `add` equals c, then shifting every entry of [i, j) by
// +add*span yields, bijectively, the subset of I_c whose left context (the
// `add` preceding windows) matches the interval's shared prefix.  We can then
// drop those positions from storage and replace them by a single pointer
// (H_offset) into the *kept* source interval, recovering them on the fly as
// C[ptr + t] + add*span.  The remaining positions of I_c are stored directly in
// one contiguous run (H_rest).  This is the two-hash-entries-per-shape scheme of
// the note.
//
// ---------------------------------------------------------------------------
// Link universe
// ---------------------------------------------------------------------------
// A hash entry H_offset = (pos, add, num) decodes to
//     { C[pos + t] + add*span : 0 <= t < num },
// where pos is the C-index of a source rank s_lo.  C stores literal positions
// in rank order, so the decode simply reads `num` consecutive *kept* entries
// and shifts them.  A link (name c, add, src = [s_lo, s_hi), covered) is
// therefore decodable iff
//   (L1) every rank in [s_lo, s_hi) is kept (otherwise C[pos+t] is some other
//        rank's position);
//   (L2) for every t, orig(s_lo+t) + add*span is the position of rank
//        covered[t], whose first symbol is c  (succ_lexpos_ valid, name == c);
//   (L3) every rank in `covered` is dropped, and covered ∪ kept(I_c) = I_c;
//   (L4) num = |covered| = s_hi - s_lo >= kMinCoverage;
//   (L5) at most one H_offset per name.
// (L1) and (L3) already force src ∩ covered = ∅.  Nothing else is required: in
// particular decoding never recurses (C holds literal positions), so no
// acyclicity constraint exists and *intra-interval* links — a source run that
// overlaps I_c itself — are perfectly legal, as are non-maximal source runs and
// runs whose internal lcp is below add+1.
//
// enumerate_candidates_ is the heuristics' target-centric view of this
// universe.  By default it emits only maximal runs with lcp >= add+1 (the
// note's lcp-interval argument, a sufficient-but-not-necessary pruning rule --
// see gcsa_cross_lcp()), which keeps the per-interval candidate list short;
// GCSA_CROSS_LCP=1 drops that extra requirement (numeric adjacency alone,
// which correctness actually needs, is enough to keep a run merged) and
// GCSA_INTRA_LINKS=0 restores the historical "source outside I_c" rule.
// enumerate_all_links_ emits the whole universe and is what
// collect_diff_problem hands to the ILP baseline, so the ILP optimum is a true
// lower bound for every algorithm below.
//
// Compression strategies (CompressAlgo):
//   Every algorithm below builds an initial accepted-link set (its "Phase
//   I"), then hands it to the same shared Phase II: the dirty-set greedy
//   unpin/retarget pass (run_phase2_) followed by the exact bounded-cluster
//   LNS pass (run_phase2_local_), run together via run_phase2_and_lns_. The
//   LNS only ever lowers |C| and composes with the retarget loop instead of
//   replacing it (retargeting reaches links -- composite adds -- the static
//   candidate cache does not hold, which is exactly where the LNS alone
//   gives ground on long repeats), so every algorithm gets both by default.
//   GCSA_LNS_ONLY=1 skips the greedy pass and measures the LNS alone;
//   GCSA_DISABLE_PHASE2=1 skips the greedy pass the same way, but the LNS
//   pass still runs regardless -- it isn't gated by that flag. What differs
//   between algorithms below is only Phase I: how the initial accepted set
//   gets built.
//
//   GreedySize   – size-order availability-aware greedy; pins sources forever.
//   GreedyDegree – richer-graph greedy MWIS: instead of sparsifying to one
//               candidate per word and solving that graph exactly, keep every
//               word's whole cached candidate list (up to kCandCacheDefaultCap)
//               as separate nodes, build the real cross-word conflict edges
//               among ALL of them (same source-vs-covered overlap test as
//               PseudoforestDp's ConflictGraph::dep, just against every
//               cached candidate of the source word instead of only its
//               single preferred one -- same-word candidates are mutually
//               exclusive by construction and are never materialized as
//               edges), then repeatedly accept the live candidate with the
//               highest score = coverage - sum, over each distinct other
//               word still contesting it, of that word's best *currently
//               conflicting* candidate (not that word's own top pick --
//               using the top pick would just re-derive edges the exact DP
//               below already prices for free; see the session note this
//               algorithm came out of). Lazy priority-queue greedy: a
//               candidate's true score only ever rises as neighboring words
//               get resolved, so a popped entry whose freshly-recomputed
//               score still matches what it was pushed with is provably the
//               current global best and can be accepted immediately. No
//               separate leftover pass is needed -- every cached candidate of
//               every word already gets its own turn in the same queue.
//   DepOrder  – dependency-order + un-pin / retarget: accept each interval's
//               globally-preferred candidate in reverse topological order of
//               the preference DAG (sinks first), pinning sources dependents
//               still need.
//   PseudoforestDp – "extract preference, solve exactly" DP (see
//               exact_pseudoforest_dp.pdf, "Sparsifying the graph"). Builds
//               the preference pseudoforest (out-degree<=1, no cov/size
//               threshold; a node's candidate only becomes a real
//               dependency edge when its rank range actually conflicts with
//               its source's own candidate -- see ConflictGraph::dep) once, over
//               every unresolved interval; tree components solve with a KEEP
//               vs COMPRESS recurrence, unicyclic components are solved
//               exactly by folding the cycle instead of dropping an edge to
//               force acyclicity. By default ("single-fire") this runs once,
//               then one leftover greedy sweep (sees every candidate for a
//               word, not just its single preferred one -- catches what the
//               DP's sparsification can't) before Phase II.
//   PseudoforestDpIterate – the same DP, but repeats the extract/solve step
//               to a fixed point instead of firing once -- the doc's "IDEA:
//               Iterative DP": each further round re-extracts each still-
//               unresolved node's currently-available best candidate over
//               the shrunk residual instance, before ever falling back to
//               leftover greedy. Not a strict improvement despite each round
//               being an exact solve: each round commits, in one batch
//               snapshot, to one candidate per name and can't fall back to
//               that name's second-best candidate within the round the way
//               the leftover sweep can. Measured on the E. coli benchmark
//               genomes, the end-to-end effect is a coin flip -- a handful
//               of positions either way out of millions kept, never both
//               zero and never large. Kept as its own algorithm (not a
//               hidden env-var toggle on PseudoforestDp) so it can be swept
//               and compared on equal footing with the rest.
//               GCSA_PFDP_MAX_ROUNDS caps the round count if ever needed
//               (unbounded by default).

#pragma once

#include "gapped_sa.hpp"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <string>
#include <chrono>
#include <functional>
#include <utility>
#include <limits>
#include <thread>
#include <atomic>
#include <mutex>
#include <cstring>
#include <cmath>
#include <map>

namespace gcsa {

// Progress / stage / trace logging.  Everything the compressor prints goes to
// stderr through gcsa_log(); set GCSA_QUIET=1 or call gcsa_set_quiet(true) to
// silence it (used by compare_algos so only its table is printed).
inline bool& gcsa_quiet_flag() {
    static bool q = (std::getenv("GCSA_QUIET") != nullptr);
    return q;
}
inline void gcsa_set_quiet(bool q) { gcsa_quiet_flag() = q; }

inline void gcsa_log(const char* msg) {
    if (gcsa_quiet_flag()) return;
    std::fputs(msg, stderr);
}

template <class A, class... Args>
inline void gcsa_log(const char* fmt, A a, Args... args) {
    if (gcsa_quiet_flag()) return;
    std::FILE* out = stderr;
    std::fprintf(out, fmt, a, args...);
}

// Fixed-width byte count needed to hold `v` (1..8). Shared by the serialized
// H_offset/H_rest field widths (serialize.hpp) and by min_coverage_breakeven()
// below, which needs it before any SerializedIndex exists.
inline int bytes_needed(uint64_t v) {
    int b = 1;
    while (b < 8 && v >= (uint64_t(1) << (8 * b))) ++b;
    return b;
}

// Minimum coverage to accept a compression link (all algos + Phase II).
// Preference-forest / DAG *edges* use the same floor (see DepOrder / PseudoforestDp).
// Runtime-adjustable (main.cpp's --min-coverage / GCSA_MIN_COVERAGE); the 3
// below is only the built-in default, read once before CompressedIndex::build
// so every enumerate_candidates_/enumerate_all_links_ call site (which all
// default to `kMinCoverage`) sees the same floor. Not constexpr on purpose.
inline int kMinCoverage = 3;

// Smallest H_offset coverage `num` for which replacing `num` directly-stored
// positions (num*pb bytes in the serialized C array) with a single H_offset
// reference (pb + off_add_b + off_cnt_b bytes added to H) is a net byte
// reduction. The marker byte itself is sunk -- every occurring name pays it
// whether or not it has an offset (see serialize.hpp's serialize_index) --
// so only the off_pos/off_add/off_num fields count against the num*pb
// reclaimed from C. `pb` is the serialized index's bytes/position (see
// estimate_pb below for a pre-DP estimate); `max_add` bounds off_add exactly
// since every accepted add is <= max_add by construction.
inline int min_coverage_breakeven(int pb, int max_add) {
    if (pb <= 0) return kMinCoverage;
    const int off_add_b = std::max(1, bytes_needed((uint64_t)std::max(0, max_add)));
    for (int num = 1; num <= 64; ++num) {
        const int off_cnt_b = std::max(1, bytes_needed((uint64_t)num));
        const int savings = pb * (num - 1) - off_add_b - off_cnt_b;
        if (savings > 0) return num;
    }
    return 64;  // pathological pb; shouldn't happen for realistic inputs
}

// Upper bound on the serialized index's pb (bytes/position), computable from
// the input length alone -- before building the suffix array or running any
// compression. pb = bytes_needed(max(maxpos, maxidx)) in serialize.hpp;
// maxpos is an original text position (< n) and maxidx indexes the kept-
// positions array C, whose size only shrinks once compression runs, so both
// are bounded above by n-1. May overestimate pb slightly once real
// compression shrinks |C| into a smaller byte bracket, which only makes the
// derived min_coverage_breakeven() conservative (a bit higher than the true
// post-compression breakeven), never too aggressive.
inline int estimate_pb(size_t n) {
    return bytes_needed(n > 1 ? (uint64_t)(n - 1) : 0);
}

// Phase II runs this many dirty generations unless it reaches a fixed point
// first. Callers override it per build (CompressedIndex::build's phase2_iters,
// i.e. --phase2-iters), or process-wide with GCSA_PHASE2_MAX_ITERS.
constexpr int kPhase2DefaultIters = 100;

// Above this SA length, Phase II defaults to a single dirty generation (same as
// GCSA_PHASE2_FAST) unless the caller overrides the budget. Keeps Phase II from
// dominating wall time on large DNA / skmer concatenations.
constexpr size_t kPhase2AutoFastM = 1000000;

// Max candidates retained per interval in the shared Phase II / leftover cache.
// Enumeration itself is capped to this (see emit_run_'s admission gate) using
// the same (coverage, add) ordering pick_best_ ranks by, so the true best
// candidate is always retained regardless of the cap -- it can never be
// beaten by anything still admitted. Override with GCSA_CAND_CACHE_CAP.
constexpr int kCandCacheDefaultCap = 64;

// Cluster size the LNS Phase II falls back to when a cluster's exact solve
// blows the enumeration budget (see run_phase2_local_).
constexpr int kLnsFallbackCluster = 4;

// Let the heuristics' enumerate_candidates_ see intra-interval links (source
// run overlapping I_c but never the ranks it covers) — see "Link universe".
// GCSA_INTRA_LINKS=0 restores the old "source entirely outside I_c" rule.
inline bool gcsa_intra_links() {
    static const bool v = [] {
        const char* e = std::getenv("GCSA_INTRA_LINKS");
        return !e || std::atoi(e) != 0;
    }();
    return v;
}

// Off by default (preserves current behavior): enumerate_candidates_'s source-
// run grouping normally also requires G_.lcp[...] >= add+base_depth_ between
// consecutive source ranks, on top of the numeric-adjacency requirement that
// alone is what correctness needs (see "Link universe" above -- L1-L5 never
// reference the source ranks' own shared LCP; pred_lexpos_/succ_lexpos_ are
// exact inverses, so every element the grouping loop even considers already
// has a verified-valid successor before the LCP check runs). The LCP clause
// mirrors the note's original LCP-interval argument (a sufficient condition,
// convenient to compute from the LCP array directly) but can fragment a
// single valid run into pieces too small to individually clear kMinCoverage
// when the source ranks' own prefixes happen to diverge despite each still
// mapping to the right target. GCSA_CROSS_LCP=1 drops that clause, so a run
// only needs numeric adjacency to stay merged -- this widens the *search*,
// never relaxes what's accepted (every extra candidate found is still fully
// L1-L5 valid), so self-test / round-trip correctness is unaffected either
// way. Measured effect (see compress.hpp's git history / the session that
// added this flag): strictly more compression on real genomic data across
// every algorithm tested, but NOT a strict win in general -- on a repeat-
// heavy synthetic input, Greedy and DepOrder got slightly *worse* (bigger
// merged source runs make emit_run_'s all-or-nothing availability check more
// fragile: one already-pinned/removed rank anywhere in a larger window now
// blocks the whole candidate, where a smaller LCP-bounded window would have
// left a neighboring piece usable), while PseudoforestDp came out ahead
// even there, presumably because
// Phase II's fuller unpin/retarget can route around it. Left off by default
// pending a proper sweep across bench/'s tunable-repetitiveness and pangenome
// categories, where this trade-off actually matters.
inline bool gcsa_cross_lcp() {
    static const bool v = [] {
        const char* e = std::getenv("GCSA_CROSS_LCP");
        return e && std::atoi(e) != 0;
    }();
    return v;
}

// Parallelism for independent work (pref enum, per-root forest DP).
// Override with GCSA_THREADS=N (N=1 disables).
inline int gcsa_num_threads() {
    int n = (int)std::thread::hardware_concurrency();
    if (n < 1) n = 1;
    if (const char* e = std::getenv("GCSA_THREADS")) {
        int v = std::atoi(e);
        if (v > 0) n = v;
    }
    return n;
}

// Run fn(i) for i in [0, n). Uses std::thread (portable on Apple clang).
template <class Fn>
inline void gcsa_parallel_for(size_t n, Fn&& fn) {
    if (n == 0) return;
    int nt = gcsa_num_threads();
    if (nt <= 1 || n == 1) {
        for (size_t i = 0; i < n; ++i) fn(i);
        return;
    }
    if ((size_t)nt > n) nt = (int)n;
    std::atomic<size_t> next{0};
    auto worker = [&]() {
        for (;;) {
            size_t i = next.fetch_add(1, std::memory_order_relaxed);
            if (i >= n) break;
            fn(i);
        }
    };
    std::vector<std::thread> pool;
    pool.reserve((size_t)nt - 1);
    for (int t = 1; t < nt; ++t) pool.emplace_back(worker);
    worker();
    for (auto& th : pool) th.join();
}

enum class CompressAlgo {
    GreedySize,            // size-first, pin sources
    GreedyDegree,          // richer-candidate-set greedy MWIS, degree/conflict-aware
    DepOrder,              // dependency-order + un-pin / retarget
    PseudoforestDp,        // exact DP on the full out-degree<=1 preference graph (single-fire)
    PseudoforestDpIterate, // same DP, iterated to a fixed point (see compress_pseudoforest_dp_)
};

inline const char* algo_name(CompressAlgo a) {
    switch (a) {
        case CompressAlgo::GreedySize:            return "greedy-size";
        case CompressAlgo::GreedyDegree:          return "greedy-degree";
        case CompressAlgo::DepOrder:              return "dep-order";
        case CompressAlgo::PseudoforestDp:        return "pseudoforest-dp";
        case CompressAlgo::PseudoforestDpIterate: return "pseudoforest-dp-iterate";
    }
    return "?";
}

// Positive-integer env override, `def` when unset/invalid.
inline int gcsa_env_int(const char* name, int def) {
    if (const char* e = std::getenv(name)) {
        int v = std::atoi(e);
        if (v > 0) return v;
    }
    return def;
}

struct HashEntry {
    // H_offset = (pos, add, num): num positions recovered as C[pos + t] + add*span
    bool     has_offset = false;
    uint64_t off_pos = 0;
    uint32_t off_add = 0;
    uint32_t off_num = 0;
    // H_rest = (pos, num): num positions stored directly at C[pos + t]
    bool     has_rest = false;
    uint64_t rest_pos = 0;
    uint32_t rest_num = 0;
};

class CompressedIndex {
public:
    // phase2_iters: Phase II dirty-generation budget; 0 = unspecified, which
    // falls back to GCSA_PHASE2_MAX_ITERS and then kPhase2DefaultIters.
    void build(const Shape& shape, std::string text,
               int max_add = 8, CompressAlgo algo = CompressAlgo::GreedySize,
               int phase2_iters = 0) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        auto t0 = Clock::now();
        G_ = build_gapped_sa(shape, std::move(text));
        auto t1 = Clock::now();
        regular_mode_ = false;
        span_ = shape.span;
        add_stride_ = span_;
        base_depth_ = 1;
        max_add_ = std::max(1, max_add);
        algo_ = algo;
        phase2_iters_ = std::max(0, phase2_iters);
        compress_();
        auto t2 = Clock::now();
        if (timing) {
            auto ms = [](Clock::time_point a, Clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count();
            };
            gcsa_log(
                "[timing] gapped_sa=%.1fms  compress(%s)=%.1fms  m=%zu intervals~kmers later\n",
                ms(t0, t1), algo_name(algo_), ms(t1, t2), G_.m());
        }
    }

    // Regular (ungapped) suffix array mode: no LexText is built, k-mer
    // occurrences are grouped by LCP >= k (character-level) instead of
    // LCP >= 1 (symbol-level). Functionally equivalent to build(Shape::parse
    // (string(k,'#')), ...) -- same positions_of/locate results -- but via a
    // cheaper construction path (see build_plain_sa in gapped_sa.hpp).
    void build(int k, std::string text,
               int max_add = 8, CompressAlgo algo = CompressAlgo::GreedySize,
               int phase2_iters = 0) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        auto t0 = Clock::now();
        G_ = build_plain_sa(k, std::move(text));
        auto t1 = Clock::now();
        regular_mode_ = true;
        span_ = k;          // window width in characters, used by locate()
        add_stride_ = 1;    // one unit of `add` = one character
        base_depth_ = k;    // full k-mer match = k characters, not 1 symbol
        max_add_ = std::max(1, max_add);
        algo_ = algo;
        phase2_iters_ = std::max(0, phase2_iters);
        compress_();
        auto t2 = Clock::now();
        if (timing) {
            auto ms = [](Clock::time_point a, Clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count();
            };
            gcsa_log(
                "[timing] plain_sa=%.1fms  compress(%s)=%.1fms  m=%zu intervals~kmers later\n",
                ms(t0, t1), algo_name(algo_), ms(t1, t2), G_.m());
        }
    }

    bool regular_mode() const { return regular_mode_; }
    int add_stride() const { return add_stride_; }
    uint64_t name_of_rank(int32_t r) const { return name_of_rank_(r); }

    CompressAlgo algo() const { return algo_; }
    const GappedSA& gsa() const { return G_; }
    const std::vector<int64_t>& compressed_positions() const { return C_; }
    size_t num_kmers() const { return table_.size(); }
    size_t stored_positions() const { return C_.size(); }
    size_t total_positions() const { return G_.m(); }

    // Return all text positions whose gapped k-mer name is `name` (a set).
    std::vector<int64_t> positions_of(uint64_t name) const {
        std::vector<int64_t> out;
        auto it = table_.find(name);
        if (it == table_.end()) return out;
        const HashEntry& e = it->second;
        if (e.has_offset)
            for (uint32_t t = 0; t < e.off_num; ++t)
                out.push_back(C_[e.off_pos + t] + (int64_t)e.off_add * add_stride_);
        if (e.has_rest)
            for (uint32_t t = 0; t < e.rest_num; ++t)
                out.push_back(C_[e.rest_pos + t]);
        return out;
    }

    std::vector<int64_t> locate(const std::string& query) const {
        std::vector<int64_t> hits;
        int last = (int)query.size() - span_;
        for (int i = 0; i <= last; ++i) {
            uint64_t name = name_at(G_.shape, query, (size_t)i);
            for (int64_t p : positions_of(name)) hits.push_back(p);
        }
        std::sort(hits.begin(), hits.end());
        hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
        return hits;
    }

    const std::unordered_map<uint64_t, HashEntry>& table() const { return table_; }

    bool    rank_kept(int32_t rank) const { return !removed_[rank]; }
    int64_t rank_c_index(int32_t rank) const { return rank_to_C_[rank]; }

    struct Interval { uint64_t name; int32_t lo, hi; };
    struct Candidate {
        int32_t target_lo = 0, target_hi = 0;
        uint64_t name = 0;
        int add = 0;
        int32_t src_lo = 0, src_hi = 0;
        std::vector<int32_t> covered;   // ranks of I_c, in source order
        int coverage() const { return (int)covered.size(); }
    };

    // Intervals + the candidate universe (kMinCoverage, no availability
    // filter).  Builds the gapped SA but does not compress / fill C_ or the
    // hash table.  `universe` records which enumeration produced `candidates`.
    struct DiffProblem {
        size_t m = 0;
        std::string universe = "full";
        std::vector<Interval> intervals;
        std::vector<Candidate> candidates;
    };

    // GCSA_LINK_UNIVERSE selects the candidate set:
    //   full    – every decodable link incl. sub-runs (default; exact optimum)
    //   maximal – maximal source runs only (no sub-runs)
    //   legacy  – enumerate_candidates_, i.e. what the heuristics see
    DiffProblem collect_diff_problem(const Shape& shape, std::string text,
                                     int max_add = 8) {
        G_ = build_gapped_sa(shape, std::move(text));
        regular_mode_ = false;
        span_ = shape.span;
        add_stride_ = span_;
        base_depth_ = 1;
        max_add_ = std::max(1, max_add);
        const size_t m = G_.m();
        rank_of_.assign(m, 0);
        for (size_t r = 0; r < m; ++r) rank_of_[G_.sa[r]] = (int32_t)r;
        // enumerate_candidates_ may consult removed_/pin_count_ only when
        // require_avail=true; still initialize so the object is consistent.
        removed_.assign(m, 0);
        removed_words_.assign((m + 63) / 64, 0);
        pin_count_.assign(m, 0);
        pin_owner_.assign(m, 0);
        multi_pin_owners_.clear();
        kept_count_ = m;

        DiffProblem P;
        P.m = m;
        for (size_t r = 0; r < m; ) {
            uint64_t name = name_of_rank_((int32_t)r);
            size_t r2 = r + 1;
            if (regular_mode_) {
                while (r2 < m && G_.lcp[r2] >= base_depth_) ++r2;
            } else {
                while (r2 < m && G_.first_symbol((int32_t)r2) == name) ++r2;
            }
            P.intervals.push_back({name, (int32_t)r, (int32_t)r2});
            r = r2;
        }
        const char* env = std::getenv("GCSA_LINK_UNIVERSE");
        P.universe = env ? env : "full";
        if (P.universe == "legacy") {
            for (const auto& iv : P.intervals) {
                auto cs = enumerate_candidates_(iv.name, iv.lo, iv.hi, /*avail=*/false);
                for (auto& c : cs) P.candidates.push_back(std::move(c));
            }
        } else {
            P.candidates = enumerate_all_links_(P.intervals,
                                                /*sub_runs=*/P.universe != "maximal");
        }
        return P;
    }

    // Materialize C_ / the hash table from an explicit link set, e.g. an ILP
    // solution.  Requires a preceding collect_diff_problem() on this object.
    // Returns false if the links are not simultaneously feasible (two offsets
    // for one name, or a source rank dropped by another link).
    bool apply_links(const std::vector<Interval>& intervals,
                     const std::vector<Candidate>& links) {
        const size_t m = G_.m();
        removed_.assign(m, 0);
        removed_words_.assign((m + 63) / 64, 0);
        pin_count_.assign(m, 0);
        pin_owner_.assign(m, 0);
        multi_pin_owners_.clear();
        kept_count_ = m;
        table_.clear();
        C_.clear();

        std::unordered_map<uint64_t, Candidate> accepted;
        for (const Candidate& link : links) {
            if (accepted.count(link.name)) return false;
            Candidate c = link;
            accept_(c, accepted);
        }
        for (const auto& kv : accepted)
            if (!run_kept_(kv.second.src_lo, kv.second.src_hi)) return false;

        finalize_(intervals, accepted);
        return true;
    }

    // Diagnostic bounds on |C| -- NOT part of the compressed index itself,
    // just a report. Reuses the gapped SA already built by build() (no
    // rebuild, unlike collect_diff_problem()). Two numbers:
    //
    //   floor_positions -- positions whose k-mer occurs < kMinCoverage times.
    //     Provably un-compressible by *any* algorithm, heuristic or exact:
    //     (L3) `covered` is a subset of the word's own interval I_c, and (L4)
    //     needs |covered| >= kMinCoverage, so |I_c| < kMinCoverage means no
    //     candidate for that word can ever exist, full stop. This is the
    //     algorithm-independent floor.
    //
    //   sum_best_coverage -- Sigma, over eligible words, of
    //     best_candidate_(word, avail=false).coverage() -- i.e. the exact
    //     same enumerate_candidates_ call every algorithm below uses for its
    //     own Phase-I preference (pseudoforest-dp/dep-order/greedy-size/greedy-degree
    //     all call this, at this same require_avail=false, before any
    //     conflict is resolved). (L5) caps every word at one H_offset, so no
    //     feasible solution built from THIS candidate search can ever delete
    //     more from a word than its own best candidate's coverage; summing
    //     those per-word caps relaxes only the *cross-word* conflicts (two
    //     words wanting overlapping source/covered ranges), which is exactly
    //     what every algorithm's Phase II then goes and resolves.
    //   i.e.  floor_positions <= m - sum_best_coverage <= any real run's stored_positions().
    //
    //   This is deliberately NOT the same quantity as the (removed) full-
    //   link-universe bound this function used to compute. That version used
    //   enumerate_all_links_ -- the ILP baseline's candidate universe, a
    //   strict superset of what any heuristic here can ever see -- so it
    //   bounded the true theoretical optimum, independent of which algorithm
    //   or candidate search produced it. It cost O(L^2) per maximal same-add
    //   run of length L, which is why --bound was opt-in in the first place. This
    //   version bounds something narrower but arguably more useful day to
    //   day: not "what's theoretically possible for the compression scheme,"
    //   but "what's the ceiling THIS run's actual candidate search could
    //   ever have reached, if Phase II resolved every conflict perfectly."
    //   Because it goes through the real enumerate_candidates_, it
    //   automatically reflects this run's GCSA_CROSS_LCP / GCSA_INTRA_LINKS /
    //   --min-coverage settings, and it's cheap: no O(L^2) sub-run
    //   enumeration, just one already-fast per-interval search (~35ms over
    //   15k eligible words on ecoli_cft073, vs. ~10s for the old universe).
    //   If you want the true, algorithm-agnostic theoretical ceiling instead,
    //   collect_diff_problem()'s "full" universe (what ilp_baseline uses) is
    //   still the reference for that -- it's just no longer what --bound
    //   reports by default.
    struct BoundReport {
        size_t m = 0;
        bool cross_lcp = false;             // this run's GCSA_CROSS_LCP setting
        bool intra_links = true;            // this run's GCSA_INTRA_LINKS setting
        size_t eligible_words = 0;          // words with |I_c| >= kMinCoverage
        size_t eligible_positions = 0;      // Sigma |I_c| over those words
        size_t words_with_candidate = 0;    // ... of which have >=1 valid link
        size_t floor_positions = 0;         // m - eligible_positions
        size_t sum_best_coverage = 0;       // Sigma best_candidate_(...).coverage()
        size_t lower_bound_C = 0;           // m - sum_best_coverage
    };

    BoundReport compute_bound_report() const {
        const size_t m = G_.m();
        std::vector<Interval> intervals;
        intervals.reserve(m / 4 + 1);
        for (size_t r = 0; r < m; ) {
            uint64_t name = name_of_rank_((int32_t)r);
            size_t r2 = r + 1;
            if (regular_mode_) {
                while (r2 < m && G_.lcp[r2] >= base_depth_) ++r2;
            } else {
                while (r2 < m && G_.first_symbol((int32_t)r2) == name) ++r2;
            }
            intervals.push_back({name, (int32_t)r, (int32_t)r2});
            r = r2;
        }

        BoundReport R;
        R.m = m;
        R.cross_lcp = gcsa_cross_lcp();
        R.intra_links = gcsa_intra_links();
        long long sum_best = 0;
        for (const auto& iv : intervals) {
            const int isize = iv.hi - iv.lo;
            if (isize >= kMinCoverage) {
                R.eligible_words++;
                R.eligible_positions += (size_t)isize;
                // Same call (name, lo, hi, require_avail=false) every
                // algorithm's own preference/DAG-edge computation makes --
                // see e.g. compress_pseudoforest_dp_'s forest-build step.
                int cov = best_candidate_(iv.name, iv.lo, iv.hi, /*require_avail=*/false).coverage();
                if (cov > 0) {
                    R.words_with_candidate++;
                    sum_best += cov;
                }
            }
        }
        R.floor_positions = m - R.eligible_positions;
        R.sum_best_coverage = (size_t)sum_best;
        R.lower_bound_C = m - (size_t)sum_best;
        return R;
    }

private:
    GappedSA G_;
    int span_ = 1;
    int max_add_ = 8;
    CompressAlgo algo_ = CompressAlgo::GreedySize;
    int phase2_iters_ = 0;  // 0 = unspecified; see run_phase2_ for precedence

    // Regular (ungapped) SA mode: no LexText, no mod-residue grouping.
    //   add_stride_ : characters shifted per unit of `add` (= shape.span in
    //                 LexText mode, since each lextext symbol packs
    //                 shape.span characters; = 1 in regular mode, since each
    //                 symbol is a single character).
    //   base_depth_ : symbols/characters required for a full name match
    //                 (= 1 in LexText mode -- DisLex's whole point is that
    //                 depth-1 already means "same gapped k-mer"; = k in
    //                 regular mode, where depth is measured in characters).
    bool regular_mode_ = false;
    int add_stride_ = 1;
    int base_depth_ = 1;

    std::vector<int64_t> C_;
    std::vector<int64_t> rank_to_C_;
    std::unordered_map<uint64_t, HashEntry> table_;

    std::vector<int32_t> rank_of_;
    std::vector<uint8_t> removed_;
    std::vector<uint64_t> removed_words_;
    std::vector<int32_t> pin_count_;   // #accepted candidates using this rank as source
    std::vector<uint64_t> pin_owner_;  // one owner name when pin_count==1 (undef if 0)
    std::unordered_map<int32_t, std::vector<uint64_t>> multi_pin_owners_; // owners when pin_count>1
    size_t kept_count_ = 0;            // #ranks with removed_==0
    int last_pin_path_ = 0;            // 0=none 1=fast 2=multi (debug/timing)

    // Undo record for a successful try_accept_with_retarget_ (surgical, not O(n)).
    struct RetargetUndo {
        bool had_old_self = false;
        Candidate old_self;
        std::vector<Candidate> old_deps;   // deps before retarget
        Candidate applied;
        std::vector<uint64_t> retargeted;  // names of old_deps / new deps
    };

    int64_t pred_lexpos_(int32_t r, int add) const {
        int64_t x = G_.sa[r];
        int64_t y = x - add;
        if (y < 0) return -1;
        if (G_.lex2orig[y] != G_.lex2orig[x] - (int64_t)add * add_stride_) return -1;
        return y;
    }

    // Inverse of pred_lexpos_: lextext position add symbols after sa[r], if the
    // original-text shift is exactly +add*add_stride_ (same residue class /
    // window in LexText mode; always true in regular mode, where every
    // position is its own "window").
    int64_t succ_lexpos_(int32_t r, int add) const {
        int64_t x = G_.sa[r];
        int64_t y = x + add;
        if (y >= (int64_t)G_.m()) return -1;
        if (G_.lex2orig[y] != G_.lex2orig[x] + (int64_t)add * add_stride_) return -1;
        return y;
    }

    bool run_kept_(int32_t a, int32_t b) const {
        if (a >= b) return true;
        if (removed_words_.empty() || (size_t)b > removed_words_.size() * 64) {
            for (int32_t r = a; r < b; ++r) if (removed_[r]) return false;
            return true;
        }
        size_t w_start = (size_t)a >> 6;
        size_t w_end   = (size_t)(b - 1) >> 6;
        uint64_t mask_start = ~0ULL << (a & 63);
        uint64_t mask_end   = ~0ULL >> (63 - ((b - 1) & 63));
        if (w_start == w_end) {
            return (removed_words_[w_start] & (mask_start & mask_end)) == 0;
        }
        if (removed_words_[w_start] & mask_start) return false;
        for (size_t w = w_start + 1; w < w_end; ++w) {
            if (removed_words_[w]) return false;
        }
        if (removed_words_[w_end] & mask_end) return false;
        return true;
    }

    // Enumerate all viable source runs for I_c = [lo,hi). Optionally require
    // sources kept and covered ranks free of pins (availability filter).
    // If only_add >= 1, enumerate that single add only (used by best_add1).
    // Never emits candidates with coverage < kMinCoverage (min_cov may raise it).
    std::vector<Candidate> enumerate_candidates_(uint64_t name, int32_t lo, int32_t hi,
                                                  bool require_avail,
                                                  int only_add = -1,
                                                  int min_cov = kMinCoverage) const {
        std::vector<Candidate> out;
        const int add_lo = (only_add >= 1) ? only_add : 1;
        const int add_hi = (only_add >= 1) ? only_add : max_add_;
        const int cov_floor = std::max(kMinCoverage, min_cov);
        const bool intra = gcsa_intra_links();
        const bool cross_lcp = gcsa_cross_lcp();
        // `out` is built and kept sorted coverage-desc, capped at `cap` --
        // emit_run_ skips materializing (and emit_intra_runs_'s scan skips
        // growing) any candidate that can't beat the current worst kept one,
        // instead of building every viable candidate and trimming after the
        // fact. See emit_run_ for why the cap check sits after availability.
        const int cap = gcsa_env_int("GCSA_CAND_CACHE_CAP", kCandCacheDefaultCap);
        for (int add = add_lo; add <= add_hi; ++add) {
            std::vector<std::pair<int32_t,int32_t>> pr;
            pr.reserve(hi - lo);
            for (int32_t r = lo; r < hi; ++r) {
                int64_t y = pred_lexpos_(r, add);
                if (y >= 0) pr.push_back({rank_of_[y], r});
            }
            if (pr.empty()) continue;
            std::sort(pr.begin(), pr.end());
            size_t i = 0;
            while (i < pr.size()) {
                size_t j = i + 1;
                // GCSA_CROSS_LCP=1: numeric adjacency alone is what L1-L5
                // require (see gcsa_cross_lcp() above); the lcp clause is an
                // extra, non-required restriction kept as the default.
                while (j < pr.size()
                       && pr[j].first == pr[j-1].first + 1
                       && (cross_lcp || G_.lcp[pr[j].first] >= add + base_depth_)) ++j;
                const int32_t s_lo = pr[i].first, s_hi = pr[j-1].first + 1;
                if (s_hi <= lo || s_lo >= hi)
                    emit_run_(name, lo, hi, add, pr, i, j, cov_floor,
                              require_avail, cap, out);
                else if (intra)
                    emit_intra_runs_(name, lo, hi, add, pr, i, j, cov_floor,
                                     require_avail, cap, out);
                i = j;
            }
        }
        return out;
    }

    // Turn the source sub-run pr[i..j) into a candidate for I_c = [lo,hi).
    // `out` is kept sorted coverage-desc (ties by add-desc) and capped at
    // `cap` entries -- the shared cache's cap (kCandCacheDefaultCap /
    // GCSA_CAND_CACHE_CAP), enforced here rather than by a separate trim
    // pass, so a candidate that can't make the cut is never materialized.
    // The cap check is a pure comparison against the current out.back(); it
    // never mutates `out` itself (only a confirmed insert below does), so
    // its position relative to the availability filter doesn't affect
    // correctness either way. It's placed first here -- an O(1) check --
    // to short-circuit before the O(window length) availability scan
    // whenever a window can't possibly beat the current worst kept one.
    void emit_run_(uint64_t name, int32_t lo, int32_t hi, int add,
                   const std::vector<std::pair<int32_t,int32_t>>& pr,
                   size_t i, size_t j, int cov_floor, bool require_avail,
                   int cap, std::vector<Candidate>& out) const {
        const int cov = (int)(j - i);
        if (cov < cov_floor) return;
        if ((int)out.size() >= cap) {
            const Candidate& worst = out.back();
            if (cov < worst.coverage() || (cov == worst.coverage() && add <= worst.add))
                return;
        }
        const int32_t s_lo = pr[i].first, s_hi = pr[j-1].first + 1;
        if (require_avail) {
            if (!run_kept_(s_lo, s_hi)) return;
            for (size_t t = i; t < j; ++t) {
                int32_t r = pr[t].second;
                if (removed_[r] || pin_count_[r] > 0) return;
            }
        }
        Candidate c;
        c.target_lo = lo; c.target_hi = hi; c.name = name;
        c.add = add; c.src_lo = s_lo; c.src_hi = s_hi;
        c.covered.reserve((size_t)cov);
        for (size_t t = i; t < j; ++t) c.covered.push_back(pr[t].second);
        auto pos = std::upper_bound(out.begin(), out.end(), c,
            [](const Candidate& a, const Candidate& b) {
                if (a.coverage() != b.coverage()) return a.coverage() > b.coverage();
                return a.add > b.add;
            });
        out.insert(pos, std::move(c));
        if ((int)out.size() > cap) out.pop_back();
    }

    // A source run reaching into I_c is still decodable as long as it holds
    // none of the ranks it covers (see "Link universe"), so instead of dropping
    // it we emit its maximal such sub-runs. Two-pointer scan in run-relative
    // coordinates: the window [p,q) is bad iff some t in [p,q) covers a rank in
    // [p,q), i.e. iff it contains both endpoints of one of the `forbid` pairs.
    void emit_intra_runs_(uint64_t name, int32_t lo, int32_t hi, int add,
                          const std::vector<std::pair<int32_t,int32_t>>& pr,
                          size_t i, size_t j, int cov_floor, bool require_avail,
                          int cap, std::vector<Candidate>& out) const {
        const int32_t len = (int32_t)(j - i);
        const int32_t s0 = pr[i].first;
        // forbid[b] = smallest window start that still traps a pair ending at b.
        std::vector<int32_t> forbid((size_t)len, -1);
        for (int32_t t = 0; t < len; ++t) {
            const int32_t u = pr[i + (size_t)t].second - s0;
            if (u < 0 || u >= len) continue;
            const int32_t a = std::min(t, u), b = std::max(t, u);
            forbid[(size_t)b] = std::max(forbid[(size_t)b], a + 1);
        }
        int32_t p = 0;
        for (int32_t q = 1; q <= len; ++q) {
            p = std::max(p, forbid[(size_t)(q - 1)]);
            // Only emit windows that cannot grow to the right without moving p;
            // p is non-decreasing, so no emitted window contains another.
            if (q == len || forbid[(size_t)q] > p)
                emit_run_(name, lo, hi, add, pr, i + (size_t)p, i + (size_t)q,
                          cov_floor, require_avail, cap, out);
        }
    }

    // Every decodable link of the whole instance (see "Link universe"), the
    // superset the ILP baseline optimizes over.  Source-centric: for a fixed
    // add the ranks split into maximal runs of consecutive ranks whose
    // add-successor exists and carries the same k-mer name; each sub-run of
    // length >= kMinCoverage that contains none of the ranks it covers is a
    // candidate.  No lcp >= add+1 and no maximality filter — those are the
    // heuristics' pruning rules, not correctness ones.
    //
    // Cost is O(m) per add plus O(L^2) per maximal run of length L (only the
    // ILP, i.e. tiny instances, uses this); sub_runs=false keeps it O(m).
    std::vector<Candidate> enumerate_all_links_(const std::vector<Interval>& intervals,
                                                bool sub_runs) const {
        const int32_t m = (int32_t)G_.m();
        std::unordered_map<uint64_t, Interval> by_name;
        by_name.reserve(intervals.size() * 2);
        for (const auto& iv : intervals) by_name[iv.name] = iv;

        std::vector<Candidate> out;
        std::vector<int32_t> tgt((size_t)m, -1);      // add-successor rank, -1 if none
        std::vector<int32_t> src_of((size_t)m, -1);   // inverse of tgt, per run
        for (int add = 1; add <= max_add_; ++add) {
            for (int32_t r = 0; r < m; ++r) {
                int64_t y = succ_lexpos_(r, add);
                tgt[(size_t)r] = (y < 0) ? -1 : rank_of_[(size_t)y];
            }
            int32_t a = 0;
            while (a < m) {
                if (tgt[(size_t)a] < 0) { ++a; continue; }
                const uint64_t name = name_of_rank_(tgt[(size_t)a]);
                int32_t b = a + 1;
                while (b < m && tgt[(size_t)b] >= 0
                       && name_of_rank_(tgt[(size_t)b]) == name) ++b;
                if (b - a >= kMinCoverage) {
                    for (int32_t t = a; t < b; ++t) src_of[(size_t)tgt[(size_t)t]] = t;
                    const Interval& iv = by_name[name];
                    const int32_t p_hi = sub_runs ? b - kMinCoverage : a;
                    for (int32_t p = a; p <= p_hi; ++p) {
                        // Grow the source run one rank at a time; `clash` marks
                        // the first q where [p,q] contains one of its own
                        // covered ranks (monotone in q, so we can stop there).
                        bool clash = false;
                        for (int32_t q = p; q < b; ++q) {
                            if (tgt[(size_t)q] >= p && tgt[(size_t)q] <= q) clash = true;
                            int32_t w = src_of[(size_t)q];
                            if (w >= p && w < q) clash = true;
                            if (clash) break;
                            const int32_t len = q + 1 - p;
                            if (len < kMinCoverage) continue;
                            if (!sub_runs && q + 1 != b) continue;
                            Candidate c;
                            c.name = name;
                            c.add = add;
                            c.src_lo = p;
                            c.src_hi = q + 1;
                            c.target_lo = iv.lo;
                            c.target_hi = iv.hi;
                            c.covered.reserve((size_t)len);
                            for (int32_t t = p; t <= q; ++t)
                                c.covered.push_back(tgt[(size_t)t]);
                            out.push_back(std::move(c));
                        }
                    }
                    for (int32_t t = a; t < b; ++t) src_of[(size_t)tgt[(size_t)t]] = -1;
                }
                a = b;
            }
        }
        return out;
    }

    static Candidate pick_best_(const std::vector<Candidate>& cs) {
        Candidate best;
        for (const auto& c : cs) {
            if (c.coverage() > best.coverage()
                || (c.coverage() == best.coverage() && c.add > best.add))
                best = c;
        }
        return best;
    }

    Candidate best_candidate_(uint64_t name, int32_t lo, int32_t hi, bool require_avail) const {
        return pick_best_(enumerate_candidates_(name, lo, hi, require_avail));
    }

    // Availability filter matching enumerate_candidates_(..., require_avail=true).
    bool cand_available_(const Candidate& c) const {
        if (!run_kept_(c.src_lo, c.src_hi)) return false;
        for (int32_t r : c.covered)
            if (removed_[r] || pin_count_[r] > 0) return false;
        return true;
    }

    // First viable candidate from a coverage-desc / add-desc sorted list.
    Candidate best_from_sorted_cache_(const std::vector<Candidate>& cands,
                                      bool require_avail) const {
        for (const auto& c : cands) {
            if (c.coverage() < kMinCoverage) continue;
            if (require_avail && !cand_available_(c)) continue;
            return c;
        }
        return Candidate{};
    }

    // The k-mer / gapped-k-mer name occupying SA rank r. In LexText mode this
    // is a cheap array lookup (the name was precomputed into G_.lex). In
    // regular mode G_.lex holds raw characters, not names, so the name is
    // recomputed on demand from the text -- called only once per interval
    // (not once per rank), so this stays cheap in aggregate.
    uint64_t name_of_rank_(int32_t r) const {
        return regular_mode_ ? name_at(G_.shape, G_.text, (size_t)G_.orig_pos(r))
                              : G_.first_symbol(r);
    }

    // ---- accept / revoke helpers ------------------------------------------
    void accept_(Candidate& c, std::unordered_map<uint64_t, Candidate>& accepted) {
        for (int32_t r = c.src_lo; r < c.src_hi; ++r) {
            if (pin_count_[r] == 0) {
                pin_owner_[r] = c.name;
            } else if (pin_count_[r] == 1) {
                multi_pin_owners_[r] = {pin_owner_[r], c.name};
            } else {
                multi_pin_owners_[r].push_back(c.name);
            }
            ++pin_count_[r];
        }
        for (int32_t r : c.covered) {
            if (!removed_[r]) {
                removed_[r] = 1;
                if (!removed_words_.empty()) removed_words_[r >> 6] |= (1ULL << (r & 63));
                --kept_count_;
            } else {
                removed_[r] = 1;
                if (!removed_words_.empty()) removed_words_[r >> 6] |= (1ULL << (r & 63));
            }
        }
        accepted[c.name] = std::move(c);
    }

    void revoke_(uint64_t name, std::unordered_map<uint64_t, Candidate>& accepted) {
        auto it = accepted.find(name);
        if (it == accepted.end()) return;
        Candidate& c = it->second;
        for (int32_t r = c.src_lo; r < c.src_hi; ++r) {
            --pin_count_[r];
            if (pin_count_[r] == 0) {
                pin_owner_[r] = 0;
                multi_pin_owners_.erase(r);
            } else {
                auto mit = multi_pin_owners_.find(r);
                if (mit != multi_pin_owners_.end()) {
                    auto& vec = mit->second;
                    vec.erase(std::remove(vec.begin(), vec.end(), name), vec.end());
                    if (pin_count_[r] == 1) {
                        if (!vec.empty()) pin_owner_[r] = vec.front();
                        multi_pin_owners_.erase(mit);
                    }
                } else if (pin_owner_[r] == name) {
                    pin_owner_[r] = 0;
                }
            }
        }
        for (int32_t r : c.covered) {
            if (removed_[r]) {
                removed_[r] = 0;
                if (!removed_words_.empty()) removed_words_[r >> 6] &= ~(1ULL << (r & 63));
                ++kept_count_;
            }
        }
        accepted.erase(it);
    }

    void undo_retarget_(RetargetUndo& u,
                        std::unordered_map<uint64_t, Candidate>& accepted) {
        revoke_(u.applied.name, accepted);
        for (uint64_t nm : u.retargeted) revoke_(nm, accepted);
        for (auto& d : u.old_deps) accept_(d, accepted);
        if (u.had_old_self) accept_(u.old_self, accepted);
    }

    // Can we retarget dependent D through new candidate W_cand?
    // D.src must be a contiguous sub-run of W_cand.covered; then D can point at
    // the corresponding sub-run of W_cand.src with add' = D.add + W_cand.add.
    static bool can_retarget_(const Candidate& dep, const Candidate& via,
                              int32_t& new_src_lo, int32_t& new_src_hi) {
        // Find dep.src_lo as an element of via.covered.
        auto& cov = via.covered;
        auto it = std::find(cov.begin(), cov.end(), dep.src_lo);
        if (it == cov.end()) return false;
        size_t t0 = (size_t)(it - cov.begin());
        size_t len = (size_t)(dep.src_hi - dep.src_lo);
        if (t0 + len > cov.size()) return false;
        // Contiguous in covered[] and equal to [src_lo, src_hi) as ranks.
        for (size_t t = 0; t < len; ++t) {
            if (cov[t0 + t] != dep.src_lo + (int32_t)t) return false;
        }
        new_src_lo = via.src_lo + (int32_t)t0;
        new_src_hi = new_src_lo + (int32_t)len;
        return true;
    }

    // Try to accept `cand` for its target, retargeting any dependents that
    // currently pin ranks inside cand.covered. Returns true on success.
    // If undo_out != nullptr, fills a surgical undo record (no O(n) snapshots).
    bool try_accept_with_retarget_(const Candidate& cand,
                                   std::unordered_map<uint64_t, Candidate>& accepted,
                                   const std::vector<Interval>& /*intervals*/,
                                   RetargetUndo* undo_out = nullptr) {
        // Covered ranks must not already be removed.
        for (int32_t r : cand.covered) if (removed_[r]) return false;
        // Source must be kept.
        if (!run_kept_(cand.src_lo, cand.src_hi)) return false;

        // Collect dependents whose source intersects cand.covered.
        std::vector<int32_t> pinned;
        pinned.reserve(8);
        bool multi_pin = false;
        for (int32_t r : cand.covered) {
            if (pin_count_[r] == 0) continue;
            pinned.push_back(r);
            if (pin_count_[r] > 1) multi_pin = true;
        }
        last_pin_path_ = pinned.empty() ? 0 : (multi_pin ? 2 : 1);

        std::vector<uint64_t> to_retarget;
        if (!pinned.empty() && !multi_pin) {
            // Fast path: each pinned covered rank has a unique owner.
            std::unordered_set<uint64_t> deps;
            for (int32_t r : pinned) deps.insert(pin_owner_[r]);
            to_retarget.assign(deps.begin(), deps.end());
            for (uint64_t nm : to_retarget) {
                auto it = accepted.find(nm);
                if (it == accepted.end()) return false;
                int32_t nlo, nhi;
                if (!can_retarget_(it->second, cand, nlo, nhi)) return false;
            }
        } else if (!pinned.empty()) {
            // Multi-pin: look up owners via pin_owner_ / multi_pin_owners_.
            std::unordered_set<uint64_t> deps;
            for (int32_t r : pinned) {
                if (pin_count_[r] == 1) {
                    deps.insert(pin_owner_[r]);
                } else if (pin_count_[r] > 1) {
                    auto mit = multi_pin_owners_.find(r);
                    if (mit != multi_pin_owners_.end()) {
                        for (uint64_t nm : mit->second) deps.insert(nm);
                    } else {
                        for (const auto& kv : accepted) {
                            if (r >= kv.second.src_lo && r < kv.second.src_hi) {
                                deps.insert(kv.first);
                            }
                        }
                    }
                }
            }
            // Every pin on a covered rank must come from a collected dep.
            for (int32_t r : pinned) {
                int seen = 0;
                for (uint64_t nm : deps) {
                    auto it = accepted.find(nm);
                    if (it == accepted.end()) continue;
                    const Candidate& d = it->second;
                    if (r >= d.src_lo && r < d.src_hi) ++seen;
                }
                if (seen != pin_count_[r]) return false;
            }
            to_retarget.assign(deps.begin(), deps.end());
            for (uint64_t nm : to_retarget) {
                auto it = accepted.find(nm);
                if (it == accepted.end()) return false;
                int32_t nlo, nhi;
                if (!can_retarget_(it->second, cand, nlo, nhi)) return false;
            }
        }

        // Snapshot surgical undo state before mutating.
        if (undo_out) {
            undo_out->had_old_self = false;
            undo_out->old_deps.clear();
            undo_out->retargeted = to_retarget;
            auto it_old = accepted.find(cand.name);
            if (it_old != accepted.end()) {
                undo_out->had_old_self = true;
                undo_out->old_self = it_old->second;
            }
            for (uint64_t nm : to_retarget)
                undo_out->old_deps.push_back(accepted[nm]);
            undo_out->applied = cand;
        }

        // If this name already has an offset, revoke it first (we'll replace).
        if (accepted.count(cand.name)) revoke_(cand.name, accepted);

        // Retarget dependents: revoke + re-accept with updated src/add.
        std::vector<Candidate> retargeted;
        retargeted.reserve(to_retarget.size());
        for (uint64_t nm : to_retarget) {
            Candidate d = accepted[nm];
            int32_t nlo, nhi;
            can_retarget_(d, cand, nlo, nhi);
            revoke_(nm, accepted);
            d.src_lo = nlo;
            d.src_hi = nhi;
            d.add += cand.add;
            retargeted.push_back(std::move(d));
        }

        Candidate applied = cand;
        accept_(applied, accepted);
        for (auto& d : retargeted) accept_(d, accepted);
        if (undo_out) undo_out->applied = accepted[cand.name];
        return true;
    }

    // Phase II iteration budget, in precedence order: --phase2-iters (the
    // caller's phase2_iters) if > 0, else GCSA_PHASE2_MAX_ITERS, else
    // GCSA_PHASE2_FAST / auto-large-m (1 gen), else kPhase2DefaultIters.
    // Shared by both Phase II implementations.
    int phase2_budget_(const char*& src) const {
        int iters = kPhase2DefaultIters;
        src = "default";
        if (G_.m() >= kPhase2AutoFastM) {
            iters = 1;
            src = "auto-large-m";
        }
        if (const char* env = std::getenv("GCSA_PHASE2_FAST")) {
            if (std::atoi(env) != 0) {
                iters = 1;
                src = "GCSA_PHASE2_FAST";
            }
        }
        if (const char* env = std::getenv("GCSA_PHASE2_MAX_ITERS")) {
            int v = std::atoi(env);
            if (v > 0) { iters = v; src = "GCSA_PHASE2_MAX_ITERS"; }
        }
        if (phase2_iters_ > 0) {
            iters = phase2_iters_;
            src = "--phase2-iters";
        }
        return iters;
    }

    // Truncate a coverage-desc candidate list for the shared cache.
    static void trim_cand_cache_(std::vector<Candidate>& cands) {
        int cap = gcsa_env_int("GCSA_CAND_CACHE_CAP", kCandCacheDefaultCap);
        if (cap > 0 && (int)cands.size() > cap) cands.resize((size_t)cap);
    }

    // ---- Phase II (local): exact cluster large-neighborhood search ---------
    //
    // Model. |C| = m - sum of the coverages of the accepted links, so we are
    // maximizing total coverage subject to exactly two constraints:
    //   (F1) at most one link per name (one H_offset per hash entry), and
    //   (F2) a chosen link's source ranks are all kept, i.e. no *other* chosen
    //        link covers them.
    // There is no acyclicity requirement (C holds literal positions, H_offset
    // never recurses), so this is a maximum-weight set packing, and crucially
    // (F2) is a *pairwise* condition between links.
    //
    // Neighborhood. Freeze the links of all names outside a small cluster S and
    // re-solve S exactly. Two names are dependent iff a candidate of one draws
    // its source from the other's interval — that is the only way (F2) can bind
    // them. Note that *sharing* a source is deliberately not a dependency:
    // sources are shared freely (pin_count_ is a count, not a lock), so two
    // links reading the same ranks never conflict.
    //
    // Composition. Covered ranks of a link always lie in its own interval, so a
    // frozen link can never drop a rank inside S — the cluster fully owns the
    // keep/drop decisions for its own ranks. The interface to the frozen part is
    // therefore just: (i) don't drop a rank some frozen link uses as a source,
    // and (ii) don't source from a rank some frozen link drops. After revoking
    // the cluster's own links, removed_ / pin_count_ describe exactly the frozen
    // part, so (i)+(ii) is precisely the existing cand_available_ predicate.
    // Internal conflicts are then the pairwise (F2) checks during enumeration.
    //
    // Since the incumbent assignment is itself feasible and is seeded as the
    // initial best, an accepted cluster solve can never increase |C|.
    void run_phase2_local_(const std::vector<Interval>& intervals,
                           std::unordered_map<uint64_t, Candidate>& accepted,
                           const char* label, bool timing,
                           std::unordered_map<uint64_t, std::vector<Candidate>>* precomputed) {
        using Clock = std::chrono::steady_clock;
        const bool trace = (std::getenv("GCSA_TRACE_LNS") != nullptr);
        const char* iters_src = "default";
        const int max_iters = phase2_budget_(iters_src);
        // 0 = pick from the measured dependency-graph density once it is built.
        int cluster_cap = gcsa_env_int("GCSA_LNS_CLUSTER", 0);
        const int opt_cap = gcsa_env_int("GCSA_LNS_OPTS", 12);
        const int deg_cap = gcsa_env_int("GCSA_LNS_DEGREE", 16);
        long node_cap = gcsa_env_int("GCSA_LNS_NODES", 0);
        auto t0 = Clock::now();
        gcsa_log("[%s] Phase II: cluster LNS...\n", label);

        std::unordered_map<uint64_t, std::vector<Candidate>> local_cache;
        auto& cand_cache = precomputed ? *precomputed : local_cache;

        // Seeds in size order, matching the greedy Phase II's bias to big wins.
        std::vector<const Interval*> order;
        for (const auto& iv : intervals) if (iv.hi - iv.lo > 1) order.push_back(&iv);
        std::sort(order.begin(), order.end(),
                  [](const Interval* a, const Interval* b){
                      return (a->hi - a->lo) > (b->hi - b->lo);
                  });
        for (const Interval* ivp : order) {
            if (cand_cache.count(ivp->name)) continue;
            auto cands = enumerate_candidates_(ivp->name, ivp->lo, ivp->hi, /*avail=*/false);
            std::sort(cands.begin(), cands.end(),
                      [](const Candidate& a, const Candidate& b){
                          if (a.coverage() != b.coverage()) return a.coverage() > b.coverage();
                          return a.add > b.add;
                      });
            cand_cache.emplace(ivp->name, std::move(cands));
        }

        // Dependency graph over names: an edge for every "c can source from d".
        std::unordered_map<uint64_t, std::vector<uint64_t>> adj;
        adj.reserve(order.size() * 2);
        auto add_edge = [&](uint64_t a, uint64_t b) {
            if (a == b) return;
            auto& va = adj[a];
            if (std::find(va.begin(), va.end(), b) == va.end()) va.push_back(b);
            auto& vb = adj[b];
            if (std::find(vb.begin(), vb.end(), a) == vb.end()) vb.push_back(a);
        };
        for (const Interval* ivp : order) {
            const auto& cands = cand_cache[ivp->name];
            for (const Candidate& c : cands) {
                uint64_t last = 0;
                bool have_last = false;
                for (int32_t r = c.src_lo; r < c.src_hi; ++r) {
                    uint64_t s = name_of_rank_(r);
                    if (have_last && s == last) continue;
                    last = s; have_last = true;
                    add_edge(ivp->name, s);
                }
            }
        }
        for (auto& kv : adj)
            if ((int)kv.second.size() > deg_cap) kv.second.resize((size_t)deg_cap);

        // Best coverage any single link could ever give a name (cand_cache is
        // coverage-desc). Summed over a cluster this upper-bounds every feasible
        // assignment, so it cheaply rules out clusters that cannot improve
        // before we pay for revoking them and filtering their options.
        // A bigger --max-add densifies this graph, so a fixed-size BFS ball
        // covers a progressively smaller share of each name's real neighborhood.
        // Size the cluster from the measured mean degree instead of a constant,
        // and grow the enumeration budget with it: a cluster that overruns the
        // budget falls back to 4 names, which costs more than the larger cluster
        // ever gained. Both remain overridable.
        // Both are quality/cost dials rather than a free win: raising them
        // together helps at every --max-add, but costs 30-80x runtime, so the
        // defaults stay put and GCSA_LNS_AUTO=1 opts in.
        size_t deg_sum = 0;
        for (const auto& kv : adj) deg_sum += kv.second.size();
        const double mean_deg =
            adj.empty() ? 0.0 : (double)deg_sum / (double)adj.size();
        const bool auto_size = (gcsa_env_int("GCSA_LNS_AUTO", 0) != 0);
        if (cluster_cap <= 0)
            cluster_cap = auto_size
                ? std::min(16, std::max(8, (int)std::lround(2.0 * mean_deg)))
                : 8;
        if (node_cap <= 0)
            node_cap = 20000L * (1 + 3 * (cluster_cap - 8));

        std::unordered_map<uint64_t, int> static_cap;
        static_cap.reserve(cand_cache.size() * 2);
        for (const auto& kv : cand_cache)
            static_cap[kv.first] = kv.second.empty() ? 0 : kv.second.front().coverage();
        auto t_adj = Clock::now();

        // Internal (F2) check between two links of a cluster. A link's source is
        // a contiguous rank range and everything it covers lies in its own
        // interval, so the ranges almost never meet and the test is O(1); only
        // on overlap do we look at individual covered ranks.
        auto conflicts = [](const Candidate& a, const Candidate& b) {
            if (a.src_lo < b.target_hi && b.target_lo < a.src_hi)
                for (int32_t r : b.covered)
                    if (r >= a.src_lo && r < a.src_hi) return true;
            if (b.src_lo < a.target_hi && a.target_lo < b.src_hi)
                for (int32_t r : a.covered)
                    if (r >= b.src_lo && r < b.src_hi) return true;
            return false;
        };

        std::vector<uint64_t> dirty;
        dirty.reserve(order.size());
        for (const Interval* ivp : order) dirty.push_back(ivp->name);
        std::unordered_map<uint64_t, const Interval*> by_name;
        by_name.reserve(intervals.size() * 2);
        for (const auto& iv : intervals) by_name[iv.name] = &iv;

        size_t clusters = 0, improved = 0, aborted = 0, skipped = 0, gain_total = 0;
        int iter = 0;
        const char* stop = "fixed-point";
        for (; iter < max_iters && !dirty.empty(); ++iter) {
            std::vector<uint64_t> next;
            std::unordered_set<uint64_t> next_seen;
            std::unordered_set<uint64_t> done;
            for (uint64_t seed : dirty) {
                if (!done.insert(seed).second) continue;

                // --- grow the cluster: BFS over the dependency graph ---------
                std::vector<uint64_t> cluster{seed};
                for (size_t qi = 0;
                     qi < cluster.size() && (int)cluster.size() < cluster_cap; ++qi) {
                    auto ait = adj.find(cluster[qi]);
                    if (ait == adj.end()) continue;
                    for (uint64_t nb : ait->second) {
                        if ((int)cluster.size() >= cluster_cap) break;
                        if (!by_name.count(nb) || !cand_cache.count(nb)) continue;
                        if (std::find(cluster.begin(), cluster.end(), nb) != cluster.end())
                            continue;
                        cluster.push_back(nb);
                    }
                }
                int ub = 0, incumbent = 0;
                for (uint64_t nm : cluster) {
                    auto sc = static_cap.find(nm);
                    if (sc != static_cap.end()) ub += sc->second;
                    auto ac = accepted.find(nm);
                    if (ac != accepted.end()) incumbent += ac->second.coverage();
                }
                if (ub <= incumbent) { ++skipped; continue; }

                // Solve S exactly; returns the coverage gain, or -1 if the node
                // cap was hit (state restored, cluster left untouched).
                auto solve_cluster = [&](const std::vector<uint64_t>& S) -> int {
                ++clusters;

                // --- the cluster's current links -----------------------------
                std::vector<Candidate> cur(S.size());
                std::vector<uint8_t> had(S.size(), 0);
                int cur_total = 0;
                for (size_t i = 0; i < S.size(); ++i) {
                    auto it = accepted.find(S[i]);
                    if (it == accepted.end()) continue;
                    had[i] = 1;
                    cur[i] = it->second;
                    cur_total += cur[i].coverage();
                }

                // Availability against the frozen part, evaluated *without*
                // revoking anything: the solve usually finds no improvement, and
                // revoke_/accept_ round trips are the dominant cost on large
                // inputs. Revoking would restore the ranks the cluster's own
                // links cover and release the pins they hold, so subtract both
                // from removed_ / pin_count_ on the fly.
                std::vector<int32_t> cl_cov;
                for (size_t i = 0; i < S.size(); ++i)
                    if (had[i])
                        cl_cov.insert(cl_cov.end(), cur[i].covered.begin(),
                                      cur[i].covered.end());
                std::sort(cl_cov.begin(), cl_cov.end());
                auto kept_after_revoke = [&](int32_t r) {
                    return !removed_[r]
                        || std::binary_search(cl_cov.begin(), cl_cov.end(), r);
                };
                auto external_pins = [&](int32_t r) {
                    int n = pin_count_[r];
                    for (size_t i = 0; i < S.size(); ++i)
                        if (had[i] && r >= cur[i].src_lo && r < cur[i].src_hi) --n;
                    return n;
                };
                auto avail = [&](const Candidate& c) {
                    for (int32_t r = c.src_lo; r < c.src_hi; ++r)
                        if (!kept_after_revoke(r)) return false;
                    for (int32_t r : c.covered)
                        if (!kept_after_revoke(r) || external_pins(r) > 0) return false;
                    return true;
                };

                // --- options per member, filtered against the frozen part ----
                std::vector<std::vector<const Candidate*>> opts(S.size());
                for (size_t i = 0; i < S.size(); ++i) {
                    const auto& cands = cand_cache[S[i]];
                    for (const Candidate& c : cands) {
                        if ((int)opts[i].size() >= opt_cap) break;
                        if (c.coverage() < kMinCoverage) continue;
                        if (!avail(c)) continue;
                        opts[i].push_back(&c);
                    }
                    // The incumbent must stay expressible even if the cap or the
                    // cache would hide it, else "never worse" is not guaranteed.
                    if (had[i]) {
                        bool found = false;
                        for (const Candidate* c : opts[i])
                            if (c->add == cur[i].add && c->src_lo == cur[i].src_lo
                                && c->src_hi == cur[i].src_hi) { found = true; break; }
                        if (!found) opts[i].push_back(&cur[i]);
                    }
                }

                // Members ordered by best-option coverage: fail fast, bound hard.
                std::vector<int> ord(S.size());
                for (size_t i = 0; i < S.size(); ++i) ord[i] = (int)i;
                std::vector<int> cap_cov(S.size(), 0);
                for (size_t i = 0; i < S.size(); ++i)
                    for (const Candidate* c : opts[i])
                        cap_cov[i] = std::max(cap_cov[i], c->coverage());
                std::sort(ord.begin(), ord.end(), [&](int a, int b){
                    return cap_cov[(size_t)a] > cap_cov[(size_t)b];
                });
                std::vector<int> suffix(S.size() + 1, 0);
                for (size_t k = S.size(); k-- > 0; )
                    suffix[k] = suffix[k + 1] + cap_cov[(size_t)ord[k]];

                // --- exact solve: DFS over one option (or none) per member ---
                std::vector<const Candidate*> pick(S.size(), nullptr);
                std::vector<const Candidate*> best(S.size(), nullptr);
                for (size_t i = 0; i < S.size(); ++i)
                    if (had[i]) best[i] = &cur[i];
                int best_total = cur_total;
                long nodes = 0;
                bool over_budget = false;

                std::vector<const Candidate*> live;
                live.reserve(S.size());
                std::function<void(size_t, int)> dfs = [&](size_t k, int total) {
                    if (over_budget) return;
                    if (++nodes > node_cap) { over_budget = true; return; }
                    if (total + suffix[k] <= best_total) return;
                    if (k == S.size()) {
                        best_total = total;
                        best = pick;
                        return;
                    }
                    const size_t i = (size_t)ord[k];
                    for (const Candidate* c : opts[i]) {
                        bool clash = false;
                        for (const Candidate* x : live)
                            if (conflicts(*c, *x)) { clash = true; break; }
                        if (clash) continue;
                        live.push_back(c);
                        pick[i] = c;
                        dfs(k + 1, total + c->coverage());
                        pick[i] = nullptr;
                        live.pop_back();
                        if (over_budget) return;
                    }
                    dfs(k + 1, total);  // leave this name uncompressed
                };
                dfs(0, 0);

                if (over_budget) return -1;
                if (best_total == cur_total) return 0;  // incumbent stands

                // --- apply: only now do we touch the global state ------------
                for (size_t i = 0; i < S.size(); ++i)
                    if (had[i]) revoke_(S[i], accepted);
                for (size_t i = 0; i < S.size(); ++i) {
                    if (!best[i]) continue;
                    Candidate c = *best[i];
                    accept_(c, accepted);
                }
                if (trace) {
                    gcsa_log("  LNS cluster seed=%s size=%zu cov %d -> %d\n",
                             name_to_string(G_.shape, seed).c_str(), S.size(),
                             cur_total, best_total);
                }
                return best_total - cur_total;
                };

                // A cluster too tangled to enumerate degrades to a smaller one
                // rather than being skipped: aborting outright loses the
                // improvements the sub-cluster would still have found.
                int gain = solve_cluster(cluster);
                if (gain < 0 && (int)cluster.size() > kLnsFallbackCluster) {
                    ++aborted;
                    cluster.resize((size_t)kLnsFallbackCluster);
                    gain = solve_cluster(cluster);
                }
                if (gain < 0) { ++aborted; continue; }

                if (gain > 0) {
                    ++improved;
                    gain_total += (size_t)gain;
                    for (uint64_t nm : cluster) {
                        if (next_seen.insert(nm).second) next.push_back(nm);
                        auto ait = adj.find(nm);
                        if (ait == adj.end()) continue;
                        for (uint64_t nb : ait->second)
                            if (cand_cache.count(nb) && next_seen.insert(nb).second)
                                next.push_back(nb);
                    }
                }
            }
            dirty.swap(next);
        }
        if (iter >= max_iters && !dirty.empty()) stop = "max-iters";

        if (timing) {
            auto ms = [](Clock::time_point a, Clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count();
            };
            gcsa_log("[timing] %s Phase II(LNS) graph=%.1fms total=%.1fms\n"
                     "         gen %d/%d (%s, stop=%s) solved=%zu skipped=%zu "
                     "improved=%zu aborted=%zu coverage_gain=%zu kept=%zu\n"
                     "         mean_degree=%.2f cluster_cap=%d node_cap=%ld\n",
                     label, ms(t0, t_adj), ms(t0, Clock::now()),
                     iter, max_iters, iters_src, stop, clusters, skipped,
                     improved, aborted, gain_total, kept_count_,
                     mean_deg, cluster_cap, node_cap);
        }
    }

    // Dirty-set unpin / retarget: try higher-coverage candidates (including
    // former sources), retargeting dependents with transitive add when needed.
    // Shared by every algorithm's Phase II (see run_phase2_and_lns_ below,
    // which pairs this with the cluster LNS pass).
    //
    // Candidates with require_avail=false depend only on SA/LCP structure, so we
    // enumerate once per interval and reuse across generations. Optional
    // `precomputed` cache (e.g. from an algorithm's own preference build)
    // avoids re-enum.
    //
    // Work queue (generation dirty-set):
    //   seed dirty = all |I|>1; process in size order.
    //   On any improve in a generation, next_dirty = intervals from this
    //   generation that can still beat their coverage (pin-blocked try_accept
    //   failures, exhausted tries with best-cached > cur, or suboptimal
    //   accepts). Avoids re-scanning skip_best / saturated names.
    //   Each interval is processed at most once per generation; generations
    //   stop when next_dirty is empty (fixed point) or the iteration budget is
    //   spent.
    //
    // Phase II runs a fixed number of dirty generations. The budget is, in
    // precedence order: the caller's phase2_iters (--phase2-iters) if > 0, else
    // GCSA_PHASE2_MAX_ITERS=K, else kPhase2DefaultIters (100). Reaching a fixed
    // point stops earlier; a budget below the fixed point leaves |C| larger.
    //
    // Adaptive early-stop (cuts long REPLACE-only plateaus), off by default:
    //   GCSA_PHASE2_STALL=S     – consecutive generations with kept-drop < G
    //                             before stopping. Default 0 = disabled.
    //   GCSA_PHASE2_MIN_GAIN=G  – min total kept-drop (|C| reduction) per dirty
    //                             generation to count as progress (default 1).
    //
    // Env GCSA_DISABLE_PHASE2 (optional): skip the greedy retarget pass (the
    // LNS pass in run_phase2_and_lns_ still runs regardless), so any
    // algorithm can be compared on its Phase I / DP output alone, without
    // the shared retarget sweep smoothing over differences between them.
    void run_phase2_(const std::vector<Interval>& intervals,
                     std::unordered_map<uint64_t, Candidate>& accepted,
                     const char* label,
                     bool trace,
                     bool timing,
                     std::unordered_map<uint64_t, std::vector<Candidate>>* precomputed = nullptr,
                     double phase1_ms = 0.0) {
        using Clock = std::chrono::steady_clock;
        if (std::getenv("GCSA_DISABLE_PHASE2") != nullptr) {
            std::fprintf(stderr, "[%s] Phase II: disabled (GCSA_DISABLE_PHASE2)\n", label);
            return;
        }
        const char* iters_src = "default";
        const int max_iters = phase2_budget_(iters_src);
        int min_gain = 1;
        if (const char* env = std::getenv("GCSA_PHASE2_MIN_GAIN")) {
            int v = std::atoi(env);
            if (v >= 0) min_gain = v;
        }
        int stall_limit = 0;
        if (const char* env = std::getenv("GCSA_PHASE2_STALL")) {
            int v = std::atoi(env);
            if (v >= 0) stall_limit = v;
        }
        const bool adaptive = (stall_limit > 0);

        double time_ratio = 1.0;
        if (const char* env = std::getenv("GCSA_PHASE2_TIME_RATIO")) {
            double v = std::atof(env);
            if (v > 0.0) time_ratio = v;
        }
        double min_ms = 50.0;
        if (const char* env = std::getenv("GCSA_PHASE2_MIN_MS")) {
            double v = std::atof(env);
            if (v >= 0.0) min_ms = v;
        }

        bool fast_mode = (G_.m() >= kPhase2AutoFastM);
        if (const char* env = std::getenv("GCSA_PHASE2_FAST")) {
            fast_mode = (std::atoi(env) != 0);
        }

        // Always time-budget against Phase I when its duration is known; FAST
        // alone still applies the min_ms floor. Caps Phase II ≈ Phase I wall.
        const bool time_budgeted = (fast_mode || phase1_ms > 0.0);
        const double time_limit_ms = time_budgeted
            ? std::max(min_ms, (phase1_ms > 0.0 ? phase1_ms : min_ms) * time_ratio)
            : 0.0;

        gcsa_log("[%s] Phase II: unpin/retarget...\n", label);
        if (timing) {
            if (precomputed) {
                gcsa_log("[timing] %s Phase II precomputed cache size=%zu\n",
                         label, precomputed->size());
            }
            if (time_budgeted) {
                gcsa_log("[timing] %s Phase II time_budget=%.1fms "
                         "(phase1=%.1fms ratio=%.2f fast=%d)\n",
                         label, time_limit_ms, phase1_ms, time_ratio,
                         (int)fast_mode);
            }
        }
        auto t0 = Clock::now();

        // Size-order once; static candidate lists (avail ignored) cached once.
        // |I|<=2 can never meet kMinCoverage, so omit them from the dirty seed.
        std::vector<const Interval*> order;
        for (const auto& iv : intervals) if (iv.hi - iv.lo > 2) order.push_back(&iv);
        std::sort(order.begin(), order.end(),
                  [](const Interval* a, const Interval* b){
                      return (a->hi-a->lo) > (b->hi-b->lo);
                  });
        const int n_ord = (int)order.size();

        std::unordered_map<uint64_t, std::vector<Candidate>> local_cache;
        auto& cand_cache = precomputed ? *precomputed : local_cache;
        size_t phase2_enum = 0, phase2_cache_hits = 0, phase2_skip_sat = 0;
        size_t phase2_skip_best = 0, phase2_tries = 0;
        size_t phase2_fail = 0, phase2_improve = 0;
        size_t phase2_fast_pin = 0, phase2_multi_pin = 0, phase2_no_pin = 0;
        size_t phase2_dirty_marks = 0;
        double ms_enum = 0, ms_try = 0;

        // Locked for the whole call, not just the cache lookup: the
        // per-generation parallel evaluate pass below calls this
        // concurrently (gcsa_parallel_for(dirty.size(), ...) -> get_cands),
        // and unordered_map modification -- the insert on a miss -- is a
        // data race across threads even for distinct keys; only individual
        // element *references* are guaranteed stable across later
        // inserts/rehashes, concurrent modification of the container itself
        // is not. Returning a reference (not a copy) matters: callers keep
        // raw pointers into the returned vector's elements past this call
        // (see "proposed[idx] = {oi, &c}" in the parallel pass below, read
        // back in the sequential decide loop after this generation's
        // parallel pass has fully returned) -- those pointers are only
        // valid because they point into cand_cache's own long-lived
        // storage, not a temporary. Misses are expected to be rare by the
        // time Phase II's generation loop runs (cand_cache is pre-warmed
        // either by the one-time pass above or, for pseudoforest-dp, by the
        // caller's own parallel candidate-graph build), so serializing the
        // miss path too -- rather than dropping the lock around
        // enumerate_candidates_ -- is an acceptable, correctness-first
        // tradeoff.
        std::mutex cand_mu;
        auto get_cands = [&](const Interval& iv) -> const std::vector<Candidate>& {
            std::lock_guard<std::mutex> lock(cand_mu);
            auto it = cand_cache.find(iv.name);
            if (it != cand_cache.end()) {
                ++phase2_cache_hits;
                return it->second;
            }
            auto te0 = Clock::now();
            auto cands = enumerate_candidates_(iv.name, iv.lo, iv.hi,
                                               /*require_avail=*/false);
            if (timing) ms_enum += std::chrono::duration<double, std::milli>(Clock::now() - te0).count();
            ++phase2_enum;
            std::sort(cands.begin(), cands.end(), [](const Candidate& a, const Candidate& b){
                if (a.coverage() != b.coverage()) return a.coverage() > b.coverage();
                return a.add > b.add;
            });
            trim_cand_cache_(cands);
            auto& slot = cand_cache[iv.name];
            slot = std::move(cands);
            return slot;
        };

        auto cov_of = [&](uint64_t name) -> int {
            auto it = accepted.find(name);
            return it == accepted.end() ? 0 : it->second.coverage();
        };

        if (cand_cache.size() < (size_t)n_ord) {
            std::vector<std::vector<Candidate>> temp_cands((size_t)n_ord);
            gcsa_parallel_for((size_t)n_ord, [&](size_t i) {
                const Interval* ivp = order[i];
                if (cand_cache.find(ivp->name) != cand_cache.end()) return;
                auto cands = enumerate_candidates_(ivp->name, ivp->lo, ivp->hi, /*require_avail=*/false);
                std::sort(cands.begin(), cands.end(), [](const Candidate& a, const Candidate& b){
                    if (a.coverage() != b.coverage()) return a.coverage() > b.coverage();
                    return a.add > b.add;
                });
                temp_cands[i] = std::move(cands);
            });
            for (size_t i = 0; i < (size_t)n_ord; ++i) {
                if (!temp_cands[i].empty() || cand_cache.find(order[i]->name) == cand_cache.end()) {
                    cand_cache.emplace(order[i]->name, std::move(temp_cands[i]));
                }
            }
        }

        // Seed: dirty queue with candidates having potential coverage gain.
        std::vector<int> dirty;
        dirty.reserve(n_ord);
        for (int i = 0; i < n_ord; ++i) {
            const Interval* ivp = order[i];
            int cur_cov = cov_of(ivp->name);
            if (cur_cov >= ivp->hi - ivp->lo) continue;
            const auto& cands = get_cands(*ivp);
            if (cands.empty() || cands.front().coverage() <= cur_cov) continue;
            dirty.push_back(i);
        }
        std::vector<char> in_next(n_ord, 0);
        std::vector<int> next_dirty;
        next_dirty.reserve((size_t)n_ord);
        std::vector<int> pin_blocked;
        pin_blocked.reserve(64);
        std::vector<char> in_pin_blocked(n_ord, 0);

        auto mark_next = [&](int oi) {
            if (oi < 0 || oi >= n_ord || in_next[oi]) return;
            const Interval* iv = order[oi];
            int cv = cov_of(iv->name);
            if (cv >= iv->hi - iv->lo) return;
            const auto& cands = get_cands(*iv);
            if (cands.empty() || cands.front().coverage() <= cv) return;
            in_next[oi] = 1;
            next_dirty.push_back(oi);
            ++phase2_dirty_marks;
        };

        int guard = 0;
        int stall = 0;
        int last_gen_kept_drop = 0;
        const char* stop_reason = "fixed-point";
        while (!dirty.empty()) {
            if (guard >= max_iters) {
                stop_reason = "max-iters";
                break;
            }
            ++guard;
            next_dirty.clear();
            std::fill(in_next.begin(), in_next.end(), 0);
            for (int oi : pin_blocked) in_pin_blocked[oi] = 0;
            pin_blocked.clear();
            bool improved_gen = false;
            int gen_kept_drop = 0;
            std::vector<int> still_open; // best-cached > cur after this visit
            still_open.reserve(64);

            // order[] is size-desc; sorting dirty indices restores size order.
            std::sort(dirty.begin(), dirty.end());

            bool hit_time_limit = false;

            struct Proposed {
                int oi = -1;
                const Candidate* cand = nullptr;
            };
            std::vector<Proposed> proposed(dirty.size());
            std::atomic<size_t> eval_counter{0};

            gcsa_parallel_for(dirty.size(), [&](size_t idx) {
                size_t cnt = eval_counter.fetch_add(1, std::memory_order_relaxed);
                if (time_budgeted && (cnt & 1023) == 0) {
                    double elapsed = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
                    if (elapsed >= time_limit_ms) {
                        hit_time_limit = true;
                        return;
                    }
                }
                if (hit_time_limit) return;

                int oi = dirty[idx];
                const Interval* ivp = order[oi];
                const int isize = ivp->hi - ivp->lo;
                int cur_cov = cov_of(ivp->name);
                if (cur_cov >= isize) return;

                const auto& cands = get_cands(*ivp);
                if (cands.empty() || cands.front().coverage() <= cur_cov) return;

                for (const auto& c : cands) {
                    if (c.coverage() < kMinCoverage) continue;
                    if (c.coverage() <= cur_cov) break;
                    if (!run_kept_(c.src_lo, c.src_hi)) continue;
                    bool ok = true;
                    for (int32_t r : c.covered) {
                        if (removed_[r]) { ok = false; break; }
                    }
                    if (!ok) continue;

                    // Quick retarget check for pinned covered ranks
                    bool pin_ok = true;
                    for (int32_t r : c.covered) {
                        if (pin_count_[r] > 0) {
                            if (pin_count_[r] == 1) {
                                uint64_t owner = pin_owner_[r];
                                auto it = accepted.find(owner);
                                if (it == accepted.end()) { pin_ok = false; break; }
                                int32_t nlo, nhi;
                                if (!can_retarget_(it->second, c, nlo, nhi)) { pin_ok = false; break; }
                            } else {
                                auto mit = multi_pin_owners_.find(r);
                                if (mit != multi_pin_owners_.end()) {
                                    for (uint64_t owner : mit->second) {
                                        auto it = accepted.find(owner);
                                        if (it == accepted.end()) { pin_ok = false; break; }
                                        int32_t nlo, nhi;
                                        if (!can_retarget_(it->second, c, nlo, nhi)) { pin_ok = false; break; }
                                    }
                                } else {
                                    pin_ok = false; break;
                                }
                            }
                        }
                    }
                    if (!pin_ok) continue;

                    proposed[idx] = {oi, &c};
                    break;
                }
            });

            if (hit_time_limit) {
                stop_reason = "phase1-time-budget";
                break;
            }

            for (size_t idx = 0; idx < dirty.size(); ++idx) {
                if (time_budgeted && (phase2_tries & 63) == 0 && phase2_tries > 0) {
                    double elapsed = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
                    if (elapsed >= time_limit_ms) {
                        stop_reason = "phase1-time-budget";
                        hit_time_limit = true;
                        break;
                    }
                }
                int oi = dirty[idx];
                const Interval* ivp = order[oi];
                const int isize = ivp->hi - ivp->lo;
                int cur_cov = cov_of(ivp->name);
                if (cur_cov >= isize) { ++phase2_skip_sat; continue; }

                const auto& cands = get_cands(*ivp);
                if (cands.empty() || cands.front().coverage() <= cur_cov) {
                    ++phase2_skip_best;
                    continue;
                }
                const int best_cov = cands.front().coverage();

                int got_cov = cur_cov;
                const Candidate* prop_cand = proposed[idx].cand;

                auto try_cand = [&](const Candidate& c) -> bool {
                    if (c.coverage() < kMinCoverage) return false;
                    if (c.coverage() <= cur_cov) return false;
                    if (!run_kept_(c.src_lo, c.src_hi)) return false;
                    for (int32_t r : c.covered) if (removed_[r]) return false;

                    size_t before = kept_count_;
                    ++phase2_tries;
                    auto tt0 = Clock::now();
                    bool accepted_ok = try_accept_with_retarget_(c, accepted, intervals,
                                                                /*undo_out=*/nullptr);
                    if (timing) ms_try += std::chrono::duration<double, std::milli>(Clock::now() - tt0).count();
                    if (last_pin_path_ == 0) ++phase2_no_pin;
                    else if (last_pin_path_ == 1) ++phase2_fast_pin;
                    else ++phase2_multi_pin;
                    if (!accepted_ok) {
                        ++phase2_fail;
                        if (last_pin_path_ != 0 && !in_pin_blocked[oi]) {
                            in_pin_blocked[oi] = 1;
                            pin_blocked.push_back(oi);
                        }
                        return false;
                    }
                    size_t after = kept_count_;
                    if (after < before) gen_kept_drop += (int)(before - after);
                    if (trace) {
                        if (after < before)
                            gcsa_log(
                                "II.%d IMPROVE I_%s src=I_%s[%d,%d) add=%d cov=%d  kept %zu->%zu\n",
                                guard, name_to_string(G_.shape, ivp->name).c_str(),
                                name_to_string(G_.shape, name_of_rank_(c.src_lo)).c_str(),
                                c.src_lo, c.src_hi, c.add, c.coverage(), before, after);
                        else
                            gcsa_log(
                                "II.%d REPLACE I_%s cov %d->%d (kept unchanged %zu)\n",
                                guard, name_to_string(G_.shape, ivp->name).c_str(),
                                cur_cov, c.coverage(), after);
                    }
                    improved_gen = true;
                    ++phase2_improve;
                    got_cov = c.coverage();
                    return true;
                };

                bool applied = false;
                if (prop_cand) {
                    applied = try_cand(*prop_cand);
                }
                if (!applied) {
                    for (const auto& c : cands) {
                        if (&c == prop_cand) continue;
                        if (c.coverage() <= cur_cov) break;
                        if (try_cand(c)) { applied = true; break; }
                    }
                }

                // Still room vs cached best (pin-blocked / source not kept yet /
                // accepted a suboptimal cand) — retry next generation if anyone
                // improved (pins / kept sets may have moved).
                if (best_cov > got_cov) still_open.push_back(oi);
            }

            if (hit_time_limit) break;
            if (time_budgeted) {
                double elapsed = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
                if (elapsed >= time_limit_ms) {
                    stop_reason = "phase1-time-budget";
                    break;
                }
            }

            last_gen_kept_drop = gen_kept_drop;
            if (improved_gen) {
                for (int oi : pin_blocked) mark_next(oi);
                for (int oi : still_open) mark_next(oi);
            }
            dirty.swap(next_dirty);

            // Adaptive stall: generations that still churn (dirty non-empty next)
            // but save fewer than min_gain kept positions. Pure REPLACE plateaus
            // (coverage up, |C| unchanged) count toward the stall streak.
            if (adaptive && !dirty.empty()) {
                if (gen_kept_drop < min_gain) {
                    if (++stall >= stall_limit) {
                        stop_reason = "adaptive-stall";
                        break;
                    }
                } else {
                    stall = 0;
                }
            } else {
                stall = 0;
            }
        }
        auto t1 = Clock::now();
        gcsa_log("[%s] Phase II done\n", label);
        if (trace) {
            gcsa_log("after Phase II: kept=%zu accepted=%zu\n",
                     kept_count_, accepted.size());
        }
        if (timing) {
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            gcsa_log(
                "[timing] %s Phase II: %.1fms  gen %d/%d (%s, stop=%s) "
                "enums=%zu cache_hits=%zu "
                "skip_sat=%zu skip_best=%zu tries=%zu "
                "fail=%zu improve=%zu dirty_marks=%zu "
                "pin(none/fast/multi)=%zu/%zu/%zu "
                "enum=%.1fms try=%.1fms",
                label, ms, guard, max_iters, iters_src, stop_reason,
                phase2_enum, phase2_cache_hits,
                phase2_skip_sat, phase2_skip_best, phase2_tries,
                phase2_fail, phase2_improve, phase2_dirty_marks,
                phase2_no_pin, phase2_fast_pin, phase2_multi_pin,
                ms_enum, ms_try);
            if (std::strcmp(stop_reason, "adaptive-stall") == 0) {
                gcsa_log(
                    "  (kept_drop=%d < min_gain=%d for stall=%d/%d)",
                    last_gen_kept_drop, min_gain, stall, stall_limit);
            } else if (std::strcmp(stop_reason, "max-iters") == 0) {
                gcsa_log("  (generation budget exhausted)");
            } else if (std::strcmp(stop_reason, "phase1-time-budget") == 0) {
                gcsa_log("  (hit Phase I time budget %.1fms)", time_limit_ms);
            }
            if (phase1_ms > 0.0) {
                gcsa_log("  phase2/phase1=%.2f", ms / phase1_ms);
            }
            gcsa_log("\n");
        }
    }

    // Phase II, generalized to every remaining algorithm: the dirty-set
    // greedy unpin/retarget (run_phase2_) followed by the exact bounded-
    // cluster LNS (run_phase2_local_). This used to be tree-dp3's pipeline
    // alone; every algorithm now gets both passes, since the LNS only ever
    // lowers |C| and composes with the retarget loop rather than replacing
    // it (retargeting reaches links -- composite adds -- that the static
    // candidate cache does not hold, which is exactly where the LNS alone
    // gives ground on long repeats). GCSA_LNS_ONLY=1 skips the greedy pass
    // and measures the LNS alone; GCSA_DISABLE_PHASE2=1 (checked inside
    // run_phase2_) skips the greedy pass the same way, but the LNS pass
    // still runs regardless -- it isn't gated by that flag.
    void run_phase2_and_lns_(const std::vector<Interval>& intervals,
                             std::unordered_map<uint64_t, Candidate>& accepted,
                             const char* label, bool trace, bool timing,
                             std::unordered_map<uint64_t, std::vector<Candidate>>* cand_cache,
                             double phase1_ms) {
        if (gcsa_env_int("GCSA_LNS_ONLY", 0) == 0)
            run_phase2_(intervals, accepted, label, trace, timing, cand_cache, phase1_ms);
        run_phase2_local_(intervals, accepted, label, timing, cand_cache);
    }

    // ---- algorithms -------------------------------------------------------
    void compress_greedy_size_(const std::vector<Interval>& intervals,
                          std::unordered_map<uint64_t, Candidate>& accepted) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        const bool trace = (std::getenv("GCSA_TRACE_GREEDY_SIZE") != nullptr);
        auto t0 = Clock::now();
        // m0 = |C| before this algorithm's own Phase I; pct_kept expresses
        // later |C| values as a percentage of that baseline.
        const size_t m0 = kept_count_;
        auto pct_kept = [&](size_t kept) {
            return 100.0 * (double)kept / (double)std::max<size_t>(1, m0);
        };
        gcsa_log("[greedy-size] size-order accept...\n");
        std::vector<const Interval*> order;
        for (const auto& iv : intervals) if (iv.hi - iv.lo > 1) order.push_back(&iv);
        std::sort(order.begin(), order.end(),
                  [](const Interval* a, const Interval* b){ return (a->hi-a->lo) > (b->hi-b->lo); });
        int step = 0;
        if (trace) {
            gcsa_log("greedy order:");
            for (auto* ivp : order)
                gcsa_log(" %s(size=%d)",
                         name_to_string(G_.shape, ivp->name).c_str(), ivp->hi - ivp->lo);
            gcsa_log("\n");
        }
        // Static candidate lists for all |I|>1, reused by Phase II + LNS
        // below (the same pattern every other algorithm follows). Phase I
        // here only ever needs each interval's single best candidate, so
        // this cache is left to be filled lazily by whichever of Phase II /
        // LNS touches a name first.
        std::unordered_map<uint64_t, std::vector<Candidate>> cand_cache;
        for (const Interval* ivp : order) {
            Candidate c = best_candidate_(ivp->name, ivp->lo, ivp->hi, /*avail=*/true);
            ++step;
            if (trace) {
                gcsa_log("step %d: I_%s [%d,%d) size=%d",
                         step, name_to_string(G_.shape, ivp->name).c_str(),
                         ivp->lo, ivp->hi, ivp->hi - ivp->lo);
                if (c.coverage() < kMinCoverage)
                    gcsa_log(" -> SKIP\n");
                else
                    gcsa_log(
                        " -> ACCEPT src=I_%s[%d,%d) add=%d cov=%d covered={%s}\n",
                        name_to_string(G_.shape, name_of_rank_(c.src_lo)).c_str(),
                        c.src_lo, c.src_hi, c.add, c.coverage(),
                        [&]{
                            std::string s;
                            for (size_t i = 0; i < c.covered.size(); ++i) {
                                if (i) s += ',';
                                s += std::to_string(c.covered[i]);
                            }
                            return s;
                        }().c_str());
            }
            if (c.coverage() < kMinCoverage) continue;
            accept_(c, accepted);
        }
        gcsa_log("[greedy-size] size-order accept done\n");
        const size_t after_greedy = kept_count_;
        gcsa_log("[greedy-size] size-order accept: |C| %zu -> %zu  (%.1f%% of original kept, -%zu)\n",
                 m0, after_greedy, pct_kept(after_greedy), m0 - after_greedy);

        auto t_p2_0 = Clock::now();
        double phase1_ms = std::chrono::duration<double, std::milli>(t_p2_0 - t0).count();
        run_phase2_and_lns_(intervals, accepted, "greedy-size", trace, timing, &cand_cache, phase1_ms);
        auto t_p2_1 = Clock::now();
        const size_t after_phase2 = kept_count_;
        gcsa_log("[greedy-size] phase II:  |C| %zu -> %zu  (%.1f%% of original kept, -%zu)\n",
                 after_greedy, after_phase2, pct_kept(after_phase2),
                 after_greedy - after_phase2);
        gcsa_log(
            "[greedy-size] success: original |C|=%zu -- size-order accept=%.1f%%, "
            "+phase II=%.1f%% (of original kept)\n",
            m0, pct_kept(after_greedy), pct_kept(after_phase2));
        if (timing) {
            auto ms = [](Clock::time_point a, Clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count();
            };
            gcsa_log("[timing] greedy-size: accept=%.1fms phase2=%.1fms total=%.1fms  "
                     "(#I>1=%zu accepted=%zu kept=%zu)\n",
                     phase1_ms, ms(t_p2_0, t_p2_1), ms(t0, Clock::now()),
                     order.size(), accepted.size(), kept_count_);
        }
    }

    // ------------------------------------------------------------------
    // GreedyDegree building blocks
    // ------------------------------------------------------------------
    // A richer, degree-aware greedy MWIS -- see the class comment above
    // (search "GreedyDegree") for the design rationale and the session note
    // it came out of. Three pieces:
    //   1. flatten every word's cached (up to kCandCacheDefaultCap) candidate
    //      list into one flat pool -- every candidate is its own node, not
    //      just each word's single best one.
    //   2. build the real cross-word conflict edges among pool entries.
    //   3. lazily greedy-pop by score = coverage - sum-of-best-conflicting-
    //      candidate-per-distinct-word, accepting through the same accept_()/
    //      cand_available_() machinery every other algorithm uses, so
    //      correctness never depends on this function getting the score
    //      heuristic right -- only compression quality does.

    // Score of pool[i]: its own coverage minus, for every *distinct* other
    // word still contesting it (i.e. not yet in `resolved`), that word's
    // best *currently conflicting* candidate -- not that word's own top
    // pick. Deduping by word matters: a word can have several cached
    // candidates all conflicting with pool[i] (it was only ever going to
    // use one of them), so summing all of them would overstate what
    // accepting pool[i] actually costs. Using each word's *global* best
    // instead of the *locally conflicting* best would collapse this into
    // exactly the pseudoforest-dp preference graph -- redundant with what
    // the exact DP there already prices for free (see the session note).
    long long greedy_degree_score_(int i, const std::vector<Candidate>& pool,
                                   const std::vector<std::vector<int>>& adjacency,
                                   const std::unordered_set<uint64_t>& resolved) const {
        std::unordered_map<uint64_t, int> best_per_word;
        for (int j : adjacency[i]) {
            uint64_t w = pool[j].name;
            if (resolved.count(w)) continue;  // that word's fate is already sealed elsewhere
            int cov = pool[j].coverage();
            auto it = best_per_word.find(w);
            if (it == best_per_word.end() || it->second < cov) best_per_word[w] = cov;
        }
        long long penalty = 0;
        for (auto& kv : best_per_word) penalty += kv.second;
        return (long long)pool[i].coverage() - penalty;
    }

    void compress_greedy_degree_(const std::vector<Interval>& intervals,
                                 std::unordered_map<uint64_t, Candidate>& accepted) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        const bool trace = (std::getenv("GCSA_TRACE_GREEDY_DEGREE") != nullptr);
        const char* label = "greedy-degree";
        auto t_all = Clock::now();
        const size_t m0 = kept_count_;
        auto pct_kept = [&](size_t kept) {
            return 100.0 * (double)kept / (double)std::max<size_t>(1, m0);
        };

        // ---- Step 1: every word's cached candidate list (reused by Phase II) ----
        gcsa_log("[%s] candidate enumeration...\n", label);
        auto t0 = Clock::now();
        std::unordered_map<uint64_t, std::vector<Candidate>> cand_cache;
        {
            std::vector<std::vector<Candidate>> temp(intervals.size());
            gcsa_parallel_for(intervals.size(), [&](size_t i) {
                const auto& iv = intervals[i];
                // |I|<=2 can never reach kMinCoverage=3 -- same reasoning as
                // every other algorithm's identical filter.
                if (iv.hi - iv.lo <= 2) return;
                temp[i] = enumerate_candidates_(iv.name, iv.lo, iv.hi, /*avail=*/false);
            });
            for (size_t i = 0; i < intervals.size(); ++i)
                if (!temp[i].empty()) cand_cache.emplace(intervals[i].name, std::move(temp[i]));
        }
        auto t1 = Clock::now();

        // ---- Step 2: flatten into one pool, indexed per word ----
        std::vector<Candidate> pool;
        std::unordered_map<uint64_t, std::vector<int>> word_to_pool;
        for (auto& kv : cand_cache) {
            auto& idxs = word_to_pool[kv.first];
            idxs.reserve(kv.second.size());
            for (const auto& c : kv.second) {
                idxs.push_back((int)pool.size());
                pool.push_back(c);  // copy -- cand_cache's originals stay put for Phase II
            }
        }
        const int N = (int)pool.size();
        gcsa_log("[%s] candidate pool: %d candidates across %zu words\n",
                 label, N, word_to_pool.size());

        // ---- Step 3: cross-word conflict edges only (same-word candidates
        // are mutually exclusive by construction -- see `resolved` below --
        // and materializing that as edges would be both redundant and, at
        // up to kCandCacheDefaultCap^2/2 pairs per word, needlessly large) ----
        auto t2 = Clock::now();
        std::vector<std::vector<int>> adjacency(N);
        for (int i = 0; i < N; ++i) {
            const Candidate& ci = pool[i];
            uint64_t src_word = name_of_rank_(ci.src_lo);
            if (src_word == ci.name) continue;  // self-reference: same-word, not a real edge
            auto it = word_to_pool.find(src_word);
            if (it == word_to_pool.end()) continue;  // source word has no candidate of its own
            for (int j : it->second) {
                const Candidate& cj = pool[j];
                // Does ci's source intersect cj's covered set? Same overlap
                // test as PseudoforestDp's ConflictGraph::dep (covered is
                // sorted ascending -- order-preserving reconstruction --
                // same reasoning as there), just run against every cached
                // candidate of the source word instead of only its single
                // preferred one.
                auto lo = std::lower_bound(cj.covered.begin(), cj.covered.end(), ci.src_lo);
                if (lo != cj.covered.end() && *lo < ci.src_hi) {
                    adjacency[i].push_back(j);
                    adjacency[j].push_back(i);
                }
            }
        }
        for (auto& adj : adjacency) {
            std::sort(adj.begin(), adj.end());
            adj.erase(std::unique(adj.begin(), adj.end()), adj.end());
        }
        auto t3 = Clock::now();
        size_t total_edges = 0;
        for (auto& adj : adjacency) total_edges += adj.size();
        gcsa_log("[%s] conflict graph: %d nodes, %zu edge-endpoints\n", label, N, total_edges);

        // ---- Step 4: seed the heap ----
        std::unordered_set<uint64_t> resolved;  // words already decided (accepted or exhausted)
        using HeapEntry = std::pair<long long, int>;  // (score, pool index)
        std::priority_queue<HeapEntry> heap;
        for (int i = 0; i < N; ++i)
            heap.push({greedy_degree_score_(i, pool, adjacency, resolved), i});
        auto t4 = Clock::now();

        // ---- Step 5: lazy greedy pop ----
        // Monotonicity: pool[i]'s true score only ever *increases* over time
        // (as contesting words get resolved and drop out of its penalty
        // sum), so every entry still sitting in the heap is an upper bound
        // on its own true current score. When the popped max's freshly
        // recomputed score still equals what it was pushed with, nothing
        // else in the heap can possibly beat it right now -- safe to accept
        // immediately. Standard lazy-priority-queue pattern (decrease-key
        // via re-push instead of an in-place update).
        size_t n_accepted = 0, n_dead = 0, n_stale_repush = 0;
        while (!heap.empty()) {
            long long score = heap.top().first;
            int i = heap.top().second;
            heap.pop();
            if (resolved.count(pool[i].name)) continue;  // word already decided by a different candidate
            long long fresh = greedy_degree_score_(i, pool, adjacency, resolved);
            if (fresh != score) { heap.push({fresh, i}); ++n_stale_repush; continue; }
            if (!cand_available_(pool[i])) { ++n_dead; continue; }  // source/covered consumed meanwhile
            Candidate c = pool[i];
            if (trace) {
                gcsa_log("  ACCEPT %s <- %s add=%d cov=%d score=%lld\n",
                    name_to_string(G_.shape, c.name).c_str(),
                    name_to_string(G_.shape, name_of_rank_(c.src_lo)).c_str(),
                    c.add, c.coverage(), score);
            }
            accept_(c, accepted);
            resolved.insert(pool[i].name);
            ++n_accepted;
        }
        auto t5 = Clock::now();
        gcsa_log("[%s] greedy pop done: accepted=%zu dead=%zu stale-repush=%zu\n",
                 label, n_accepted, n_dead, n_stale_repush);
        const size_t after_greedy = kept_count_;
        gcsa_log("[%s] greedy: |C| %zu -> %zu  (%.1f%% of original kept, -%zu)\n",
                 label, m0, after_greedy, pct_kept(after_greedy), m0 - after_greedy);

        // No separate leftover pass: every cached candidate of every word
        // already got a turn in the same heap (a word only ever leaves the
        // loop unresolved if none of its up-to-cap candidates ever cleared
        // both a fresh-score check and cand_available_, i.e. it genuinely
        // has nothing left, exactly the leftover-pass exit condition too --
        // see run_pseudoforest_leftover_'s comment for that same floor).
        auto t_p2_0 = Clock::now();
        double phase1_ms = std::chrono::duration<double, std::milli>(t5 - t0).count();
        run_phase2_and_lns_(intervals, accepted, label, trace, timing, &cand_cache, phase1_ms);
        auto t_p2_1 = Clock::now();
        const size_t after_phase2 = kept_count_;
        gcsa_log("[%s] phase II:  |C| %zu -> %zu  (%.1f%% of original kept, -%zu)\n",
                 label, after_greedy, after_phase2, pct_kept(after_phase2),
                 after_greedy - after_phase2);
        gcsa_log(
            "[%s] success: original |C|=%zu -- greedy=%.1f%%, +phase II=%.1f%% (of original kept)\n",
            label, m0, pct_kept(after_greedy), pct_kept(after_phase2));

        if (timing) {
            auto ms = [](Clock::time_point a, Clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count();
            };
            const double phase2_total = ms(t_p2_0, t_p2_1);
            gcsa_log(
                "[timing] %s enum=%.1fms flatten=%.1fms graph=%.1fms heap_init=%.1fms "
                "pop=%.1fms phase2=%.1fms total=%.1fms  pool=%d edges=%zu "
                "accepted=%zu kept=%zu\n",
                label, ms(t0, t1), ms(t1, t2), ms(t2, t3), ms(t3, t4), ms(t4, t5),
                phase2_total, ms(t_all, Clock::now()), N, total_edges,
                accepted.size(), kept_count_);
        }
    }

    void compress_dep_order_(const std::vector<Interval>& intervals,
                             std::unordered_map<uint64_t, Candidate>& accepted) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        auto t_all = Clock::now();
        // m0 = |C| before this algorithm's own Phase I; pct_kept
        // expresses later |C| values as a percentage of that baseline.
        const size_t m0 = kept_count_;
        auto pct_kept = [&](size_t kept) {
            return 100.0 * (double)kept / (double)std::max<size_t>(1, m0);
        };
        // Preferred candidate per interval (ignore availability — global best).
        struct Pref { const Interval* iv; Candidate cand; uint64_t src_name; };
        std::vector<Pref> prefs;
        std::unordered_map<uint64_t, int> id_of;   // k-mer name -> node id
        for (size_t i = 0; i < intervals.size(); ++i)
            id_of[intervals[i].name] = (int)i;
        const int N = (int)intervals.size();

        // Preference relation: only intervals with |I_c| > 2 participate
        // (cardinality-2 intervals are ignored for prefs / DAG / Phase I).
        // Static cand lists for all |I|>1 are reused by Phase II.
        std::unordered_map<uint64_t, std::vector<Candidate>> cand_cache;
        gcsa_log("[dep-order] preference / DAG build...\n");
        auto t0 = Clock::now();
        std::vector<std::vector<Candidate>> temp_cands(intervals.size());
        std::vector<Pref> temp_prefs(intervals.size());
        std::vector<char> has_pref(intervals.size(), 0);

        gcsa_parallel_for(intervals.size(), [&](size_t i) {
            const auto& iv = intervals[i];
            if (iv.hi - iv.lo <= 1) return;
            auto cands = enumerate_candidates_(iv.name, iv.lo, iv.hi, /*avail=*/false);
            if (iv.hi - iv.lo > 2) {
                Candidate c = pick_best_(cands);
                if (c.coverage() >= kMinCoverage) {
                    uint64_t src_name = name_of_rank_(c.src_lo);
                    temp_prefs[i] = Pref{&iv, std::move(c), src_name};
                    has_pref[i] = 1;
                }
            }
            std::sort(cands.begin(), cands.end(), [](const Candidate& a, const Candidate& b){
                if (a.coverage() != b.coverage()) return a.coverage() > b.coverage();
                return a.add > b.add;
            });
            temp_cands[i] = std::move(cands);
        });

        for (size_t i = 0; i < intervals.size(); ++i) {
            if (has_pref[i]) prefs.push_back(std::move(temp_prefs[i]));
            if (intervals[i].hi - intervals[i].lo > 1) {
                cand_cache.emplace(intervals[i].name, std::move(temp_cands[i]));
            }
        }
        auto t1 = Clock::now();

        std::unordered_map<uint64_t, Pref*> pref_map;
        pref_map.reserve(prefs.size());
        for (auto& p : prefs) pref_map[p.iv->name] = &p;

        // DAG edge: src_name -> target_name (target depends on source).
        // Only introduce an edge when the preferred candidate has cov >= kMinCoverage
        // (same floor as Phase I/II acceptance; weaker prefs are not recorded).
        std::vector<std::vector<int>> outs(N), ins(N);
        std::vector<int> indeg(N, 0);
        auto pref_of = [&](uint64_t name) -> Pref* {
            auto it = pref_map.find(name);
            return it != pref_map.end() ? it->second : nullptr;
        };
        for (auto& p : prefs) {
            if (p.cand.coverage() < kMinCoverage) continue;   // no preference-DAG edge below floor
            int t = id_of[p.iv->name];
            auto it = id_of.find(p.src_name);
            if (it == id_of.end()) continue;
            int s = it->second;
            if (s == t) continue;                  // self-loop (shouldn't after disjoint)
            outs[s].push_back(t);
            ins[t].push_back(s);
            ++indeg[t];
        }

        const bool trace = (std::getenv("GCSA_TRACE_DEP") != nullptr);
        if (trace) {
            gcsa_log("=== preferred ===\n");
            for (auto& p : prefs)
                gcsa_log("  %s -> prefers %s add=%d cov=%d src[%d,%d)\n",
                    name_to_string(G_.shape, p.iv->name).c_str(),
                    name_to_string(G_.shape, p.src_name).c_str(),
                    p.cand.add, p.cand.coverage(), p.cand.src_lo, p.cand.src_hi);
            gcsa_log("=== DAG edges (src -> target, cov>=%d) ===\n", kMinCoverage);
            for (auto& p : prefs) {
                if (p.cand.coverage() < kMinCoverage) continue;
                if (p.iv->name == p.src_name) continue;
                gcsa_log("  %s -> %s\n",
                    name_to_string(G_.shape, p.src_name).c_str(),
                    name_to_string(G_.shape, p.iv->name).c_str());
            }
        }

        // Kahn topological order; sources before dependents.
        // We *accept* in reverse topo (sinks / dependents first), then improve.
        std::queue<int> q;
        for (int i = 0; i < N; ++i) if (indeg[i] == 0) q.push(i);
        std::vector<int> topo;
        topo.reserve(N);
        std::vector<int> indeg2 = indeg;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            topo.push_back(u);
            for (int v : outs[u]) if (--indeg2[v] == 0) q.push(v);
        }
        // Cycles (if any): append remaining nodes by interval size desc.
        if ((int)topo.size() < N) {
            std::vector<char> seen(N, 0);
            for (int u : topo) seen[u] = 1;
            std::vector<int> rest;
            for (int i = 0; i < N; ++i) if (!seen[i]) rest.push_back(i);
            std::sort(rest.begin(), rest.end(), [&](int a, int b){
                return (intervals[a].hi-intervals[a].lo) > (intervals[b].hi-intervals[b].lo);
            });
            for (int u : rest) topo.push_back(u);
        }
        gcsa_log("[dep-order] preference / DAG done\n");

        if (trace) {
            gcsa_log("=== Phase I accept order (reverse topo, sinks first) ===\n  ");
            for (int k = N - 1; k >= 0; --k) {
                const auto& iv = intervals[topo[k]];
                if (iv.hi - iv.lo <= 2) continue;
                if (!pref_of(iv.name)) continue;
                gcsa_log("%s ", name_to_string(G_.shape, iv.name).c_str());
            }
            gcsa_log("\n");
        }

        // Pass 1: accept preferred candidates in reverse topo (sinks first),
        // with availability. This pins sources that dependents need.
        gcsa_log("[dep-order] Phase I: reverse-topo accept...\n");
        auto t2 = Clock::now();
        int step = 0;
        for (int k = N - 1; k >= 0; --k) {
            uint64_t name = intervals[topo[k]].name;
            Pref* p = pref_of(name);
            if (!p) continue;
            // Recompute under current availability (source may already be gone).
            Candidate c = best_candidate_(p->iv->name, p->iv->lo, p->iv->hi, /*avail=*/true);
            ++step;
            if (trace) {
                gcsa_log("I.%d I_%s", step, name_to_string(G_.shape, name).c_str());
                if (c.coverage() < kMinCoverage) gcsa_log(" -> SKIP\n");
                else gcsa_log(" -> ACCEPT src=I_%s[%d,%d) add=%d cov=%d\n",
                    name_to_string(G_.shape, name_of_rank_(c.src_lo)).c_str(),
                    c.src_lo, c.src_hi, c.add, c.coverage());
            }
            if (c.coverage() < kMinCoverage) continue;
            accept_(c, accepted);
        }
        auto t3 = Clock::now();
        gcsa_log("[dep-order] Phase I done\n");
        if (trace) {
            size_t kept = 0; for (auto x : removed_) if (!x) ++kept;
            gcsa_log("after Phase I: kept=%zu accepted=%zu\n", kept, accepted.size());
        }
        const size_t after_phase1 = kept_count_;
        gcsa_log("[dep-order] Phase I:   |C| %zu -> %zu  (%.1f%% of original kept, -%zu)\n",
                 m0, after_phase1, pct_kept(after_phase1), m0 - after_phase1);

        double phase1_ms = std::chrono::duration<double, std::milli>(t3 - t0).count();
        run_phase2_and_lns_(intervals, accepted, "dep-order", trace, timing, &cand_cache, phase1_ms);
        auto t4 = Clock::now();
        const size_t after_phase2 = kept_count_;
        gcsa_log("[dep-order] phase II:  |C| %zu -> %zu  (%.1f%% of original kept, -%zu)\n",
                 after_phase1, after_phase2, pct_kept(after_phase2),
                 after_phase1 - after_phase2);
        gcsa_log(
            "[dep-order] success: original |C|=%zu -- Phase I=%.1f%%, "
            "+phase II=%.1f%% (of original kept)\n",
            m0, pct_kept(after_phase1), pct_kept(after_phase2));
        if (timing) {
            auto ms = [](Clock::time_point a, Clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count();
            };
            gcsa_log(
                "[timing] dep-order prefs=%.1fms phaseI=%.1fms phaseII=%.1fms total=%.1fms\n"
                "         N=%d prefs=%zu\n",
                ms(t0, t1), ms(t2, t3), ms(t3, t4), ms(t_all, t4),
                N, prefs.size());
        }
    }

    // ------------------------------------------------------------------
    // PseudoforestDp building blocks
    // ------------------------------------------------------------------
    // The three pieces below implement one "extract preference, solve
    // exactly" round (exact_pseudoforest_dp.pdf, "Sparsifying the graph" /
    // "MWIS on pseudoforests"):
    //   1. build_pseudoforest_graph_    -- extract the preference pseudoforest
    //   2. solve_pseudoforest_dp_       -- solve it exactly (KEEP vs COMPRESS)
    //   3. accept_pseudoforest_chosen_  -- materialize the DP's choices
    // compress_pseudoforest_dp_ further below just calls these in a loop,
    // each round over whatever intervals earlier rounds left unresolved --
    // see "IDEA: Iterative DP" in the note.

    // A single round's preference pseudoforest: every interval with
    // |I_c|>2 that isn't already resolved gets, at most, one candidate it
    // would use to compress (`cand` below) -- and, separately, at most one
    // REAL dependency (`dep` below) on whichever interval that candidate's
    // source actually collides with. These are deliberately not the same
    // map: every node with a candidate gets an entry in `cand`, but only a
    // node whose candidate's source range genuinely overlaps its source's
    // own chosen candidate's covered set gets an entry in `dep` -- a node
    // can have a candidate and no dependency at all. That's exactly what
    // always happens for a self-reference (L1+L3 guarantee its source can
    // never overlap its own covered set), but it's not special-cased: any
    // node whose candidate doesn't actually conflict with anything ends up
    // in `cand` only, self-reference or not.
    //
    // Only `dep` entries are real edges: they're the only things cycle
    // detection chases and the only things registered in `children` below.
    // A candidate's source always lies entirely inside a single interval,
    // since LCP>=add+1>=2 forces every source rank to share its own first
    // symbol (see "Sparsifying the graph"), so out-degree in `dep` is <=1
    // per node -- that's what makes this a pseudoforest: every connected
    // component of real dependencies is either a tree (rooted at a node
    // with no dependency) or unicyclic (exactly one cycle, with trees
    // hanging off each cycle node) -- never anything worse. No edge is ever
    // dropped to force acyclicity: unicyclic components are
    // handled by folding the cycle instead (see solve_pseudoforest_dp_), so
    // no cov/size threshold is needed either -- every still-unresolved
    // interval participates.
    //
    // Splitting "has a candidate" from "has a real dependency" like this
    // (instead of one record with a conflicts flag riding along on it)
    // means there's no runtime "is this actually a conflict" check left to
    // do once the graph is built: presence in `dep` *is* the conflict, and
    // solve_pseudoforest_dp_ never needs to ask again.
    struct ConflictGraph {
        // Every interval in the instance, including already-resolved /
        // |I|<=2 ones -- needed because those can still be valid *sources*
        // even though they're never DP-decidable nodes this round.
        std::unordered_map<uint64_t, const Interval*> by_name;
        // This round's DP-decidable intervals are exactly the keys of
        // `cand` below -- a node only matters to the DP once it actually
        // has something to compress via, so there's no separate "nodes"
        // set: a node with no candidate can never be anyone's real
        // dependency (dep is only ever keyed from `cand`) and can never
        // itself become COMPRESS, so tracking it here would be inert.
        std::unordered_map<uint64_t, Candidate> cand;       // node -> its one chosen candidate, if any
        std::unordered_map<uint64_t, uint64_t> dep;         // node -> the node it has a REAL dependency on (conflicting)
        std::unordered_map<uint64_t, std::vector<uint64_t>> children;  // src -> dependents (cycle-internal edge excluded)
        std::vector<std::vector<uint64_t>> cycles;          // each cycle as an ordered node list
        std::unordered_map<uint64_t, uint64_t> cycle_pred;  // cyc[i] -> cyc[i-1] (the excluded edge)
    };

    // Extract one round's preference pseudoforest. `intervals` is always the
    // whole instance (interval bounds never change); `accepted` marks which
    // names earlier rounds already resolved, so they're excluded from
    // `cand` here -- they remain valid, fixed sources (see
    // solve_pseudoforest_dp_) but are never reconsidered for compression.
    // require_avail=false (round 0, nothing accepted yet) picks each node's
    // globally-best candidate, exactly like the original single-pass
    // algorithm; require_avail=true (every later round) picks the best
    // candidate that is still available given what earlier rounds
    // pinned/removed -- the "extract a second set of preferences that do
    // not conflict with the fixed solution" step from the note.
    // `cand_cache` accumulates every node's full (avail=false)
    // coverage-sorted candidate list across rounds, for reuse by the
    // leftover pass and Phase II: a name's avail=false list only depends on
    // static SA/LCP structure, so it's computed once (whichever round first
    // visits that name) and reused by every later round instead of
    // re-enumerated.
    ConflictGraph build_pseudoforest_graph_(
            const std::vector<Interval>& intervals,
            const std::unordered_map<uint64_t, Candidate>& accepted,
            bool require_avail,
            std::unordered_map<uint64_t, std::vector<Candidate>>& cand_cache,
            bool rank_aware = true,
            // Set (if non-null) right after the per-interval candidate-scoring
            // loop below, before the dependency/cycle-detection pass -- lets
            // callers split "cand" (independent per interval, parallelized via
            // GCSA_THREADS like DepOrder's preference enum) from
            // "graph" (the dependency-map + cycle-detection walk, inherently
            // sequential) for GCSA_TIMING, without this function owning a
            // `timing` flag itself.
            std::chrono::steady_clock::time_point* t_cand_end = nullptr) const {
        ConflictGraph graph;
        graph.by_name.reserve(intervals.size() * 2);
        for (const auto& iv : intervals) graph.by_name[iv.name] = &iv;

        // Independent per interval -- same shape as DepOrder's
        // parallelized preference enum, now parallelized the same way via
        // GCSA_THREADS. cand_cache may already hold iv.name's entry from an
        // earlier round (see the class comment above build_pseudoforest_
        // graph_): unlike DepOrder's enumeration, which only ever runs
        // once per build against an empty cache, this function is called once
        // per pseudoforest-dp round, so a later round can hit a prior round's
        // cache. Both the lookup and the eventual insert touch that same
        // shared map, so both are locked -- an unlocked find() racing an
        // emplace() (from this or another interval) is a data race even
        // though every interval's own key is unique, because emplace can
        // trigger a rehash that invalidates every other thread's read.
        std::mutex cand_mu;
        gcsa_parallel_for(intervals.size(), [&](size_t i) {
            const Interval& iv = intervals[i];
            // |I|<=2 can never meet kMinCoverage (coverage is a subset of the
            // interval's own rows, so coverage<=|I|; the default kMinCoverage=3
            // makes |I|<=2 structurally unsatisfiable) -- same reasoning
            // compress_tree_dp_ uses for its own |I|>2 filter. Skipping here
            // avoids a wasted enumerate_candidates_ call and cand_cache entry
            // for intervals that can never produce a usable candidate; it does
            // not change which candidates get accepted.
            if (iv.hi - iv.lo <= 2) return;
            if (accepted.count(iv.name)) return;  // already resolved: fixed, not a DP node this round

            std::vector<Candidate> fresh;
            bool cached = false;
            {
                std::lock_guard<std::mutex> lock(cand_mu);
                auto cit = cand_cache.find(iv.name);
                if (cit != cand_cache.end()) {
                    fresh = cit->second;  // copy out: can't keep a ref past the unlock
                    cached = true;
                }
            }
            if (!cached) {
                fresh = enumerate_candidates_(iv.name, iv.lo, iv.hi, /*avail=*/false);
                // Stable, not std::sort: pick_best_ (the original, unrefactored
                // preference pick) scanned enumerate_candidates_'s output in its
                // native order and kept the first candidate to strictly improve
                // on (coverage, add); an unstable sort can reorder same-
                // (coverage, add) candidates arbitrarily and silently pick a
                // different -- if equally "good" by that metric -- source run,
                // which can still cascade into a different DP outcome.
                std::stable_sort(fresh.begin(), fresh.end(), [](const Candidate& a, const Candidate& b){
                    if (a.coverage() != b.coverage()) return a.coverage() > b.coverage();
                    return a.add > b.add;
                });
            }

            Candidate best = best_from_sorted_cache_(fresh, require_avail);

            std::lock_guard<std::mutex> lock(cand_mu);
            if (best.coverage() >= 2) {
                uint64_t src = name_of_rank_(best.src_lo);
                // Self-references (src == iv.name) are allowed here: nothing
                // downstream needs to know a candidate is a self-reference
                // as such -- it's just one instance of a node with no real
                // dependency (see the dep pass below), handled by the same
                // general rule as any other dependency-free candidate.
                if (graph.by_name.count(src)) graph.cand.emplace(iv.name, std::move(best));
            }
            if (!cached) cand_cache.emplace(iv.name, std::move(fresh));
        });
        if (t_cand_end) *t_cand_end = std::chrono::steady_clock::now();

        // ---- which nodes have a REAL (rank-level conflict) dependency -----
        // v has a real dependency on src = name_of_rank_(v's candidate's
        // src_lo) only if src's own candidate would, upon compressing,
        // actually delete rows that overlap v's specific source range --
        // both v's and src's candidates are already fixed at this point, so
        // this is knowable now, without running the DP. A node with a
        // candidate but no entry in `dep` (self-references always
        // included, by L1+L3) has nothing tying it to another node: it's
        // solved as an independent root, its candidate existing only as
        // something to compress via, never as a constraint.
        for (auto& kv : graph.cand) {
            uint64_t v = kv.first;
            const Candidate& vc = kv.second;
            uint64_t src = name_of_rank_(vc.src_lo);
            if (!rank_aware) { graph.dep.emplace(v, src); continue; }  // legacy: every candidate is a real dependency
            auto sit = graph.cand.find(src);
            if (sit == graph.cand.end()) continue;  // src has no candidate of its own: can't compress, no dependency
            const std::vector<int32_t>& covered = sit->second.covered;
            auto lo = std::lower_bound(covered.begin(), covered.end(), vc.src_lo);
            if (lo != covered.end() && *lo < vc.src_hi) graph.dep.emplace(v, src);
        }

        // ---- pure graph decomposition: find every cycle up front ----------
        // (independent of the DP; identifies exactly which single incoming
        // edge per cycle node is "cycle-internal" so it can be excluded from
        // the generic children[] map used by solve_pseudoforest_dp_ below --
        // that edge is instead handled explicitly by its fold_cycle scenario
        // walk, never both, which would double count it. A cycle can only
        // involve nodes with a real dependency in `dep`; every other name --
        // including any node whose candidate turned out not to conflict --
        // has no dep entry, so a chain through it always terminates instead
        // of looping; it's just another independent root. Only nodes in
        // `cand` can ever start or continue such a chain -- dep is keyed
        // exclusively from `cand` -- so walking graph.cand's keys covers
        // every possible cycle.)
        enum class Color : uint8_t { White, Gray, Black };
        std::unordered_map<uint64_t, Color> color;
        for (auto& kv : graph.cand) color[kv.first] = Color::White;

        for (auto& kv : graph.cand) {
            uint64_t start = kv.first;
            if (color[start] != Color::White) continue;
            std::vector<uint64_t> path;
            uint64_t v = start;
            bool hit_cycle = false;
            uint64_t cycle_entry = 0;
            while (true) {
                if (color[v] == Color::Black) break;
                if (color[v] == Color::Gray) { hit_cycle = true; cycle_entry = v; break; }
                color[v] = Color::Gray;
                path.push_back(v);
                auto dit = graph.dep.find(v);
                if (dit == graph.dep.end()) break;  // no real dependency: chain ends here
                v = dit->second;
            }
            if (hit_cycle) {
                auto it = std::find(path.begin(), path.end(), cycle_entry);
                std::vector<uint64_t> cyc(it, path.end());
                int L = (int)cyc.size();
                for (int i = 0; i < L; ++i)
                    graph.cycle_pred[cyc[i]] = cyc[(i + L - 1) % L];
                graph.cycles.push_back(std::move(cyc));
            }
            for (uint64_t u : path) color[u] = Color::Black;
        }

        // children[w] = dependents of w, EXCLUDING the one cycle-internal
        // edge per cycle node (handled explicitly by fold_cycle instead). w
        // itself may or may not have a candidate of its own (already
        // accepted, |I|<=2, or -- legacy mode only -- below-threshold
        // coverage) -- see solve_pseudoforest_dp_.
        for (auto& kv : graph.dep) {
            uint64_t v = kv.first, src = kv.second;
            auto it = graph.cycle_pred.find(src);
            if (it != graph.cycle_pred.end() && it->second == v) continue;
            graph.children[src].push_back(v);
        }

        return graph;
    }

    // GCSA_TRACE_PFDP dump of one round's graph.
    void trace_pseudoforest_graph_(const ConflictGraph& graph, size_t round) const {
        std::fprintf(stderr, "=== pseudoforest-dp preference graph (round %zu) ===\n", round);
        for (auto& kv : graph.cand) {
            auto dit = graph.dep.find(kv.first);
            if (dit != graph.dep.end()) {
                std::fprintf(stderr, "  %s -> %s (cov=%d)\n",
                    name_to_string(G_.shape, kv.first).c_str(),
                    name_to_string(G_.shape, dit->second).c_str(),
                    kv.second.coverage());
            } else {
                std::fprintf(stderr, "  %s (no real dependency, cov=%d)\n",
                    name_to_string(G_.shape, kv.first).c_str(),
                    kv.second.coverage());
            }
        }
        std::fprintf(stderr, "cycles found: %zu\n", graph.cycles.size());
        for (auto& cyc : graph.cycles) {
            std::fprintf(stderr, "  [");
            for (size_t i = 0; i < cyc.size(); ++i)
                std::fprintf(stderr, "%s%s", i ? "," : "",
                             name_to_string(G_.shape, cyc[i]).c_str());
            std::fprintf(stderr, "]\n");
        }
    }

    // Exact MWIS ("KEEP vs COMPRESS" per node) on one preference pseudoforest
    // (exact_pseudoforest_dp.pdf, "MWIS on pseudoforests"): tree components
    // use the tree recurrence; each unicyclic component's cycle is folded
    // into two forced scenarios (anchor kept vs anchor compressed) and the
    // cheaper one wins. Any name with no entry in graph.cand (already
    // accepted by an earlier round, |I|<=2, or -- legacy mode only --
    // below-threshold coverage) is never DP-decidable, so the recurrence
    // naturally treats it as permanently KEEP: it's never reconsidered,
    // only ever used as a fixed, always-available source for this round's
    // real nodes. Returns the node -> candidate map to accept.
    //
    // Every entry in graph.children/graph.cycles exists only because
    // build_pseudoforest_graph_ found a real, rank-conflicting dependency
    // (see ConflictGraph::dep) -- a node whose candidate didn't actually
    // conflict with anything never got registered as anyone's dependent.
    // So any child reached via children_of(), or any node walked as part of
    // a cycle here, is unconditionally blocked (src_ok=false) whenever its
    // dependency compresses: there's no "maybe it doesn't actually
    // conflict" case left to check at solve time -- that used to be a
    // rank-aware lookup here (ranges_conflict), but once the graph itself
    // stopped registering non-conflicting candidates as dependencies at
    // all, the check became tautologically true everywhere it was still
    // invoked, so it's gone. `rank_aware` only affects how
    // build_pseudoforest_graph_ populates `dep` in the first place.
    std::unordered_map<uint64_t, Candidate> solve_pseudoforest_dp_(
            const ConflictGraph& graph, bool rank_aware, bool trace) const {
        (void)rank_aware;
        auto isize = [&](uint64_t n) -> int {
            auto it = graph.by_name.find(n);
            return it == graph.by_name.end() ? 0 : (it->second->hi - it->second->lo);
        };
        auto has_pref = [&](uint64_t n) { return graph.cand.count(n) != 0; };
        auto children_of = [&](uint64_t n) -> const std::vector<uint64_t>& {
            static const std::vector<uint64_t> kEmpty;
            auto it = graph.children.find(n);
            return it == graph.children.end() ? kEmpty : it->second;
        };

        // ---- KEEP / COMPRESS DP --------------------------------------------
        struct Cell { long long cost; bool compress; };
        std::map<std::pair<uint64_t, int>, Cell> memo;

        std::function<Cell(uint64_t, bool)> solve = [&](uint64_t v, bool src_ok) -> Cell {
            auto key = std::make_pair(v, src_ok ? 1 : 0);
            auto it = memo.find(key);
            if (it != memo.end()) return it->second;

            const auto& ch = children_of(v);
            long long cost_keep = isize(v);
            for (uint64_t w : ch) cost_keep += solve(w, /*src_ok=*/true).cost;

            long long cost_comp = std::numeric_limits<long long>::max();
            bool can_comp = false;
            if (src_ok && has_pref(v)) {
                can_comp = true;
                long long cc = isize(v) - graph.cand.at(v).coverage();
                for (uint64_t w : ch)
                    cc += solve(w, /*src_ok=*/false).cost;  // w is only here because it really conflicts
                cost_comp = cc;
            }

            Cell cell = (can_comp && cost_comp < cost_keep) ? Cell{cost_comp, true}
                                                             : Cell{cost_keep, false};
            memo[key] = cell;
            return cell;
        };

        std::unordered_map<uint64_t, Candidate> chosen;  // node -> compressed-via candidate

        std::function<void(uint64_t, bool)> apply = [&](uint64_t v, bool src_ok) {
            Cell cell = solve(v, src_ok);
            if (cell.compress) {
                chosen[v] = graph.cand.at(v);
                for (uint64_t w : children_of(v)) apply(w, /*src_ok=*/false);
            } else {
                for (uint64_t w : children_of(v)) apply(w, /*src_ok=*/true);
            }
        };

        // Fold a unicyclic component's cycle cyc[0]->cyc[1]->...->cyc[L-1]->cyc[0]
        // (cyc[i]'s real dependency is cyc[i+1 mod L]) into two forced
        // sub-cases. cyc[0]'s own compress-ability requires walking all the
        // way around and back to itself -- a circular dependency broken by
        // fixing cyc[0]'s status externally instead of deriving it.
        auto fold_cycle = [&](const std::vector<uint64_t>& cyc) {
            int L = (int)cyc.size();

            // Scenario A: cyc[0] forced KEEP (always valid on its own).
            long long costA = isize(cyc[0]);
            for (uint64_t w : children_of(cyc[0])) costA += solve(w, true).cost;
            std::vector<char> statusA(L, 1);
            {
                bool avail = true;   // cyc[0] is available
                for (int i = L - 1; i >= 1; --i) {
                    Cell cell = solve(cyc[i], /*src_ok=*/avail);
                    costA += cell.cost;
                    statusA[i] = cell.compress ? 0 : 1;
                    avail = (statusA[i] == 1);
                }
            }

            // Scenario B: cyc[0] forced COMPRESS. Walk the rest of the ring
            // from cyc[L-1] down to cyc[1] with cyc[0] unavailable as a
            // source. Every cyc[i] here -- cyc[1] included, the other end of
            // the fixed cyc[0]->cyc[1] dependency -- is only in this cycle
            // at all because it's a real, conflicting dependency, so cyc[0]
            // compressing unconditionally blocks cyc[1]; the rest of the
            // ring is gated normally by its own dependency (cyc[2], or
            // cyc[0] again when L==2) via the same `avail` chain.
            long long costB = isize(cyc[0]) - graph.cand.at(cyc[0]).coverage();
            for (uint64_t w : children_of(cyc[0])) costB += solve(w, /*src_ok=*/false).cost;
            std::vector<char> statusB(L, 1);
            statusB[0] = 0;
            {
                bool avail = false;  // cyc[0] compresses: unavailable to cyc[L-1]
                for (int i = L - 1; i >= 1; --i) {
                    bool src_ok_i = (i == 1) ? false : avail;
                    Cell cell = solve(cyc[i], src_ok_i);
                    costB += cell.cost;
                    statusB[i] = cell.compress ? 0 : 1;
                    avail = (statusB[i] == 1);
                }
            }

            bool use_a = (costA <= costB);
            const std::vector<char>& status = use_a ? statusA : statusB;
            if (trace) {
                std::fprintf(stderr,
                    "cycle(len=%d) anchor-keep=%lld anchor-compress=%lld -> %s\n",
                    L, costA, costB, use_a ? "keep anchor" : "compress anchor");
            }

            for (int i = 0; i < L; ++i) {
                uint64_t v = cyc[i];
                if (status[i] == 0) {
                    chosen[v] = graph.cand.at(v);
                    for (uint64_t w : children_of(v)) apply(w, /*src_ok=*/false);
                } else {
                    for (uint64_t w : children_of(v)) apply(w, /*src_ok=*/true);
                }
            }
        };

        // Entry points: nodes the DP starts walking down from with
        // src_ok=true. Two kinds -- a fixed source some node in this round
        // depends on but that has no candidate of its own (already accepted
        // by an earlier round, |I|<=2, or -- legacy mode only -- below-
        // threshold coverage), and a node with a candidate but no real
        // dependency (a genuine root). The first kind is permanently-KEEP
        // as far as this round is concerned, and matters only because a
        // still-unresolved node's best available candidate can perfectly
        // well source from an interval that has no candidate of its own --
        // without it, every node hanging off such a source would silently
        // never get visited.
        std::unordered_set<uint64_t> entry_points;
        for (auto& kv : graph.children)
            if (!has_pref(kv.first)) entry_points.insert(kv.first);
        // A node with a candidate but no real dependency (self-references
        // included, always dependency-free by L1+L3) was never registered
        // as anyone's dependent -- it never made it into graph.children at
        // all -- so it would otherwise never be visited by the top-down
        // apply() walk below. It must be entered explicitly, same as any
        // other root.
        for (auto& kv : graph.cand)
            if (!graph.dep.count(kv.first)) entry_points.insert(kv.first);

        for (uint64_t v : entry_points) apply(v, /*src_ok=*/true);
        for (auto& cyc : graph.cycles) fold_cycle(cyc);

        return chosen;
    }

    // Materialize one round's DP choices: accept each chosen candidate
    // against the *current* removed_/pin_count_ state. By DP construction a
    // compress choice is only ever taken under src_ok=true (source already
    // decided KEPT/fixed), so accept order across `chosen` cannot conflict
    // and no dependency ordering is required; the availability re-check
    // below is defensive insurance, matching the other algorithms.
    void accept_pseudoforest_chosen_(
            const std::unordered_map<uint64_t, Candidate>& chosen,
            const std::unordered_map<uint64_t, const Interval*>& by_name,
            std::unordered_map<uint64_t, Candidate>& accepted,
            bool trace) {
        for (auto& kv : chosen) {
            uint64_t v = kv.first;
            Candidate c = kv.second;
            bool ok = run_kept_(c.src_lo, c.src_hi);
            for (int32_t rr : c.covered)
                if (removed_[rr] || pin_count_[rr] > 0) ok = false;
            if (!ok) {
                auto it = by_name.find(v);
                if (it == by_name.end()) continue;
                const Interval* iv = it->second;
                c = best_candidate_(iv->name, iv->lo, iv->hi, /*avail=*/true);
            }
            if (c.coverage() >= 2) {
                if (trace) {
                    std::fprintf(stderr, "  ACCEPT %s <- %s add=%d cov=%d\n",
                        name_to_string(G_.shape, v).c_str(),
                        name_to_string(G_.shape, name_of_rank_(c.src_lo)).c_str(),
                        c.add, c.coverage());
                }
                accept_(c, accepted);
            }
        }
    }

    // Greedy fallback pass: for every still-unresolved interval, accept its
    // best AVAILABLE candidate outright (best_candidate_ / enumerate_
    // candidates_, not the preference graph). Unlike the pseudoforest DP,
    // this sees an interval's *entire* candidate universe -- every add,
    // every source range, self-referencing included -- not just the single
    // highest-coverage candidate build_pseudoforest_graph_ extracted as
    // "the" preference for that word. That sparsification (one candidate
    // per word, needed to keep the graph a pseudoforest at all) is what
    // this sweep exists to backstop: a word's second-best candidate can
    // easily be the one that's actually still available after the DP's
    // choices land, and no number of DP rounds over the single-preference
    // graph will ever surface it; only this unrestricted, live-availability
    // sweep does. Processed size-descending, same bias as the DP's own root
    // order. Returns how many intervals it resolved, so the caller can tell
    // whether a round made progress.
    size_t run_pseudoforest_leftover_(const std::vector<Interval>& intervals,
                                      std::unordered_map<uint64_t, Candidate>& accepted,
                                      bool trace) {
        std::vector<const Interval*> leftover;
        for (const auto& iv : intervals) {
            if (iv.hi - iv.lo <= 1) continue;
            if (accepted.count(iv.name)) continue;
            leftover.push_back(&iv);
        }
        std::sort(leftover.begin(), leftover.end(),
                  [](const Interval* a, const Interval* b) {
                      return (a->hi - a->lo) > (b->hi - b->lo);
                  });
        size_t n = 0;
        for (const Interval* ivp : leftover) {
            Candidate c = best_candidate_(ivp->name, ivp->lo, ivp->hi, /*avail=*/true);
            if (c.coverage() < 2) continue;
            if (trace) {
                std::fprintf(stderr, "  LEFTOVER-ACCEPT %s <- %s add=%d cov=%d\n",
                    name_to_string(G_.shape, ivp->name).c_str(),
                    name_to_string(G_.shape, name_of_rank_(c.src_lo)).c_str(),
                    c.add, c.coverage());
            }
            accept_(c, accepted);
            ++n;
        }
        return n;
    }

    // "Extract preference, solve exactly" DP (exact_pseudoforest_dp.pdf,
    // "Sparsifying the graph"): one preference-graph extraction over every
    // unresolved interval, one exact DP solve, then (by default) straight to
    // a greedy leftover sweep that mops up whatever the DP's
    // one-candidate-per-word sparsification couldn't see (self-links
    // included -- see run_pseudoforest_leftover_), then Phase II. This is
    // the "single-fire" mode and mirrors the original algorithm exactly.
    //
    // GCSA_PFDP_ITERATE=1 turns on the doc's "IDEA: Iterative DP": instead
    // of falling straight to leftover after the first DP round, repeat
    // extract/solve to a fixed point first. Every round after the first
    // re-extracts the preference pseudoforest with require_avail=true --
    // i.e. each still-unresolved node's *currently available* best
    // candidate, given every earlier round's accepted set -- over just the
    // residual (shrunk) instance, and solves it exactly again, same as
    // round 0. Leftover still runs exactly once, after the loop reaches a
    // fixed point (a round whose DP accepts nothing) or hits
    // GCSA_PFDP_MAX_ROUNDS -- it remains the only backstop for whatever no
    // round's one-candidate-per-word sparsification could even pose as a DP
    // node (see run_pseudoforest_leftover_). The loop needs no artificial
    // cap to terminate: a productive round accepts at least one more name,
    // which strictly shrinks the residual instance, so it's bounded by the
    // interval count regardless; GCSA_PFDP_MAX_ROUNDS is only there to bail
    // out early if that's ever desired.
    //
    // NOT a strict improvement, despite every individual round being an
    // exact MWIS solve: each round's DP commits, in one batch snapshot, to
    // exactly one candidate per still-unresolved node and can never
    // reconsider a *different* candidate for that same node within the
    // round, even when the chosen one loses to a conflict. Plain leftover
    // greedy, by contrast, re-derives "the best still-available candidate"
    // for one node at a time against always-current state (accept_'s
    // effects from every prior node already applied), so it can fall back
    // to a node's second-best candidate mid-pass in a way no DP round can
    // until an entirely new round starts. Measured on the E. coli benchmark
    // genomes (`bench/`'s pool, several shapes and max-add values): the
    // net effect through leftover (and again through Phase II) is a coin
    // flip -- typically zero, otherwise a handful of |C| positions either
    // way out of several million kept, never consistently one direction.
    // Off by default pending a wider real-data sweep; build_pseudoforest_
    // graph_'s require_avail parameter and the round-oriented framing below
    // existed for exactly this experiment.
    void compress_pseudoforest_dp_(const std::vector<Interval>& intervals,
                                   std::unordered_map<uint64_t, Candidate>& accepted) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        const bool trace = (std::getenv("GCSA_TRACE_PFDP") != nullptr);
        const bool rank_aware = [] {
            const char* e = std::getenv("GCSA_PFDP_RANK_AWARE");
            return !e || std::atoi(e) != 0;
        }();
        // Iteration is primarily selected via --algo pseudoforest-dp-iterate
        // (CompressAlgo::PseudoforestDpIterate); GCSA_PFDP_ITERATE=1 remains
        // as a legacy override that also turns iteration on for plain
        // --algo pseudoforest-dp.
        const bool iterate = (algo_ == CompressAlgo::PseudoforestDpIterate) || [] {
            const char* e = std::getenv("GCSA_PFDP_ITERATE");
            return e && std::atoi(e) != 0;
        }();
        // 0 == unbounded (run to a fixed point; see comment above for why
        // that's already guaranteed to terminate without a cap).
        const size_t max_rounds = [] () -> size_t {
            const char* e = std::getenv("GCSA_PFDP_MAX_ROUNDS");
            if (!e) return 0;
            long v = std::atol(e);
            return v > 0 ? (size_t)v : 0;
        }();
        const char* label = algo_name(algo_);
        auto t_all = Clock::now();
        const int nthreads = gcsa_num_threads();
        // m0 = |C| before this algorithm's own Phase I; pct_kept
        // expresses later |C| values as a percentage of that baseline.
        const size_t m0 = kept_count_;
        auto pct_kept = [&](size_t kept) {
            return 100.0 * (double)kept / (double)std::max<size_t>(1, m0);
        };

        std::unordered_map<uint64_t, std::vector<Candidate>> cand_cache;  // reused by Phase II

        // Wall-clock start of Phase I proper (graph build through leftover),
        // for run_phase2_'s time-budget calc. Captured once here (not per
        // round) so it
        // covers every round when iterating, the same way it covered the
        // single round before.
        auto t_phase1_start = Clock::now();
        // Per-phase ms, summed across every round, for GCSA_TIMING.
        double cand_ms = 0.0, graph_ms = 0.0, dp_ms = 0.0, accept_ms = 0.0;
        size_t total_dp_nodes = 0, total_cycles = 0, total_chosen = 0;
        size_t rounds = 0;

        for (;;) {
            // Round 0 matches the original algorithm's single pass exactly
            // (require_avail=false: each node's globally-best candidate,
            // nothing accepted yet). Every later round asks for each
            // still-unresolved node's best candidate that is still
            // available given what earlier rounds already accepted -- the
            // "extract a second set of preferences that do not conflict
            // with the fixed solution" step from the doc.
            const bool require_avail = (rounds > 0);
            gcsa_log("[%s] round %zu: preference graph build...\n", label, rounds);
            auto t_graph0 = Clock::now();
            Clock::time_point t_cand_end = t_graph0;
            ConflictGraph graph = build_pseudoforest_graph_(intervals, accepted,
                                                        require_avail, cand_cache,
                                                        rank_aware, &t_cand_end);
            auto t_graph1 = Clock::now();
            if (trace) trace_pseudoforest_graph_(graph, rounds);
            cand_ms += std::chrono::duration<double, std::milli>(t_cand_end - t_graph0).count();
            graph_ms += std::chrono::duration<double, std::milli>(t_graph1 - t_cand_end).count();

            if (graph.cand.empty()) {
                gcsa_log("[%s] round %zu: no preferences to extract\n", label, rounds);
                ++rounds;
                break;
            }

            gcsa_log("[%s] round %zu: DP on %zu nodes...\n", label, rounds, graph.cand.size());
            auto t_dp0 = Clock::now();
            std::unordered_map<uint64_t, Candidate> chosen =
                solve_pseudoforest_dp_(graph, rank_aware, trace);
            auto t_dp1 = Clock::now();
            dp_ms += std::chrono::duration<double, std::milli>(t_dp1 - t_dp0).count();
            const size_t n_chosen = chosen.size();
            gcsa_log("[%s] round %zu: DP done (cand=%zu cycles=%zu chosen=%zu)\n",
                     label, rounds, graph.cand.size(), graph.cycles.size(), n_chosen);
            total_dp_nodes += graph.cand.size();
            total_cycles += graph.cycles.size();
            total_chosen += n_chosen;

            if (n_chosen) {
                auto t_acc0 = Clock::now();
                accept_pseudoforest_chosen_(chosen, graph.by_name, accepted, trace);
                auto t_acc1 = Clock::now();
                accept_ms += std::chrono::duration<double, std::milli>(t_acc1 - t_acc0).count();
            }
            ++rounds;

            if (!iterate) break;               // single-fire: stop after round 0
            if (n_chosen == 0) break;           // fixed point: nothing left for another round to find
            if (max_rounds && rounds >= max_rounds) {
                gcsa_log("[%s] GCSA_PFDP_MAX_ROUNDS=%zu reached, stopping\n", label, max_rounds);
                break;
            }
        }

        const size_t after_dp = kept_count_;
        gcsa_log(
            "[%s] DP alone (%zu round%s): |C| %zu -> %zu  (%.1f%% of original kept, -%zu)\n",
            label, rounds, rounds == 1 ? "" : "s", m0, after_dp, pct_kept(after_dp), m0 - after_dp);

        gcsa_log("[%s] leftover greedy...\n", label);
        auto t_left0 = Clock::now();
        size_t n_leftover = run_pseudoforest_leftover_(intervals, accepted, trace);
        auto t_left1 = Clock::now();
        const size_t after_leftover = kept_count_;
        gcsa_log(
            "[%s] leftover:  |C| %zu -> %zu  (%.1f%% of original kept, -%zu)\n",
            label, after_dp, after_leftover, pct_kept(after_leftover), after_dp - after_leftover);

        // Phase II + LNS: same combined pass every algorithm now runs (see
        // run_phase2_and_lns_). Pass the real Phase I wall time (graph
        // build through leftover, every round included) for run_phase2_'s
        // time-budget calc.
        auto t_p2_0 = Clock::now();
        double phase1_ms = std::chrono::duration<double, std::milli>(t_left1 - t_phase1_start).count();
        run_phase2_and_lns_(intervals, accepted, label, trace, timing, &cand_cache, phase1_ms);
        auto t_p2_1 = Clock::now();
        const size_t after_phase2 = kept_count_;
        gcsa_log(
            "[%s] phase II:  |C| %zu -> %zu  (%.1f%% of original kept, -%zu)\n",
            label, after_leftover, after_phase2, pct_kept(after_phase2),
            after_leftover - after_phase2);
        gcsa_log(
            "[%s] success: original |C|=%zu -- DP alone=%.1f%%, "
            "+leftover=%.1f%%, +phase II=%.1f%% (of original kept)\n",
            label, m0, pct_kept(after_dp), pct_kept(after_leftover), pct_kept(after_phase2));

        if (timing) {
            auto ms = [](Clock::time_point a, Clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count();
            };
            const double phase2_total = ms(t_p2_0, t_p2_1);
            const double leftover_total = ms(t_left0, t_left1);
            gcsa_log(
                "[timing] %s rounds=%zu cand=%.1fms graph=%.1fms dp=%.1fms accept=%.1fms "
                "leftover=%.1fms phase2=%.1fms total=%.1fms "
                "phase2/phase1=%.2f\n"
                "         dp_nodes=%zu cycles=%zu chosen=%zu n_leftover=%zu "
                "accepted=%zu kept=%zu threads=%d\n",
                label, rounds, cand_ms, graph_ms,
                dp_ms, accept_ms, leftover_total,
                phase2_total, ms(t_all, Clock::now()),
                phase2_total / std::max(1.0, phase1_ms),
                total_dp_nodes, total_cycles, total_chosen, n_leftover,
                accepted.size(), kept_count_, nthreads);
        }
    }

    void compress_() {
        const size_t m = G_.m();
        rank_of_.assign(m, 0);
        for (size_t r = 0; r < m; ++r) rank_of_[G_.sa[r]] = (int32_t)r;
        removed_.assign(m, 0);
        removed_words_.assign((m + 63) / 64, 0);
        pin_count_.assign(m, 0);
        pin_owner_.assign(m, 0);
        multi_pin_owners_.clear();
        kept_count_ = m;
        table_.clear();
        C_.clear();

        gcsa_log("[%s] compressing...\n", algo_name(algo_));

        std::vector<Interval> intervals;
        for (size_t r = 0; r < m; ) {
            uint64_t name = name_of_rank_((int32_t)r);
            size_t r2 = r + 1;
            if (regular_mode_) {
                while (r2 < m && G_.lcp[r2] >= base_depth_) ++r2;
            } else {
                while (r2 < m && G_.first_symbol((int32_t)r2) == name) ++r2;
            }
            intervals.push_back({name, (int32_t)r, (int32_t)r2});
            r = r2;
        }

        std::unordered_map<uint64_t, Candidate> accepted;
        if (algo_ == CompressAlgo::DepOrder)
            compress_dep_order_(intervals, accepted);
        else if (algo_ == CompressAlgo::PseudoforestDp ||
                 algo_ == CompressAlgo::PseudoforestDpIterate)
            compress_pseudoforest_dp_(intervals, accepted);
        else if (algo_ == CompressAlgo::GreedyDegree)
            compress_greedy_degree_(intervals, accepted);
        else
            compress_greedy_size_(intervals, accepted);

        gcsa_log("[%s] finalize: build C / hash table...\n", algo_name(algo_));
        finalize_(intervals, accepted);
        gcsa_log("[%s] finalize done\n", algo_name(algo_));
    }

    // Build C (kept positions in rank order) and the two-entries-per-name hash
    // table from an accepted link set. Shared by every algorithm (compress_()
    // calls this once, after whichever compress_<algo>_ ran) -- neither this
    // function nor any compress_<algo>_'s own [timing] line covered it before;
    // its cost only ever showed up folded into CompressedIndex::build's outer
    // "compress(<algo>)=" total. Two passes:
    //   1. compact:  C_ in rank order. Sequential by construction -- a kept
    //      rank's C_ index depends on how many earlier ranks also survived.
    //   2. table:    one HashEntry per interval (k-mer). Independent per
    //      interval once rank_to_C_ exists from pass 1 -- a parallelization
    //      candidate, not parallelized yet.
    void finalize_(const std::vector<Interval>& intervals,
                   const std::unordered_map<uint64_t, Candidate>& accepted) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        auto t0 = Clock::now();

        const size_t m = G_.m();
        rank_to_C_.assign(m, -1);
        C_.reserve(m);
        for (size_t r = 0; r < m; ++r) {
            if (!removed_[r]) {
                rank_to_C_[r] = (int64_t)C_.size();
                C_.push_back(G_.orig_pos((int32_t)r));
            }
        }
        auto t1 = Clock::now();

        for (const auto& iv : intervals) {
            HashEntry e;
            auto it = accepted.find(iv.name);
            if (it != accepted.end()) {
                const Candidate& c = it->second;
                e.has_offset = true;
                e.off_pos = (uint64_t)rank_to_C_[c.src_lo];
                e.off_add = (uint32_t)c.add;
                e.off_num = (uint32_t)c.covered.size();
            }
            int64_t first = -1; uint32_t cnt = 0;
            for (int32_t r = iv.lo; r < iv.hi; ++r) {
                if (removed_[r]) continue;
                if (first < 0) first = rank_to_C_[r];
                ++cnt;
            }
            if (cnt > 0) { e.has_rest = true; e.rest_pos = (uint64_t)first; e.rest_num = cnt; }
            table_.emplace(iv.name, e);
        }
        auto t2 = Clock::now();

        if (timing) {
            auto ms = [](Clock::time_point a, Clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count();
            };
            gcsa_log(
                "[timing] %s finalize: compact=%.1fms table=%.1fms total=%.1fms "
                "m=%zu intervals=%zu kept=%zu\n",
                algo_name(algo_), ms(t0, t1), ms(t1, t2), ms(t0, t2),
                m, intervals.size(), C_.size());
        }
    }
};

} // namespace gcsa
