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
// universe.  It emits only maximal runs with lcp >= add+1 (the note's
// lcp-interval argument), which keeps the per-interval candidate list short;
// GCSA_INTRA_LINKS=0 additionally restores the historical "source outside I_c"
// rule.  enumerate_all_links_ emits the whole universe and is what
// collect_diff_problem hands to the ILP baseline, so the ILP optimum is a true
// lower bound for every algorithm below.
//
// Compression strategies (CompressAlgo):
//   Greedy    – size-order availability-aware greedy; pins sources forever.
//   DepOrder  – dependency-order + un-pin / retarget.
//   TreeDp    – DP on the preference forest (strong cov>=kMinCoverage edges, |I_c|>2):
//               for each source hub, choose KEEP vs COMPRESS knowing how
//               dependents' costs change; then Phase II unpin/retarget.
//   TreeDp2   – same as TreeDp but without Phase II (forest DP + accept + leftover).
//   PseudoforestDp – DP on the *whole* preference graph (out-degree<=1, no
//               cov/size threshold): tree components use the same KEEP vs
//               COMPRESS recurrence; unicyclic components (e.g. a mutual
//               pair) are solved exactly by folding the cycle instead of
//               dropping an edge to force acyclicity; then Phase II.

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

// Minimum coverage to accept a compression link (all algos + Phase II).
// Preference-forest / DAG *edges* use the same floor (see DepOrder / TreeDp).
constexpr int kMinCoverage = 3;

// Phase II runs this many dirty generations unless it reaches a fixed point
// first. Callers override it per build (CompressedIndex::build's phase2_iters,
// i.e. --phase2-iters), or process-wide with GCSA_PHASE2_MAX_ITERS.
constexpr int kPhase2DefaultIters = 100;

// Above this SA length, Phase II defaults to a single dirty generation (same as
// GCSA_PHASE2_FAST) unless the caller overrides the budget. Keeps Phase II from
// dominating wall time on large DNA / skmer concatenations.
constexpr size_t kPhase2AutoFastM = 1000000;

// Max candidates retained per interval in the shared Phase II / leftover cache.
// Preference pick runs on the full enum first; the cache only needs coverage-
// descending heads for try_accept. Override with GCSA_CAND_CACHE_CAP.
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
    Greedy,        // size-first, pin sources
    DepOrder,      // dependency-order + un-pin / retarget
    TreeDp,        // preference-forest DP (KEEP vs COMPRESS per hub) + Phase II
    TreeDp2,       // same as TreeDp without Phase II
    TreeDp3,       // same as TreeDp with the cluster-LNS Phase II
    TreeDp4,       // same as TreeDp with value-based cycle repair in Phase I
    PseudoforestDp, // exact DP on the full out-degree<=1 preference graph
};

inline const char* algo_name(CompressAlgo a) {
    switch (a) {
        case CompressAlgo::Greedy:        return "greedy";
        case CompressAlgo::DepOrder:      return "dep-order";
        case CompressAlgo::TreeDp:        return "tree-dp";
        case CompressAlgo::TreeDp2:       return "tree-dp2";
        case CompressAlgo::TreeDp3:       return "tree-dp3";
        case CompressAlgo::TreeDp4:       return "tree-dp4";
        case CompressAlgo::PseudoforestDp: return "pseudoforest-dp";
    }
    return "?";
}

// What compress_tree_dp_ runs after the forest DP + leftover greedy.
enum class Phase2Mode {
    None,    // tree-dp2
    Greedy,  // tree-dp: dirty-set unpin/retarget
    Local,   // tree-dp3: exact cluster large-neighborhood search
};

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
               int max_add = 8, CompressAlgo algo = CompressAlgo::Greedy,
               int phase2_iters = 0) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        auto t0 = Clock::now();
        G_ = build_gapped_sa(shape, std::move(text));
        auto t1 = Clock::now();
        span_ = shape.span;
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
                out.push_back(C_[e.off_pos + t] + (int64_t)e.off_add * span_);
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
        span_ = shape.span;
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
            uint64_t name = G_.first_symbol((int32_t)r);
            size_t r2 = r + 1;
            while (r2 < m && G_.first_symbol((int32_t)r2) == name) ++r2;
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

private:
    GappedSA G_;
    int span_ = 1;
    int max_add_ = 8;
    CompressAlgo algo_ = CompressAlgo::Greedy;
    int phase2_iters_ = 0;  // 0 = unspecified; see run_phase2_ for precedence

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
        if (G_.lex2orig[y] != G_.lex2orig[x] - (int64_t)add * span_) return -1;
        return y;
    }

    // Inverse of pred_lexpos_: lextext position add symbols after sa[r], if the
    // original-text shift is exactly +add*span (same residue class / window).
    int64_t succ_lexpos_(int32_t r, int add) const {
        int64_t x = G_.sa[r];
        int64_t y = x + add;
        if (y >= (int64_t)G_.m()) return -1;
        if (G_.lex2orig[y] != G_.lex2orig[x] + (int64_t)add * span_) return -1;
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
                while (j < pr.size()
                       && pr[j].first == pr[j-1].first + 1
                       && G_.lcp[pr[j].first] >= add + 1) ++j;
                const int32_t s_lo = pr[i].first, s_hi = pr[j-1].first + 1;
                if (s_hi <= lo || s_lo >= hi)
                    emit_run_(name, lo, hi, add, pr, i, j, cov_floor,
                              require_avail, out);
                else if (intra)
                    emit_intra_runs_(name, lo, hi, add, pr, i, j, cov_floor,
                                     require_avail, out);
                i = j;
            }
        }
        return out;
    }

    // Turn the source sub-run pr[i..j) into a candidate for I_c = [lo,hi).
    void emit_run_(uint64_t name, int32_t lo, int32_t hi, int add,
                   const std::vector<std::pair<int32_t,int32_t>>& pr,
                   size_t i, size_t j, int cov_floor, bool require_avail,
                   std::vector<Candidate>& out) const {
        const int cov = (int)(j - i);
        if (cov < cov_floor) return;
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
        out.push_back(std::move(c));
    }

    // A source run reaching into I_c is still decodable as long as it holds
    // none of the ranks it covers (see "Link universe"), so instead of dropping
    // it we emit its maximal such sub-runs. Two-pointer scan in run-relative
    // coordinates: the window [p,q) is bad iff some t in [p,q) covers a rank in
    // [p,q), i.e. iff it contains both endpoints of one of the `forbid` pairs.
    void emit_intra_runs_(uint64_t name, int32_t lo, int32_t hi, int add,
                          const std::vector<std::pair<int32_t,int32_t>>& pr,
                          size_t i, size_t j, int cov_floor, bool require_avail,
                          std::vector<Candidate>& out) const {
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
                          cov_floor, require_avail, out);
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
                const uint64_t name = G_.first_symbol(tgt[(size_t)a]);
                int32_t b = a + 1;
                while (b < m && tgt[(size_t)b] >= 0
                       && G_.first_symbol(tgt[(size_t)b]) == name) ++b;
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

    // First-symbol name of an SA rank.
    uint64_t name_of_rank_(int32_t r) const { return G_.first_symbol(r); }

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
    // Shared by DepOrder Phase II and TreeDp Phase II.
    //
    // Candidates with require_avail=false depend only on SA/LCP structure, so we
    // enumerate once per interval and reuse across generations. Optional
    // `precomputed` cache (e.g. from TreeDp preference build) avoids re-enum.
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
    // Env GCSA_DISABLE_PHASE2 (optional): skip Phase II entirely, so callers
    // (DepOrder, TreeDp, PseudoforestDp) can be compared on their Phase I /
    // DP output alone, without the shared retarget sweep smoothing over
    // differences between them.
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

        auto get_cands = [&](const Interval& iv) -> const std::vector<Candidate>& {
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

    // ---- algorithms -------------------------------------------------------
    void compress_greedy_(const std::vector<Interval>& intervals,
                          std::unordered_map<uint64_t, Candidate>& accepted) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        auto t0 = Clock::now();
        gcsa_log("[greedy] size-order accept...\n");
        std::vector<const Interval*> order;
        for (const auto& iv : intervals) if (iv.hi - iv.lo > 1) order.push_back(&iv);
        std::sort(order.begin(), order.end(),
                  [](const Interval* a, const Interval* b){ return (a->hi-a->lo) > (b->hi-b->lo); });
        // Optional trace: set GCSA_TRACE_GREEDY=1
        const bool trace = (std::getenv("GCSA_TRACE_GREEDY") != nullptr);
        int step = 0;
        if (trace) {
            gcsa_log("greedy order:");
            for (auto* ivp : order)
                gcsa_log(" %s(size=%d)",
                         name_to_string(G_.shape, ivp->name).c_str(), ivp->hi - ivp->lo);
            gcsa_log("\n");
        }
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
        gcsa_log("[greedy] size-order accept done\n");
        if (timing) {
            double ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
            gcsa_log("[timing] greedy: %.1fms  (#I>1=%zu accepted=%zu)\n",
                     ms, order.size(), accepted.size());
        }
    }

    void compress_dep_order_(const std::vector<Interval>& intervals,
                             std::unordered_map<uint64_t, Candidate>& accepted) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        auto t_all = Clock::now();
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

        double phase1_ms = std::chrono::duration<double, std::milli>(t3 - t0).count();
        run_phase2_(intervals, accepted, "dep-order", trace, timing, &cand_cache, phase1_ms);
        auto t4 = Clock::now();
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

    // Preference-forest DP: each strong preference edge (cov>=kMinCoverage, |I|>2) makes
    // the target a child of its preferred source.  For every node we choose
    // KEEP (all ranks stored, children may compress against us) vs COMPRESS
    // (drop cov ranks via the preferred candidate; children lose that source).
    // Exact for |C| on this forest model; leftover names get an avail. greedy.
    // phase2: which Phase II follows — none (TreeDp2), the dirty-set
    // unpin/retarget (TreeDp), or the exact cluster LNS (TreeDp3).
    //
    // Hot-path notes (large N):
    //   - Preference enum is independent per interval → GCSA_THREADS parallel.
    //   - cand_cache reused for leftover + Phase II (no re-enum).
    //   - Forest DP uses dense node ids + vector memo; roots are independent.
    void compress_tree_dp_(const std::vector<Interval>& intervals,
                           std::unordered_map<uint64_t, Candidate>& accepted,
                           Phase2Mode phase2 = Phase2Mode::Greedy,
                           bool cycle_repair = false) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        const bool trace = (std::getenv("GCSA_TRACE_DP") != nullptr);
        const char* label = algo_name(algo_);
        auto t_all = Clock::now();
        const int nthreads = gcsa_num_threads();

        std::unordered_map<uint64_t, const Interval*> by_name;
        by_name.reserve(intervals.size() * 2);
        for (const auto& iv : intervals) by_name[iv.name] = &iv;

        struct Pref { const Interval* iv; Candidate cand; uint64_t src_name; };
        std::unordered_map<uint64_t, Pref> prefs;   // target -> preferred
        // Static (avail=false) candidate lists — reused by leftover + Phase II.
        // Prefs use pick_best_ on enum order (stable tie-break); Phase II wants
        // coverage-desc order, so we sort only after recording the preferred.
        std::unordered_map<uint64_t, std::vector<Candidate>> cand_cache;

        gcsa_log("[%s] preference forest build...\n", label);
        if (timing)
            gcsa_log("[timing] %s pref threads=%d\n", label, nthreads);
        auto t_pref0 = Clock::now();

        std::vector<const Interval*> big_ivs;
        big_ivs.reserve(intervals.size());
        for (const auto& iv : intervals)
            if (iv.hi - iv.lo > 1) big_ivs.push_back(&iv);

        // Stream into prefs/cand_cache under a mutex instead of retaining a
        // PrefWork[] the size of |{I : |I|>1}| — that doubled peak RAM on large
        // skmer concatenations (hundreds of thousands of intervals) and OOM'd.
        prefs.reserve(big_ivs.size());
        cand_cache.reserve(big_ivs.size() * 2);
        std::mutex pref_mu;
        gcsa_parallel_for(big_ivs.size(), [&](size_t i) {
            const Interval& iv = *big_ivs[i];
            auto cands = enumerate_candidates_(iv.name, iv.lo, iv.hi, /*avail=*/false);
            bool has_pref = false;
            Candidate best;
            uint64_t src_name = 0;
            if (iv.hi - iv.lo > 2) {
                best = pick_best_(cands);
                if (best.coverage() >= kMinCoverage) {
                    has_pref = true;
                    src_name = name_of_rank_(best.src_lo);
                }
            }
            std::sort(cands.begin(), cands.end(), [](const Candidate& a, const Candidate& b){
                if (a.coverage() != b.coverage()) return a.coverage() > b.coverage();
                return a.add > b.add;
            });
            trim_cand_cache_(cands);
            std::lock_guard<std::mutex> lock(pref_mu);
            if (has_pref)
                prefs.emplace(iv.name, Pref{&iv, std::move(best), src_name});
            cand_cache.emplace(iv.name, std::move(cands));
        });
        auto t_pref1 = Clock::now();

        // Forest edges: only cov>=kMinCoverage (same threshold as DepOrder DAG edges).
        std::unordered_map<uint64_t, std::vector<uint64_t>> children;
        std::unordered_map<uint64_t, uint64_t> parent;
        // Cycle check: t->s would cycle iff s already has t as an ancestor
        // (parent-chain walk). Equivalent to the old downward DFS reaches(t,s).
        auto has_ancestor = [&](uint64_t s, uint64_t t) {
            for (uint64_t u = s; ;) {
                if (u == t) return true;
                auto it = parent.find(u);
                if (it == parent.end()) return false;
                u = it->second;
            }
        };

        auto t_forest0 = Clock::now();
        // Deterministic edge order (name asc) so parallel pref build does not
        // change which cycle-breaking edges survive.
        std::vector<uint64_t> pref_names;
        pref_names.reserve(prefs.size());
        for (auto& kv : prefs) pref_names.push_back(kv.first);
        std::sort(pref_names.begin(), pref_names.end());
        size_t cycles_seen = 0, cycles_repaired = 0, cycle_gain = 0;
        const int cycle_min_gain = gcsa_env_int("GCSA_CYCLE_MIN_GAIN", 1);
        for (uint64_t t : pref_names) {
            auto& p = prefs[t];
            uint64_t s = p.src_name;
            if (p.cand.coverage() < kMinCoverage) continue;
            if (s == t) continue;
            if (!by_name.count(s)) continue;
            if (!has_ancestor(s, t)) {
                parent[t] = s;
                children[s].push_back(t);
                continue;
            }
            // t -> s would close a cycle. Which edge of that cycle dies is
            // otherwise decided by name order, not by value: drop the
            // least-valuable one instead. The node whose edge is dropped is the
            // one that cannot compress via its preference, so the loss is that
            // edge's coverage. Adding t -> s creates exactly this one cycle, so
            // removing any single edge of it leaves the forest acyclic.
            ++cycles_seen;
            if (!cycle_repair) continue;
            uint64_t victim = t;
            int victim_cov = p.cand.coverage();
            for (uint64_t u = s; u != t; u = parent[u]) {
                auto pit = prefs.find(u);
                if (pit == prefs.end()) { victim = t; break; }
                int cov = pit->second.cand.coverage();
                if (cov < victim_cov) { victim_cov = cov; victim = u; }
            }
            // Rewiring mid-chain is more disruptive than dropping t's edge (t is
            // still a leaf here), so only do it when the edge we save is
            // meaningfully better than the one we sacrifice.
            if (victim == t) continue;  // name order already picked the cheapest
            if (p.cand.coverage() - victim_cov < cycle_min_gain) continue;
            uint64_t vp = parent[victim];
            auto& sib = children[vp];
            sib.erase(std::remove(sib.begin(), sib.end(), victim), sib.end());
            parent.erase(victim);
            parent[t] = s;
            children[s].push_back(t);
            ++cycles_repaired;
            cycle_gain += (size_t)(p.cand.coverage() - victim_cov);
        }

        std::unordered_set<uint64_t> nodes;
        for (auto& kv : prefs) nodes.insert(kv.first);
        for (auto& kv : children) {
            nodes.insert(kv.first);
            for (uint64_t c : kv.second) nodes.insert(c);
        }
        std::vector<uint64_t> roots;
        for (uint64_t n : nodes)
            if (!parent.count(n)) roots.push_back(n);
        std::sort(roots.begin(), roots.end());
        auto t_forest1 = Clock::now();

        if (trace) {
            gcsa_log("=== %s forest (src -> dependents) ===\n", label);
            for (auto& kv : children) {
                gcsa_log("  %s ->",
                         name_to_string(G_.shape, kv.first).c_str());
                for (uint64_t c : kv.second)
                    gcsa_log(" %s",
                             name_to_string(G_.shape, c).c_str());
                gcsa_log("\n");
            }
            gcsa_log("roots:");
            for (uint64_t r : roots)
                gcsa_log(" %s", name_to_string(G_.shape, r).c_str());
            gcsa_log("\n");
        }

        // Dense forest arrays for DP.
        std::vector<uint64_t> node_of;
        node_of.reserve(nodes.size());
        for (uint64_t n : nodes) node_of.push_back(n);
        std::sort(node_of.begin(), node_of.end());
        std::unordered_map<uint64_t, int> id_of;
        id_of.reserve(node_of.size() * 2);
        for (size_t i = 0; i < node_of.size(); ++i) id_of[node_of[i]] = (int)i;
        const int NN = (int)node_of.size();

        std::vector<std::vector<int>> ch(NN);
        std::vector<int> isize_arr(NN, 0);
        std::vector<int> pref_cov(NN, -1);          // -1 = cannot compress via pref
        std::vector<const Candidate*> pref_ptr(NN, nullptr);
        for (int i = 0; i < NN; ++i) {
            auto it = by_name.find(node_of[i]);
            if (it != by_name.end())
                isize_arr[i] = it->second->hi - it->second->lo;
            auto pit = prefs.find(node_of[i]);
            if (pit != prefs.end() && pit->second.cand.coverage() >= kMinCoverage) {
                pref_cov[i] = pit->second.cand.coverage();
                pref_ptr[i] = &pit->second.cand;
            }
        }
        for (auto& kv : children) {
            int p = id_of[kv.first];
            for (uint64_t c : kv.second) ch[p].push_back(id_of[c]);
        }

        struct Cell { int cost; bool compress; };
        // Per-root DP (disjoint trees) — parallelize over roots.
        std::vector<std::vector<std::pair<int, Candidate>>> chosen_by_root(roots.size());

        gcsa_log("[%s] DP on forest...\n", label);
        auto t_dp0 = Clock::now();
        gcsa_parallel_for(roots.size(), [&](size_t ri) {
            const int root = id_of.at(roots[ri]);
            // memo[2*v + src_ok]: cost < 0 => empty
            std::vector<Cell> memo((size_t)NN * 2, Cell{-1, false});

            std::function<Cell(int, bool)> solve = [&](int v, bool src_ok) -> Cell {
                Cell& slot = memo[(size_t)v * 2 + (src_ok ? 1 : 0)];
                if (slot.cost >= 0) return slot;

                const auto& kids = ch[v];
                int cost_keep = isize_arr[v];
                for (int w : kids) cost_keep += solve(w, /*src_ok=*/true).cost;

                int cost_comp = std::numeric_limits<int>::max();
                bool can_comp = false;
                if (src_ok && pref_cov[v] >= kMinCoverage) {
                    can_comp = true;
                    int cc = isize_arr[v] - pref_cov[v];
                    for (int w : kids) cc += solve(w, /*src_ok=*/false).cost;
                    cost_comp = cc;
                }

                if (can_comp && cost_comp < cost_keep)
                    slot = {cost_comp, true};
                else
                    slot = {cost_keep, false};
                return slot;
            };

            auto& local = chosen_by_root[ri];
            std::function<void(int, bool)> apply = [&](int v, bool src_ok) {
                Cell cell = solve(v, src_ok);
                if (cell.compress) {
                    local.push_back({v, *pref_ptr[v]});
                    for (int w : ch[v]) apply(w, /*src_ok=*/false);
                } else {
                    for (int w : ch[v]) apply(w, /*src_ok=*/true);
                }
            };

            // Root: KEEP, or COMPRESS if preferred source is outside this tree.
            const auto& kids = ch[root];
            int cost_keep = isize_arr[root];
            for (int w : kids) cost_keep += solve(w, true).cost;

            int cost_comp = std::numeric_limits<int>::max();
            bool can_comp = false;
            if (pref_cov[root] >= kMinCoverage) {
                auto pit = prefs.find(node_of[root]);
                uint64_t s = pit->second.src_name;
                if (s != node_of[root] && !has_ancestor(s, node_of[root])) {
                    can_comp = true;
                    int cc = isize_arr[root] - pref_cov[root];
                    for (int w : kids) cc += solve(w, false).cost;
                    cost_comp = cc;
                }
            }

            if (trace) {
                // fprintf is not great under parallel; only print when single-threaded.
                if (nthreads <= 1)
                    gcsa_log("root %s: keep=%d comp=%s\n",
                             name_to_string(G_.shape, node_of[root]).c_str(), cost_keep,
                             can_comp ? std::to_string(cost_comp).c_str() : "n/a");
            }

            if (can_comp && cost_comp < cost_keep) {
                local.push_back({root, *pref_ptr[root]});
                for (int w : kids) apply(w, false);
            } else {
                for (int w : kids) apply(w, true);
            }
        });

        std::unordered_map<uint64_t, Candidate> chosen;
        chosen.reserve(nodes.size());
        for (auto& vec : chosen_by_root)
            for (auto& p : vec)
                chosen.emplace(node_of[p.first], std::move(p.second));
        auto t_dp1 = Clock::now();
        gcsa_log("[%s] DP done\n", label);

        // Materialize: accept deepest dependents first (pins hubs before hubs decide).
        gcsa_log("[%s] accept chosen...\n", label);
        auto t_acc0 = Clock::now();
        const std::vector<uint64_t> kEmptyKids;
        auto children_of = [&](uint64_t v) -> const std::vector<uint64_t>& {
            auto it = children.find(v);
            return it == children.end() ? kEmptyKids : it->second;
        };
        std::function<void(uint64_t)> accept_subtree = [&](uint64_t v) {
            for (uint64_t w : children_of(v)) accept_subtree(w);
            auto it = chosen.find(v);
            if (it == chosen.end()) return;
            Candidate c = it->second;
            bool ok = cand_available_(c);
            if (!ok) {
                auto cit = cand_cache.find(v);
                if (cit != cand_cache.end())
                    c = best_from_sorted_cache_(cit->second, /*require_avail=*/true);
                else {
                    const Interval* iv = by_name[v];
                    c = best_candidate_(iv->name, iv->lo, iv->hi, /*avail=*/true);
                }
            }
            if (c.coverage() >= kMinCoverage) {
                if (trace) {
                    gcsa_log("  ACCEPT %s <- %s add=%d cov=%d\n",
                        name_to_string(G_.shape, v).c_str(),
                        name_to_string(G_.shape, name_of_rank_(c.src_lo)).c_str(),
                        c.add, c.coverage());
                }
                accept_(c, accepted);
            }
        };
        for (uint64_t r : roots) accept_subtree(r);
        auto t_acc1 = Clock::now();

        // Leftover intervals: greedy with availability, using static cand_cache.
        gcsa_log("[%s] leftover greedy...\n", label);
        auto t_left0 = Clock::now();
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
        for (const Interval* ivp : leftover) {
            Candidate c;
            auto cit = cand_cache.find(ivp->name);
            if (cit != cand_cache.end())
                c = best_from_sorted_cache_(cit->second, /*require_avail=*/true);
            else
                c = best_candidate_(ivp->name, ivp->lo, ivp->hi, /*avail=*/true);
            if (c.coverage() < kMinCoverage) continue;
            accept_(c, accepted);
        }
        auto t_left1 = Clock::now();
        gcsa_log("[%s] leftover done\n", label);

        // Phase II (reuses cand_cache either way).
        auto t_p2_0 = Clock::now();
        double phase1_ms = std::chrono::duration<double, std::milli>(t_left1 - t_pref0).count();
        if (phase2 == Phase2Mode::Greedy) {
            run_phase2_(intervals, accepted, label, trace, timing, &cand_cache, phase1_ms);
        } else if (phase2 == Phase2Mode::Local) {
            // The cluster solve only ever lowers |C|, so it composes with the
            // retarget loop instead of replacing it: retargeting reaches links
            // (composite adds) the static candidate cache does not contain,
            // which is exactly where the LNS alone gives ground on long
            // repeats. GCSA_LNS_ONLY=1 measures the LNS on its own.
            if (gcsa_env_int("GCSA_LNS_ONLY", 0) == 0)
                run_phase2_(intervals, accepted, label, trace, timing, &cand_cache, phase1_ms);
            run_phase2_local_(intervals, accepted, label, timing, &cand_cache);
        }
        auto t_p2_1 = Clock::now();

        if (timing) {
            auto ms = [](Clock::time_point a, Clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count();
            };
            const double phase1_total = ms(t_pref0, t_left1);
            const double phase2_total = ms(t_p2_0, t_p2_1);
            gcsa_log(
                "[timing] %s pref=%.1fms forest=%.1fms dp=%.1fms accept=%.1fms "
                "leftover=%.1fms phase2=%.1fms total=%.1fms "
                "phase2/phase1=%.2f\n"
                "         roots=%zu forest_nodes=%zu prefs=%zu #I>1=%zu "
                "accepted=%zu kept=%zu threads=%d "
                "cycles=%zu repaired=%zu edge_gain=%zu\n",
                label, ms(t_pref0, t_pref1), ms(t_forest0, t_forest1),
                ms(t_dp0, t_dp1), ms(t_acc0, t_acc1), ms(t_left0, t_left1),
                phase2_total, ms(t_all, Clock::now()),
                phase2_total / std::max(1.0, phase1_total),
                roots.size(), nodes.size(), prefs.size(), cand_cache.size(),
                accepted.size(), kept_count_, nthreads,
                cycles_seen, cycles_repaired, cycle_gain);
        }
    }

    // Pseudoforest DP: every interval with |I_c|>=2 and a preferred candidate
    // of coverage>=2 gets exactly one outgoing "prefers" edge, to whichever
    // shape holds its source (a candidate's source always lies entirely
    // inside a single shape's interval, since LCP>=add+1>=2 forces every
    // source rank to share its own first symbol). A graph with out-degree<=1
    // per node is a pseudoforest: every connected component is either a tree
    // (rooted at a node with no preference) or unicyclic (exactly one cycle,
    // with trees hanging off each cycle node) -- never anything worse.
    //
    // Unlike TreeDp, no edge is ever dropped to force acyclicity. Tree
    // components use the same KEEP/COMPRESS recurrence as TreeDp; unicyclic
    // components are solved exactly by folding the cycle into two forced
    // sub-cases (anchor node forced KEEP, or forced COMPRESS with its
    // successor forced KEEP) and taking whichever is cheaper. Because
    // cycles are solved rather than avoided, no cov>=3 or |I_c|>2 threshold
    // is needed either -- every compressible interval participates.
    void compress_pseudoforest_dp_(const std::vector<Interval>& intervals,
                                   std::unordered_map<uint64_t, Candidate>& accepted) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        const bool trace = (std::getenv("GCSA_TRACE_PFDP") != nullptr);
        auto t_all = Clock::now();

        std::unordered_map<uint64_t, const Interval*> by_name;
        for (const auto& iv : intervals) by_name[iv.name] = &iv;

        struct Pref { Candidate cand; uint64_t src_name; };
        std::unordered_map<uint64_t, Pref> prefs;   // node -> its unique preference
        std::unordered_map<uint64_t, std::vector<Candidate>> cand_cache;  // reused by Phase II

        std::fprintf(stderr, "[pseudoforest-dp] preference graph build...\n");
        for (const auto& iv : intervals) {
            if (iv.hi - iv.lo < 2) continue;
            auto cands = enumerate_candidates_(iv.name, iv.lo, iv.hi, /*avail=*/false);
            Candidate best = pick_best_(cands);
            if (best.coverage() >= 2) {
                uint64_t src = name_of_rank_(best.src_lo);
                if (src != iv.name && by_name.count(src))
                    prefs.emplace(iv.name, Pref{best, src});
            }
            std::sort(cands.begin(), cands.end(), [](const Candidate& a, const Candidate& b){
                if (a.coverage() != b.coverage()) return a.coverage() > b.coverage();
                return a.add > b.add;
            });
            cand_cache.emplace(iv.name, std::move(cands));
        }

        std::unordered_set<uint64_t> all_nodes;
        for (const auto& iv : intervals) if (iv.hi - iv.lo >= 2) all_nodes.insert(iv.name);

        auto isize = [&](uint64_t n) -> int {
            auto it = by_name.find(n);
            return it == by_name.end() ? 0 : (it->second->hi - it->second->lo);
        };
        auto has_pref = [&](uint64_t n) { return prefs.count(n) != 0; };
        auto pref_src = [&](uint64_t n) { return prefs.at(n).src_name; };

        // ---- pure graph decomposition: find every cycle up front ----------
        // (independent of the DP; identifies exactly which single incoming
        // edge per cycle node is "cycle-internal" so it can be excluded from
        // the generic children[] map used by solve()/apply() below -- that
        // edge is instead handled explicitly by fold_cycle's scenario walk,
        // never both, which would double count it.)
        enum class Color : uint8_t { White, Gray, Black };
        std::unordered_map<uint64_t, Color> color;
        for (uint64_t n : all_nodes) color[n] = Color::White;
        std::vector<std::vector<uint64_t>> cycles;
        std::unordered_map<uint64_t, uint64_t> cycle_pred;  // cyc[i] -> cyc[i-1] (to exclude)

        for (uint64_t start : all_nodes) {
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
                if (!has_pref(v)) break;   // tree root: end of chain
                v = pref_src(v);
            }
            if (hit_cycle) {
                auto it = std::find(path.begin(), path.end(), cycle_entry);
                std::vector<uint64_t> cyc(it, path.end());
                int L = (int)cyc.size();
                for (int i = 0; i < L; ++i)
                    cycle_pred[cyc[i]] = cyc[(i + L - 1) % L];
                cycles.push_back(std::move(cyc));
            }
            for (uint64_t u : path) color[u] = Color::Black;
        }

        // children[w] = dependents of w, EXCLUDING the one cycle-internal
        // edge per cycle node (handled explicitly by fold_cycle instead).
        std::unordered_map<uint64_t, std::vector<uint64_t>> children;
        for (auto& kv : prefs) {
            uint64_t v = kv.first, src = kv.second.src_name;
            auto it = cycle_pred.find(src);
            if (it != cycle_pred.end() && it->second == v) continue;
            children[src].push_back(v);
        }

        if (trace) {
            std::fprintf(stderr, "=== pseudoforest-dp preference graph ===\n");
            for (uint64_t v : all_nodes) {
                if (!has_pref(v)) continue;
                std::fprintf(stderr, "  %s -> %s (cov=%d)\n",
                    name_to_string(G_.shape, v).c_str(),
                    name_to_string(G_.shape, pref_src(v)).c_str(),
                    prefs.at(v).cand.coverage());
            }
            std::fprintf(stderr, "cycles found: %zu\n", cycles.size());
            for (auto& cyc : cycles) {
                std::fprintf(stderr, "  [");
                for (size_t i = 0; i < cyc.size(); ++i)
                    std::fprintf(stderr, "%s%s", i ? "," : "",
                                 name_to_string(G_.shape, cyc[i]).c_str());
                std::fprintf(stderr, "]\n");
            }
        }

        // ---- KEEP / COMPRESS DP (identical recurrence to TreeDp) ----------
        struct Cell { long long cost; bool compress; };
        std::map<std::pair<uint64_t, int>, Cell> memo;

        std::function<Cell(uint64_t, bool)> solve = [&](uint64_t v, bool src_ok) -> Cell {
            auto key = std::make_pair(v, src_ok ? 1 : 0);
            auto it = memo.find(key);
            if (it != memo.end()) return it->second;

            const auto& ch = children[v];
            long long cost_keep = isize(v);
            for (uint64_t w : ch) cost_keep += solve(w, /*src_ok=*/true).cost;

            long long cost_comp = std::numeric_limits<long long>::max();
            bool can_comp = false;
            if (src_ok && has_pref(v)) {
                can_comp = true;
                long long cc = isize(v) - prefs.at(v).cand.coverage();
                for (uint64_t w : ch) cc += solve(w, /*src_ok=*/false).cost;
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
                chosen[v] = prefs.at(v).cand;
                for (uint64_t w : children[v]) apply(w, /*src_ok=*/false);
            } else {
                for (uint64_t w : children[v]) apply(w, /*src_ok=*/true);
            }
        };

        // Fold a unicyclic component's cycle cyc[0]->cyc[1]->...->cyc[L-1]->cyc[0]
        // (cyc[i]'s preferred source is cyc[i+1 mod L]) into two forced
        // sub-cases. cyc[0]'s own compress-ability requires walking all the
        // way around and back to itself -- a circular dependency broken by
        // fixing cyc[0]'s status externally instead of deriving it.
        auto fold_cycle = [&](const std::vector<uint64_t>& cyc) {
            int L = (int)cyc.size();

            // Scenario A: cyc[0] forced KEEP (always valid on its own).
            long long costA = isize(cyc[0]);
            for (uint64_t w : children[cyc[0]]) costA += solve(w, true).cost;
            std::vector<char> statusA(L, 1);
            {
                bool avail = true;   // cyc[0] is available
                for (int i = L - 1; i >= 1; --i) {
                    Cell cell = solve(cyc[i], avail);
                    costA += cell.cost;
                    statusA[i] = cell.compress ? 0 : 1;
                    avail = (statusA[i] == 1);
                }
            }

            // Scenario B: cyc[0] forced COMPRESS -> cyc[1] forced KEEP (must
            // be physically stored for cyc[0] to point at). Walk the rest of
            // the cycle from cyc[L-1] with cyc[0] unavailable.
            long long costB = isize(cyc[0]) - prefs.at(cyc[0]).cand.coverage();
            for (uint64_t w : children[cyc[0]]) costB += solve(w, false).cost;
            {
                long long cyc1_keep = isize(cyc[1]);
                for (uint64_t w : children[cyc[1]]) cyc1_keep += solve(w, true).cost;
                costB += cyc1_keep;
            }
            std::vector<char> statusB(L, 1);
            statusB[0] = 0;
            {
                bool avail = false;  // cyc[0] unavailable
                for (int i = L - 1; i >= 2; --i) {
                    Cell cell = solve(cyc[i], avail);
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
                    chosen[v] = prefs.at(v).cand;
                    for (uint64_t w : children[v]) apply(w, /*src_ok=*/false);
                } else {
                    for (uint64_t w : children[v]) apply(w, /*src_ok=*/true);
                }
            }
        };

        std::fprintf(stderr, "[pseudoforest-dp] DP on components...\n");
        size_t n_roots = 0;
        for (uint64_t v : all_nodes) {
            if (has_pref(v)) continue;
            apply(v, /*src_ok=*/true);
            ++n_roots;
        }
        for (auto& cyc : cycles) fold_cycle(cyc);
        std::fprintf(stderr, "[pseudoforest-dp] DP done (roots=%zu cycles=%zu)\n",
                     n_roots, cycles.size());

        // iter_dp groundwork: how many shapes with a preferred candidate
        // (a `prefs` entry) actually got compressed via it (`chosen`)? Every
        // prefs entry is visited exactly once by apply()/fold_cycle() above,
        // and only ever lands in `chosen` when the DP picked compress=true
        // for it, so prefs.size()-chosen.size() is exactly the number of
        // shapes whose preferred candidate lost to KEEP in this pseudoforest
        // -- the population an iterative second round would re-examine.
        // Reported before Materialize/leftover/Phase II so it isn't muddied
        // by their separate (non-preferred-candidate) fixups.
        {
            const size_t n_prefs = prefs.size();
            const size_t n_orphaned = n_prefs - chosen.size();
            std::fprintf(stderr,
                "[pseudoforest-dp] preferred candidates: total=%zu accepted_by_dp=%zu "
                "orphaned=%zu (%.1f%% orphaned)\n",
                n_prefs, chosen.size(), n_orphaned,
                n_prefs ? 100.0 * (double)n_orphaned / (double)n_prefs : 0.0);
        }

        // Materialize. Every accepted candidate's source shape is, by DP
        // construction, never itself compressed in the same solution (a
        // compress choice is only ever taken under src_ok=true, i.e. when
        // the source was already decided KEPT), so accept order across
        // `chosen` cannot conflict -- no dependency ordering is required.
        std::fprintf(stderr, "[pseudoforest-dp] accept chosen...\n");
        for (auto& kv : chosen) {
            uint64_t v = kv.first;
            Candidate c = kv.second;
            bool ok = run_kept_(c.src_lo, c.src_hi);
            for (int32_t rr : c.covered)
                if (removed_[rr] || pin_count_[rr] > 0) ok = false;
            if (!ok) {
                const Interval* iv = by_name[v];
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

        // Leftover: should be empty by construction (every |I_c|>=2 interval
        // is either a root or belongs to exactly one component above); kept
        // only as defensive insurance, matching the other algorithms.
        std::fprintf(stderr, "[pseudoforest-dp] leftover greedy...\n");
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
        for (const Interval* ivp : leftover) {
            Candidate c = best_candidate_(ivp->name, ivp->lo, ivp->hi, /*avail=*/true);
            if (c.coverage() < 2) continue;
            accept_(c, accepted);
        }
        std::fprintf(stderr, "[pseudoforest-dp] leftover done\n");

        // Phase II: same unpin/retarget fixed point as DepOrder/TreeDp.
        run_phase2_(intervals, accepted, "pseudoforest-dp", trace, timing, &cand_cache);

        if (timing) {
            double ms = std::chrono::duration<double, std::milli>(Clock::now() - t_all).count();
            std::fprintf(stderr,
                "[timing] pseudoforest-dp total: %.1fms  roots=%zu cycles=%zu "
                "nodes=%zu accepted=%zu kept=%zu\n",
                ms, n_roots, cycles.size(), all_nodes.size(), accepted.size(), kept_count_);
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
            uint64_t name = G_.first_symbol((int32_t)r);
            size_t r2 = r + 1;
            while (r2 < m && G_.first_symbol((int32_t)r2) == name) ++r2;
            intervals.push_back({name, (int32_t)r, (int32_t)r2});
            r = r2;
        }

        std::unordered_map<uint64_t, Candidate> accepted;
        if (algo_ == CompressAlgo::DepOrder)
            compress_dep_order_(intervals, accepted);
        else if (algo_ == CompressAlgo::TreeDp)
            compress_tree_dp_(intervals, accepted, Phase2Mode::Greedy);
        else if (algo_ == CompressAlgo::TreeDp2)
            compress_tree_dp_(intervals, accepted, Phase2Mode::None);
        else if (algo_ == CompressAlgo::TreeDp3)
            compress_tree_dp_(intervals, accepted, Phase2Mode::Local);
        else if (algo_ == CompressAlgo::TreeDp4)
            compress_tree_dp_(intervals, accepted, Phase2Mode::Greedy,
                              /*cycle_repair=*/true);
        else if (algo_ == CompressAlgo::PseudoforestDp)
            compress_pseudoforest_dp_(intervals, accepted);
        else
            compress_greedy_(intervals, accepted);

        gcsa_log("[%s] finalize: build C / hash table...\n", algo_name(algo_));
        finalize_(intervals, accepted);
        gcsa_log("[%s] finalize done\n", algo_name(algo_));
    }

    // Build C (kept positions in rank order) and the two-entries-per-name hash
    // table from an accepted link set.
    void finalize_(const std::vector<Interval>& intervals,
                   const std::unordered_map<uint64_t, Candidate>& accepted) {
        const size_t m = G_.m();
        rank_to_C_.assign(m, -1);
        C_.reserve(m);
        for (size_t r = 0; r < m; ++r) {
            if (!removed_[r]) {
                rank_to_C_[r] = (int64_t)C_.size();
                C_.push_back(G_.orig_pos((int32_t)r));
            }
        }
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
    }
};

} // namespace gcsa
