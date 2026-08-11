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
// Compression strategies (CompressAlgo):
//   Greedy    – size-order availability-aware greedy; pins sources forever.
//   DepOrder  – dependency-order + un-pin / retarget.
//   GreedyDfs – pick the add=+1 source hub of greatest total coverage; DFS to
//               its add=+1 targets (|I|>2), then their add=+1 targets, etc.,
//               always pointing at the DFS root with cumulative add (1,2,3,...).
//   TreeDp    – DP on the preference forest (strong cov>=kMinCoverage edges, |I_c|>2):
//               for each source hub, choose KEEP vs COMPRESS knowing how
//               dependents' costs change; then Phase II unpin/retarget.

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
#include <map>
#include <utility>
#include <limits>

namespace gcsa {

// Minimum coverage to accept a compression link (all algos + Phase II).
// Preference-forest / DAG *edges* use the same floor (see DepOrder / TreeDp).
constexpr int kMinCoverage = 3;

enum class CompressAlgo {
    Greedy,    // size-first, pin sources
    DepOrder,  // dependency-order + un-pin / retarget
    GreedyDfs, // DFS from best add=+1 hub, cumulative add to root
    TreeDp,    // preference-forest DP (KEEP vs COMPRESS per hub)
};

inline const char* algo_name(CompressAlgo a) {
    switch (a) {
        case CompressAlgo::Greedy:   return "greedy";
        case CompressAlgo::DepOrder: return "dep-order";
        case CompressAlgo::GreedyDfs: return "greedy-dfs";
        case CompressAlgo::TreeDp:   return "tree-dp";
    }
    return "?";
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
    void build(const Shape& shape, std::string text,
               int max_add = 8, CompressAlgo algo = CompressAlgo::Greedy) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        auto t0 = Clock::now();
        G_ = build_gapped_sa(shape, std::move(text));
        auto t1 = Clock::now();
        span_ = shape.span;
        max_add_ = std::max(1, max_add);
        algo_ = algo;
        compress_();
        auto t2 = Clock::now();
        if (timing) {
            auto ms = [](Clock::time_point a, Clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count();
            };
            std::fprintf(stderr,
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

private:
    GappedSA G_;
    int span_ = 1;
    int max_add_ = 8;
    CompressAlgo algo_ = CompressAlgo::Greedy;

    std::vector<int64_t> C_;
    std::vector<int64_t> rank_to_C_;
    std::unordered_map<uint64_t, HashEntry> table_;

    std::vector<int32_t> rank_of_;
    std::vector<uint8_t> removed_;
    std::vector<int32_t> pin_count_;   // #accepted candidates using this rank as source
    std::vector<uint64_t> pin_owner_;  // one owner name when pin_count==1 (undef if 0)
    size_t kept_count_ = 0;            // #ranks with removed_==0
    int last_pin_path_ = 0;            // 0=none 1=fast 2=multi (debug/timing)

    struct Interval { uint64_t name; int32_t lo, hi; };
    struct Candidate {
        int32_t target_lo = 0, target_hi = 0;
        uint64_t name = 0;
        int add = 0;
        int32_t src_lo = 0, src_hi = 0;
        std::vector<int32_t> covered;   // ranks of I_c, in source order
        int coverage() const { return (int)covered.size(); }
    };

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

    bool run_kept_(int32_t a, int32_t b) const {
        for (int32_t r = a; r < b; ++r) if (removed_[r]) return false;
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
                int32_t s_lo = pr[i].first, s_hi = pr[j-1].first + 1;
                bool disjoint = (s_hi <= lo || s_lo >= hi);
                int cov = (int)(j - i);
                if (disjoint && cov >= cov_floor) {
                    bool ok = true;
                    if (require_avail) {
                        ok = run_kept_(s_lo, s_hi);
                        for (size_t t = i; t < j && ok; ++t) {
                            int32_t r = pr[t].second;
                            if (removed_[r] || pin_count_[r] > 0) ok = false;
                        }
                    }
                    if (ok) {
                        Candidate c;
                        c.target_lo = lo; c.target_hi = hi; c.name = name;
                        c.add = add; c.src_lo = s_lo; c.src_hi = s_hi;
                        c.covered.reserve(cov);
                        for (size_t t = i; t < j; ++t) c.covered.push_back(pr[t].second);
                        out.push_back(std::move(c));
                    }
                }
                i = j;
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

    // First-symbol name of an SA rank.
    uint64_t name_of_rank_(int32_t r) const { return G_.first_symbol(r); }

    // ---- accept / revoke helpers ------------------------------------------
    void accept_(Candidate& c, std::unordered_map<uint64_t, Candidate>& accepted) {
        for (int32_t r = c.src_lo; r < c.src_hi; ++r) {
            if (pin_count_[r] == 0) pin_owner_[r] = c.name;
            ++pin_count_[r];
        }
        for (int32_t r : c.covered) {
            if (!removed_[r]) {
                removed_[r] = 1;
                --kept_count_;
            } else {
                removed_[r] = 1;
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
            } else if (pin_owner_[r] == name) {
                // Rare multi-pin: pick any remaining owner of r.
                pin_owner_[r] = 0;
                for (const auto& kv : accepted) {
                    if (kv.first == name) continue;
                    if (r >= kv.second.src_lo && r < kv.second.src_hi) {
                        pin_owner_[r] = kv.first;
                        break;
                    }
                }
            }
        }
        for (int32_t r : c.covered) {
            if (removed_[r]) {
                removed_[r] = 0;
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
            // Multi-pin: scan accepted once against the (usually small) pinned list.
            std::unordered_set<uint64_t> deps;
            for (auto& kv : accepted) {
                const Candidate& d = kv.second;
                for (int32_t r : pinned) {
                    if (r >= d.src_lo && r < d.src_hi) {
                        deps.insert(kv.first);
                        break;
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
    //   stop when next_dirty is empty (fixed point) or max_iters is hit.
    //
    // Env GCSA_PHASE2_MAX_ITERS=K (optional): hard-cap dirty generations.
    // Default is min(|I|+5, 32). A cap may leave |C| slightly larger.
    void run_phase2_(const std::vector<Interval>& intervals,
                     std::unordered_map<uint64_t, Candidate>& accepted,
                     const char* label,
                     bool trace,
                     bool timing,
                     std::unordered_map<uint64_t, std::vector<Candidate>>* precomputed = nullptr) {
        using Clock = std::chrono::steady_clock;
        const int N = (int)intervals.size();
        int max_iters = std::min(N + 5, 32);
        if (const char* env = std::getenv("GCSA_PHASE2_MAX_ITERS")) {
            int v = std::atoi(env);
            if (v > 0) max_iters = v;
        }
        std::fprintf(stderr, "[%s] Phase II: unpin/retarget...\n", label);
        if (timing && precomputed) {
            std::fprintf(stderr, "[timing] %s Phase II precomputed cache size=%zu\n",
                         label, precomputed->size());
        }
        auto t0 = Clock::now();

        // Size-order once; static candidate lists (avail ignored) cached once.
        std::vector<const Interval*> order;
        for (const auto& iv : intervals) if (iv.hi - iv.lo > 1) order.push_back(&iv);
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
            auto& slot = cand_cache[iv.name];
            slot = std::move(cands);
            return slot;
        };

        auto cov_of = [&](uint64_t name) -> int {
            auto it = accepted.find(name);
            return it == accepted.end() ? 0 : it->second.coverage();
        };

        // Seed: all |I|>1 (saturated / skip_best exit cheaply on first touch).
        std::vector<int> dirty(n_ord);
        for (int i = 0; i < n_ord; ++i) dirty[i] = i;
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
            in_next[oi] = 1;
            next_dirty.push_back(oi);
            ++phase2_dirty_marks;
        };

        int guard = 0;
        while (!dirty.empty() && guard++ < max_iters) {
            next_dirty.clear();
            std::fill(in_next.begin(), in_next.end(), 0);
            for (int oi : pin_blocked) in_pin_blocked[oi] = 0;
            pin_blocked.clear();
            bool improved_gen = false;
            std::vector<int> still_open; // best-cached > cur after this visit
            still_open.reserve(64);

            // order[] is size-desc; sorting dirty indices restores size order.
            std::sort(dirty.begin(), dirty.end());

            for (int oi : dirty) {
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
                for (const auto& c : cands) {
                    if (c.coverage() < kMinCoverage) continue;
                    if (c.coverage() <= cur_cov) break;
                    if (!run_kept_(c.src_lo, c.src_hi)) continue;
                    bool ok = true;
                    for (int32_t r : c.covered) if (removed_[r]) { ok = false; break; }
                    if (!ok) continue;

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
                        continue;
                    }
                    size_t after = kept_count_;
                    if (trace) {
                        if (after < before)
                            std::fprintf(stderr,
                                "II.%d IMPROVE I_%s src=I_%s[%d,%d) add=%d cov=%d  kept %zu->%zu\n",
                                guard, name_to_string(G_.shape, ivp->name).c_str(),
                                name_to_string(G_.shape, name_of_rank_(c.src_lo)).c_str(),
                                c.src_lo, c.src_hi, c.add, c.coverage(), before, after);
                        else
                            std::fprintf(stderr,
                                "II.%d REPLACE I_%s cov %d->%d (kept unchanged %zu)\n",
                                guard, name_to_string(G_.shape, ivp->name).c_str(),
                                cur_cov, c.coverage(), after);
                    }
                    improved_gen = true;
                    ++phase2_improve;
                    got_cov = c.coverage();
                    break;
                }
                // Still room vs cached best (pin-blocked / source not kept yet /
                // accepted a suboptimal cand) — retry next generation if anyone
                // improved (pins / kept sets may have moved).
                if (best_cov > got_cov) still_open.push_back(oi);
            }

            if (improved_gen) {
                for (int oi : pin_blocked) mark_next(oi);
                for (int oi : still_open) mark_next(oi);
            }
            dirty.swap(next_dirty);
        }
        auto t1 = Clock::now();
        std::fprintf(stderr, "[%s] Phase II done\n", label);
        if (trace) {
            std::fprintf(stderr, "after Phase II: kept=%zu accepted=%zu\n",
                         kept_count_, accepted.size());
        }
        if (timing) {
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            std::fprintf(stderr,
                "[timing] %s Phase II: %.1fms  iters=%d enums=%zu cache_hits=%zu "
                "skip_sat=%zu skip_best=%zu tries=%zu "
                "fail=%zu improve=%zu dirty_marks=%zu "
                "pin(none/fast/multi)=%zu/%zu/%zu "
                "enum=%.1fms try=%.1fms\n",
                label, ms, guard, phase2_enum, phase2_cache_hits,
                phase2_skip_sat, phase2_skip_best, phase2_tries,
                phase2_fail, phase2_improve, phase2_dirty_marks,
                phase2_no_pin, phase2_fast_pin, phase2_multi_pin,
                ms_enum, ms_try);
        }
    }

    // ---- algorithms -------------------------------------------------------
    void compress_greedy_(const std::vector<Interval>& intervals,
                          std::unordered_map<uint64_t, Candidate>& accepted) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        auto t0 = Clock::now();
        std::fprintf(stderr, "[greedy] size-order accept...\n");
        std::vector<const Interval*> order;
        for (const auto& iv : intervals) if (iv.hi - iv.lo > 1) order.push_back(&iv);
        std::sort(order.begin(), order.end(),
                  [](const Interval* a, const Interval* b){ return (a->hi-a->lo) > (b->hi-b->lo); });
        // Optional trace: set GCSA_TRACE_GREEDY=1
        const bool trace = (std::getenv("GCSA_TRACE_GREEDY") != nullptr);
        int step = 0;
        if (trace) {
            std::fprintf(stderr, "greedy order:");
            for (auto* ivp : order)
                std::fprintf(stderr, " %s(size=%d)",
                             name_to_string(G_.shape, ivp->name).c_str(), ivp->hi - ivp->lo);
            std::fprintf(stderr, "\n");
        }
        for (const Interval* ivp : order) {
            Candidate c = best_candidate_(ivp->name, ivp->lo, ivp->hi, /*avail=*/true);
            ++step;
            if (trace) {
                std::fprintf(stderr, "step %d: I_%s [%d,%d) size=%d",
                             step, name_to_string(G_.shape, ivp->name).c_str(),
                             ivp->lo, ivp->hi, ivp->hi - ivp->lo);
                if (c.coverage() < kMinCoverage)
                    std::fprintf(stderr, " -> SKIP\n");
                else
                    std::fprintf(stderr,
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
        std::fprintf(stderr, "[greedy] size-order accept done\n");
        if (timing) {
            double ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
            std::fprintf(stderr, "[timing] greedy: %.1fms  (#I>1=%zu accepted=%zu)\n",
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
        std::fprintf(stderr, "[dep-order] preference / DAG build...\n");
        auto t0 = Clock::now();
        for (const auto& iv : intervals) {
            if (iv.hi - iv.lo <= 1) continue;
            auto cands = enumerate_candidates_(iv.name, iv.lo, iv.hi, /*avail=*/false);
            if (iv.hi - iv.lo > 2) {
                Candidate c = pick_best_(cands);
                if (c.coverage() >= kMinCoverage) {
                    uint64_t src_name = name_of_rank_(c.src_lo);
                    prefs.push_back({&iv, std::move(c), src_name});
                }
            }
            std::sort(cands.begin(), cands.end(), [](const Candidate& a, const Candidate& b){
                if (a.coverage() != b.coverage()) return a.coverage() > b.coverage();
                return a.add > b.add;
            });
            cand_cache.emplace(iv.name, std::move(cands));
        }
        auto t1 = Clock::now();

        // DAG edge: src_name -> target_name (target depends on source).
        // Only introduce an edge when the preferred candidate has cov >= kMinCoverage
        // (same floor as Phase I/II acceptance; weaker prefs are not recorded).
        std::vector<std::vector<int>> outs(N), ins(N);
        std::vector<int> indeg(N, 0);
        auto pref_of = [&](uint64_t name) -> Pref* {
            for (auto& p : prefs) if (p.iv->name == name) return &p;
            return nullptr;
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
            std::fprintf(stderr, "=== preferred ===\n");
            for (auto& p : prefs)
                std::fprintf(stderr, "  %s -> prefers %s add=%d cov=%d src[%d,%d)\n",
                    name_to_string(G_.shape, p.iv->name).c_str(),
                    name_to_string(G_.shape, p.src_name).c_str(),
                    p.cand.add, p.cand.coverage(), p.cand.src_lo, p.cand.src_hi);
            std::fprintf(stderr, "=== DAG edges (src -> target, cov>=%d) ===\n", kMinCoverage);
            for (auto& p : prefs) {
                if (p.cand.coverage() < kMinCoverage) continue;
                if (p.iv->name == p.src_name) continue;
                std::fprintf(stderr, "  %s -> %s\n",
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
        std::fprintf(stderr, "[dep-order] preference / DAG done\n");

        if (trace) {
            std::fprintf(stderr, "=== Phase I accept order (reverse topo, sinks first) ===\n  ");
            for (int k = N - 1; k >= 0; --k) {
                const auto& iv = intervals[topo[k]];
                if (iv.hi - iv.lo <= 2) continue;
                if (!pref_of(iv.name)) continue;
                std::fprintf(stderr, "%s ", name_to_string(G_.shape, iv.name).c_str());
            }
            std::fprintf(stderr, "\n");
        }

        // Pass 1: accept preferred candidates in reverse topo (sinks first),
        // with availability. This pins sources that dependents need.
        std::fprintf(stderr, "[dep-order] Phase I: reverse-topo accept...\n");
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
                std::fprintf(stderr, "I.%d I_%s", step, name_to_string(G_.shape, name).c_str());
                if (c.coverage() < kMinCoverage) std::fprintf(stderr, " -> SKIP\n");
                else std::fprintf(stderr, " -> ACCEPT src=I_%s[%d,%d) add=%d cov=%d\n",
                    name_to_string(G_.shape, name_of_rank_(c.src_lo)).c_str(),
                    c.src_lo, c.src_hi, c.add, c.coverage());
            }
            if (c.coverage() < kMinCoverage) continue;
            accept_(c, accepted);
        }
        auto t3 = Clock::now();
        std::fprintf(stderr, "[dep-order] Phase I done\n");
        if (trace) {
            size_t kept = 0; for (auto x : removed_) if (!x) ++kept;
            std::fprintf(stderr, "after Phase I: kept=%zu accepted=%zu\n", kept, accepted.size());
        }

        run_phase2_(intervals, accepted, "dep-order", trace, timing, &cand_cache);
        auto t4 = Clock::now();
        if (timing) {
            auto ms = [](Clock::time_point a, Clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count();
            };
            std::fprintf(stderr,
                "[timing] dep-order prefs=%.1fms phaseI=%.1fms phaseII=%.1fms total=%.1fms\n"
                "         N=%d prefs=%zu\n",
                ms(t0, t1), ms(t2, t3), ms(t3, t4), ms(t_all, t4),
                N, prefs.size());
        }
    }

    // Best add=1 candidate for target whose source lies entirely inside source_iv.
    Candidate best_add1_from_(const Interval& target, const Interval& source,
                              bool require_avail) const {
        Candidate best;
        // Only enumerate add=1 — callers never accept other adds here.
        auto cands = enumerate_candidates_(target.name, target.lo, target.hi,
                                           require_avail, /*only_add=*/1);
        for (auto& c : cands) {
            if (c.src_lo < source.lo || c.src_hi > source.hi) continue;
            if (c.coverage() > best.coverage()
                || (c.coverage() == best.coverage() && c.src_lo < best.src_lo))
                best = c;
        }
        return best;
    }

    // Greedy-DFS: root = interval maximizing total add=+1 outbound coverage to
    // targets with |I|>2.  DFS along add=+1 links; every reached target points
    // at the DFS root with cumulative add (depth).  Repeat on leftovers, then
    // a final availability-aware greedy sweep.
    void compress_greedy_dfs_(const std::vector<Interval>& intervals,
                              std::unordered_map<uint64_t, Candidate>& accepted) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        auto t_all = Clock::now();
        double score_ms = 0, grow_ms = 0;
        size_t score_calls = 0, trees = 0, enum_calls = 0;

        const bool trace = (std::getenv("GCSA_TRACE_DFS") != nullptr);
        std::unordered_map<uint64_t, const Interval*> by_name;
        for (const auto& iv : intervals) by_name[iv.name] = &iv;

        auto isize = [](const Interval& iv) { return iv.hi - iv.lo; };

        // Outbound add=+1 score for a prospective root.
        auto root_score = [&](const Interval& R) -> int {
            auto ts = Clock::now();
            int score = 0;
            for (const auto& T : intervals) {
                if (T.name == R.name || isize(T) <= 2) continue;
                if (accepted.count(T.name)) continue;
                ++enum_calls;
                Candidate c = best_add1_from_(T, R, /*avail=*/false);
                if (c.coverage() >= kMinCoverage) score += c.coverage();
            }
            score_ms += std::chrono::duration<double, std::milli>(Clock::now() - ts).count();
            ++score_calls;
            return score;
        };

        // Compose child --add1--> parent into child --add'--> root via parent's
        // accepted-to-root candidate. Returns empty (cov 0) on failure.
        auto compose_to_root = [&](const Candidate& link_add1,
                                   const Candidate& parent_to_root) -> Candidate {
            Candidate out;
            int32_t nlo = 0, nhi = 0;
            // Treat link_add1 as a "dependent" whose source is in parent;
            // parent_to_root.covered lists the parent ranks recovered from root.
            if (!can_retarget_(link_add1, parent_to_root, nlo, nhi)) return out;
            out = link_add1;
            out.src_lo = nlo;
            out.src_hi = nhi;
            out.add = parent_to_root.add + 1;
            if (out.add > max_add_) { out.covered.clear(); return out; }
            return out;
        };

        std::unordered_set<uint64_t> is_root;  // roots stay stored (no offset)

        // Grow as many DFS trees as profitable.
        std::fprintf(stderr, "[greedy-dfs] DFS trees (score+grow)...\n");
        while (true) {
            const Interval* root = nullptr;
            int best = 0;
            for (const auto& R : intervals) {
                if (isize(R) <= 2) continue;
                if (accepted.count(R.name) || is_root.count(R.name)) continue;
                // Root ranks must still be kept.
                bool kept = true;
                for (int32_t r = R.lo; r < R.hi; ++r) if (removed_[r]) { kept = false; break; }
                if (!kept) continue;
                int sc = root_score(R);
                if (sc > best) { best = sc; root = &R; }
            }
            if (!root || best <= 0) break;

            auto tg0 = Clock::now();
            ++trees;
            is_root.insert(root->name);
            if (trace)
                std::fprintf(stderr, "DFS root=I_%s size=%d outbound_cov=%d\n",
                    name_to_string(G_.shape, root->name).c_str(), isize(*root), best);

            // Identity map for the root: covered = all its ranks, add=0, src=itself.
            Candidate root_id;
            root_id.name = root->name;
            root_id.target_lo = root->lo;
            root_id.target_hi = root->hi;
            root_id.add = 0;
            root_id.src_lo = root->lo;
            root_id.src_hi = root->hi;
            root_id.covered.clear();
            for (int32_t r = root->lo; r < root->hi; ++r) root_id.covered.push_back(r);

            // parent_to_root[name] for nodes in the tree (including root).
            std::unordered_map<uint64_t, Candidate> to_root;
            to_root[root->name] = root_id;

            struct Frame { uint64_t name; };
            std::vector<Frame> stack;

            // Seed stack with add=+1 children of the root (|I|>2).
            {
                std::vector<Candidate> kids;
                for (const auto& T : intervals) {
                    if (T.name == root->name || isize(T) <= 2) continue;
                    if (accepted.count(T.name)) continue;
                    ++enum_calls;
                    Candidate link = best_add1_from_(T, *root, /*avail=*/true);
                    if (link.coverage() < kMinCoverage) continue;
                    kids.push_back(std::move(link));
                }
                std::sort(kids.begin(), kids.end(), [](const Candidate& a, const Candidate& b){
                    return a.coverage() > b.coverage();
                });
                for (auto& link : kids) {
                    // Directly to root with add=1 (source already inside root).
                    if (!run_kept_(link.src_lo, link.src_hi)) continue;
                    bool ok = true;
                    for (int32_t r : link.covered) if (removed_[r] || pin_count_[r] > 0) { ok = false; break; }
                    if (!ok) continue;
                    accept_(link, accepted);
                    to_root[link.name] = link;
                    stack.push_back({link.name});
                    if (trace)
                        std::fprintf(stderr, "  depth1 I_%s <- root add=1 cov=%d\n",
                            name_to_string(G_.shape, link.name).c_str(), link.coverage());
                }
            }

            // DFS: from each node, attach add=+1 children pointing at root with add+1.
            while (!stack.empty()) {
                Frame fr = stack.back();
                stack.pop_back();
                const Interval* P = by_name[fr.name];
                auto itP = to_root.find(fr.name);
                if (!P || itP == to_root.end()) continue;
                const Candidate& parent_tr = itP->second;

                std::vector<Candidate> kids;
                for (const auto& U : intervals) {
                    if (U.name == fr.name || U.name == root->name) continue;
                    if (isize(U) <= 2) continue;
                    if (accepted.count(U.name) || is_root.count(U.name)) continue;
                    // Link U <- P with add=1 (ignore avail; compose handles pins).
                    ++enum_calls;
                    Candidate link = best_add1_from_(U, *P, /*avail=*/false);
                    if (link.coverage() < kMinCoverage) continue;
                    kids.push_back(std::move(link));
                }
                std::sort(kids.begin(), kids.end(), [](const Candidate& a, const Candidate& b){
                    return a.coverage() > b.coverage();
                });

                for (auto& link : kids) {
                    Candidate composed = compose_to_root(link, parent_tr);
                    if (composed.coverage() < kMinCoverage) continue;
                    // Availability on composed source (subset of root) and covered.
                    if (!run_kept_(composed.src_lo, composed.src_hi)) continue;
                    bool ok = true;
                    for (int32_t r : composed.covered)
                        if (removed_[r] || pin_count_[r] > 0) { ok = false; break; }
                    if (!ok) continue;
                    accept_(composed, accepted);
                    to_root[composed.name] = composed;
                    stack.push_back({composed.name});
                    if (trace)
                        std::fprintf(stderr, "  depth%d I_%s <- root add=%d cov=%d (via I_%s)\n",
                            composed.add,
                            name_to_string(G_.shape, composed.name).c_str(),
                            composed.add, composed.coverage(),
                            name_to_string(G_.shape, fr.name).c_str());
                }
            }
            grow_ms += std::chrono::duration<double, std::milli>(Clock::now() - tg0).count();
        }
        std::fprintf(stderr, "[greedy-dfs] DFS trees done\n");

        // Leftover sweep: classic size-greedy on remaining intervals.
        std::fprintf(stderr, "[greedy-dfs] leftover greedy sweep...\n");
        auto tl0 = Clock::now();
        std::vector<const Interval*> order;
        for (const auto& iv : intervals)
            if (isize(iv) > 1 && !accepted.count(iv.name)) order.push_back(&iv);
        std::sort(order.begin(), order.end(),
                  [](const Interval* a, const Interval* b){ return (a->hi-a->lo) > (b->hi-b->lo); });
        for (const Interval* ivp : order) {
            Candidate c = best_candidate_(ivp->name, ivp->lo, ivp->hi, /*avail=*/true);
            if (c.coverage() < kMinCoverage) continue;
            accept_(c, accepted);
            if (trace)
                std::fprintf(stderr, "leftover I_%s <- I_%s add=%d cov=%d\n",
                    name_to_string(G_.shape, ivp->name).c_str(),
                    name_to_string(G_.shape, name_of_rank_(c.src_lo)).c_str(),
                    c.add, c.coverage());
        }
        auto tl1 = Clock::now();
        std::fprintf(stderr, "[greedy-dfs] leftover sweep done\n");
        if (trace) {
            size_t kept = 0; for (auto x : removed_) if (!x) ++kept;
            std::fprintf(stderr, "greedy-dfs done: kept=%zu accepted=%zu roots=%zu\n",
                         kept, accepted.size(), is_root.size());
        }
        if (timing) {
            double leftover_ms = std::chrono::duration<double, std::milli>(tl1 - tl0).count();
            double total_ms = std::chrono::duration<double, std::milli>(tl1 - t_all).count();
            std::fprintf(stderr,
                "[timing] greedy-dfs score=%.1fms grow=%.1fms leftover=%.1fms total=%.1fms\n"
                "         trees=%zu score_calls=%zu best_add1_enums=%zu #I=%zu\n",
                score_ms, grow_ms, leftover_ms, total_ms,
                trees, score_calls, enum_calls, intervals.size());
        }
    }

    // Preference-forest DP: each strong preference edge (cov>=kMinCoverage, |I|>2) makes
    // the target a child of its preferred source.  For every node we choose
    // KEEP (all ranks stored, children may compress against us) vs COMPRESS
    // (drop cov ranks via the preferred candidate; children lose that source).
    // Exact for |C| on this forest model; leftover names get an avail. greedy.
    void compress_tree_dp_(const std::vector<Interval>& intervals,
                           std::unordered_map<uint64_t, Candidate>& accepted) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        const bool trace = (std::getenv("GCSA_TRACE_DP") != nullptr);
        auto t_all = Clock::now();

        std::unordered_map<uint64_t, const Interval*> by_name;
        for (const auto& iv : intervals) by_name[iv.name] = &iv;

        struct Pref { const Interval* iv; Candidate cand; uint64_t src_name; };
        std::unordered_map<uint64_t, Pref> prefs;   // target -> preferred
        // Static (avail=false) candidate lists — reused by Phase II.
        // Prefs use pick_best_ on enum order (stable tie-break); Phase II wants
        // coverage-desc order, so we sort only after recording the preferred.
        std::unordered_map<uint64_t, std::vector<Candidate>> cand_cache;

        std::fprintf(stderr, "[tree-dp] preference forest build...\n");
        for (const auto& iv : intervals) {
            if (iv.hi - iv.lo <= 1) continue;
            auto cands = enumerate_candidates_(iv.name, iv.lo, iv.hi, /*avail=*/false);
            if (iv.hi - iv.lo > 2) {
                Candidate best = pick_best_(cands);
                if (best.coverage() >= kMinCoverage) {
                    uint64_t src = name_of_rank_(best.src_lo);
                    prefs.emplace(iv.name, Pref{&iv, std::move(best), src});
                }
            }
            std::sort(cands.begin(), cands.end(), [](const Candidate& a, const Candidate& b){
                if (a.coverage() != b.coverage()) return a.coverage() > b.coverage();
                return a.add > b.add;
            });
            cand_cache.emplace(iv.name, std::move(cands));
        }

        // Forest edges: only cov>=kMinCoverage (same threshold as DepOrder DAG edges).
        std::unordered_map<uint64_t, std::vector<uint64_t>> children;
        std::unordered_map<uint64_t, uint64_t> parent;
        auto reaches = [&](uint64_t from, uint64_t to) {
            std::unordered_set<uint64_t> seen;
            std::vector<uint64_t> st = {from};
            while (!st.empty()) {
                uint64_t u = st.back(); st.pop_back();
                if (!seen.insert(u).second) continue;
                if (u == to) return true;
                auto it = children.find(u);
                if (it == children.end()) continue;
                for (uint64_t v : it->second) st.push_back(v);
            }
            return false;
        };

        for (auto& kv : prefs) {
            uint64_t t = kv.first;
            uint64_t s = kv.second.src_name;
            if (kv.second.cand.coverage() < kMinCoverage) continue;
            if (s == t) continue;
            if (!by_name.count(s)) continue;
            // Skip edge if it would create a cycle in the forest.
            if (reaches(t, s)) continue;
            parent[t] = s;
            children[s].push_back(t);
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

        if (trace) {
            std::fprintf(stderr, "=== tree-dp forest (src -> dependents) ===\n");
            for (auto& kv : children) {
                std::fprintf(stderr, "  %s ->",
                             name_to_string(G_.shape, kv.first).c_str());
                for (uint64_t c : kv.second)
                    std::fprintf(stderr, " %s",
                                 name_to_string(G_.shape, c).c_str());
                std::fprintf(stderr, "\n");
            }
            std::fprintf(stderr, "roots:");
            for (uint64_t r : roots)
                std::fprintf(stderr, " %s", name_to_string(G_.shape, r).c_str());
            std::fprintf(stderr, "\n");
        }

        auto isize = [&](uint64_t n) -> int {
            auto it = by_name.find(n);
            if (it == by_name.end()) return 0;
            return it->second->hi - it->second->lo;
        };

        // Memo: (node, source_available_for_this_node) -> (cost, compress?)
        struct Cell { int cost; bool compress; };
        std::map<std::pair<uint64_t, int>, Cell> memo;

        std::function<Cell(uint64_t, bool)> solve = [&](uint64_t v, bool src_ok) -> Cell {
            auto key = std::make_pair(v, src_ok ? 1 : 0);
            auto it = memo.find(key);
            if (it != memo.end()) return it->second;

            const auto& ch = children[v];
            // KEEP v: all ranks stored; children may compress against v.
            int cost_keep = isize(v);
            for (uint64_t w : ch) cost_keep += solve(w, /*src_ok=*/true).cost;

            // COMPRESS v via preferred candidate, only if that source is available.
            int cost_comp = std::numeric_limits<int>::max();
            bool can_comp = false;
            if (src_ok) {
                auto pit = prefs.find(v);
                if (pit != prefs.end() && pit->second.cand.coverage() >= kMinCoverage) {
                    can_comp = true;
                    int cc = isize(v) - pit->second.cand.coverage();
                    for (uint64_t w : ch) cc += solve(w, /*src_ok=*/false).cost;
                    cost_comp = cc;
                }
            }

            Cell cell;
            if (can_comp && cost_comp < cost_keep) {
                cell = {cost_comp, true};
            } else {
                cell = {cost_keep, false};
            }
            memo[key] = cell;
            return cell;
        };

        // Chosen compressions (target name -> candidate).
        std::unordered_map<uint64_t, Candidate> chosen;

        std::function<void(uint64_t, bool)> apply = [&](uint64_t v, bool src_ok) {
            Cell cell = solve(v, src_ok);
            if (cell.compress) {
                chosen[v] = prefs[v].cand;
                for (uint64_t w : children[v]) apply(w, /*src_ok=*/false);
            } else {
                for (uint64_t w : children[v]) apply(w, /*src_ok=*/true);
            }
        };

        std::fprintf(stderr, "[tree-dp] DP on forest...\n");
        for (uint64_t r : roots) {
            // Root: KEEP, or COMPRESS if it has a preferred source outside its subtree.
            int cost_keep = isize(r);
            for (uint64_t w : children[r]) cost_keep += solve(w, true).cost;

            int cost_comp = std::numeric_limits<int>::max();
            bool can_comp = false;
            auto pit = prefs.find(r);
            if (pit != prefs.end() && pit->second.cand.coverage() >= kMinCoverage) {
                uint64_t s = pit->second.src_name;
                // Source must not lie in this preference subtree.
                if (s != r && !reaches(r, s)) {
                    can_comp = true;
                    int cc = isize(r) - pit->second.cand.coverage();
                    for (uint64_t w : children[r]) cc += solve(w, false).cost;
                    cost_comp = cc;
                }
            }

            if (trace) {
                std::fprintf(stderr, "root %s: keep=%d comp=%s\n",
                             name_to_string(G_.shape, r).c_str(), cost_keep,
                             can_comp ? std::to_string(cost_comp).c_str() : "n/a");
            }

            if (can_comp && cost_comp < cost_keep) {
                chosen[r] = pit->second.cand;
                for (uint64_t w : children[r]) apply(w, false);
            } else {
                for (uint64_t w : children[r]) apply(w, true);
            }
        }
        std::fprintf(stderr, "[tree-dp] DP done\n");

        // Materialize: accept deepest dependents first (pins hubs before hubs decide).
        std::fprintf(stderr, "[tree-dp] accept chosen...\n");
        std::function<void(uint64_t)> accept_subtree = [&](uint64_t v) {
            for (uint64_t w : children[v]) accept_subtree(w);
            auto it = chosen.find(v);
            if (it == chosen.end()) return;
            // Re-validate under availability; fall back to best available.
            Candidate c = it->second;
            bool ok = run_kept_(c.src_lo, c.src_hi);
            for (int32_t rr : c.covered)
                if (removed_[rr] || pin_count_[rr] > 0) ok = false;
            if (!ok) {
                const Interval* iv = by_name[v];
                c = best_candidate_(iv->name, iv->lo, iv->hi, /*avail=*/true);
            }
            if (c.coverage() >= kMinCoverage) {
                if (trace) {
                    std::fprintf(stderr, "  ACCEPT %s <- %s add=%d cov=%d\n",
                        name_to_string(G_.shape, v).c_str(),
                        name_to_string(G_.shape, name_of_rank_(c.src_lo)).c_str(),
                        c.add, c.coverage());
                }
                accept_(c, accepted);
            }
        };
        for (uint64_t r : roots) accept_subtree(r);

        // Leftover intervals (not decided / not compressed): greedy with availability.
        std::fprintf(stderr, "[tree-dp] leftover greedy...\n");
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
            if (c.coverage() < kMinCoverage) continue;
            accept_(c, accepted);
        }
        std::fprintf(stderr, "[tree-dp] leftover done\n");

        // Phase II: same unpin/retarget fixed point as DepOrder (reuse cand_cache).
        run_phase2_(intervals, accepted, "tree-dp", trace, timing, &cand_cache);

        if (timing) {
            double ms = std::chrono::duration<double, std::milli>(Clock::now() - t_all).count();
            std::fprintf(stderr,
                "[timing] tree-dp total: %.1fms  roots=%zu forest_nodes=%zu accepted=%zu kept=%zu\n",
                ms, roots.size(), nodes.size(), accepted.size(), kept_count_);
        }
    }

    void compress_() {
        const size_t m = G_.m();
        rank_of_.assign(m, 0);
        for (size_t r = 0; r < m; ++r) rank_of_[G_.sa[r]] = (int32_t)r;
        removed_.assign(m, 0);
        pin_count_.assign(m, 0);
        pin_owner_.assign(m, 0);
        kept_count_ = m;
        table_.clear();
        C_.clear();

        std::fprintf(stderr, "[%s] compressing...\n", algo_name(algo_));

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
        else if (algo_ == CompressAlgo::GreedyDfs)
            compress_greedy_dfs_(intervals, accepted);
        else if (algo_ == CompressAlgo::TreeDp)
            compress_tree_dp_(intervals, accepted);
        else
            compress_greedy_(intervals, accepted);

        // Build C and hash table.
        std::fprintf(stderr, "[%s] finalize: build C / hash table...\n", algo_name(algo_));
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
        std::fprintf(stderr, "[%s] finalize done\n", algo_name(algo_));
    }
};

} // namespace gcsa
