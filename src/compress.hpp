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
// Two source-selection strategies are implemented (CompressAlgo):
//   Greedy   – size-order availability-aware greedy (v1); pins sources forever.
//   DepOrder – compress sink k-mers first; then compress former sources via a
//              deeper LCP-interval, retargeting dependents' (pos,add) through
//              the new source (transitive add). Iterates to a fixed point.

#pragma once

#include "gapped_sa.hpp"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <string>

namespace gcsa {

enum class CompressAlgo {
    Greedy,    // size-first, pin sources
    DepOrder,  // dependency-order + un-pin / retarget
};

inline const char* algo_name(CompressAlgo a) {
    switch (a) {
        case CompressAlgo::Greedy:   return "greedy";
        case CompressAlgo::DepOrder: return "dep-order";
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
        G_ = build_gapped_sa(shape, std::move(text));
        span_ = shape.span;
        max_add_ = std::max(1, max_add);
        algo_ = algo;
        compress_();
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

    struct Interval { uint64_t name; int32_t lo, hi; };
    struct Candidate {
        int32_t target_lo = 0, target_hi = 0;
        uint64_t name = 0;
        int add = 0;
        int32_t src_lo = 0, src_hi = 0;
        std::vector<int32_t> covered;   // ranks of I_c, in source order
        int coverage() const { return (int)covered.size(); }
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
    std::vector<Candidate> enumerate_candidates_(uint64_t name, int32_t lo, int32_t hi,
                                                  bool require_avail) const {
        std::vector<Candidate> out;
        for (int add = 1; add <= max_add_; ++add) {
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
                if (disjoint && cov >= 2) {
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
        for (int32_t r = c.src_lo; r < c.src_hi; ++r) ++pin_count_[r];
        for (int32_t r : c.covered) removed_[r] = 1;
        accepted[c.name] = c;
    }

    void revoke_(uint64_t name, std::unordered_map<uint64_t, Candidate>& accepted) {
        auto it = accepted.find(name);
        if (it == accepted.end()) return;
        Candidate& c = it->second;
        for (int32_t r = c.src_lo; r < c.src_hi; ++r) --pin_count_[r];
        for (int32_t r : c.covered) removed_[r] = 0;
        accepted.erase(it);
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
    bool try_accept_with_retarget_(Candidate cand,
                                   std::unordered_map<uint64_t, Candidate>& accepted,
                                   const std::vector<Interval>& /*intervals*/) {
        // Covered ranks must not already be removed.
        for (int32_t r : cand.covered) if (removed_[r]) return false;
        // Source must be kept.
        if (!run_kept_(cand.src_lo, cand.src_hi)) return false;
        // Source must be disjoint from target (already true from enum).
        // Find dependents whose source intersects cand.covered.
        std::vector<uint64_t> to_retarget;
        for (auto& kv : accepted) {
            const Candidate& d = kv.second;
            // Intersection of [d.src_lo,d.src_hi) with cand.covered (as a set of ranks)?
            // Fast path: any pin in cand.covered means some dependent uses it.
            bool hits = false;
            for (int32_t r : cand.covered) if (pin_count_[r] > 0) {
                // Check if this dependent's source contains r.
                if (r >= d.src_lo && r < d.src_hi) { hits = true; break; }
            }
            if (!hits) continue;
            int32_t nlo, nhi;
            if (!can_retarget_(d, cand, nlo, nhi)) return false;
            to_retarget.push_back(kv.first);
        }

        // Also: covered ranks that are pinned must ALL belong to retargetable deps
        // (no leftover pin from a dep we can't retarget — already enforced above).
        for (int32_t r : cand.covered) {
            if (pin_count_[r] == 0) continue;
            // Every pin on r must be from a dep in to_retarget.
            // pin_count only tells cardinality; verify each accepted dep that
            // uses r is in to_retarget.
            for (auto& kv : accepted) {
                const Candidate& d = kv.second;
                if (r >= d.src_lo && r < d.src_hi) {
                    if (std::find(to_retarget.begin(), to_retarget.end(), kv.first)
                        == to_retarget.end())
                        return false;
                }
            }
        }

        // If this name already has an offset, revoke it first (we'll replace).
        if (accepted.count(cand.name)) revoke_(cand.name, accepted);

        // Retarget dependents: revoke + re-accept with updated src/add.
        std::vector<Candidate> retargeted;
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

        accept_(cand, accepted);
        for (auto& d : retargeted) accept_(d, accepted);
        return true;
    }

    // ---- algorithms -------------------------------------------------------
    void compress_greedy_(const std::vector<Interval>& intervals,
                          std::unordered_map<uint64_t, Candidate>& accepted) {
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
                if (c.coverage() < 2)
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
            if (c.coverage() < 2) continue;
            accept_(c, accepted);
        }
    }

    void compress_dep_order_(const std::vector<Interval>& intervals,
                             std::unordered_map<uint64_t, Candidate>& accepted) {
        // Preferred candidate per interval (ignore availability — global best).
        struct Pref { const Interval* iv; Candidate cand; uint64_t src_name; };
        std::vector<Pref> prefs;
        std::unordered_map<uint64_t, int> id_of;   // k-mer name -> node id
        for (size_t i = 0; i < intervals.size(); ++i)
            id_of[intervals[i].name] = (int)i;
        const int N = (int)intervals.size();

        for (const auto& iv : intervals) {
            if (iv.hi - iv.lo <= 1) continue;
            Candidate c = best_candidate_(iv.name, iv.lo, iv.hi, /*avail=*/false);
            if (c.coverage() < 2) continue;
            uint64_t src_name = name_of_rank_(c.src_lo);
            prefs.push_back({&iv, std::move(c), src_name});
        }

        // DAG edge: src_name -> target_name (target depends on source).
        std::vector<std::vector<int>> outs(N), ins(N);
        std::vector<int> indeg(N, 0);
        auto pref_of = [&](uint64_t name) -> Pref* {
            for (auto& p : prefs) if (p.iv->name == name) return &p;
            return nullptr;
        };
        for (auto& p : prefs) {
            int t = id_of[p.iv->name];
            auto it = id_of.find(p.src_name);
            if (it == id_of.end()) continue;
            int s = it->second;
            if (s == t) continue;                  // self-loop (shouldn't after disjoint)
            outs[s].push_back(t);
            ins[t].push_back(s);
            ++indeg[t];
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

        // Pass 1: accept preferred candidates in reverse topo (sinks first),
        // with availability. This pins sources that dependents need.
        for (int k = N - 1; k >= 0; --k) {
            uint64_t name = intervals[topo[k]].name;
            Pref* p = pref_of(name);
            if (!p) continue;
            // Recompute under current availability (source may already be gone).
            Candidate c = best_candidate_(p->iv->name, p->iv->lo, p->iv->hi, /*avail=*/true);
            if (c.coverage() < 2) continue;
            accept_(c, accepted);
        }

        // Pass 2: fixed-point un-pin — try to compress intervals further
        // (including former sources) via deeper candidates, retargeting deps.
        bool changed = true;
        int guard = 0;
        while (changed && guard++ < N + 5) {
            changed = false;
            // Try larger intervals first for bigger wins.
            std::vector<const Interval*> order;
            for (const auto& iv : intervals) if (iv.hi - iv.lo > 1) order.push_back(&iv);
            std::sort(order.begin(), order.end(),
                      [](const Interval* a, const Interval* b){
                          return (a->hi-a->lo) > (b->hi-b->lo);
                      });
            for (const Interval* ivp : order) {
                // Enumerate available candidates that may require retargeting:
                // temporarily ignore pin_count on covered (enumerate without
                // avail), then filter with try_accept_with_retarget_.
                auto cands = enumerate_candidates_(ivp->name, ivp->lo, ivp->hi,
                                                   /*require_avail=*/false);
                // Prefer higher coverage, then deeper add.
                std::sort(cands.begin(), cands.end(), [](const Candidate& a, const Candidate& b){
                    if (a.coverage() != b.coverage()) return a.coverage() > b.coverage();
                    return a.add > b.add;
                });
                int cur_cov = 0;
                auto it = accepted.find(ivp->name);
                if (it != accepted.end()) cur_cov = it->second.coverage();

                for (auto& c : cands) {
                    if (c.coverage() <= cur_cov) break;   // sorted; no better left
                    // Source must currently be kept (not removed).
                    if (!run_kept_(c.src_lo, c.src_hi)) continue;
                    // Covered must not include already-removed ranks.
                    bool ok = true;
                    for (int32_t r : c.covered) if (removed_[r]) { ok = false; break; }
                    if (!ok) continue;
                    // Snapshot to allow revert on failure — try_accept mutates.
                    // (try_accept_with_retarget_ is transactional on failure.)
                    // Count kept positions before.
                    size_t before = 0;
                    for (uint8_t x : removed_) if (!x) ++before;

                    // Need a deep copy of accepted for rollback.
                    auto acc_snap = accepted;
                    auto rem_snap = removed_;
                    auto pin_snap = pin_count_;

                    if (!try_accept_with_retarget_(c, accepted, intervals)) {
                        accepted = std::move(acc_snap);
                        removed_ = std::move(rem_snap);
                        pin_count_ = std::move(pin_snap);
                        continue;
                    }
                    size_t after = 0;
                    for (uint8_t x : removed_) if (!x) ++after;
                    if (after < before) {
                        changed = true;
                        break;   // re-sort / rescan from largest
                    }
                    // No improvement (e.g. replaced with equal) — keep if more
                    // coverage on this name even if global kept count same.
                    if (c.coverage() > cur_cov) { changed = true; break; }
                    // Roll back if no gain.
                    accepted = std::move(acc_snap);
                    removed_ = std::move(rem_snap);
                    pin_count_ = std::move(pin_snap);
                }
            }
        }
    }

    void compress_() {
        const size_t m = G_.m();
        rank_of_.assign(m, 0);
        for (size_t r = 0; r < m; ++r) rank_of_[G_.sa[r]] = (int32_t)r;
        removed_.assign(m, 0);
        pin_count_.assign(m, 0);
        table_.clear();
        C_.clear();

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
        else
            compress_greedy_(intervals, accepted);

        // Build C and hash table.
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
