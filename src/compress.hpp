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

namespace gcsa {

enum class CompressAlgo {
    Greedy,    // size-first, pin sources
    DepOrder,  // dependency-order + un-pin / retarget
    GreedyDfs, // DFS from best add=+1 hub, cumulative add to root
};

inline const char* algo_name(CompressAlgo a) {
    switch (a) {
        case CompressAlgo::Greedy:   return "greedy";
        case CompressAlgo::DepOrder: return "dep-order";
        case CompressAlgo::GreedyDfs: return "greedy-dfs";
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

    // ---- algorithms -------------------------------------------------------
    void compress_greedy_(const std::vector<Interval>& intervals,
                          std::unordered_map<uint64_t, Candidate>& accepted) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        auto t0 = Clock::now();
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

        auto t0 = Clock::now();
        for (const auto& iv : intervals) {
            if (iv.hi - iv.lo <= 1) continue;
            Candidate c = best_candidate_(iv.name, iv.lo, iv.hi, /*avail=*/false);
            if (c.coverage() < 2) continue;
            uint64_t src_name = name_of_rank_(c.src_lo);
            prefs.push_back({&iv, std::move(c), src_name});
        }
        auto t1 = Clock::now();

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

        const bool trace = (std::getenv("GCSA_TRACE_DEP") != nullptr);
        if (trace) {
            std::fprintf(stderr, "=== preferred ===\n");
            for (auto& p : prefs)
                std::fprintf(stderr, "  %s -> prefers %s add=%d cov=%d src[%d,%d)\n",
                    name_to_string(G_.shape, p.iv->name).c_str(),
                    name_to_string(G_.shape, p.src_name).c_str(),
                    p.cand.add, p.cand.coverage(), p.cand.src_lo, p.cand.src_hi);
            std::fprintf(stderr, "=== DAG edges (src -> target) ===\n");
            for (auto& p : prefs) {
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

        if (trace) {
            std::fprintf(stderr, "=== Phase I accept order (reverse topo, sinks first) ===\n  ");
            for (int k = N - 1; k >= 0; --k) {
                const auto& iv = intervals[topo[k]];
                if (iv.hi - iv.lo <= 1) continue;
                if (!pref_of(iv.name)) continue;
                std::fprintf(stderr, "%s ", name_to_string(G_.shape, iv.name).c_str());
            }
            std::fprintf(stderr, "\n");
        }

        // Pass 1: accept preferred candidates in reverse topo (sinks first),
        // with availability. This pins sources that dependents need.
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
                if (c.coverage() < 2) std::fprintf(stderr, " -> SKIP\n");
                else std::fprintf(stderr, " -> ACCEPT src=I_%s[%d,%d) add=%d cov=%d\n",
                    name_to_string(G_.shape, name_of_rank_(c.src_lo)).c_str(),
                    c.src_lo, c.src_hi, c.add, c.coverage());
            }
            if (c.coverage() < 2) continue;
            accept_(c, accepted);
        }
        auto t3 = Clock::now();
        if (trace) {
            size_t kept = 0; for (auto x : removed_) if (!x) ++kept;
            std::fprintf(stderr, "after Phase I: kept=%zu accepted=%zu\n", kept, accepted.size());
        }

        // Pass 2: fixed-point un-pin — try to compress intervals further
        // (including former sources) via deeper candidates, retargeting deps.
        bool changed = true;
        int guard = 0;
        size_t phase2_enum = 0, phase2_tries = 0;
        size_t phase2_fail = 0, phase2_improve = 0;
        size_t phase2_fast_pin = 0, phase2_multi_pin = 0, phase2_no_pin = 0;
        double ms_enum = 0, ms_try = 0;
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
                auto te0 = Clock::now();
                auto cands = enumerate_candidates_(ivp->name, ivp->lo, ivp->hi,
                                                   /*require_avail=*/false);
                if (timing) ms_enum += std::chrono::duration<double, std::milli>(Clock::now() - te0).count();
                ++phase2_enum;
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

                    // try_accept is transactional on failure (no mutation).
                    // Successful tries always keep: we only consider cov > cur_cov,
                    // so the old O(n) snapshot / surgical-undo rollback path is dead.
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
                    changed = true;
                    ++phase2_improve;
                    break;   // continue pass; while-loop will rescan if needed
                }
            }
        }
        auto t4 = Clock::now();
        if (trace) {
            size_t kept = 0; for (auto x : removed_) if (!x) ++kept;
            std::fprintf(stderr, "after Phase II: kept=%zu accepted=%zu\n", kept, accepted.size());
        }
        if (timing) {
            auto ms = [](Clock::time_point a, Clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count();
            };
            std::fprintf(stderr,
                "[timing] dep-order prefs=%.1fms phaseI=%.1fms phaseII=%.1fms total=%.1fms\n"
                "         N=%d prefs=%zu phaseII_iters=%d enums=%zu tries=%zu\n"
                "         fail=%zu improve=%zu pin_path(none/fast/multi)=%zu/%zu/%zu\n"
                "         phaseII detail: enum=%.1fms try=%.1fms\n",
                ms(t0, t1), ms(t2, t3), ms(t3, t4), ms(t_all, t4),
                N, prefs.size(), guard, phase2_enum, phase2_tries,
                phase2_fail, phase2_improve,
                phase2_no_pin, phase2_fast_pin, phase2_multi_pin,
                ms_enum, ms_try);
        }
    }

    // Best add=1 candidate for target whose source lies entirely inside source_iv.
    Candidate best_add1_from_(const Interval& target, const Interval& source,
                              bool require_avail) const {
        Candidate best;
        auto cands = enumerate_candidates_(target.name, target.lo, target.hi, require_avail);
        for (auto& c : cands) {
            if (c.add != 1) continue;
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
                if (c.coverage() >= 2) score += c.coverage();
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
                    if (link.coverage() < 2) continue;
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
                    if (link.coverage() < 2) continue;
                    kids.push_back(std::move(link));
                }
                std::sort(kids.begin(), kids.end(), [](const Candidate& a, const Candidate& b){
                    return a.coverage() > b.coverage();
                });

                for (auto& link : kids) {
                    Candidate composed = compose_to_root(link, parent_tr);
                    if (composed.coverage() < 2) continue;
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

        // Leftover sweep: classic size-greedy on remaining intervals.
        auto tl0 = Clock::now();
        std::vector<const Interval*> order;
        for (const auto& iv : intervals)
            if (isize(iv) > 1 && !accepted.count(iv.name)) order.push_back(&iv);
        std::sort(order.begin(), order.end(),
                  [](const Interval* a, const Interval* b){ return (a->hi-a->lo) > (b->hi-b->lo); });
        for (const Interval* ivp : order) {
            Candidate c = best_candidate_(ivp->name, ivp->lo, ivp->hi, /*avail=*/true);
            if (c.coverage() < 2) continue;
            accept_(c, accepted);
            if (trace)
                std::fprintf(stderr, "leftover I_%s <- I_%s add=%d cov=%d\n",
                    name_to_string(G_.shape, ivp->name).c_str(),
                    name_to_string(G_.shape, name_of_rank_(c.src_lo)).c_str(),
                    c.add, c.coverage());
        }
        auto tl1 = Clock::now();
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
