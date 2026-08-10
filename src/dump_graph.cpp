// dump_graph.cpp — print preferred candidates and the dep-order DAG
// for a text/shape (default: the note's #.# example).

#include "compress.hpp"
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <vector>

using namespace gcsa;

// Replicate the candidate enumeration used inside CompressedIndex by building
// a thin wrapper that exposes the same logic via a public dump helper.
// We duplicate the minimal computation here for clarity (same formulas).

struct Cand {
    uint64_t name;
    std::string kmer;
    int add;
    int32_t src_lo, src_hi;
    uint64_t src_name;
    std::string src_kmer;
    int coverage;
    std::vector<int32_t> covered;
};

static int64_t pred_lexpos(const GappedSA& G, int32_t r, int add) {
    int64_t x = G.sa[r];
    int64_t y = x - add;
    if (y < 0) return -1;
    if (G.lex2orig[y] != G.lex2orig[x] - (int64_t)add * G.shape.span) return -1;
    return y;
}

static std::vector<Cand> best_per_interval(const GappedSA& G, int max_add) {
    const size_t m = G.m();
    std::vector<int32_t> rank_of(m);
    for (size_t r = 0; r < m; ++r) rank_of[G.sa[r]] = (int32_t)r;

    std::vector<Cand> prefs;
    for (size_t r = 0; r < m; ) {
        uint64_t name = G.first_symbol((int32_t)r);
        size_t r2 = r + 1;
        while (r2 < m && G.first_symbol((int32_t)r2) == name) ++r2;
        int32_t lo = (int32_t)r, hi = (int32_t)r2;
        // Preference relation: only |I_c| > 2 (match DepOrder).
        if (hi - lo <= 2) { r = r2; continue; }

        Cand best{};
        best.coverage = 0;
        for (int add = 1; add <= max_add; ++add) {
            std::vector<std::pair<int32_t,int32_t>> pr;
            for (int32_t i = lo; i < hi; ++i) {
                int64_t y = pred_lexpos(G, i, add);
                if (y >= 0) pr.push_back({rank_of[y], i});
            }
            std::sort(pr.begin(), pr.end());
            size_t i = 0;
            while (i < pr.size()) {
                size_t j = i + 1;
                while (j < pr.size()
                       && pr[j].first == pr[j-1].first + 1
                       && G.lcp[pr[j].first] >= add + 1) ++j;
                int32_t s_lo = pr[i].first, s_hi = pr[j-1].first + 1;
                bool disjoint = (s_hi <= lo || s_lo >= hi);
                int cov = (int)(j - i);
                if (disjoint && cov >= 2
                    && (cov > best.coverage || (cov == best.coverage && add > best.add))) {
                    best.name = name;
                    best.kmer = name_to_string(G.shape, name);
                    best.add = add;
                    best.src_lo = s_lo;
                    best.src_hi = s_hi;
                    best.src_name = G.first_symbol(s_lo);
                    best.src_kmer = name_to_string(G.shape, best.src_name);
                    best.coverage = cov;
                    best.covered.clear();
                    for (size_t t = i; t < j; ++t) best.covered.push_back(pr[t].second);
                }
                i = j;
            }
        }
        if (best.coverage >= 2) prefs.push_back(best);
        r = r2;
    }
    return prefs;
}

int main(int argc, char** argv) {
    std::string text = "ACGTCTTAAACCCTCGTCTTAAACCCAACGTCTTAAACCC";
    std::string shape = "#.#";
    int max_add = 8;
    if (argc >= 2) text = argv[1];
    if (argc >= 3) shape = argv[2];
    if (argc >= 4) max_add = std::stoi(argv[3]);

    Shape sh = Shape::parse(shape);
    GappedSA G = build_gapped_sa(sh, text);
    auto prefs = best_per_interval(G, max_add);

    std::cout << "text  = " << text << "\n";
    std::cout << "shape = " << shape << "  m=" << G.m() << "\n\n";

    std::cout << "=== Preferred candidates (ignore availability) ===\n";
    std::cout << std::left
              << std::setw(6) << "kmer"
              << std::setw(8) << "|I_c|"
              << std::setw(6) << "add"
              << std::setw(8) << "cov"
              << std::setw(10) << "src_kmer"
              << std::setw(14) << "src_ranks"
              << "covered_ranks\n";

    // interval sizes
    std::map<uint64_t, int> isize;
    for (size_t r = 0; r < G.m(); ) {
        uint64_t name = G.first_symbol((int32_t)r);
        size_t r2 = r + 1;
        while (r2 < G.m() && G.first_symbol((int32_t)r2) == name) ++r2;
        isize[name] = (int)(r2 - r);
        r = r2;
    }

    for (auto& c : prefs) {
        std::cout << std::setw(6) << c.kmer
                  << std::setw(8) << isize[c.name]
                  << std::setw(6) << c.add
                  << std::setw(8) << c.coverage
                  << std::setw(10) << c.src_kmer
                  << "[" << c.src_lo << "," << c.src_hi << ")   ";
        std::cout << "[";
        for (size_t i = 0; i < c.covered.size(); ++i)
            std::cout << c.covered[i] << (i + 1 < c.covered.size() ? "," : "");
        std::cout << "]\n";
    }

    std::cout << "\n=== Dependency DAG  (edge: source_kmer → target_kmer, cov>=3) ===\n";
    std::cout << "Meaning: target prefers a differential source inside source_kmer.\n"
              << "Edges only for preferred candidates with coverage >= 3.\n\n";

    // Collect all nodes that appear
    std::map<std::string, std::vector<std::string>> outs;
    std::map<std::string, int> indeg;
    for (auto& c : prefs) {
        indeg[c.kmer]; // ensure node
        indeg[c.src_kmer];
        if (c.coverage < 3) continue;  // no preference-DAG edge for cov < 3
        if (c.kmer == c.src_kmer) continue;
        outs[c.src_kmer].push_back(c.kmer + "(add=" + std::to_string(c.add)
                                   + ",cov=" + std::to_string(c.coverage) + ")");
        indeg[c.kmer]++;
        // also ensure src has indeg entry counted only for targets
    }
    // recompute indeg properly
    indeg.clear();
    for (auto& c : prefs) {
        indeg[c.kmer] = 0;
        indeg[c.src_kmer] = 0;
    }
    for (auto& c : prefs) {
        if (c.coverage < 3) continue;
        if (c.kmer == c.src_kmer) continue;
        indeg[c.kmer]++;
    }

    for (auto& kv : outs) {
        std::cout << "  " << kv.first << "  →  ";
        for (size_t i = 0; i < kv.second.size(); ++i)
            std::cout << kv.second[i] << (i + 1 < kv.second.size() ? ", " : "");
        std::cout << "\n";
    }

    // Nodes with no outgoing preferred edge
    std::cout << "\nNodes with no preferred outgoing edge (never used as preferred source):\n  ";
    {
        std::map<std::string, bool> has_out;
        for (auto& kv : outs) has_out[kv.first] = true;
        bool any = false;
        for (auto& c : prefs) {
            if (!has_out[c.kmer] && !has_out[c.src_kmer]) { /* skip */ }
        }
        // list all k-mers in SA
        for (auto& kv : isize) {
            std::string k = name_to_string(sh, kv.first);
            if (!has_out[k]) { std::cout << k << " "; any = true; }
        }
        if (!any) std::cout << "(none)";
        std::cout << "\n";
    }

    std::cout << "\nIndegree (how many names prefer this node as source):\n";
    // Actually indeg on target = number of preferred edges INTO target (always 0 or 1)
    // Better: count how often a kmer is chosen as src
    std::map<std::string, int> as_src;
    for (auto& c : prefs) as_src[c.src_kmer]++;
    for (auto& kv : as_src)
        std::cout << "  " << kv.first << " is preferred source for " << kv.second << " name(s)\n";

    // Topo / reverse topo
    std::cout << "\n=== Mermaid diagram ===\n";
    std::cout << "```mermaid\nflowchart LR\n";
    for (auto& c : prefs) {
        if (c.coverage < 3) continue;
        if (c.kmer == c.src_kmer) continue;
        std::cout << "  " << c.src_kmer << " -->|add=" << c.add
                  << " cov=" << c.coverage << "| " << c.kmer << "\n";
    }
    std::cout << "```\n";

    return 0;
}
