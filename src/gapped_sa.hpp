// gapped_sa.hpp
//
// Construction of the gapped suffix array for an arbitrary shape via the DisLex
// transform (Horton 2008).  The suffix sorting is done by a linear-time SA-IS
// for integer alphabets (see sais.hpp).
//
// DisLex in one paragraph.  For a shape of span s we form, for every text start
// position p, the gapped k-mer name at p (see shape.hpp).  We then lay these
// names out grouped by p mod s (mod 0 group first, then mod 1, ...).  The suffix
// array of this integer string ("lextext") is exactly the gapped suffix array:
// consecutive symbols within a mod group are the names at p, p+s, p+2s, ... i.e.
// exactly the gapped suffix read at p.  Following the note we do NOT insert
// separators between groups; boundary windows read into '$' padding, and the
// unique all-'$' name (= 0) at the very end acts as the global terminator.
//
// The lextext names are large integers (up to ~2^62), so we sort them with an
// integer-alphabet SA-IS.  Before sorting we order-preservingly remap the
// distinct names to compact ranks 0..k-1 (sort-unique + lower_bound) so the
// alphabet fed to SA-IS is dense and small; the original arithmetic names are
// kept in G.lex for hashing/lookup.  The LCP array is computed in *symbol* space
// with Kasai's algorithm over G.lex.

#pragma once

#include "shape.hpp"
#include "sais.hpp"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>

namespace gcsa {

struct GappedSA {
    Shape shape;
    std::string text;                 // cleaned DNA text (no padding)
    size_t n = 0;                     // text length
    int alphabet_size = 0;            // number of distinct lextext names (SA-IS alphabet)

    std::vector<uint64_t> lex;        // lextext names, length m = n+1 (DisLex order)
    std::vector<int64_t>  lex2orig;   // lextext position -> original text start position
    std::vector<int32_t>  sa;         // gapped SA: sa[r] = lextext position, length m
    std::vector<int32_t>  lcp;        // lcp[r] = lcp(sa[r-1], sa[r]) in symbols; lcp[0]=0

    size_t m() const { return lex.size(); }

    // Original text position of the r-th suffix in the gapped SA.
    int64_t orig_pos(int32_t rank) const { return lex2orig[sa[rank]]; }

    // First symbol (= gapped k-mer name) of the r-th suffix.
    uint64_t first_symbol(int32_t rank) const { return lex[sa[rank]]; }
};

// Kasai's algorithm: LCP array (in units of `lex` symbols) for a suffix array
// `sa` of the integer string `lex`. lcp[r] = shared-prefix length (in `lex`
// symbols) between the suffixes at sa[r-1] and sa[r]; lcp[0] = 0. Shared by
// both build_gapped_sa (symbol space = gapped k-mer names) and build_plain_sa
// (symbol space = raw characters) -- the algorithm itself doesn't care what
// the symbols mean.
inline std::vector<int32_t> kasai_lcp(const std::vector<uint64_t>& lex,
                                       const std::vector<int32_t>& sa) {
    const size_t m = lex.size();
    std::vector<int32_t> rank(m);
    for (size_t r = 0; r < m; ++r) rank[sa[r]] = (int32_t)r;
    std::vector<int32_t> lcp(m, 0);
    int h = 0;
    for (size_t i = 0; i < m; ++i) {
        if (rank[i] > 0) {
            size_t j = sa[rank[i] - 1];
            while (i + h < m && j + h < m && lex[i + h] == lex[j + h]) ++h;
            lcp[rank[i]] = h;
            if (h > 0) --h;
        } else {
            h = 0;
        }
    }
    return lcp;
}

// Plain (ungapped) suffix array of the raw text, built directly over the
// small {$,A,C,G,T} alphabet -- no per-position k-mer name computation, no
// mod-residue grouping, no alphabet remap. `k` only affects G.shape (used for
// naming k-mers of that width downstream, e.g. interval identity in
// compress.hpp); the SA/LCP themselves don't depend on k at all. This isolates
// the cost of DisLex's per-position name computation + big-alphabet sort from
// suffix-sorting itself: for a contiguous shape "######" (span==weight==k),
// build_gapped_sa produces the same total suffix order as this function.
inline GappedSA build_plain_sa(int k, std::string text) {
    using Clock = std::chrono::steady_clock;
    const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
    auto ms = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };

    GappedSA G;
    G.shape = Shape::parse(std::string((size_t)k, '#'));
    G.text  = std::move(text);
    G.n     = G.text.size();

    // 1) Map characters directly to a small dense alphabet: $ = 0 (unique
    //    sentinel at position n), A/C/G/T = 1..4. No name_at() calls.
    auto t0 = Clock::now();
    const size_t m = G.n + 1;
    G.lex.resize(m);
    G.lex2orig.resize(m);
    for (size_t i = 0; i < m; ++i) {
        G.lex[i] = (i < G.n) ? (uint64_t)(base_value(G.text[i]) + 1) : 0;
        G.lex2orig[i] = (int64_t)i;   // identity: no mod-grouping in plain mode
    }
    G.alphabet_size = SIGMA + 1;
    auto t1 = Clock::now();

    // 2) Suffix-sort the small-alphabet character string directly with SA-IS.
    std::vector<int32_t> chars(m);
    for (size_t i = 0; i < m; ++i) chars[i] = (int32_t)G.lex[i];
    G.sa = sais_int(chars, G.alphabet_size);
    auto t2 = Clock::now();

    // 3) Kasai LCP in character space.
    G.lcp = kasai_lcp(G.lex, G.sa);
    auto t3 = Clock::now();

    if (timing) {
        std::fprintf(stderr,
            "[timing] plain charmap=%.1fms  SA-IS=%.1fms  LCP=%.1fms  m=%zu sigma=%d\n",
            ms(t0, t1), ms(t1, t2), ms(t2, t3), m, G.alphabet_size);
    }
    return G;
}

inline GappedSA build_gapped_sa(const Shape& shape, std::string text) {
    using Clock = std::chrono::steady_clock;
    const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
    auto ms = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };

    GappedSA G;
    G.shape = shape;
    G.text  = std::move(text);
    G.n     = G.text.size();

    // 1) Build the lextext in DisLex (grouped-by-residue) order.
    auto t0 = Clock::now();
    const int s = shape.span;
    G.lex.reserve(G.n + 1);
    G.lex2orig.reserve(G.n + 1);
    for (int mod = 0; mod < s; ++mod) {
        for (long p = mod; p <= (long)G.n; p += s) {
            G.lex.push_back(name_at(shape, G.text, (size_t)p));
            G.lex2orig.push_back(p);
        }
    }
    const size_t m = G.lex.size();
    auto t1 = Clock::now();

    // 2) Order-preservingly remap the distinct names to compact ranks 0..k-1.
    //    (This keeps the alphabet passed to SA-IS dense; G.lex keeps the
    //    original arithmetic names for hashing/lookup.)
    std::vector<uint64_t> uniq(G.lex);
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    G.alphabet_size = (int)uniq.size();

    std::vector<int32_t> remapped(m);
    for (size_t i = 0; i < m; ++i)
        remapped[i] = (int32_t)(std::lower_bound(uniq.begin(), uniq.end(), G.lex[i]) - uniq.begin());
    auto t2 = Clock::now();

    // 3) Suffix-sort the integer string directly with SA-IS.
    G.sa = sais_int(remapped, G.alphabet_size);
    auto t3 = Clock::now();

    // 4) Kasai LCP in symbol space.
    G.lcp = kasai_lcp(G.lex, G.sa);
    auto t4 = Clock::now();
    if (timing) {
        std::fprintf(stderr,
            "[timing] DisLex+name=%.1fms  remap=%.1fms  SA-IS=%.1fms  LCP=%.1fms  m=%zu sigma=%d\n",
            ms(t0, t1), ms(t1, t2), ms(t2, t3), ms(t3, t4), m, G.alphabet_size);
    }
    return G;
}

} // namespace gcsa
