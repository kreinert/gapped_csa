// sais.hpp
//
// Linear-time suffix array construction for INTEGER alphabets via SA-IS
// (induced sorting; Nong, Zhang & Chan 2009).  This replaces the earlier
// fixed-width byte encoding + libdivsufsort + boundary-subsampling trick, which
// paid a b-times memory blow-up to squeeze a large integer alphabet through a
// byte-only sorter.  SA-IS sorts the integer string directly.
//
// Convention.  We use the standard SA-IS convention that the end of the string
// is the smallest possible suffix (equivalently: a unique smallest sentinel is
// appended).  This is exactly divsufsort's implicit ordering, so `sais_int`
// returns the SAME suffix array the byte-based path produced (the paper-table
// validation is the ground truth for this).
//
// Indices are int32_t here (matching the rest of the project).  The core works
// on `int`; to move to 64-bit lengths, change the core's index type / the
// `Index` alias and `sais_int`'s return type to int64_t -- the algorithm is
// otherwise identical.

#pragma once

#include <cstdint>
#include <vector>
#include <algorithm>

namespace gcsa {
namespace sais_detail {

using Idx = int;   // internal index type for the core (switch to long for >2^31)

inline bool is_lms(const std::vector<char>& t, Idx i) {
    return i > 0 && t[i] && !t[i - 1];        // S-type preceded by L-type
}

// Fill bucket boundaries.  end=false -> bucket heads, end=true -> bucket ends.
inline void get_buckets(const std::vector<Idx>& s, std::vector<Idx>& bkt,
                        Idx n, Idx K, bool end) {
    std::fill(bkt.begin(), bkt.begin() + K, 0);
    for (Idx i = 0; i < n; ++i) ++bkt[s[i]];
    Idx sum = 0;
    for (Idx c = 0; c < K; ++c) { sum += bkt[c]; bkt[c] = end ? sum : sum - bkt[c]; }
}

// Induce L-type suffixes from already-placed (LMS or sorted-LMS) suffixes.
inline void induce_L(const std::vector<Idx>& s, std::vector<Idx>& sa,
                     std::vector<Idx>& bkt, const std::vector<char>& t, Idx n, Idx K) {
    get_buckets(s, bkt, n, K, false);         // heads
    for (Idx i = 0; i < n; ++i) {
        Idx j = sa[i] - 1;
        if (sa[i] > 0 && !t[j]) sa[bkt[s[j]]++] = j;
    }
}

// Induce S-type suffixes (scan right to left).
inline void induce_S(const std::vector<Idx>& s, std::vector<Idx>& sa,
                     std::vector<Idx>& bkt, const std::vector<char>& t, Idx n, Idx K) {
    get_buckets(s, bkt, n, K, true);          // ends
    for (Idx i = n - 1; i >= 0; --i) {
        Idx j = sa[i] - 1;
        if (sa[i] > 0 && t[j]) sa[--bkt[s[j]]] = j;
    }
}

inline bool lms_equal(const std::vector<Idx>& s, const std::vector<char>& t,
                      Idx a, Idx b) {
    for (Idx d = 0; ; ++d) {
        bool aL = (d > 0) && is_lms(t, a + d);
        bool bL = (d > 0) && is_lms(t, b + d);
        if (aL && bL) return true;                       // both LMS substrings end here
        if (aL != bL) return false;                      // different lengths
        if (s[a + d] != s[b + d] || t[a + d] != t[b + d]) return false;
    }
}

// Core SA-IS.  s: string of length n over [0,K); s[n-1] must be the unique
// smallest sentinel.  Writes the length-n suffix array into sa.
inline void sais_core(const std::vector<Idx>& s, std::vector<Idx>& sa, Idx n, Idx K) {
    if (n == 1) { sa[0] = 0; return; }

    // 1) classify L / S types (t[i]=1 means S-type).
    std::vector<char> t(n);
    t[n - 1] = 1;
    if (n >= 2) t[n - 2] = 0;
    for (Idx i = n - 3; i >= 0; --i)
        t[i] = (s[i] < s[i + 1]) || (s[i] == s[i + 1] && t[i + 1]);

    std::vector<Idx> bkt(K);

    // 2) place LMS suffixes at their bucket ends, then induce.
    get_buckets(s, bkt, n, K, true);
    std::fill(sa.begin(), sa.begin() + n, (Idx)-1);
    for (Idx i = 1; i < n; ++i)
        if (is_lms(t, i)) sa[--bkt[s[i]]] = i;
    induce_L(s, sa, bkt, t, n, K);
    induce_S(s, sa, bkt, t, n, K);

    // 3) collect sorted LMS positions into the front of sa.
    Idx n1 = 0;
    for (Idx i = 0; i < n; ++i)
        if (is_lms(t, sa[i])) sa[n1++] = sa[i];

    // 4) name the LMS substrings.
    std::vector<Idx> lms_name(n, -1);
    Idx name = 0, prev = -1;
    for (Idx i = 0; i < n1; ++i) {
        Idx pos = sa[i];
        bool diff = (prev < 0) || !lms_equal(s, t, prev, pos);
        if (diff) { ++name; prev = pos; }
        lms_name[pos] = name - 1;
    }

    // 5) build reduced string s1 (LMS positions in text order).
    std::vector<Idx> lms_order;
    lms_order.reserve(n1);
    for (Idx i = 1; i < n; ++i) if (is_lms(t, i)) lms_order.push_back(i);
    std::vector<Idx> s1(n1);
    for (Idx i = 0; i < n1; ++i) s1[i] = lms_name[lms_order[i]];

    // 6) recurse (or invert directly when all names are unique).
    std::vector<Idx> sa1(n1);
    if (name < n1) {
        sais_core(s1, sa1, n1, name);
    } else {
        for (Idx i = 0; i < n1; ++i) sa1[s1[i]] = i;
    }

    // 7) induce the final SA from the sorted LMS suffixes.
    get_buckets(s, bkt, n, K, true);
    std::fill(sa.begin(), sa.begin() + n, (Idx)-1);
    for (Idx i = n1 - 1; i >= 0; --i) {
        Idx pos = lms_order[sa1[i]];
        sa[--bkt[s[pos]]] = pos;
    }
    induce_L(s, sa, bkt, t, n, K);
    induce_S(s, sa, bkt, t, n, K);
}

} // namespace sais_detail

// Suffix array of an integer string `s` (values in [0, alphabet_size)) using the
// end-of-string-is-smallest convention (identical ordering to libdivsufsort).
// A unique smallest sentinel is handled internally, so `s` needs no sentinel.
// Returns a suffix array of length s.size().
inline std::vector<int32_t> sais_int(const std::vector<int32_t>& s, int alphabet_size) {
    using sais_detail::Idx;
    const Idx n = (Idx)s.size();
    if (n == 0) return {};

    // Shift values up by one and append the 0 sentinel (unique smallest).
    std::vector<Idx> str((size_t)n + 1);
    for (Idx i = 0; i < n; ++i) str[i] = (Idx)s[i] + 1;
    str[n] = 0;
    const Idx K = alphabet_size + 1;

    std::vector<Idx> sa((size_t)n + 1);
    sais_detail::sais_core(str, sa, n + 1, K);

    // sa[0] is the sentinel (position n); drop it. The rest is the SA of s.
    std::vector<int32_t> out((size_t)n);
    for (Idx i = 0; i < n; ++i) out[(size_t)i] = (int32_t)sa[(size_t)i + 1];
    return out;
}

} // namespace gcsa
