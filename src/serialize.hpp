// serialize.hpp
//
// Byte-level serialization of the compressed index exactly as described in the
// note: the compressed suffix array C is stored as fixed-width positions, and
// the hash table is stored in a byte array H with a select-supported bitvector
// S_H marking the start of each entry.
//
// Marker byte layout (MSB first), matching the note's worked example where the
// entry (2000000,16,300),(3000000000,24) starts with the byte 10100100:
//
//   bit7      : H_offset present (1) or nil (0)
//   bit6-bit5 : (offset count bytes) - 1        [1..4 bytes]
//   bit4-bit3 : (offset addition bytes) - 1     [1..4 bytes]
//   bit2      : H_rest present (1) or nil (0)
//   bit1-bit0 : (rest count bytes) - 1          [1..4 bytes]
//
// Data following the marker (only present parts):
//   H_offset : pos (pb bytes) | add (add-bytes) | num (count-bytes)
//   H_rest   : pos (pb bytes) | num (count-bytes)
//
// `pb` is the fixed number of bytes for a position/index (ceil(bits(|C|,n)/8)).
// `add` and `num` are limited to 4 bytes (2^32), per the note.

#pragma once

#include "compress.hpp"
#include <cstdint>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>

namespace gcsa {

// bytes_needed() now lives in compress.hpp (compress.hpp is included above)
// -- min_coverage_breakeven() needs it before any SerializedIndex exists.

inline void put_be(std::vector<uint8_t>& out, uint64_t v, int nbytes) {
    for (int k = 0; k < nbytes; ++k)
        out.push_back((uint8_t)((v >> (8 * (nbytes - 1 - k))) & 0xFF));
}
inline uint64_t get_be(const uint8_t* p, int nbytes) {
    uint64_t v = 0;
    for (int k = 0; k < nbytes; ++k) v = (v << 8) | p[k];
    return v;
}

// Elias-Fano cost, in bits, of `n` monotone values drawn from [0, u). Shared by
// every EF-priced component of the size report (S_H, the name index, and the
// uncompressed baseline's copies of both) so they are all costed identically.
//
// Split each value at bit `ell`. The low `ell` bits go into a fixed-width array
// (n*ell bits). The high parts are non-decreasing and bounded by u/2^ell, so
// they go into a bucket bitvector carrying one 1 per value and one 0 per bucket
// (n + ceil(u/2^ell) bits). That is what makes select(i) = x_i O(1):
//
//     x_i = ((select1(high, i) - i) << ell) | low[i]
//
// since select1(high, i) - i counts the zeros before the i-th one, i.e. the
// bucket index = the high part. Predecessor / membership (name -> rank) goes the
// other way via select0 plus a scan inside one bucket: O(1) expected, not worst
// case. Total space is within ~0.5 bit/value of log2(C(u,n)), the information-
// theoretic floor for storing the set at all.
//
// `ell` is chosen to minimize the total: one more low bit costs n bits and
// halves the bucket count, so the optimum sits at ell ~ log2(u/n) and checking
// floor(log2(u/n)) with one step either side is exact.
//
// NOTE (changed Sep 18 2026): this used to return n*(2 + floor(log2(u/n))),
// i.e. it assumed the bucket array was exactly 2n bits. With the floor choice of
// ell the bucket array is n + ceil(u/2^ell), which runs up to 3n, so the old
// formula was only tight when u/n sat near a power of two. It understated the
// human frequent-k-mer run's S_H (u/n = 6.2 = 2^2.64, the worst case) by 14%
// and `names` (u/n = 2^15.27) by 1.2%. Both sides of the uncompressed-baseline
// comparison carried the same understatement, so ratios from older logs are
// still sound (88.45% -> 88.51% on that run); absolute S_H/names byte counts
// from before this date are not comparable with ones printed after it.
//
// Not counted, on either side: the o(n) rank/select support over the bucket
// bitvector, ~0.1-0.25 bits/value for a practical sd_vector-style structure.
inline size_t ef_bits(size_t n, uint64_t u) {
    if (n == 0) return 0;
    if (u < n) u = n;                       // degenerate; keep the split well-defined
    const uint64_t q = u / n;               // >= 1, and floor(log2(q)) == floor(log2(u/n))
    int base = 0;
    while (base < 62 && (uint64_t(1) << (base + 1)) <= q) ++base;

    size_t best = SIZE_MAX;
    for (int ell = std::max(0, base - 1); ell <= base + 1; ++ell) {
        const uint64_t buckets = (u + ((uint64_t(1) << ell) - 1)) >> ell;   // ceil(u / 2^ell)
        const size_t bits = n * (size_t)ell + n + (size_t)buckets;
        if (bits < best) best = bits;
    }
    return best;
}

// Serialized field widths of one hash-table entry. Factored out so that the
// three places that care about the encoding -- emit_entry (writes the bytes),
// entry_bytes (counts them without writing), and positions_of (reads them back
// off the marker) -- cannot drift apart.
struct EntryWidths {
    int off_add_b = 1, off_cnt_b = 1, rest_cnt_b = 1;
};

inline EntryWidths entry_widths(const HashEntry& e) {
    EntryWidths w;
    if (e.has_offset) {
        w.off_add_b = std::max(1, bytes_needed(e.off_add));
        w.off_cnt_b = std::max(1, bytes_needed(e.off_num));
        if (w.off_add_b > 4 || w.off_cnt_b > 4) throw std::runtime_error("add/num exceed 4 bytes");
    }
    if (e.has_rest) {
        w.rest_cnt_b = std::max(1, bytes_needed(e.rest_num));
        if (w.rest_cnt_b > 4) throw std::runtime_error("rest num exceeds 4 bytes");
    }
    return w;
}

// Bytes one entry occupies in H: the marker byte plus whichever of
// H_offset = (pos, add, num) / H_rest = (pos, num) are present. Must stay equal
// to what emit_entry() appends -- both derive their widths from entry_widths(),
// and the self-test path checks the two agree byte for byte.
inline size_t entry_bytes(const HashEntry& e, int pb) {
    const EntryWidths w = entry_widths(e);
    size_t b = 1;  // marker
    if (e.has_offset) b += (size_t)pb + w.off_add_b + w.off_cnt_b;
    if (e.has_rest)   b += (size_t)pb + w.rest_cnt_b;
    return b;
}

// Composition of a hash table and where its H bytes actually go. Two questions
// the single `H=` figure cannot answer: how many k-mers ended up with an
// H_offset at all (i.e. how much of the table compression actually touched),
// and how H divides between the marker bytes, the position fields and the
// count/add fields. The position fields are the interesting ones -- SA
// intervals are contiguous and in name order, so `rest_pos` is in principle
// derivable rather than stored, and rest_pos_bytes is exactly what dropping it
// would reclaim.
//
// Invariants, checked in the self-test path:
//   total_bytes()    == the serialized index's bytes_H()
//   offset_positions == m - |C|     (every position not kept is reconstructed)
//   rest_positions   == |C|         (every kept position is named exactly once)
struct HashStats {
    size_t entries = 0;
    size_t rest_only = 0, offset_only = 0, both = 0, empty = 0;
    // H byte breakdown; sums to total_bytes()
    size_t marker_bytes = 0;
    size_t off_pos_bytes = 0, off_add_bytes = 0, off_num_bytes = 0;
    size_t rest_pos_bytes = 0, rest_num_bytes = 0;
    // link statistics
    size_t offset_positions = 0;   // sum of off_num
    size_t rest_positions = 0;     // sum of rest_num

    size_t offsets() const { return offset_only + both; }
    size_t total_bytes() const {
        return marker_bytes + off_pos_bytes + off_add_bytes + off_num_bytes
             + rest_pos_bytes + rest_num_bytes;
    }
    double avg_coverage() const {
        return offsets() ? (double)offset_positions / (double)offsets() : 0.0;
    }
};

// Fold one entry into a HashStats. Field widths come from entry_widths(), the
// same source emit_entry() writes from, so the breakdown cannot disagree with
// the bytes actually emitted.
inline void accumulate_entry(HashStats& S, const HashEntry& e, int pb) {
    const EntryWidths w = entry_widths(e);
    ++S.entries;
    S.marker_bytes += 1;
    if (e.has_offset) {
        S.off_pos_bytes += (size_t)pb;
        S.off_add_bytes += (size_t)w.off_add_b;
        S.off_num_bytes += (size_t)w.off_cnt_b;
        S.offset_positions += e.off_num;
    }
    if (e.has_rest) {
        S.rest_pos_bytes += (size_t)pb;
        S.rest_num_bytes += (size_t)w.rest_cnt_b;
        S.rest_positions += e.rest_num;
    }
    if (e.has_offset && e.has_rest) ++S.both;
    else if (e.has_offset)          ++S.offset_only;
    else if (e.has_rest)            ++S.rest_only;
    else                            ++S.empty;
}

inline HashStats hash_stats(const CompressedIndex& idx, int pb) {
    HashStats S;
    for (const auto& kv : idx.table()) accumulate_entry(S, kv.second, pb);
    return S;
}

// Append one entry to H: marker byte (layout at the top of this file) followed
// by the present parts. The single writer of the encoding -- serialize_index()
// and the uncompressed baseline both go through it.
inline void emit_entry(std::vector<uint8_t>& H, const HashEntry& e, int pb) {
    const EntryWidths w = entry_widths(e);
    uint8_t marker = 0;
    if (e.has_offset) marker |= 0x80;
    marker |= (uint8_t)((w.off_cnt_b - 1) & 0x3) << 5;
    marker |= (uint8_t)((w.off_add_b - 1) & 0x3) << 3;
    if (e.has_rest) marker |= 0x04;
    marker |= (uint8_t)((w.rest_cnt_b - 1) & 0x3);
    H.push_back(marker);
    if (e.has_offset) {
        put_be(H, e.off_pos, pb);
        put_be(H, e.off_add, w.off_add_b);
        put_be(H, e.off_num, w.off_cnt_b);
    }
    if (e.has_rest) {
        put_be(H, e.rest_pos, pb);
        put_be(H, e.rest_num, w.rest_cnt_b);
    }
}

// A minimal select-supporting bitvector (dense sampling). Marks entry starts in
// H; select(i) returns the bit position (= byte offset) of the i-th entry.
// For a production build this is where an Elias-Fano vector would go.
struct SelectBitvector {
    std::vector<uint64_t> starts;              // sorted set bit positions
    size_t nbits = 0;
    uint64_t select(size_t i) const { return starts[i]; }
    size_t count() const { return starts.size(); }
    // Heuristic serialized size of an EF representation of these positions:
    // n*(2 + ceil(log2(u/n))) bits.  (We store `starts` in memory for speed.)
    size_t ef_bits() const {
        size_t n = starts.size(); if (n == 0) return 0;
        uint64_t u = nbits ? nbits : (starts.back() + 1);
        return gcsa::ef_bits(n, u);
    }
};

struct SerializedIndex {
    int pb = 4;                                 // bytes per position/index
    int span = 1;
    std::vector<uint8_t>  C_bytes;              // |C| positions, pb bytes each
    size_t                C_count = 0;
    std::vector<uint8_t>  H;                    // hash-table bytes
    SelectBitvector       S_H;                  // entry starts in H
    std::vector<uint64_t> names_sorted;         // occurring k-mer names, ascending

    // ---- sizes ---------------------------------------------------------------
    size_t bytes_C()  const { return C_bytes.size(); }
    size_t bytes_H()  const { return H.size(); }
    size_t bytes_SH() const { return (S_H.ef_bits() + 7) / 8; }
    // name index: an EF/rank structure over names; approximate with an EF cost.
    size_t bytes_names() const {
        size_t n = names_sorted.size(); if (!n) return 0;
        return (gcsa::ef_bits(n, names_sorted.back() + 1) + 7) / 8;
    }
    size_t total_bytes() const { return bytes_C() + bytes_H() + bytes_SH() + bytes_names(); }

    // ---- lookup --------------------------------------------------------------
    int64_t position_at(size_t idx) const { return (int64_t)get_be(&C_bytes[idx * pb], pb); }

    // rank of name among occurring names (or -1 if absent).
    long name_rank(uint64_t name) const {
        auto it = std::lower_bound(names_sorted.begin(), names_sorted.end(), name);
        if (it == names_sorted.end() || *it != name) return -1;
        return (long)(it - names_sorted.begin());
    }

    // Decode positions for a k-mer directly from the serialized bytes.
    std::vector<int64_t> positions_of(uint64_t name) const {
        std::vector<int64_t> out;
        long i = name_rank(name);
        if (i < 0) return out;
        const uint8_t* p = &H[S_H.select((size_t)i)];
        uint8_t marker = *p++;
        bool has_off = marker & 0x80;
        int  off_cnt_b = ((marker >> 5) & 0x3) + 1;
        int  off_add_b = ((marker >> 3) & 0x3) + 1;
        bool has_rest = marker & 0x04;
        int  rest_cnt_b = (marker & 0x3) + 1;
        if (has_off) {
            uint64_t pos = get_be(p, pb);          p += pb;
            uint64_t add = get_be(p, off_add_b);   p += off_add_b;
            uint64_t num = get_be(p, off_cnt_b);   p += off_cnt_b;
            for (uint64_t t = 0; t < num; ++t)
                out.push_back(position_at(pos + t) + (int64_t)add * span);
        }
        if (has_rest) {
            uint64_t pos = get_be(p, pb);          p += pb;
            uint64_t num = get_be(p, rest_cnt_b);  p += rest_cnt_b;
            for (uint64_t t = 0; t < num; ++t)
                out.push_back(position_at(pos + t));
        }
        return out;
    }
};

inline SerializedIndex serialize_index(const CompressedIndex& idx) {
    SerializedIndex S;
    // Not shape.span: that's the k-mer/window width, which diverges from the
    // add-shift stride in regular (ungapped) mode (window=k, stride=1). Pull
    // the stride straight from the index so both modes decode correctly.
    S.span = idx.add_stride();

    const std::vector<int64_t>& C = idx.compressed_positions();
    S.C_count = C.size();
    uint64_t maxpos = C.empty() ? 0 : (uint64_t)*std::max_element(C.begin(), C.end());
    uint64_t maxidx = C.empty() ? 0 : (C.size() - 1);
    S.pb = bytes_needed(std::max(maxpos, maxidx));
    for (int64_t v : C) put_be(S.C_bytes, (uint64_t)v, S.pb);

    // Entries ordered by k-mer name.
    std::map<uint64_t, HashEntry> ordered(idx.table().begin(), idx.table().end());
    for (auto& kv : ordered) {
        const HashEntry& e = kv.second;
        S.names_sorted.push_back(kv.first);
        S.S_H.starts.push_back(S.H.size());          // entry start (byte offset)

        emit_entry(S.H, e, S.pb);
    }
    S.S_H.nbits = S.H.size();
    return S;
}

// ---------------------------------------------------------------------------
// The uncompressed baseline
//
// The size line's "full SA" figure compares this index against a bare suffix
// array -- but a bare SA is not a structure that answers our query. It has no
// k-mer -> interval map, so locating a k-mer in it costs a binary search over
// the text, not constant time. The like-for-like competitor is the
// *uncompressed* index with the same query interface:
//
//     the full suffix array (all m positions, nothing shared)
//   + a hash table in which every k-mer's entry is a single
//     H_rest = (pos, num) pointing at that k-mer's own SA interval I_c.
//
// No H_offset ever appears there, because nothing is reconstructed from
// anything else. Everything else is deliberately held fixed -- same marker-byte
// layout, same variable-width count fields, same EF cost model for S_H and the
// name index, same pb rule -- so the gap between the two totals is exactly what
// differential compression buys and costs: C shrinks from m to |C| positions,
// while H grows by one H_offset (pb + add + num bytes) per compressed k-mer.
//
// uncompressed_baseline_sizes() is size-only and never materializes the m*pb
// baseline position array (2.2 GB on the human frequent-k-mer input): it is one
// pass over the hash table, because a k-mer's interval size |I_c| is recoverable
// from its entry as off_num + rest_num -- exactly the count positions_of()
// returns. serialize_uncompressed_baseline() below builds the thing for real,
// for round-trip validation on small inputs.
// ---------------------------------------------------------------------------
struct BaselineSizes {
    int    pb = 4;            // bytes per position/index in the baseline
    size_t entries = 0;       // hash-table entries (= distinct k-mers)
    size_t positions = 0;     // m -- the full SA keeps every entry
    size_t bytes_C = 0;       // the full SA itself: m * pb
    size_t bytes_H = 0;
    size_t bytes_SH = 0;
    size_t bytes_names = 0;
    HashStats stats;          // same breakdown as the index's, for comparison
    size_t total_bytes() const { return bytes_C + bytes_H + bytes_SH + bytes_names; }
};

inline BaselineSizes uncompressed_baseline_sizes(const CompressedIndex& idx) {
    BaselineSizes B;
    const size_t m = idx.total_positions();
    const size_t n = idx.gsa().n;
    B.positions = m;
    // serialize_index()'s pb rule with the baseline's own maxima: its position
    // array is indexed over all m entries rather than |C|, and it stores the
    // same text positions. In regular mode lex2orig is the identity, so m-1 is
    // itself a stored position; in gapped mode positions are < n <= m-1.
    B.pb = bytes_needed(std::max<uint64_t>(m ? m - 1 : 0, n ? n - 1 : 0));
    B.bytes_C = m * (size_t)B.pb;

    uint64_t max_name = 0;
    for (const auto& kv : idx.table()) {
        const HashEntry& e = kv.second;
        HashEntry base;                            // rest-only: I_c stored directly
        base.has_rest = true;
        base.rest_pos = 0;                         // fixed width pb; the value is irrelevant here
        base.rest_num = e.off_num + e.rest_num;    // = |I_c|
        accumulate_entry(B.stats, base, B.pb);
        max_name = std::max(max_name, kv.first);
        ++B.entries;
    }
    B.bytes_H = B.stats.total_bytes();
    B.bytes_SH    = (ef_bits(B.entries, B.bytes_H) + 7) / 8;
    B.bytes_names = (ef_bits(B.entries, max_name + 1) + 7) / 8;
    return B;
}

// The same baseline, actually built and encoded through the identical writer.
// Costs m*pb bytes of memory plus a pass over the SA, so it is meant for the
// self-test path and small inputs -- use uncompressed_baseline_sizes() for the
// size report on real inputs. Returning a SerializedIndex lets the caller
// round-trip it: the baseline has to answer the same queries this index does,
// or it is not a fair competitor.
inline SerializedIndex serialize_uncompressed_baseline(const CompressedIndex& idx) {
    const GappedSA& G = idx.gsa();
    const size_t m = G.m();
    SerializedIndex S;
    S.span = idx.add_stride();
    S.C_count = m;

    uint64_t maxpos = 0;
    for (size_t r = 0; r < m; ++r)
        maxpos = std::max(maxpos, (uint64_t)G.orig_pos((int32_t)r));
    S.pb = bytes_needed(std::max<uint64_t>(maxpos, m ? m - 1 : 0));
    S.C_bytes.reserve(m * (size_t)S.pb);
    for (size_t r = 0; r < m; ++r)
        put_be(S.C_bytes, (uint64_t)G.orig_pos((int32_t)r), S.pb);

    // One rest-only entry per SA interval. sa_intervals() already yields them
    // in name order, which is the order serialize_index() stores entries in.
    for (const auto& iv : idx.sa_intervals()) {
        HashEntry e;
        e.has_rest = true;
        e.rest_pos = (uint64_t)iv.lo;
        e.rest_num = (uint32_t)(iv.hi - iv.lo);
        S.names_sorted.push_back(iv.name);
        S.S_H.starts.push_back(S.H.size());
        emit_entry(S.H, e, S.pb);
    }
    S.S_H.nbits = S.H.size();
    return S;
}

} // namespace gcsa
