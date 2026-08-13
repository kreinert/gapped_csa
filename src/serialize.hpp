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
        int lo = 0; while ((uint64_t(1) << (lo + 1)) * n <= u) ++lo;
        return n * (size_t)(2 + lo);
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
        uint64_t u = names_sorted.back() + 1;
        int lo = 0; while ((uint64_t(1) << (lo + 1)) * n <= u) ++lo;
        return (n * (size_t)(2 + lo) + 7) / 8;
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
    S.span = idx.gsa().shape.span;

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

        int off_add_b = 1, off_cnt_b = 1, rest_cnt_b = 1;
        if (e.has_offset) {
            off_add_b = std::max(1, bytes_needed(e.off_add));
            off_cnt_b = std::max(1, bytes_needed(e.off_num));
            if (off_add_b > 4 || off_cnt_b > 4) throw std::runtime_error("add/num exceed 4 bytes");
        }
        if (e.has_rest) {
            rest_cnt_b = std::max(1, bytes_needed(e.rest_num));
            if (rest_cnt_b > 4) throw std::runtime_error("rest num exceeds 4 bytes");
        }
        uint8_t marker = 0;
        if (e.has_offset) marker |= 0x80;
        marker |= (uint8_t)((off_cnt_b - 1) & 0x3) << 5;
        marker |= (uint8_t)((off_add_b - 1) & 0x3) << 3;
        if (e.has_rest) marker |= 0x04;
        marker |= (uint8_t)((rest_cnt_b - 1) & 0x3);
        S.H.push_back(marker);

        if (e.has_offset) {
            put_be(S.H, e.off_pos, S.pb);
            put_be(S.H, e.off_add, off_add_b);
            put_be(S.H, e.off_num, off_cnt_b);
        }
        if (e.has_rest) {
            put_be(S.H, e.rest_pos, S.pb);
            put_be(S.H, e.rest_num, rest_cnt_b);
        }
    }
    S.S_H.nbits = S.H.size();
    return S;
}

} // namespace gcsa
