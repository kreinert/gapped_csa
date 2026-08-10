// shape.hpp
//
// Gapped shape parsing and the order-preserving arithmetic naming of gapped
// k-mers described in the "Compressed SA for gapped shapes" note.
//
// A *shape* is a string over {'#','.'}, e.g. "#.#".  The '#' positions are the
// care positions that are actually read; the number of '#' is the weight w and
// the length of the string is the span s.
//
// Naming.  Each gapped k-mer is a sequence of w characters drawn from
// {$, A, C, G, T} with the constraint that '$' (the text terminator) only ever
// occurs as a *suffix* of the k-mer (it comes from padding past the end of the
// text).  We assign to every such valid k-mer its lexicographic rank, treating
// '$' < A < C < G < T.  This is an order-preserving name and it reproduces the
// arithmetic codes of the note (AA=2, C$=6, TT=20, ...), so we can use it
// directly as the integer symbol for the DisLex lextext without a separate rank
// determination step.
//
// Derivation of the closed form.  Let g(i) be the number of valid k-mer
// suffixes of length (w-i): g(w)=1 and g(i)=1+sigma*g(i+1) (either the position
// is '$' and everything after is '$' -> the "+1", or it is one of sigma bases
// followed by any valid suffix).  Then the name of a k-mer c_0..c_{w-1} is the
// sum over positions, stopping at the first '$', of  1 + v_i * g(i+1)  where
// v_i in {0..sigma-1} is the 0-based base value (A=0,C=1,G=2,T=3).

#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <stdexcept>

namespace gcsa {

constexpr int SIGMA = 4;                 // DNA alphabet size (A,C,G,T)
constexpr int8_t DOLLAR = -1;            // sentinel value used inside a window

// Map a DNA character to 0..3, or DOLLAR for the '$' sentinel.
// Any non-ACGT character is treated as 'A' (matches the existing tools here).
inline int8_t base_value(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        case '$':           return DOLLAR;
        default:            return 0;
    }
}

inline char base_char(int v) {
    switch (v) { case 0: return 'A'; case 1: return 'C'; case 2: return 'G'; case 3: return 'T'; }
    return '$';
}

struct Shape {
    std::string pattern;          // e.g. "#.#"
    std::vector<int> care;        // 0-based offsets of '#'
    int span = 0;                 // pattern length
    int weight = 0;               // number of '#'
    std::vector<uint64_t> g;      // g[0..weight], g[weight]=1

    static Shape parse(const std::string& p) {
        Shape s;
        s.pattern = p;
        for (int i = 0; i < (int)p.size(); ++i) {
            if (p[i] == '#') s.care.push_back(i);
            else if (p[i] != '.') throw std::runtime_error("invalid shape char (use # and .)");
        }
        if (s.care.empty()) throw std::runtime_error("shape must contain at least one '#'");
        s.span = (int)p.size();
        s.weight = (int)s.care.size();
        // Precompute g(i).
        s.g.assign(s.weight + 1, 0);
        s.g[s.weight] = 1;
        for (int i = s.weight - 1; i >= 0; --i)
            s.g[i] = 1 + (uint64_t)SIGMA * s.g[i + 1];
        return s;
    }

    // Largest possible name + 1 (i.e. the size of the name space) = g[0].
    uint64_t name_space() const { return g[0]; }

    // Compute the order-preserving name of the window whose care characters are
    // given as base values (0..3) or DOLLAR.  vals must have length == weight.
    uint64_t name_from_values(const int8_t* vals) const {
        uint64_t name = 0;
        for (int i = 0; i < weight; ++i) {
            if (vals[i] == DOLLAR) break;            // '$' contributes 0 and stops
            name += 1 + (uint64_t)vals[i] * g[i + 1];
        }
        return name;
    }
};

// Convenience: name a gapped k-mer directly from a text and a start position.
// `text` may be shorter than pos+span; positions past the end are treated as '$'.
inline uint64_t name_at(const Shape& shape, const std::string& text, size_t pos) {
    std::vector<int8_t> vals(shape.weight);
    for (int i = 0; i < shape.weight; ++i) {
        size_t idx = pos + shape.care[i];
        char c = (idx < text.size()) ? text[idx] : '$';
        vals[i] = base_value(c);
    }
    return shape.name_from_values(vals.data());
}

// Human-readable spelling of a gapped k-mer name (inverse of the naming), using
// '$' for terminator positions.  Useful for debugging / printing tables.
inline std::string name_to_string(const Shape& shape, uint64_t name) {
    std::string out;
    for (int i = 0; i < shape.weight; ++i) {
        if (name == 0) { out.push_back('$'); continue; }
        // At position i, value 0 is '$' (contributes 0). Values 1..sigma map to
        // bases with contribution 1 + (v-1)*g[i+1]. Recover v.
        uint64_t block = shape.g[i + 1];
        // name = 1 + base*block + rest, base in 0..sigma-1
        uint64_t t = name - 1;
        uint64_t base = t / block;
        out.push_back(base_char((int)base));
        name = t - base * block;
    }
    return out;
}

} // namespace gcsa
