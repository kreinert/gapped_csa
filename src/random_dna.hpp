// random_dna.hpp
//
// Shared helper for generating a random DNA string. Pulled out of
// bench_repetition.cpp, compare_algos.cpp, simulate_repeats.cpp, and
// find_minimal_gap.cpp, which had each grown their own (identical, modulo
// the parameter type) copy of this function.

#pragma once

#include <random>
#include <string>

namespace gcsa {

inline std::string random_dna(size_t len, std::mt19937_64& rng) {
    static const char bases[] = "ACGT";
    std::uniform_int_distribution<int> d(0, 3);
    std::string s(len, 'A');
    for (size_t i = 0; i < len; ++i) s[i] = bases[d(rng)];
    return s;
}

}  // namespace gcsa
