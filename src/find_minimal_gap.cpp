// find_minimal_gap.cpp
// Search for a short DNA string + shape where dep-order |C| < greedy |C|.

#include "compress.hpp"
#include "random_dna.hpp"
#include <algorithm>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

using namespace gcsa;

static size_t stored(const std::string& text, const std::string& shape,
                     CompressAlgo algo, int max_add = 8) {
    CompressedIndex idx;
    idx.build(Shape::parse(shape), text, max_add, algo);
    return idx.stored_positions();
}

static bool correct(const std::string& text, const std::string& shape, CompressAlgo algo) {
    Shape sh = Shape::parse(shape);
    CompressedIndex idx;
    idx.build(sh, text, 8, algo);
    for (long p = 0; p <= (long)text.size(); ++p) {
        uint64_t nm = name_at(sh, text, (size_t)p);
        // spot-check later via full brute in main when found
        (void)nm;
    }
    // full check
    std::map<uint64_t, std::vector<int64_t>> truth;
    for (long p = 0; p <= (long)text.size(); ++p)
        truth[name_at(sh, text, (size_t)p)].push_back(p);
    for (auto& kv : truth) {
        auto got = idx.positions_of(kv.first);
        std::sort(got.begin(), got.end());
        auto exp = kv.second; std::sort(exp.begin(), exp.end());
        if (got != exp) return false;
    }
    return true;
}

int main() {
    std::vector<std::string> shapes = {
        "#.#", "##.#", "#.##", "#..#", "##.##", "###", "####", "#.#.#",
        "##", "#.#.", "###.#", "#.###"
    };

    std::mt19937_64 rng(42);
    int best_n = 1e9;
    std::string best_text, best_shape;
    size_t best_g = 0, best_d = 0;

    // 1) random short strings
    for (int n = 8; n <= 40; ++n) {
        for (int trial = 0; trial < 2000; ++trial) {
            std::string t = random_dna(n, rng);
            for (auto& sh : shapes) {
                if ((int)sh.size() > n) continue;
                size_t g = stored(t, sh, CompressAlgo::Greedy);
                size_t d = stored(t, sh, CompressAlgo::DepOrder);
                if (d < g && (int)t.size() < best_n) {
                    if (!correct(t, sh, CompressAlgo::Greedy)) continue;
                    if (!correct(t, sh, CompressAlgo::DepOrder)) continue;
                    best_n = (int)t.size();
                    best_text = t; best_shape = sh;
                    best_g = g; best_d = d;
                    std::cout << "FOUND n=" << best_n
                              << " shape=" << sh
                              << " greedy=" << g << " dep=" << d
                              << "  text=" << t << "\n";
                    std::cout.flush();
                }
            }
        }
        if (best_n <= 16) break; // good enough minimal
    }

    // 2) repetitive: motif * reps
    for (int mlen = 4; mlen <= 12 && best_n > 12; ++mlen) {
        for (int reps = 2; reps <= 8; ++reps) {
            for (int trial = 0; trial < 400; ++trial) {
                std::string motif = random_dna(mlen, rng);
                std::string t;
                for (int i = 0; i < reps; ++i) t += motif;
                for (auto& sh : shapes) {
                    if ((int)sh.size() > (int)t.size()) continue;
                    size_t g = stored(t, sh, CompressAlgo::Greedy);
                    size_t d = stored(t, sh, CompressAlgo::DepOrder);
                    if (d < g && (int)t.size() < best_n) {
                        if (!correct(t, sh, CompressAlgo::Greedy)) continue;
                        if (!correct(t, sh, CompressAlgo::DepOrder)) continue;
                        best_n = (int)t.size();
                        best_text = t; best_shape = sh;
                        best_g = g; best_d = d;
                        std::cout << "FOUND n=" << best_n
                                  << " shape=" << sh
                                  << " greedy=" << g << " dep=" << d
                                  << "  text=" << t
                                  << "  (motif=" << motif << " x" << reps << ")\n";
                        std::cout.flush();
                    }
                }
            }
        }
    }

    if (best_text.empty()) {
        std::cout << "No gap found.\n";
        return 1;
    }
    std::cout << "\n=== Minimal example ===\n"
              << "text  = " << best_text << "\n"
              << "shape = " << best_shape << "\n"
              << "greedy |C|    = " << best_g << "\n"
              << "dep-order |C| = " << best_d << "\n"
              << "saved         = " << (best_g - best_d) << "\n";
    return 0;
}
