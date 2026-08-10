// validate.cpp — reproduce the tables from the note for the #.# example and
// allow ad-hoc (text, shape) experiments.
//
// Usage:
//   ./validate                       # runs the note's #.# example
//   ./validate "<TEXT>" "<SHAPE>"    # custom text and shape

#include "gapped_sa.hpp"
#include <cstdio>
#include <iostream>
#include <iomanip>

using namespace gcsa;

static void print_table(const GappedSA& G, int max_suffix_syms = 12) {
    std::cout << "alphabet_size=" << G.alphabet_size
              << "  m=" << G.m() << "  name_space=" << G.shape.name_space() << "\n\n";
    std::cout << std::right
              << std::setw(4) << "Rk" << " | "
              << std::setw(6) << "lexSA" << " | "
              << std::setw(6) << "orig" << " | "
              << std::setw(6) << "code" << " | "
              << std::setw(6) << "kmer" << " | "
              << std::setw(4) << "lcp" << " | suffix(codes)\n";
    for (size_t r = 0; r < G.m(); ++r) {
        int32_t lexpos = G.sa[r];
        std::cout << std::setw(4) << r << " | "
                  << std::setw(6) << lexpos << " | "
                  << std::setw(6) << G.orig_pos((int32_t)r) << " | "
                  << std::setw(6) << G.first_symbol((int32_t)r) << " | "
                  << std::setw(6) << name_to_string(G.shape, G.first_symbol((int32_t)r)) << " | "
                  << std::setw(4) << G.lcp[r] << " | ";
        for (int k = 0; k < max_suffix_syms && lexpos + k < (int)G.m(); ++k)
            std::cout << G.lex[lexpos + k] << ' ';
        std::cout << "\n";
    }
}

int main(int argc, char** argv) {
    std::string text  = "aCGTCTTAAACCCtCGTCTTAAACCCaaCGTCTTAAACCC";
    std::string shape = "#.#";
    if (argc >= 2) text  = argv[1];
    if (argc >= 3) shape = argv[2];

    // Normalise text to ACGT.
    std::string clean;
    for (char c : text) if (std::isalpha((unsigned char)c)) clean.push_back(base_char(base_value(c)));

    Shape sh = Shape::parse(shape);
    std::cout << "text  = " << clean << "  (n=" << clean.size() << ")\n";
    std::cout << "shape = " << shape << "  span=" << sh.span
              << " weight=" << sh.weight << "\n\n";

    // Show the naming table for all names that occur.
    GappedSA G = build_gapped_sa(sh, clean);
    print_table(G);
    return 0;
}
