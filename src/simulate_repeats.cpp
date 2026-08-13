// simulate_repeats.cpp
//
// Write a DNA FASTA consisting of x identical concatenated copies of a random
// motif of length y.
//
// Usage:
//   ./simulate_repeats -x <reps> -y <motif_len> [-o out.fasta] [--seed S]
//                      [--header NAME] [--line-width W]
//
// Example:
//   ./simulate_repeats -x 50 -y 200 -o repeats.fasta
//   ./simulate_repeats -x 100 -y 200 --seed 1 -o /tmp/rep.fasta

#include "random_dna.hpp"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " -x <reps> -y <motif_len> [-o out.fasta] [--seed S]\n"
        << "       [--header NAME] [--line-width W]\n"
        << "\n"
        << "  -x            number of concatenated motif copies (required, >= 1)\n"
        << "  -y            motif length in bp (required, >= 1)\n"
        << "  -o            output FASTA path (default: stdout)\n"
        << "  --seed S      RNG seed (default: random_device)\n"
        << "  --header NAME FASTA header without '>' (default: repeats_x{X}_y{Y})\n"
        << "  --line-width W  wrap sequence every W columns; 0 = no wrap (default: 60)\n";
}

static void write_fasta(std::ostream& out, const std::string& header,
                        const std::string& seq, int line_width) {
    out << '>' << header << '\n';
    if (line_width <= 0) {
        out << seq << '\n';
        return;
    }
    for (size_t i = 0; i < seq.size(); i += (size_t)line_width)
        out << seq.substr(i, (size_t)line_width) << '\n';
}

int main(int argc, char** argv) {
    long reps = -1;
    long motif_len = -1;
    std::string out_path;
    std::string header;
    int line_width = 60;
    bool have_seed = false;
    uint64_t seed = 0;

    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto need = [&](const char* n) -> std::string {
            if (++i >= argc) {
                std::cerr << "Missing value for " << n << "\n";
                usage(argv[0]);
                std::exit(1);
            }
            return argv[i];
        };
        if (k == "-x")              reps = std::stol(need("-x"));
        else if (k == "-y")         motif_len = std::stol(need("-y"));
        else if (k == "-o")         out_path = need("-o");
        else if (k == "--seed")   { seed = std::stoull(need("--seed")); have_seed = true; }
        else if (k == "--header")   header = need("--header");
        else if (k == "--line-width") line_width = std::stoi(need("--line-width"));
        else if (k == "-h" || k == "--help") { usage(argv[0]); return 0; }
        else {
            std::cerr << "Unknown argument: " << k << "\n";
            usage(argv[0]);
            return 1;
        }
    }

    if (reps < 1 || motif_len < 1) {
        std::cerr << "Both -x (>=1) and -y (>=1) are required.\n";
        usage(argv[0]);
        return 1;
    }

    if (!have_seed) seed = std::random_device{}();
    std::mt19937_64 rng(seed);

    std::string motif = gcsa::random_dna((size_t)motif_len, rng);
    std::string seq;
    seq.reserve((size_t)reps * (size_t)motif_len);
    for (long i = 0; i < reps; ++i) seq += motif;

    if (header.empty())
        header = "repeats_x" + std::to_string(reps) + "_y" + std::to_string(motif_len)
                 + " motif_len=" + std::to_string(motif_len)
                 + " reps=" + std::to_string(reps)
                 + " total=" + std::to_string(seq.size())
                 + " seed=" + std::to_string(seed);

    if (out_path.empty()) {
        write_fasta(std::cout, header, seq, line_width);
    } else {
        std::ofstream out(out_path);
        if (!out) {
            std::cerr << "Cannot open '" << out_path << "' for writing\n";
            return 1;
        }
        write_fasta(out, header, seq, line_width);
        std::cerr << "Wrote " << seq.size() << " bp ("
                  << reps << " x " << motif_len
                  << ") to " << out_path
                  << "  seed=" << seed << "\n";
    }
    return 0;
}
