// simulate_random.cpp
//
// Write a uniform-random DNA FASTA. This is the "category A" benchmark
// input from the benchmark design doc: no repeat structure at all, so
// keep% should land near 100% and the gzip ratio near 1 -- the sanity-check
// row that everything else is compared against.
//
// Usage:
//   ./simulate_random -n <len> [-o out.fasta] [--seed S]
//                      [--header NAME] [--line-width W]
//
// Example:
//   ./simulate_random -n 1000000 --seed 1 -o random_1e6.fasta

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
        << " -n <len> [-o out.fasta] [--seed S] [--header NAME] [--line-width W]\n"
        << "\n"
        << "  -n            sequence length in bp (required, >= 1)\n"
        << "  -o            output FASTA path (default: stdout)\n"
        << "  --seed S      RNG seed (default: random_device)\n"
        << "  --header NAME FASTA header without '>' (default: random_n{N}_seed{S})\n"
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
    long len = -1;
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
        if (k == "-n")               len = std::stol(need("-n"));
        else if (k == "-o")          out_path = need("-o");
        else if (k == "--seed")    { seed = std::stoull(need("--seed")); have_seed = true; }
        else if (k == "--header")    header = need("--header");
        else if (k == "--line-width") line_width = std::stoi(need("--line-width"));
        else if (k == "-h" || k == "--help") { usage(argv[0]); return 0; }
        else {
            std::cerr << "Unknown argument: " << k << "\n";
            usage(argv[0]);
            return 1;
        }
    }

    if (len < 1) {
        std::cerr << "-n (>=1) is required.\n";
        usage(argv[0]);
        return 1;
    }

    if (!have_seed) seed = std::random_device{}();
    std::mt19937_64 rng(seed);

    std::string seq = gcsa::random_dna((size_t)len, rng);

    if (header.empty())
        header = "random_n" + std::to_string(len) + " len=" + std::to_string(len)
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
        std::cerr << "Wrote " << seq.size() << " bp to " << out_path
                  << "  seed=" << seed << "\n";
    }
    return 0;
}
