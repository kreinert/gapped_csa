// simulate_pangenome.cpp
//
// Write a DNA FASTA consisting of `n` concatenated near-identical genome
// copies: a reference (either loaded from a FASTA or freshly randomized),
// followed by n-1 point-mutated copies at a chosen divergence rate. This is
// "category E" in the benchmark design doc -- the email's example is 1, 2,
// 4, 8 (near-)identical E. coli strains; divergence is the second knob kept
// separate and fixed so strain count is the only thing varying in a sweep.
//
// Usage:
//   ./simulate_pangenome -r <ref.fasta> -n <strains> --divergence D
//                         [-o out.fasta] [--seed S] [--header NAME]
//                         [--line-width W]
//   ./simulate_pangenome -y <ref_len> -n <strains> --divergence D  # no
//                         reference file: generate a random one of length y
//
// Divergence D is the per-base substitution probability applied
// independently to each of the n-1 non-reference copies (so all copies
// differ from the reference at the same expected rate, but not from each
// other in a correlated way -- this matches "strains diverged from a common
// ancestor" better than "each strain diverged from the previous one").
//
// Example:
//   ./simulate_pangenome -r ecoli_k12.fasta -n 4 --divergence 0.01 --seed 1 -o pangenome_ecoli_n4.fasta
//   ./simulate_pangenome -y 200000 -n 8 --divergence 0.01 --seed 1 -o /tmp/pg8.fasta

#include "random_dna.hpp"

#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " (-r <ref.fasta> | -y <ref_len>) -n <strains> --divergence D\n"
        << "       [-o out.fasta] [--seed S] [--header NAME] [--line-width W]\n"
        << "\n"
        << "  -r <ref.fasta>   reference genome to diverge from (mutually\n"
        << "                   exclusive with -y)\n"
        << "  -y <ref_len>     generate a random reference of this length instead\n"
        << "                   of loading one (mutually exclusive with -r)\n"
        << "  -n <strains>     total genome copies in the output, including the\n"
        << "                   reference itself (required, >= 1)\n"
        << "  --divergence D   per-base substitution probability for each of the\n"
        << "                   n-1 non-reference copies, 0 <= D <= 1 (required)\n"
        << "  -o               output FASTA path (default: stdout)\n"
        << "  --seed S         RNG seed (default: random_device)\n"
        << "  --header NAME    FASTA header without '>' (default: derived)\n"
        << "  --line-width W   wrap sequence every W columns; 0 = no wrap (default: 60)\n";
}

static std::string load_dna(const std::string& path) {
    std::ifstream in(path);
    if (!in) { std::cerr << "cannot open " << path << "\n"; std::exit(1); }
    std::string out, line;
    static const char bases[] = "ACGT";
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') continue;
        for (char c : line) {
            if (!std::isalpha((unsigned char)c)) continue;
            char u = (char)std::toupper((unsigned char)c);
            out.push_back((u == 'A' || u == 'C' || u == 'G' || u == 'T') ? u : bases[0]);
        }
    }
    return out;
}

// Mutate `n` bases of `s`, uniformly at random, each to one of the other 3
// bases (never a silent self-substitution, so the expected divergence is
// exactly the requested rate rather than 3/4 of it).
static std::string mutate(const std::string& s, double rate, std::mt19937_64& rng) {
    static const char bases[] = "ACGT";
    std::string out = s;
    std::binomial_distribution<size_t> nmut_dist(s.size(), rate);
    size_t nmut = nmut_dist(rng);
    std::uniform_int_distribution<size_t> pos_dist(0, s.empty() ? 0 : s.size() - 1);
    std::uniform_int_distribution<int> base_dist(0, 2);
    // Positions are drawn with replacement, so on rare collisions a position
    // gets mutated more than once; the realized divergence is then a hair
    // below the requested rate. Fine for a benchmark knob, not exact stats.
    for (size_t i = 0; i < nmut && !s.empty(); ++i) {
        size_t p = pos_dist(rng);
        char cur = out[p];
        char alts[3]; int k = 0;
        for (char b : bases) if (b != cur) alts[k++] = b;
        out[p] = alts[base_dist(rng)];
    }
    return out;
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
    std::string ref_path, out_path, header;
    long ref_len = -1;
    long n_strains = -1;
    double divergence = -1.0;
    int line_width = 60;
    bool have_seed = false;
    uint64_t seed = 0;

    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto need = [&](const char* n) -> std::string {
            if (++i >= argc) { std::cerr << "Missing value for " << n << "\n"; usage(argv[0]); std::exit(1); }
            return argv[i];
        };
        if (k == "-r")                    ref_path = need("-r");
        else if (k == "-y")                ref_len = std::stol(need("-y"));
        else if (k == "-n")                n_strains = std::stol(need("-n"));
        else if (k == "--divergence")      divergence = std::stod(need("--divergence"));
        else if (k == "-o")                out_path = need("-o");
        else if (k == "--seed")          { seed = std::stoull(need("--seed")); have_seed = true; }
        else if (k == "--header")          header = need("--header");
        else if (k == "--line-width")      line_width = std::stoi(need("--line-width"));
        else if (k == "-h" || k == "--help") { usage(argv[0]); return 0; }
        else { std::cerr << "Unknown argument: " << k << "\n"; usage(argv[0]); return 1; }
    }

    if (n_strains < 1) { std::cerr << "-n (>=1) is required.\n"; usage(argv[0]); return 1; }
    if (divergence < 0.0 || divergence > 1.0) {
        std::cerr << "--divergence is required, 0 <= D <= 1.\n"; usage(argv[0]); return 1;
    }
    if (ref_path.empty() == (ref_len < 0)) {
        // Both or neither given; exactly one of -r / -y is required.
        // (ref_path.empty() && ref_len<0) -> neither; (!ref_path.empty() && ref_len>=0) -> both.
        std::cerr << "Exactly one of -r <ref.fasta> or -y <ref_len> is required.\n";
        usage(argv[0]);
        return 1;
    }

    if (!have_seed) seed = std::random_device{}();
    std::mt19937_64 rng(seed);

    std::string ref = ref_path.empty() ? gcsa::random_dna((size_t)ref_len, rng)
                                        : load_dna(ref_path);
    if (ref.empty()) { std::cerr << "reference sequence is empty\n"; return 1; }

    std::string seq;
    seq.reserve(ref.size() * (size_t)n_strains);
    seq += ref;
    for (long s = 1; s < n_strains; ++s)
        seq += mutate(ref, divergence, rng);

    if (header.empty())
        header = "pangenome_n" + std::to_string(n_strains)
                 + " ref_len=" + std::to_string(ref.size())
                 + " divergence=" + std::to_string(divergence)
                 + " strains=" + std::to_string(n_strains)
                 + " total=" + std::to_string(seq.size())
                 + " seed=" + std::to_string(seed)
                 + (ref_path.empty() ? " ref=random" : (" ref=" + ref_path));

    if (out_path.empty()) {
        write_fasta(std::cout, header, seq, line_width);
    } else {
        std::ofstream out(out_path);
        if (!out) { std::cerr << "Cannot open '" << out_path << "' for writing\n"; return 1; }
        write_fasta(out, header, seq, line_width);
        std::cerr << "Wrote " << seq.size() << " bp (" << n_strains << " strains x "
                  << ref.size() << " bp, divergence=" << divergence << ") to "
                  << out_path << "  seed=" << seed << "\n";
    }
    return 0;
}
