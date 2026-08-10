// compare_algos.cpp
//
// Side-by-side comparison of CompressAlgo::Greedy vs CompressAlgo::DepOrder
// on the note's #.# example and on the repetition suite (weight-30 shapes).
//
// Build:  make compare_algos
// Usage:  ./compare_algos [--min-rep N] [--max-rep N] [--step S] [--seed S]
//                         [--motif-len L] [--max-add A]

#include "compress.hpp"
#include "serialize.hpp"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

using namespace gcsa;
using Clock = std::chrono::steady_clock;

static std::string random_dna(size_t len, std::mt19937_64& rng) {
    static const char bases[] = "ACGT";
    std::uniform_int_distribution<int> d(0, 3);
    std::string s(len, 'A');
    for (size_t i = 0; i < len; ++i) s[i] = bases[d(rng)];
    return s;
}

static std::string concat_copies(const std::string& motif, int reps) {
    std::string t;
    t.reserve(motif.size() * (size_t)reps);
    for (int i = 0; i < reps; ++i) t += motif;
    return t;
}

static std::string cares(int n) { return std::string((size_t)n, '#'); }
static std::string gaps(int n)  { return std::string((size_t)n, '.'); }

struct ShapeSpec { std::string pattern, label; };

static std::vector<ShapeSpec> make_shapes() {
    std::vector<ShapeSpec> v;
    v.push_back({cares(30), "30cont"});
    for (int g : {1, 2, 4, 8})
        v.push_back({cares(20) + gaps(g) + cares(10), "20+g" + std::to_string(g) + "+10"});
    for (int g : {1, 2, 4, 8})
        v.push_back({cares(10) + gaps(g) + cares(20), "10+g" + std::to_string(g) + "+20"});
    for (int g : {1, 2, 4})
        v.push_back({cares(5) + gaps(g) + cares(20) + gaps(g) + cares(5),
                     "5+g" + std::to_string(g) + "+20+g" + std::to_string(g) + "+5"});
    {
        std::string left, right;
        for (int i = 0; i < 5; ++i) { left += "#."; right += ".#"; }
        v.push_back({left + cares(20) + right, "sparse5+20+sparse5"});
    }
    return v;
}

static std::map<uint64_t, std::vector<int64_t>>
brute_positions(const Shape& sh, const std::string& text) {
    std::map<uint64_t, std::vector<int64_t>> m;
    for (long p = 0; p <= (long)text.size(); ++p)
        m[name_at(sh, text, (size_t)p)].push_back(p);
    return m;
}

struct Result {
    size_t m = 0, C = 0, bytes = 0;
    double keep_pct = 0, size_pct = 0, build_ms = 0;
    bool ok = false, rt = false;
};

static Result run_algo(const Shape& sh, const std::string& text,
                       int max_add, CompressAlgo algo) {
    Result r;
    auto t0 = Clock::now();
    CompressedIndex idx;
    idx.build(sh, text, max_add, algo);
    r.build_ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    r.m = idx.total_positions();
    r.C = idx.stored_positions();
    r.keep_pct = 100.0 * r.C / std::max<size_t>(1, r.m);

    auto truth = brute_positions(sh, text);
    r.ok = true;
    for (auto& kv : truth) {
        auto got = idx.positions_of(kv.first);
        std::sort(got.begin(), got.end());
        auto exp = kv.second; std::sort(exp.begin(), exp.end());
        if (got != exp) { r.ok = false; break; }
    }
    SerializedIndex S = serialize_index(idx);
    r.rt = true;
    for (auto& kv : truth) {
        auto got = S.positions_of(kv.first);
        std::sort(got.begin(), got.end());
        auto exp = kv.second; std::sort(exp.begin(), exp.end());
        if (got != exp) { r.rt = false; break; }
    }
    r.bytes = S.total_bytes();
    r.size_pct = 100.0 * r.bytes / std::max<size_t>(1, r.m * (size_t)S.pb);
    return r;
}

static void print_pair(const std::string& tag, const Result& g, const Result& d) {
    double dC = (g.C == 0) ? 0.0 : 100.0 * (1.0 - (double)d.C / (double)g.C);
    std::cout << std::left << std::setw(28) << tag
              << std::setw(8) << g.C
              << std::setw(8) << std::fixed << std::setprecision(1) << g.keep_pct
              << std::setw(8) << d.C
              << std::setw(8) << d.keep_pct
              << std::setw(9) << std::setprecision(1) << dC
              << std::setw(6) << ((g.ok && g.rt) ? "ok" : "FAIL")
              << std::setw(6) << ((d.ok && d.rt) ? "ok" : "FAIL")
              << std::setw(8) << std::setprecision(1) << g.build_ms
              << std::setw(8) << d.build_ms
              << "\n";
}

int main(int argc, char** argv) {
    uint64_t seed = 1;
    size_t motif_len = 200;
    int min_rep = 10, max_rep = 100, step = 10;
    int max_add = 8;

    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto need = [&](const char* n) {
            if (++i >= argc) { std::cerr << "missing " << n << "\n"; std::exit(1); }
            return std::string(argv[i]);
        };
        if      (k == "--seed")      seed = std::stoull(need("--seed"));
        else if (k == "--motif-len") motif_len = std::stoul(need("--motif-len"));
        else if (k == "--min-rep")   min_rep = std::stoi(need("--min-rep"));
        else if (k == "--max-rep")   max_rep = std::stoi(need("--max-rep"));
        else if (k == "--step")      step = std::stoi(need("--step"));
        else if (k == "--max-add")   max_add = std::stoi(need("--max-add"));
        else {
            std::cerr << "Usage: " << argv[0]
                      << " [--seed S] [--motif-len L] [--min-rep N] [--max-rep N]\n"
                      << "       [--step S] [--max-add A]\n";
            return 1;
        }
    }

    std::cout << "Comparing greedy vs dep-order (max_add=" << max_add << ")\n\n";
    std::cout << std::left
              << std::setw(28) << "case"
              << std::setw(8)  << "g|C|"
              << std::setw(8)  << "g%"
              << std::setw(8)  << "d|C|"
              << std::setw(8)  << "d%"
              << std::setw(9)  << "dC%↓"
              << std::setw(6)  << "g"
              << std::setw(6)  << "d"
              << std::setw(8)  << "g_ms"
              << std::setw(8)  << "d_ms"
              << "\n";

    int fails = 0;
    int dep_better = 0, greedy_better = 0, tie = 0;
    long long sum_gC = 0, sum_dC = 0;

    // --- Note's #.# example ------------------------------------------------
    {
        std::string text = "ACGTCTTAAACCCTCGTCTTAAACCCAACGTCTTAAACCC";
        Shape sh = Shape::parse("#.#");
        Result g = run_algo(sh, text, max_add, CompressAlgo::Greedy);
        Result d = run_algo(sh, text, max_add, CompressAlgo::DepOrder);
        print_pair("note #.# (n=40)", g, d);
        if (!g.ok || !g.rt || !d.ok || !d.rt) ++fails;
        sum_gC += (long long)g.C; sum_dC += (long long)d.C;
        if (d.C < g.C) ++dep_better; else if (g.C < d.C) ++greedy_better; else ++tie;
    }

    // --- A few small shapes on a short repetitive text ---------------------
    {
        std::mt19937_64 rng(seed);
        std::string motif = random_dna(200, rng);
        std::string text = concat_copies(motif, 20);
        for (const char* pat : {"#####", "##.##", "####.####", "#.#.#.#"}) {
            Shape sh = Shape::parse(pat);
            Result g = run_algo(sh, text, max_add, CompressAlgo::Greedy);
            Result d = run_algo(sh, text, max_add, CompressAlgo::DepOrder);
            print_pair(std::string("rep20x200 ") + pat, g, d);
            if (!g.ok || !g.rt || !d.ok || !d.rt) ++fails;
            sum_gC += (long long)g.C; sum_dC += (long long)d.C;
            if (d.C < g.C) ++dep_better; else if (g.C < d.C) ++greedy_better; else ++tie;
        }
    }

    std::cout << "\n";

    // --- Full weight-30 repetition suite -----------------------------------
    std::mt19937_64 rng(seed);
    std::string motif = random_dna(motif_len, rng);
    auto shapes = make_shapes();
    std::cout << "Weight-30 suite: motif_len=" << motif_len
              << " reps=" << min_rep << ".." << max_rep << " step " << step
              << " shapes=" << shapes.size() << "\n";

    for (int reps = min_rep; reps <= max_rep; reps += step) {
        std::string text = concat_copies(motif, reps);
        for (auto& spec : shapes) {
            Shape sh = Shape::parse(spec.pattern);
            Result g = run_algo(sh, text, max_add, CompressAlgo::Greedy);
            Result d = run_algo(sh, text, max_add, CompressAlgo::DepOrder);
            std::string tag = "x" + std::to_string(reps) + " " + spec.label;
            print_pair(tag, g, d);
            std::cout.flush();
            if (!g.ok || !g.rt || !d.ok || !d.rt) ++fails;
            sum_gC += (long long)g.C; sum_dC += (long long)d.C;
            if (d.C < g.C) ++dep_better;
            else if (g.C < d.C) ++greedy_better;
            else ++tie;
        }
    }

    int total = dep_better + greedy_better + tie;
    double overall = (sum_gC == 0) ? 0.0 : 100.0 * (1.0 - (double)sum_dC / (double)sum_gC);
    std::cout << "\n=== Summary ===\n"
              << "configs: " << total
              << "  dep-order better: " << dep_better
              << "  greedy better: " << greedy_better
              << "  tie: " << tie << "\n"
              << "total |C| greedy=" << sum_gC
              << "  dep-order=" << sum_dC
              << "  overall reduction: " << std::fixed << std::setprecision(2)
              << overall << "%\n"
              << "correctness failures: " << fails << "\n";
    return fails ? 1 : 0;
}
