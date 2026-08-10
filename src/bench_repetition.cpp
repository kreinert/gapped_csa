// bench_repetition.cpp
//
// Standalone benchmark for the compressed gapped SA on highly repetitive DNA:
//   - draw a random motif of length 200
//   - concatenate it R times for R in {10,20,...,100} (configurable)
//   - for various shapes of weight 30 that each contain a run of 20 consecutive
//     care positions ('#'), measure:
//       * correctness (positions_of vs brute force; serialized round-trip)
//       * compression (|C|/m and serialized bytes vs full SA)
//       * build / query wall time
//
// Build:  make bench_repetition
// Usage:  ./bench_repetition [--seed S] [--motif-len L] [--min-rep N] [--max-rep N]
//                            [--step S] [--max-add A] [--queries Q]

#include "compress.hpp"
#include "serialize.hpp"

#include <chrono>
#include <cmath>
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

// Shape generators: weight 30, containing a contiguous run of 20 '#'.
// Remaining 10 cares are placed with the given gap layout.
struct ShapeSpec {
    std::string pattern;
    std::string label;   // short description
};

static std::string cares(int n) { return std::string((size_t)n, '#'); }
static std::string gaps(int n)  { return std::string((size_t)n, '.'); }

static std::vector<ShapeSpec> make_shapes() {
    // All have weight 30 and contain "####################" (20 cares).
    std::vector<ShapeSpec> v;
    // Fully contiguous weight-30 k-mer (contains 20 consecutive).
    v.push_back({cares(30), "30cont"});
    // 20-run then 10 cares, separated by g = 1,2,4,8.
    for (int g : {1, 2, 4, 8})
        v.push_back({cares(20) + gaps(g) + cares(10), "20+g" + std::to_string(g) + "+10"});
    // 10 cares, gap, then 20-run.
    for (int g : {1, 2, 4, 8})
        v.push_back({cares(10) + gaps(g) + cares(20), "10+g" + std::to_string(g) + "+20"});
    // 5 + gap + 20-run + gap + 5, symmetric.
    for (int g : {1, 2, 4})
        v.push_back({cares(5) + gaps(g) + cares(20) + gaps(g) + cares(5),
                     "5+g" + std::to_string(g) + "+20+g" + std::to_string(g) + "+5"});
    // 20-run flanked by sparse cares: (#.){5} + #{20} + (.#){5} → weight 30.
    {
        std::string left, right;
        for (int i = 0; i < 5; ++i) { left += "#."; right += ".#"; }
        v.push_back({left + cares(20) + right, "sparse5+20+sparse5"});
    }
    // Verify invariants.
    for (auto& s : v) {
        int w = 0, run = 0, maxrun = 0;
        for (char c : s.pattern) {
            if (c == '#') { ++w; ++run; maxrun = std::max(maxrun, run); }
            else run = 0;
        }
        if (w != 30 || maxrun < 20)
            throw std::runtime_error("bad shape " + s.label + ": w=" + std::to_string(w)
                                     + " maxrun=" + std::to_string(maxrun));
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

struct BenchRow {
    int reps = 0;
    size_t n = 0;
    std::string label;
    std::string pattern;
    int span = 0;
    int weight = 0;
    size_t m = 0;
    size_t C_size = 0;
    size_t distinct = 0;
    double keep_pct = 0;      // 100 * |C|/m
    size_t bytes_full_sa = 0;
    size_t bytes_total = 0;
    double size_pct = 0;      // 100 * total / full_sa
    bool correct = false;
    bool roundtrip = false;
    double build_ms = 0;
    double query_ms = 0;
    size_t query_hits = 0;
};

static BenchRow run_one(const std::string& text, int reps,
                        const ShapeSpec& spec, int max_add, int nqueries,
                        std::mt19937_64& rng) {
    BenchRow r;
    r.reps = reps;
    r.n = text.size();
    r.label = spec.label;
    r.pattern = spec.pattern;

    Shape sh = Shape::parse(spec.pattern);
    r.span = sh.span;
    r.weight = 30;

    auto t0 = Clock::now();
    CompressedIndex idx;
    idx.build(sh, text, max_add);
    auto t1 = Clock::now();
    r.build_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    r.m = idx.total_positions();
    r.C_size = idx.stored_positions();
    r.distinct = idx.num_kmers();
    r.keep_pct = 100.0 * r.C_size / std::max<size_t>(1, r.m);

    // Correctness vs brute force.
    auto truth = brute_positions(sh, text);
    r.correct = true;
    for (auto& kv : truth) {
        auto got = idx.positions_of(kv.first);
        std::sort(got.begin(), got.end());
        auto exp = kv.second;
        std::sort(exp.begin(), exp.end());
        if (got != exp) { r.correct = false; break; }
    }

    SerializedIndex S = serialize_index(idx);
    r.roundtrip = true;
    for (auto& kv : truth) {
        auto got = S.positions_of(kv.first);
        std::sort(got.begin(), got.end());
        auto exp = kv.second;
        std::sort(exp.begin(), exp.end());
        if (got != exp) { r.roundtrip = false; break; }
    }
    r.bytes_full_sa = r.m * (size_t)S.pb;
    r.bytes_total = S.total_bytes();
    r.size_pct = 100.0 * r.bytes_total / std::max<size_t>(1, r.bytes_full_sa);

    // Query: random windows of length = span drawn from the text.
    r.query_hits = 0;
    std::uniform_int_distribution<size_t> pos_dist(
        0, text.size() > (size_t)sh.span ? text.size() - (size_t)sh.span : 0);
    auto tq0 = Clock::now();
    for (int q = 0; q < nqueries; ++q) {
        size_t p = pos_dist(rng);
        std::string qstr = text.substr(p, (size_t)sh.span);
        auto hits = idx.locate(qstr);
        r.query_hits += hits.size();
    }
    auto tq1 = Clock::now();
    r.query_ms = std::chrono::duration<double, std::milli>(tq1 - tq0).count();

    return r;
}

static void print_header() {
    std::cout << std::left
              << std::setw(5)  << "reps"
              << std::setw(7)  << "n"
              << std::setw(22) << "shape"
              << std::setw(6)  << "span"
              << std::setw(8)  << "m"
              << std::setw(8)  << "|C|"
              << std::setw(8)  << "keep%"
              << std::setw(8)  << "kmers"
              << std::setw(10) << "bytes"
              << std::setw(8)  << "size%"
              << std::setw(6)  << "ok"
              << std::setw(6)  << "rt"
              << std::setw(10) << "build_ms"
              << std::setw(10) << "query_ms"
              << "\n";
}

static void print_row(const BenchRow& r) {
    std::cout << std::left
              << std::setw(5)  << r.reps
              << std::setw(7)  << r.n
              << std::setw(22) << r.label
              << std::setw(6)  << r.span
              << std::setw(8)  << r.m
              << std::setw(8)  << r.C_size
              << std::setw(8)  << std::fixed << std::setprecision(1) << r.keep_pct
              << std::setw(8)  << r.distinct
              << std::setw(10) << r.bytes_total
              << std::setw(8)  << std::fixed << std::setprecision(1) << r.size_pct
              << std::setw(6)  << (r.correct ? "PASS" : "FAIL")
              << std::setw(6)  << (r.roundtrip ? "PASS" : "FAIL")
              << std::setw(10) << std::fixed << std::setprecision(1) << r.build_ms
              << std::setw(10) << std::fixed << std::setprecision(1) << r.query_ms
              << "\n";
}

int main(int argc, char** argv) {
    uint64_t seed = 1;
    size_t motif_len = 200;
    int min_rep = 10, max_rep = 100, step = 10;
    int max_add = 8;
    int nqueries = 100;

    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto need = [&](const char* n) {
            if (++i >= argc) { std::cerr << "missing value for " << n << "\n"; std::exit(1); }
            return std::string(argv[i]);
        };
        if      (k == "--seed")      seed = std::stoull(need("--seed"));
        else if (k == "--motif-len") motif_len = std::stoul(need("--motif-len"));
        else if (k == "--min-rep")   min_rep = std::stoi(need("--min-rep"));
        else if (k == "--max-rep")   max_rep = std::stoi(need("--max-rep"));
        else if (k == "--step")      step = std::stoi(need("--step"));
        else if (k == "--max-add")   max_add = std::stoi(need("--max-add"));
        else if (k == "--queries")   nqueries = std::stoi(need("--queries"));
        else {
            std::cerr << "Usage: " << argv[0]
                      << " [--seed S] [--motif-len L] [--min-rep N] [--max-rep N]\n"
                      << "       [--step S] [--max-add A] [--queries Q]\n";
            return 1;
        }
    }

    std::mt19937_64 rng(seed);
    std::string motif = random_dna(motif_len, rng);
    auto shapes = make_shapes();

    std::cout << "motif length = " << motif_len
              << "  seed = " << seed
              << "  reps = " << min_rep << ".." << max_rep << " step " << step
              << "  shapes = " << shapes.size()
              << " (weight 30, each with ≥20 consecutive #)\n";
    std::cout << "shapes:\n";
    for (auto& s : shapes)
        std::cout << "  " << std::setw(22) << std::left << s.label
                  << " span=" << s.pattern.size() << "  " << s.pattern << "\n";
    std::cout << "\n";

    print_header();

    int fails = 0;
    std::vector<BenchRow> rows;
    for (int reps = min_rep; reps <= max_rep; reps += step) {
        std::string text = concat_copies(motif, reps);
        for (auto& spec : shapes) {
            BenchRow r = run_one(text, reps, spec, max_add, nqueries, rng);
            print_row(r);
            std::cout.flush();
            if (!r.correct || !r.roundtrip) ++fails;
            rows.push_back(std::move(r));
        }
        std::cout << "\n";
    }

    // Summary: mean keep% and size% by shape (over all reps), and overall.
    std::cout << "=== Summary by shape (mean over reps) ===\n";
    std::cout << std::left
              << std::setw(22) << "shape"
              << std::setw(10) << "mean_keep%"
              << std::setw(10) << "mean_size%"
              << std::setw(10) << "mean_build"
              << "all_ok\n";
    for (auto& spec : shapes) {
        double sk = 0, ss = 0, sb = 0;
        int cnt = 0;
        bool all_ok = true;
        for (auto& r : rows) if (r.label == spec.label) {
            sk += r.keep_pct; ss += r.size_pct; sb += r.build_ms; ++cnt;
            if (!r.correct || !r.roundtrip) all_ok = false;
        }
        if (!cnt) continue;
        std::cout << std::setw(22) << spec.label
                  << std::setw(10) << std::fixed << std::setprecision(1) << (sk / cnt)
                  << std::setw(10) << (ss / cnt)
                  << std::setw(10) << (sb / cnt)
                  << (all_ok ? "PASS" : "FAIL") << "\n";
    }

    std::cout << "\nOverall: " << rows.size() << " configs, "
              << fails << " failures.\n";
    return fails ? 1 : 0;
}
