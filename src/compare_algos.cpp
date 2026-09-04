// compare_algos.cpp
//
// Side-by-side comparison of CompressAlgo::GreedySize, GreedyDegree, DepOrder,
// PseudoforestDp, PseudoforestDpIterate on the note's #.# example, the
// GCCTTTAAAG×3 demo, short repetitive DNA, and the weight-30 repetition
// suite.
//
// Build:  make compare_algos
// Usage:  ./compare_algos [--min-rep N] [--max-rep N] [--step S] [--seed S]
//                         [--motif-len L] [--max-add A] [--verbose]
//
// Per-algorithm stage/progress logging is suppressed by default so that only
// the table and the summary are printed; pass --verbose to keep it.

#include "compress.hpp"
#include "random_dna.hpp"
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
                       int max_add, CompressAlgo algo, int phase2_iters) {
    Result r;
    auto t0 = Clock::now();
    CompressedIndex idx;
    idx.build(sh, text, max_add, algo, phase2_iters);
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

static void print_header() {
    std::cout << std::left
              << std::setw(28) << "case"
              << std::setw(7)  << "g|C|"
              << std::setw(6)  << "g%"
              << std::setw(7)  << "gd|C|"
              << std::setw(6)  << "gd%"
              << std::setw(7)  << "d|C|"
              << std::setw(6)  << "d%"
              << std::setw(7)  << "p|C|"
              << std::setw(6)  << "p%"
              << std::setw(7)  << "pi|C|"
              << std::setw(6)  << "pi%"
              << std::setw(20) << "best"
              << std::setw(4)  << "g"
              << std::setw(4)  << "gd"
              << std::setw(4)  << "d"
              << std::setw(4)  << "p"
              << std::setw(4)  << "pi"
              << "\n";
}

// g=greedy-size, gd=greedy-degree, d=dep-order, p=pseudoforest-dp,
// pi=pseudoforest-dp-iterate
static void print_row(const std::string& tag,
                      const Result& g, const Result& gd, const Result& d,
                      const Result& p, const Result& pi) {
    size_t bC = std::min({g.C, gd.C, d.C, p.C, pi.C});
    std::string best;
    if (g.C == bC) best += (best.empty() ? "" : "/") + std::string("g");
    if (gd.C == bC) best += (best.empty() ? "" : "/") + std::string("gd");
    if (d.C == bC) best += (best.empty() ? "" : "/") + std::string("d");
    if (p.C == bC) best += (best.empty() ? "" : "/") + std::string("p");
    if (pi.C == bC) best += (best.empty() ? "" : "/") + std::string("pi");
    if (best == "g/gd/d/p/pi") best = "tie";

    auto ok = [](const Result& r) { return (r.ok && r.rt) ? "ok" : "FAIL"; };
    std::cout << std::left << std::setw(28) << tag
              << std::setw(7) << g.C
              << std::setw(6) << std::fixed << std::setprecision(1) << g.keep_pct
              << std::setw(7) << gd.C
              << std::setw(6) << gd.keep_pct
              << std::setw(7) << d.C
              << std::setw(6) << d.keep_pct
              << std::setw(7) << p.C
              << std::setw(6) << p.keep_pct
              << std::setw(7) << pi.C
              << std::setw(6) << pi.keep_pct
              << std::setw(20) << best
              << std::setw(4) << ok(g)
              << std::setw(4) << ok(gd)
              << std::setw(4) << ok(d)
              << std::setw(4) << ok(p)
              << std::setw(4) << ok(pi)
              << "\n";
}

struct Totals {
    long long sum_g = 0, sum_gd = 0, sum_d = 0, sum_p = 0, sum_pi = 0;
    int win_g = 0, win_gd = 0, win_d = 0, win_p = 0, win_pi = 0, tie = 0;
    int fails = 0, n = 0;

    void add(const Result& g, const Result& gd, const Result& d,
             const Result& p, const Result& pi) {
        ++n;
        sum_g += (long long)g.C;
        sum_gd += (long long)gd.C;
        sum_d += (long long)d.C;
        sum_p += (long long)p.C;
        sum_pi += (long long)pi.C;
        if (!g.ok || !g.rt || !gd.ok || !gd.rt || !d.ok || !d.rt ||
            !p.ok || !p.rt || !pi.ok || !pi.rt) ++fails;
        size_t b = std::min({g.C, gd.C, d.C, p.C, pi.C});
        int winners = (g.C == b) + (gd.C == b) + (d.C == b) + (p.C == b) + (pi.C == b);
        if (winners > 1) ++tie;
        else if (g.C == b) ++win_g;
        else if (gd.C == b) ++win_gd;
        else if (d.C == b) ++win_d;
        else if (p.C == b) ++win_p;
        else ++win_pi;
    }
};

int main(int argc, char** argv) {
    uint64_t seed = 1;
    size_t motif_len = 200;
    int min_rep = 10, max_rep = 100, step = 10;
    int max_add = 8;
    int phase2_iters = 0;  // 0 = unset: GCSA_PHASE2_MAX_ITERS, else the default
    bool verbose = false;

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
        else if (k == "--phase2-iters") {
            phase2_iters = std::stoi(need("--phase2-iters"));
            if (phase2_iters < 1) {
                std::cerr << "--phase2-iters must be >= 1\n";
                return 1;
            }
        }
        else if (k == "--verbose")   verbose = true;
        else {
            std::cerr << "Usage: " << argv[0]
                      << " [--seed S] [--motif-len L] [--min-rep N] [--max-rep N]\n"
                      << "       [--step S] [--max-add A] [--phase2-iters N] [--verbose]\n";
            return 1;
        }
    }

    gcsa_set_quiet(!verbose);

    std::cout << "Comparing greedy-size | greedy-degree | dep-order | pseudoforest-dp | pseudoforest-dp-iterate"
              << "  (max_add=" << max_add << ")\n"
              << "columns: g=greedy-size  gd=greedy-degree  d=dep-order"
              << "  p=pseudoforest-dp  pi=pseudoforest-dp-iterate\n\n";
    print_header();

    Totals tot;

    auto run_case = [&](const std::string& tag, const Shape& sh, const std::string& text) {
        Result g = run_algo(sh, text, max_add, CompressAlgo::GreedySize, phase2_iters);
        Result gd = run_algo(sh, text, max_add, CompressAlgo::GreedyDegree, phase2_iters);
        Result d = run_algo(sh, text, max_add, CompressAlgo::DepOrder, phase2_iters);
        Result p = run_algo(sh, text, max_add, CompressAlgo::PseudoforestDp, phase2_iters);
        Result pi = run_algo(sh, text, max_add, CompressAlgo::PseudoforestDpIterate, phase2_iters);
        print_row(tag, g, gd, d, p, pi);
        std::cout.flush();
        tot.add(g, gd, d, p, pi);
    };

    // --- Hand examples -----------------------------------------------------
    run_case("note #.# (n=40)",
             Shape::parse("#.#"),
             "ACGTCTTAAACCCTCGTCTTAAACCCAACGTCTTAAACCC");
    run_case("GCCTTTAAAG x3 #.#",
             Shape::parse("#.#"),
             "GCCTTTAAAGGCCTTTAAAGGCCTTTAAAG");
    run_case("GGGCGGCGGC ##",
             Shape::parse("##"),
             "GGGCGGCGGC");
    run_case("ac12 #.#",
             Shape::parse("#.#"),
             "ACACACACACAC");
    run_case("acgt20 #.#",
             Shape::parse("#.#"),
             "ACGTACGTACGTACGTACGT");

    // --- Cyclic preference structures --------------------------------------
    // Each name points at one preferred source, so the preference graph has
    // out-degree 1 and its cycles are what PseudoforestDp's exact unicyclic
    // solve (and, upstream, DepOrder's dependency order) has to handle
    // directly rather than approximate. These are periodic by construction:
    // a period-p repeat makes the p gapped k-mers derive from each other
    // cyclically, which stresses that handling at scale.
    auto rep = [](const std::string& motif, size_t n) {
        std::string s;
        s.reserve(n + motif.size());
        while (s.size() < n) s += motif;
        s.resize(n);
        return s;
    };
    const Shape sh_gap = Shape::parse("#.#");
    const Shape sh_adj = Shape::parse("##");
    run_case("cyc ac48 #.#", sh_gap, rep("AC", 48));
    run_case("cyc taaa48 #.#", sh_gap, rep("TAAA", 48));
    run_case("cyc gaca120 #.#", sh_gap, rep("GACA", 120));
    run_case("cyc taaa400 #.#", sh_gap, rep("TAAA", 400));
    run_case("cyc tgttct120 ##", sh_adj, rep("TGTTCT", 120));
    run_case("cyc acgt48 #.#", sh_gap, rep("ACGT", 48));
    run_case("cyc gtac48 #.#", sh_gap, rep("GTAC", 48));
    run_case("cyc gtac400 #.#", sh_gap, rep("GTAC", 400));

    // --- Short repetitive DNA ----------------------------------------------
    {
        std::mt19937_64 rng(seed);
        std::string motif = random_dna(200, rng);
        std::string text = concat_copies(motif, 20);
        for (const char* pat : {"#####", "##.##", "####.####", "#.#.#.#"})
            run_case(std::string("rep20x200 ") + pat, Shape::parse(pat), text);
    }

    std::cout << "\n";

    // --- Weight-30 suite ---------------------------------------------------
    std::mt19937_64 rng(seed);
    std::string motif = random_dna(motif_len, rng);
    auto shapes = make_shapes();
    std::cout << "Weight-30 suite: motif_len=" << motif_len
              << " reps=" << min_rep << ".." << max_rep << " step " << step
              << " shapes=" << shapes.size() << "\n";
    print_header();

    for (int reps = min_rep; reps <= max_rep; reps += step) {
        std::string text = concat_copies(motif, reps);
        for (auto& spec : shapes) {
            run_case("x" + std::to_string(reps) + " " + spec.label,
                     Shape::parse(spec.pattern), text);
        }
    }

    std::cout << "\n=== Summary ===\n"
              << "configs: " << tot.n << "\n"
              << "unique wins:  greedy-size=" << tot.win_g
              << "  greedy-degree=" << tot.win_gd
              << "  dep-order=" << tot.win_d
              << "  pseudoforest-dp=" << tot.win_p
              << "  pseudoforest-dp-iterate=" << tot.win_pi
              << "  ties/shared=" << tot.tie << "\n"
              << "total |C|:  greedy-size=" << tot.sum_g
              << "  greedy-degree=" << tot.sum_gd
              << "  dep-order=" << tot.sum_d
              << "  pseudoforest-dp=" << tot.sum_p
              << "  pseudoforest-dp-iterate=" << tot.sum_pi << "\n"
              << "vs greedy-size:" << std::fixed << std::setprecision(2);
    // Stored positions relative to greedy-size: "fewer" is better, "more" is worse.
    auto vs_greedy = [&](const char* name, long long a) {
        double p = (tot.sum_g == 0)
                       ? 0.0
                       : 100.0 * (1.0 - (double)a / (double)tot.sum_g);
        std::cout << "  " << name << " " << (p < 0 ? -p : p) << "% "
                  << (p < 0 ? "more" : "fewer");
    };
    vs_greedy("greedy-degree", tot.sum_gd);
    vs_greedy("dep-order", tot.sum_d);
    vs_greedy("pseudoforest-dp", tot.sum_p);
    vs_greedy("pseudoforest-dp-iterate", tot.sum_pi);
    std::cout << "\ncorrectness failures: " << tot.fails << "\n";
    return tot.fails ? 1 : 0;
}
