// compare_algos.cpp
//
// Side-by-side comparison of CompressAlgo::Greedy, DepOrder, TreeDp, TreeDp2,
// TreeDp3, TreeDp4, PseudoforestDp on the note's #.# example, the
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
              << std::setw(7)  << "d|C|"
              << std::setw(6)  << "d%"
              << std::setw(7)  << "t|C|"
              << std::setw(6)  << "t%"
              << std::setw(7)  << "t2|C|"
              << std::setw(6)  << "t2%"
              << std::setw(7)  << "t3|C|"
              << std::setw(6)  << "t3%"
              << std::setw(7)  << "t4|C|"
              << std::setw(6)  << "t4%"
              << std::setw(7)  << "p|C|"
              << std::setw(6)  << "p%"
              << std::setw(16) << "best"
              << std::setw(4)  << "g"
              << std::setw(4)  << "d"
              << std::setw(4)  << "t"
              << std::setw(4)  << "t2"
              << std::setw(4)  << "t3"
              << std::setw(4)  << "t4"
              << std::setw(4)  << "p"
              << "\n";
}

// g=greedy, d=dep-order, t=tree-dp, t2=tree-dp2, t3=tree-dp3, t4=tree-dp4,
// p=pseudoforest-dp
static void print_row(const std::string& tag,
                      const Result& g, const Result& d, const Result& t,
                      const Result& t2, const Result& t3, const Result& t4,
                      const Result& p) {
    size_t bC = std::min({g.C, d.C, t.C, t2.C, t3.C, t4.C, p.C});
    std::string best;
    if (g.C == bC) best += (best.empty() ? "" : "/") + std::string("g");
    if (d.C == bC) best += (best.empty() ? "" : "/") + std::string("d");
    if (t.C == bC) best += (best.empty() ? "" : "/") + std::string("t");
    if (t2.C == bC) best += (best.empty() ? "" : "/") + std::string("t2");
    if (t3.C == bC) best += (best.empty() ? "" : "/") + std::string("t3");
    if (t4.C == bC) best += (best.empty() ? "" : "/") + std::string("t4");
    if (p.C == bC) best += (best.empty() ? "" : "/") + std::string("p");
    if (best == "g/d/t/t2/t3/t4/p") best = "tie";

    auto ok = [](const Result& r) { return (r.ok && r.rt) ? "ok" : "FAIL"; };
    std::cout << std::left << std::setw(28) << tag
              << std::setw(7) << g.C
              << std::setw(6) << std::fixed << std::setprecision(1) << g.keep_pct
              << std::setw(7) << d.C
              << std::setw(6) << d.keep_pct
              << std::setw(7) << t.C
              << std::setw(6) << t.keep_pct
              << std::setw(7) << t2.C
              << std::setw(6) << t2.keep_pct
              << std::setw(7) << t3.C
              << std::setw(6) << t3.keep_pct
              << std::setw(7) << t4.C
              << std::setw(6) << t4.keep_pct
              << std::setw(7) << p.C
              << std::setw(6) << p.keep_pct
              << std::setw(16) << best
              << std::setw(4) << ok(g)
              << std::setw(4) << ok(d)
              << std::setw(4) << ok(t)
              << std::setw(4) << ok(t2)
              << std::setw(4) << ok(t3)
              << std::setw(4) << ok(t4)
              << std::setw(4) << ok(p)
              << "\n";
}

struct Totals {
    long long sum_g = 0, sum_d = 0, sum_t = 0, sum_t2 = 0, sum_t3 = 0, sum_t4 = 0, sum_p = 0;
    int win_g = 0, win_d = 0, win_t = 0, win_t2 = 0, win_t3 = 0, win_t4 = 0, win_p = 0, tie = 0;
    int fails = 0, n = 0;

    void add(const Result& g, const Result& d, const Result& t,
             const Result& t2, const Result& t3, const Result& t4,
             const Result& p) {
        ++n;
        sum_g += (long long)g.C;
        sum_d += (long long)d.C;
        sum_t += (long long)t.C;
        sum_t2 += (long long)t2.C;
        sum_t3 += (long long)t3.C;
        sum_t4 += (long long)t4.C;
        sum_p += (long long)p.C;
        if (!g.ok || !g.rt || !d.ok || !d.rt ||
            !t.ok || !t.rt || !t2.ok || !t2.rt || !t3.ok || !t3.rt ||
            !t4.ok || !t4.rt || !p.ok || !p.rt) ++fails;
        size_t b = std::min({g.C, d.C, t.C, t2.C, t3.C, t4.C, p.C});
        int winners = (g.C == b) + (d.C == b) + (t.C == b) + (t2.C == b)
                    + (t3.C == b) + (t4.C == b) + (p.C == b);
        if (winners > 1) ++tie;
        else if (g.C == b) ++win_g;
        else if (d.C == b) ++win_d;
        else if (t.C == b) ++win_t;
        else if (t2.C == b) ++win_t2;
        else if (t3.C == b) ++win_t3;
        else if (t4.C == b) ++win_t4;
        else ++win_p;
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

    std::cout << "Comparing greedy | dep-order | tree-dp | tree-dp2 | tree-dp3 | tree-dp4 | pseudoforest-dp"
              << "  (max_add=" << max_add << ")\n"
              << "columns: g=greedy  d=dep-order  t=tree-dp  t2=tree-dp2"
              << "  t3=tree-dp3  t4=tree-dp4  p=pseudoforest-dp\n\n";
    print_header();

    Totals tot;

    auto run_case = [&](const std::string& tag, const Shape& sh, const std::string& text) {
        Result g = run_algo(sh, text, max_add, CompressAlgo::Greedy, phase2_iters);
        Result d = run_algo(sh, text, max_add, CompressAlgo::DepOrder, phase2_iters);
        Result t = run_algo(sh, text, max_add, CompressAlgo::TreeDp, phase2_iters);
        Result t2 = run_algo(sh, text, max_add, CompressAlgo::TreeDp2, phase2_iters);
        Result t3 = run_algo(sh, text, max_add, CompressAlgo::TreeDp3, phase2_iters);
        Result t4 = run_algo(sh, text, max_add, CompressAlgo::TreeDp4, phase2_iters);
        Result p = run_algo(sh, text, max_add, CompressAlgo::PseudoforestDp, phase2_iters);
        print_row(tag, g, d, t, t2, t3, t4, p);
        std::cout.flush();
        tot.add(g, d, t, t2, t3, t4, p);
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
    // Known tree-dp counterexamples (ILP optimum |C|=9 resp. 11).
    run_case("ac12 #.#",
             Shape::parse("#.#"),
             "ACACACACACAC");
    run_case("acgt20 #.#",
             Shape::parse("#.#"),
             "ACGTACGTACGTACGTACGT");

    // --- Cyclic preference structures --------------------------------------
    // Each name points at one preferred source, so the preference graph has
    // out-degree 1 and its cycles are what the forest build must break. It
    // breaks them at the member processed last -- the largest packed k-mer name
    // -- which is an artifact of iteration order, not of value; tree-dp4 drops
    // the cycle's cheapest edge instead. These are the families where that
    // choice is load-bearing. They are periodic by construction: a period-p
    // repeat makes the p gapped k-mers derive from each other cyclically.
    //
    // The rows below deliberately cover both signs. Repair maximizes the forest
    // DP's own objective, but the DP only feeds Phase II, and re-rooting the
    // chain can leave the retarget loop with less to work with -- so the same
    // move that wins 97 on tgac-style motifs loses 97 on gtac ones.
    auto rep = [](const std::string& motif, size_t n) {
        std::string s;
        s.reserve(n + motif.size());
        while (s.size() < n) s += motif;
        s.resize(n);
        return s;
    };
    const Shape sh_gap = Shape::parse("#.#");
    const Shape sh_adj = Shape::parse("##");
    // Period 2: the ac12 counterexample at scale. One 2-cycle, coverages 3 vs
    // 4; the gain stays 1 however long the text gets.
    run_case("cyc ac48 #.#", sh_gap, rep("AC", 48));
    // Period 4 with motif[1]==motif[3]: the odd positions collapse to a single
    // name, leaving one cycle whose repair unlocks the whole chain. Nominal
    // edge gain is 1, but end to end it halves |C| and hits the ILP optimum.
    run_case("cyc taaa48 #.#", sh_gap, rep("TAAA", 48));
    run_case("cyc gaca120 #.#", sh_gap, rep("GACA", 120));
    run_case("cyc taaa400 #.#", sh_gap, rep("TAAA", 400));
    // Two interacting cycles, both repaired.
    run_case("cyc tgttct120 ##", sh_adj, rep("TGTTCT", 120));
    // Boundary: a period-4 forest *path*, no cycle to repair (dep-order and
    // tree-dp3 reach the optimum here; cycle repair cannot).
    run_case("cyc acgt48 #.#", sh_gap, rep("ACGT", 48));
    // Counter-family: a rotation of the same cyclic sequence, where plain
    // tree-dp is already optimal and repairing the cycle breaks it.
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
              << "unique wins:  greedy=" << tot.win_g
              << "  dep-order=" << tot.win_d
              << "  tree-dp=" << tot.win_t
              << "  tree-dp2=" << tot.win_t2
              << "  tree-dp3=" << tot.win_t3
              << "  tree-dp4=" << tot.win_t4
              << "  pseudoforest-dp=" << tot.win_p
              << "  ties/shared=" << tot.tie << "\n"
              << "total |C|:  greedy=" << tot.sum_g
              << "  dep-order=" << tot.sum_d
              << "  tree-dp=" << tot.sum_t
              << "  tree-dp2=" << tot.sum_t2
              << "  tree-dp3=" << tot.sum_t3
              << "  tree-dp4=" << tot.sum_t4
              << "  pseudoforest-dp=" << tot.sum_p << "\n"
              << "vs greedy:" << std::fixed << std::setprecision(2);
    // Stored positions relative to greedy: "fewer" is better, "more" is worse.
    auto vs_greedy = [&](const char* name, long long a) {
        double p = (tot.sum_g == 0)
                       ? 0.0
                       : 100.0 * (1.0 - (double)a / (double)tot.sum_g);
        std::cout << "  " << name << " " << (p < 0 ? -p : p) << "% "
                  << (p < 0 ? "more" : "fewer");
    };
    vs_greedy("dep-order", tot.sum_d);
    vs_greedy("tree-dp", tot.sum_t);
    vs_greedy("tree-dp2", tot.sum_t2);
    vs_greedy("tree-dp3", tot.sum_t3);
    vs_greedy("tree-dp4", tot.sum_t4);
    vs_greedy("pseudoforest-dp", tot.sum_p);
    std::cout << "\ncorrectness failures: " << tot.fails << "\n";
    return tot.fails ? 1 : 0;
}
