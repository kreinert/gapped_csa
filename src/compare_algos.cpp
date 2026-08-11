// compare_algos.cpp
//
// Side-by-side comparison of CompressAlgo::Greedy, DepOrder, GreedyDfs,
// TreeDp, TreeDp2 on the note's #.# example, the GCCTTTAAAG×3 demo, short
// repetitive DNA, and the weight-30 repetition suite.
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

static void print_header() {
    std::cout << std::left
              << std::setw(28) << "case"
              << std::setw(7)  << "g|C|"
              << std::setw(6)  << "g%"
              << std::setw(7)  << "d|C|"
              << std::setw(6)  << "d%"
              << std::setw(7)  << "f|C|"
              << std::setw(6)  << "f%"
              << std::setw(7)  << "t|C|"
              << std::setw(6)  << "t%"
              << std::setw(7)  << "t2|C|"
              << std::setw(6)  << "t2%"
              << std::setw(12) << "best"
              << std::setw(4)  << "g"
              << std::setw(4)  << "d"
              << std::setw(4)  << "f"
              << std::setw(4)  << "t"
              << std::setw(4)  << "t2"
              << "\n";
}

// g=greedy, d=dep-order, f=greedy-dfs, t=tree-dp, t2=tree-dp2 (no Phase II)
static void print_row(const std::string& tag,
                      const Result& g, const Result& d,
                      const Result& f, const Result& t, const Result& t2) {
    size_t bC = std::min({g.C, d.C, f.C, t.C, t2.C});
    std::string best;
    if (g.C == bC) best += (best.empty() ? "" : "/") + std::string("g");
    if (d.C == bC) best += (best.empty() ? "" : "/") + std::string("d");
    if (f.C == bC) best += (best.empty() ? "" : "/") + std::string("f");
    if (t.C == bC) best += (best.empty() ? "" : "/") + std::string("t");
    if (t2.C == bC) best += (best.empty() ? "" : "/") + std::string("t2");
    if (best == "g/d/f/t/t2") best = "tie";

    auto ok = [](const Result& r) { return (r.ok && r.rt) ? "ok" : "FAIL"; };
    std::cout << std::left << std::setw(28) << tag
              << std::setw(7) << g.C
              << std::setw(6) << std::fixed << std::setprecision(1) << g.keep_pct
              << std::setw(7) << d.C
              << std::setw(6) << d.keep_pct
              << std::setw(7) << f.C
              << std::setw(6) << f.keep_pct
              << std::setw(7) << t.C
              << std::setw(6) << t.keep_pct
              << std::setw(7) << t2.C
              << std::setw(6) << t2.keep_pct
              << std::setw(12) << best
              << std::setw(4) << ok(g)
              << std::setw(4) << ok(d)
              << std::setw(4) << ok(f)
              << std::setw(4) << ok(t)
              << std::setw(4) << ok(t2)
              << "\n";
}

struct Totals {
    long long sum_g = 0, sum_d = 0, sum_f = 0, sum_t = 0, sum_t2 = 0;
    int win_g = 0, win_d = 0, win_f = 0, win_t = 0, win_t2 = 0, tie = 0;
    int fails = 0, n = 0;

    void add(const Result& g, const Result& d, const Result& f,
             const Result& t, const Result& t2) {
        ++n;
        sum_g += (long long)g.C;
        sum_d += (long long)d.C;
        sum_f += (long long)f.C;
        sum_t += (long long)t.C;
        sum_t2 += (long long)t2.C;
        if (!g.ok || !g.rt || !d.ok || !d.rt || !f.ok || !f.rt ||
            !t.ok || !t.rt || !t2.ok || !t2.rt) ++fails;
        size_t b = std::min({g.C, d.C, f.C, t.C, t2.C});
        int winners = (g.C == b) + (d.C == b) + (f.C == b) + (t.C == b) + (t2.C == b);
        if (winners > 1) ++tie;
        else if (g.C == b) ++win_g;
        else if (d.C == b) ++win_d;
        else if (f.C == b) ++win_f;
        else if (t.C == b) ++win_t;
        else ++win_t2;
    }
};

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

    std::cout << "Comparing greedy | dep-order | greedy-dfs | tree-dp | tree-dp2"
              << "  (max_add=" << max_add << ")\n"
              << "columns: g=greedy  d=dep-order  f=greedy-dfs  t=tree-dp"
              << "  t2=tree-dp2 (no Phase II)\n\n";
    print_header();

    Totals tot;

    auto run_case = [&](const std::string& tag, const Shape& sh, const std::string& text) {
        Result g = run_algo(sh, text, max_add, CompressAlgo::Greedy);
        Result d = run_algo(sh, text, max_add, CompressAlgo::DepOrder);
        Result f = run_algo(sh, text, max_add, CompressAlgo::GreedyDfs);
        Result t = run_algo(sh, text, max_add, CompressAlgo::TreeDp);
        Result t2 = run_algo(sh, text, max_add, CompressAlgo::TreeDp2);
        print_row(tag, g, d, f, t, t2);
        std::cout.flush();
        tot.add(g, d, f, t, t2);
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

    auto pct = [](long long a, long long b) -> double {
        if (b == 0) return 0.0;
        return 100.0 * (1.0 - (double)a / (double)b);
    };
    std::cout << "\n=== Summary ===\n"
              << "configs: " << tot.n << "\n"
              << "unique wins:  greedy=" << tot.win_g
              << "  dep-order=" << tot.win_d
              << "  greedy-dfs=" << tot.win_f
              << "  tree-dp=" << tot.win_t
              << "  tree-dp2=" << tot.win_t2
              << "  ties/shared=" << tot.tie << "\n"
              << "total |C|:  greedy=" << tot.sum_g
              << "  dep-order=" << tot.sum_d
              << "  greedy-dfs=" << tot.sum_f
              << "  tree-dp=" << tot.sum_t
              << "  tree-dp2=" << tot.sum_t2 << "\n"
              << "vs greedy:  dep-order " << std::fixed << std::setprecision(2)
              << pct(tot.sum_d, tot.sum_g) << "% fewer"
              << "  greedy-dfs " << pct(tot.sum_f, tot.sum_g) << "% fewer"
              << "  tree-dp " << pct(tot.sum_t, tot.sum_g) << "% fewer"
              << "  tree-dp2 " << pct(tot.sum_t2, tot.sum_g) << "% fewer\n"
              << "correctness failures: " << tot.fails << "\n";
    return tot.fails ? 1 : 0;
}
