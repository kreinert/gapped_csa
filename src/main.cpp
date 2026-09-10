// gcsa — build the compressed gapped-shape position index, validate it, print
// the hash table, and answer queries.
//
// Usage:
//   ./gcsa                                  # run the note's #.# example + self-test
//   ./gcsa -g <fasta> -s <shape> [-q <query>] [--table] [--max-add N]
//                                           [--algo A] [--phase2-iters N]
//   ./gcsa -g <fasta> -s <shape> -r <reads.fasta>   # locate reads (like kmer/locate)
//   ./gcsa --help                           # full option list

#include "compress.hpp"
#include "serialize.hpp"
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <sstream>
#include <unordered_map>

using namespace gcsa;

static std::string load_dna(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string out, line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') continue;
        for (char c : line)
            if (std::isalpha((unsigned char)c)) out.push_back(base_char(base_value(c)));
    }
    return out;
}

struct FastaRec { std::string header, seq; };
static std::vector<FastaRec> parse_fasta(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::vector<FastaRec> v; std::string line; FastaRec cur;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') { if (!cur.header.empty()||!cur.seq.empty()) v.push_back(std::move(cur)); cur = {line.substr(1), ""}; }
        else { for (char c: line) if (std::isalpha((unsigned char)c)) cur.seq.push_back(base_char(base_value(c))); }
    }
    if (!cur.header.empty()||!cur.seq.empty()) v.push_back(std::move(cur));
    return v;
}

// Brute-force ground truth: name -> set of text positions (0..n).
static std::map<uint64_t, std::vector<int64_t>> brute_positions(const Shape& sh, const std::string& text) {
    std::map<uint64_t, std::vector<int64_t>> m;
    for (long p = 0; p <= (long)text.size(); ++p)
        m[name_at(sh, text, (size_t)p)].push_back(p);
    return m;
}

static bool self_test(const CompressedIndex& idx, const Shape& sh, const std::string& text, bool verbose) {
    auto truth = brute_positions(sh, text);
    bool ok = true;
    for (auto& kv : truth) {
        auto got = idx.positions_of(kv.first);
        std::sort(got.begin(), got.end());
        auto exp = kv.second;
        std::sort(exp.begin(), exp.end());
        if (got != exp) {
            ok = false;
            if (verbose) {
                std::cout << "MISMATCH name=" << kv.first << " (" << name_to_string(sh, kv.first) << ")\n";
                std::cout << "  expected:"; for (auto x: exp) std::cout << ' ' << x; std::cout << "\n";
                std::cout << "  got     :"; for (auto x: got) std::cout << ' ' << x; std::cout << "\n";
            }
        }
    }
    return ok;
}

// Print the full gapped suffix array (with shapes/codes and LCP) and, side by
// side, the compressed SA: kept rows show their index in C, removed rows show x.
static void print_sa_tables(const CompressedIndex& idx, int max_suffix_syms = 14) {
    const GappedSA& G = idx.gsa();
    const Shape& sh = G.shape;

    std::cout << "\n=== Full " << (idx.regular_mode() ? "plain" : "gapped")
              << " suffix array (shape " << sh.pattern << ") ===\n";
    std::cout << std::right
              << std::setw(4) << "Rk" << " | " << std::setw(6) << "lexSA" << " | "
              << std::setw(5) << "orig" << " | " << std::setw(6) << "code" << " | "
              << std::setw(6) << "kmer" << " | " << std::setw(4) << "lcp" << " | suffix(codes)\n";
    for (size_t r = 0; r < G.m(); ++r) {
        int32_t lexpos = G.sa[r];
        uint64_t nm = idx.name_of_rank((int32_t)r);
        std::cout << std::setw(4) << r << " | " << std::setw(6) << lexpos << " | "
                  << std::setw(5) << G.orig_pos((int32_t)r) << " | "
                  << std::setw(6) << nm << " | "
                  << std::setw(6) << name_to_string(sh, nm) << " | "
                  << std::setw(4) << G.lcp[r] << " | ";
        for (int k = 0; k < max_suffix_syms && lexpos + k < (int)G.m(); ++k)
            std::cout << G.lex[lexpos + k] << ' ';
        std::cout << "\n";
    }

    std::cout << "\n=== Compressed suffix array (x = removed / recovered via a pointer) ===\n";
    std::cout << std::right
              << std::setw(4) << "Rk" << " | " << std::setw(7) << "C-idx" << " | "
              << std::setw(5) << "orig" << " | " << std::setw(6) << "kmer" << " | "
              << std::setw(4) << "lcp" << " | stored?\n";
    for (size_t r = 0; r < G.m(); ++r) {
        bool kept = idx.rank_kept((int32_t)r);
        std::cout << std::setw(4) << r << " | ";
        if (kept) std::cout << std::setw(7) << idx.rank_c_index((int32_t)r);
        else      std::cout << std::setw(7) << "x";
        std::cout << " | " << std::setw(5) << G.orig_pos((int32_t)r) << " | "
                  << std::setw(6) << name_to_string(sh, idx.name_of_rank((int32_t)r)) << " | "
                  << std::setw(4) << G.lcp[r] << " | "
                  << (kept ? "C[" + std::to_string(idx.rank_c_index((int32_t)r)) + "] = "
                             + std::to_string(G.orig_pos((int32_t)r))
                           : std::string("(pointer)")) << "\n";
    }
    std::cout << "\ncompressed position array C (kept orig positions, rank order):\n  [";
    const auto& C = idx.compressed_positions();
    for (size_t i = 0; i < C.size(); ++i) std::cout << C[i] << (i + 1 < C.size() ? ", " : "");
    std::cout << "]\n";
}

static void print_hash_table(const CompressedIndex& idx) {
    const Shape& sh = idx.gsa().shape;
    // order by k-mer name
    std::map<uint64_t, HashEntry> ordered(idx.table().begin(), idx.table().end());
    std::cout << "\n" << std::left
              << std::setw(6) << "kmer" << std::setw(8) << "name"
              << std::setw(22) << "H_offset(pos,add,num)"
              << "H_rest(pos,num)\n";
    for (auto& kv : ordered) {
        const HashEntry& e = kv.second;
        std::ostringstream o, r;
        if (e.has_offset) o << "(" << e.off_pos << "," << e.off_add << "," << e.off_num << ")"; else o << "nil";
        if (e.has_rest)   r << "(" << e.rest_pos << "," << e.rest_num << ")";              else r << "nil";
        std::cout << std::setw(6) << name_to_string(sh, kv.first)
                  << std::setw(8) << kv.first
                  << std::setw(22) << o.str()
                  << r.str() << "\n";
    }
}

static void print_usage(const char* argv0) {
    std::cout
        << "Usage:\n"
        << "  " << argv0 << "                                 # note's #.# example + self-test\n"
        << "  " << argv0 << " -g <fasta> -s <shape> [options]\n"
        << "  " << argv0 << " -g <fasta> -s <shape> -r <reads.fasta>   # locate reads\n"
        << "\nOptions:\n"
        << "  -g <fasta>          genome FASTA (default: the note's example text)\n"
        << "  -s <shape>          gapped shape, e.g. \"#.#\" (default: #.#); or a bare\n"
        << "                      positive integer k (e.g. \"6\") for a REGULAR\n"
        << "                      (ungapped) suffix array compressed to support k-mer\n"
        << "                      lookup: no LexText is built, compression intervals\n"
        << "                      are LCP>=k (characters) instead of LCP>=1 (symbols).\n"
        << "                      Same positions_of()/locate() results as the\n"
        << "                      equivalent contiguous shape (k '#'s), different\n"
        << "                      (cheaper) construction path.\n"
        << "  --shape=<value>     long-form alias for -s\n"
        << "  -q <query>          query string to locate\n"
        << "  -r <reads.fasta>    locate every read in a FASTA\n"
        << "  --algo <name>       greedy-size|greedy-degree|dep-order|pseudoforest-dp|\n"
        << "                      pseudoforest-dp-iterate (default: greedy-size)\n"
        << "  --max-add N         max differential offset (default: 8)\n"
        << "  --phase2-iters N    Phase II dirty-generation budget (default: "
        << kPhase2DefaultIters << ").\n"
        << "                      Phase II stops earlier at a fixed point. Applies to\n"
        << "                      every --algo. Precedence: this flag >\n"
        << "                      GCSA_PHASE2_MAX_ITERS > default.\n"
        << "                      Phase II's two passes switch off independently, for\n"
        << "                      any --algo: GCSA_DISABLE_RETARGET=1 skips the greedy\n"
        << "                      unpin/retarget sweep, GCSA_DISABLE_LNS=1 skips the\n"
        << "                      exact cluster LNS, both together leave Phase I's\n"
        << "                      output as-is (e.g. to measure the forest DP alone).\n"
        << "                      GCSA_LNS_ONLY / GCSA_DISABLE_PHASE2 are deprecated\n"
        << "                      spellings of GCSA_DISABLE_RETARGET.\n"
        << "  --min-coverage N|auto  Minimum rows a candidate must reroute to be\n"
        << "                      accepted (default: " << kMinCoverage << "). \"auto\" derives the\n"
        << "                      break-even coverage from the input length and\n"
        << "                      --max-add instead (see min_coverage_breakeven() in\n"
        << "                      compress.hpp). Precedence: this flag >\n"
        << "                      GCSA_MIN_COVERAGE > default.\n"
        << "  --table             print the hash table\n"
        << "  --sa                print the full gapped suffix array\n"
        << "  --bound             print lower bounds on |C| (structural floor + a\n"
        << "                      per-word best-candidate-ignoring-conflicts bound built\n"
        << "                      from THIS run's own enumerate_candidates_ search --\n"
        << "                      i.e. what this exact algorithm/settings combo could\n"
        << "                      reach if Phase II resolved every conflict perfectly;\n"
        << "                      see CompressedIndex::compute_bound_report in\n"
        << "                      compress.hpp). Cheap (no O(L^2) sub-run enumeration);\n"
        << "                      still opt-in since it's not part of the index build.\n"
        << "  -h, --help          this message\n";
}

int main(int argc, char** argv) {
    std::string genome_path, shape_str = "#.#", query, reads_path;
    bool want_table = false, want_sa = false, want_bound = false;
    int max_add = 8;
    int phase2_iters = 0;  // 0 = unset: GCSA_PHASE2_MAX_ITERS, else the default
    CompressAlgo algo = CompressAlgo::GreedySize;
    std::string min_cov_arg;  // "" = unset -> GCSA_MIN_COVERAGE -> built-in default

    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto need = [&](const char* n){ if (++i>=argc){ std::cerr<<"missing value for "<<n<<"\n"; std::exit(1);} return std::string(argv[i]); };
        if      (k == "-g") genome_path = need("-g");
        else if (k == "-s") shape_str   = need("-s");
        else if (k.rfind("--shape=", 0) == 0) shape_str = k.substr(8);
        else if (k == "-q") query       = need("-q");
        else if (k == "-r") reads_path  = need("-r");
        else if (k == "--table") want_table = true;
        else if (k == "--sa") want_sa = true;
        else if (k == "--bound") want_bound = true;
        else if (k == "--max-add") max_add = std::stoi(need("--max-add"));
        else if (k == "--phase2-iters") {
            phase2_iters = std::stoi(need("--phase2-iters"));
            if (phase2_iters < 1) {
                std::cerr << "--phase2-iters must be >= 1\n";
                return 1;
            }
        }
        else if (k == "--min-coverage") min_cov_arg = need("--min-coverage");
        else if (k == "--help" || k == "-h") { print_usage(argv[0]); return 0; }
        else if (k == "--algo") {
            std::string a = need("--algo");
            if (a == "greedy-size" || a == "greedy" || a == "gs") algo = CompressAlgo::GreedySize;
            else if (a == "greedy-degree" || a == "gd") algo = CompressAlgo::GreedyDegree;
            else if (a == "dep-order" || a == "dep") algo = CompressAlgo::DepOrder;
            else if (a == "pseudoforest-dp" || a == "pf-dp" || a == "pfdp") algo = CompressAlgo::PseudoforestDp;
            else if (a == "pseudoforest-dp-iterate" || a == "pf-dp-iterate" || a == "pfdp-iterate" || a == "pfdpi") algo = CompressAlgo::PseudoforestDpIterate;
            else {
                std::cerr << "--algo must be greedy-size|greedy-degree|dep-order|pseudoforest-dp|pseudoforest-dp-iterate\n";
                return 1;
            }
        }
        else { std::cerr << "unknown arg: " << k << "\n"; print_usage(argv[0]); return 1; }
    }

    bool is_example = genome_path.empty();
    std::string text;
    if (is_example) {
        text = "aCGTCTTAAACCCtCGTCTTAAACCCaaCGTCTTAAACCC";
        std::string clean; for (char c: text) if (std::isalpha((unsigned char)c)) clean.push_back(base_char(base_value(c)));
        text = clean;
        want_table = true;
        want_sa = true;
        if (query.empty()) query = "CTTAAC";
        std::cout << "=== Note's #.# example ===\ntext = " << text << "\n";
        // Verify the marker-byte encoding against the note's worked example:
        // (2000000,16,300),(3000000000,24) must start with 10100100.
        {
            HashEntry e; e.has_offset = true; e.off_pos = 2000000; e.off_add = 16; e.off_num = 300;
            e.has_rest = true; e.rest_pos = 3000000000ULL; e.rest_num = 24;
            int off_add_b = bytes_needed(e.off_add), off_cnt_b = bytes_needed(e.off_num), rest_cnt_b = bytes_needed(e.rest_num);
            uint8_t marker = 0x80 | ((uint8_t)(off_cnt_b-1)<<5) | ((uint8_t)(off_add_b-1)<<3) | 0x04 | (uint8_t)(rest_cnt_b-1);
            std::cout << "marker-byte check for (2000000,16,300),(3000000000,24) = ";
            for (int b = 7; b >= 0; --b) std::cout << ((marker>>b)&1);
            std::cout << (marker == 0b10100100 ? "  (== 10100100, OK)\n" : "  (MISMATCH!)\n");
        }
    } else {
        text = load_dna(genome_path);
    }

    // A bare positive integer for -s/--shape means "regular (ungapped) suffix
    // array, compressed to support k-mer lookup" instead of a #/. pattern: no
    // LexText is built, and compression intervals are LCP>=k instead of
    // LCP>=1 (see build(int, ...) in compress.hpp / build_plain_sa in
    // gapped_sa.hpp). The Shape below is only used for naming/self-test here
    // -- it's the same contiguous shape build(int,...) uses internally, so
    // it's consistent with the index by construction.
    bool regular = !shape_str.empty() &&
        shape_str.find_first_not_of("0123456789") == std::string::npos;
    int regular_k = 0;
    if (regular) {
        regular_k = std::stoi(shape_str);
        if (regular_k < 1) {
            std::cerr << "-s/--shape: k must be >= 1\n";
            return 1;
        }
    }
    Shape sh = regular ? Shape::parse(std::string((size_t)regular_k, '#'))
                        : Shape::parse(shape_str);
    if (regular)
        std::cout << "mode=regular (no LexText) k=" << regular_k
                  << " n=" << text.size() << "  algo=" << algo_name(algo) << "\n";
    else
        std::cout << "shape=" << shape_str << " span=" << sh.span << " weight=" << sh.weight
                  << " n=" << text.size() << "  algo=" << algo_name(algo) << "\n";

    // Resolve the min-coverage floor before building: it gates candidate
    // enumeration inside idx.build(), so it must be set first. Precedence:
    // --min-coverage > GCSA_MIN_COVERAGE > built-in default (kMinCoverage's
    // compiled-in value). Either source may be "auto" to derive the
    // break-even coverage from the input length (pre-DP pb estimate) and
    // --max-add instead of the fixed default.
    {
        std::string src = min_cov_arg;
        const char* origin = "--min-coverage";
        if (src.empty()) {
            if (const char* env = std::getenv("GCSA_MIN_COVERAGE")) { src = env; origin = "GCSA_MIN_COVERAGE"; }
        }
        if (!src.empty()) {
            if (src == "auto") {
                const int pb_est = estimate_pb(text.size());
                const int mc = min_coverage_breakeven(pb_est, max_add);
                std::cout << "min-coverage: auto (" << origin << ") pb_estimate=" << pb_est
                          << " max_add=" << max_add << " -> min_coverage=" << mc << "\n";
                kMinCoverage = mc;
            } else {
                int mc = std::stoi(src);
                if (mc < 1) {
                    std::cerr << "--min-coverage must be >= 1 (or \"auto\")\n";
                    return 1;
                }
                std::cout << "min-coverage: " << mc << " (" << origin << ")\n";
                kMinCoverage = mc;
            }
        }
    }

    CompressedIndex idx;
    if (regular) idx.build(regular_k, text, max_add, algo, phase2_iters);
    else         idx.build(sh, text, max_add, algo, phase2_iters);

    std::cout << "distinct k-mers      : " << idx.num_kmers() << "\n";
    std::cout << "SA entries (m)       : " << idx.total_positions() << "\n";
    std::cout << "stored positions (|C|): " << idx.stored_positions()
              << "  (" << std::fixed << std::setprecision(1)
              << 100.0 * idx.stored_positions() / std::max<size_t>(1, idx.total_positions())
              << "% of full SA)\n";

    // --bound: lower bounds on |C|, computed off the already-built gapped SA
    // (no rebuild). See CompressedIndex::compute_bound_report in compress.hpp
    // for the derivation of both numbers. The second one is deliberately
    // algorithm-consistent: it's built from this exact run's enumerate_candidates_
    // search (same GCSA_CROSS_LCP / GCSA_INTRA_LINKS / --min-coverage this
    // build used), not the ILP baseline's full link universe -- see the
    // comment on compute_bound_report for why, and for how to get the other
    // (algorithm-agnostic, more expensive) number if you want it instead.
    if (want_bound) {
        using Clock = std::chrono::steady_clock;
        const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
        auto t0 = Clock::now();
        CompressedIndex::BoundReport B = idx.compute_bound_report();
        auto t1 = Clock::now();
        const size_t m = std::max<size_t>(1, idx.total_positions());
        std::cout << "bound: floor(structural, |I_c|<" << kMinCoverage << ")=" << B.floor_positions
                  << "  (" << std::fixed << std::setprecision(3)
                  << 100.0 * B.floor_positions / m << "% of m)\n";
        std::cout << "       reachable(best-candidate-per-word via enumerate_candidates_,"
                     " ignoring conflicts, cross_lcp=" << (B.cross_lcp ? 1 : 0)
                  << " intra_links=" << (B.intra_links ? 1 : 0) << ")=" << B.lower_bound_C
                  << "  (" << std::fixed << std::setprecision(3)
                  << 100.0 * B.lower_bound_C / m << "% of m)"
                  << "  [floor <= this <= |C| from any run with these same settings]\n";
        std::cout << "       eligible words (|I_c|>=" << kMinCoverage << "): " << B.eligible_words
                  << "  positions=" << B.eligible_positions
                  << "  words_with_no_candidate=" << (B.eligible_words - B.words_with_candidate)
                  << "  sum_best_coverage=" << B.sum_best_coverage << "\n";
        if (timing) {
            std::cout << "[timing] bound="
                      << std::chrono::duration<double, std::milli>(t1 - t0).count()
                      << "ms\n";
        }
    }

    // Brute-force self-test + round-trip walk every text position. On large N
    // that dominates wall time / RAM; GCSA_SKIP_SELFTEST=1 keeps compress+timing.
    const bool skip_selftest = [] {
        const char* e = std::getenv("GCSA_SKIP_SELFTEST");
        return e && std::atoi(e) != 0;
    }();

    if (skip_selftest) {
        std::cout << "self-test (vs brute force): SKIPPED (GCSA_SKIP_SELFTEST)\n";
        std::cout << "serialized round-trip     : SKIPPED (GCSA_SKIP_SELFTEST)\n";
        SerializedIndex S = serialize_index(idx);
        double full_sa = (double)idx.total_positions() * S.pb;
        std::cout << "pb(bytes/pos)=" << S.pb
                  << "  |C|*pb=" << S.bytes_C()
                  << "  H=" << S.bytes_H()
                  << "  S_H=" << S.bytes_SH()
                  << "  names=" << S.bytes_names()
                  << "  total=" << S.total_bytes()
                  << "  (full SA would be " << (size_t)full_sa << " bytes; "
                  << std::fixed << std::setprecision(1)
                  << 100.0 * S.total_bytes() / std::max(1.0, full_sa) << "%)\n";
    } else {
        bool ok = self_test(idx, sh, text, is_example);
        std::cout << "self-test (vs brute force): " << (ok ? "PASS" : "FAIL") << "\n";

        // Byte-level serialization + round-trip check + size accounting.
        SerializedIndex S = serialize_index(idx);
        bool rt_ok = true;
        {
            auto truth = brute_positions(sh, text);
            for (auto& kv : truth) {
                auto got = S.positions_of(kv.first);
                std::sort(got.begin(), got.end());
                auto exp = kv.second; std::sort(exp.begin(), exp.end());
                if (got != exp) { rt_ok = false; break; }
            }
        }
        std::cout << "serialized round-trip     : " << (rt_ok ? "PASS" : "FAIL") << "\n";
        double full_sa = (double)idx.total_positions() * S.pb;   // uncompressed SA bytes
        std::cout << "pb(bytes/pos)=" << S.pb
                  << "  |C|*pb=" << S.bytes_C()
                  << "  H=" << S.bytes_H()
                  << "  S_H=" << S.bytes_SH()
                  << "  names=" << S.bytes_names()
                  << "  total=" << S.total_bytes()
                  << "  (full SA would be " << (size_t)full_sa << " bytes; "
                  << std::fixed << std::setprecision(1)
                  << 100.0 * S.total_bytes() / std::max(1.0, full_sa) << "%)\n";
    }

    if (want_sa) print_sa_tables(idx);
    if (want_table) print_hash_table(idx);

    if (!query.empty()) {
        std::cout << "\nquery Q=" << query << "\n";
        for (int i = 0; i + sh.span <= (int)query.size(); ++i) {
            uint64_t nm = name_at(sh, query, i);
            auto pos = idx.positions_of(nm);
            std::sort(pos.begin(), pos.end());
            std::cout << "  window@" << i << " kmer=" << name_to_string(sh, nm)
                      << " (name " << nm << ") -> positions:";
            for (auto p : pos) std::cout << ' ' << p;
            std::cout << "\n";
        }
        auto all = idx.locate(query);
        std::cout << "  locate(all windows, sorted unique):";
        for (auto p : all) std::cout << ' ' << p;
        std::cout << "\n";
    }

    if (!reads_path.empty()) {
        auto reads = parse_fasta(reads_path);
        long long total_hits = 0;
        for (auto& r : reads) total_hits += (long long)idx.locate(r.seq).size();
        std::cout << "\nreads processed: " << reads.size()
                  << "  total hit positions: " << total_hits << "\n";
    }
    return 0;
}
