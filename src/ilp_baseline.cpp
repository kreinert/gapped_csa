// ilp_baseline.cpp — exact |C| via ILP over the full differential-link universe
// of compress.hpp (kMinCoverage=3, max_add default 8): see "Link universe"
// there. Because that universe is a superset of what every heuristic can build,
// the optimum reported here is a valid lower bound for all of them. Emits
// CPLEX .lp; solves with cbc/glpsol if present, else brute-force for tiny
// instances.
//
// Build:  make ilp_baseline
// Usage:  ./ilp_baseline "<TEXT>" "<SHAPE>" [--max-add 8] [--lp-out path.lp]
//                        [--universe full|maximal|legacy]

#include "compress.hpp"

#include <unistd.h>
#include <fcntl.h>

#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace gcsa;
using Clock = std::chrono::steady_clock;

struct IlpSolution {
    bool solved = false;
    std::string method;          // cbc / glpsol / brute / none
    int64_t opt_C = -1;
    std::vector<uint8_t> y;      // selected candidates
    std::vector<uint8_t> x;      // kept ranks
    double solve_ms = 0;
    int n_vars = 0;
    int n_cons = 0;
    std::string message;
};

static size_t heuristic_C(const Shape& sh, const std::string& text,
                          int max_add, CompressAlgo algo) {
    // compress.hpp logs progress to stderr; silence for clean CLI output.
    fflush(stderr);
    int saved = dup(STDERR_FILENO);
    int nullfd = open("/dev/null", O_WRONLY);
    if (saved >= 0 && nullfd >= 0) dup2(nullfd, STDERR_FILENO);
    if (nullfd >= 0) close(nullfd);

    CompressedIndex idx;
    idx.build(sh, text, max_add, algo);
    size_t C = idx.stored_positions();

    if (saved >= 0) {
        fflush(stderr);
        dup2(saved, STDERR_FILENO);
        close(saved);
    }
    return C;
}

static std::string find_on_path(const char* name) {
    const char* path = std::getenv("PATH");
    if (!path) return {};
    std::stringstream ss(path);
    std::string dir;
    while (std::getline(ss, dir, ':')) {
        if (dir.empty()) continue;
        std::string p = dir + "/" + name;
        if (::access(p.c_str(), X_OK) == 0) return p;
    }
    return {};
}

// Write CPLEX LP for the differential source-selection ILP.
static void write_lp(const std::string& path,
                     const CompressedIndex::DiffProblem& P,
                     int& n_vars, int& n_cons) {
    const size_t m = P.m;
    const size_t E = P.candidates.size();
    n_vars = (int)(m + E);
    n_cons = 0;

    std::ofstream out(path);
    // CBC's LP reader rejects CPLEX-style '\' comment lines ("Unknown image").
    // Keep the file comment-free so both cbc and glpsol can load it.
    out << "Minimize\n obj:";
    for (size_t r = 0; r < m; ++r) {
        if (r % 16 == 0) out << "\n ";
        out << " + x" << r;
    }
    out << "\nSubject To\n";

    // (1) ≤1 candidate per target name
    std::map<uint64_t, std::vector<size_t>> by_name;
    for (size_t e = 0; e < E; ++e)
        by_name[P.candidates[e].name].push_back(e);
    for (const auto& kv : by_name) {
        out << " one_" << kv.first << ":";
        for (size_t e : kv.second) out << " + y" << e;
        out << " <= 1\n";
        ++n_cons;
    }

    // (2) y_e=1 & r covered ⇒ x_r=0   →  x_r + y_e <= 1
    for (size_t e = 0; e < E; ++e) {
        for (int32_t r : P.candidates[e].covered) {
            out << " drop_e" << e << "_r" << r << ": x" << r << " + y" << e << " <= 1\n";
            ++n_cons;
        }
    }

    // (3) y_e=1 & r in source ⇒ x_r=1  →  x_r - y_e >= 0
    for (size_t e = 0; e < E; ++e) {
        const auto& c = P.candidates[e];
        for (int32_t r = c.src_lo; r < c.src_hi; ++r) {
            out << " pin_e" << e << "_r" << r << ": x" << r << " - y" << e << " >= 0\n";
            ++n_cons;
        }
    }

    // (4) no cover ⇒ keep: x_r + sum_{e: r covered} y_e >= 1
    std::vector<std::vector<size_t>> covers(m);
    for (size_t e = 0; e < E; ++e)
        for (int32_t r : P.candidates[e].covered)
            covers[(size_t)r].push_back(e);
    for (size_t r = 0; r < m; ++r) {
        out << " keep_r" << r << ": x" << r;
        for (size_t e : covers[r]) out << " + y" << e;
        out << " >= 1\n";
        ++n_cons;
    }

    out << "Bounds\nBinary\n";
    for (size_t r = 0; r < m; ++r) {
        if (r % 16 == 0) out << " ";
        out << " x" << r;
        if (r % 16 == 15) out << "\n";
    }
    out << "\n";
    for (size_t e = 0; e < E; ++e) {
        if (e % 16 == 0) out << " ";
        out << " y" << e;
        if (e % 16 == 15) out << "\n";
    }
    out << "\nEnd\n";
}

static bool parse_cbc_sol(const std::string& path, size_t m, size_t E,
                          IlpSolution& sol) {
    std::ifstream in(path);
    if (!in) return false;
    sol.x.assign(m, 0);
    sol.y.assign(E, 0);
    std::string line;
    bool seen_obj = false;
    while (std::getline(in, line)) {
        // CBC solu: "Optimal - objective value N" or "Integer ..."
        if (line.find("objective value") != std::string::npos) {
            auto pos = line.find_last_of(" \t");
            if (pos != std::string::npos) {
                sol.opt_C = (int64_t)std::llround(std::stod(line.substr(pos + 1)));
                seen_obj = true;
            }
        }
        // Variable lines: "      0 x12               1" or "x12 1"
        std::istringstream iss(line);
        std::string a, b, c;
        if (!(iss >> a)) continue;
        // Skip index if numeric
        std::string name, val;
        if (!a.empty() && std::isdigit((unsigned char)a[0])) {
            if (!(iss >> name >> val)) continue;
        } else {
            name = a;
            if (!(iss >> val)) continue;
        }
        double v = 0;
        try { v = std::stod(val); } catch (...) { continue; }
        if (v < 0.5) continue;
        if (name.size() >= 2 && name[0] == 'x') {
            size_t r = (size_t)std::stoul(name.substr(1));
            if (r < m) sol.x[r] = 1;
        } else if (name.size() >= 2 && name[0] == 'y') {
            size_t e = (size_t)std::stoul(name.substr(1));
            if (e < E) sol.y[e] = 1;
        }
    }
    if (!seen_obj && sol.opt_C < 0) {
        int64_t sum = 0;
        for (auto v : sol.x) sum += v;
        sol.opt_C = sum;
    }
    return seen_obj || sol.opt_C >= 0;
}

static bool parse_glpsol_sol(const std::string& path, size_t m, size_t E,
                             IlpSolution& sol) {
    // glpsol -o plain text. Columns section lines look like:
    //   "     1 x0           *              1             0             1"
    //   No.  name  [*]  Activity  Lower  Upper
    // Must take Activity (first number after name/*), not the trailing Upper
    // bound — otherwise every binary var with UB=1 is read as selected.
    std::ifstream in(path);
    if (!in) return false;
    sol.x.assign(m, 0);
    sol.y.assign(E, 0);
    std::string line;
    bool got = false;
    bool in_columns = false;
    while (std::getline(in, line)) {
        if (line.find("Objective:") != std::string::npos) {
            auto eq = line.find('=');
            if (eq != std::string::npos) {
                std::istringstream iss(line.substr(eq + 1));
                double v;
                if (iss >> v) { sol.opt_C = (int64_t)std::llround(v); got = true; }
            }
        }
        if (line.find("Column name") != std::string::npos) {
            in_columns = true;
            continue;
        }
        if (!in_columns) continue;
        if (line.find("-----") != std::string::npos) continue;
        if (line.empty() || line.find("Integer") != std::string::npos
            || line.find("KKT.") != std::string::npos
            || line.find("End of") != std::string::npos) {
            if (!line.empty()) in_columns = false;
            continue;
        }
        std::istringstream iss(line);
        std::string tok;
        std::vector<std::string> toks;
        while (iss >> tok) toks.push_back(tok);
        if (toks.size() < 2) continue;

        size_t i = 0;
        // Optional leading column index
        if (!toks[i].empty() && std::isdigit((unsigned char)toks[i][0])) ++i;
        if (i >= toks.size()) continue;
        const std::string& name = toks[i++];
        if (name.empty() || (name[0] != 'x' && name[0] != 'y')) continue;
        // Optional '*' status marker
        if (i < toks.size() && toks[i] == "*") ++i;
        if (i >= toks.size()) continue;
        double v = 0;
        try { v = std::stod(toks[i]); } catch (...) { continue; }
        if (v < 0.5) continue;
        if (name[0] == 'x') {
            size_t r = (size_t)std::stoul(name.substr(1));
            if (r < m) sol.x[r] = 1;
        } else {
            size_t e = (size_t)std::stoul(name.substr(1));
            if (e < E) sol.y[e] = 1;
        }
    }
    if (!got) {
        int64_t sum = 0;
        for (auto v : sol.x) sum += v;
        if (sum > 0) { sol.opt_C = sum; got = true; }
    }
    return got;
}

// Evaluate |C| for a 0/1 selection of candidates (at most one per name).
static int64_t eval_keep(const CompressedIndex::DiffProblem& P,
                         const std::vector<uint8_t>& y_sel,
                         std::vector<uint8_t>* x_out = nullptr) {
    const size_t m = P.m;
    std::vector<uint8_t> covered(m, 0);
    std::vector<uint8_t> pinned(m, 0);
    std::set<uint64_t> used_name;
    for (size_t e = 0; e < P.candidates.size(); ++e) {
        if (!y_sel[e]) continue;
        const auto& c = P.candidates[e];
        if (!used_name.insert(c.name).second) return -1; // >1 per name
        for (int32_t r : c.covered) covered[(size_t)r] = 1;
        for (int32_t r = c.src_lo; r < c.src_hi; ++r) pinned[(size_t)r] = 1;
    }
    // Feasibility: source ranks must not be covered (dropped).
    for (size_t r = 0; r < m; ++r)
        if (pinned[r] && covered[r]) return -1;

    std::vector<uint8_t> x(m, 0);
    int64_t sum = 0;
    for (size_t r = 0; r < m; ++r) {
        // keep iff not covered (constraint 4 + 2)
        x[r] = covered[r] ? 0 : 1;
        // constraint 3 already enforced by pinned⊆keep when covered∩pinned=∅
        sum += x[r];
    }
    if (x_out) *x_out = std::move(x);
    return sum;
}

static bool self_check(const CompressedIndex::DiffProblem& P, const IlpSolution& sol,
                       std::string& err) {
    if (!sol.solved || sol.y.size() != P.candidates.size() || sol.x.size() != P.m) {
        err = "solution size mismatch";
        return false;
    }
    // (1) at most one y per name
    std::map<uint64_t, int> cnt;
    for (size_t e = 0; e < sol.y.size(); ++e)
        if (sol.y[e]) ++cnt[P.candidates[e].name];
    for (auto& kv : cnt) if (kv.second > 1) {
        err = "more than one candidate for a name";
        return false;
    }
    for (size_t e = 0; e < sol.y.size(); ++e) {
        if (!sol.y[e]) continue;
        const auto& c = P.candidates[e];
        for (int32_t r : c.covered)
            if (sol.x[(size_t)r]) { err = "covered rank kept"; return false; }
        for (int32_t r = c.src_lo; r < c.src_hi; ++r)
            if (!sol.x[(size_t)r]) { err = "source rank not kept"; return false; }
    }
    // (4) + derive keep from covers
    std::vector<uint8_t> any(P.m, 0);
    for (size_t e = 0; e < sol.y.size(); ++e)
        if (sol.y[e])
            for (int32_t r : P.candidates[e].covered) any[(size_t)r] = 1;
    int64_t sum = 0;
    for (size_t r = 0; r < P.m; ++r) {
        if (!any[r] && !sol.x[r]) { err = "uncovered rank dropped"; return false; }
        if (any[r] && sol.x[r]) { err = "covered rank kept (4)"; return false; }
        sum += sol.x[r];
    }
    if (sum != sol.opt_C) {
        err = "opt_C != sum x_r";
        return false;
    }
    // Cross-check eval_keep
    std::vector<uint8_t> x2;
    int64_t eC = eval_keep(P, sol.y, &x2);
    if (eC < 0 || eC != sol.opt_C || x2 != sol.x) {
        err = "eval_keep inconsistency";
        return false;
    }
    return true;
}

// End-to-end check: materialize the ILP solution as a real index and decode
// every k-mer against brute force. The universe drops the note's lcp >= add+1
// rule, so this is what actually proves the extra links are decodable.
static bool decode_check(CompressedIndex& helper,
                         const CompressedIndex::DiffProblem& P,
                         const IlpSolution& sol,
                         const Shape& sh, const std::string& text,
                         size_t& C_out, std::string& err) {
    std::vector<CompressedIndex::Candidate> chosen;
    for (size_t e = 0; e < sol.y.size(); ++e)
        if (sol.y[e]) chosen.push_back(P.candidates[e]);
    if (!helper.apply_links(P.intervals, chosen)) {
        err = "link set not simultaneously feasible";
        return false;
    }
    C_out = helper.stored_positions();
    if ((int64_t)C_out != sol.opt_C) {
        err = "|C| after apply_links (" + std::to_string(C_out)
              + ") != opt_C (" + std::to_string(sol.opt_C) + ")";
        return false;
    }
    std::map<uint64_t, std::vector<int64_t>> truth;
    for (long p = 0; p <= (long)text.size(); ++p)
        truth[name_at(sh, text, (size_t)p)].push_back(p);
    for (auto& kv : truth) {
        auto got = helper.positions_of(kv.first);
        std::sort(got.begin(), got.end());
        auto exp = kv.second;
        std::sort(exp.begin(), exp.end());
        if (got != exp) {
            err = "positions_of mismatch for " + name_to_string(sh, kv.first);
            return false;
        }
    }
    return true;
}

static IlpSolution brute_force(const CompressedIndex::DiffProblem& P) {
    IlpSolution sol;
    sol.method = "brute";
    const size_t E = P.candidates.size();
    // Group candidates by name; choose at most one index per group (+ none).
    std::map<uint64_t, std::vector<size_t>> by_name;
    for (size_t e = 0; e < E; ++e)
        by_name[P.candidates[e].name].push_back(e);
    std::vector<std::vector<size_t>> groups;
    groups.reserve(by_name.size());
    for (auto& kv : by_name) groups.push_back(std::move(kv.second));

    // Product of (1+|g|) choices; abort if too large.
    uint64_t combos = 1;
    for (auto& g : groups) {
        uint64_t f = 1 + (uint64_t)g.size();
        if (combos > (uint64_t)1e7 / f) {
            sol.message = "brute-force too large (" + std::to_string(groups.size())
                          + " names, " + std::to_string(E) + " cands)";
            return sol;
        }
        combos *= f;
    }

    auto t0 = Clock::now();
    sol.y.assign(E, 0);
    sol.x.assign(P.m, 1);
    sol.opt_C = (int64_t)P.m;
    std::vector<int> choice(groups.size(), -1); // -1 = none, else index into group

    std::function<void(size_t)> rec = [&](size_t gi) {
        if (gi == groups.size()) {
            std::vector<uint8_t> y(E, 0);
            for (size_t i = 0; i < groups.size(); ++i)
                if (choice[i] >= 0) y[groups[i][(size_t)choice[i]]] = 1;
            std::vector<uint8_t> x;
            int64_t C = eval_keep(P, y, &x);
            if (C >= 0 && C < sol.opt_C) {
                sol.opt_C = C;
                sol.y = std::move(y);
                sol.x = std::move(x);
            }
            return;
        }
        choice[gi] = -1;
        rec(gi + 1);
        for (size_t k = 0; k < groups[gi].size(); ++k) {
            choice[gi] = (int)k;
            rec(gi + 1);
        }
    };
    rec(0);
    sol.solve_ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    sol.solved = true;
    sol.n_vars = (int)(P.m + E);
    return sol;
}

static IlpSolution solve_with_external(const std::string& lp_path,
                                       const CompressedIndex::DiffProblem& P,
                                       int n_vars, int n_cons) {
    IlpSolution sol;
    sol.n_vars = n_vars;
    sol.n_cons = n_cons;
    const size_t m = P.m, E = P.candidates.size();

    std::string cbc = find_on_path("cbc");
    std::string glp = find_on_path("glpsol");

    auto run = [&](const std::string& cmd, const std::string& sol_path,
                   bool (*parse)(const std::string&, size_t, size_t, IlpSolution&),
                   const char* method) -> bool {
        auto t0 = Clock::now();
        int rc = std::system(cmd.c_str());
        sol.solve_ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
        (void)rc;
        if (!parse(sol_path, m, E, sol)) return false;
        // Prefer sum of x if objective parse flaky
        int64_t sum = 0;
        for (auto v : sol.x) sum += v;
        if (sum > 0) sol.opt_C = sum;
        sol.solved = true;
        sol.method = method;
        return true;
    };

    if (!cbc.empty()) {
        std::string sol_path = lp_path + ".cbc.sol";
        // Quiet solve; write solution file
        std::string cmd = "\"" + cbc + "\" \"" + lp_path
            + "\" -solve -solu \"" + sol_path + "\" >/dev/null 2>&1";
        if (run(cmd, sol_path, parse_cbc_sol, "cbc")) return sol;
        sol.message = "cbc ran but solution parse failed";
    }
    if (!glp.empty()) {
        std::string sol_path = lp_path + ".glpk.sol";
        std::string cmd = "\"" + glp + "\" --lp \"" + lp_path
            + "\" -o \"" + sol_path + "\" >/dev/null 2>&1";
        if (run(cmd, sol_path, parse_glpsol_sol, "glpsol")) return sol;
        if (sol.message.empty())
            sol.message = "glpsol ran but solution parse failed";
    }

    // Fallback brute for small instances
    IlpSolution br = brute_force(P);
    br.n_vars = n_vars;
    br.n_cons = n_cons;
    if (br.solved) return br;

    sol.method = "none";
    sol.solved = false;
    if (sol.message.empty())
        sol.message = "no MIP solver on PATH (tried cbc, glpsol); "
                      "brute-force infeasible; wrote LP only";
    if (!br.message.empty())
        sol.message += "; " + br.message;
    return sol;
}

static void usage(const char* argv0) {
    std::cerr << "Usage: " << argv0
              << " \"<TEXT>\" \"<SHAPE>\" [--max-add 8] [--lp-out path.lp]\n"
              << "       [--universe full|maximal|legacy]   (= GCSA_LINK_UNIVERSE)\n"
              << "  full    every decodable link incl. non-maximal source runs;\n"
              << "          the only setting that is a true lower bound for all\n"
              << "          algorithms, but O(L^2) links per run of length L\n"
              << "  maximal maximal source runs only\n"
              << "  legacy  the heuristics' enumerate_candidates_ view\n";
}

int main(int argc, char** argv) {
    if (argc < 3) { usage(argv[0]); return 1; }
    std::string text = argv[1];
    std::string shape_str = argv[2];
    int max_add = 8;
    std::string lp_out;

    for (int i = 3; i < argc; ++i) {
        std::string k = argv[i];
        auto need = [&](const char* flag) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "missing value for " << flag << "\n";
                std::exit(1);
            }
            return argv[++i];
        };
        if (k == "--max-add") max_add = std::stoi(need("--max-add"));
        else if (k == "--lp-out") lp_out = need("--lp-out");
        else if (k == "--universe") ::setenv("GCSA_LINK_UNIVERSE",
                                             need("--universe").c_str(), 1);
        else {
            std::cerr << "unknown arg: " << k << "\n";
            usage(argv[0]);
            return 1;
        }
    }

    Shape sh = Shape::parse(shape_str);
    CompressedIndex helper;
    auto P = helper.collect_diff_problem(sh, text, max_add);

    if (lp_out.empty()) {
        // Sanitize shape only — keep the ".lp" suffix so CBC auto-detects LP format.
        // (Replacing '.' in the whole string used to yield "ilp_____lp" for "#.#".)
        std::string stem = "ilp_" + shape_str;
        for (char& c : stem)
            if (c == '#' || c == '.' || c == '/') c = '_';
        lp_out = stem + ".lp";
    }

    int n_vars = 0, n_cons = 0;
    write_lp(lp_out, P, n_vars, n_cons);

    std::cout << "=== ILP baseline ===\n"
              << "text_len=" << text.size()
              << " shape=" << shape_str
              << " max_add=" << max_add
              << " kMinCoverage=" << kMinCoverage << "\n"
              << "m=" << P.m
              << " intervals=" << P.intervals.size()
              << " universe=" << P.universe
              << " candidates=" << P.candidates.size() << "\n"
              << "LP written: " << lp_out
              << "  (#vars=" << n_vars << " #constraints=" << n_cons << ")\n";

    // Heuristics
    size_t Cg = heuristic_C(sh, text, max_add, CompressAlgo::Greedy);
    size_t Cd = heuristic_C(sh, text, max_add, CompressAlgo::DepOrder);
    size_t Ct = heuristic_C(sh, text, max_add, CompressAlgo::TreeDp);
    size_t Ct3 = heuristic_C(sh, text, max_add, CompressAlgo::TreeDp3);
    size_t Ct4 = heuristic_C(sh, text, max_add, CompressAlgo::TreeDp4);

    IlpSolution sol = solve_with_external(lp_out, P, n_vars, n_cons);

    if (!sol.solved) {
        std::cout << "solver: NOT FOUND / failed\n"
                  << "  " << sol.message << "\n"
                  << "heuristic |C|: greedy=" << Cg
                  << " dep-order=" << Cd
                  << " tree-dp=" << Ct
                  << " tree-dp3=" << Ct3
                  << " tree-dp4=" << Ct4 << "\n"
                  << "Open " << lp_out << " with cbc/glpsol/gurobi to get optimal |C|.\n";
        return 2;
    }

    std::string err;
    bool ok = self_check(P, sol, err);
    std::string derr;
    size_t decoded_C = 0;
    bool dec_ok = decode_check(helper, P, sol, sh, text, decoded_C, derr);
    size_t n_sel = 0;
    for (auto v : sol.y) n_sel += v;

    size_t best_h = std::min({Cg, Cd, Ct, Ct3, Ct4});
    bool bounds = (sol.opt_C >= 0 && (size_t)sol.opt_C <= best_h);

    auto gap = [&](size_t h) -> double {
        if (sol.opt_C <= 0) return 0;
        return 100.0 * ((double)h - (double)sol.opt_C) / (double)sol.opt_C;
    };

    std::cout << "solver: " << sol.method
              << "  time_ms=" << sol.solve_ms
              << "  #vars=" << sol.n_vars
              << "  #constraints=" << sol.n_cons << "\n"
              << "optimal |C| = " << sol.opt_C
              << "  (selected " << n_sel << " candidates)\n"
              << "heuristic |C|: greedy=" << Cg
              << " dep-order=" << Cd
              << " tree-dp=" << Ct
              << " tree-dp3=" << Ct3
              << " tree-dp4=" << Ct4 << "\n"
              << "gap vs opt (%): greedy=" << gap(Cg)
              << " dep-order=" << gap(Cd)
              << " tree-dp=" << gap(Ct)
              << " tree-dp3=" << gap(Ct3)
              << " tree-dp4=" << gap(Ct4) << "\n"
              << "self-check: " << (ok ? "OK" : ("FAIL (" + err + ")")) << "\n"
              << "decode-check: " << (dec_ok ? "OK" : ("FAIL (" + derr + ")"))
              << "  (|C|=" << decoded_C << ")\n"
              << "bounds all algos: " << (bounds ? "YES" : "NO")
              << "  (opt=" << sol.opt_C << " min heuristic=" << best_h << ")\n";

    return (ok && dec_ok && bounds) ? 0 : 3;
}
