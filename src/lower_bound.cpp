// lower_bound.cpp -- provable lower bounds on |C| that scale to the full
// benchmark inputs, plus an exporter for the compact exact ILP
// (tools/exact_ilp.py).
//
// The problem, as fixed by compress.hpp's "Link universe" (L1)-(L5) and solved
// exactly by ilp_baseline: choose at most one link per k-mer c; a link at add
// a reads a run of >= kMinCoverage consecutive SA ranks, all kept, whose
// a-successors all carry c, and drops those successors.  Minimise the kept
// count |C|.  Every bound below is the optimum of a relaxation of that problem
// restricted to adds 1..A (A = --max-add), so each is <= OPT(A) and hence <=
// the |C| of any run whose links all have add <= A.  (Phase II's retarget pass
// can compose adds beyond A -- see the `add > cap` check in tools/; a run that
// does so is outside OPT(A)'s feasible set and must be compared with a bound
// computed at a larger --max-add.)
//
//   LB0  words with < kMinCoverage occurrences can never be covered.
//        m - sum_{|I_c| >= 3} |I_c|.
//   LB1  one H_offset per word: a word loses at most cap_c positions, where
//        cap_c is its longest source run over adds 1..A (the maximal runs of
//        enumerate_all_links_, no availability, no clash test).
//        m - sum_c cap_c.  Implies LB0.
//   DOM  a covered position p needs a *kept* position p - a*span (a <= A) in
//        its own residue chain, and only adds at which p is coverable at all
//        count.  On each chain that is a covering problem on a path, solved
//        exactly by a sliding-window DP.  Strong on repeats (where LB1 says
//        "almost everything can go"), useless on unique sequence.
//   LAG  Lagrangian relaxation combining both: dualise the per-word caps of
//        LB1 into the chain DP of DOM.  Any multipliers lambda in [0,1] give a
//        valid bound; lambda = 1 reproduces LB1 and lambda = 0 reproduces DOM
//        exactly, and projected subgradient ascent then improves on both.
//
// Usage:
//   ./lower_bound -g genome.fasta -s SHAPE [--max-add 8[,16,32...]] [--ub U] [--iters 300]
//                 [--sa-load blob.gsa | --sa-cache DIR] [--dump-runs runs.bin] [--text ACGT...]
// Prints one `key=value` line per requested max-add (the SA is built once; the
// add scan is cumulative, so 8,16,32 costs the same as 32 alone plus the DPs).
// --dump-runs writes the runs for adds 1..max(--max-add) for tools/exact_ilp.py.

#include "compress.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace gcsa;
using Clock = std::chrono::steady_clock;

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

static double ms_since(Clock::time_point t) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t).count();
}

int main(int argc, char** argv) {
    std::string genome, shape_str, sa_load, sa_cache, dump_runs, text_arg;
    std::vector<int> As = {8};
    long long ub = -1;
    int iters = 300;
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto need = [&]() -> std::string {
            if (i + 1 >= argc) { std::cerr << "missing value for " << k << "\n"; std::exit(1); }
            return argv[++i];
        };
        if (k == "-g") genome = need();
        else if (k == "-s" || k == "--shape") shape_str = need();
        else if (k == "--max-add") {
            As.clear();
            std::string v = need(), tok;
            for (char ch : v + ",") { if (ch == ',') { if (!tok.empty()) As.push_back(std::stoi(tok)); tok.clear(); } else tok.push_back(ch); }
            std::sort(As.begin(), As.end());
        }
        else if (k == "--ub") ub = std::stoll(need());
        else if (k == "--iters") iters = std::stoi(need());
        else if (k == "--sa-load") sa_load = need();
        else if (k == "--sa-cache") sa_cache = need();
        else if (k == "--dump-runs") dump_runs = need();
        else if (k == "--text") text_arg = need();
        else if (k == "--min-coverage") kMinCoverage = std::stoi(need());
        else { std::cerr << "unknown arg " << k << "\n"; return 1; }
    }
    if (shape_str.empty()) { std::cerr << "need -s SHAPE\n"; return 1; }
    Shape sh = Shape::parse(shape_str);
    const int MC = kMinCoverage;
    const int A = As.back();   // scan adds 1..max(As); report a bound at every requested A

    auto t0 = Clock::now();
    GappedSA G;
    if (!sa_cache.empty() && sa_load.empty()) {
        // same content-addressed blob name gcsa --sa-cache uses; never builds/writes
        if (genome.empty()) { std::cerr << "--sa-cache needs -g to derive the blob name\n"; return 1; }
        sa_load = sa_cache + "/" + gsa_blob_name(gcsa_text_hash(load_dna(genome)), sh.pattern);
    }
    if (!sa_load.empty()) {
        std::string err;
        if (!gsa_load(G, sa_load, &sh.pattern, &err)) { std::cerr << "--sa-load: " << err << "\n"; return 1; }
    } else {
        std::string text = !text_arg.empty() ? text_arg : load_dna(genome);
        if (!text_arg.empty()) {
            std::string clean;
            for (char c : text) if (std::isalpha((unsigned char)c)) clean.push_back(base_char(base_value(c)));
            text = clean;
        }
        G = build_gapped_sa(sh, std::move(text));
    }
    const double t_sa = ms_since(t0);
    std::vector<int32_t>().swap(G.lcp);   // not needed

    const size_t m = G.m();
    t0 = Clock::now();
    std::vector<int32_t> rank_of(m);
    for (size_t r = 0; r < m; ++r) rank_of[(size_t)G.sa[r]] = (int32_t)r;

    // word id per rank (interval index in SA order = name order)
    std::vector<int32_t> word_of_rank(m);
    std::vector<int64_t> isize;   // |I_c|
    {
        int32_t w = -1;
        uint64_t prev = 0;
        for (size_t r = 0; r < m; ++r) {
            uint64_t nm = G.lex[(size_t)G.sa[r]];
            if (r == 0 || nm != prev) { ++w; isize.push_back(0); prev = nm; }
            word_of_rank[r] = w;
            ++isize[(size_t)w];
        }
    }
    const size_t K = isize.size();
    std::vector<uint64_t>().swap(G.lex);
    std::string().swap(G.text);

    std::vector<int32_t> word_of_lex(m);
    for (size_t r = 0; r < m; ++r) word_of_lex[(size_t)G.sa[r]] = word_of_rank[r];
    const int64_t span = sh.span;
    // chain boundaries (lextext is grouped by residue: one block per residue)
    std::vector<size_t> chain_start;
    for (size_t i = 0; i < m; ++i)
        if (i == 0 || G.lex2orig[i] != G.lex2orig[i - 1] + span) chain_start.push_back(i);
    chain_start.push_back(m);

    int64_t elig = 0;
    for (size_t c = 0; c < K; ++c) if (isize[c] >= MC) elig += isize[c];
    const int64_t LB0 = (int64_t)m - elig;

    // succ in lextext: y = x + a valid iff lex2orig[y] == lex2orig[x] + a*span
    auto succ = [&](int64_t x, int a) -> int64_t {
        int64_t y = x + a;
        if (y >= (int64_t)m) return -1;
        if (G.lex2orig[(size_t)y] != G.lex2orig[(size_t)x] + (int64_t)a * span) return -1;
        return y;
    };

    std::vector<int64_t> cap(K, 0);
    std::vector<uint16_t> h(m, 0);              // max coverable add, per lextext pos (as target)
    std::ofstream dump;
    if (!dump_runs.empty()) {
        dump.open(dump_runs, std::ios::binary);
        int64_t hdr[3] = {(int64_t)m, (int64_t)K, (int64_t)A};
        dump.write((const char*)hdr, sizeof(hdr));
    }

    // ---- chain DP ---------------------------------------------------------
    // cost(keep)=1, cost(cover p)=lambda[word(p)]; covered p (local index k)
    // needs a kept j with k - h_k <= j < k.  Adds the min cost to `total` and
    // fills covcount (covered positions per word) from the argmin.
    std::vector<int32_t> pred(m);
    std::vector<double> f(m), S(m);
    std::vector<int64_t> covcount(K);
    auto solve = [&](const std::vector<float>& lam, double& total) {
        std::fill(covcount.begin(), covcount.end(), 0);
        total = 0;
        for (size_t ci = 0; ci + 1 < chain_start.size(); ++ci) {
            const size_t b0 = chain_start[ci], b1 = chain_start[ci + 1];
            const int64_t L = (int64_t)(b1 - b0);
            double acc = 0;
            for (int64_t k = 0; k < L; ++k) { acc += lam[(size_t)word_of_lex[b0 + k]]; S[b0 + k] = acc; }
            auto Sv = [&](int64_t j) -> double { return j < 0 ? 0.0 : S[b0 + (size_t)j]; };
            auto fv = [&](int64_t j) -> double { return j < 0 ? 0.0 : f[b0 + (size_t)j]; };
            std::deque<int64_t> dq;   // candidate j (incl. -1), increasing key fv(j)-Sv(j)
            int64_t lo = -1;
            dq.push_back(-1);
            for (int64_t k = 0; k <= L; ++k) {
                if (k >= 1) {
                    const int64_t g = (k - 1) - (int64_t)h[b0 + (size_t)(k - 1)];
                    lo = std::max(lo, g);
                    const double key = fv(k - 1) - Sv(k - 1);
                    while (!dq.empty() && fv(dq.back()) - Sv(dq.back()) >= key) dq.pop_back();
                    dq.push_back(k - 1);
                }
                while (!dq.empty() && dq.front() < lo) dq.pop_front();
                const int64_t j = dq.front();
                const double val = fv(j) - Sv(j) + Sv(k - 1);
                if (k < L) { f[b0 + (size_t)k] = 1.0 + val; pred[b0 + (size_t)k] = (int32_t)j; }
                else {
                    total += val;
                    int64_t cur = L, pj = j;
                    while (true) {
                        for (int64_t t = pj + 1; t < cur; ++t) ++covcount[(size_t)word_of_lex[b0 + (size_t)t]];
                        if (pj < 0) break;
                        cur = pj;
                        pj = pred[b0 + (size_t)pj];
                    }
                }
            }
        }
    };

    auto evaluate = [&](int Acur, int64_t runs_so_far) {
        auto te = Clock::now();
        int64_t sumcap = 0;
        for (size_t c = 0; c < K; ++c) sumcap += cap[c];
        const int64_t LB1 = (int64_t)m - sumcap;
        std::vector<float> lam(K, 0.0f);
        double tot = 0;
        solve(lam, tot);
        const double DOM = tot;   // lambda = 0
        std::fill(lam.begin(), lam.end(), 1.0f);
        solve(lam, tot);
        double L_at1 = tot;       // lambda = 1 must reproduce LB1
        for (size_t c = 0; c < K; ++c) L_at1 -= (double)cap[c];
        double best = std::max(DOM, L_at1);
        if (DOM > L_at1) std::fill(lam.begin(), lam.end(), 0.0f);
        double theta = 1.0;
        int since_improve = 0, it_done = 0;
        const double target = (ub > 0) ? (double)ub : (double)m;
        for (int it = 0; it < iters; ++it) {
            solve(lam, tot);
            double Lval = tot, gnorm = 0;
            for (size_t c = 0; c < K; ++c) {
                Lval -= (double)lam[c] * (double)cap[c];
                double g = (double)covcount[c] - (double)cap[c];
                if ((lam[c] <= 0 && g < 0) || (lam[c] >= 1 && g > 0)) continue;
                gnorm += g * g;
            }
            ++it_done;
            if (Lval > best + 1e-9) { best = Lval; since_improve = 0; }
            else if (++since_improve >= 15) { theta *= 0.5; since_improve = 0; }
            if (gnorm < 1e-12 || theta < 1e-4) break;
            const double step = theta * std::max(1.0, target - Lval) / gnorm;
            for (size_t c = 0; c < K; ++c) {
                double v = (double)lam[c] + step * ((double)covcount[c] - (double)cap[c]);
                lam[c] = (float)std::min(1.0, std::max(0.0, v));
            }
        }
        const int64_t LAG = (int64_t)std::ceil(best - 1e-6);
        const int64_t DOMi = (int64_t)std::ceil(DOM - 1e-6);
        int64_t coverable = 0;
        for (size_t i = 0; i < m; ++i) coverable += (h[i] > 0);
        std::cout << "m=" << m << " K=" << K << " max_add=" << Acur << " min_cov=" << MC
                  << " LB0=" << LB0 << " LB1=" << LB1 << " DOM=" << DOMi << " LAG=" << LAG
                  << " LB=" << std::max({LB0, LB1, DOMi, LAG})
                  << " lam1_check=" << (int64_t)std::llround(L_at1)
                  << " coverable=" << coverable << " runs=" << runs_so_far
                  << " lag_iters=" << it_done
                  << " t_sa_ms=" << (int64_t)t_sa << " t_scan_ms=" << (int64_t)ms_since(t0)
                  << " t_lag_ms=" << (int64_t)ms_since(te) << std::endl;
    };

    std::vector<int32_t> tgt(m);
    int64_t runs_total = 0;
    size_t next_report = 0;
    for (int a = 1; a <= A; ++a) {
        for (size_t r = 0; r < m; ++r) {
            int64_t y = succ(G.sa[r], a);
            tgt[r] = (y < 0) ? -1 : rank_of[(size_t)y];
        }
        size_t s = 0;
        while (s < m) {
            if (tgt[s] < 0) { ++s; continue; }
            const int32_t w = word_of_rank[(size_t)tgt[s]];
            size_t e = s + 1;
            while (e < m && tgt[e] >= 0 && word_of_rank[(size_t)tgt[e]] == w) ++e;
            const int64_t L = (int64_t)(e - s);
            if (L >= MC) {
                ++runs_total;
                cap[(size_t)w] = std::max(cap[(size_t)w], std::min<int64_t>(L, isize[(size_t)w]));
                for (size_t r = s; r < e; ++r) {
                    int64_t y = G.sa[(size_t)tgt[r]];
                    if (h[(size_t)y] < a) h[(size_t)y] = (uint16_t)a;
                }
                if (dump) {
                    int64_t rec[4] = {(int64_t)a, (int64_t)w, (int64_t)s, L};
                    dump.write((const char*)rec, sizeof(rec));
                    for (size_t r = s; r < e; ++r) dump.write((const char*)&tgt[r], sizeof(int32_t));
                }
            }
            s = e;
        }
        while (next_report < As.size() && As[next_report] == a) { evaluate(a, runs_total); ++next_report; }
    }
    if (dump) {
        int64_t endrec[4] = {0, 0, 0, 0};
        dump.write((const char*)endrec, sizeof(endrec));
        dump.write((const char*)word_of_rank.data(), (std::streamsize)(m * sizeof(int32_t)));
        dump.close();
    }
    return 0;
}
