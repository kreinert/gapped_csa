// link_stats.cpp -- build one index exactly as gcsa does (same env vars) and
// report |C| plus what the final link set looks like, in particular the
// largest add any link uses. Phase II's retarget pass composes adds
// (d.add += cand.add), so a run at --max-add A can end with links of add > A;
// a lower bound is only valid for that run if computed at an add >= max_add_used.
//
//   ./link_stats FASTA|BLOB.gsa SHAPE|K MAX_ADD ALGO     (ALGO: gs|gd|dep|pfdp|pfdpi|rundp)
// SHAPE as a bare integer k runs gcsa's regular (plain-SA) mode.
// A path ending in .gsa is loaded as a prepared SA blob (gcsa --sa-cache / run_suite.py).
// prints: C= m= links= max_add_used= links_over_cap= positions_covered_by_over_cap=
#include "compress.hpp"
#include <fstream>
#include <iostream>
using namespace gcsa;
static std::string load_dna(const std::string& path) {
    std::ifstream in(path); std::string out, line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') continue;
        for (char c : line) if (std::isalpha((unsigned char)c)) out.push_back(base_char(base_value(c)));
    }
    return out;
}
int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "usage: link_stats FASTA|BLOB.gsa SHAPE MAX_ADD gs|gd|dep|pfdp|pfdpi|rundp\n"; return 1; }
    const std::string in = argv[1]; const std::string shs = argv[2]; int A = std::stoi(argv[3]);
    // a bare integer k = gcsa's regular mode (plain SA, one add = one character)
    const bool regular = shs.find_first_not_of("0123456789") == std::string::npos;
    Shape sh = Shape::parse(regular ? std::string((size_t)std::stoi(shs), '#') : shs);
    const bool is_blob = in.size() > 4 && in.compare(in.size() - 4, 4, ".gsa") == 0;
    std::string al = argv[4];
    CompressAlgo algo = al == "dep" ? CompressAlgo::DepOrder : al == "gs" ? CompressAlgo::GreedySize
                      : al == "rundp" ? CompressAlgo::RunDp
                      : al == "gd" ? CompressAlgo::GreedyDegree : al == "pfdp" ? CompressAlgo::PseudoforestDp
                      : CompressAlgo::PseudoforestDpIterate;
    gcsa_set_quiet(true);
    CompressedIndex idx;
    if (is_blob) {
        // blobs are keyed on the pattern, or "plain" for the (k-independent) plain SA
        const std::string key = regular ? std::string(kGsaPlainKey) : sh.pattern;
        GappedSA G; std::string err;
        if (!gsa_load(G, in, &key, &err)) { std::cerr << "blob: " << err << "\n"; return 1; }
        if (regular) idx.build_from_plain_sa(std::move(G), std::stoi(shs), A, algo);
        else         idx.build_from_gsa(std::move(G), A, algo);
    } else if (regular) {
        idx.build(std::stoi(shs), load_dna(in), A, algo);
    } else {
        idx.build(sh, load_dna(in), A, algo);
    }
    int mx = 0; size_t over = 0, over_pos = 0, nlinks = 0;
    for (auto& kv : idx.table()) {
        if (!kv.second.has_offset) continue;
        ++nlinks; mx = std::max<int>(mx, (int)kv.second.off_add);
        if ((int)kv.second.off_add > A) { ++over; over_pos += kv.second.off_num; }
    }
    std::cout << "C=" << idx.stored_positions() << " m=" << idx.total_positions() << " links=" << nlinks
              << " max_add_used=" << mx << " links_over_cap=" << over
              << " positions_covered_by_over_cap=" << over_pos << "\n";
}
