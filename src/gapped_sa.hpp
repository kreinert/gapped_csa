// gapped_sa.hpp
//
// Construction of the gapped suffix array for an arbitrary shape via the DisLex
// transform (Horton 2008).  The suffix sorting is done by a linear-time SA-IS
// for integer alphabets (see sais.hpp).
//
// DisLex in one paragraph.  For a shape of span s we form, for every text start
// position p, the gapped k-mer name at p (see shape.hpp).  We then lay these
// names out grouped by p mod s (mod 0 group first, then mod 1, ...).  The suffix
// array of this integer string ("lextext") is exactly the gapped suffix array:
// consecutive symbols within a mod group are the names at p, p+s, p+2s, ... i.e.
// exactly the gapped suffix read at p.  Following the note we do NOT insert
// separators between groups; boundary windows read into '$' padding, and the
// unique all-'$' name (= 0) at the very end acts as the global terminator.
//
// The lextext names are large integers (up to ~2^62), so we sort them with an
// integer-alphabet SA-IS.  Before sorting we order-preservingly remap the
// distinct names to compact ranks 0..k-1 (sort-unique + lower_bound) so the
// alphabet fed to SA-IS is dense and small; the original arithmetic names are
// kept in G.lex for hashing/lookup.  The LCP array is computed in *symbol* space
// with Kasai's algorithm over G.lex.

#pragma once

#include "shape.hpp"
#include "sais.hpp"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <unistd.h>

namespace gcsa {

struct GappedSA {
    Shape shape;
    std::string text;                 // cleaned DNA text (no padding)
    size_t n = 0;                     // text length
    int alphabet_size = 0;            // number of distinct lextext names (SA-IS alphabet)

    std::vector<uint64_t> lex;        // lextext names, length m = n+1 (DisLex order)
    std::vector<int64_t>  lex2orig;   // lextext position -> original text start position
    std::vector<int32_t>  sa;         // gapped SA: sa[r] = lextext position, length m
    std::vector<int32_t>  lcp;        // lcp[r] = lcp(sa[r-1], sa[r]) in symbols; lcp[0]=0

    size_t m() const { return lex.size(); }

    // Original text position of the r-th suffix in the gapped SA.
    int64_t orig_pos(int32_t rank) const { return lex2orig[sa[rank]]; }

    // First symbol (= gapped k-mer name) of the r-th suffix.
    uint64_t first_symbol(int32_t rank) const { return lex[sa[rank]]; }
};

// Kasai's algorithm: LCP array (in units of `lex` symbols) for a suffix array
// `sa` of the integer string `lex`. lcp[r] = shared-prefix length (in `lex`
// symbols) between the suffixes at sa[r-1] and sa[r]; lcp[0] = 0. Shared by
// both build_gapped_sa (symbol space = gapped k-mer names) and build_plain_sa
// (symbol space = raw characters) -- the algorithm itself doesn't care what
// the symbols mean.
inline std::vector<int32_t> kasai_lcp(const std::vector<uint64_t>& lex,
                                       const std::vector<int32_t>& sa) {
    const size_t m = lex.size();
    std::vector<int32_t> rank(m);
    for (size_t r = 0; r < m; ++r) rank[sa[r]] = (int32_t)r;
    std::vector<int32_t> lcp(m, 0);
    int h = 0;
    for (size_t i = 0; i < m; ++i) {
        if (rank[i] > 0) {
            size_t j = sa[rank[i] - 1];
            while (i + h < m && j + h < m && lex[i + h] == lex[j + h]) ++h;
            lcp[rank[i]] = h;
            if (h > 0) --h;
        } else {
            h = 0;
        }
    }
    return lcp;
}

// Plain (ungapped) suffix array of the raw text, built directly over the
// small {$,A,C,G,T} alphabet -- no per-position k-mer name computation, no
// mod-residue grouping, no alphabet remap. `k` only affects G.shape (used for
// naming k-mers of that width downstream, e.g. interval identity in
// compress.hpp); the SA/LCP themselves don't depend on k at all. This isolates
// the cost of DisLex's per-position name computation + big-alphabet sort from
// suffix-sorting itself: for a contiguous shape "######" (span==weight==k),
// build_gapped_sa produces the same total suffix order as this function.
inline GappedSA build_plain_sa(int k, std::string text) {
    using Clock = std::chrono::steady_clock;
    const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
    auto ms = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };

    GappedSA G;
    G.shape = Shape::parse(std::string((size_t)k, '#'));
    G.text  = std::move(text);
    G.n     = G.text.size();

    // 1) Map characters directly to a small dense alphabet: $ = 0 (unique
    //    sentinel at position n), A/C/G/T = 1..4. No name_at() calls.
    auto t0 = Clock::now();
    const size_t m = G.n + 1;
    G.lex.resize(m);
    G.lex2orig.resize(m);
    for (size_t i = 0; i < m; ++i) {
        G.lex[i] = (i < G.n) ? (uint64_t)(base_value(G.text[i]) + 1) : 0;
        G.lex2orig[i] = (int64_t)i;   // identity: no mod-grouping in plain mode
    }
    G.alphabet_size = SIGMA + 1;
    auto t1 = Clock::now();

    // 2) Suffix-sort the small-alphabet character string directly with SA-IS.
    std::vector<int32_t> chars(m);
    for (size_t i = 0; i < m; ++i) chars[i] = (int32_t)G.lex[i];
    G.sa = sais_int(chars, G.alphabet_size);
    auto t2 = Clock::now();

    // 3) Kasai LCP in character space.
    G.lcp = kasai_lcp(G.lex, G.sa);
    auto t3 = Clock::now();

    if (timing) {
        std::fprintf(stderr,
            "[timing] plain charmap=%.1fms  SA-IS=%.1fms  LCP=%.1fms  m=%zu sigma=%d\n",
            ms(t0, t1), ms(t1, t2), ms(t2, t3), m, G.alphabet_size);
    }
    return G;
}

inline GappedSA build_gapped_sa(const Shape& shape, std::string text) {
    using Clock = std::chrono::steady_clock;
    const bool timing = (std::getenv("GCSA_TIMING") != nullptr);
    auto ms = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };

    GappedSA G;
    G.shape = shape;
    G.text  = std::move(text);
    G.n     = G.text.size();

    // 1) Build the lextext in DisLex (grouped-by-residue) order.
    auto t0 = Clock::now();
    const int s = shape.span;
    G.lex.reserve(G.n + 1);
    G.lex2orig.reserve(G.n + 1);
    for (int mod = 0; mod < s; ++mod) {
        for (long p = mod; p <= (long)G.n; p += s) {
            G.lex.push_back(name_at(shape, G.text, (size_t)p));
            G.lex2orig.push_back(p);
        }
    }
    const size_t m = G.lex.size();
    auto t1 = Clock::now();

    // 2) Order-preservingly remap the distinct names to compact ranks 0..k-1.
    //    (This keeps the alphabet passed to SA-IS dense; G.lex keeps the
    //    original arithmetic names for hashing/lookup.)
    std::vector<uint64_t> uniq(G.lex);
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    G.alphabet_size = (int)uniq.size();

    std::vector<int32_t> remapped(m);
    for (size_t i = 0; i < m; ++i)
        remapped[i] = (int32_t)(std::lower_bound(uniq.begin(), uniq.end(), G.lex[i]) - uniq.begin());
    auto t2 = Clock::now();

    // 3) Suffix-sort the integer string directly with SA-IS.
    G.sa = sais_int(remapped, G.alphabet_size);
    auto t3 = Clock::now();

    // 4) Kasai LCP in symbol space.
    G.lcp = kasai_lcp(G.lex, G.sa);
    auto t4 = Clock::now();
    if (timing) {
        std::fprintf(stderr,
            "[timing] DisLex+name=%.1fms  remap=%.1fms  SA-IS=%.1fms  LCP=%.1fms  m=%zu sigma=%d\n",
            ms(t0, t1), ms(t1, t2), ms(t2, t3), ms(t3, t4), m, G.alphabet_size);
    }
    return G;
}

// ---------------------------------------------------------------------------
// Serializing a built GappedSA ("blob")
//
// Why: building the SA is 37% of a benchmark sweep's total compute and depends
// only on (text, shape) -- not on --algo, --max-add or --min-coverage -- so the
// 15 runs of one (dataset, shape) cell rebuild the identical array 15 times.
// Dumping it once and loading it back removes that, and removes the SA's
// machine-to-machine timing variance from the running-time comparison.
//
// This is a same-machine cache, not an interchange format: the arrays go down
// in host byte order and native widths. The header records enough to refuse
// anything it does not recognize, which matters more here than portability --
// silently loading the wrong SA would corrupt every downstream number without
// crashing anything.
//
// Layout:
//   magic       8   "GCSAGSA\0"
//   version     u32 format version; a bump invalidates every existing blob
//   endian      u32 0x01020304, refuses a foreign byte order
//   shape_len   u32 followed by that many bytes of the cache key: the shape
//               pattern for a gapped SA ("##.##"), or kGsaPlainKey ("plain")
//               for the plain character SA of gcsa's regular mode. The plain
//               SA does not depend on k (only the LCP >= k grouping applied
//               later does), so one plain blob per text serves every -s k, and
//               it never collides with a lextext blob of a contiguous k-mer.
//   n           u64 text length
//   m           u64 lex/sa/lcp length (== n+1)
//   alphabet    i64 alphabet_size
//   text_hash   u64 FNV-1a over the cleaned text
//   text        n bytes
//   lex         m * 8   (uint64)
//   lex2orig    m * 8   (int64)
//   sa          m * 4   (int32)
//   lcp         m * 4   (int32)
static constexpr uint32_t kGsaBlobVersion = 1;
static constexpr uint32_t kGsaBlobEndian  = 0x01020304u;

// FNV-1a. Used to name a blob after the text it was built from, so a cache
// entry follows the content rather than a filename -- the concat inputs are
// regenerated (fresh mtime, same bytes) before every sweep, so anything keyed
// on the path or timestamp would miss exactly where the cache is worth most.
inline uint64_t gcsa_fnv1a(const void* data, size_t len) {
    const unsigned char* p = (const unsigned char*)data;
    uint64_t h = 1469598103934665603ull;
    for (size_t i = 0; i < len; ++i) { h ^= p[i]; h *= 1099511628211ull; }
    return h;
}
inline uint64_t gcsa_text_hash(const std::string& t) {
    return gcsa_fnv1a(t.data(), t.size());
}

// A blob's cache key: the #/. pattern for a gapped SA, kGsaPlainKey for a plain
// SA (regular mode, any k). Kept distinct on purpose -- the two SAs index the
// same k-mers but are different arrays (lextext names vs raw characters), and
// loading one for the other would corrupt every number.
static const char* const kGsaPlainKey = "plain";
inline bool gsa_key_is_plain(const std::string& key) { return key == kGsaPlainKey; }

// Blob filename for a (text, key) pair. Both are folded in, so two shapes
// over one genome are different entries and can coexist in one directory.
inline std::string gsa_blob_name(uint64_t text_hash, const std::string& shape_pattern) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%016llx_%016llx.gsa",
                  (unsigned long long)text_hash,
                  (unsigned long long)gcsa_fnv1a(shape_pattern.data(), shape_pattern.size()));
    return std::string(buf);
}

namespace detail {
template <class T>
inline void gsa_put_vec(std::ostream& o, const std::vector<T>& v) {
    if (!v.empty()) o.write((const char*)v.data(), (std::streamsize)(v.size() * sizeof(T)));
}
template <class T>
inline bool gsa_get_vec(std::istream& i, std::vector<T>& v, size_t count) {
    v.resize(count);
    if (!count) return true;
    i.read((char*)v.data(), (std::streamsize)(count * sizeof(T)));
    return (size_t)i.gcount() == count * sizeof(T);
}
inline void gsa_put(std::ostream& o, const void* p, size_t n) { o.write((const char*)p, (std::streamsize)n); }
inline bool gsa_get(std::istream& i, void* p, size_t n) {
    i.read((char*)p, (std::streamsize)n);
    return (size_t)i.gcount() == n;
}
}  // namespace detail

// Bytes a blob for this SA will occupy -- for sizing a cache directory before
// filling it, and to validate a file's length on load.
inline uint64_t gsa_blob_bytes(uint64_t n, uint64_t m, uint32_t shape_len) {
    return 8 + 4 + 4 + 4 + (uint64_t)shape_len + 8 + 8 + 8 + 8
           + n + m * 8 + m * 8 + m * 4 + m * 4;
}

// Write `G` to `path`. The bytes go to a process-unique temp name in the same
// directory and are renamed into place, so a blob at `path` is always whole:
// a killed or out-of-disk writer leaves a .part behind rather than a truncated
// SA that the next run would load and silently compute against.
// `key` is the cache key to record; empty means the SA's own shape pattern
// (a gapped SA). A plain SA must pass kGsaPlainKey.
inline bool gsa_save(const GappedSA& G, const std::string& path, std::string* err,
                     const std::string& key = std::string()) {
    const std::string& rec = key.empty() ? G.shape.pattern : key;
    const std::string tmp = path + "." + std::to_string((long)::getpid()) + ".part";
    {
        std::ofstream o(tmp, std::ios::binary);
        if (!o) { if (err) *err = "cannot open " + tmp + " for writing"; return false; }
        const uint32_t ver = kGsaBlobVersion, endian = kGsaBlobEndian;
        const uint32_t slen = (uint32_t)rec.size();
        const uint64_t n = G.n, m = G.lex.size(), th = gcsa_text_hash(G.text);
        const int64_t alpha = G.alphabet_size;
        detail::gsa_put(o, "GCSAGSA\0", 8);
        detail::gsa_put(o, &ver, 4);
        detail::gsa_put(o, &endian, 4);
        detail::gsa_put(o, &slen, 4);
        detail::gsa_put(o, rec.data(), slen);
        detail::gsa_put(o, &n, 8);
        detail::gsa_put(o, &m, 8);
        detail::gsa_put(o, &alpha, 8);
        detail::gsa_put(o, &th, 8);
        detail::gsa_put(o, G.text.data(), G.text.size());
        detail::gsa_put_vec(o, G.lex);
        detail::gsa_put_vec(o, G.lex2orig);
        detail::gsa_put_vec(o, G.sa);
        detail::gsa_put_vec(o, G.lcp);
        o.flush();
        if (!o) { if (err) *err = "write failed (disk full?): " + tmp; std::remove(tmp.c_str()); return false; }
    }
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        if (err) *err = "cannot rename " + tmp + " -> " + path;
        std::remove(tmp.c_str());
        return false;
    }
    return true;
}

// Read a blob into `G`. Every header field is checked, and the file's length
// must equal exactly what the header implies -- a short file means a truncated
// or interrupted write, which must fail loudly rather than load a partial SA.
// `expect_shape`, when non-null, additionally requires the blob to have been
// built for that cache key (a shape pattern, or kGsaPlainKey); that
// is the check that stops a stale or mismatched cache entry from being
// compressed against the wrong pattern -- or a plain SA against a lextext one.
inline bool gsa_load(GappedSA& G, const std::string& path,
                     const std::string* expect_shape, std::string* err) {
    std::ifstream i(path, std::ios::binary);
    if (!i) { if (err) *err = "cannot open " + path; return false; }
    char magic[8];
    uint32_t ver = 0, endian = 0, slen = 0;
    if (!detail::gsa_get(i, magic, 8) || std::memcmp(magic, "GCSAGSA\0", 8) != 0) {
        { if (err) *err = path + ": not a gapped-SA blob";
          return false; }
    }
    if (!detail::gsa_get(i, &ver, 4) || ver != kGsaBlobVersion) {
        if (err) *err = path + ": format version " + std::to_string(ver) + ", this build writes "
                       + std::to_string(kGsaBlobVersion) + " -- rebuild the cache";
        return false;
    }
    if (!detail::gsa_get(i, &endian, 4) || endian != kGsaBlobEndian) {
        { if (err) *err = path + ": built on a machine with a different byte order";
          return false; }
    }
    if (!detail::gsa_get(i, &slen, 4) || slen == 0 || slen > 4096) {
        { if (err) *err = path + ": bad shape length";
          return false; }
    }
    std::string pattern(slen, '\0');
    if (!detail::gsa_get(i, &pattern[0], slen)) {
        { if (err) *err = path + ": truncated header";
          return false; }
    }
    if (expect_shape && pattern != *expect_shape) {
        auto what = [](const std::string& k) {
            return gsa_key_is_plain(k) ? std::string("a plain SA") : "shape " + k; };
        if (err) *err = path + ": blob was built for " + what(pattern) + ", not " + what(*expect_shape);
        return false;
    }
    uint64_t n = 0, m = 0, th = 0; int64_t alpha = 0;
    if (!detail::gsa_get(i, &n, 8) || !detail::gsa_get(i, &m, 8)
        || !detail::gsa_get(i, &alpha, 8) || !detail::gsa_get(i, &th, 8)) {
        { if (err) *err = path + ": truncated header";
          return false; }
    }
    if (m != n + 1) {
        { if (err) *err = path + ": inconsistent header (m != n+1)";
          return false; }
    }
    {   // a blob shorter than its own header describes is a partial write
        std::ifstream sz(path, std::ios::binary | std::ios::ate);
        const uint64_t actual = (uint64_t)sz.tellg();
        const uint64_t want = gsa_blob_bytes(n, m, slen);
        if (actual != want) {
            if (err) *err = path + ": size " + std::to_string(actual) + " bytes, header implies "
                           + std::to_string(want) + " -- truncated or corrupt";
            return false;
        }
    }
    G = GappedSA{};
    // A plain SA carries no k: the caller supplies it (build_from_plain_sa sets
    // the contiguous k-shape). Until then it holds a placeholder 1-mer shape.
    G.shape = gsa_key_is_plain(pattern) ? Shape::parse("#") : Shape::parse(pattern);
    G.n = (size_t)n;
    G.alphabet_size = (int)alpha;
    G.text.resize((size_t)n);
    if (n && !detail::gsa_get(i, &G.text[0], (size_t)n)) {
        { if (err) *err = path + ": truncated text";
          return false; }
    }
    if (!detail::gsa_get_vec(i, G.lex, (size_t)m)
        || !detail::gsa_get_vec(i, G.lex2orig, (size_t)m)
        || !detail::gsa_get_vec(i, G.sa, (size_t)m)
        || !detail::gsa_get_vec(i, G.lcp, (size_t)m)) {
        { if (err) *err = path + ": truncated arrays";
          return false; }
    }
    if (gcsa_text_hash(G.text) != th) {
        { if (err) *err = path + ": text hash mismatch -- the blob is corrupt";
          return false; }
    }
    return true;
}

} // namespace gcsa
