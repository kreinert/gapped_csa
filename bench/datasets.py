"""Single source of truth for what goes into the benchmark suite.

Each entry's `kind` decides how run_suite.py turns it into a local FASTA
path -- this is the answer to "we can't store all the datasets": nothing
here is committed. `bench/` holds only code; every byte of sequence data
either gets regenerated from a seed, downloaded on demand into --data-dir,
or is expected to already exist on disk outside the repo (matching how
gapped_csa/.gitignore already ignores *.fasta, and how data/ already sits
outside the git repo entirely as a sibling of gapped_csa/).

  "synthetic": regenerated fresh on every run via one of the simulate_*
      binaries (built by the top-level Makefile). Cheap -- it's a seeded
      RNG -- so it is never cached and never touches --data-dir.

  "fetched": downloaded once from `url` by fetch_data.py, decompressed if
      needed, cached as <data-dir>/<name>.fasta, and verified against
      `sha256` on every subsequent run so results stay reproducible without
      committing a multi-hundred-MB FASTA. Leave sha256=None until the first
      successful fetch, then copy the printed hash back in here.

  "provided": already sitting on disk outside the repo (the curated
      frequent-k-mer-region files -- see category D below). `path` is
      relative to --data-dir; run_suite.py just checks it exists and tells
      you how to get it if not.

A dataset's `args` list may reference another dataset's resolved path with
"@name" (used by the pangenome entries below to point simulate_pangenome at
the E. coli reference) -- run_suite.py substitutes it before invoking the
generator.

Where DATA_DIR/BIN_DIR/TMP_DIR actually point on your machine is
deliberately not decided here -- see config.py. Nothing in this file (or
anywhere else tracked by git) should hardcode a real filesystem path.
"""

DATASETS = [
    # --- A. random --------------------------------------------------------
    dict(name="random_1e5", category="random", kind="synthetic",
         generator="simulate_random", args=["-n", "100000", "--seed", "1"]),
    dict(name="random_1e6", category="random", kind="synthetic",
         generator="simulate_random", args=["-n", "1000000", "--seed", "1"]),
    dict(name="random_1e7", category="random", kind="synthetic",
         generator="simulate_random", args=["-n", "10000000", "--seed", "1"]),

    # --- B. repetitive (tunable %) -----------------------------------------
    dict(name="repeat_100pct_x50_y200", category="repetitive", kind="synthetic",
         generator="simulate_repeats", args=["-x", "50", "-y", "200", "--seed", "1"]),
    dict(name="repeat_50pct_x50_y200", category="repetitive", kind="synthetic",
         generator="simulate_repeats",
         args=["-x", "50", "-y", "200", "--repetitive-frac", "0.5", "--seed", "1"]),
    dict(name="repeat_10pct_x50_y200", category="repetitive", kind="synthetic",
         generator="simulate_repeats",
         args=["-x", "50", "-y", "200", "--repetitive-frac", "0.1", "--seed", "1"]),

    # --- C. real genomes -----------------------------------------------
    # URLs/accessions to fill in for real use -- see docs/ for the exact
    # commands once chosen. Left unpinned (sha256=None) deliberately: the
    # design doc flags these as needing a decision (which human chromosome,
    # whole-chromosome vs subregion) before fetch_data.py should run for
    # real. Placeholder shown for E. coli, which is small and unambiguous.
    dict(name="ecoli_k12", category="real_genome", kind="fetched",
         source="NCBI RefSeq NC_000913.3",
         url=("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
              "?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text"),
         sha256=None),

    # --- D. curated frequent-k-mer regions ----------------------------------
    # Already exist somewhere outside this repo; point --data-dir (see
    # config.py) at wherever they live to use them as-is, no fetch/generate
    # step needed.
    dict(name="cod_freq_kmers", category="freq_kmers", kind="provided",
         path="cod_freq_skmers.fasta"),
    dict(name="kestrel_freq_kmers", category="freq_kmers", kind="provided",
         path="kestrel_freq_skmers.fasta"),
    dict(name="human_freq_kmers", category="freq_kmers", kind="provided",
         path="human_freq_skmers.fasta"),

    # --- E. pangenome sweep -------------------------------------------------
    # Same reference, same divergence; strain count is the one variable the
    # email asks about (1, 2, 4, 8). "@ecoli_k12" is resolved to that
    # dataset's local path by run_suite.py before simulate_pangenome runs.
    *[dict(name=f"pangenome_ecoli_n{n}", category="pangenome", kind="synthetic",
           generator="simulate_pangenome",
           args=["-r", "@ecoli_k12", "-n", str(n), "--divergence", "0.01", "--seed", "1"])
      for n in (1, 2, 4, 8)],
]

# Shapes to sweep per dataset. Trimmed down from the weight-30 family in
# bench_repetition.cpp/compare_algos.cpp (kept in sync by hand for now --
# see the note in run_suite.py) plus one plain-kmer shape for comparability
# with ordinary k-mer indexes.
SHAPES = [
    "#" * 21,                        # k=21 contiguous
    "#" * 30,                        # 30cont
    "#" * 20 + "." * 4 + "#" * 10,    # 20+g4+10
    "#" * 10 + "." * 4 + "#" * 20,    # 10+g4+20
]

# Per your Aug-14 call: tree-dp / tree-dp2 / tree-dp3 are out of scope.
ALGOS = ["greedy", "dep-order", "tree-dp4", "pseudoforest-dp"]
MAX_ADDS = [8, 16]

# Above this input size, run_suite.py sets GCSA_SKIP_SELFTEST=1 (the
# brute-force check in main.cpp is O(n) extra work per config; fine at
# repeat/random-synthetic scale, not worth it once real genomes are in the
# mix). Correctness is still covered by bench_repetition.cpp/compare_algos.cpp
# and ilp_baseline on small inputs -- this suite is about compression/timing,
# not re-proving correctness on every run.
SKIP_SELFTEST_ABOVE_BYTES = 5_000_000
