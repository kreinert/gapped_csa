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

  "concat": resolves each entry in `refs` (also "@name" references) and
      concatenates their FASTAs, in order, into one file under --tmp-dir.
      Used for a *real* (not simulated-divergence) multi-strain pangenome
      entry -- three distinct real genomes glued together, no seed or
      divergence knob because there's nothing being simulated.

A dataset's `args` (synthetic) or `refs` (concat) list may reference another
dataset's resolved path with "@name" -- run_suite.py substitutes it before
invoking the generator or concatenating.

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
    # URLs verified live (Aug 2026) -- see docs/benchmark_design.md S3C for
    # how each was chosen and why NCBI's static FTP mirror + UCSC goldenPath
    # were used instead of Ensembl (Ensembl restructured its FTP layout this
    # month; the old species-name URLs are being retired, so anything copied
    # from an older tutorial is likely already broken). sha256 is still
    # unpinned -- run `./fetch_data.py --print-hash <name>` after the first
    # successful fetch and copy the value back in here.
    dict(name="ecoli_k12", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000005845.2 (ASM584v2), chromosome NC_000913.3",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/"
              "GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"),
         sha256="53bb6a51b6e92139ced1e38f74b7938781027c52200922ff03718c2237d23bb4"),
    dict(name="dmel_genome", category="real_genome", kind="fetched",
         source="UCSC dm6 (= Ensembl BDGP6.46 / GCA_000001215.4), whole genome, soft-masked",
         url="https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz",
         sha256="2c211a6789ebaee418ecf9df847daee6d60e14d9b3377a2e55678890a212d7a9"),
    dict(name="human_chr21", category="real_genome", kind="fetched",
         source="UCSC hg38 (= GRCh38), chromosome 21 only",
         url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz",
         sha256="35c71b68436d1a278ecb6a1e875af3ba4020738a028a7feac769a6d62790ae1f"),

    # Second and third real E. coli strains, for the "real strains" pangenome
    # entry below -- distinct pathotypes/phylogroups from K-12, so genuinely
    # different genomes rather than resequencing the same isolate.
    dict(name="ecoli_sakai", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000008865.2 (ASM886v2), O157:H7 str. Sakai",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/"
              "GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz"),
         sha256="71c2e5c364293c9ba36fc2c7acbcaa75cd6884295fe06260ba198826a8b1ddd3"),
    dict(name="ecoli_cft073", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_014262945.1 (ASM1426294v1), uropathogenic CFT073",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/262/945/"
              "GCF_014262945.1_ASM1426294v1/GCF_014262945.1_ASM1426294v1_genomic.fna.gz"),
         sha256="0aa2c239f4ec0f7bfd2dc3c6cef1ec15969c80be0393aefd1b0289c143240978"),

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

    # Real-strain counterpart (secondary/stretch goal per the design doc):
    # three genuinely distinct, independently-sequenced E. coli genomes
    # concatenated, rather than one reference plus simulated point mutations.
    # kind="concat" just cats the resolved FASTAs of `refs` together --
    # no divergence knob, no seed, because there's nothing to simulate.
    dict(name="pangenome_ecoli_real_n3", category="pangenome", kind="concat",
         refs=["@ecoli_k12", "@ecoli_sakai", "@ecoli_cft073"]),
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
