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

# Ordered pool of 16 real strains for the n=8/n=16 real-strain pangenome
# sweep below, each entry taking a prefix of this list. Order is hand-picked
# so the n=8 prefix alone already spans lab/commensal, EHEC, UPEC/ExPEC,
# EAEC and ETEC pathotypes rather than front-loading near-duplicates -- see
# the pathotype notes on each dataset in category C below.
REAL_ECOLI_STRAIN_POOL = [
    "@ecoli_k12", "@ecoli_sakai", "@ecoli_cft073", "@ecoli_hs",
    "@ecoli_042", "@ecoli_edl933", "@ecoli_e24377a", "@ecoli_uti89",
    "@ecoli_umn026", "@ecoli_apec_o1", "@ecoli_iai1", "@ecoli_536",
    "@ecoli_55989", "@ecoli_s88", "@ecoli_ed1a", "@ecoli_iai39",
]

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

    # 13 more real strains (verified live, Aug 2026, same way as the three
    # above) so the real-strain pangenome sweep below can go up to n=8/n=16
    # without ever resequencing the same isolate twice. Pathotypes are as
    # reported in Rasko et al. 2008 J. Bacteriol (EDL933/HS/E24377A/UTI89/
    # 536/042/APEC O1) and Touchon et al. 2009 PLoS Genet (IAI1/ED1a/S88/
    # UMN026/IAI39/55989) -- both are classic, widely-used E. coli diversity
    # panels, not an arbitrary grab of whatever NCBI returned first.
    dict(name="ecoli_hs", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000017765.1 (ASM1776v1), commensal strain HS",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/017/765/"
              "GCF_000017765.1_ASM1776v1/GCF_000017765.1_ASM1776v1_genomic.fna.gz"),
         sha256="d427a238a5a5a32e4cb504ab8e11169c9c2486b56905d311fe8e698e4b492b8a"),
    dict(name="ecoli_042", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000027125.1 (ASM2712v1), EAEC strain 042",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/125/"
              "GCF_000027125.1_ASM2712v1/GCF_000027125.1_ASM2712v1_genomic.fna.gz"),
         sha256="6b449c0f341ccc448d1de7ba2e61d449cf66d0c336257e4caa673d46131498c5"),
    dict(name="ecoli_edl933", category="real_genome", kind="fetched",
         source=("NCBI RefSeq GCF_000006665.1 (ASM666v1), EHEC O157:H7 str. "
                 "EDL933 -- a second, independently isolated O157:H7 lineage "
                 "alongside ecoli_sakai above"),
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/665/"
              "GCF_000006665.1_ASM666v1/GCF_000006665.1_ASM666v1_genomic.fna.gz"),
         sha256="84ae40b897e1c7651612a3f62c8d5050e604608d8e573d07c193f9c37462af23"),
    dict(name="ecoli_e24377a", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000017745.1 (ASM1774v1), ETEC strain E24377A",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/017/745/"
              "GCF_000017745.1_ASM1774v1/GCF_000017745.1_ASM1774v1_genomic.fna.gz"),
         sha256="919e37f0a296bac5d8e4b221cdc0f4514df725cca4db53d170b72a6afabfd15e"),
    dict(name="ecoli_uti89", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000013265.1 (ASM1326v1), UPEC/ExPEC strain UTI89",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/265/"
              "GCF_000013265.1_ASM1326v1/GCF_000013265.1_ASM1326v1_genomic.fna.gz"),
         sha256="81c3f1807af638bce3f95c15a8d87ee24184b02777d52fb05634e00b55429124"),
    dict(name="ecoli_umn026", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000026325.1 (ASM2632v2), ExPEC strain UMN026",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/325/"
              "GCF_000026325.1_ASM2632v2/GCF_000026325.1_ASM2632v2_genomic.fna.gz"),
         sha256="d40e06f1490c5777c379ad925caa6361e8588e9c65535437d46ca44b3bebb46b"),
    dict(name="ecoli_apec_o1", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000014845.1 (ASM1484v1), avian-pathogenic ExPEC strain APEC O1",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/014/845/"
              "GCF_000014845.1_ASM1484v1/GCF_000014845.1_ASM1484v1_genomic.fna.gz"),
         sha256="f0ab95865593beb469d9423c75a5ba49bb5f3f142f61c5d0820b9980a41903f2"),
    dict(name="ecoli_536", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000013305.1 (ASM1330v1), UPEC/ExPEC strain 536",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/305/"
              "GCF_000013305.1_ASM1330v1/GCF_000013305.1_ASM1330v1_genomic.fna.gz"),
         sha256="68319392e5ea6cb477d6a4908298e936c8f8f0d1cff85054d956c17891589534"),
    dict(name="ecoli_iai1", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000026265.1 (ASM2626v1), commensal strain IAI1",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/265/"
              "GCF_000026265.1_ASM2626v1/GCF_000026265.1_ASM2626v1_genomic.fna.gz"),
         sha256="c65a68318dc65b95141d6bdb0f7457d6b3c561ae7e8271407bc215317074a5ae"),
    dict(name="ecoli_55989", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000026245.1 (ASM2624v1), EAEC strain 55989",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/245/"
              "GCF_000026245.1_ASM2624v1/GCF_000026245.1_ASM2624v1_genomic.fna.gz"),
         sha256="bd737a35a5e4b8d731d43d5708270f60381ba226acc63ca58a1cca83ec2c371a"),
    dict(name="ecoli_s88", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000026285.1 (ASM2628v2), ExPEC strain S88 (meningitis-associated)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/285/"
              "GCF_000026285.1_ASM2628v2/GCF_000026285.1_ASM2628v2_genomic.fna.gz"),
         sha256="335dd2635bcb9cc98b29a2058f060201ee065cc49c2c877fa9b5134e4694f1f8"),
    dict(name="ecoli_ed1a", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000026305.1 (ASM2630v1), commensal strain ED1a",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/305/"
              "GCF_000026305.1_ASM2630v1/GCF_000026305.1_ASM2630v1_genomic.fna.gz"),
         sha256="cab7fda0166436e7426dbe5f4341b07b15b517a9f7cc7fc8d9bc782c47d3b29e"),
    dict(name="ecoli_iai39", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_000026345.1 (ASM2634v1), ExPEC strain IAI39",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/345/"
              "GCF_000026345.1_ASM2634v1/GCF_000026345.1_ASM2634v1_genomic.fna.gz"),
         sha256="a8384b63f8071022f01ebe2cd9ff331a04117564a8f0f04aeceba2aa4b554261"),

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
    # genuinely distinct, independently-sequenced E. coli genomes
    # concatenated, rather than one reference plus simulated point mutations.
    # kind="concat" just cats the resolved FASTAs of `refs` together --
    # no divergence knob, no seed, because there's nothing to simulate.
    dict(name="pangenome_ecoli_real_n3", category="pangenome", kind="concat",
         refs=["@ecoli_k12", "@ecoli_sakai", "@ecoli_cft073"]),

    # n=8 and n=16 prefixes of REAL_ECOLI_STRAIN_POOL above -- the real-data
    # counterpart of the simulated pangenome_ecoli_n{1,2,4,8} sweep, now that
    # 16 distinct real strains are available to draw from.
    *[dict(name=f"pangenome_ecoli_real_n{n}", category="pangenome", kind="concat",
           refs=REAL_ECOLI_STRAIN_POOL[:n])
      for n in (8, 16)],
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
