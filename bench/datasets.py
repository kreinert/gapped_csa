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
      entry -- several distinct real genomes glued together, no seed or
      divergence knob because there's nothing being simulated.

A dataset's `args` (synthetic) or `refs` (concat) list may reference another
dataset's resolved path with "@name" -- run_suite.py substitutes it before
invoking the generator or concatenating.

Where DATA_DIR/BIN_DIR/TMP_DIR actually point on your machine is
deliberately not decided here -- see config.py. Nothing in this file (or
anywhere else tracked by git) should hardcode a real filesystem path.
"""

# Ordered pool of 100 real strains for the n=1/10/100 real-strain
# pangenome sweep below, each entry taking a prefix of this list.
# Aug-28 update: this used to be 16 strains hand-ordered to span
# pathotypes/phylogroups (see git history) -- replaced wholesale with a
# 100-strain pool confirmed very similar to one another (Mete, Aug 28
# 2026), so pathotype-spanning order no longer applies. Order here is
# just the reference strain ("@ecoli_003") first, then the rest in the
# order Mete's source list gave them.
REAL_ECOLI_STRAIN_POOL = [
    "@ecoli_003", "@ecoli_001", "@ecoli_002", "@ecoli_004", "@ecoli_005",
    "@ecoli_006", "@ecoli_007", "@ecoli_008", "@ecoli_009", "@ecoli_010",
    "@ecoli_011", "@ecoli_012", "@ecoli_013", "@ecoli_014", "@ecoli_015",
    "@ecoli_016", "@ecoli_017", "@ecoli_018", "@ecoli_019", "@ecoli_020",
    "@ecoli_021", "@ecoli_022", "@ecoli_023", "@ecoli_024", "@ecoli_025",
    "@ecoli_026", "@ecoli_027", "@ecoli_028", "@ecoli_029", "@ecoli_030",
    "@ecoli_031", "@ecoli_032", "@ecoli_033", "@ecoli_034", "@ecoli_035",
    "@ecoli_036", "@ecoli_037", "@ecoli_038", "@ecoli_039", "@ecoli_040",
    "@ecoli_041", "@ecoli_042", "@ecoli_043", "@ecoli_044", "@ecoli_045",
    "@ecoli_046", "@ecoli_047", "@ecoli_048", "@ecoli_049", "@ecoli_050",
    "@ecoli_051", "@ecoli_052", "@ecoli_053", "@ecoli_054", "@ecoli_055",
    "@ecoli_056", "@ecoli_057", "@ecoli_058", "@ecoli_059", "@ecoli_060",
    "@ecoli_061", "@ecoli_062", "@ecoli_063", "@ecoli_064", "@ecoli_065",
    "@ecoli_066", "@ecoli_067", "@ecoli_068", "@ecoli_069", "@ecoli_070",
    "@ecoli_071", "@ecoli_072", "@ecoli_073", "@ecoli_074", "@ecoli_075",
    "@ecoli_076", "@ecoli_077", "@ecoli_078", "@ecoli_079", "@ecoli_080",
    "@ecoli_081", "@ecoli_082", "@ecoli_083", "@ecoli_084", "@ecoli_085",
    "@ecoli_086", "@ecoli_087", "@ecoli_088", "@ecoli_089", "@ecoli_090",
    "@ecoli_091", "@ecoli_092", "@ecoli_093", "@ecoli_094", "@ecoli_095",
    "@ecoli_096", "@ecoli_097", "@ecoli_098", "@ecoli_099", "@ecoli_100",
]

# Categories that exist only to be referenced by "@name" from other entries
# (the real_genome_pool strains feed REAL_ECOLI_STRAIN_POOL above, which
# feeds the pangenome_ecoli_real_n{1,10,100} concat entries below; ecoli_003
# is both the pool's first strain and a concat ingredient in its own right)
# -- not meant to be run through the full shape x algo x max_add grid on
# their own. run_suite.py excludes them from a default (no
# --only-category/--only-dataset) sweep; by_name resolution for "@name"
# refs is unaffected either way since it's built from the unfiltered
# DATASETS. Pass --only-category real_genome_pool (or --only-dataset
# ecoli_042 etc.) if you deliberately want to sweep one of these directly.
INGREDIENT_ONLY_CATEGORIES = {"real_genome", "real_genome_pool"}

DATASETS = [
    # --- A. random --------------------------------------------------------
    # A single 1MB uniformly-random sequence -- the "no repeat structure at
    # all" synthetic baseline.
    dict(name="random_1e6", category="random", kind="synthetic",
         generator="simulate_random", args=["-n", "1000000", "--seed", "1"]),

    # --- B. repetitive (tunable %, tunable repeat length) --------------------
    # Multiple copies of a fixed repeat unit, followed by a non-repetitive
    # random block. Two knobs are swept: repetitiveness (what fraction of
    # the output is the repeat block -- 100% / 50% / 10%, via
    # --repetitive-frac) and repeat length (-x = number of copies, -y =
    # length of each copy in bp). x*y is held at 10000bp across all three
    # (x, y) pairs below, so only the *shape* of the repeat varies (many
    # short copies vs. few long ones), not its total footprint.
    *[dict(name=f"repeat_{pct}pct_x{x}_y{y}", category="repetitive", kind="synthetic",
           generator="simulate_repeats",
           args=(["-x", str(x), "-y", str(y)]
                 + (["--repetitive-frac", str(frac)] if frac != 1.0 else [])
                 + ["--seed", "1"]))
      for x, y in ((50, 200), (100, 100), (200, 50))
      for pct, frac in ((100, 1.0), (50, 0.5), (10, 0.1))],

    # --- C. real genomes -----------------------------------------------
    # Aug-28 update: the 16 hand-picked, pathotype-diverse E. coli strains
    # previously here (K-12, Sakai, CFT073, +13 more spanning EHEC/UPEC/
    # EAEC/ETEC pathotypes) have been REMOVED and replaced wholesale with a
    # pool of 100 E. coli assemblies confirmed (by Mete, Aug 28 2026) to be
    # very similar to one another -- i.e. this batch is deliberately near-
    # duplicate rather than pathotype-diverse, which answers (differently
    # than expected) the open question in benchmark_design.md S7 about
    # whether strain selection should be phylogroup-balanced. All 100 are
    # `fetched`-kind, `sha256=None` -- this session's sandbox (both the
    # cloud container and the device_bash bridge to Mete's machine) cannot
    # reach ftp.ncbi.nlm.nih.gov at all (org egress allowlist -- confirmed
    # blocked, not just untested), so none of the 100 could actually be
    # downloaded or hashed here. Run `./fetch_data.py` on a machine with
    # normal internet access, then copy each printed sha256 back into the
    # matching entry below (same workflow the previous batch used).
    #
    # ecoli_003 (GCF_003697165.2, ASM369716v2) is the single "primary" real
    # E. coli genome -- the one entry that (a) sits in the main real_genome
    # category alongside dmel_genome/human_chr21 for the full algorithm x
    # shape x max_add grid, and (b) is the `-r` reference simulate_pangenome
    # uses for the simulated n=1/2/4/8 sweep below, taking over ecoli_k12's
    # old role in both places. Picked by accession number alone (a 2018-era,
    # low-numbered RefSeq accession vs. the 040000000-048000000 range
    # everything else in this batch was assigned in 2025-2026) as weak
    # evidence of an older/more-established assembly -- CONFIRMED after
    # fetching (Aug 28): 2 contigs (`NZ_CP033092.2` chromosome + one
    # plasmid, both "complete genome"/"complete sequence"), and it's
    # actually the E. coli type strain (DSM 30083 = JCM 1649 = ATCC 11775).
    # Good pick.
    dict(name="ecoli_003", category="real_genome", kind="fetched",
         source="NCBI RefSeq GCF_003697165.2 (ASM369716v2)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/697/165/"
              "GCF_003697165.2_ASM369716v2/GCF_003697165.2_ASM369716v2_genomic.fna.gz"),
         sha256="17124ac56df45b547706ddb0b2a56274d893022f42fe445afb5a08e19eb874b1"),

    # The other 99 strains in the 100-strain "very similar" pool -- category
    # real_genome_pool (NOT real_genome) so they feed REAL_ECOLI_STRAIN_POOL
    # below without each also multiplying the main algorithm x shape x
    # max_add grid by 99x on top of ecoli_003/dmel_genome/human_chr21.
    # Named ecoli_001..ecoli_100 by the row order of the URL list Mete
    # supplied (Aug 28 2026); ecoli_003 above is #003 of that same list,
    # pulled out and renumbered nowhere -- it just keeps its dataset name
    # and additionally plays the primary/reference role.
    dict(name="ecoli_001", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_046581325.1 (ASM4658132v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/046/581/325/"
              "GCF_046581325.1_ASM4658132v1/GCF_046581325.1_ASM4658132v1_genomic.fna.gz"),
         sha256="ad7365b41099582386a97c59b8a5b3d9114b257cfff59974c71e31c0fe626558"),
    dict(name="ecoli_002", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_040279995.2 (ASM4027999v2)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/279/995/"
              "GCF_040279995.2_ASM4027999v2/GCF_040279995.2_ASM4027999v2_genomic.fna.gz"),
         sha256="4765ba8b1abb0ee521aced530bc47682d0727741f9c6b3514ebca897fa7c965c"),
    dict(name="ecoli_004", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_046581375.1 (ASM4658137v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/046/581/375/"
              "GCF_046581375.1_ASM4658137v1/GCF_046581375.1_ASM4658137v1_genomic.fna.gz"),
         sha256="b6a9a555183553b7556521920b270875533a11c39c8304eae62dae7fffd5c148"),
    dict(name="ecoli_005", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_046603975.1 (ASM4660397v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/046/603/975/"
              "GCF_046603975.1_ASM4660397v1/GCF_046603975.1_ASM4660397v1_genomic.fna.gz"),
         sha256="26630dacc2f0936132463135c9b0afb703805f9d206bcf82ecddd3f265ad52ed"),
    dict(name="ecoli_006", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_000008865.2 (ASM886v2)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/"
              "GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz"),
         sha256="71c2e5c364293c9ba36fc2c7acbcaa75cd6884295fe06260ba198826a8b1ddd3"),
    dict(name="ecoli_007", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_046581385.1 (ASM4658138v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/046/581/385/"
              "GCF_046581385.1_ASM4658138v1/GCF_046581385.1_ASM4658138v1_genomic.fna.gz"),
         sha256="1fe8c90f82b74eb9fca5538baf46c3fd92e3efe95bba35faeaa1a63d9c788e59"),
    dict(name="ecoli_008", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_040272325.2 (ASM4027232v2)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/272/325/"
              "GCF_040272325.2_ASM4027232v2/GCF_040272325.2_ASM4027232v2_genomic.fna.gz"),
         sha256="18b4b7b3244f8b39914ade590c7c25c03f1915c2d174ad2cbeeaeddd0d151931"),
    dict(name="ecoli_009", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_047115195.1 (ASM4711519v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/115/195/"
              "GCF_047115195.1_ASM4711519v1/GCF_047115195.1_ASM4711519v1_genomic.fna.gz"),
         sha256="fffa6f2d2aa564424ae828b5b4deeca9477ac466bab0cf7978f371d65fdb281f"),
    dict(name="ecoli_010", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_047300965.1 (ASM4730096v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/300/965/"
              "GCF_047300965.1_ASM4730096v1/GCF_047300965.1_ASM4730096v1_genomic.fna.gz"),
         sha256="75c368f5ad26621728bc271f4741786aba9de165a5f5f725e03dc54019a3bdef"),
    dict(name="ecoli_011", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_046605345.1 (ASM4660534v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/046/605/345/"
              "GCF_046605345.1_ASM4660534v1/GCF_046605345.1_ASM4660534v1_genomic.fna.gz"),
         sha256="26dfaf253231bea8cd124bfb69eb5333066b2d96d13ac67ae0194930f8c8e77d"),
    dict(name="ecoli_012", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_047300905.1 (ASM4730090v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/300/905/"
              "GCF_047300905.1_ASM4730090v1/GCF_047300905.1_ASM4730090v1_genomic.fna.gz"),
         sha256="98b74de52bf472f22aa4e4f60b81f126180f02342bce94d73c6fe6eff6faec80"),
    dict(name="ecoli_013", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_046713585.1 (18AR0845)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/046/713/585/"
              "GCF_046713585.1_18AR0845/GCF_046713585.1_18AR0845_genomic.fna.gz"),
         sha256="0934c29476c3ed5610328d2725fb246a088da8e345c4228a769c4e3155eb1e32"),
    dict(name="ecoli_014", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_047456025.1 (ASM4745602v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/456/025/"
              "GCF_047456025.1_ASM4745602v1/GCF_047456025.1_ASM4745602v1_genomic.fna.gz"),
         sha256="be6adb12087ece0b1e01deae6c4a968767dbc6b53c5e39c3981ff085c9a1bb71"),
    dict(name="ecoli_015", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_047037185.1 (ASM4703718v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/037/185/"
              "GCF_047037185.1_ASM4703718v1/GCF_047037185.1_ASM4703718v1_genomic.fna.gz"),
         sha256="ff592b8cbc17c5eaeb999e1e8e95f193691e4513d04946ac1d0865dada1df566"),
    dict(name="ecoli_016", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_047551885.1 (ASM4755188v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/551/885/"
              "GCF_047551885.1_ASM4755188v1/GCF_047551885.1_ASM4755188v1_genomic.fna.gz"),
         sha256="f3f31fe812490c63a1e3090e0e371ab0e34b4ff8b1890a66492c2982ae9fad44"),
    dict(name="ecoli_017", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048020285.1 (ASM4802028v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/020/285/"
              "GCF_048020285.1_ASM4802028v1/GCF_048020285.1_ASM4802028v1_genomic.fna.gz"),
         sha256="0bdb9f1ca644840a7110812788ae15588d514a1d1dd20ce8d11231f860f52d80"),
    dict(name="ecoli_018", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048161525.1 (Br_54_kylling)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/161/525/"
              "GCF_048161525.1_Br_54_kylling/GCF_048161525.1_Br_54_kylling_genomic.fna.gz"),
         sha256="418d21dd4c768b75189f1aef18e2a1e84e011621926eadae2722f7b4e5829866"),
    dict(name="ecoli_019", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048020275.1 (ASM4802027v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/020/275/"
              "GCF_048020275.1_ASM4802027v1/GCF_048020275.1_ASM4802027v1_genomic.fna.gz"),
         sha256="bedcabd825f1e93eee37da903868c46f94b78363298ce59bf7b01213fc5b5cfd"),
    dict(name="ecoli_020", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048016545.1 (ASM4801654v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/016/545/"
              "GCF_048016545.1_ASM4801654v1/GCF_048016545.1_ASM4801654v1_genomic.fna.gz"),
         sha256="cfea459db42125a7b73b39de99404ddaa2298bdbb848678e5f7422b3b5c3bf89"),
    dict(name="ecoli_021", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048163035.1 (Br_178_kylling)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/163/035/"
              "GCF_048163035.1_Br_178_kylling/GCF_048163035.1_Br_178_kylling_genomic.fna.gz"),
         sha256="f0e1b92985f69715e103669c72ebd2c54550322516f0acf862de3f8167c7fa93"),
    dict(name="ecoli_022", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048163005.1 (Br_167_kylling)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/163/005/"
              "GCF_048163005.1_Br_167_kylling/GCF_048163005.1_Br_167_kylling_genomic.fna.gz"),
         sha256="5ece80455781b001ec6d85212259bf72a31574241f6b150f51d3af5dbc04f639"),
    dict(name="ecoli_023", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048160285.1 (ASM4816028v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/160/285/"
              "GCF_048160285.1_ASM4816028v1/GCF_048160285.1_ASM4816028v1_genomic.fna.gz"),
         sha256="0c5429cf10562ca1b1fded7165f7fd6a678adbe5e37cfd5e6598964bf5983d91"),
    dict(name="ecoli_024", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048162985.1 (Br_112_kylling)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/162/985/"
              "GCF_048162985.1_Br_112_kylling/GCF_048162985.1_Br_112_kylling_genomic.fna.gz"),
         sha256="e97ec009f190e67afa400f087ef509d47d5420796cecd71d93b37b93c81b5541"),
    dict(name="ecoli_025", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048163845.1 (A65EC)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/163/845/"
              "GCF_048163845.1_A65EC/GCF_048163845.1_A65EC_genomic.fna.gz"),
         sha256="0a099b07e0bfee90c758fab6b84b961b0c50a2443c5334b6708b30b266f61822"),
    dict(name="ecoli_026", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048163865.1 (A134_1EC)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/163/865/"
              "GCF_048163865.1_A134_1EC/GCF_048163865.1_A134_1EC_genomic.fna.gz"),
         sha256="4f0f9b2d6f749ab19b192e59fc44d6ca96a93886642babf88d42fc35d8398afa"),
    dict(name="ecoli_027", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048296705.1 (ASM4829670v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/296/705/"
              "GCF_048296705.1_ASM4829670v1/GCF_048296705.1_ASM4829670v1_genomic.fna.gz"),
         sha256="dd0fb8a8d5b993ed9a7bd25b79e4c69245418c499026b4c9e267d8102ec5e1e8"),
    dict(name="ecoli_028", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048401735.1 (ASM4840173v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/401/735/"
              "GCF_048401735.1_ASM4840173v1/GCF_048401735.1_ASM4840173v1_genomic.fna.gz"),
         sha256="bb14e42e12e443121929b4dc7337f9e66789f700cb655b081948faf0e3e0aca1"),
    dict(name="ecoli_029", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048402635.1 (ASM4840263v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/402/635/"
              "GCF_048402635.1_ASM4840263v1/GCF_048402635.1_ASM4840263v1_genomic.fna.gz"),
         sha256="d630aa0709005d7cb168725cff40be6a323c48e26f034b5f76a119748db54358"),
    dict(name="ecoli_030", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048403475.1 (ASM4840347v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/403/475/"
              "GCF_048403475.1_ASM4840347v1/GCF_048403475.1_ASM4840347v1_genomic.fna.gz"),
         sha256="d48b0691e7fbbbda0292754a648cacf6f69af6a0bd8d7213f8caa5cd03f6a059"),
    dict(name="ecoli_031", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048403405.1 (ASM4840340v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/403/405/"
              "GCF_048403405.1_ASM4840340v1/GCF_048403405.1_ASM4840340v1_genomic.fna.gz"),
         sha256="7cdbe60d8f679a09134bc7cb7860f4174b4371fd5577752a8154f4ed6d267ef8"),
    dict(name="ecoli_032", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048403485.1 (ASM4840348v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/403/485/"
              "GCF_048403485.1_ASM4840348v1/GCF_048403485.1_ASM4840348v1_genomic.fna.gz"),
         sha256="7b1bca962bf08bda0aa6e1b35b4265dd757cd4a6bb8edcb86e3f2210c4df8711"),
    dict(name="ecoli_033", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048457185.1 (ASM4845718v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/457/185/"
              "GCF_048457185.1_ASM4845718v1/GCF_048457185.1_ASM4845718v1_genomic.fna.gz"),
         sha256="c3c1146c11013292cf11d153ad5d4926f5bc3f4d83ec55574abbe875bfb725fd"),
    dict(name="ecoli_034", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048405095.1 (ASM4840509v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/405/095/"
              "GCF_048405095.1_ASM4840509v1/GCF_048405095.1_ASM4840509v1_genomic.fna.gz"),
         sha256="2abda9d277d9a4659d42fc5681223d2065528a0aebde265fa1a8b9af0dca228a"),
    dict(name="ecoli_035", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048405105.1 (ASM4840510v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/405/105/"
              "GCF_048405105.1_ASM4840510v1/GCF_048405105.1_ASM4840510v1_genomic.fna.gz"),
         sha256="a4cfc73b594b49d1dcd1f5ff240f3c9bcbe2cda67491d7f3112c50dfb24f100b"),
    dict(name="ecoli_036", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048457195.1 (ASM4845719v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/457/195/"
              "GCF_048457195.1_ASM4845719v1/GCF_048457195.1_ASM4845719v1_genomic.fna.gz"),
         sha256="580743a7daae969688c52c349c12ae595923f29ea09acd4ddb90c8bba2385bdb"),
    dict(name="ecoli_037", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048568955.1 (ASM4856895v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/568/955/"
              "GCF_048568955.1_ASM4856895v1/GCF_048568955.1_ASM4856895v1_genomic.fna.gz"),
         sha256="c0842147ba51e03c823a48a7d8c68e74f74225f1bee1f359c134ffa86281b0a9"),
    dict(name="ecoli_038", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048568935.1 (ASM4856893v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/568/935/"
              "GCF_048568935.1_ASM4856893v1/GCF_048568935.1_ASM4856893v1_genomic.fna.gz"),
         sha256="fa09ed682965c350aab0a5c4e135fdb1b987f7fe7f5528ddd29ed0675dcfcc52"),
    dict(name="ecoli_039", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048568965.1 (ASM4856896v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/568/965/"
              "GCF_048568965.1_ASM4856896v1/GCF_048568965.1_ASM4856896v1_genomic.fna.gz"),
         sha256="bf8aa2af0e5c160c786c4ce2fc1744e538d7b4dea1966c6bd295ba62cdd59b72"),
    dict(name="ecoli_040", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048568995.1 (ASM4856899v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/568/995/"
              "GCF_048568995.1_ASM4856899v1/GCF_048568995.1_ASM4856899v1_genomic.fna.gz"),
         sha256="83a037badbed363a7b5e003f31fbdaeb2d27979d8c92931d390ff55189710119"),
    dict(name="ecoli_041", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048569025.1 (ASM4856902v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/569/025/"
              "GCF_048569025.1_ASM4856902v1/GCF_048569025.1_ASM4856902v1_genomic.fna.gz"),
         sha256="d3b1d12737bf455d83135382d2cc9ee5924477613731543ed3f5fa41b8c2683a"),
    dict(name="ecoli_042", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048569015.1 (ASM4856901v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/569/015/"
              "GCF_048569015.1_ASM4856901v1/GCF_048569015.1_ASM4856901v1_genomic.fna.gz"),
         sha256="6b449c0f341ccc448d1de7ba2e61d449cf66d0c336257e4caa673d46131498c5"),
    dict(name="ecoli_043", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571005.1 (ASM4857100v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/005/"
              "GCF_048571005.1_ASM4857100v1/GCF_048571005.1_ASM4857100v1_genomic.fna.gz"),
         sha256="f89ae1338e2ad301296536ee6ce6b5a74e8b13cfe7ac9c4dbe8430fee61217c5"),
    dict(name="ecoli_044", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048569075.1 (ASM4856907v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/569/075/"
              "GCF_048569075.1_ASM4856907v1/GCF_048569075.1_ASM4856907v1_genomic.fna.gz"),
         sha256="37f771b40a796b765c49d22a53c32226cb44ae1c8c8be35037ff20cdda416699"),
    dict(name="ecoli_045", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571285.1 (ASM4857128v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/285/"
              "GCF_048571285.1_ASM4857128v1/GCF_048571285.1_ASM4857128v1_genomic.fna.gz"),
         sha256="7816d6dd3ae5cbaef6469c8d3a5f776bf80472114528e1bd9453b49bf99f5a81"),
    dict(name="ecoli_046", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571325.1 (ASM4857132v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/325/"
              "GCF_048571325.1_ASM4857132v1/GCF_048571325.1_ASM4857132v1_genomic.fna.gz"),
         sha256="d576219797d56d41ee168186551e29821d51c26c473b92bab4a97c41dd295a1d"),
    dict(name="ecoli_047", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048569065.1 (ASM4856906v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/569/065/"
              "GCF_048569065.1_ASM4856906v1/GCF_048569065.1_ASM4856906v1_genomic.fna.gz"),
         sha256="f70c10858e0d76631fc5cd99e5dbcb85b7ea1cc5e5bede899496c3f7896b9534"),
    dict(name="ecoli_048", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571335.1 (ASM4857133v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/335/"
              "GCF_048571335.1_ASM4857133v1/GCF_048571335.1_ASM4857133v1_genomic.fna.gz"),
         sha256="48ff3e98e364489e5eb502d67175c991a0572780df91835b11f3cf0511ea738f"),
    dict(name="ecoli_049", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571345.1 (ASM4857134v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/345/"
              "GCF_048571345.1_ASM4857134v1/GCF_048571345.1_ASM4857134v1_genomic.fna.gz"),
         sha256="42715923f5566f207e5a4d053c73785536355aaf5a53b1a9cd839ab8cce8dd23"),
    dict(name="ecoli_050", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571375.1 (ASM4857137v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/375/"
              "GCF_048571375.1_ASM4857137v1/GCF_048571375.1_ASM4857137v1_genomic.fna.gz"),
         sha256="fc26f0fc7896488aa043dfa6894679e8b495a04d148f2d823ed7182406a5bb8c"),
    dict(name="ecoli_051", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571365.1 (ASM4857136v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/365/"
              "GCF_048571365.1_ASM4857136v1/GCF_048571365.1_ASM4857136v1_genomic.fna.gz"),
         sha256="7c06a5bd911c64ba029795b23982f81b08f274b3405cbedc263e2f416c6e3972"),
    dict(name="ecoli_052", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571405.1 (ASM4857140v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/405/"
              "GCF_048571405.1_ASM4857140v1/GCF_048571405.1_ASM4857140v1_genomic.fna.gz"),
         sha256="1f37e5521db8a1dbdb0de24a1ca531e111a8eba763eccc2286422012d7c79a58"),
    dict(name="ecoli_053", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571355.1 (ASM4857135v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/355/"
              "GCF_048571355.1_ASM4857135v1/GCF_048571355.1_ASM4857135v1_genomic.fna.gz"),
         sha256="53a8a704a4ad98758716ccdd4fe09363dde67c8334547c362d167d5ef78855b6"),
    dict(name="ecoli_054", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571385.1 (ASM4857138v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/385/"
              "GCF_048571385.1_ASM4857138v1/GCF_048571385.1_ASM4857138v1_genomic.fna.gz"),
         sha256="a0692534778e904169feb9372e0c7b072f0229b91c02be96ad703a2126b73880"),
    dict(name="ecoli_055", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571395.1 (ASM4857139v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/395/"
              "GCF_048571395.1_ASM4857139v1/GCF_048571395.1_ASM4857139v1_genomic.fna.gz"),
         sha256="3ed3a678cdf84d81d2989227676af29d8a5ae98be125d9e1e78866ac1efe3d6b"),
    dict(name="ecoli_056", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571435.1 (ASM4857143v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/435/"
              "GCF_048571435.1_ASM4857143v1/GCF_048571435.1_ASM4857143v1_genomic.fna.gz"),
         sha256="0ff57486a488d0a8e943b46e72f2aabdc441fac55fb184bf69e4f2ac1ac89d6b"),
    dict(name="ecoli_057", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571445.1 (ASM4857144v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/445/"
              "GCF_048571445.1_ASM4857144v1/GCF_048571445.1_ASM4857144v1_genomic.fna.gz"),
         sha256="7cc33a25e4dab2898cdaf95c300f4272412adf00c7aacc72a05b3c93f0830c27"),
    dict(name="ecoli_058", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571465.1 (ASM4857146v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/465/"
              "GCF_048571465.1_ASM4857146v1/GCF_048571465.1_ASM4857146v1_genomic.fna.gz"),
         sha256="ee46dbe1b4883c55eddd56479c646220a9a0da0a7c7d2ced8bd7fdcb83dc1cf9"),
    dict(name="ecoli_059", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571535.1 (ASM4857153v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/535/"
              "GCF_048571535.1_ASM4857153v1/GCF_048571535.1_ASM4857153v1_genomic.fna.gz"),
         sha256="3742b93740f975dbf39dfde29df9b4f5f8fee9b98db71270979afc8c9ee24bef"),
    dict(name="ecoli_060", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571575.1 (ASM4857157v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/575/"
              "GCF_048571575.1_ASM4857157v1/GCF_048571575.1_ASM4857157v1_genomic.fna.gz"),
         sha256="207b74ee05f34673cf8eac00ad40ca9aaf6b38c88b2663727425585ddfd45243"),
    dict(name="ecoli_061", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571485.1 (ASM4857148v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/485/"
              "GCF_048571485.1_ASM4857148v1/GCF_048571485.1_ASM4857148v1_genomic.fna.gz"),
         sha256="df64732e1467a97732f47cb68f717fa2557fee644fb63cb4254cc00cb2715dc5"),
    dict(name="ecoli_062", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571455.1 (ASM4857145v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/455/"
              "GCF_048571455.1_ASM4857145v1/GCF_048571455.1_ASM4857145v1_genomic.fna.gz"),
         sha256="803f05eab0e30bc7121b14e104da141ccb458dacefd0c0d45988810c63950a67"),
    dict(name="ecoli_063", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571555.1 (ASM4857155v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/555/"
              "GCF_048571555.1_ASM4857155v1/GCF_048571555.1_ASM4857155v1_genomic.fna.gz"),
         sha256="758fb9bc5afff66cb7c31996df2434e7312e657ccc06ffbf55afa9723de48110"),
    dict(name="ecoli_064", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571585.1 (ASM4857158v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/585/"
              "GCF_048571585.1_ASM4857158v1/GCF_048571585.1_ASM4857158v1_genomic.fna.gz"),
         sha256="276a53f888e36a880b07b1ee328cae8f5245d90ee45c1c47393c893301f47691"),
    dict(name="ecoli_065", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048571595.1 (ASM4857159v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/571/595/"
              "GCF_048571595.1_ASM4857159v1/GCF_048571595.1_ASM4857159v1_genomic.fna.gz"),
         sha256="d57c52f4bc1da4d63701177b9e87d34c9e394616addb4b73554a699b14af28b7"),
    dict(name="ecoli_066", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048572565.1 (ASM4857256v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/572/565/"
              "GCF_048572565.1_ASM4857256v1/GCF_048572565.1_ASM4857256v1_genomic.fna.gz"),
         sha256="7684044ba7d85f76240528689856997f2d04e7209e46eace833798d12ca0347a"),
    dict(name="ecoli_067", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048572625.1 (ASM4857262v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/572/625/"
              "GCF_048572625.1_ASM4857262v1/GCF_048572625.1_ASM4857262v1_genomic.fna.gz"),
         sha256="121c3991354cc2f44f0ecb56f8cfb04b8ba8e795a531dc279785f9cedf23f23e"),
    dict(name="ecoli_068", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048572595.1 (ASM4857259v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/572/595/"
              "GCF_048572595.1_ASM4857259v1/GCF_048572595.1_ASM4857259v1_genomic.fna.gz"),
         sha256="e39988714c167bd6bfd7160821956d3a884a22826ec522835a259b1ebc913cbd"),
    dict(name="ecoli_069", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048572585.1 (ASM4857258v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/572/585/"
              "GCF_048572585.1_ASM4857258v1/GCF_048572585.1_ASM4857258v1_genomic.fna.gz"),
         sha256="5d285e9f084c683b9a596005ad816e71cd3c038ec84963d5da2555ceeaee7f9e"),
    dict(name="ecoli_070", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048572605.1 (ASM4857260v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/572/605/"
              "GCF_048572605.1_ASM4857260v1/GCF_048572605.1_ASM4857260v1_genomic.fna.gz"),
         sha256="8062e6de28c32569611648d25f8492232134da5ce4b0e2452838d65081b03d9b"),
    dict(name="ecoli_071", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048572615.1 (ASM4857261v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/572/615/"
              "GCF_048572615.1_ASM4857261v1/GCF_048572615.1_ASM4857261v1_genomic.fna.gz"),
         sha256="7aa63c3882a2df3e3ce3a5124c2b6bd24f0b32a06b81b691de8842ee7ca81bb6"),
    dict(name="ecoli_072", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048572645.1 (ASM4857264v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/572/645/"
              "GCF_048572645.1_ASM4857264v1/GCF_048572645.1_ASM4857264v1_genomic.fna.gz"),
         sha256="4edf5625d3046119471c3f5bc7d2723dc68be6126c20b08c423fedd73e07b791"),
    dict(name="ecoli_073", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048572655.1 (ASM4857265v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/572/655/"
              "GCF_048572655.1_ASM4857265v1/GCF_048572655.1_ASM4857265v1_genomic.fna.gz"),
         sha256="0d63ab2bb9a5048b6ac0fa346b2e9d7b4e0f40ff0f4b22932a417d3afeab2b8c"),
    dict(name="ecoli_074", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048572665.1 (ASM4857266v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/572/665/"
              "GCF_048572665.1_ASM4857266v1/GCF_048572665.1_ASM4857266v1_genomic.fna.gz"),
         sha256="b511b07affcce306abee2c301e3598aa0cf30c0d6833af34d792f5f919b9c0b7"),
    dict(name="ecoli_075", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048572705.1 (VTB96933v)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/572/705/"
              "GCF_048572705.1_VTB96933v/GCF_048572705.1_VTB96933v_genomic.fna.gz"),
         sha256="b77da74e8bd1f1bb880784e457310879ccf98159c5611898b404c8f569fefe27"),
    dict(name="ecoli_076", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048572675.1 (ASM4857267v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/572/675/"
              "GCF_048572675.1_ASM4857267v1/GCF_048572675.1_ASM4857267v1_genomic.fna.gz"),
         sha256="3674fe6d31715a16a15224a16f84bb83ef46a503d593fa149eb605a0dc045a51"),
    dict(name="ecoli_077", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585475.1 (ASM4858547v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/475/"
              "GCF_048585475.1_ASM4858547v1/GCF_048585475.1_ASM4858547v1_genomic.fna.gz"),
         sha256="faee44d58e669788f752f03e30444144b504bc2c708a60a7796fc17e5498a05a"),
    dict(name="ecoli_078", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585465.1 (ASM4858546v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/465/"
              "GCF_048585465.1_ASM4858546v1/GCF_048585465.1_ASM4858546v1_genomic.fna.gz"),
         sha256="c7b953e803859b1b25a59ddeaf2f66fecd7d34060c954b849f5a74b90a19d89b"),
    dict(name="ecoli_079", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048572685.1 (ASM4857268v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/572/685/"
              "GCF_048572685.1_ASM4857268v1/GCF_048572685.1_ASM4857268v1_genomic.fna.gz"),
         sha256="3cfcdd8d7176052ba78c2834a83c504bec16c0ced3c80d327eec3068c1d42f7c"),
    dict(name="ecoli_080", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585485.1 (ASM4858548v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/485/"
              "GCF_048585485.1_ASM4858548v1/GCF_048585485.1_ASM4858548v1_genomic.fna.gz"),
         sha256="7aea96151a4a29b9aebf77c2eae8b93d620429cba31dc70b025bb7dc43760f0d"),
    dict(name="ecoli_081", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585495.1 (ASM4858549v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/495/"
              "GCF_048585495.1_ASM4858549v1/GCF_048585495.1_ASM4858549v1_genomic.fna.gz"),
         sha256="84242d44535ef4615f977dda8113c86cbbcd637d05f426cc4fecd702bed1bd3b"),
    dict(name="ecoli_082", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585505.1 (ASM4858550v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/505/"
              "GCF_048585505.1_ASM4858550v1/GCF_048585505.1_ASM4858550v1_genomic.fna.gz"),
         sha256="1af2370f252fb50cb220824d9662fd0e15d817676f8099b5bc72a7f14c96bfa6"),
    dict(name="ecoli_083", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585515.1 (ASM4858551v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/515/"
              "GCF_048585515.1_ASM4858551v1/GCF_048585515.1_ASM4858551v1_genomic.fna.gz"),
         sha256="7bfbe55a2e4b412577c01d3ed21514b3a7d1d2a855190838791a084937389288"),
    dict(name="ecoli_084", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585525.1 (ASM4858552v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/525/"
              "GCF_048585525.1_ASM4858552v1/GCF_048585525.1_ASM4858552v1_genomic.fna.gz"),
         sha256="29c84e51f61dfaf8dbb2dbe780fcdd8f284dd86c4e08b29e1fd3fae0a230a652"),
    dict(name="ecoli_085", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585535.1 (ASM4858553v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/535/"
              "GCF_048585535.1_ASM4858553v1/GCF_048585535.1_ASM4858553v1_genomic.fna.gz"),
         sha256="acfa210645d9b14dbdf029baef2263b27c398640355624e0d86824e1d2421206"),
    dict(name="ecoli_086", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585545.1 (ASM4858554v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/545/"
              "GCF_048585545.1_ASM4858554v1/GCF_048585545.1_ASM4858554v1_genomic.fna.gz"),
         sha256="9fb6cff5ce1aa4b79f217cd1f0cb42ace5a1f63c1aa74a1281b81883b1735679"),
    dict(name="ecoli_087", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585555.1 (ASM4858555v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/555/"
              "GCF_048585555.1_ASM4858555v1/GCF_048585555.1_ASM4858555v1_genomic.fna.gz"),
         sha256="cc75abce743b75cf9b863edb1ddc27761fdbae719fbd6b80551fc8e8870c99e7"),
    dict(name="ecoli_088", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585565.1 (ASM4858556v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/565/"
              "GCF_048585565.1_ASM4858556v1/GCF_048585565.1_ASM4858556v1_genomic.fna.gz"),
         sha256="842fe3a8837dc73af6e94a7ef671f29fc92290a2c3a3b79b409e606549a7e2de"),
    dict(name="ecoli_089", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585575.1 (ASM4858557v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/575/"
              "GCF_048585575.1_ASM4858557v1/GCF_048585575.1_ASM4858557v1_genomic.fna.gz"),
         sha256="b3e0fced9ce573263bab21d6ff12426611cecc84b336229104fa2fa4a0984a2d"),
    dict(name="ecoli_090", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585585.1 (ASM4858558v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/585/"
              "GCF_048585585.1_ASM4858558v1/GCF_048585585.1_ASM4858558v1_genomic.fna.gz"),
         sha256="ab761cd124a2294e649f857b642f1001a861d87b49801c7ea487113025eac35b"),
    dict(name="ecoli_091", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585595.1 (ASM4858559v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/595/"
              "GCF_048585595.1_ASM4858559v1/GCF_048585595.1_ASM4858559v1_genomic.fna.gz"),
         sha256="8261a8195c8de294acaa5d0f866e6067f6b4455a3022803b0659cf80762e875e"),
    dict(name="ecoli_092", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048585605.1 (ASM4858560v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/585/605/"
              "GCF_048585605.1_ASM4858560v1/GCF_048585605.1_ASM4858560v1_genomic.fna.gz"),
         sha256="08bd6f1c9996f47bbd361aac4346aacee67b1fe4f077bd390ee066c2f2b17f16"),
    dict(name="ecoli_093", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_047782965.1 (ASM4778296v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/782/965/"
              "GCF_047782965.1_ASM4778296v1/GCF_047782965.1_ASM4778296v1_genomic.fna.gz"),
         sha256="9d1ced90656939a5d941e5857e075be8645fdfbb3333e076e84cc62538b5fa56"),
    dict(name="ecoli_094", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_047715505.1 (de_novo)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/715/505/"
              "GCF_047715505.1_de_novo/GCF_047715505.1_de_novo_genomic.fna.gz"),
         sha256="36d5bbd3381241b9ea5b42fa51687d66ea32c1a44873570fd8e617efb82eb56c"),
    dict(name="ecoli_095", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_047782975.1 (ASM4778297v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/782/975/"
              "GCF_047782975.1_ASM4778297v1/GCF_047782975.1_ASM4778297v1_genomic.fna.gz"),
         sha256="d49b19dff675e2c4155d9248d8f2e5d6191b04567927b9a7a55a3e3771154dd0"),
    dict(name="ecoli_096", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048011955.1 (ASM4801195v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/011/955/"
              "GCF_048011955.1_ASM4801195v1/GCF_048011955.1_ASM4801195v1_genomic.fna.gz"),
         sha256="77e479c7fe08ab25a8099b983a9b3eb1d1993634a916c33798d864b9960d7e32"),
    dict(name="ecoli_097", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_047782985.1 (ASM4778298v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/782/985/"
              "GCF_047782985.1_ASM4778298v1/GCF_047782985.1_ASM4778298v1_genomic.fna.gz"),
         sha256="e1a90f2726216f2fef428f358d815f7b01c6319c351713efeb6e2313532a9199"),
    dict(name="ecoli_098", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048012845.1 (ASM4801284v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/012/845/"
              "GCF_048012845.1_ASM4801284v1/GCF_048012845.1_ASM4801284v1_genomic.fna.gz"),
         sha256="fa9f951dea067c72546ac83e19a54f640bc672f6561f5cdebc7e1bfadbeb93a0"),
    dict(name="ecoli_099", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048011965.1 (ASM4801196v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/011/965/"
              "GCF_048011965.1_ASM4801196v1/GCF_048011965.1_ASM4801196v1_genomic.fna.gz"),
         sha256="ce574735c2274cfc363382317a1327dec81344c8fdf172df2edb8b24ac3d40fe"),
    dict(name="ecoli_100", category="real_genome_pool", kind="fetched",
         source="NCBI RefSeq GCF_048012795.1 (ASM4801279v1)",
         url=("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/012/795/"
              "GCF_048012795.1_ASM4801279v1/GCF_048012795.1_ASM4801279v1_genomic.fna.gz"),
         sha256="ce226c738424b8d389d6293948b948aada8d37871c2b283870a8c018d7ea2716"),

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

    # --- E. E. coli pangenome collections (1, 10, 100 genomes, ~95% identity)
    # Real, independently-sequenced strains from the 100-strain pool above
    # (REAL_ECOLI_STRAIN_POOL), concatenated -- not simulated divergence.
    # kind="concat" just cats the resolved FASTAs of `refs` together, in
    # pool order (ecoli_003 first, so every larger n's collection is a
    # superset of every smaller n's). Matches the paper's "(1) collections
    # of 1, 10, and 100 E. coli genomes sharing 95% identity" real dataset.
    *[dict(name=f"pangenome_ecoli_real_n{n}", category="pangenome", kind="concat",
           refs=REAL_ECOLI_STRAIN_POOL[:n])
      for n in (1, 10, 100)],
]

# Shapes to sweep per dataset. Trimmed down from the weight-30 family in
# bench_repetition.cpp/compare_algos.cpp (kept in sync by hand for now --
# see the note in run_suite.py) plus one plain-kmer shape for comparability
# with ordinary k-mer indexes.
SHAPES = [
    "#" * 21,                        # k=21 contiguous
    "#" * 30,                        # 30cont
    # "#" * 20 + "." * 4 + "#" * 10,    # 20+g4+10
    # "#" * 10 + "." * 4 + "#" * 20,    # 10+g4+20

    # Canonical spaced seeds from the literature (Sep-8; see the design doc
    # for the literature survey these were pulled from). Each comment gives
    # the seed's original 1/0 bit string and source.
    "###.#..#.#..##.###",              # PatternHunter optimal seed (Ma,
                                        # Tromp & Li, Bioinformatics 2002):
                                        # 111010010100110111, span18/w11,
                                        # optimized for ~70%-identity DNA
                                        # homology search.
    "###.#..##..#.#.####",             # BLASTZ "12of19" (Schwartz, Kent et
                                        # al., Genome Research 2003; still
                                        # LASTZ's default today):
                                        # 1110100110010101111, span19/w12.
    "##.##.##.##.##.#",                # NCBI discontiguous megablast
                                        # weight-11/span-16 "coding"
                                        # template (periodic "drop every
                                        # 3rd position", in contrast to the
                                        # two irregular optimized seeds
                                        # above): 1101101101101101.
    "###################...#.#.##.#.###..##.##",
                                        # Mete-supplied shape (Sep-8):
                                        # 11111111111111111110001010110101110011011,
                                        # span41/w31.
]

# Full retained set: greedy-size and greedy-degree (Phase I greedy variants),
# dep-order (dependency-order un-pin/retarget), and pseudoforest-dp /
# pseudoforest-dp-iterate (exact DP, single-fire and iterated). Every one of
# these now runs the same Phase II (dirty-set retarget) + LNS pipeline; see
# compress.hpp's top-of-file CompressAlgo comment.
ALGOS = ["greedy-size", "greedy-degree", "dep-order", "pseudoforest-dp", "pseudoforest-dp-iterate"]
MAX_ADDS = [8, 16, 32, 64, 128]

# Above this input size, run_suite.py sets GCSA_SKIP_SELFTEST=1 (the
# brute-force check in main.cpp is O(n) extra work per config; fine at
# repeat/random-synthetic scale, not worth it once real genomes are in the
# mix). Correctness is still covered by bench_repetition.cpp/compare_algos.cpp
# and ilp_baseline on small inputs -- this suite is about compression/timing,
# not re-proving correctness on every run.
SKIP_SELFTEST_ABOVE_BYTES = 5_000_000
