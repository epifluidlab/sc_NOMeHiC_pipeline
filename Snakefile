import yaml
import os

# Define config file
configfile: "configs/config.yaml"

# Set working directory
workdir: config["workdir"]

def get_fastq_prefixes(directory):
    prefixes = set()
    for filename in os.listdir(directory):
        if filename.endswith('.fastq.gz'):
        #if filename.endswith('.R2_fastqc.zip'):
            #prefix = filename.rsplit('.')[0]
            prefix = filename.rsplit('_', 2)[0]  # keep everything before the last two underscores
            prefixes.add(prefix)
    return list(prefixes)

def discover_pairs_from_trimmed(trimmed_dir, indices=None):
    """Scan `trimmed_dir` for files like {prefix}.{index}.R1_val_1.fq.gz and
    return the sorted list of (prefix, index) PAIRS that actually exist.
    The `index` is the suffix after the last `.` in the filename stem, so the
    prefix can itself contain dots. If `indices` is provided, only pairs whose
    `index` is in that set are returned. If `indices` is None, every well-
    formed `{prefix}.{index}.R1_val_1.fq.gz` file contributes a pair regardless
    of what its `index` is — useful when starting from already-demultiplexed
    trimmed fastqs and you don't want to maintain a barcode whitelist."""
    if not os.path.isdir(trimmed_dir):
        return []
    indices_set = set(indices) if indices is not None else None
    pairs = set()
    for fn in os.listdir(trimmed_dir):
        if not fn.endswith(".R1_val_1.fq.gz"):
            continue
        stem = fn[:-len(".R1_val_1.fq.gz")]
        if "." not in stem:
            continue
        prefix, _, idx = stem.rpartition(".")
        if indices_set is None or idx in indices_set:
            pairs.add((prefix, idx))
    return sorted(pairs)

# Pipeline starting mode — see configs/config.yaml for the description.
# fileindex is REQUIRED in raw mode (drives demultiplexing). In trimmed mode
# it's OPTIONAL — if not set, INDICES is auto-derived from the trimmed dir
# (every filename's suffix after the last `.` before `.R1_val_1.fq.gz`).
START_FROM = config.get("start_from", "raw")
_HAS_FILEINDEX = bool(config.get("fileindex"))
if START_FROM == "raw" and not _HAS_FILEINDEX:
    raise ValueError("start_from='raw' requires config['fileindex'] for demultiplexing.")
if _HAS_FILEINDEX:
    with open(config["fileindex"]) as f:
        INDICES = [line.strip() for line in f if line.strip()]
else:
    INDICES = []  # populated below from the trimmed dir
print("INDICES (from fileindex):", INDICES)

# Trimmed-fastq directory. Both the (raw-mode) demultiplex_fastqc_trim rule's
# trimmed outputs and the mapping rule's inputs live here, so override this
# when your already-trimmed fastqs are in a non-default folder.
TRIMMED_DATA = config.get("trimmed_data", "03.trimmed_fastq_snakemake")

if START_FROM == "raw":
    FASTQ_PREFIXES = get_fastq_prefixes(config["data"])
    CELLS = [(p, i) for p in FASTQ_PREFIXES for i in INDICES]
    SKIP_DEMUX_TRIM = False
elif START_FROM == "trimmed":
    # Three ways to determine which cells to process, in order of precedence:
    #   1. config['fastq_prefixes']      — inline list of prefixes. Cross with
    #      INDICES if fileindex is set; otherwise pair with whatever indices
    #      each listed prefix actually has on disk.
    #   2. config['fastq_prefixes_file'] — same semantics as (1) but read from
    #      a one-per-line text file.
    #   3. (default) auto-discover from TRIMMED_DATA — every (prefix, index)
    #      pair that actually exists. If fileindex is set, used as a barcode
    #      whitelist; otherwise every distinct suffix found becomes a valid
    #      index. This handles datasets where some prefixes only have a
    #      subset of barcodes, or where you mix samples from batches that
    #      used different barcode sets entirely.
    _explicit_prefixes = None
    if config.get("fastq_prefixes"):
        _explicit_prefixes = list(config["fastq_prefixes"])
    elif config.get("fastq_prefixes_file"):
        with open(config["fastq_prefixes_file"]) as f:
            _explicit_prefixes = [line.strip() for line in f if line.strip()]
    if _explicit_prefixes is not None:
        FASTQ_PREFIXES = _explicit_prefixes
        if INDICES:
            CELLS = [(p, i) for p in FASTQ_PREFIXES for i in INDICES]
        else:
            # No fileindex — pair the listed prefixes with the indices each
            # actually has on disk.
            disk_pairs = discover_pairs_from_trimmed(TRIMMED_DATA, indices=None)
            allowed = set(FASTQ_PREFIXES)
            CELLS = sorted(pi for pi in disk_pairs if pi[0] in allowed)
    else:
        # Pure auto-discover.
        CELLS = discover_pairs_from_trimmed(TRIMMED_DATA, INDICES or None)
        FASTQ_PREFIXES = sorted({p for p, _ in CELLS})
    if not INDICES:
        # In trimmed mode without fileindex, derive INDICES from the chosen pairs.
        INDICES = sorted({i for _, i in CELLS})
        print("INDICES (auto-derived from", TRIMMED_DATA, "):", INDICES)
    if not CELLS:
        raise ValueError(
            f"start_from='trimmed' but no (prefix, index) pairs resolved. Set "
            f"config['fastq_prefixes'] (list), config['fastq_prefixes_file'] "
            f"(path to one-per-line file), or place "
            f"{{prefix}}.{{index}}.R1_val_1.fq.gz files under {TRIMMED_DATA}/."
        )
    SKIP_DEMUX_TRIM = True
else:
    raise ValueError(f"Unknown start_from mode: {START_FROM!r}. Use 'raw' or 'trimmed'.")

# Global wildcard constraints — pin {index} to the resolved values so snakemake
# never picks the wrong wildcard split for files like
# `batch1_sc1.ATCACG_sorted_by_name.calmd.bam`.
wildcard_constraints:
    index = "|".join(INDICES) if INDICES else r"[^./]+",

# Shell prefix prepended to every rule's shell. Default puts the scnomehic conda
# env's bin FIRST in PATH (so `samtools` resolves to the conda-shipped 1.21
# instead of the system samtools-1.16 which links against an absent
# libcrypto.so.10 on Quest compute nodes), and appends hic_env's bin (provides
# `hicluster` for the hiccluster rules). Override via config['shell_prefix']
# to point at a different conda layout, or set it to "" to disable.
_DEFAULT_SHELL_PREFIX = (
    "export PATH=/projects/b1198/epifluidlab/yoshii/software/conda/envs/scnomehic/bin:"
    "$PATH:/projects/b1198/epifluidlab/yoshii/software/conda/envs/hic_env/bin; "
)
shell.prefix(config.get("shell_prefix", _DEFAULT_SHELL_PREFIX))

include: "rules/hiccluster.smk"
#include: "rules/GCHnorm.smk"

def _trinuc_autosome():
    """Match Bis-QC.pl's autosome selection: chr21 if the reference has it,
    else chr19 (so mm10/mouse runs work). Read from the reference .fai when
    available; otherwise just default to chr21 for backward compat with hg38."""
    fai = config.get("reference", "") + ".fa.fai"
    try:
        with open(fai) as f:
            chroms = {line.split("\t", 1)[0] for line in f if line.strip()}
        if "chr21" in chroms:
            return "chr21"
        if "chr19" in chroms:
            return "chr19"
    except (FileNotFoundError, OSError):
        pass
    return "chr21"

TRINUC_AUTOSOME = _trinuc_autosome()

def _resolve_bed_or_gz(template):
    """Return a snakemake input function that, at DAG-build time, picks the
    bistools .bed.gz output by default, falling back to a legacy plain .bed
    only if that's what's actually on disk (e.g., for older runs that haven't
    been bulk-gzip-migrated yet). `template` may end in either `.bed` or
    `.bed.gz`; this function normalizes."""
    def _inner(wildcards):
        path = template.format(**dict(wildcards.items()))
        gz = path if path.endswith(".gz") else path + ".gz"
        plain = path[:-3] if path.endswith(".gz") else path
        # Legacy on-disk plain bed → use it (lets older runs still be consumed
        # without forcing a rerun of bistools)
        if os.path.exists(plain) and not os.path.exists(gz):
            return plain
        # Default: the .gz is either already on disk or will be produced by
        # rule bistools.
        return gz
    return _inner

def _expand_cells(template, **extra):
    """Expand `template` over every (prefix, index) pair in CELLS, plus any
    extra wildcards. Like snakemake's expand() but pair-aware so we don't
    construct invalid cross-product cells."""
    if not extra:
        return [template.format(prefix=p, index=i) for p, i in CELLS]
    keys = list(extra.keys())
    vals = [extra[k] for k in keys]
    out = []
    for p, i in CELLS:
        # iterate Cartesian over extra dims
        from itertools import product
        for combo in product(*vals):
            fmt = {"prefix": p, "index": i}
            fmt.update(dict(zip(keys, combo)))
            out.append(template.format(**fmt))
    return out

# Methylation per-bin output configuration. The pipeline produces one TSV
# per cell × context × bin_definition. Bin definitions can be either:
#   - fixed bin sizes in bp (config['methylation_bin_sizes'])
#   - custom BED files of regions (config['methylation_region_beds'] is a
#     mapping label → bed path; the label becomes the filename suffix)
METHYLATION_CONTEXTS = ["GCH", "HCG"]
METHYLATION_BIN_SIZES = list(config.get("methylation_bin_sizes", [1000, 50000, 100000, 250000]))
METHYLATION_REGION_LABELS = list((config.get("methylation_region_beds") or {}).keys())

def _bin_size_to_label(n):
    """Render a bp count as a human-readable filename label: 1Mb / 5kb / 500bp.
    Only uses Mb/kb when n divides the unit evenly; otherwise falls back to bp
    so the round-trip parse is unambiguous."""
    n = int(n)
    if n >= 1_000_000 and n % 1_000_000 == 0:
        return f"{n // 1_000_000}Mb"
    if n >= 1000 and n % 1000 == 0:
        return f"{n // 1000}kb"
    return f"{n}bp"

def _bin_label_to_size(label):
    """Inverse of _bin_size_to_label: parse a 1Mb / 5kb / 500bp label back to bp."""
    if label.endswith("Mb"):
        return int(label[:-2]) * 1_000_000
    if label.endswith("kb"):
        return int(label[:-2]) * 1000
    if label.endswith("bp"):
        return int(label[:-2])
    raise ValueError(f"Unknown bin-size label: {label!r} (expected NNbp/NNkb/NNMb)")

METHYLATION_BIN_LABELS = [_bin_size_to_label(n) for n in METHYLATION_BIN_SIZES]

# Build rule_all target list. When SKIP_DEMUX_TRIM is set, drop the
# 02.fastqc outputs (which are produced by demultiplex_fastqc_trim).
_rule_all_targets = []
if not SKIP_DEMUX_TRIM:
    _rule_all_targets += _expand_cells("02.fastqc_out_snakemake/{prefix}.{index}.R1_fastqc.html")
    _rule_all_targets += _expand_cells("02.fastqc_out_snakemake/{prefix}.{index}.R1_fastqc.zip")
    _rule_all_targets += _expand_cells("02.fastqc_out_snakemake/{prefix}.{index}.R2_fastqc.html")
    _rule_all_targets += _expand_cells("02.fastqc_out_snakemake/{prefix}.{index}.R2_fastqc.zip")
# always — downstream targets, applicable in both modes
_rule_all_targets += _expand_cells("04.alignment_snakemake/{prefix}.{index}.summary.txt.gz")
_rule_all_targets += _expand_cells("07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.hist.txt.gz")
# hicluster intentionally NOT in default `rule all` — it runs only when the
# user provides `selected_cells_file` in config and invokes the
# `hicluster_selected` rule explicitly. See rules/hiccluster.smk.
_rule_all_targets += _expand_cells("07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.{context}.6plus2.bed.gz",
                                   context=METHYLATION_CONTEXTS)
_rule_all_targets += _expand_cells("08.methylation_snakemake/{prefix}.{index}.{context}.{label}_methylation.txt.gz",
                                   context=METHYLATION_CONTEXTS, label=METHYLATION_BIN_LABELS)
if METHYLATION_REGION_LABELS:
    _rule_all_targets += _expand_cells("08.methylation_snakemake/{prefix}.{index}.{context}.{label}_methylation.txt.gz",
                                       context=METHYLATION_CONTEXTS, label=METHYLATION_REGION_LABELS)

rule all:
    input: _rule_all_targets

rule demultiplex_fastqc_trim: # combination step to trigger re-runs if failed
    input:
        # Placeholders if data/fileindex are absent — the rule is unused in
        # trimmed mode but its definition still has to parse.
        r1 = os.path.join(config.get("data") or ".", "{prefix}_R1_001.fastq.gz"),
        r2 = os.path.join(config.get("data") or ".", "{prefix}_R2_001.fastq.gz"),
        index_file = config.get("fileindex") or "."
    output:
        # demultiplex output
        r1_demul = "01.demul_fastq_snakemake/{prefix}.{index}.R1.fastq.gz",
        r2_demul = "01.demul_fastq_snakemake/{prefix}.{index}.R2.fastq.gz",

        # first fastqc output
        r1_out_html = "02.fastqc_out_snakemake/{prefix}.{index}.R1_fastqc.html",
        r1_out_zip = "02.fastqc_out_snakemake/{prefix}.{index}.R1_fastqc.zip",
        r2_out_html = "02.fastqc_out_snakemake/{prefix}.{index}.R2_fastqc.html",
        r2_out_zip = "02.fastqc_out_snakemake/{prefix}.{index}.R2_fastqc.zip",

        # trimming and second fastqc output
        r1_out_v1_html = "02.fastqc_out_snakemake/{prefix}.{index}.R1_val_1_fastqc.html",
        r1_out_v1_zip = "02.fastqc_out_snakemake/{prefix}.{index}.R1_val_1_fastqc.zip",
        r1_out_v2_html = "02.fastqc_out_snakemake/{prefix}.{index}.R2_val_2_fastqc.html",
        r1_out_v2_zip = "02.fastqc_out_snakemake/{prefix}.{index}.R2_val_2_fastqc.zip",

        r1_out_trimmed = f"{TRIMMED_DATA}/{{prefix}}.{{index}}.R1_val_1.fq.gz",
        r1_out_report  = f"{TRIMMED_DATA}/{{prefix}}.{{index}}.R1.fastq.gz_trimming_report.txt",
        r2_out_trimmed = f"{TRIMMED_DATA}/{{prefix}}.{{index}}.R2_val_2.fq.gz",
        r2_out_report  = f"{TRIMMED_DATA}/{{prefix}}.{{index}}.R2.fastq.gz_trimming_report.txt"
    threads: 8
    resources:
        mem_mb=32000
    benchmark: "benchmarks/00.demultiplex/demul.{prefix}.{index}_benchmark.txt"
    params:
        out_prefix = "01.demul_fastq_snakemake/{prefix}",
        outdir = "02.fastqc_out_snakemake",
        bisulfitehic = config["bisulfitehic"],

        # trimming params
        outdir_fq = TRIMMED_DATA,
        clip_r1 = config["clip_r1"],
        clip_r2 = config["clip_r2"],
        three_prime_clip_r1 = config["three_prime_clip_r1"],
        three_prime_clip_r2 = config["three_prime_clip_r2"]
    log:
        "logs/00.demultiplex_snakemake/demultiplex_{prefix}_{index}.log"
    shell: # temporary sleep and touch() solution
        """
        echo $(pwd)
        module load java
        perl {params.bisulfitehic}/src/perl/demultiplex_ecker_scWGBS.pl \
        {params.out_prefix} {input.index_file} {input.r1} {input.r2} \
        2> {log}
        sleep 2

        fastqc --outdir {params.outdir} -t {threads} {output.r1_demul} {output.r2_demul} \
        2> {log}
        echo "First fastqc done"
        sleep 2

        trim_galore --paired_end \
        --clip_R1 {params.clip_r1} \
        --clip_R2 {params.clip_r2} \
        --three_prime_clip_R1 {params.three_prime_clip_r1} \
        --three_prime_clip_R2 {params.three_prime_clip_r2} \
        -j {threads} \
        -o {params.outdir_fq} --gzip --fastqc --fastqc_args "--outdir {params.outdir} -t 32" {output.r1_demul} {output.r2_demul} \
        2> {log}
        echo "Trimming done"
        """

rule mapping:
    input:
        r1_trimmed = f"{TRIMMED_DATA}/{{prefix}}.{{index}}.R1_val_1.fq.gz",
        r2_trimmed = f"{TRIMMED_DATA}/{{prefix}}.{{index}}.R2_val_2.fq.gz"
    output:
        bam = "04.alignment_snakemake/{prefix}.{index}.bam"
    threads: 8
    resources:
        mem_mb=16000
    benchmark: "benchmarks/02.mapping/alignment.{prefix}.{index}_benchmark.txt"
    params:
        reference = config["reference"],
        picard = config["picard"],
        bisulfitehic = config["bisulfitehic"],
        # Bhmem -enzymeList. If config['enzyme_list'] is set use that;
        # else fall back to the legacy {reference}.DpnII.span_region.bedgraph
        # convention (works for the hg38 setup, fails for mm10 where the file
        # is named differently).
        enzyme_list = config.get("enzyme_list") or f'{config["reference"]}.DpnII.span_region.bedgraph'
    log:
        "logs/02.alignment_snakemake/realignment.{prefix}.{index}.log"
    shell: # how to load jdk here without messing up -> ok it works, wondering if this is best way to go about this
        """
        module load java/jdk-17.0.2+8
        picard={params.picard}

        java -Xmx15G -Djava.library.path={params.bisulfitehic}/jbwa/jbwa-1.0.0/src/main/native \
        -cp "{params.bisulfitehic}/target/bisulfitehic-0.38-jar-with-dependencies.jar:{params.bisulfitehic}/jbwa/jbwa-1.0.0/jbwa.jar" \
        main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem {params.reference}.fa \
        {output.bam} {input.r1_trimmed} {input.r2_trimmed} \
        -t {threads} -rgId {wildcards.prefix}.{wildcards.index} -rgSm scNOMeHiC -nonDirectional -pbat -buffer 1000 -enzymeList {params.enzyme_list} -outputMateDiffChr \
        > {log} 2>&1
        """

rule bamprocess:
    input:
        bam = "04.alignment_snakemake/{prefix}.{index}.bam"
    output:
        # kept output
        calmd_bam = "04.alignment_snakemake/{prefix}.{index}.calmd.bam",
        calmd_bam_bai = "04.alignment_snakemake/{prefix}.{index}.calmd.bam.bai",
        sorted_bam = "04.alignment_snakemake/{prefix}.{index}_sorted_by_name.calmd.bam",
        filtered_bam = "04.alignment_snakemake/{prefix}.{index}.calmd.sorted_by_name.filtered.bam",
        bed = "04.alignment_snakemake/{prefix}.{index}.bed.gz"

        # temporary output (maybe could optimize to use temp() here instead of temp in the shell command?)
    threads: 16
    resources:
        mem_mb=64000
    benchmark: "benchmarks/03.bamprocess/bamprocess.{prefix}.{index}_benchmark.txt"
    params:
        reference = config["reference"],
        picard = config["picard"],
        shortcut = "04.alignment_snakemake/{prefix}.{index}"
    log:
        "logs/03.bamprocess_snakemake/bamprocess.{prefix}.{index}.log"
    shell: # if doesn't work in a force rerun, remove threads and maybe 2>/dev/null
        """
        module load java/jdk-17.0.2+8
	picard={params.picard}

        samtools sort --threads {threads} -T {params.shortcut}.tmp -n {input.bam} | \
        samtools fixmate -m --threads {threads} - - | \
        samtools sort --threads {threads} -T {params.shortcut}.cor - | \
        samtools markdup -T {params.shortcut}.mdups --threads {threads} - - | \
        samtools calmd --threads {threads} -b - {params.reference}.fa 2>/dev/null > {output.calmd_bam} 

        samtools index -@ {threads} {output.calmd_bam}

        samtools sort -@ {threads} -n -o {output.sorted_bam} {output.calmd_bam}
        samtools view -@ {threads} -f 2 -F 780 -q 30 -b {output.sorted_bam} > {output.filtered_bam}

        bedtools bamtobed -bedpe -mate1 -i {output.filtered_bam} | gzip -nc > {output.bed} 2> {log}
        """

rule qc:
    input:
        sorted_bam = "04.alignment_snakemake/{prefix}.{index}_sorted_by_name.calmd.bam"
    output:
        summary = "04.alignment_snakemake/{prefix}.{index}.summary.txt.gz"
    threads: 8
    resources:
        mem_mb=16000
    benchmark: "benchmarks/qc/qc.{prefix}.{index}_benchmark.txt"
    params:
        bisulfitehic = config["bisulfitehic"],
        summary_tmp = "04.alignment_snakemake/{prefix}.{index}.summary.txt"
    log:
        "logs/04.qc_snakemake/qc.{prefix}.{index}.log"
    shell: # mh_reads_summary.v2.py writes plain text; gzip after to honor .gz output
        """
        python {params.bisulfitehic}/src/python/mh_reads_summary.v2.py \
        --in_cram {input.sorted_bam} --out_summary {params.summary_tmp} \
        2> {log}
        gzip -f {params.summary_tmp}
        """

rule bisqc:
    input:
        calmd_bam = "04.alignment_snakemake/{prefix}.{index}.calmd.bam"
    output:
        # qc output — .txt files post-gzipped after Bis-QC.pl finishes
        qc_hist = "07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.hist.txt.gz",
        qc_autosome = f"07.bistools_snakemake/qc/{{prefix}}.{{index}}/{{prefix}}.{{index}}.calmd.trinuc_methy.{TRINUC_AUTOSOME}.txt.gz",
        qc_chrM = "07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.trinuc_methy.chrM.txt.gz",
        qc_conv_plot = "07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.WCH.bisuflite_conv_distribution_plot.pdf",
        qc_bias_plot = "07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.WCH.methy_bias_plot.pdf",
        qc_cycle = "07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.WCH.methy.cycle.txt.gz"
    threads: 1
    resources:
        mem_mb=10000
    params:
        reference = config["reference"],
        vcf = config["variants"],
        indir = "04.alignment_snakemake/{prefix}.{index}",
        qc_outdir = "07.bistools_snakemake/qc/{prefix}.{index}/",
        methylation_outdir = "07.bistools_snakemake/methylation/{prefix}.{index}/",
        WCG_outdir = "07.bistools_snakemake/WCG/{prefix}.{index}/",
        picard = config["picard"],
        bistools = config["bistools"]
    log:
        "logs/07.bisqc_snakemake/bisqc.{prefix}.{index}"
    shell:
        # BisSNP-1.0.1 / GATK 3.8 require Java 8 for reflective walker discovery
        # (QuickMethylationLevel etc are not findable under Java 17).
        # BISTOOLS env var must be EXPORTED so methylation_bias_plot.pl (which
        # does `my $bistools_path = \`echo $BISTOOLS\`;`) can read it.
        """
        export picard={params.picard}
        module load java/jdk1.8.0_191
        module load libpng
        export BISTOOLS={params.bistools}

        perl {params.bistools}/Bis-QC/Bis-QC.pl \
        --QC_mode 1 \
        --disable_enzyme_eff_check \
        --disable_coverage_check \
        --pattern WCH \
        --nt {threads} \
        --mem 10 \
        --genome {params.reference}.fa \
        --dbsnp {params.vcf} \
        --bistools_path {params.bistools} {input.calmd_bam} \
        > {log} 2>&1
        cp {params.indir}*.txt {params.qc_outdir}
        cp {params.indir}*.pdf {params.qc_outdir}
        gzip -f {params.qc_outdir}*.txt
        """

rule bistools:
    input:
        calmd_bam = "04.alignment_snakemake/{prefix}.{index}.calmd.bam"
    output:
        # raw VCFs left uncompressed (downstream tools expect plain .vcf);
        # per-read tables + per-base BEDs are gzipped post-hoc (text, ~5x smaller).
        cyt_vcf       = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.vcf",
        snp_vcf       = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.snp.vcf",
        gch_per_read  = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.GCH.txt.gz",
        hcg_per_read  = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.HCG.txt.gz",
        # post-processed filtered + sorted VCF and per-base BEDs (consumed downstream)
        cyt_filt_vcf       = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.vcf",
        cyt_filt_summary   = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.vcf.cpgSummary.txt.gz",
        gch_6plus2         = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.GCH.6plus2.bed.gz",
        hcg_6plus2         = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.HCG.6plus2.bed.gz",
        gch_strand_6plus2  = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.GCH.strand.6plus2.bed.gz",
        hcg_strand_6plus2  = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.HCG.strand.6plus2.bed.gz"
    threads: 1
    resources:
        mem_mb=30000
    params:
        reference  = config["reference"],
        vcf        = config["variants"],
        bistools   = config["bistools"],
        bissnp_jar = config["bissnp_jar"],
        outdir     = "07.bistools_snakemake/methylation/{prefix}.{index}",
        sample     = "{prefix}.{index}"
    log:
        "logs/07.bistools_snakemake/bistools.{prefix}.{index}.log"
    shell:
        # bissnp_nomehic_usage.pl writes outputs in cwd, named from --prefix.
        # cd into the output dir, then call with the BAM/ref/vcf as absolute paths.
        # After the tool runs, gzip the big text outputs (.bed and per-read .txt)
        # to save ~5x space. VCFs stay uncompressed for downstream compatibility.
        """
        mkdir -p {params.outdir}
        BAM_ABS=$(readlink -f {input.calmd_bam})
        LOG_ABS=$(readlink -f {log} 2>/dev/null || echo "$PWD/{log}")
        cd {params.outdir}
        perl {params.bistools}/Bis-SNP/bissnp_nomehic_usage.pl \
            --prefix {params.sample} --mem 20 \
            {params.bissnp_jar} \
            ${{BAM_ABS}} {params.reference}.fa {params.vcf} \
            > ${{LOG_ABS}} 2>&1
        gzip -f *.6plus2.bed *.strand.6plus2.bed *.GCH.txt *.HCG.txt *.cpgSummary.txt 2>/dev/null || true
        """

rule rerun_bistools_all:
    input:
        _expand_cells("07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.GCH.6plus2.bed.gz")
        + _expand_cells("07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.HCG.6plus2.bed.gz")

_METHYL_BED_TEMPLATE = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.{context}.6plus2.bed"

rule methylation:
    # Fixed-bin-size methylation: one TSV.GZ per cell × context × bin_size.
    # Input picks .bed.gz if it exists on disk, else .bed (plain).
    # {label} constrained to NNbp/NNkb/NNMb so this rule doesn't collide
    # with methylation_region (non-digit {label}).
    # calcmethylation.py uses pandas to_csv which auto-gzips when the output
    # path ends in .gz.
    input:
        bed = _resolve_bed_or_gz(_METHYL_BED_TEMPLATE)
    output:
        methylation = "08.methylation_snakemake/{prefix}.{index}.{context}.{label}_methylation.txt.gz"
    wildcard_constraints:
        context = "GCH|HCG",
        label = r"\d+(bp|kb|Mb)"
    threads: 1
    resources:
        mem_mb=8000
    params:
        scripts = config["scripts"],
        reference = config["reference"],
        bin_size_bp = lambda wildcards: _bin_label_to_size(wildcards.label)
    log:
        "logs/08.methylation_snakemake/methylation.{prefix}.{index}.{context}.{label}.log"
    shell:
        """
        python {params.scripts}/calcmethylation.py \
            --outfile {output.methylation} \
            --bed_path {input.bed} \
            --context {wildcards.context} \
            --bin_size {params.bin_size_bp} \
            --chrom_size_filepath {params.reference}.chrom.sizes \
            2> {log}
        """


rule methylation_region:
    # Custom-region methylation: one TSV.GZ per cell × context × region_label.
    # Input picks .bed.gz if present, else .bed. {label} constrained to
    # start with a non-digit so this rule doesn't collide with methylation.
    input:
        bed = _resolve_bed_or_gz(_METHYL_BED_TEMPLATE)
    output:
        methylation = "08.methylation_snakemake/{prefix}.{index}.{context}.{label}_methylation.txt.gz"
    wildcard_constraints:
        context = "GCH|HCG",
        label = r"[A-Za-z][A-Za-z0-9_\-]*"
    threads: 1
    resources:
        mem_mb=8000
    params:
        scripts = config["scripts"],
        region_bed = lambda wildcards: (config.get("methylation_region_beds") or {})[wildcards.label]
    log:
        "logs/08.methylation_snakemake/methylation.{prefix}.{index}.{context}.{label}.log"
    shell:
        """
        python {params.scripts}/calcmethylation.py \
            --outfile {output.methylation} \
            --bed_path {input.bed} \
            --context {wildcards.context} \
            --bin_bed {params.region_bed} \
            2> {log}
        """
