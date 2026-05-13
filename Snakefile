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

# Pipeline starting mode — see configs/config.yaml for the description.
START_FROM = config.get("start_from", "raw")
if START_FROM == "raw":
    FASTQ_PREFIXES = get_fastq_prefixes(config["data"])
    SKIP_DEMUX_TRIM = False
elif START_FROM == "trimmed":
    FASTQ_PREFIXES = list(config.get("fastq_prefixes") or [])
    if not FASTQ_PREFIXES:
        raise ValueError("start_from='trimmed' requires `fastq_prefixes` to be set in config (a list of cell prefixes).")
    SKIP_DEMUX_TRIM = True
else:
    raise ValueError(f"Unknown start_from mode: {START_FROM!r}. Use 'raw' or 'trimmed'.")
with open(config["fileindex"]) as f:
    INDICES = [line.strip() for line in f]
print(INDICES)
include: "rules/hiccluster.smk"
#include: "rules/GCHnorm.smk"

# Build rule_all target list. When SKIP_DEMUX_TRIM is set, drop the
# 02.fastqc outputs (which are produced by demultiplex_fastqc_trim).
_rule_all_targets = []
if not SKIP_DEMUX_TRIM:
    _rule_all_targets += expand("02.fastqc_out_snakemake/{prefix}.{index}.R1_fastqc.html", prefix=FASTQ_PREFIXES, index=INDICES)
    _rule_all_targets += expand("02.fastqc_out_snakemake/{prefix}.{index}.R1_fastqc.zip",  prefix=FASTQ_PREFIXES, index=INDICES)
    _rule_all_targets += expand("02.fastqc_out_snakemake/{prefix}.{index}.R2_fastqc.html", prefix=FASTQ_PREFIXES, index=INDICES)
    _rule_all_targets += expand("02.fastqc_out_snakemake/{prefix}.{index}.R2_fastqc.zip",  prefix=FASTQ_PREFIXES, index=INDICES)
# always — downstream targets, applicable in both modes
_rule_all_targets += expand("04.alignment_snakemake/{prefix}.{index}.summary.txt", prefix=FASTQ_PREFIXES, index=INDICES)
_rule_all_targets += expand("07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.hist.txt", prefix=FASTQ_PREFIXES, index=INDICES)
_rule_all_targets += ["06.hiccluster_snakemake/hicluster_250kb_embed_dir/all_merged.pad1_std1_rp0.5_sqrtvc.svd20.hdf5"]
_rule_all_targets += expand("07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.HCG.6plus2.bed", prefix=FASTQ_PREFIXES, index=INDICES)
_rule_all_targets += expand("08.methylation_snakemake/{prefix}.{index}_methylation.txt", prefix=FASTQ_PREFIXES, index=INDICES)

rule all:
    input: _rule_all_targets

rule demultiplex_fastqc_trim: # combination step to trigger re-runs if failed
    input:
        r1 = os.path.join(config["data"], "{prefix}_R1_001.fastq.gz"),
        r2 = os.path.join(config["data"], "{prefix}_R2_001.fastq.gz"),
        index_file = config["fileindex"]
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

        r1_out_trimmed = "03.trimmed_fastq_snakemake/{prefix}.{index}.R1_val_1.fq.gz",
        r1_out_report = "03.trimmed_fastq_snakemake/{prefix}.{index}.R1.fastq.gz_trimming_report.txt",
        r2_out_trimmed = "03.trimmed_fastq_snakemake/{prefix}.{index}.R2_val_2.fq.gz",
        r2_out_report = "03.trimmed_fastq_snakemake/{prefix}.{index}.R2.fastq.gz_trimming_report.txt"
    threads: 8
    resources:
        mem_mb=32000
    benchmark: "benchmarks/00.demultiplex/demul.{prefix}.{index}_benchmark.txt"
    params:
        out_prefix = "01.demul_fastq_snakemake/{prefix}",
        outdir = "02.fastqc_out_snakemake",
        bisulfitehic = config["bisulfitehic"],

        # trimming params
        outdir_fq = "03.trimmed_fastq_snakemake",
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
        r1_trimmed = "03.trimmed_fastq_snakemake/{prefix}.{index}.R1_val_1.fq.gz",
        r2_trimmed = "03.trimmed_fastq_snakemake/{prefix}.{index}.R2_val_2.fq.gz"
    output:
        bam = "04.alignment_snakemake/{prefix}.{index}.bam"
    threads: 8
    resources:
        mem_mb=16000
    benchmark: "benchmarks/02.mapping/alignment.{prefix}.{index}_benchmark.txt"
    params:
        reference = config["reference"],
        picard = config["picard"],
        bisulfitehic = config["bisulfitehic"]
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
        -t {threads} -rgId {wildcards.prefix}.{wildcards.index} -rgSm scNOMeHiC -nonDirectional -pbat -buffer 1000 -enzymeList {params.reference}.DpnII.span_region.bedgraph -outputMateDiffChr \
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
        summary = "04.alignment_snakemake/{prefix}.{index}.summary.txt"
    threads: 8
    resources:
        mem_mb=16000
    benchmark: "benchmarks/qc/qc.{prefix}.{index}_benchmark.txt"
    params:
        bisulfitehic = config["bisulfitehic"]
    log:
        "logs/04.qc_snakemake/qc.{prefix}.{index}.log"
    shell: # is this sort redundant?
        # samtools sort -@ 8 -n -o ${inputfolder}/${sample}_sorted_by_name.calmd.bam ${inputfolder}/${sample}.calmd.bam
        """
        python {params.bisulfitehic}/src/python/mh_reads_summary.v2.py \
        --in_cram {input.sorted_bam} --out_summary {output.summary} \
        2> {log}
        """
        # rm ${inputfolder}/${sample}_sorted_by_name.calmd.bam

rule bisqc:
    input:
        calmd_bam = "04.alignment_snakemake/{prefix}.{index}.calmd.bam"
    output:
        # qc output
        qc_hist = "07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.hist.txt",
        qc_chr21 = "07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.trinuc_methy.chr21.txt",
        qc_chrM = "07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.trinuc_methy.chrM.txt",
        qc_conv_plot = "07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.WCH.bisuflite_conv_distribution_plot.pdf",
        qc_bias_plot = "07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.WCH.methy_bias_plot.pdf",
        qc_cycle = "07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.WCH.methy.cycle.txt"
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
        """

rule bistools:
    input:
        calmd_bam = "04.alignment_snakemake/{prefix}.{index}.calmd.bam"
    output:
        # raw VCFs + per-read methylation tables produced directly by BisulfiteGenotyper
        cyt_vcf       = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.vcf",
        snp_vcf       = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.snp.vcf",
        gch_per_read  = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.GCH.txt",
        hcg_per_read  = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.HCG.txt",
        # post-processed filtered + sorted VCF and per-base BEDs (consumed downstream)
        cyt_filt_vcf       = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.vcf",
        cyt_filt_summary   = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.vcf.cpgSummary.txt",
        gch_6plus2         = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.GCH.6plus2.bed",
        hcg_6plus2         = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.HCG.6plus2.bed",
        gch_strand_6plus2  = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.GCH.strand.6plus2.bed",
        hcg_strand_6plus2  = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.HCG.strand.6plus2.bed"
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
        """

rule rerun_bistools_all:
    input:
        expand("07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.GCH.6plus2.bed",
               prefix=FASTQ_PREFIXES, index=INDICES) +
        expand("07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.HCG.6plus2.bed",
               prefix=FASTQ_PREFIXES, index=INDICES)

rule methylation:
    input:
        meth_HCG_6plus2 = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.cyt.filtered.sort.HCG.6plus2.bed"
    output:
        methylation = "08.methylation_snakemake/{prefix}.{index}_methylation.txt"
    threads: 8
    resources:
        mem_mb=16000
    params:
        sample = "{prefix}.{index}",
        scripts = config["scripts"],
        reference = config["reference"]
    log:
        "logs/08.methylation_snakemake/methylation.{prefix}.{index}.log"
    shell:
        """
        python {params.scripts}/calcmethylation.py \
        --chrom_size_filepath {params.reference}.chrom.sizes \
        --sample_prefix {params.sample} \
        --outfolder 08.methylation_snakemake \
        --bin_size 1000000 \
        2> {log}
        """

        # python calcmethylation.py --chrom_size_filepath reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.chrom.sizes --sample_prefix ${sample} --outfolder 08.methylation --bin_size 25000
