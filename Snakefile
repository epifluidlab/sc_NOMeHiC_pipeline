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
            prefix = filename.rsplit('_', 2)[0]  # keep everything before the last two underscores
            prefixes.add(prefix)
    return list(prefixes)

FASTQ_PREFIXES = get_fastq_prefixes(config["data"])

with open(config["fileindex"]) as f:
    INDICES = [line.strip() for line in f]

include: "rules/hiccluster.smk"
include: "rules/GCHnorm.smk"

rule all: # does this need only the end output or every single one??
    input:
        # 01.fastqc output
        expand("02.fastqc_out_snakemake/{prefix}.{index}.R1_fastqc.html", prefix=FASTQ_PREFIXES, index=INDICES),
        expand("02.fastqc_out_snakemake/{prefix}.{index}.R1_fastqc.zip", prefix=FASTQ_PREFIXES, index=INDICES),
        
        expand("02.fastqc_out_snakemake/{prefix}.{index}.R2_fastqc.html", prefix=FASTQ_PREFIXES, index=INDICES),
        expand("02.fastqc_out_snakemake/{prefix}.{index}.R2_fastqc.zip", prefix=FASTQ_PREFIXES, index=INDICES),

        # 04.qc output GETS MOVED BY 07.BISTOOLS
        expand("04.alignment_snakemake/{prefix}.{index}.summary.txt", prefix = FASTQ_PREFIXES, index = INDICES),

        # 04b.bisqc output
        expand("07.bistools_snakemake/qc/{prefix}.{index}/{prefix}.{index}.calmd.hist.txt", prefix = FASTQ_PREFIXES, index = INDICES),

        # 06.hiccluster output
        ## generate-matrix output (single cell dependent)
        # expand("06.hiccluster_snakemake/hicluster_250kb_raw_dir/{prefix}.{index}.generatematrix.done", prefix = FASTQ_PREFIXES, index = INDICES),
        ## merge-chromosome-data output (perhaps for later step?)
        # "06.hiccluster_snakemake/hicluster_250kb_embed_dir/all_merged.pad1_std1_rp0.5_sqrtvc.svd20.hdf5",

        # 07.bistools output
        expand("07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.BisSNP-0.90.WCG.coverage.bedgraph", prefix = FASTQ_PREFIXES, index = INDICES),

        # 08.methylation output
        expand("08.methylation_snakemake/{prefix}.{index}_methylation.txt", prefix = FASTQ_PREFIXES, index = INDICES),

        # 09.GCHnorm output
        expand("09.GCHnorm_snakemake/methytab2pbetabinom/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.5kb_interval.pbetabinom.txt", prefix = FASTQ_PREFIXES, index = INDICES)

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
        """
        picard={params.picard}
        module load java/jdk-17.0.2+8
        module load libpng
        BISTOOLS={params.bistools}

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
        # methylation output
        meth_GCH_bissnp = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.BisSNP-0.90.GCH.coverage.bedgraph",
        meth_HCG_bissnp = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.BisSNP-0.90.HCG.coverage.bedgraph",
        meth_GCH_6plus2 = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.6plus2.bed",
        meth_GCH_bedgraph = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.bedgraph",
        meth_GCH_bw = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.bw",
        meth_HCG_6plus2 = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.HCG.6plus2.bed",
        meth_HCG_bedgraph = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.HCG.bedgraph",
        meth_HCG_bw = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.HCG.bw",
        meth_vcf = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.vcf",
        meth_cpgSummary = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.vcf.cpgSummary.txt",
        meth_idx = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.vcf.idx",
        meth_raw_vcf = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.raw.sort.vcf",
        meth_raw_idx = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.raw.sort.vcf.idx",
        meth_raw_MethySummarizeList = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.raw.vcf.MethySummarizeList.txt",
        meth_snp_vcf = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.snp.filtered.sort.vcf",
        meth_snp_cpgSummary = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.snp.filtered.sort.vcf.cpgSummary.txt",
        meth_snp_idx = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.snp.filtered.sort.vcf.idx",
        meth_snp_raw_vcf = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.snp.raw.sort.vcf",
        meth_snp_raw_idx = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.snp.raw.sort.vcf.idx",

        # WCG output
        WCG_coverage = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.BisSNP-0.90.WCG.coverage.bedgraph",
        WCG_filtered_vcf = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.vcf",
        WCG_filtered_cpgSummary = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.vcf.cpgSummary.txt",
        WCG_filtered_idx = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.vcf.idx",
        WCG_filtered_6plus2 = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.WCG.6plus2.bed",
        WCG_filtered_bedgraph = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.WCG.bedgraph",
        WCG_filtered_bw = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.WCG.bw",
        WCG_raw = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.raw.sort.vcf",
        WCG_raw_idx = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.raw.sort.vcf.idx",
        WCG_raw_MethySummarizeList = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.raw.vcf.MethySummarizeList.txt",
        WCG_snp = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.snp.filtered.sort.vcf",
        WCG_snp_cpgSummary = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.snp.filtered.sort.vcf.cpgSummary.txt",
        WCG_snp_idx = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.snp.filtered.sort.vcf.idx",
        WCG_snp_raw = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.snp.raw.sort.vcf",
        WCG_snp_raw_idx = "07.bistools_snakemake/WCG/{prefix}.{index}/{prefix}.{index}.calmd.snp.raw.sort.vcf.idx"
    threads: 16
    resources:
        mem_mb=64000
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
        "logs/07.bistools_snakemake/bistools.{prefix}.{index}"
    shell:
        """
        picard={params.picard}
        module load java/jdk-17.0.2+8
        module load libpng
        BISTOOLS={params.bistools}

        perl {params.bistools}/Bis-SNP/bissnp_easy_usage.pl -use_bad_mates \
        --bistools_path {params.bistools} --nomeseq --lowCov --mmq 30 \
        --nt 1 --mem 10 \
        {params.bistools}/Bis-SNP/BisSNP-0.90.jar \
        {input.calmd_bam} {params.reference}.fa {params.vcf} > {log} 2>&1
        mv {params.indir}*vcf* {params.methylation_outdir}
        mv {params.indir}*.bedgraph {params.methylation_outdir}
        mv {params.indir}*.bed {params.methylation_outdir}
        mv {params.indir}*.bw {params.methylation_outdir}

        perl {params.bistools}/Bis-SNP/bissnp_easy_usage.pl --use_bad_mates \
        --bistools_path {params.bistools} --cytosines WCG,2 --outMode 2 --allC --mmq 30 \
        --nt 12 --mem 60 \
        {params.bistools}/Bis-SNP/BisSNP-0.90.jar \
        {input.calmd_bam} {params.reference}.fa {params.vcf} >> {log} 2>&1
        mv {params.indir}*vcf* {params.WCG_outdir}
        mv {params.indir}*.bedgraph {params.WCG_outdir}
        mv {params.indir}*.bed {params.WCG_outdir}
        mv {params.indir}*.bw {params.WCG_outdir}
        """

rule methylation:
    input:
        meth_HCG_6plus2 = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.HCG.6plus2.bed"
    output:
        methylation = "08.methylation_snakemake/{prefix}.{index}_methylation.txt"
    threads: 8
    resources:
        mem_mb=16000
    params:
        sample = "{prefix}.{index}",
        scripts = config["scripts"]
    log:
        "logs/08.methylation_snakemake/methylation.{prefix}.{index}.log"
    shell:
        """
        python {params.scripts}/calcmethylation.py \
        --chrom_size_filepath reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.chrom.sizes \
        --sample_prefix {params.sample} \
        --outfolder 08.methylation_snakemake \
        --bin_size 1000000 \
        2> {log}
        """

        # python calcmethylation.py --chrom_size_filepath reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.chrom.sizes --sample_prefix ${sample} --outfolder 08.methylation --bin_size 25000
