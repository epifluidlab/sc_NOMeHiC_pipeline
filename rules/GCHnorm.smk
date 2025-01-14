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

with open('00.raw_fastq/index.txt') as f:
    INDICES = [line.strip() for line in f]

# rule all:
#     input:
#         expand("09.GCHnorm_snakemake/methytab2pbetabinom/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.5kb_interval.pbetabinom.txt", prefix = FASTQ_PREFIXES, index = INDICES)

rule bed2bigwig:
    input:
        bed = "07.bistools_snakemake/methylation/{prefix}.{index}/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.6plus2.bed"
    output:
        methy_bw = "09.GCHnorm_snakemake/bed2bigwig/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.methy.bw",
        cov_bw = "09.GCHnorm_snakemake/bed2bigwig/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.cov.bw",
        m_count_bw = "09.GCHnorm_snakemake/bed2bigwig/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.methy_count.bw"
    threads: 1
    resources:
        mem_mb = 16000
    params:
        out_prefix = "09.GCHnorm_snakemake/bed2bigwig/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH",
        scripts = config["scripts"]
    log:
        "logs/09.GCHnorm/bed2bw/bed2bw.{prefix}.{index}.log"
    shell:
        """
        perl {params.scripts}/bed6plus2bw.pl {params.out_prefix} {input.bed} \
        2> {log}
        """

rule bw2tab:
    input:
        counts_bw = "09.GCHnorm_snakemake/bed2bigwig/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.methy_count.bw",
        coverage_bw = "09.GCHnorm_snakemake/bed2bigwig/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.cov.bw"
    output:
        # 5kb tabs
        counts_5kb_tab = "09.GCHnorm_snakemake/bw2tab/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.methy_count.5kb_interval.tab",
        coverage_5kb_tab = "09.GCHnorm_snakemake/bw2tab/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.cov.5kb_interval.tab",

        # 250kb tabs
        # counts_250kb_tab = "09.GCHnorm_snakemake/bed2bigwig/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.methy_count.250kb_interval.tab",
        # coverage_250kb_tab = "09.GCHnorm_snakemake/bed2bigwig/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.cov.250kb_interval.tab"
    threads: 4
    resources:
        mem_mb=32000
    params:
        hg38_5kb = config["hg38_5kb_no_dark_regions"],
        # hg38_250kb = config["hg38_250kb_no_dark_regions"]
    log:
        "logs/09.GCHnorm/bw2tab/bw2tab.{prefix}.{index}.log"
    shell:
        """
        echo "5kb"
        bigWigAverageOverBed \
        {input.counts_bw} \
        {params.hg38_5kb} \
        {output.counts_5kb_tab} \
        2> {log}
        
        bigWigAverageOverBed \
        {input.coverage_bw} \
        {params.hg38_5kb} \
        {output.coverage_5kb_tab} \
        2>> {log}
        """
        # echo "250kb"
        # bigWigAverageOverBed \
        # {input.counts_bw} \
        # {params.hg38_250kb} \
        # {output.counts_250kb_tab}
        
        # bigWigAverageOverBed \
        # {input.coverage_bw} \
        # {params.hg38_250kb} \
        # {output.coverage_250kb_tab}

rule methytab2pbetabinom:
    input:
        counts_5kb_tab = "09.GCHnorm_snakemake/bw2tab/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.methy_count.5kb_interval.tab",
        coverage_5kb_tab = "09.GCHnorm_snakemake/bw2tab/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.cov.5kb_interval.tab",
        # counts_250kb_tab = "09.GCHnorm_snakemake/bed2bigwig/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.methy_count.250kb_interval.tab",
        # coverage_250kb_tab = "09.GCHnorm_snakemake/bed2bigwig/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.cov.250kb_interval.tab"
    output:
        pbetabinom_5kb = "09.GCHnorm_snakemake/methytab2pbetabinom/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.5kb_interval.pbetabinom.txt",
        # pbetabinom_250kb = "09.GCHnorm_snakemake/methytab2pbetabinom/{prefix}.{index}.calmd.cytosine.filtered.sort.GCH.250kb_interval.pbetabinom.txt"
    threads: 1
    resources:
        mem_mb=16000
    params:
        wd = config["workdir"],
        scripts = config["scripts"]
    log:
        "logs/09.GCHnorm/tab2pbeta/tab2pbeta.{prefix}.{index}.log"
    shell:
        """
        R --no-restore --no-save --args \
        wd={params.wd} \
        file1={input.counts_5kb_tab} \
        file2={input.coverage_5kb_tab} \
        output={output.pbetabinom_5kb} \
        minCT=1 \
        logp=T \
        < {params.scripts}/methy_tab_to_pbetabinom.R \
        2> {log}
        """
        # R --no-restore --no-save --args \
        # wd={params.wd} \
        # file1={input.counts_250kb_tab} \
        # file2={input.coverage_250kb_tab} \
        # output={output.pbetabinom_250kb} \
        # minCT=1 \
        # logp=T \
        # < {params.scripts}/methy_tab_to_pbetabinom.R