import yaml
import os
import glob

# Define config file
configfile: "configs/config.yaml"

# Set working directory
workdir: config["workdir"]

# FASTQ_PREFIXES, INDICES, and CELLS come from the main Snakefile (the include
# happens after they are set). SAMPLES here is derived from CELLS so we only
# carry the (prefix, index) pairs that actually exist on disk.
SAMPLES = [f"{p}.{i}" for p, i in CELLS]

# CHROM = [str(c) for c in range(1, 23)] + ['EBV', 'M', 'Un', 'X', 'Y']
CHROM = [str(c) for c in range(1, 23)]

def schic_exclude_missing_chrom():
    schic_passed = []
    schic_excluded = []
    for sample in SAMPLES:
        chrom_exists_sum = 0
        for chr in CHROM:
            if os.path.exists(f"06.hiccluster_snakemake/hicluster_250kb_raw_dir/chr{chr}/{sample}_chr{chr}.txt"):
                chrom_exists_sum+=1
        if chrom_exists_sum==len(CHROM):
            schic_passed.append(sample)
        else:
            schic_excluded.append(sample)
    return schic_passed, schic_excluded

# keep commented out if you want to call it from the main snakefile
# rule all:
#     input:
#         # mergechrom
#         "06.hiccluster_snakemake/hicluster_250kb_embed_dir/all_merged.pad1_std1_rp0.5_sqrtvc.svd20.hdf5"

rule scbam2hic:
    input:
        calmd_bam = "04.alignment_snakemake/{prefix}.{index}.calmd.bam"
    output:
        good_reads = "06.hiccluster_snakemake/{prefix}.{index}.good_reads.bam",
        hic = "06.hiccluster_snakemake/{prefix}.{index}.hic.txt.gz"
    threads: 5
    params:
        reference = config["reference"],
        restriction_sites = config["restriction_sites"],
        bisulfitehic = config["bisulfitehic"]
    log:
        "logs/06.hiccluster_snakemake/scbam2hic/scbam2hic.{prefix}.{index}.log"
    shell:
        ###generate contact matrix file for each single-cell
        # sam2juicer_new.py writes Hi-C pairs to stdout; pipe through gzip
        # to save ~5x (these files are ~500MB-1GB plain text per cell).
        """
        samtools view -bh -q 30 -f 1 -F 1804 {input.calmd_bam} > {output.good_reads} && \
        python {params.bisulfitehic}/src/python/sam2juicer_new.py \
        -s {output.good_reads} -f {params.restriction_sites} 2> {log} | gzip -nc > {output.hic}
        """
        # maybe use GCA restriction sites?

rule hicprocess:
    # preprocess.hg38.sh uses `zcat -f` so it transparently handles both
    # .hic.txt.gz (current scbam2hic output) and legacy plain .hic.txt.
    input:
        hic = "06.hiccluster_snakemake/{prefix}.{index}.hic.txt.gz"
    output:
        hic_matrix = "06.hiccluster_snakemake/{prefix}.{index}.hic_matrix.txt.gz"
    threads: 1
    params:
        scripts = config["scripts"]
    log:
        "logs/06.hiccluster_snakemake/hicprocess/hicprocess.{prefix}.{index}.log"
    shell:
        """
        {params.scripts}/preprocess.hg38.sh \
        {input.hic} {output.hic_matrix} \
        2> {log}
        """

rule generatematrix:
    input:
        hic_matrix = "06.hiccluster_snakemake/{prefix}.{index}.hic_matrix.txt.gz"
    output:
        directory("06.hiccluster_snakemake/hicluster_250kb_raw_dir/{prefix}.{index}"),
        donefile = touch("06.hiccluster_snakemake/hicluster_250kb_raw_dir/{prefix}.{index}.generatematrix.done")
    conda:
        "../envs/schicluster_test.yaml"
    threads: 1
    params:
        reference = config["reference"],
        outdir = "06.hiccluster_snakemake/hicluster_250kb_raw_dir/{prefix}.{index}/",
        shortcut = "{prefix}.{index}"
    log:
        "logs/06.hiccluster_snakemake/generatematrix/generatematrix.{prefix}.{index}.log"
    shell: 
        ###hicluster generatematrix-cell
        """
        hicluster generatematrix-cell \
        --infile {input.hic_matrix} --outdir {params.outdir} \
        --chrom_file {params.reference}.chrom.sizes \
        --res 250000 --cell {params.shortcut} --chr1 1 --pos1 2 --chr2 5 --pos2 6 \
        2> {log}

        cp -r {params.outdir}* 06.hiccluster_snakemake/hicluster_250kb_raw_dir/
        """

checkpoint excludemissing:
    input:
        [f"06.hiccluster_snakemake/hicluster_250kb_raw_dir/{p}.{i}.generatematrix.done" for p, i in CELLS]
    output:
        checkpoint_done = touch("06.hiccluster_snakemake/hicluster_check/check_missing_chrom.done"),
        passed = "06.hiccluster_snakemake/hicluster_check/passed_cells.txt",
        excluded = "06.hiccluster_snakemake/hicluster_check/excluded_cells.txt"
    run:
        PASSED, EXCLUDED = schic_exclude_missing_chrom()

        # write to file
        with open("06.hiccluster_snakemake/hicluster_check/passed_cells.txt", "w") as f:
            for cell in PASSED:
                f.write(cell + "\n")

        # write to file CHANGE THIS TO THE NAME THAT'S NEEDED
        with open("06.hiccluster_snakemake/hicluster_check/excluded_cells.txt", "w") as f:
            for cell in EXCLUDED:
                f.write(cell + "\n")
        
        print("Excluded cells:", EXCLUDED)

def get_passed_samples(wildcards):
    checkpoint_output = checkpoints.excludemissing.get().output.passed
    with open(checkpoint_output) as f:
        passed_samples = [line.strip() for line in f]
    return passed_samples

rule imputecell: # USE DIFFERENT WILDCARDS FOR THIS RULE
    input:
        hic_matrix = "06.hiccluster_snakemake/{sample}.hic_matrix.txt.gz",
        matrixdone = "06.hiccluster_snakemake/hicluster_250kb_raw_dir/{sample}.generatematrix.done",
        checkpoint_done = "06.hiccluster_snakemake/hicluster_check/check_missing_chrom.done"
    output:
        imputed_cells = "06.hiccluster_snakemake/hicluster_250kb_impute_dir/chr{chr}/{sample}_chr{chr}_pad1_std1_rp0.5_sqrtvc.hdf5"
    conda:
        "../envs/schicluster_test.yaml"
    threads: 1
    params:
        reference = config["reference"],
        indir = "06.hiccluster_snakemake/hicluster_250kb_raw_dir/chr{chr}/",
        outdir = "06.hiccluster_snakemake/hicluster_250kb_impute_dir/chr{chr}/",
        shortcut = "{sample}"
    log:
        "logs/06.hiccluster_snakemake/imputecell/imputecell_{sample}/imputecell.{sample}.chr{chr}.log"
    shell:
        ###hicluster impute-cell
        """
        hicluster impute-cell \
        --indir {params.indir} \
        --outdir {params.outdir} \
        --cell {params.shortcut} --chrom chr{wildcards.chr} --res 250000 \
        --chrom_file {params.reference}.chrom.sizes 2> {log}
        """

rule imputelist:
    input:
        lambda wildcards: [
            f"06.hiccluster_snakemake/hicluster_250kb_impute_dir/chr{wildcards.chr}/{sample}_chr{wildcards.chr}_pad1_std1_rp0.5_sqrtvc.hdf5"
            for sample in [line.strip() for line in open(checkpoints.excludemissing.get().output.passed)]
        ]
    output:
        file_list = "06.hiccluster_snakemake/hicluster_250kb_chr{chr}_impute_file_list.txt"
    conda:
        "../envs/schicluster_test.yaml"
    threads: 1
    params:
        imputed_cells = "06.hiccluster_snakemake/hicluster_250kb_impute_dir/chr{chr}/*chr{chr}_pad1_std1_rp0.5_sqrtvc.hdf5"
    log:
        "logs/06.hiccluster_snakemake/imputelist/imputelist.chr{chr}.log"
    shell:
        """
        ls {params.imputed_cells} > {output.file_list} 2> {log}
        """

rule concatcells:
    input:
        file_list = "06.hiccluster_snakemake/hicluster_250kb_chr{chr}_impute_file_list.txt"
    output:
        npz = "06.hiccluster_snakemake/hicluster_250kb_embed_dir/pad1_std1_rp0.5_sqrtvc_chr{chr}.npz",
        svd50_npy = "06.hiccluster_snakemake/hicluster_250kb_embed_dir/pad1_std1_rp0.5_sqrtvc_chr{chr}.svd50.npy"
    conda:
        "../envs/schicluster_test.yaml"
    threads: 1
    params:
        outprefix = "06.hiccluster_snakemake/hicluster_250kb_embed_dir/pad1_std1_rp0.5_sqrtvc_chr{chr}"
    shell:
        ###concatnate to form a matrix
        # concatenate cells for each chromosome
        """
        hicluster embed-concatcell-chr \
        --cell_list {input.file_list} \
        --outprefix {params.outprefix} --res 250000
        """

rule mergechrom:
    input:
        expand("06.hiccluster_snakemake/hicluster_250kb_embed_dir/pad1_std1_rp0.5_sqrtvc_chr{chr}.svd50.npy", chr = CHROM)
    output:
        embed_file_list = "06.hiccluster_snakemake/hicluster_250kb_embed_file_list.txt",
        all_merged = "06.hiccluster_snakemake/hicluster_250kb_embed_dir/all_merged.pad1_std1_rp0.5_sqrtvc.svd20.hdf5"
    conda:
        "../envs/schicluster_test.yaml"
    threads: 1
    params:
        outprefix = "06.hiccluster_snakemake/hicluster_250kb_embed_dir/all_merged.pad1_std1_rp0.5_sqrtvc",
        concat_files = "06.hiccluster_snakemake/hicluster_250kb_embed_dir/pad1_std1_rp0.5_sqrtvc_*npy"
    log:
        "logs/06.hiccluster_snakemake/mergechrom/mergechrom.log"
    shell:
        # merge chromosomes (hicluster_100kb_embed_dir/pad1_std1_rp0.5_sqrtvc.svd50.hdf5)
        """
        ls {params.concat_files} > {output.embed_file_list} 2> {log}
        hicluster embed-mergechr --embed_list {output.embed_file_list} --outprefix {params.outprefix} 2>> {log}
        """
