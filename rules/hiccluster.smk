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

CHROM = [str(c) for c in range(1, 23)]

# ── selected_cells_file ──────────────────────────────────────────────────
# Hicluster's expensive imputation/embedding stage runs ONLY on cells listed
# in this file. The file is a strict whitelist: one "{prefix}.{index}" per
# line. Lines blank or starting with `#` are skipped.
#
# When the user invokes `snakemake hicluster_selected`, this list controls
# which cells participate in concatcells/mergechrom. Cells named here that
# aren't in the discovered CELLS pool (no trimmed fastq pair) trigger a hard
# error at DAG-build time.
#
# Cells that ARE in CELLS but lack upstream outputs (calmd.bam, hic.txt.gz,
# etc.) are NOT silently skipped — snakemake's normal dependency cascade
# triggers the upstream rules to run for them.
SELECTED_CELLS_FILE = config.get("selected_cells_file")


def _load_selected_samples():
    """Strict parse of selected_cells_file → list of "prefix.index" strings.
    Returns [] when the config option is unset. Fails loudly on a missing
    file or unknown cell names so DAG construction never silently drops
    work the user asked for.

    A whitelisted cell is considered valid if EITHER:
      (a) it's in CELLS (trimmed fastq pair present → full pipeline can run), or
      (b) its `04.alignment_snakemake/{cell}.calmd.bam` exists on disk
          (hicluster's chain starts from the BAM, so the trimmed fastq isn't
          strictly required — useful when /scratch was purged but BAMs
          survive in the project workdir).
    A cell that satisfies neither is reported as unknown."""
    if not SELECTED_CELLS_FILE:
        return []
    if not os.path.exists(SELECTED_CELLS_FILE):
        raise FileNotFoundError(
            f"selected_cells_file does not exist: {SELECTED_CELLS_FILE}"
        )
    with open(SELECTED_CELLS_FILE) as fh:
        samples = []
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            samples.append(line)
    cells_set = set(SAMPLES)
    unknown = []
    for s in samples:
        if s in cells_set:
            continue
        if os.path.exists(f"04.alignment_snakemake/{s}.calmd.bam"):
            continue
        unknown.append(s)
    if unknown:
        raise ValueError(
            f"selected_cells_file ({SELECTED_CELLS_FILE}) lists "
            f"{len(unknown)} cell(s) with no upstream inputs available "
            f"(neither a trimmed fastq pair nor a .calmd.bam). "
            f"First few: {unknown[:5]}"
        )
    if not samples:
        raise ValueError(
            f"selected_cells_file ({SELECTED_CELLS_FILE}) is empty or "
            f"contains only comments."
        )
    return samples


SELECTED_SAMPLES = _load_selected_samples()

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

# Selected-cells whitelist is parsed at module load (see SELECTED_SAMPLES
# above); we materialize it to disk so it shows up as a real Snakemake target
# and downstream rules have a clean file-based dependency, but the contents
# come from the user-provided config file, not from a runtime check.
rule write_selected_cells:
    output:
        selected = "06.hiccluster_snakemake/hicluster_check/selected_cells.txt"
    run:
        os.makedirs(os.path.dirname(output.selected), exist_ok=True)
        with open(output.selected, "w") as fh:
            for s in SELECTED_SAMPLES:
                fh.write(s + "\n")


def _selected_imputed_hdf5s(wildcards):
    """Inputs for imputelist: every selected cell's imputed hdf5 for this
    chromosome. Uses SELECTED_SAMPLES (parsed at DAG-build time from the
    user's selected_cells_file) — no runtime checkpoint needed."""
    return [
        f"06.hiccluster_snakemake/hicluster_250kb_impute_dir/chr{wildcards.chr}/"
        f"{sample}_chr{wildcards.chr}_pad1_std1_rp0.5_sqrtvc.hdf5"
        for sample in SELECTED_SAMPLES
    ]


rule imputecell:
    input:
        hic_matrix = "06.hiccluster_snakemake/{sample}.hic_matrix.txt.gz",
        matrixdone = "06.hiccluster_snakemake/hicluster_250kb_raw_dir/{sample}.generatematrix.done"
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
        _selected_imputed_hdf5s
    output:
        file_list = "06.hiccluster_snakemake/hicluster_250kb_chr{chr}_impute_file_list.txt"
    threads: 1
    log:
        "logs/06.hiccluster_snakemake/imputelist/imputelist.chr{chr}.log"
    run:
        # Write the (DAG-build-time) list of imputed hdf5s for this chrom.
        # Snakemake guarantees they all exist by the time this rule runs.
        with open(output.file_list, "w") as fh:
            for p in input:
                fh.write(p + "\n")

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


# Friendly entry point. After populating `selected_cells_file` in config,
# trigger the whole hicluster cascade for the selected cells with:
#
#     snakemake hicluster_selected --configfile config.yaml --profile slurm
#
# Snakemake then runs scbam2hic → hicprocess → generatematrix → imputecell
# → imputelist → concatcells → mergechrom for every cell in the whitelist
# (and skips anything already done).
rule hicluster_selected:
    input:
        "06.hiccluster_snakemake/hicluster_250kb_embed_dir/all_merged.pad1_std1_rp0.5_sqrtvc.svd20.hdf5",
        "06.hiccluster_snakemake/hicluster_check/selected_cells.txt"
