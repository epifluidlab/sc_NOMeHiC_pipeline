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

# ── autosome list ────────────────────────────────────────────────────────
# hicluster's per-chromosome imputation/embedding runs over every chromosome
# listed here. Override via config["hicluster_chromosomes"] (list of strings
# without "chr" prefix, e.g. ["1", "2", ..., "19"]). Default: auto-detect
# autosomes from {reference}.chrom.sizes — picks every "chr<N>" with integer N,
# so mm10 produces 1..19 and hg38 produces 1..22 without any config change.
def _autosomes_from_chrom_sizes(path):
    import re
    autosomes = set()
    with open(path) as fh:
        for line in fh:
            chrom = line.split('\t', 1)[0].strip()
            m = re.match(r'^chr([0-9]+)$', chrom)
            if m:
                autosomes.add(int(m.group(1)))
    if not autosomes:
        raise ValueError(f"No 'chr<N>' autosomes found in {path}")
    return [str(i) for i in sorted(autosomes)]

_CHROM_CFG = config.get("hicluster_chromosomes")
if _CHROM_CFG:
    CHROM = [str(c) for c in _CHROM_CFG]
else:
    CHROM = _autosomes_from_chrom_sizes(config["reference"] + ".chrom.sizes")

# ── hicluster resolutions ────────────────────────────────────────────────
# Hicluster's imputation+embedding is run at every resolution listed here.
# Override via config["hicluster_resolutions"] (list of int bp values).
# Default produces both 100kb and 250kb embeddings in one invocation.
def _res_label(bp):
    bp = int(bp)
    if bp >= 1_000_000 and bp % 1_000_000 == 0:
        return f"{bp // 1_000_000}Mb"
    if bp >= 1_000 and bp % 1_000 == 0:
        return f"{bp // 1_000}kb"
    return f"{bp}bp"

RESOLUTIONS_BP = [int(x) for x in config.get("hicluster_resolutions", [100000, 250000])]
RES_LABELS = [_res_label(b) for b in RESOLUTIONS_BP]
RES_BP_BY_LABEL = dict(zip(RES_LABELS, RESOLUTIONS_BP))

# ── cell-set label ─────────────────────────────────────────────────────────
# Optional suffix that separates the CELL-SET-DEPENDENT combining outputs
# (imputelist file_list, concatcells/mergechrom embed_dir, selected_cells.txt)
# of runs that share the SAME workdir but use different cell whitelists. The
# per-cell outputs (scbam2hic, hicprocess, generatematrix raw_dir, imputecell
# hdf5s) are keyed by cell name and stay UNsuffixed so they're reused across
# cell-sets. Empty (default) -> current unsuffixed paths (backward compatible).
# Set via config["cellset_label"] (e.g. "all3", "batch23").
CELLSET_LABEL = str(config.get("cellset_label", "") or "")
_CS = f"_{CELLSET_LABEL}" if CELLSET_LABEL else ""

wildcard_constraints:
    res_label = "|".join(RES_LABELS)

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
    threads: 16
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
        samtools view -@ {threads} -bh -q 30 -f 1 -F 1804 {input.calmd_bam} > {output.good_reads} && \
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
        directory("06.hiccluster_snakemake/hicluster_{res_label}_raw_dir/{prefix}.{index}"),
        donefile = touch("06.hiccluster_snakemake/hicluster_{res_label}_raw_dir/{prefix}.{index}.generatematrix.done")
    conda:
        "../envs/schicluster_test.yaml"
    threads: 1
    params:
        reference = config["reference"],
        outdir = "06.hiccluster_snakemake/hicluster_{res_label}_raw_dir/{prefix}.{index}/",
        rawdir = "06.hiccluster_snakemake/hicluster_{res_label}_raw_dir/",
        res_bp = lambda w: RES_BP_BY_LABEL[w.res_label],
        shortcut = "{prefix}.{index}"
    log:
        "logs/06.hiccluster_snakemake/generatematrix/generatematrix.{prefix}.{index}.{res_label}.log"
    shell:
        ###hicluster generatematrix-cell
        """
        hicluster generatematrix-cell \
        --infile {input.hic_matrix} --outdir {params.outdir} \
        --chrom_file {params.reference}.chrom.sizes \
        --res {params.res_bp} --cell {params.shortcut} --chr1 1 --pos1 2 --chr2 5 --pos2 6 \
        2> {log}

        cp -r {params.outdir}* {params.rawdir}
        """

# Selected-cells whitelist is parsed at module load (see SELECTED_SAMPLES
# above); we materialize it to disk so it shows up as a real Snakemake target
# and downstream rules have a clean file-based dependency, but the contents
# come from the user-provided config file, not from a runtime check.
rule write_selected_cells:
    output:
        selected = f"06.hiccluster_snakemake/hicluster_check/selected_cells{_CS}.txt"
    run:
        os.makedirs(os.path.dirname(output.selected), exist_ok=True)
        with open(output.selected, "w") as fh:
            for s in SELECTED_SAMPLES:
                fh.write(s + "\n")


def _selected_imputed_hdf5_paths(wildcards):
    """Per-cell imputed hdf5 paths for this (res, chrom). Pure path-list helper
    — used both as DAG inputs (default) and as runtime params (when
    concat_from_existing_imputes is set)."""
    return [
        f"06.hiccluster_snakemake/hicluster_{wildcards.res_label}_impute_dir/chr{wildcards.chr}/"
        f"{sample}_chr{wildcards.chr}_pad1_std1_rp0.5_sqrtvc.hdf5"
        for sample in SELECTED_SAMPLES
    ]


# Config flag: when True, the imputelist rule treats the per-cell hdf5s as
# already-on-disk inputs (not as DAG-tracked outputs), bypassing the slow
# imputecell rule entirely. Used to re-run concatcells + mergechrom for a new
# cell list when the per-cell impute hdf5s already exist (e.g., a previous run
# at a different resolution / cell whitelist).
_CONCAT_FROM_EXISTING = bool(config.get("concat_from_existing_imputes", False))


def _selected_imputed_hdf5s(wildcards):
    """Inputs for imputelist: every selected cell's imputed hdf5 for this
    (resolution, chromosome). Uses SELECTED_SAMPLES (parsed at DAG-build time
    from the user's selected_cells_file) — no runtime checkpoint needed.

    When config.concat_from_existing_imputes is set, returns [] so snakemake's
    DAG won't try to (re-)build the hdf5s — caller asserts they already exist."""
    if _CONCAT_FROM_EXISTING:
        return []
    return _selected_imputed_hdf5_paths(wildcards)


rule imputecell:
    input:
        hic_matrix = "06.hiccluster_snakemake/{sample}.hic_matrix.txt.gz",
        matrixdone = "06.hiccluster_snakemake/hicluster_{res_label}_raw_dir/{sample}.generatematrix.done"
    output:
        imputed_cells = "06.hiccluster_snakemake/hicluster_{res_label}_impute_dir/chr{chr}/{sample}_chr{chr}_pad1_std1_rp0.5_sqrtvc.hdf5"
    conda:
        "../envs/schicluster_test.yaml"
    threads: 1
    params:
        reference = config["reference"],
        indir = "06.hiccluster_snakemake/hicluster_{res_label}_raw_dir/chr{chr}/",
        outdir = "06.hiccluster_snakemake/hicluster_{res_label}_impute_dir/chr{chr}/",
        res_bp = lambda w: RES_BP_BY_LABEL[w.res_label],
        shortcut = "{sample}"
    log:
        "logs/06.hiccluster_snakemake/imputecell/imputecell_{sample}/imputecell.{sample}.{res_label}.chr{chr}.log"
    shell:
        ###hicluster impute-cell
        """
        hicluster impute-cell \
        --indir {params.indir} \
        --outdir {params.outdir} \
        --cell {params.shortcut} --chrom chr{wildcards.chr} --res {params.res_bp} \
        --chrom_file {params.reference}.chrom.sizes 2> {log}
        """

rule imputelist:
    input:
        _selected_imputed_hdf5s
    output:
        file_list = f"06.hiccluster_snakemake/hicluster_{{res_label}}_chr{{chr}}_impute_file_list{_CS}.txt"
    threads: 1
    params:
        # ALWAYS the per-cell hdf5 paths (regardless of concat_from_existing_imputes).
        # In default mode, these match `input` 1-for-1. In concat-only mode, `input`
        # is empty but params.paths still has the full path list to write.
        paths = _selected_imputed_hdf5_paths,
    log:
        "logs/06.hiccluster_snakemake/imputelist/imputelist.{res_label}.chr{chr}.log"
    run:
        # Write the (DAG-build-time) list of imputed hdf5s for this (res, chrom).
        # Default mode: snakemake's DAG already guaranteed each path exists.
        # concat_from_existing_imputes mode: caller asserts they exist; we
        # fail loudly here if any are missing rather than concatcells later.
        import os
        missing = [p for p in params.paths if not os.path.exists(p)]
        if missing:
            raise FileNotFoundError(
                f"imputelist ({wildcards.res_label} chr{wildcards.chr}): "
                f"{len(missing)} of {len(params.paths)} expected per-cell hdf5s "
                f"missing. First few: {missing[:5]}"
            )
        with open(output.file_list, "w") as fh:
            for p in params.paths:
                fh.write(p + "\n")

rule concatcells:
    input:
        file_list = f"06.hiccluster_snakemake/hicluster_{{res_label}}_chr{{chr}}_impute_file_list{_CS}.txt"
    output:
        npz = f"06.hiccluster_snakemake/hicluster_{{res_label}}_embed_dir{_CS}/pad1_std1_rp0.5_sqrtvc_chr{{chr}}.npz",
        svd50_npy = f"06.hiccluster_snakemake/hicluster_{{res_label}}_embed_dir{_CS}/pad1_std1_rp0.5_sqrtvc_chr{{chr}}.svd50.npy"
    conda:
        "../envs/schicluster_test.yaml"
    threads: 16
    # Concatenating 1000+ cells of a large chromosome at 100kb resolution then
    # running SVD is memory-heavy; the snakemake-profile default (1000 MB)
    # was triggering OOM kills on chr10 100kb.
    resources:
        mem_mb=128000
    params:
        outprefix = f"06.hiccluster_snakemake/hicluster_{{res_label}}_embed_dir{_CS}/pad1_std1_rp0.5_sqrtvc_chr{{chr}}",
        res_bp = lambda w: RES_BP_BY_LABEL[w.res_label]
    shell:
        ###concatenate to form a matrix
        # concatenate cells for each chromosome
        """
        export OMP_NUM_THREADS={threads}
        export MKL_NUM_THREADS={threads}
        export OPENBLAS_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        hicluster embed-concatcell-chr \
        --cell_list {input.file_list} \
        --outprefix {params.outprefix} --res {params.res_bp}
        """

rule mergechrom:
    input:
        expand(f"06.hiccluster_snakemake/hicluster_{{{{res_label}}}}_embed_dir{_CS}/pad1_std1_rp0.5_sqrtvc_chr{{chr}}.svd50.npy", chr = CHROM)
    output:
        embed_file_list = f"06.hiccluster_snakemake/hicluster_{{res_label}}_embed_file_list{_CS}.txt",
        all_merged = f"06.hiccluster_snakemake/hicluster_{{res_label}}_embed_dir{_CS}/all_merged.pad1_std1_rp0.5_sqrtvc.svd20.hdf5"
    conda:
        "../envs/schicluster_test.yaml"
    threads: 1
    params:
        outprefix = f"06.hiccluster_snakemake/hicluster_{{res_label}}_embed_dir{_CS}/all_merged.pad1_std1_rp0.5_sqrtvc",
        concat_files = f"06.hiccluster_snakemake/hicluster_{{res_label}}_embed_dir{_CS}/pad1_std1_rp0.5_sqrtvc_*npy"
    log:
        "logs/06.hiccluster_snakemake/mergechrom/mergechrom.{res_label}.log"
    shell:
        # merge chromosomes per resolution
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
        expand(
            f"06.hiccluster_snakemake/hicluster_{{res_label}}_embed_dir{_CS}/all_merged.pad1_std1_rp0.5_sqrtvc.svd20.hdf5",
            res_label=RES_LABELS,
        ),
        f"06.hiccluster_snakemake/hicluster_check/selected_cells{_CS}.txt"
