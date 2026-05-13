# Per-project launch template

Copy `config.yaml` and `run.sh` from this directory into a new project
directory (one per sample set), edit the values in `config.yaml`, and submit.

## Quick start

```bash
PROJECT=/projects/b1198/epifluidlab/yoshii/runs/2026-05-13_cellline_batch5
mkdir -p ${PROJECT}
cp /gpfs/projects/b1198/epifluidlab/yoshii/software/sc_NOMeHiC_pipeline/profiles/template/{config.yaml,run.sh} ${PROJECT}/

# Edit only the values you care about (workdir, start_from, inputs).
# Everything else has sensible defaults pre-filled.
$EDITOR ${PROJECT}/config.yaml

# Stage your inputs into ${PROJECT}, then submit:
sbatch ${PROJECT}/run.sh
```

## What each file is

- **config.yaml** — self-contained pipeline config. Pass this as the ONLY
  `--configfile` to snakemake; no need to also reference the repo's
  `configs/config.yaml`. Fields are grouped:
  1. **Project-specific** (edit per run): workdir, start_from, data/trimmed_data, fileindex, fastq_prefixes(_file)
  2. **Shared paths** (rarely change): reference, picard, bisulfitehic, bistools, etc.
  3. **Tuning** (rarely change): trimming params, qc thresholds

- **run.sh** — SLURM batch script. Single allocation that runs `snakemake -j 16`
  internally. Resolves its own project dir, cd's into it, activates the
  `scnomehic` conda env (plus appends `hic_env`'s bin for `hicluster`), loads
  Java 17, and runs the workflow against `config.yaml` in the same dir.

## Two starting modes

**Raw fastqs** (default — runs demultiplex + fastqc + trim, then mapping onward):

```yaml
workdir: "/projects/.../runs/myrun"
start_from: "raw"
data: "00.raw_data/"
fileindex: "00.raw_data/index.txt"
# trimmed_data defaults to "03.trimmed_fastq_snakemake"
```

Then stage raw inputs to `${PROJECT}/00.raw_data/{prefix}_R{1,2}_001.fastq.gz`
and `${PROJECT}/00.raw_data/index.txt`.

**Already-trimmed fastqs** (skip demux + fastqc + trim):

```yaml
workdir: "/projects/.../runs/myrun"
start_from: "trimmed"
fileindex: "00.raw_data/index.txt"   # still needed for INDICES
trimmed_data: "/scratch/.../cellline_fastq"   # wherever the trimmed fastqs live
# auto-discover scans trimmed_data; or set fastq_prefixes / fastq_prefixes_file
```

Then ensure your trimmed inputs exist at
`${trimmed_data}/{prefix}.{index}.R1_val_1.fq.gz` (and `R2_val_2.fq.gz`)
for each index in `fileindex`.

## Alternative: per-rule SLURM via the slurm profile

For very large runs (many cells × many heavy rules), the
`profiles/slurm/config.yaml` profile lets each rule submit its own SLURM
job — the `snakemake` driver runs from the login node. Substitute this for
the `snakemake` invocation in `run.sh` (and DON'T sbatch the whole thing;
run inside tmux/screen):

```bash
PIPELINE=/gpfs/projects/b1198/epifluidlab/yoshii/software/sc_NOMeHiC_pipeline
snakemake \
  -s ${PIPELINE}/Snakefile \
  --profile ${PIPELINE}/profiles/slurm \
  --configfile config.yaml \
  -j 100   # max simultaneous SLURM jobs
```
