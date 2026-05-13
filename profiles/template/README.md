# Per-project launch template

Copy `config.yaml` and `run.sh` from this directory into a new project
directory (one per sample set), edit the values, and submit.

## Quick start

```bash
PROJECT=/projects/b1198/epifluidlab/yoshii/runs/2026-05-13_cellline_batch5
mkdir -p ${PROJECT}
cp /gpfs/projects/b1198/epifluidlab/yoshii/software/sc_NOMeHiC_pipeline/profiles/template/{config.yaml,run.sh} ${PROJECT}/

# Edit the project's config.yaml — set workdir, start_from, etc.
$EDITOR ${PROJECT}/config.yaml

# Stage your inputs into ${PROJECT}:
#   - start_from: raw    → ${PROJECT}/00.raw_data/{prefix}_R{1,2}_001.fastq.gz
#   - start_from: trimmed → ${PROJECT}/03.trimmed_fastq_snakemake/{prefix}.{index}.R{1,2}_val_{1,2}.fq.gz
# plus ${PROJECT}/00.raw_data/index.txt (one cell index per line)

# Submit:
sbatch ${PROJECT}/run.sh
```

## What's in each file

- **config.yaml** — only the project-specific values. Everything you don't
  override (reference paths, software paths, qc thresholds, trimming params)
  is inherited from `configs/config.yaml` in the pipeline repo. Snakemake
  layers the two `--configfile` args; your project file wins on overlapping
  keys.

- **run.sh** — SLURM batch script. Single allocation that runs `snakemake -j
  16` internally. Resolves the project directory from the script's own path,
  cd's into it, activates the scnomehic conda env (plus appends hic_env's
  bin for `hicluster`), loads Java 17, and runs the workflow.

## Alternative: per-rule SLURM via the slurm profile

For very large runs (many cells × many heavy rules), the
`profiles/slurm/config.yaml` profile lets each rule submit its own SLURM
job — the `snakemake` driver process orchestrates from the login node.
Substitute this for the `snakemake` invocation in `run.sh` (and DON'T
sbatch the whole thing; run it inside tmux/screen):

```bash
snakemake \
  -s ${PIPELINE}/Snakefile \
  --profile ${PIPELINE}/profiles/slurm \
  --configfile ${PIPELINE}/configs/config.yaml \
  --configfile config.yaml \
  -j 100   # max simultaneous SLURM jobs
```
