#!/bin/bash
# Per-project pipeline launcher (per-rule SLURM pattern).
#
# Drop this file into your project directory (next to that project's
# config.yaml) and launch from the login node (inside tmux/screen, since
# the driver process must stay alive):
#   bash snakemake_command.sh
#
# Snakemake runs on the login node and submits each rule as its own SLURM
# job via the slurm profile. Use this when you have many cells × heavy
# rules and want each to size its own SLURM allocation. For smaller runs,
# prefer the single-sbatch run.sh template alongside this file.

# Activate scnomehic so the snakemake binary and snakemake-executor-plugin-slurm
# are importable. Idempotent — safe if already active. Without this you'll get
# `argument --executor/-e: invalid choice: 'slurm'`.
source /projects/b1198/epifluidlab/yoshii/software/conda/etc/profile.d/conda.sh
conda activate scnomehic

PIPELINE=/gpfs/projects/b1198/epifluidlab/yoshii/software/sc_NOMeHiC_pipeline

snakemake \
  -s ${PIPELINE}/Snakefile \
  --profile ${PIPELINE}/profiles/slurm \
  --configfile config.yaml \
  -j 1000
