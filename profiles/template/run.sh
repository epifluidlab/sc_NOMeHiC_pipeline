#!/bin/bash
# Per-project pipeline launcher (single-sbatch pattern).
#
# Drop this file into your project directory (next to that project's
# config.yaml) and submit with:
#   sbatch run.sh
#
# It runs snakemake inside one SLURM allocation that parallelizes rules via
# `-j N` internally. Simpler than the per-rule `--profile slurm` pattern;
# fine for runs that fit in a single multi-core node.

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name=scnomehic
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err

set -e

# Locations — adjust PIPELINE if the pipeline repo moves.
PIPELINE=/gpfs/projects/b1198/epifluidlab/yoshii/software/sc_NOMeHiC_pipeline

# This project's directory (the one containing this run.sh + config.yaml).
PROJECT="$(dirname "$(readlink -f "$0")")"
cd "${PROJECT}"
mkdir -p logs

# Conda env: scnomehic provides samtools/bedtools/perl/python/etc.
# hic_env's bin appended so `hicluster` resolves for hiccluster rules.
source /projects/b1198/epifluidlab/yoshii/software/conda/etc/profile.d/conda.sh
conda activate scnomehic
export PATH="${PATH}:/projects/b1198/epifluidlab/yoshii/software/conda/envs/hic_env/bin"

# Java 17 for bhmem (mapping). bisqc and bistools rules internally use Java 8
# via the wrappers' hardcoded path, so this module load is fine for both.
module load java/jdk-17.0.2+8

# Snakemake invocation. `--configfile ${PIPELINE}/configs/config.yaml` first
# loads the repo defaults; `--configfile config.yaml` then overlays this
# project's values.
snakemake \
  -s ${PIPELINE}/Snakefile \
  --configfile ${PIPELINE}/configs/config.yaml \
  --configfile config.yaml \
  -j ${SLURM_CPUS_PER_TASK:-16} \
  -p --keep-going
