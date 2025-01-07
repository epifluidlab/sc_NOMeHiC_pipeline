# scNOMeHiC Data Processing Pipeline (Snakemake)
scNOMeHiC is a method that allows for simultaneous capture of chromatin accessibility, DNA methylation, and 3D genome interactions within singular molecules at a single cell resolution via combination of HiC, NOMeSeq, and BisulfiteSeq techniques. This scNOMeHiC pipeline, developed by Yaping Liu and Yoshii Ma, processes the paired read FASTQ files into GCH and HCG methylation calls, as well as single cell HiC embeddings to be used to single cell clustering and analysis. It has been built in Snakemake 8.25 to consolidate each processing step and allow for easy cluster execution and resource management.

## Prerequisites
Please download/build the following modules from their respective online resources:

- FastQC: https://github.com/s-andrews/FastQC
- Juicer: https://github.com/aidenlab/juicer?tab=readme-ov-file
- Picard: https://github.com/broadinstitute/picard
- TrimGalore: https://github.com/FelixKrueger/TrimGalore
- BedTools: https://github.com/arq5x/bedtools2
- Mustache: https://github.com/mustache/mustache
- Bismark: https://github.com/FelixKrueger/Bismark

NOTE: The following prerequisites have a chance of errors in being added to $PATH or being properly utilized when creating the conda environment from the environment.yaml file. If errors persist when running the pipeline after creating the conda environment, see resources below:

- SamTools: https://github.com/samtools/samtools?tab=readme-ov-file
- Cutadapt: https://cutadapt.readthedocs.io/en/stable/

## Creating the Conda Environment/Installing Snakemake
We must first create the scNOMeHiC conda environment in order to both install Snakemake and scNOMeHiC dependencies, and run the pipeline

### Step 1: Install Miniforge
It is recommended to install Miniforge as it includes the Mamba package manager

Please open a terminal and run the following on the platform of your choice:

**Linux or Windows Subsystem for Linux**
```
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

**Mac with x86_64 architecture**
```
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh -o Miniforge3-MacOSX-x86_64.sh
bash Miniforge3-MacOSX-x86_64.sh
```

**Mac with ARM/M1 architecture**
```
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh -o Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh
```

### Step 2: Create your workflow directory
Choose a desired location for all your data files, scNOMeHiC workflow files, and processed data files to be stored. Make sure to have plenty of space and genomic data can be large.
```
mkdir <path/to/folder>
cd <path/to/folder>
```

Clone the files from this directory into your workflow folder
```
git clone https://github.com/yoshihiko1218/sc_nomehic_pipeline_snakemake .
```

### Step 3: Create Snakemake scNOMeHiC environment for required software
First, activate the conda base environment
```
conda activate base
```

Then, using the environment.yaml file within the envs folder downloaded from this repository, we can create our conda environment and all required software (besides prerequisites)
```
mamba env create --name scnomehic --file envs/environment.yaml
```

If you don't have Mamba due to using a different conda distribution than Miniforge, run the following command to get mamba
```
conda install -n base -c conda-forge mamba
```

And run the above `mamba env create...` command again.

### Step 4: Activating the environment and start executing the pipeline!
Activate the environment using
```
conda activate scnomehic
```

Your working directory may change, so `cd` back into it
```
cd <path/to/folder>
```

As this pipeline is built in Snakemake, it is with Snakemake commands, the most useful of which will be detailed in the next portion of this readme. To get a full, comprehensive list and usage of snakemake commands, you can run
```
snakemake --help
```
Or read the official snakemake documentation: https://snakemake.readthedocs.io/

To deactivate the environment (not now, since we will run the pipeline), run
```
conda deactivate
```

## Download reference genome of choice, and generate Bisulfite genome from reference genome
(Unfinished part of readme here, to be filled with Bismark's steps)

## Running the pipeline with Snakemake commands
Snakemake is a workflow management tool that links outputs to inputs from one step to another in processes called "rules", in order to work **backwards** from a desired output to create a file pipeline that leads to it. For a more detailed description of how Snakemake operates, read the documentation (https://snakemake.readthedocs.io/). Here, we will only provide a brief list of commands useful to running the pipeline

### Setting up the Config file
Open the file `configs/config.yaml` within the working directory with a text editor of your choice. If you wish to continue working from the terminal, run
```
nano configs/config.yaml
```

Here, **these arguments** MUST be modified to proper file locations/paths for the pipeline to run:

- workdir: change to your workflow folder
- reference: change to your reference file PREFIX before ".fa"
- variants: .vcf file corresponding to reference genome of choice
- restriction_sites: DpnII file
- picard: path to your picard folder
- bisulfitehic: path to bisulfitehic folder
- bistools: path to bistools folder
- scripts: path to scripts folder
- data: path to folder with raw FASTQ files
- fileindex: path to file with single cell indices

Here is example of a the fleshed out arguments within the config file
```
# working directory
workdir: "/eric/projects/scnomehic_opt"

# reference genome
reference: "reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set"
variants: "reference/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
restriction_sites: "reference/hg38/hg38_DpnII.txt"

# software
picard: "/eric/software/picard/picard.jar"
bisulfitehic: "software/bisulfitehic"
bistools: "software/Bis-tools"
scripts: "scripts"
```

The further config file arguments are for pipeline step parameters, and can be edited at one's own discretion

### Running the Snakemake pipeline locally
To perform a dry-run on the entire pipeline (without running processes, see if the pipeline is able to process all steps properly with current code), run:
```
snakemake -np
```

If this is successful, you can simply run the command
```
snakemake
```

To run the pipeline locally. Essentially, the command `snakemake` reads from the file `Snakefile`, which executes for the output within the `rule all:` at the top of the file (take a look to understand the mechanisms of snakemake). 

### Running the Snakemake pipeline on a cluster
Please see `profiles/slurm/config.yaml` and edit the arguments under `default-resources` to change the SLURM partition, account, wall time, and max number of jobs submittable to the SLURM. Once that is done, you can run:
```
snakemake --profile profiles/slurm
```

To execute the pipeline on your cluster

### Tracing errors
If errors arise, read the terminal output for said errors, or the most recent dated snakemake log in `.snakemake/logs/<YYYY-DD-MM...snakemake.log>` to trace error messages or to be redirected to job-specific log files within the folder `logs`

Happy snakemaking!
