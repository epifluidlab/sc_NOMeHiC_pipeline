import yaml
import os
import glob

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

SAMPLES = [f"{prefix}.{index}" for prefix in FASTQ_PREFIXES for index in INDICES]

rule all:
    input:
        "10.qc_analysis/qc_summary.txt"

rule getSCprefixes:
    output:
        sc_prefixes = "sc_prefixes.txt"
    run:
        with open('sc_prefixes.txt', 'w') as f:
            for cell in SAMPLES:
                f.write(f"{cell}\n")

# maybe add hic qc step before here

rule qc:
    input:
        sc_prefixes = "sc_prefixes.txt"
    output:
        qc_summary = "10.qc_analysis/qc_summary.txt",
        qc_summary_figures = "10.qc_analysis/qc_summary_figures.pdf",
        utilized_cell = "10.qc_analysis/utilized_cell.txt"
    params:
        scripts = config["scripts"],
        fragnum_threshold = config["fragnum_threshold"],
        noncpg_threshold = config["noncpg_threshold"],
        GCH_threshold = config["GCH_threshold"],
        HCG_threshold = config["HCG_threshold"],
        transrate_threshold = config["transrate_threshold"],
        cis1kbnum_threshold = config["cis1kbnum_threshold"]
    shell:
        """
        python {params.scripts}/qc_analysis.py \
        --prefix_file {input.sc_prefixes} \
        --output_file {output.qc_summary} \
        --summary_pdf {output.qc_summary_figures} \
        --utilized_cell_output {output.utilized_cell} \
        --fragnum_threshold {params.fragnum_threshold} \
        --noncpg_threshold {params.noncpg_threshold} \
        --GCH_threshold {params.GCH_threshold} \
        --HCG_threshold {params.HCG_threshold} \
        --transrate_threshold {params.transrate_threshold} \
        --cis1kbnum_threshold {params.cis1kbnum_threshold}
        """