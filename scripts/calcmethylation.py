import pandas as pd
import numpy as np
import argparse
import os
import sys

def read_chrom_sizes(filepath):
    chromsize = pd.read_csv(filepath, sep='\t', header=None, names=['chr', 'size'])
    chromsize = chromsize[chromsize['chr'].isin([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'])]
    return chromsize

def generate_methyldf(chromsize, bin_size):
    methyldf = pd.DataFrame(columns=['chr', 'start', 'end'])
    for i in chromsize['chr'].unique().tolist():
        size = chromsize[chromsize['chr'] == i]['size'].values[0]
        start = np.arange(0, size, bin_size)
        end = np.arange(bin_size, size + bin_size, bin_size)
        for j in range(len(start)):
            methyldf.loc[i + ':' + (start[j] + 1).astype(str) + '-' + end[j].astype(str)] = [i, start[j] + 1, end[j]]
    return methyldf

def compute_methylation_levels(methyldf, prefix):
    methyllevel = []
    for i in range(len(methyldf.index.tolist())):
        chr = methyldf.iloc[i]['chr']
        start = methyldf.iloc[i]['start']
        end = methyldf.iloc[i]['end']
        filename = f'07.bistools_snakemake/methylation/{prefix}/{prefix}.cyt.filtered.sort.HCG.6plus2.bed'
        bed = pd.read_csv(filename, sep='\t', header=None, skiprows=[0], names=['chr', 'start', 'end', 'name', 'ha', 'direction', 'rate', 'num'])
        bed = bed[bed['chr'].isin([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'])]
        subbed = bed[(bed['chr'] == chr) & (bed['start'] >= start) & (bed['end'] <= end)]
        allpoints = subbed['num'].sum()
        methylpoints = (subbed['rate'] / 100 * subbed['num']).sum()
        methyllevel.append(methylpoints / allpoints if allpoints > 0 else np.nan)
    methyldf[prefix] = methyllevel
    return methyldf

def main():
    parser = argparse.ArgumentParser(description='Compute Methylation Levels')
    parser.add_argument('--chrom_size_filepath', required=True, help='File path to chromosome size file')
    parser.add_argument('--sample_prefix', required=True, help='Prefix path for sample files')
    parser.add_argument('--outfolder', required=True, help='Output folder for the resulting CSV file')
    parser.add_argument('--bin_size', type=int, default=1000000, help='Bin size for methylation calculation')

    args = parser.parse_args()

    outfile_path = os.path.join(args.outfolder, f'{args.sample_prefix}_methylation.txt')
    # Ensure the output directory exist
    os.makedirs(args.outfolder, exist_ok=True)

    if os.path.exists(outfile_path) == True:
        sys.exit()

    # Process the chromosome sizes
    chromsize = read_chrom_sizes(args.chrom_size_filepath)

    # Generate the base methylation DataFrame
    methyldf = generate_methyldf(chromsize, args.bin_size)

    # Compute methylation levels for the specified cell
    methyldf = compute_methylation_levels(methyldf, args.sample_prefix)

    # Save the DataFrame to the specified output folder
    methyldf.to_csv(outfile_path, sep='\t')

if __name__ == '__main__':
    main()
