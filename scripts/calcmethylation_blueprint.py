import pandas as pd
import numpy as np
import argparse
import os
import sys

def read_chrom_sizes(filepath='/home/jmj7858/epifluidlab/sc_nomehic/reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.chrom.sizes'):
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

def compute_methylation_levels(methyldf, prefix, merged_df, bin_size=25000):
    # Adding bin information to merged_df
    merged_df['bin'] = (merged_df['start'] // bin_size) * bin_size + 1
    merged_df['bin_end'] = merged_df['bin'] + bin_size - 1

    # Group by chromosome and bin
    grouped = merged_df.groupby(['chr', 'bin'])

    # Compute methylation levels
    methylation_data = grouped.apply(lambda x: (x['methylation_rate'] / 100 * x['read_number']).sum() / x['read_number'].sum() if x['read_number'].sum() > 0 else np.nan)

    # Reset index to merge with methyldf
    methylation_data = methylation_data.reset_index()
    methylation_data.columns = ['chr', 'start', prefix]

    # Merge with methyldf
    methyldf = methyldf.merge(methylation_data, how='left', on=['chr', 'start'])
    return methyldf

def main():
    parser = argparse.ArgumentParser(description='Compute Methylation Levels')
    parser.add_argument('--chrom_size_filepath', required=True, help='File path to chromosome size file')
    parser.add_argument('--sample_prefix', required=True, help='Prefix path for sample files')
    parser.add_argument('--outfolder', required=True, help='Output folder for the resulting CSV file')
    parser.add_argument('--bin_size', type=int, default=1000000, help='Bin size for methylation calculation')

    args = parser.parse_args()

    bin_size = args.bin_size

    outfile_path = os.path.join(args.outfolder, f'{args.sample_prefix}_methylation_hcg_{args.bin_size}.txt')
    # Ensure the output directory exists
    os.makedirs(args.outfolder, exist_ok=True)

    if os.path.exists(outfile_path):
        sys.exit()

    # Process the chromosome sizes
    chromsize = read_chrom_sizes(args.chrom_size_filepath)

    # Generate the base methylation DataFrame
    methyldf = generate_methyldf(chromsize, bin_size)
    
    depth = pd.read_csv('blueprint/' + args.sample_prefix + '_sequencing_depth.bed', sep='\t', header=None, names=['chr', 'start', 'end', 'id', 'read_number'])
    methylation = pd.read_csv('blueprint/' + args.sample_prefix + '_methylation_signal.bed', sep='\t', header=None, names=['chr', 'start', 'end', 'id', 'methylation_rate'])
    merged_df = pd.merge(depth, methylation, on=['chr', 'start', 'end', 'id'])
    merged_df = merged_df[merged_df['chr'].isin([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'])]
    
    HCG_loc_pos = pd.read_csv('/home/jmj7858/epifluidlab/sc_nomehic/analysis/HCG_motif_seqkit_all_pos.hg38.bed',sep='\t',header=None,names=['chr','start','end','seq','score','strand'])
    mask = merged_df.set_index(['chr', 'start']).index.isin(HCG_loc_pos.set_index(['chr', 'end']).index)
    merged_df_hcg = merged_df[mask]
    # calculate the percentage of HCG sites and only keep two decimal places
    print(f"Percentage of HCG sites: {merged_df_hcg.shape[0]/merged_df.shape[0]*100:.2f}%")

    # Compute methylation levels for the specified cell
    methyldf = compute_methylation_levels(methyldf, args.sample_prefix, merged_df_hcg, bin_size=bin_size)

    # Save the DataFrame to the specified output folder
    methyldf.to_csv(outfile_path, sep='\t', index=False)

if __name__ == '__main__':
    main()
