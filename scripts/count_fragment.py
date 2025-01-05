import pandas as pd
import numpy as np
import argparse

def read_chrom_sizes(filepath):
    chromsize = pd.read_csv(filepath, sep='\t', header=None, names=['chr', 'size'])
    chromsize = chromsize[chromsize['chr'].isin([f'chr{i}' for i in range(1, 24)] + ['chrX', 'chrY'])]
    return chromsize

def read_bed(filepath, binsize=1000000):
    bed_df = pd.read_csv(filepath, sep='\t', header=None, names=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2'])
    bed_df = bed_df[bed_df['chrom1'].isin(['chr'+str(i) for i in range(1, 24)]+['chrX','chrY'])]
    bed_df['midpoint'] = (bed_df[['start1', 'start2', 'end1', 'end2']].min(axis=1) + bed_df[['start1', 'start2', 'end1', 'end2']].max(axis=1)) // 2
    bed_df['interval'] = pd.IntervalIndex.from_arrays(bed_df['midpoint']//binsize*binsize+1, (bed_df['midpoint']//binsize+1)*binsize, closed='right')
    return bed_df

def generate_bindf(chromsize, binsize=1000000):
    bindf = pd.DataFrame({
        'chr': np.repeat(chromsize['chr'].values, (chromsize['size'] // binsize + 1)),
        'start': np.hstack([np.arange(0, size, binsize) + 1 for size in chromsize['size']]),
        'end': np.hstack([np.arange(binsize, size + binsize, binsize) for size in chromsize['size']])
    })
    bindf['interval'] = pd.IntervalIndex.from_arrays(bindf['start'], bindf['end'], closed='right')
    return bindf

def count_bin(bindf, bed_df):
    # Merge dataframes on intervals
    result = pd.merge(bed_df, bindf, left_on=['chrom1', 'interval'], right_on=['chr', 'interval'])
    count_series = result.groupby(['chr', 'start', 'end']).size()
    
    return count_series

def read_prefixes(prefix_file):
    with open(prefix_file, 'r') as file:
        prefixes = [line.strip() for line in file.readlines()]
    return prefixes

def main():
    parser = argparse.ArgumentParser(description='Count Fragments in Bins')
    parser.add_argument('--chrom_size_filepath', required=True, help='File path to chromosome size file')
    parser.add_argument('--prefix_file', required=True, help='File path to the file with BED file prefixes')
    parser.add_argument('--bed_directory', required=True, help='Directory where BED files are stored')
    parser.add_argument('--binsize', type=int, default=1000000, help='Size of bins in base pairs')
    parser.add_argument('--outfile', required=True, help='Output file for the resulting CSV')
    
    args = parser.parse_args()
    
    chromsize = read_chrom_sizes(args.chrom_size_filepath)    
    bindf = generate_bindf(chromsize,binsize=args.binsize)
    result_df = bindf.set_index(['chr', 'start', 'end'])
    
    bed_files = read_prefixes(args.prefix_file)
    
    for prefix in bed_files:
        bed_df = read_bed(args.bed_directory + '/' + prefix + '.bed', binsize=args.binsize)
        count_series = count_bin(bindf, bed_df)
        result_df[prefix] = count_series

    result_df.fillna(0, inplace=True)
    result_df.drop(columns=['interval'], inplace=True)
    result_df.to_csv(args.outfile, sep='\t')

if __name__ == "__main__":
    main()
