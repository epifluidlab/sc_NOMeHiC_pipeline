"""Compute per-bin methylation level for one cell, one context, one bin size.

Reads the per-base BisSNP output
  07.bistools_snakemake/methylation/{prefix}/{prefix}.cyt.filtered.sort.{context}.6plus2.bed
and writes a per-bin TSV
  {outfolder}/{prefix}.{context}.{bin_size}bp_methylation.txt

Each row is one genome-wide bin (chrN:start-end) with the weighted methylation
level (sum(rate/100 * num_reads) / sum(num_reads)) over CpG/GCH sites that fall
within it.
"""

import pandas as pd
import numpy as np
import argparse
import os
import sys


AUTOSOMES_AND_SEX = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']


def read_chrom_sizes(filepath):
    chromsize = pd.read_csv(filepath, sep='\t', header=None, names=['chr', 'size'])
    chromsize = chromsize[chromsize['chr'].isin(AUTOSOMES_AND_SEX)]
    return chromsize


def generate_methyldf(chromsize, bin_size):
    rows = []
    for chrom, size in chromsize[['chr', 'size']].itertuples(index=False):
        starts = np.arange(0, size, bin_size, dtype=np.int64)
        ends = np.minimum(starts + bin_size, size)
        for s, e in zip(starts, ends):
            rows.append((f'{chrom}:{s + 1}-{e}', chrom, int(s + 1), int(e)))
    methyldf = pd.DataFrame(rows, columns=['bin_id', 'chr', 'start', 'end']).set_index('bin_id')
    return methyldf


def read_bed(filename):
    """Read a `.cyt.filtered.sort.{context}.6plus2.bed` once. Filter to autosomes + XY."""
    bed = pd.read_csv(
        filename, sep='\t', header=None, skiprows=[0],
        names=['chr', 'start', 'end', 'name', 'ha', 'direction', 'rate', 'num']
    )
    return bed[bed['chr'].isin(AUTOSOMES_AND_SEX)].reset_index(drop=True)


def compute_methylation_levels(methyldf, bed):
    """Aggregate weighted methylation per bin in `methyldf`. Reads `bed` once,
    groups its rows by chromosome, then for each bin slices the chrom group by
    [start, end] coordinate range. O(n_bins + n_sites) per chromosome rather
    than O(n_bins * n_sites)."""
    methylevel = np.full(len(methyldf), np.nan, dtype=np.float64)
    by_chrom = {chrom: g[['start', 'num', 'rate']].sort_values('start').reset_index(drop=True)
                for chrom, g in bed.groupby('chr')}
    for idx, (_, chrom, start, end) in enumerate(methyldf.reset_index().itertuples(index=False)):
        g = by_chrom.get(chrom)
        if g is None or g.empty:
            continue
        # Bed rows are CpG/GCH sites at single positions. The BisSNP 6plus2 BED
        # stores them as half-open intervals where `start` is 0-based; methyldf
        # bins are 1-based inclusive. A site is "in bin" iff its `start` (BED
        # 0-based) >= bin.start - 1 AND `end` (BED 0-based exclusive) <= bin.end.
        sub = g[(g['start'] >= start - 1) & (g['start'] < end)]
        n = sub['num'].sum()
        if n > 0:
            methylevel[idx] = (sub['rate'] / 100.0 * sub['num']).sum() / n
    methyldf = methyldf.copy()
    methyldf['methylation_level'] = methylevel
    return methyldf


def main():
    parser = argparse.ArgumentParser(description='Compute per-bin methylation level for one cell + context.')
    parser.add_argument('--chrom_size_filepath', required=True, help='Path to chromosome size file (chr\\tsize)')
    parser.add_argument('--sample_prefix', required=True, help='Cell prefix (e.g. batch5_scD_1.CGATGT)')
    parser.add_argument('--context', required=True, choices=['GCH', 'HCG'],
                        help='Cytosine context to aggregate (GCH or HCG)')
    parser.add_argument('--outfolder', required=True, help='Output folder for the per-bin TSV')
    parser.add_argument('--bin_size', type=int, default=1000000,
                        help='Bin size in bp for methylation calculation (default 1Mb)')
    parser.add_argument('--bed_root', default='07.bistools_snakemake/methylation',
                        help='Root dir for the per-cell BisSNP output')

    args = parser.parse_args()

    bed_path = os.path.join(
        args.bed_root,
        args.sample_prefix,
        f'{args.sample_prefix}.cyt.filtered.sort.{args.context}.6plus2.bed',
    )
    outfile_path = os.path.join(
        args.outfolder,
        f'{args.sample_prefix}.{args.context}.{args.bin_size}bp_methylation.txt',
    )

    os.makedirs(args.outfolder, exist_ok=True)
    if os.path.exists(outfile_path):
        # Snakemake will recompute mtimes; respect already-present output.
        sys.exit()

    if not os.path.exists(bed_path):
        sys.exit(f"ERROR: bed not found: {bed_path}")

    chromsize = read_chrom_sizes(args.chrom_size_filepath)
    methyldf = generate_methyldf(chromsize, args.bin_size)
    bed = read_bed(bed_path)
    methyldf = compute_methylation_levels(methyldf, bed)
    methyldf.to_csv(outfile_path, sep='\t')


if __name__ == '__main__':
    main()
