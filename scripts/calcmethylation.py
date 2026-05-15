"""Compute per-bin methylation level for one cell, one context.

Reads the per-base BisSNP output
  {bed_root}/{prefix}/{prefix}.cyt.filtered.sort.{context}.6plus2.bed
and writes a per-bin TSV to --outfile.

Bin definitions come from one of two sources:
  --bin_size N  (with --chrom_size_filepath)
      Fixed-size tiling across autosomes + chrX + chrY.
  --bin_bed FILE
      Custom BED of regions. The first three cols are chr/start/end (BED
      0-based half-open). A fourth col (if present) is used as the bin name
      (e.g. gene_tss BEDs naming each region by gene); otherwise the bin_id
      is the genomic range `chr:start+1-end`. Header rows (where col 2 is
      not an integer) and `#`-comment lines are skipped.

Output schema (TSV with header):
  bin_id  chr  start  end  methylated_reads  total_reads  methylation_level
"""

import pandas as pd
import numpy as np
import argparse
import gzip
import os
import sys


def _open_text(path):
    """Open `path` for text reading, transparently handling .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


AUTOSOMES_AND_SEX = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']


def read_chrom_sizes(filepath):
    chromsize = pd.read_csv(filepath, sep='\t', header=None, names=['chr', 'size'])
    chromsize = chromsize[chromsize['chr'].isin(AUTOSOMES_AND_SEX)]
    return chromsize


def methyldf_from_binsize(chromsize, bin_size):
    rows = []
    for chrom, size in chromsize[['chr', 'size']].itertuples(index=False):
        starts = np.arange(0, size, bin_size, dtype=np.int64)
        ends = np.minimum(starts + bin_size, size)
        for s, e in zip(starts, ends):
            rows.append((f'{chrom}:{s + 1}-{e}', chrom, int(s + 1), int(e)))
    return pd.DataFrame(rows, columns=['bin_id', 'chr', 'start', 'end']).set_index('bin_id')


def methyldf_from_bin_bed(path):
    """Parse a BED file into bins. Tolerates `#` comments, header rows, and
    transparent gzip (`.bed.gz`)."""
    rows = []
    seen_ids = {}
    with _open_text(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#') or line.startswith('track'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                # header row (e.g. `chr\tstart\tend\tgene_name`)
                continue
            chrom = parts[0]
            name = parts[3] if len(parts) >= 4 and parts[3] not in ('', '.') else f'{chrom}:{start + 1}-{end}'
            # Ensure bin_id is unique — append :N suffix on collision
            key = name
            n = seen_ids.get(key, 0)
            if n:
                name = f'{key}#{n}'
            seen_ids[key] = n + 1
            rows.append((name, chrom, start + 1, end))
    if not rows:
        raise ValueError(f"No usable bins parsed from {path}")
    return pd.DataFrame(rows, columns=['bin_id', 'chr', 'start', 'end']).set_index('bin_id')


def read_bed(filename):
    """Read the per-base BisSNP 6plus2 BED once. Filter to autosomes + XY."""
    bed = pd.read_csv(
        filename, sep='\t', header=None, skiprows=[0],
        names=['chr', 'start', 'end', 'name', 'ha', 'direction', 'rate', 'num']
    )
    return bed[bed['chr'].isin(AUTOSOMES_AND_SEX)].reset_index(drop=True)


def compute_methylation_levels(methyldf, bed):
    """Aggregate per bin in `methyldf`. Adds methylated_reads, total_reads,
    methylation_level columns. O(n_bins + n_sites) per chromosome."""
    n_bins = len(methyldf)
    methylated_reads = np.zeros(n_bins, dtype=np.int64)
    total_reads = np.zeros(n_bins, dtype=np.int64)
    methylation_level = np.full(n_bins, np.nan, dtype=np.float64)
    by_chrom = {chrom: g[['start', 'num', 'rate']].sort_values('start').reset_index(drop=True)
                for chrom, g in bed.groupby('chr')}
    for idx, (_, chrom, start, end) in enumerate(methyldf.reset_index().itertuples(index=False)):
        g = by_chrom.get(chrom)
        if g is None or g.empty:
            continue
        # methyldf bins are 1-based inclusive; bed start is 0-based.
        sub = g[(g['start'] >= start - 1) & (g['start'] < end)]
        n = int(sub['num'].sum())
        if n > 0:
            meth_float = float((sub['rate'] / 100.0 * sub['num']).sum())
            total_reads[idx] = n
            methylated_reads[idx] = int(round(meth_float))
            methylation_level[idx] = meth_float / n
    out = methyldf.copy()
    out['methylated_reads'] = methylated_reads
    out['total_reads'] = total_reads
    out['methylation_level'] = methylation_level
    return out


def main():
    parser = argparse.ArgumentParser(description='Compute per-bin methylation for one cell + context.')
    parser.add_argument('--outfile', required=True, help='Output TSV path')
    parser.add_argument('--sample_prefix', help='Cell prefix (e.g. batch5_scD_1.CGATGT); used '
                        'to construct the BED path when --bed_path is omitted')
    parser.add_argument('--context', choices=['GCH', 'HCG'],
                        help='Cytosine context to aggregate; required for path construction '
                        'when --bed_path is omitted')
    parser.add_argument('--bed_path',
                        help='Explicit path to the per-base BisSNP 6plus2 BED (plain or .gz). '
                        'Overrides the --bed_root/--sample_prefix/--context construction.')
    parser.add_argument('--bed_root', default='07.bistools_snakemake/methylation',
                        help='Root dir for the per-cell BisSNP 6plus2 BED (used when '
                        '--bed_path is not given)')
    # Bin definition (one required)
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument('--bin_size', type=int, help='Fixed bin size in bp (requires --chrom_size_filepath)')
    grp.add_argument('--bin_bed', help='Path to a BED file defining custom regions')
    parser.add_argument('--chrom_size_filepath', help='Required when --bin_size is used')

    args = parser.parse_args()

    if args.bin_size is not None and not args.chrom_size_filepath:
        parser.error("--bin_size requires --chrom_size_filepath")

    # Resolve the site BED path. Explicit --bed_path wins; otherwise build it
    # from prefix + context and try .bed then .bed.gz.
    if args.bed_path:
        site_bed_path = args.bed_path
        if not os.path.exists(site_bed_path):
            sys.exit(f"ERROR: site bed not found: {site_bed_path}")
    else:
        if not args.sample_prefix or not args.context:
            parser.error("Either --bed_path OR (--sample_prefix AND --context) must be given.")
        plain = os.path.join(
            args.bed_root, args.sample_prefix,
            f'{args.sample_prefix}.cyt.filtered.sort.{args.context}.6plus2.bed',
        )
        if os.path.exists(plain):
            site_bed_path = plain
        elif os.path.exists(plain + ".gz"):
            site_bed_path = plain + ".gz"
        else:
            sys.exit(f"ERROR: site bed not found: {plain} (also tried {plain}.gz)")

    if args.bin_size is not None:
        chromsize = read_chrom_sizes(args.chrom_size_filepath)
        methyldf = methyldf_from_binsize(chromsize, args.bin_size)
    else:
        methyldf = methyldf_from_bin_bed(args.bin_bed)

    bed = read_bed(site_bed_path)
    out = compute_methylation_levels(methyldf, bed)

    os.makedirs(os.path.dirname(args.outfile) or '.', exist_ok=True)
    out.to_csv(args.outfile, sep='\t')


if __name__ == '__main__':
    main()
