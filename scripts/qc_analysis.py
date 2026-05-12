import argparse
from pandas.errors import EmptyDataError
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm

def count_rows(filename):
    with open(filename, 'r') as file:
        row_count = sum(1 for _ in file)
    return row_count

# Function to plot distribution of each single ones
def plot_hist_by_condition(df, column, title, xlabel, pdf_pages, bins=50):
    # Determine the common bins for all conditions
    all_data = df[column]
    min_value, max_value = all_data.min(), all_data.max()
    bins = np.linspace(min_value, max_value, bins+1)  # 100 bins

    # Colors for each condition
    colors = {'positive': 'red', 'negative': 'blue', 'sample': 'green'}

    # Plotting
    plt.figure(figsize=(10, 6))
    maxcounts = 0

    for category in df['condition'].unique():
        subset = df[df['condition'] == category][column]
        
        # Calculate the histogram data for this subset
        counts, _ = np.histogram(subset, bins=bins)
        maxcounts = max(maxcounts, max(counts))
        
        # Calculate the width for each bar
        width = np.diff(bins)[0]
        
        # Plot the bars
        plt.bar(bins[:-1], counts, width=width, alpha=0.6, label=f'Condition {category}', color=colors.get(category, 'gray'), align='edge')

    # Plotted a vertical dotted line for the location of bottom 10% sample value location and top 10% sample value location
    numbers_array = np.array(all_data)
    sorted_numbers = np.sort(numbers_array)
    percentile_index_top = int(len(sorted_numbers) * 0.9)
    percentile_index_bottom = int(len(sorted_numbers) * 0.1)
    top_10_percent_value_top = sorted_numbers[percentile_index_top]
    top_10_percent_value_bottom = sorted_numbers[percentile_index_bottom]

    plt.vlines(top_10_percent_value_top, 0, 200, linestyles='dotted', colors='green', label='Top 10%')
    plt.vlines(top_10_percent_value_bottom, 0, 200, linestyles='dotted', colors='red', label='Bottom 10%')
    
    plt.ylim(0, maxcounts + 0.1 * maxcounts)
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    plt.title(title)

    # Save the current figure to the PDF
    pdf_pages.savefig()

    # Close the current figure
    plt.close()

def main(args):
    # Create list of cell names from the prefix file passed via command-line argument
    with open(args.prefix_file, 'r') as f:
        prefixes = f.read().rstrip().split('\n')
        indexes = [x.split('_')[0] + '_' + x.split('_')[1] + '.' + x.split('.')[1] for x in prefixes]

    qcdf = pd.DataFrame(index=indexes)

    noncpg = []
    mapq30rate = []
    transnum = []
    cis1kbnum = []
    transrate = []
    cis1kbrate = []
    endo = []
    exo =  []
    rownumGCH = []
    rownumHCG = []
    uniquemap = []
    fragnum = []
    control = []

    # Initialize tqdm progress bar for processing prefixes
    for i in tqdm(range(len(prefixes)), desc="Processing prefixes", unit="prefix"):
        # Assign condition based on provided control arguments, otherwise default to 'sample'
        print(indexes[i])
        print(args.negative_controls)
        if args.positive_controls and indexes[i] in args.positive_controls:
            control.append('positive')
        elif args.negative_controls and indexes[i] in args.negative_controls:
            control.append('negative')
        else:
            control.append('sample')

        filename = f'07.bistools_snakemake/qc/{prefixes[i]}/{prefixes[i]}.calmd.trinuc_methy.chr21.txt'
        try:
            df = pd.read_csv(filename, sep='\t', header=None)
            noncpg.append(float(df[df[0] == 'ACT:'][2].values[0].split('%')[0]))
            endo.append(float(df[df[0] == 'ACG:'][2].values[0].split('%')[0]))
            exo.append(float(df[df[0] == 'GCT:'][2].values[0].split('%')[0]))
        except FileNotFoundError:
            print(filename + ' not found')
            noncpg.append(0)
            endo.append(0)
            exo.append(0)
        
        filename = f'04.alignment_snakemake/{prefixes[i]}.summary.txt'
        try:
            df = pd.read_csv(filename, sep='\t', header=None)
            uniquemap.append(float(df[df[0] == 'UniqMapped:'][2].values[0].split('%')[0]))
            mapq30rate.append(float(df[df[0] == 'UniqMappedMapQ30:'][2].values[0].split('%')[0]))
            fragnum.append(int(df[df[0] == 'TotalFragments:'][1].values[0]))
            cis1kbnum.append(int(df[df[0] == 'UniqMappedMapQ30NoPcrCisMore1kb:'][1].values[0]))
            transnum.append(int(df[df[0] == 'UniqMappedMapQ30NoPcrTrans:'][1].values[0]))
            cis1kbrate.append(float(df[df[0] == 'UniqMappedMapQ30NoPcrCisMore1kb:'][2].values[0].split('%')[0]))
            transrate.append(float(df[df[0] == 'UniqMappedMapQ30NoPcrTrans:'][2].values[0].split('%')[0]))
        except EmptyDataError:
            print(filename+' not found')
            mapq30rate.append(0)
            cis1kbnum.append(0)
            transnum.append(0)
            cis1kbrate.append(0)
            transrate.append(0)
            fragnum.append(0)
            uniquemap.append(0)

        try:
            filename = f'07.bistools_snakemake/methylation/{prefixes[i]}/{prefixes[i]}.calmd.cytosine.filtered.sort.GCH.6plus2.bed'
            rownumGCH.append(count_rows(filename)-1)
        except FileNotFoundError:
            print(filename+' not found')
            rownumGCH.append(0)
        try:
            filename = f'07.bistools_snakemake/methylation/{prefixes[i]}/{prefixes[i]}.calmd.cytosine.filtered.sort.HCG.6plus2.bed'
            rownumHCG.append(count_rows(filename)-1)
        except FileNotFoundError:
            print(filename+' not found')
            rownumHCG.append(0)

    qcdf['prefix'] = prefixes
    qcdf['fragnum'] = fragnum
    qcdf['uniquemap'] = uniquemap
    qcdf['mapq30rate'] = mapq30rate
    qcdf['transnum'] = transnum
    qcdf['cis1kbnum'] = cis1kbnum
    qcdf['transrate'] = transrate
    qcdf['cis1kbrate'] = cis1kbrate
    qcdf['noncpg'] = noncpg
    qcdf['endo'] = endo
    qcdf['exo'] = exo
    qcdf['GCH'] = rownumGCH
    qcdf['HCG'] = rownumHCG
    qcdf['condition'] = control

    # Save the summary to the output file passed as an argument
    qcdf.to_csv(args.output_file, sep='\t')

    # Create a PDF to store all figures
    with PdfPages(args.summary_pdf) as pdf_pages:
        plot_hist_by_condition(qcdf, 'fragnum', 'Distribution of raw fragment number', '# of fragments', pdf_pages)
        plot_hist_by_condition(qcdf, 'uniquemap', 'Distribution of unique mapping rate', 'Mapping rate(%)', pdf_pages)
        plot_hist_by_condition(qcdf, 'mapq30rate', 'Distribution of rate of fragments with mapq>30', 'Rate of fragments with MAPQ>30 (%)', pdf_pages)
        plot_hist_by_condition(qcdf, 'transnum', 'Distribution of number of trans interaction for Hi-C', '# of trans interaction', pdf_pages)
        plot_hist_by_condition(qcdf, 'cis1kbnum', 'Distribution of number of cis interaction with 1kb away for Hi-C', '# of cis interaction > 1kb distance', pdf_pages)
        plot_hist_by_condition(qcdf, 'transrate', 'Distribution of rate of trans interaction for Hi-C', 'Rate of trans interaction', pdf_pages)
        plot_hist_by_condition(qcdf, 'cis1kbrate', 'Distribution of rate of cis interaction with 1kb away for Hi-C', 'Rate of cis interaction > 1kb distance', pdf_pages)
        plot_hist_by_condition(qcdf, 'noncpg', 'Distribution of non-CpG methylation rate', 'Non-CpG methylation rate(%)', pdf_pages)
        plot_hist_by_condition(qcdf, 'endo', 'Distribution of endogenous methylation rate', 'Endogenous methylation rate(%)', pdf_pages)
        plot_hist_by_condition(qcdf, 'exo', 'Distribution of exogenous methylation rate', 'Exogenous methylation rate(%)', pdf_pages)
        plot_hist_by_condition(qcdf, 'GCH', 'Distribution of number of GCH sites', 'Number of GCH sites', pdf_pages)
        plot_hist_by_condition(qcdf, 'HCG', 'Distribution of number of HCG sites', 'Number of HCG sites', pdf_pages)

    # Filter and get sample names
    qcdf_filtered = qcdf[(qcdf['condition'] == 'sample') & 
                         (qcdf['fragnum'] > args.fragnum_threshold) & 
                         (qcdf['noncpg'] < args.noncpg_threshold) & 
                         (qcdf['GCH'] > args.GCH_threshold) & 
                         (qcdf['HCG'] > args.HCG_threshold) & 
                         (qcdf['transrate'] > args.transrate_threshold) & 
                         (qcdf['cis1kbnum'] > args.cis1kbnum_threshold)]
    
    utilized_cell = qcdf_filtered['prefix'].tolist()

    # Save the utilized cell list to the output file passed via argument
    with open(args.utilized_cell_output, 'w') as file: 
        for item in utilized_cell:
            file.write("%s\n" % item)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process QC data and plot distributions")
    parser.add_argument('--prefix_file', required=True, help='File containing prefixes')
    parser.add_argument('--output_file', required=True, help='Output file for QC summary')
    parser.add_argument('--summary_pdf', required=True, help='Output PDF file for summary figures')
    parser.add_argument('--utilized_cell_output', required=True, help='Output file for utilized cells')
    parser.add_argument('--positive_controls', nargs='+', help='List of positive control prefixes (optional)')
    parser.add_argument('--negative_controls', nargs='+', help='List of negative control prefixes (optional)')
    parser.add_argument('--fragnum_threshold', type=int, default=800000, help='Fragment number threshold')
    parser.add_argument('--noncpg_threshold', type=float, default=5, help='Non-CpG methylation rate threshold')
    parser.add_argument('--GCH_threshold', type=int, default=1000000, help='GCH site count threshold')
    parser.add_argument('--HCG_threshold', type=int, default=150000, help='HCG site count threshold')
    parser.add_argument('--transrate_threshold', type=float, default=0.7, help='Trans interaction rate threshold')
    parser.add_argument('--cis1kbnum_threshold', type=int, default=2000, help='Cis interaction > 1kb distance threshold')

    # Parse the arguments and pass them to the main function
    args = parser.parse_args()
    main(args)
