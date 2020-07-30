import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb

import preprocess
import visualize
import epochs
import haplotypes

os.chdir('..')

# snv_data_dir = 'data/per_chr_singleton'
snv_data_file = 'data/bxd_singletons_reformatted.csv'
meta_data = 'data/strain_summary.csv'
ht_data = 'data/hmm_haplotypes'

results_df = 'out/dfs/'
results_figs = 'out/figs/'

# raw_summary, filtered_summary, muts_by_strains = preprocess.main(snv_data_file)

# raw_summary.to_csv(results_df + 'raw_singletons')
# filtered_summary.to_csv(results_df + 'filtered_singletons')
# muts_by_strains.to_csv(results_df + 'mutations_by_strains')

raw_df = pd.read_csv(results_df+'raw_singletons', index_col=0)
filtered_df = pd.read_csv(results_df+'filtered_singletons', index_col=0)
mut_strain_df = pd.read_csv(results_df+'mutations_by_strains', index_col=[0, 1, 2, 3], header=[0, 1])

# visualize.epoch_bar_charts(mut_strain_df)

# snv_df, bl_df, dba_df = visualize.mutation_spectrum_barchartsv1(mutations_strains,
# show=True, save=False, results_dir=results_figs)

# snv_frac_per_strain, snv_frac_strain_avg, ht_snv_frac_strain_avg, snv_tot_frac, ht_snv_tot_frac = \
#    visualize.mutation_spectrum_barcharts(mutation_strain_df, show=False, save=False, results_dir=results_figs)

meta_data, gen_seq = epochs.main(meta_data, mut_strain_df)

#haplotypes.load_ht_data(ht_data)
