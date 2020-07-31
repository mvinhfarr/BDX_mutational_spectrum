import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb

import preprocess
import visualize
import meta
import haplotypes


snv_data_dir = 'data/per_chr_singleton'
snv_data_file = 'data/bxd_singletons_reformatted.csv'
meta_data = 'data/strain_summary.csv'
ht_data = 'data/hmm_haplotypes'

results_df = 'out/dfs/'
results_figs = 'out/figs/'

# load raw data; filter; format
raw_summary, filtered_summary, muts_by_strains = preprocess.main(snv_data_file)

raw_summary.to_csv(results_df + 'raw_singletons')
filtered_summary.to_csv(results_df + 'filtered_singletons')
muts_by_strains.to_csv(results_df + 'mutations_by_strains')

raw_df = pd.read_csv(results_df+'raw_singletons', index_col=0)
filtered_df = pd.read_csv(results_df+'filtered_singletons', index_col=0)
mut_strain_df = pd.read_csv(results_df+'mutations_by_strains', index_col=[0, 1, 2, 3], header=[0, 1])

# load meta data; extract generation at sequence and epoch
meta_data, gen_seq, epochs = meta.main(meta_data, mut_strain_df)

meta_data.to_csv(results_df+'meta_data')
gen_seq.to_csv(results_df+'gen_at_seq')
epochs.to_csv(results_df+'epochs')

meta_df = pd.read_csv(results_df+'meta_data', index_col=0, header=0)
generations_df = pd.read_csv(results_df+'gen_at_seq', index_col=0, header=0)
epoch_df = pd.read_csv(results_df+'epochs', index_col=0, header=0)

# visualize.epoch_bar_charts(mut_strain_df)

# snv_df, bl_df, dba_df = visualize.mutation_spectrum_barchartsv1(mutations_strains,
# show=True, save=False, results_dir=results_figs)

# snv_frac_per_strain, snv_frac_strain_avg, ht_snv_frac_strain_avg, snv_tot_frac, ht_snv_tot_frac = \
#    visualize.mutation_spectrum_barcharts(mutation_strain_df, show=False, save=False, results_dir=results_figs)

# haplotypes.main(ht_data, filtered_df)
