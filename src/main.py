import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
from scipy import stats

import preprocess
import visualize
import meta
import haplotypes
import per_chrom
import stat_tests

# snv_data_dir = 'data/per_chr_singleton/'
snv_data_file = 'data/bxd_singletons_reformatted.csv'
meta_data_file = 'data/strain_summary.csv'
ht_data_dir = 'data/hmm_haplotypes'
ht_dict_dir = 'ht_dict/'

results_df = 'out/dfs/'
results_figs = 'out/figs/'

# # load raw data; filter; format; mutation spectra counts
# raw_summary, filtered_summary, muts_by_strains = preprocess.main(snv_data_file)
#
# raw_summary.to_csv(results_df + 'raw_singletons')
# filtered_summary.to_csv(results_df + 'filtered_singletons')
# muts_by_strains.to_csv(results_df + 'mutations_by_strains')

raw_df = pd.read_csv(results_df + 'raw_singletons', index_col=0)
filtered_df = pd.read_csv(results_df + 'filtered_singletons', index_col=0)
mut_strain_df = pd.read_csv(results_df + 'mutations_by_strains', index_col=[0, 1, 2, 3], header=0)

# # load meta data; filter; extract generation at sequence and epoch
# meta_data, gen_seq, epochs = meta.main(meta_data_file, mut_strain_df)
#
# meta_data.to_csv(results_df+'meta_data')
# gen_seq.to_csv(results_df+'gen_at_seq')
# epochs.to_csv(results_df+'epochs')

meta_df = pd.read_csv(results_df + 'meta_data', index_col=0, header=0)
gens_df = pd.read_csv(results_df + 'gen_at_seq', index_col=0, header=0)
epoch_df = pd.read_csv(results_df + 'epochs', index_col=0, header=0)

'''
# snv_df -> count and fraction of total for each the 6 snv types
# bl_df/dba_df -> count and fract of per haplotype total
snv_df, bl_df, dba_df = visualize.mutation_spectrum_barcharts_simple(mut_strain_df, show=False, save=False,
                                                                     results_dir=results_figs)

# snv_frac_per_strain -> divide mut_strain_df by per strain totals
# snv_frac_strain_avg -> collapse into snv and find average strain
# ht_snv_frac_strain_avg -> collapse into snv by ht and find average strain
# snv_tot_frac -> collapse counts across strains (same as snv_df)
# ht_snv_tot_frac -> collapse counts across strains but divide by ht (same as bl_df/dba_df)
snv_frac_per_strain, snv_frac_strain_avg, ht_snv_frac_strain_avg, snv_tot_frac, ht_snv_tot_frac = \
    visualize.mutation_spectrum_barcharts(mut_strain_df, show=False, save=False, results_dir=results_figs)

muts_per_strain, muts_per_strain_per_gen = visualize.strain_distrb(mut_strain_df, epoch_df, gens_df,
                                                                   show=False, save=False, results_dir=results_figs)

visualize.other_bar_charts(mut_strain_df, show=False)

# NOT WORKING
# visualize.epoch_bar_charts(mut_strain_df)
'''
# ht_dict, d2_frac_df, muts_per_chrom = haplotypes.main(ht_data_dir, filtered_df)


# for key, df in ht_dict.items():
#     df.to_csv(results_df + ht_dict_dir + key)
# d2_frac_df.to_csv(results_df + 'd2_frac_per_chrom')
# muts_per_chrom.to_csv(results_df + 'muts_per_chrom')

ht_strain_dict = {}
for f in os.listdir(results_df+ht_dict_dir):
    ht_strain_dict[f] = pd.read_csv(results_df+ht_dict_dir+f, index_col=0, header=0)
d2_frac_df = pd.read_csv(results_df+'d2_frac_per_chrom', index_col=0, header=0)
muts_per_chrom = pd.read_csv(results_df+'muts_per_chrom', index_col=[0, 1], header=0)

# # chr_windows = haplotypes.chrom_ht_windows(ht_strain_dict, filtered_df)
# chr_windows = haplotypes.chrom_ht_windows(ht_strain_dict, filtered_df, num_win=10)
#
# filtered_df, chr_windows = haplotypes.chrom_win_muts(filtered_df, chr_windows)
#
# filtered_df.to_csv(results_df + 'filtered_df')
# chr_windows.to_csv(results_df + 'chrom_windows')

filtered_df = pd.read_csv(results_df + 'filtered_df')
chr_windows = pd.read_csv(results_df + 'chrom_windows')

chr_windows = stat_tests.chr_windows_chi_sq_test(chr_windows)

chr_windows = chr_windows.set_index(['chrom', 'window'])

df1 = chr_windows.b6_muts_per_bp / chr_windows.d2_muts_per_bp
df1 = df1.unstack(level='chrom')

df2 = chr_windows.b6_bp / chr_windows.d2_bp
df2 = df2.unstack(level='chrom')

df3 = chr_windows.b6_muts / chr_windows.d2_muts
df3 = df3.unstack(level='chrom')

fig, ax = plt.subplots(2, 2)

sb.heatmap(df1, ax=ax[0][0], cmap='bwr', cbar=True, square=True, vmin=0.5, vmax=1.5, xticklabels=True, yticklabels=True)
sb.heatmap(df2, ax=ax[0][1], cmap='bwr', cbar=True, square=True, center=1, xticklabels=True, yticklabels=True)
sb.heatmap(df3, ax=ax[1][0], cmap='bwr', cbar=True, square=True, vmin=0.5, vmax=1.5, xticklabels=True, yticklabels=True)

# # needs to be redone with expected probability = 1
# ab_p_vals = stat_tests.ab_binomial_test(filtered_df)
# visualize.ab_scores(ab_p_vals, filtered_df)

# muts_per_chrom_per_gen = visualize.mutation_rate(muts_per_chrom, epoch_df, gens_df,
#                                                  show=False, save=True, results_dir=results_figs)

# p_vals = mutspec_stats.mut_spec_chi_sq(mut_strain_df)
# mut_fracs, ht_ratio = visualize.mutation_spectrum_heatmap(mut_strain_df, show=True, save=False,
#                                                           results_dir=results_figs, vrange=0.5)

# chrom_spect_dict = per_chrom.per_chrom_mut_spectrum(filtered_df)
# chrom_ratios, chrom_snv_ratios = per_chrom.plot(chrom_spect_dict, show=True, save=True, save_dir=results_figs)
