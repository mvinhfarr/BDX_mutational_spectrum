import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb

import preprocess
import visualize
import meta
import haplotypes
import per_chrom

# snv_data_file = 'data/bxd_singletons_reformatted.csv'
meta_data_file = 'data/strain_summary.csv'
ht_data_dir = 'data/hmm_haplotypes'

main_dir = 'out/dfs/'
# het_df = 'out/het/dfs/'
het_figs = 'out/het/figs/all_hets/'

raw_df = pd.read_csv(main_dir + 'raw_singletons', index_col=0)

filtered_summary = raw_df[~raw_df['kmer'].str.contains('indel')]
filtered_summary = filtered_summary[filtered_summary['gt'] == 1]
filtered_df = filtered_summary.drop('BXD087/RwwJ', axis='index', inplace=False)

mut_strain_df = preprocess.mutations_by_strains_df(filtered_df)

# load meta data; filter; extract generation at sequence and epoch
meta_df, gens_df, epoch_df = meta.main(meta_data_file, mut_strain_df)

# snv_df, bl_df, dba_df = visualize.mutation_spectrum_barcharts_simple(mut_strain_df, show=False, save=True,
#                                                                      results_dir=het_figs)
#
# snv_frac_per_strain, snv_frac_strain_avg, ht_snv_frac_strain_avg, snv_tot_frac, ht_snv_tot_frac = \
#     visualize.mutation_spectrum_barcharts(mut_strain_df, show=False, save=True, results_dir=het_figs)

muts_per_strain, muts_per_strain_per_gen = visualize.strain_distrb(mut_strain_df, epoch_df, gens_df,
                                                                   show=False, save=False, results_dir=het_figs)

# visualize.other_bar_charts(mut_strain_df, show=False)

# ht_dict, d2_frac_df, muts_per_chrom = haplotypes.main(ht_data_dir, filtered_df)
#
# muts_per_chrom_per_gen = visualize.mutation_rate(muts_per_chrom, epoch_df, gens_df,
#                                                  show=False, save=True, results_dir=het_figs)

ht_ratio = visualize.mutation_spectrum_heatmap(mut_strain_df, show=False, save=False, results_dir=het_figs)

# chrom_spect_dict = per_chrom.per_chrom_mut_spectrum(filtered_df)
# chrom_ratios, chrom_snv_ratios = per_chrom.plot(chrom_spect_dict, show=False, save=True, save_dir=het_figs)

notable = muts_per_strain[muts_per_strain > 500]
others = muts_per_strain[muts_per_strain <= 500]

ab_scores = filtered_df.ab
notable_ab_scores = ab_scores.loc[notable.index]
other_ab_scores = ab_scores.loc[others.index]
