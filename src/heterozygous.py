import os
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sb

import preprocess
import visualize
import meta
import haplotypes
import per_chrom
import stat_tests

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

# ht_ratio = visualize.mutation_spectrum_heatmap(mut_strain_df, show=False, save=False, results_dir=het_figs)

# chrom_spect_dict = per_chrom.per_chrom_mut_spectrum(filtered_df)
# chrom_ratios, chrom_snv_ratios = per_chrom.plot(chrom_spect_dict, show=False, save=True, save_dir=het_figs)

notable = muts_per_strain[muts_per_strain > 500]
others = muts_per_strain[muts_per_strain <= 500]

ab_scores = filtered_df.ab
notable_ab_scores = ab_scores.loc[notable.index]
other_ab_scores = ab_scores.loc[others.index]

ab_p_vals = stat_tests.ab_binomial_test(filtered_df)

fig, ax = plt.subplots(3,3)
ax[0][0].hist(ab_scores)
ax[1][0].hist(ab_p_vals.binom_p)
ax[1][1].hist(ab_p_vals.ztest_p)
ax[0][2].hist(filtered_df.dp, range=(0, 30))
ax[1][2].boxplot(filtered_df.dp)
ax[0][1].hist(filtered_df.dp, bins=50, range=(0,100))
ax[2][0].hist(ab_p_vals.binom_p, bins=50)
ax[2][0].hist(ab_p_vals.ztest_p, bins=50, alpha=0.5)
ax[2][1].hist(ab_p_vals.binom_p, bins=25)
ax[2][2].hist(ab_p_vals[ab_p_vals.binom_p <= 0.05].binom_p, bins=20, alpha=0.50, range=(0, 1))
ax[2][2].hist(ab_p_vals[ab_p_vals.binom_p >= 0.05].binom_p, bins=20, alpha=0.50)

shaky = ab_p_vals[(ab_p_vals.binom_p >= 0.05) & (ab_p_vals.binom_p <= 0.4)]
fig, ax = plt.subplots()
ax.scatter(x=(ab_p_vals.x/ab_p_vals.n), y=(ab_p_vals.binom_p))

fig, ax = plt.subplots(2, 2)
ax[0][0].scatter(ab_p_vals.n, ab_p_vals.x/ab_p_vals.n)
ax[0][1].scatter(ab_p_vals.n, ab_p_vals.binom_p)

temp_ab = ab_p_vals.set_index('strain')
df1 = temp_ab.loc[notable.index]
df2 = temp_ab.loc[others.index]

ax[1][0].scatter(df1.n, df1.x/df1.n, alpha=0.5)
ax[1][0].scatter(df2.n, df2.x/df2.n, alpha=0.5)
ax[1][1].scatter(df1.x/df1.n, df1.binom_p, alpha=0.5)
ax[1][1].scatter(df2.x/df2.n, df2.binom_p, alpha=0.5)