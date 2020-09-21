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


def mut_distance(muts, strain_index):
    per_strain_mut_distance = {}
    per_strain_mut_loc = {}

    for strain in strain_index:
        strain_df = muts.loc[strain]
        per_chrom_mut_distance = {}
        per_chrom_mut_loc = {}
        for chrom in strain_df.chrom.unique():
            df = strain_df.loc[strain_df.chrom == chrom]
            mut_locs = df.start
            mut_distances = mut_locs[1:] - mut_locs[:-1]

            per_chrom_mut_distance[chrom] = mut_distances.to_numpy(dtype=int)
            per_chrom_mut_loc[chrom] = mut_locs.to_numpy(dtype=int)

        per_strain_mut_distance[strain] = per_chrom_mut_distance
        per_strain_mut_loc[strain] = per_chrom_mut_loc

    return per_strain_mut_distance, per_strain_mut_loc


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

# ht_ratio = visualize.mutation_spectrum_heatmap(mut_strain_df, show=False, save=False, results_dir=het_figs)

# visualize.other_bar_charts(mut_strain_df, show=False)

# ht_dict, d2_frac_df, muts_per_chrom = haplotypes.main(ht_data_dir, filtered_df)
#
# muts_per_chrom_per_gen = visualize.mutation_rate(muts_per_chrom, epoch_df, gens_df,
#                                                  show=False, save=True, results_dir=het_figs)
# chrom_spect_dict = per_chrom.per_chrom_mut_spectrum(filtered_df)
# chrom_ratios, chrom_snv_ratios = per_chrom.plot(chrom_spect_dict, show=False, save=True, save_dir=het_figs)

notable = muts_per_strain[muts_per_strain > 500]
others = muts_per_strain[muts_per_strain <= 500]

ab_scores = filtered_df.ab
notable_ab_scores = ab_scores.loc[notable.index]
other_ab_scores = ab_scores.loc[others.index]

# ab_p_vals = stat_tests.ab_binomial_test(filtered_df)

# shaky = ab_p_vals[(ab_p_vals.binom_p >= 0.05) & (ab_p_vals.binom_p <= 0.4)]
#
# temp_ab = ab_p_vals.set_index('strain')
# df1 = temp_ab.loc[notable.index]
# df2 = temp_ab.loc[others.index]

# visualize.ab_scores(ab_p_vals, filtered_df)

# ax[1][0].scatter(df1.n, df1.x/df1.n, alpha=0.5)
# ax[1][0].scatter(df2.n, df2.x/df2.n, alpha=0.5)
# ax[1][1].scatter(df2.x/df2.n, df2.binom_p, alpha=0.5)
# ax[1][1].scatter(df1.x/df1.n, df1.binom_p, alpha=0.5)

# phastcons = filtered_df.phastCons.copy(deep=True)
# phastcons = pd.to_numeric(phastcons, errors='coerce')
# phastcons.dropna(inplace=True)
#
# hist, bins = np.histogram(phastcons, bins=20, density=True)
# values = np.cumsum((0.05 * hist))
#
# homo_phastcons = pd.read_csv(main_dir+'filtered_singletons', index_col=0)
# homo_phastcons = homo_phastcons.phastCons.copy(deep=True)
# homo_phastcons = pd.to_numeric(homo_phastcons, errors='coerce')
# homo_phastcons.dropna(inplace=True)
#
# hist2, bins2 = np.histogram(homo_phastcons, bins=20, density=True)
# values2 = np.cumsum(0.05*hist2)
#
# fig, ax = plt.subplots()
#
# ax.plot(bins[1:], values, 'o-', label='heterozygous')
# ax.plot(bins2[1:], values2, 'o-', label='homozygous')
# ax.set_ylim(0, 1)
# ax.legend()

mut_distances, mut_locs = mut_distance(filtered_df, epoch_df.index)


fig, ax = plt.subplots()

mut_loc_df = pd.DataFrame(mut_locs)

for strain, dict in mut_locs.items():
    i = 0
    for chrom, arr in dict.items():
        ax.scatter(arr, np.full_like(arr, i))
        i += 1

fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()

i = 0
for chrom in filtered_df.chrom.unique():
    if chrom == 'chr15':
        j=0
    df = filtered_df.loc[filtered_df.chrom == chrom]
    ax2.hist(df.start, bins=10000, range=(0, 2e8), alpha=0.5)

    arr = df.start
    ax3.scatter(arr, np.full_like(arr, i), label=chrom)
    i += 1
