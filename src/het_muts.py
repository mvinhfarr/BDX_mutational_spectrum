import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb

import preprocess
import visualize
import meta
import haplotypes
import per_chrom

snv_data_dir = 'data/per_chr_singleton'
snv_data_file = 'data/bxd_singletons_reformatted.csv'
meta_data_file = 'data/strain_summary.csv'
ht_data_dir = 'data/hmm_haplotypes'

results_df = 'out/dfs/'
het_df = 'out/het/dfs/'
het_figs = 'out/het/figs/'

raw_df = pd.read_csv(results_df + 'raw_singletons', index_col=0)

filtered_summary = raw_df[~raw_df['kmer'].str.contains('indel')]
filtered_summary = filtered_summary[filtered_summary['gt'] == 1]
filtered_summary.drop('BXD087/RwwJ', axis='index', inplace=True)

muts_strains = preprocess.mutations_by_strains_df(filtered_summary)

filtered_summary.to_csv(het_df + 'filtered_singletons')
muts_strains.to_csv(het_df + 'mutations_by_strains')

filtered_df = pd.read_csv(het_df + 'filtered_singletons', index_col=0)
mut_strain_df = pd.read_csv(het_df + 'mutations_by_strains', index_col=[0, 1, 2, 3], header=0)

# # load meta data; filter; extract generation at sequence and epoch
# meta_data, gen_seq, epochs = meta.main(meta_data_file, mut_strain_df)
#
# meta_data.to_csv(het_df+'meta_data')
# gen_seq.to_csv(het_df+'gen_at_seq')
# epochs.to_csv(het_df+'epochs')

meta_df = pd.read_csv(het_df + 'meta_data', index_col=0, header=0)
gens_df = pd.read_csv(het_df + 'gen_at_seq', index_col=0, header=0)
epoch_df = pd.read_csv(het_df + 'epochs', index_col=0, header=0)


# snv_df -> count and fraction of total for each the 6 snv types
# bl_df/dba_df -> count and fract of per haplotype total
snv_df, bl_df, dba_df = visualize.mutation_spectrum_barcharts_simple(mut_strain_df, show=False, save=False,
                                                                     results_dir=het_figs)

# snv_frac_per_strain -> divide mut_strain_df by per strain totals
# snv_frac_strain_avg -> collapse into snv and find average strain
# ht_snv_frac_strain_avg -> collapse into snv by ht and find average strain
# snv_tot_frac -> collapse counts across strains (same as snv_df)
# ht_snv_tot_frac -> collapse counts across strains but divide by ht (same as bl_df/dba_df)
snv_frac_per_strain, snv_frac_strain_avg, ht_snv_frac_strain_avg, snv_tot_frac, ht_snv_tot_frac = \
    visualize.mutation_spectrum_barcharts(mut_strain_df, show=False, save=False, results_dir=het_figs)

muts_per_strain, muts_per_strain_per_gen = visualize.strain_distrb(mut_strain_df, epoch_df, gens_df,
                                                                   show=False, save=False, results_dir=het_figs)

visualize.other_bar_charts(mut_strain_df, show=False)

ht_dict, d2_frac_df, muts_per_chrom = haplotypes.main(ht_data_dir, filtered_df)

muts_per_chrom_per_gen = visualize.mutation_rate(muts_per_chrom, epoch_df, gens_df,
                                                 show=False, save=False, results_dir=het_figs)

ht_ratio = visualize.mutation_spectrum_heatmap(mut_strain_df, show=False, save=False, results_dir=het_figs)

chrom_spect_dict = per_chrom.per_chrom_mut_spectrum(filtered_df)
chrom_ratios, chrom_snv_ratios = per_chrom.plot(chrom_spect_dict, show=True, save=False, save_dir=het_figs)

