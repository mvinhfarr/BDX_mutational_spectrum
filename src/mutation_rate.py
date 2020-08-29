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

df_dir = 'out/dfs/'

raw_df = pd.read_csv(df_dir + 'raw_singletons', index_col=0)

# create two filtered dfs: homozygous only and homo+het muts
homo_filtered = preprocess.filter_raw_data(raw_df, remove_hetero=True)
combined_filtered = preprocess.filter_raw_data(raw_df, remove_hetero=False)
# combined_filtered.drop('BXD087/RwwJ', axis='index', inplace=True)

homo_muts = preprocess.mutations_by_strains_df(homo_filtered)
combined_muts = preprocess.mutations_by_strains_df(combined_filtered)

homo_meta_data, homo_gens, homo_epochs = meta.main(meta_data_file, homo_muts)
combined_meta_data, combined_gens, combined_epochs = meta.main(meta_data_file, combined_muts)

homo_muts_per_strain, homo_muts_per_strain_per_gen = visualize.strain_distrb(homo_muts, homo_epochs, homo_gens)
combined_muts_per_strain, combined_muts_per_strain_per_gen = visualize.strain_distrb(combined_muts, combined_epochs, combined_gens)
