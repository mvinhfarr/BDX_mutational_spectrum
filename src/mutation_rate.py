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

raw_df = pd.read_csv(df_dir + 'raw_singletons', index_col=None)

# create two filtered dfs: homozygous only and homo+het muts
homo_filtered = preprocess.filter_raw_data(raw_df, remove_hetero=True)
combined_filtered = preprocess.filter_raw_data(raw_df, remove_hetero=False)
combined_filtered.drop('BXD087/RwwJ', axis='index', inplace=True)

homo_muts = preprocess.mutations_by_strains_df(homo_filtered)
combined_muts = preprocess.mutations_by_strains_df(combined_filtered)





