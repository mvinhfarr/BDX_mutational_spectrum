import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
from scipy import stats
from sklearn.decomposition import PCA

import preprocess
import visualize
import meta
import haplotypes
import per_chrom
import stat_tests

snv_data_file = 'data/bxd_singletons_reformatted.csv'
meta_data_file = 'data/strain_summary.csv'
ht_data_dir = 'data/hmm_haplotypes'
ht_dict_dir = 'ht_dict/'

results_df = 'out/dfs/'
results_figs = 'out/figs/'


def per_strain_chrom_windows(filtered_df, ht_dict, chr_windows):
    strains = filtered_df.bxd_strain.unique()
    # chroms = chr_windows.chrom.unique()

    strain_chr_win = pd.DataFrame(index=strains, columns=[chr_windows.chrom, chr_windows.window])
    strain_chr_win.fillna(0, inplace=True)
    # df.drop('chrY', axis=1, level=0, inplace=True)

    for key, df in ht_dict.items():
        if key in strains:
            for index, row in df.iterrows():
                # window with row.start
                start_win = chr_windows[(row.chrom == chr_windows.chrom) & (row.start >= chr_windows.start) & (row.start <= chr_windows.end)]
                # win with row.end
                end_win = chr_windows[(row.chrom == chr_windows.chrom) & (row.end >= chr_windows.start) & (row.end <= chr_windows.end)]

                # only count b6 bp
                if row.ht == 0:
                    if start_win.index == end_win.index:
                        strain_chr_win.loc[key, (row.chrom, start_win.window)] += row.num_bp
                    elif start_win.index == end_win.index - 1:
                        strain_chr_win.loc[key, (row.chrom, start_win.window)] += (start_win.end - row.start).values
                        strain_chr_win.loc[key, (row.chrom, end_win.window)] += (row.end - end_win.start).values
                    else:
                        mid_wins = chr_windows.iloc[start_win.index[0]+1:end_win.index[0], :]
                        strain_chr_win.loc[key, (row.chrom, start_win.window)] += (start_win.end - row.start).values
                        strain_chr_win.loc[key, (row.chrom, end_win.window)] += (row.end - end_win.start).values
                        strain_chr_win.loc[key, (row.chrom, mid_wins.window)] += mid_wins.loc[:, 'step'].values
                elif row.ht == 1:
                    if start_win.index == end_win.index:
                        strain_chr_win.loc[key, (row.chrom, start_win.window)] -= row.num_bp
                    elif start_win.index == end_win.index - 1:
                        strain_chr_win.loc[key, (row.chrom, start_win.window)] -= (start_win.end - row.start).values
                        strain_chr_win.loc[key, (row.chrom, end_win.window)] -= (row.end - end_win.start).values
                    else:
                        mid_wins = chr_windows.iloc[start_win.index[0]+1:end_win.index[0], :]
                        strain_chr_win.loc[key, (row.chrom, start_win.window)] -= (start_win.end - row.start).values
                        strain_chr_win.loc[key, (row.chrom, end_win.window)] -= (row.end - end_win.start).values
                        strain_chr_win.loc[key, (row.chrom, mid_wins.window)] -= mid_wins.loc[:, 'step'].values
    return strain_chr_win


# def strain_chr_win_pca(df):


if __name__ == '__main__':
    filtered_df = pd.read_csv(results_df + 'filtered_df', index_col=0)
    chr_windows = pd.read_csv(results_df + 'chrom_windows', index_col=0)

    meta_df = pd.read_csv(results_df + 'meta_data', index_col=0, header=0)
    gens_df = pd.read_csv(results_df + 'gen_at_seq', index_col=0, header=0)
    epoch_df = pd.read_csv(results_df + 'epochs', index_col=0, header=0)

    ht_dict = {}
    for f in os.listdir(results_df + ht_dict_dir):
        ht_dict[f] = pd.read_csv(results_df + ht_dict_dir + f, index_col=0, header=0)

    # # for later
    # formatted_col_names = ['{}_{}'.format(c, i) for c, i in zip(chr_windows.chrom, chr_windows.window)]
    # formatted_col_names = [col for col in formatted_col_names if 'chrY' not in col]

    strain_chr_win = per_strain_chrom_windows(filtered_df, ht_dict, chr_windows)
