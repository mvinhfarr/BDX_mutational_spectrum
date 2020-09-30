import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
from scipy import stats
from sklearn.preprocessing import StandardScaler
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

    multiindex = pd.MultiIndex.from_product([strains, ['b6', 'd2']], names=['strain', 'ht'])
    multicols = pd.MultiIndex.from_frame(chr_windows.loc[:, ['chrom', 'window']])
    strain_chr_win = pd.DataFrame(0, index=multiindex, columns=multicols)

    for key, df in ht_dict.items():
        if key in strains:
            for index, row in df.iterrows():
                temp = chr_windows[chr_windows.chrom == row.chrom]
                # window with row.start
                start_win = temp[(row.start >= temp.start) & (row.start <= temp.end)]
                # win with row.end
                end_win = temp[(row.end >= temp.start) & (row.end <= temp.end)]

                row_ht = 'b6' if row.ht == 0 else 'd2'
                if start_win.index == end_win.index:
                    strain_chr_win.loc[(key, row_ht), (row.chrom, start_win.window)] += row.num_bp
                elif start_win.index == end_win.index - 1:
                    strain_chr_win.loc[(key, row_ht), (row.chrom, start_win.window)] += (start_win.end - row.start).values
                    strain_chr_win.loc[(key, row_ht), (row.chrom, end_win.window)] += (row.end - end_win.start).values
                else:
                    mid_wins = chr_windows.iloc[start_win.index[0]+1:end_win.index[0], :]
                    strain_chr_win.loc[(key, row_ht), (row.chrom, start_win.window)] += (start_win.end - row.start).values
                    strain_chr_win.loc[(key, row_ht), (row.chrom, end_win.window)] += (row.end - end_win.start).values
                    strain_chr_win.loc[(key, row_ht), (row.chrom, mid_wins.window)] += mid_wins.loc[:, 'step'].values

    return strain_chr_win


# def strain_chr_win_pca(df):


if __name__ == '__main__':
    filtered_df = pd.read_csv(results_df + 'filtered_singletons', index_col=0)

    meta_df = pd.read_csv(results_df + 'meta_data', index_col=0, header=0)
    gens_df = pd.read_csv(results_df + 'gen_at_seq', index_col=0, header=0)
    epoch_df = pd.read_csv(results_df + 'epochs', index_col=0, header=0)

    ht_dict = {}
    for f in os.listdir(results_df + ht_dict_dir):
        ht_dict[f] = pd.read_csv(results_df + ht_dict_dir + f, index_col=0, header=0)

    # chr_windows = pd.read_csv(results_df + 'chrom_windows', index_col=0)
    chr_windows = haplotypes.chrom_ht_windows(ht_dict, filtered_df, win_size=10e6)
    # chr_windows = haplotypes.chrom_ht_windows(ht_dict, filtered_df, num_win=10)
    filtered_df, chr_wndows = haplotypes.chrom_win_muts(filtered_df, chr_windows)
    print('created chr_windows, size {}, {}'.format(chr_windows.shape[0], chr_windows.shape[1]))
    # # for later
    # formatted_col_names = ['{}_{}'.format(c, i) for c, i in zip(chr_windows.chrom, chr_windows.window)]
    # formatted_col_names = [col for col in formatted_col_names if 'chrY' not in col]

    strain_chr_win = per_strain_chrom_windows(filtered_df, ht_dict, chr_windows)
    strain_b6_frac = strain_chr_win.xs('b6', level='ht', axis=0) / strain_chr_win.sum(level='strain', axis=0)
    print('created strain_b6_frac, size {}, {}'.format(strain_b6_frac.shape[0], strain_b6_frac.shape[1]))

    x = strain_b6_frac.values
    # x = strain_b6_frac.drop('chrY', axis=1, level=0)
    x = StandardScaler().fit_transform(x)

    print('StandardScaler applied')

    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(x)

    pc_df = pd.DataFrame(principal_components, index=strain_b6_frac.index, columns=['pc1', 'pc2'])

    explained_variance = pca.explained_variance_ratio_

    print('Explained variation: pc1 = {}, pc2 = {}'.format(explained_variance[0], explained_variance[1]))

    epoch_hue = filtered_df.loc[:, ['bxd_strain', 'epoch']].drop_duplicates()
    epoch_hue.set_index('bxd_strain', inplace=True)
    pc_df = pd.concat([pc_df, epoch_hue], axis=1)

    pc_df.loc['BXD32_TyJ_0415', 'epoch'] = 1

    fig, ax = plt.subplots(figsize=(10, 8))

    sb.scatterplot(x='pc1', y='pc2', hue='epoch', data=pc_df, ax=ax, palette='tab10')

    ax.set_xlabel('Principal Component 1    {:.2f}% var. expl.'.format(explained_variance[0] * 100))
    ax.set_ylabel('Principal Component 2    {:.2f}% var. expl.'.format(explained_variance[1] * 100))
    ax.legend()

    # plt.savefig('out/figs/presentation/strain_similarity_pca')
