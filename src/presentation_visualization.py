import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import preprocess
import meta
import haplotypes
import per_chrom
import stat_tests

params = {'legend.fontsize': 'x-small',
          'legend.framealpha': 0.5,
          'legend.title_fontsize': 'small',
          'xtick.labelsize': 'small',
          'ytick.labelsize': 'small',
          'axes.titlesize': 'large',
          'axes.labelsize': 'medium',
          'figure.autolayout': True}
plt.rcParams.update(params)


def mutation_rate(muts, gens, epochs):
    muts_per_strain = muts.sum(axis=0)
    muts_per_strain_per_gen = muts_per_strain / gens.gen

    muts_per_strain.rename('muts', inplace=True)
    muts_per_strain_per_gen.rename('muts_per_gen', inplace=True)

    df = pd.concat([muts_per_strain, muts_per_strain_per_gen, epochs], axis=1)

    fig1, ax1 = plt.subplots(figsize=(7, 6))

    sb.boxplot(x='epoch', y='muts', data=df, color='white', fliersize=0, ax=ax1)
    sb.stripplot(x='epoch', y='muts', data=df, palette='tab10', size=4, ax=ax1)

    ax1.set_xlabel('Epoch')
    ax1.set_ylabel('# of Mutations')
    ax1.set_title('Number of Mutations per Strain')

    plt.tight_layout()
    # plt.savefig('out/figs/presentation/het_muts_per_strain')

    fig2, ax2 = plt.subplots(figsize=(7, 6))

    sb.boxplot(x='epoch', y='muts_per_gen', data=df, color='white', fliersize=0, ax=ax2)
    sb.stripplot(x='epoch', y='muts_per_gen', data=df, palette='tab10', size=4, ax=ax2)

    ax2.set_xlabel('Epoch')
    ax2.set_ylabel('# muts / # gens')
    ax2.set_title('Number of Mutations per Strain per Generation')

    plt.tight_layout()
    # plt.savefig('out/figs/presentation/het_muts_per_strain_per_gen')

    return df


def mut_spectrum(muts):
    strain_mut_frac = muts.sum(axis=0, level=0) / muts.sum(axis=0)
    strain_mut_frac = strain_mut_frac.stack().reset_index()
    strain_mut_frac.rename(columns={'level_1': 'strain', 0: 'frac'}, inplace=True)

    strain_per_ht_mut_frac = muts.sum(axis=0, level=[0, -1]) / muts.sum(axis=0, level=-1)
    strain_per_ht_mut_frac = strain_per_ht_mut_frac.stack().reset_index()
    strain_per_ht_mut_frac.rename(columns={'level_2': 'strain', 0: 'frac'}, inplace=True)

    df = pd.concat([strain_per_ht_mut_frac, strain_mut_frac], axis=0)
    df.fillna('Combined', inplace=True)

    fig, ax = plt.subplots(figsize=(7, 6))

    # sb.barplot(x='snv', y='frac', data=strain_mut_frac, ci=None, color='white', edgecolor='black', linewidth=1.5, ax=ax)
    # sb.barplot(x='snv', y='frac', hue='ht', data=strain_per_ht_mut_frac, ci='sd', errwidth=1, alpha=0.8, ax=ax)


    # sb.boxplot(x='snv', y='frac', data=strain_mut_frac, color='white', ax=ax)
    sb.boxplot(x='snv', y='frac', hue='ht', data=df, ax=ax)

    ax.set_ylabel('Fraction')
    ax.set_xlabel('Mutation Type')
    ax.set_title('Per haplotype mutation spectrum vs overall')
    ax.legend(title='Haplotype')

    plt.tight_layout()

    # plt.savefig('out/figs/presentation/mut_spectrum')

    return strain_mut_frac, strain_per_ht_mut_frac


def mut_spec_heatmap(muts, vrange=0.15):
    p_vals = stat_tests.mut_spec_chi_sq(muts)

    muts = muts.sum(axis=1).unstack(level='ht')
    # sum of all mutations by haplotype
    ht_tot = muts.sum(axis=0)
    # mutation spectrum as fraction of per haplotype mutations i.e. BL6 and DBA sum to 1
    mut_spec_frac = muts / ht_tot

    ratio_props = mut_spec_frac['BL'] / mut_spec_frac['DBA']
    ratio_props = ratio_props.unstack(level=2)

    fig, ax = plt.subplots(figsize=(3, 6))

    annot = p_vals.stack()
    arr = annot.to_numpy()
    arr = np.where(arr < 0.05, '.', arr)
    annot = pd.Series(arr, index=annot.index).unstack(level=-1)
    sb.heatmap(ratio_props, ax=ax, cmap='bwr', cbar=True, square=True,
               vmin=1-vrange, vmax=1+vrange, xticklabels=True, yticklabels=True)
    sb.heatmap(p_vals, ax=ax, cbar=False, square=True, alpha=0.0, mask=p_vals > 0.05, annot=annot, fmt='',
               annot_kws={'color': 'black', 'size': 50, 'va': 'baseline'})
    ax.hlines(range(4, 96, 4), *ax.get_xlim())

    ax.set_title('BL6 vs. D2')
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    ax.set_ylabel('5\'')
    ax.set_yticklabels(ratio_props.index.get_level_values(1))

    plt.tight_layout()
    # plt.savefig('out/figs/presentation/mut_ratio_heatmap')

    return mut_spec_frac, ratio_props, p_vals


def mut_spec_pca(muts, epochs, by_ht=False):
    mut_spec_frac = muts.sum(axis=0, level=[0, 1, 2]) / muts.sum(axis=0)

    arr = mut_spec_frac.T

    x = arr
    x = StandardScaler().fit_transform(x)

    pca = PCA(n_components=2)
    pricipal_components = pca.fit_transform(x)

    pc_df = pd.DataFrame(pricipal_components, index=muts.columns, columns=['pc1', 'pc2'])
    df = pd.concat([pc_df, epochs], axis=1)

    explained_variance = pca.explained_variance_ratio_

    print('Explained vairation: pc1 = {}, pc2 = {}'.format(explained_variance[0], explained_variance[1]))

    fig, ax = plt.subplots()

    plt.xlabel('Principal Component 1 {} var. expl.'.format(explained_variance[0]))
    plt.ylabel('Principal Component 2 {} var. expl.'.format(explained_variance[1]))

    sb.scatterplot(x='pc1', y='pc2', hue='epoch', data=df, palette='tab10', ax=ax)

    return x, df, explained_variance


def het_mut_spec(homo_frac, het_frac):
    homo_frac['gt'] = 'homo'
    het_frac['gt'] = 'het'

    df = pd.concat([homo_frac, het_frac], axis=0)

    fig, ax = plt.subplots()

    sb.boxplot(x='snv', y='frac', hue='gt', data=df, ax=ax)

    ax.set_ylabel('Fraction')
    ax.set_xlabel('Mutation Type')
    ax.set_title('Heterozygous vs. Homozygous')
    ax.legend(title='Genotype')

    plt.tight_layout()

    # plt.savefig('out/figs/presentation/het_v_homo_mut_spectrum')


def ab_dp_distrb(homo_df, het_df):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))

    ax1.hist(homo_df.ab, bins=50, range=(0, 1), label='homo', alpha=0.75)
    ax1.hist(het_df.ab, bins=50, range=(0, 1), label='het', alpha=0.75)
    ax2.hist(homo_df.dp, bins=50, range=(0, 100), label='homo', alpha=0.75)
    ax2.hist(het_df.dp, bins=50, range=(0, 100), label='het', alpha=0.75)

    ax1.set_title('Allele Balance Distribution')
    ax1.set_xlabel('allele balance')
    ax1.set_ylabel('count')
    ax1.legend()

    ax2.set_title('Depth Distribution')
    ax2.set_xlabel('depth')
    ax2.set_ylabel('count')
    ax2.legend()

    plt.tight_layout()
    # plt.savefig('out/figs/presentation/ab_dp_distrb')


def ab_dp_scatter(p_vals):
    fig, ax = plt.subplots(figsize=(10, 10))

    ax.scatter(p_vals.n, p_vals.x / p_vals.n)

    ax.set_title('Depth vs. Allele Balance')
    ax.set_xlabel('depth')
    ax.set_ylabel('allele balance')

    plt.tight_layout()
    # plt.savefig('out/figs/presentation/ab_dp_distrb')


def ab_binom_pval_distrb(p_vals):
    fig, ax = plt.subplots(figsize=(5, 5))

    ax.hist(p_vals.binom_p, bins=20)

    ax.set_title('Binomial Test p-value Distribution')
    ax.set_xlabel('p-val')
    ax.set_ylabel('count')

    plt.tight_layout()
    # plt.savefig('out/figs/presentation/ab_binom_pval_distrb')


def ab_binom_pval_scatter(p_vals):
    fig, ax = plt.subplots(figsize=(10, 8))

    ax.scatter(x=(p_vals.x / p_vals.n), y=p_vals.binom_p)

    ax.set_xlim(0, 1)
    ax.set_title('Allele Balance vs Binomial Test p-val')
    ax.set_xlabel('allele balance')
    ax.set_ylabel('p-val')

    plt.tight_layout()
    # plt.savefig('out/figs/presentation/ab_binom_pval_scatter')


def phastcons_distrb(homo_df, het_df):
    homo_phastcons = homo_df.phastCons.copy(deep=True)
    homo_phastcons = pd.to_numeric(homo_phastcons, errors='coerce')
    homo_phastcons.dropna(inplace=True)

    hist1, bins1 = np.histogram(homo_phastcons, bins=20, density=True)
    values1 = np.cumsum((0.05 * hist1))

    het_phastcons = het_df.phastCons.copy(deep=True)
    het_phastcons = pd.to_numeric(het_phastcons, errors='coerce')
    het_phastcons.dropna(inplace=True)

    hist2, bins2 = np.histogram(het_phastcons, bins=20, density=True)
    values2 = np.cumsum((0.05 * hist2))

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.plot(bins1[1:], values1, 'o-', label='homozygous')
    ax.plot(bins2[1:], values2, 'o-', label='heterozygous')

    ax.set_ylim(0, 1.05)
    # ax.set_title()
    ax.set_xlabel('PhastCons Score')
    ax.set_ylabel('Cumulative Fraction')
    ax.legend(loc='upper left')

    plt.tight_layout()
    # plt.savefig('out/figs/presentation/phastcons_distrb')


def chr_window_bp_heatmap(chr_win_df):
    df = chr_win_df.copy(deep=True)
    df.set_index(['chrom', 'window'], inplace=True)

    df = df.b6_bp / df.d2_bp
    df = df.unstack(level='chrom')
    # df = df.loc[:, chr_win_df.chrom.unique()]
    df = df.reindex(columns=(chr_win_df.chrom.unique()))

    fig, ax = plt.subplots(figsize=(9, 5))
    sb.heatmap(df, cmap='bwr', cbar=True, square=True, vmin=0, vmax=2,
               ax=ax, xticklabels=True, yticklabels=True)
    ax.set_xlabel('chromosome')
    ax.set_ylabel('window')

    plt.tight_layout()
    # plt.savefig('out/figs/presentation/chr_win_num_bp_heatmap')


def chr_window_mut_per_bp_heatmap(chr_win_df, alpha_level=0.001):
    reindex_df = chr_win_df.copy(deep=True)
    reindex_df.set_index(['chrom', 'window'], inplace=True)

    df = reindex_df.b6_muts_per_bp / reindex_df.d2_muts_per_bp
    df = df.unstack(level='chrom')
    # df = df.loc[:, chr_win_df.chrom.unique()]
    df = df.reindex(columns=(chr_win_df.chrom.unique()))

    p_vals = reindex_df.chi2_pval
    p_vals.fillna(1, inplace=True)
    p_vals = p_vals.unstack(level='chrom')
    p_vals = p_vals.reindex(columns=(chr_win_df.chrom.unique()))

    arr = p_vals.to_numpy()
    arr = np.where(arr < alpha_level, '.', arr)
    annot = pd.DataFrame(arr, index=p_vals.index, columns=p_vals.columns)

    fig, ax = plt.subplots(figsize=(9, 5))
    sb.heatmap(df, cmap='bwr', cbar=True, square=True, vmin=0.5, vmax=1.5,
               ax=ax, xticklabels=True, yticklabels=True)
    sb.heatmap(p_vals, cbar=False, square=True, alpha=0.0, mask=p_vals >= alpha_level, annot=annot, fmt='',
               ax=ax, annot_kws={'color': 'black', 'size': 50, 'va': 'baseline'})
    ax.set_xlabel('chromosome')
    ax.set_ylabel('window')

    plt.tight_layout()
    # plt.savefig('out/figs/presentation/chr_win_mut_per_bp_heatmap')


if __name__ == '__main__':
    snv_data_file = 'data/bxd_singletons_reformatted.csv'
    meta_data_file = 'data/strain_summary.csv'
    ht_data_dir = 'data/hmm_haplotypes'
    ht_dict_dir = 'ht_dict/'

    results_df = 'out/dfs/'

    # Load homozygous data
    raw_df = pd.read_csv(results_df + 'raw_singletons', index_col=0)
    filtered_df = pd.read_csv(results_df + 'filtered_singletons', index_col=0)
    mut_strain_df = pd.read_csv(results_df + 'mutations_by_strains', index_col=[0, 1, 2, 3], header=0)

    meta_df = pd.read_csv(results_df + 'meta_data', index_col=0, header=0)
    gens_df = pd.read_csv(results_df + 'gen_at_seq', index_col=0, header=0)
    epoch_df = pd.read_csv(results_df + 'epochs', index_col=0, header=0)

    # Load heterozygous data
    het_filtered_df = preprocess.filter_raw_data(raw_df, remove_hetero=False, remove_homo=True)
    het_filtered_df = het_filtered_df.drop('BXD087/RwwJ', axis='index', inplace=False)
    het_mut_strain_df = preprocess.mutations_by_strains_df(het_filtered_df)

    het_meta_df, het_gens_df, het_epoch_df = meta.main(meta_data_file, het_mut_strain_df)

    '''
    # Load both homo & heterozygous data
    combined_filtered_df = preprocess.filter_raw_data(raw_df, remove_hetero=False, remove_homo=False)
    combined_filtered_df = combined_filtered_df.drop('BXD087/RwwJ', axis='index', inplace=False)
    combined_mut_strain_df = preprocess.mutations_by_strains_df(combined_filtered_df)
    
    combined_meta_df, combined_gens_df, combined_epoch_df = meta.main(meta_data_file, combined_mut_strain_df)
    '''

    # mut_strain_df.drop('BXD194/RwwJ', axis=1, inplace=True)

    '''HOMO Plots'''
    # mut_rate_df = mutation_rate(mut_strain_df, gens_df, epoch_df)
    # mut_frac, ht_mut_frac = mut_spectrum(mut_strain_df)
    mut_3mer_frac, ht_mut_ratio, heatmap_p_vals = mut_spec_heatmap(mut_strain_df, vrange=0.25)

    # # x_df, pca_dc, expl_var = mut_spec_pca(mut_strain_df, epoch_df)

    ''' AND FOR HETEROZYGOUS
    # muts_per_chrom_per_gen = visualize.mutation_rate(muts_per_chrom, epoch_df, gens_df,
    #                                                  show=False, save=True, results_dir=results_figs)

    # chrom_spect_dict = per_chrom.per_chrom_mut_spectrum(filtered_df)
    # chrom_ratios, chrom_snv_ratios = per_chrom.plot(chrom_spect_dict, show=True, save=True, save_dir=results_figs)
    '''

    # ab_p_vals = stat_tests.ab_binomial_test(filtered_df)
    # # ab_dp_distrb(filtered_df, het_filtered_df)
    # ab_dp_scatter(ab_p_vals)
    # ab_binom_pval_distrb(ab_p_vals)
    # ab_binom_pval_scatter(ab_p_vals)

    # ht_strain_dict = {}
    # for f in os.listdir(results_df + ht_dict_dir):
    #     ht_strain_dict[f] = pd.read_csv(results_df + ht_dict_dir + f, index_col=0, header=0)
    # chr_windows = haplotypes.chrom_ht_windows(ht_strain_dict, filtered_df, num_win=10)
    # filtered_df, chr_windows = haplotypes.chrom_win_muts(filtered_df, chr_windows)

    # autosome_chr_windows = chr_windows[(chr_windows.chrom != 'chrX') & (chr_windows.chrom != 'chrY')]
    # autosome_chr_windows = stat_tests.chr_windows_chi2_elevated_bp_ht_adjust(autosome_chr_windows)

    # chr_window_bp_heatmap(chr_windows)
    # chr_window_mut_per_bp_heatmap(autosome_chr_windows)

    '''HET Plots'''
    # het_mut_rate_df = mutation_rate(het_mut_strain_df, het_gens_df, het_epoch_df)
    # het_mut_frac, het_ht_mut_frac = mut_spectrum(het_mut_strain_df)
    # het_mut_spec(mut_frac, het_mut_frac)
    # het_mut_3mer_frac, het_ht_mut_ratio, het_heatmap_p_vals = mut_spec_heatmap(het_mut_strain_df, vrange=0.25)

    # het_ab_p_vals = stat_tests.ab_binomial_test(het_filtered_df)
    # ab_dp_distrb(filtered_df, het_filtered_df)
    # ab_dp_scatter(het_ab_p_vals)
    # ab_binom_pval_distrb(het_ab_p_vals)
    # ab_binom_pval_scatter(het_ab_p_vals)

    # phastcons_distrb(filtered_df, het_filtered_df)

    # combined_mut_rate_df = mutation_rate(combined_mut_strain_df, combined_gens_df, combined_epoch_df)
    # combined_mut_3mer_frac, combined_ht_mut_ratio, combined_heatmap_p_vals = \
    #     mut_spec_heatmap(combined_mut_strain_df, vrange=0.25)

    plt.show()
