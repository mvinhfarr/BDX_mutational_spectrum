import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb

import preprocess
import visualize


def per_chrom_mut_spectrum(df):
    per_chrom_muts = {}

    chroms = df.chrom.unique()

    for chr in chroms:
        chr_df = df[df.chrom == chr]
        chr_mut_strain = preprocess.mutations_by_strains_df(chr_df)

        per_chrom_muts[chr] = chr_mut_strain

    return per_chrom_muts


def plot(per_chrom_muts_dict, show=True, save=False, save_dir=None):
    # fig, ax = plt.subplots(3, 7)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,   10))

    chrom_mut_fracs = {}
    chrom_ht_ratios = {}

    for chr, df in per_chrom_muts_dict.items():
        if chr == 'chr12': continue
        mut_frac, ht_ratio = visualize.mutation_spectrum_heatmap(df, per_chrom=True)

        chrom_mut_fracs[chr] = mut_frac
        chrom_ht_ratios[chr] = ht_ratio

    ht_ratio_df = pd.concat(chrom_ht_ratios, axis=1)
    ht_ratio_df.sort_index(axis=1, level=0, inplace=True)

    sb.heatmap(ht_ratio_df.unstack(level=-1), ax=ax1, cmap='bwr', cbar=True, square=True,
               vmin=0.5, vmax=1.5, xticklabels=True, yticklabels=True)

    ax1.hlines(range(4, 96, 4), *ax1.get_xlim())
    ax1.vlines(range(4, 96, 4), *ax1.get_ylim())
    ax1.set_title('Ratio of per Chrom SNV Proportions between BL6 & D2')

    chrom_mut_fracs_df = pd.concat(chrom_mut_fracs, axis=1)
    chrom_mut_fracs_df = chrom_mut_fracs_df.sum(axis=0, level=0)
    chrom_snv_ratio = chrom_mut_fracs_df.xs('BL', axis=1, level='ht') / chrom_mut_fracs_df.xs('DBA', axis=1, level='ht')
    chrom_snv_ratio.sort_index(axis=1, inplace=True)

    sb.heatmap(chrom_snv_ratio, ax=ax2, cmap='bwr', cbar=True, square=True,
               vmin=0.5, vmax=1.5, xticklabels=True, yticklabels=True)

    # ax2.hlines(range(4, 96, 4), *ax2.get_xlim())
    ax2.set_title('Ratio of per Chrom SNV Proportions between BL6 & D2')

    plt.tight_layout()

    if save:
        plt.savefig(save_dir+'per_chrom_ht_ratio_heatmap.pdf')
    if show:
        plt.show()

    return ht_ratio_df, chrom_snv_ratio
