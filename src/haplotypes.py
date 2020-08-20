import os
import math
import numpy as np
import pandas as pd


def load_ht_data(ht_dir, filtered_df):
    strain_ht_dict = {}
    for f in os.listdir(ht_dir):
        path = os.path.join(ht_dir, f)
        ht_temp = pd.read_csv(path, index_col=None, header=None)
        ht_temp.columns = ['chrom', 'start', 'end', 'ht']
        name = '_'.join(f.strip('.csv').split('_')[:-1])
        # strain = filtered_df[filtered_df.bxd_strain == name].index.unique()
        # if len(strain) == 0: strain = name
        # else: strain = strain[0]
        # strain_ht_dict[strain] = ht_temp
        strain_ht_dict[name] = ht_temp
    # df = pd.concat(strain_ht_dict, axis=0)
    return strain_ht_dict


def calc_chrom_ht_frac(ht_dict):
    d2_chrom_frac_dict = {}

    for strain, df in ht_dict.items():
        df['num_bp'] = df.end - df.start

        chroms = df['chrom'].unique()
        d2_chrom_frac = pd.Series(index=chroms, name='d2_frac')

        for c in chroms:
            # bl6_sum = df.loc[(df.chrom == c) & (df.ht == 0), 'num_bp'].sum()
            d2_sum = df.loc[(df.chrom == c) & (df.ht == 1), 'num_bp'].sum()
            c_tot = df.loc[df.chrom == c, 'num_bp'].sum()
            d2_chrom_frac.loc[c] = d2_sum / c_tot

        d2_chrom_frac_dict[strain] = d2_chrom_frac

    d2_frac_df = pd.concat(d2_chrom_frac_dict, axis=1)

    return d2_frac_df


def mutations_per_chrom_per_strain(ht_df, filtered_df):
    # chrom_ht_index = pd.MultiIndex.from_product([['bl', 'dba'], ht_df.index], names=['chrom', 'ht'])
    muts_per_chrom = {}

    for strain in filtered_df.index.unique():
        df = filtered_df.loc[strain, :]
        bl_counts = df[df.haplotype == 0].chrom.value_counts()
        dba_counts = df[df.haplotype == 1].chrom.value_counts()
        strain_muts_per_chrom = pd.concat([bl_counts, dba_counts], axis=1, keys=['bl', 'dba'])

        muts_per_chrom[strain] = strain_muts_per_chrom

    muts_per_chrom = pd.concat(muts_per_chrom, axis=1)
    muts_per_chrom.fillna(0, inplace=True)

    # DO I ACTUALLY NEED THIS????
    muts_per_chrom = muts_per_chrom.reindex(ht_df.index, axis='index')

    return muts_per_chrom.stack(level=1)


def chrom_mut_distrb(ht_dict, filtered_df, num_win=10):
    summary_strain = ht_dict['B6_0333_no_gap']
    summary_strain.set_index('chrom')

    # list of chromosomes
    chroms = summary_strain.chrom.unique()
    # list of strains with data
    strains = filtered_df.bxd_strain.unique()

    #
    chr_windows = {}

    for chr in chroms:
        win_size = math.ceil(summary_strain.loc[chr].num_bp / num_win)
        chr_start = summary_strain.loc[chr].start
        chr_end = summary_strain.loc[chr].end

        chr_win_df = pd.DataFrame(index=range(num_win), columns=['start', 'end', 'b6_bp', 'd2_bp', 'crossovers'])
        chr_win_df.start = range(chr_start, chr_end, win_size)

        for key, df in ht_dict.items():
            if key in strains:
                df = df[df.chrom == chr]
                for index, row in df.itterrows():
                    window_cross_occurs_in = (row.end - chr_start) // win_size

        '''
        
        for i in range(num_win):
            win_start = chr_start + (i * win_size)
            win_end = win_start + win_size

            # number of B6 base pairs
            bl_bp = 0
            # number of D2 base pairs
            dba_bp = 0

            for key, df in ht_dict.items():
                if key in strains:
                    df = df[df.chrom == chr]
                    for index, row in df.itterrows():
                        '''


def main(ht_dir, filtered_df):
    ht_dict = load_ht_data(ht_dir, filtered_df)
    d2_frac_df = calc_chrom_ht_frac(ht_dict)

    muts_per_chrom = mutations_per_chrom_per_strain(d2_frac_df, filtered_df)

    return ht_dict, d2_frac_df, muts_per_chrom
    # print(ht_dict)
    # print(d2_frac_df)
