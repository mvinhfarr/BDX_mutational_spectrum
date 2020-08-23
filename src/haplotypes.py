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


def chrom_ht_windows(ht_dict, filtered_df, win_size=15e6, num_win=None):
    summary_strain = ht_dict['B6_0333_no_gap']
    summary_strain.set_index('chrom', inplace=True)

    # list of chromosomes
    chroms = summary_strain.index
    # list of strains with data
    strains = filtered_df.bxd_strain.unique()

    # # testing
    # chroms = ['chr1']
    # strains = ['BXD001_TyJ_0361']

    chr_windows = {}

    for chrm in chroms:
        chr_start = summary_strain.loc[chrm].start
        chr_end = summary_strain.loc[chrm].end

        chr_win_df = pd.DataFrame(columns=['start', 'end', 'step', 'b6_bp', 'd2_bp', 'crossovers'])

        if num_win is not None:
            win_arr, win_size = np.linspace(chr_start, chr_end+1, num_win+1, endpoint=True, dtype=int, retstep=True)
        else:
            win_arr = np.arange(chr_start, chr_end, win_size, dtype=int)
            win_arr = np.append(win_arr, chr_end+1)

        chr_win_df.start = win_arr[:-1]
        chr_win_df.end = win_arr[1:] - 1
        chr_win_df.step = chr_win_df.end - chr_win_df.start
        chr_win_df.b6_bp = 0
        chr_win_df.d2_bp = 0
        chr_win_df.crossovers = 0

        for key, df in ht_dict.items():
            if key in strains:
                df = df[df.chrom == chrm]
                for index, row in df.iterrows():
                    # window with row.start
                    start_win = chr_win_df[(row.start >= chr_win_df.start) & (row.start <= chr_win_df.end)].index.values[0]
                    # win with row.end
                    end_win = chr_win_df[(row.end >= chr_win_df.start) & (row.end <= chr_win_df.end)].index.values[0]

                    row_ht = 'b6_bp' if row.ht == 0 else 'd2_bp'

                    if start_win == end_win:
                        chr_win_df.loc[start_win, row_ht] += row.num_bp
                    elif start_win == end_win-1:
                        chr_win_df.loc[start_win, row_ht] += chr_win_df.loc[start_win, 'end'] - row.start
                        chr_win_df.loc[end_win, row_ht] += row.end - chr_win_df.loc[end_win, 'start']
                    else:
                        mid_wins = list(range(start_win + 1, end_win))
                        chr_win_df.loc[start_win, row_ht] += chr_win_df.loc[start_win, 'end'] - row.start
                        chr_win_df.loc[end_win, row_ht] += row.end - chr_win_df.loc[end_win, 'start']
                        chr_win_df.loc[mid_wins, row_ht] += chr_win_df.loc[mid_wins, 'step']

                    if row.end != chr_win_df.iloc[-1].end:
                        chr_win_df.loc[end_win, 'crossovers'] += 1

        chr_windows[chrm] = chr_win_df

    df = pd.concat(chr_windows, axis=0)
    df.index.set_names(['chrom', 'window'], inplace=True)
    df.reset_index(inplace=True)
    return df


def chrom_win_muts(df, filtered_muts):
    df['b6_muts'] = 0
    df['d2_muts'] = 0

    for mut in filtered_muts.itertuples():
        window = df[(mut.chrom == df.chrom) & (mut.start >= df.start) & (mut.end <= df.end)].index
        ht = 'b6_muts' if mut.haplotype == 0 else 'd2_muts'
        df.loc[window, ht] += 1

    return df


def main(ht_dir, filtered_df):
    ht_dict = load_ht_data(ht_dir, filtered_df)
    d2_frac_df = calc_chrom_ht_frac(ht_dict)

    muts_per_chrom = mutations_per_chrom_per_strain(d2_frac_df, filtered_df)

    # chr_windows = chrom_ht_windows(ht_dict, filtered_df)

    return ht_dict, d2_frac_df, muts_per_chrom
