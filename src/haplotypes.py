import os
import numpy as np
import pandas as pd


def load_ht_data(ht_dir):
    strain_ht_dict = {}
    for f in os.listdir(ht_dir):
        path = os.path.join(ht_dir, f)
        ht_temp = pd.read_csv(path, index_col=None, header=None)
        ht_temp.columns = ['chrom', 'start', 'end', 'ht']
        strain_ht_dict[f.strip('.csv')] = ht_temp

    # df = pd.concat(strain_ht_dict, axis=0)
    return strain_ht_dict


def calc_chrom_ht_frac(ht_dict):
    d2_chrom_frac_dict = {}

    for strain, df in ht_dict.items():
        df['num_bp'] = df.end - df.start

        chroms = df['chrom'].unique()
        d2_chrom_frac = pd.Series(index=chroms, name='d2_frac')

        for c in chroms:
            bl6_sum = df.loc[(df.chrom == c) & (df.ht == 0), 'num_bp'].sum()
            d2_sum = df.loc[(df.chrom == c) & (df.ht == 1), 'num_bp'].sum()
            c_tot = bl6_sum + d2_sum
            d2_chrom_frac.loc[c] = d2_sum / c_tot

        d2_chrom_frac_dict[strain] = d2_chrom_frac

    d2_frac_df = pd.concat(d2_chrom_frac_dict, axis=1)

    return d2_frac_df


def mutations_per_chrom_per_strain(ht_df, mut_df):
    muts_per_chrom = mut_df.unstack()
    print(muts_per_chrom)


def main(ht_dir, filtered_df):
    ht_dict = load_ht_data(ht_dir)
    d2_frac_df = calc_chrom_ht_frac(ht_dict)
    mutations_per_chrom_per_strain(d2_frac_df, filtered_df)

    # print(ht_dict)
    # print(d2_frac_df)
