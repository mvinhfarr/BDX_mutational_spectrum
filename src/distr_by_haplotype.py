import os
import sys
import datetime
import dill
import pandas as pd
import numpy as np
import matplotlib


def load_data(dir):
    chr_singleton_dir = dir

    chr_singletons = []

    for f_name in os.listdir(chr_singleton_dir):
        path = os.path.join(chr_singleton_dir, f_name)
        chrom_mut_spec = pd.read_csv(path, header=0, index_col=None)
        chr_singletons.append(chrom_mut_spec)

    raw_mut_spec = pd.concat(chr_singletons, axis=0)
    return raw_mut_spec


def filter_raw_data(raw_summary, remove_hetero=True):
    # df = pd.DataFrame(raw_summary.loc[:, ['chrom', 'bxd_strain', 'k']])
    df = raw_summary.copy(deep=True)

    # remove indel mutations
    # df = df.loc(~df['kmer'].str.contains('indel'))
    df = df[~df['kmer'].str.contains('del')]
    # remove heterozygous mutations
    if remove_hetero:
        df = df[df['gt'] == 2]

    df.drop(['start', 'end', 'gt', 'phastCons'], axis=1)

    return df


if __name__ == '__main__':
    os.chdir('..')
    data_dir = 'data/per_chr_singleton'

    raw_singleton_summary = load_data(data_dir)
    filtered_singletons = filter_raw_data(raw_singleton_summary)

    filtered_singletons.to_csv('out/filtered_singletons', sep='\t')
