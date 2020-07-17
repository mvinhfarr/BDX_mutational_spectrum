import os
import sys
import datetime
import dill
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


def load_data(data):
    chr_singleton_dir = data

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

    df = df.drop(['start', 'end', 'gt', 'dp', 'ab', 'phastCons'], axis=1)

    return df


def get_complement(seq):
    base_map = str.maketrans('ATGC>', 'TACG>')

    return seq.translate(base_map)


def classify_snv(kmer, include_flank_nuc=True):
    snv = [kmer[1], kmer[-2]]

    '''
    if snv[0] == 'G' or snv[0] == 'T':
        kmer = get_complement(kmer)
        snv = [kmer[1], kmer[-2]]
    '''

    if include_flank_nuc:
        five_nuc = kmer[0]
        three_nuc = kmer[2]
        return {'5-base': five_nuc, '3-base': three_nuc, 'snv': snv}
    else:
        return {'snv': snv}


def convert_kmers(df):
    kmer_df = pd.DataFrame(df['kmer'].apply(lambda k: classify_snv(k, include_flank_nuc=True)).values.tolist(), index=df.index)
    df = df.drop('kmer', axis=1)
    converted = pd.concat([df, kmer_df], axis=1)
    return converted


def visualize(df):
    # distribution of snv's
    snv_list = ['A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T']
    plt.bar()


def main():
    data_dir = 'data/per_chr_singleton'

    raw_singleton_summary = load_data(data_dir)
    filtered_singletons = filter_raw_data(raw_singleton_summary)
    collapsed_complements = convert_kmers(filtered_singletons)


    # collapsed_complements.to_csv('out/converted_kmers', sep='\t')


if __name__ == '__main__':
    os.chdir('..')
    main()
