import datetime
import os
import sys
import dill

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt


def load_data(data_dir):
    chr_singletons = []

    for f_name in os.listdir(data_dir):
        path = os.path.join(data_dir, f_name)
        chr_mut_spec = pd.read_csv(path, header=0, index_col=None)
        chr_singletons.append(chr_mut_spec)

    raw_mut_spec = pd.concat(chr_singletons, axis=0)
    return raw_mut_spec


def filter_raw_data(raw_summary, remove_hetero=True):
    df = pd.DataFrame(raw_summary.loc[:, ['chrom', 'bxd_strain', 'haplotype', 'kmer', 'gt']])
    # df = raw_summary.copy(deep=True)

    # remove indel mutations
    df = df[~df['kmer'].str.contains('indel')]

    # remove heterozygous mutations
    if remove_hetero:
        df = df[df['gt'] == 2]

    df = df.drop(['gt'], axis=1)

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
        return {'5-base': five_nuc, '3-base': three_nuc, 'snv': '>'.join(snv)}
    else:
        return {'snv': '>'.join(snv)}


def convert_kmers(df):
    kmer_df = pd.DataFrame(df['kmer'].apply(lambda k: classify_snv(k, include_flank_nuc=True)).values.tolist(), index=df.index)
    df = df.drop('kmer', axis=1)
    converted = pd.concat([df, kmer_df], axis=1)
    return converted


def visualize(df):
    snv_list = ['A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T']

    # distribution of snv's
    snv_counts = df['snv'].value_counts().reindex(index=snv_list)
    plt.bar(snv_list, snv_counts)
    plt.savefig('out/snv_distr.pdf')

    # distribution by haplotype
    bl6j_df = df[df['haplotype'] == 0]
    dba2j_df = df[df['haplotype'] == 1]

    bl6j_snv_counts = bl6j_df['snv'].value_counts().reindex(index=snv_list)
    dba2j_snv_counts = dba2j_df['snv'].value_counts().reindex(index=snv_list)

    x = np.arange(len(snv_list))
    width = 0.35

    fig, ax = plt.subplots()
    bl6j_bar = ax.bar(x - width/2, bl6j_snv_counts, width=width, label='C57BL/6J')
    dba2j_bar = ax.bar(x + width / 2, dba2j_snv_counts, width=width, label='DBA/2J')

    ax.set_ylabel('Counts')
    ax.set_title('Mutation Spectrum Counts by Haplotype')
    ax.set_xticklabels(snv_list)
    ax.legend()

    fig.tight_layout()
    # fig.savefig('out/snv_distr_by_haplo_new.pdf')

    plt.show()
