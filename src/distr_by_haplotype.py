import datetime
import os
import sys
import dill

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt

import seaborn as sb


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


def classify_kmer(kmer, include_flank_nuc=True):
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
    kmer_df = pd.DataFrame(df['kmer'].apply(lambda k: classify_kmer(k, include_flank_nuc=True)).values.tolist(), index=df.index)
    df = df.drop('kmer', axis=1)
    converted = pd.concat([df, kmer_df], axis=1)
    return converted


def visualize_mutation_spectrum(df):
    snv_list = ['A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T']

    # distribution of snv's
    snv_counts = df['snv'].value_counts().reindex(index=snv_list)

    fig1= plt.figure(1)

    plt.bar(snv_list, snv_counts)

    plt.ylabel('Counts')
    plt.xlabel('Mutation Type')
    plt.title('Mutation Spectrum Counts')

    # distribution by haplotype
    bl6j_df = df[df['haplotype'] == 0]
    dba2j_df = df[df['haplotype'] == 1]

    bl6j_snv_counts = bl6j_df['snv'].value_counts().reindex(index=snv_list)
    dba2j_snv_counts = dba2j_df['snv'].value_counts().reindex(index=snv_list)

    x = np.arange(len(snv_list))
    width = 0.35

    fig2, ax2 = plt.subplots()
    ax2.bar(x, bl6j_snv_counts, width=width, align='edge', label='C57BL/6J')
    ax2.bar(x, dba2j_snv_counts, width=-width, align='edge', label='DBA/2J')

    ax2.set_ylabel('Counts')
    ax2.set_title('Mutation Spectrum Counts by Haplotype')
    ax2.set_xticklabels(snv_list)
    ax2.set_xlabel('Mutation Type')
    ax2.legend()

    fig1.tight_layout()
    fig2.tight_layout()
    fig1.savefig('out/snv_distr.pdf')
    fig2.savefig('out/snv_distr_by_haplo_new.pdf')

    plt.show()


def visualize_strain_distrb(df):
    strain_counts = df['bxd_strain'].value_counts()
    fig = plt.figure(1)
    plt.bar(strain_counts.index, strain_counts)
    plt.xticks([])
    plt.xlabel('Strains')
    plt.ylabel('# Mutations')
    plt.title('Number of Mutations per Strains')

    plt.tight_layout()
    # plt.savefig('out/numSNV_per_strain')

    plt.show()
    print(strain_counts)


def convert_into_proportion(df):
    multi_index_names = [['C>A', 'C>G', 'C>T', 'A>G', 'A>C', 'A>T'], ['A', 'T', 'C', 'G']]
    propsdf = pd.DataFrame(index=pd.MultiIndex.from_product(multi_index_names, names=('SNV', '5\'')), columns=['A', 'T', 'C', 'G'], dtype=float)

    snv_counts = df['kmer'].value_counts()
    tot = snv_counts.sum()
    snv_props = snv_counts/tot

    for index, row in snv_props.iteritems():
        kmer = classify_kmer(index)
        propsdf.loc[(kmer['snv'], kmer['5-base']), kmer['3-base']] = row

    return propsdf


def mutation_spectrum_heatmap(df):
    b_df = df[df['haplotype'] == 0]
    d_df = df[df['haplotype'] == 1]

    b_props = convert_into_proportion(b_df)
    d_props = convert_into_proportion(d_df)

    ratio_props = b_props/d_props

    print(ratio_props)

    # tempdf = np.random.random((16,16))

    ax = sb.heatmap(ratio_props, cmap='bwr', cbar=True, square=True)
    plt.title('Ratio of proportions of SNVs between BL6: DBA')
    plt.tight_layout()

    plt.savefig('out/ratio_props_heatmap')
    plt.show()


