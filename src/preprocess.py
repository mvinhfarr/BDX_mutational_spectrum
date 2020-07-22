import os

import pandas as pd


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
    base_map = str.maketrans('ATGC', 'TACG')
    return seq.translate(base_map)


def classify_kmer(kmer, strip_snv=False):
    snv = [kmer[1], kmer[-2]]

    ''' Finds complement
    if snv[0] == 'G' or snv[0] == 'T':
        kmer = get_complement(kmer)
        snv = [kmer[1], kmer[-2]]
    '''

    if not strip_snv:
        five_nuc = kmer[0]
        three_nuc = kmer[2]
        return {'5-base': five_nuc, '3-base': three_nuc, 'snv': '>'.join(snv)}
    else:
        return {'snv': '>'.join(snv)}


def format_kmers(df, strip_snv=False):
    kmer_df = pd.DataFrame(df['kmer'].apply(lambda k: classify_kmer(k, strip_snv=strip_snv)).values.tolist(),
                           index=df.index)
    df = df.drop('kmer', axis=1)
    converted = pd.concat([df, kmer_df], axis=1)
    return converted


def convert_into_fraction(df):
    snv_counts = df['kmer'].value_counts()
    tot = snv_counts.sum()
    snv_fracs = snv_counts / tot

    '''
    for index, row in snv_props.iteritems():
        kmer = classify_kmer(index)
        props_df.loc[(kmer['snv'], kmer['5-base'], kmer['3-base']), 'Frequency'] = row
    '''

    return snv_fracs


def mutations_by_strains_df(filtered_df):
    strains = filtered_df['bxd_strain'].value_counts()
    mut_strain = []

    for index, row in strains.iteritems():
        per_strain_df = filtered_df[filtered_df['bxd_strain'] == index]

        bl_strain_df = per_strain_df[per_strain_df['haplotype'] == 0]
        dba_strain_df = per_strain_df[per_strain_df['haplotype'] == 1]

        bl_snv_counts = bl_strain_df['kmer'].value_counts()
        dba_snv_counts = dba_strain_df['kmer'].value_counts()

        strain_snv_counts = pd.concat([bl_snv_counts, dba_snv_counts], keys=['BL', 'DBA'],
                                      names=['haplotype', 'kmer'])

        # strain_snv_counts = per_strain_df['kmer'].value_counts()
        strain_snv_counts.rename(index, inplace=True)

        # strain_fracs_df = convert_into_fraction(per_strain_df)
        # strain_fracs_df.rename(index, inplace=True)

        mut_strain.append(strain_snv_counts)

    mut_strain_df = pd.concat(mut_strain, axis=1)
    mut_strain_df.fillna(0, inplace=True)

    kmer_index = mut_strain_df.index.get_level_values('kmer')
    haplo_index = mut_strain_df.index.get_level_values('haplotype')
    multi_index = [kmer_index.str[1] + '>' + kmer_index.str[-2], kmer_index.str[0], kmer_index.str[-1], haplo_index]
    mut_strain_df.index = pd.MultiIndex.from_arrays(multi_index, names=['snv', '5\'', '3\'', 'ht'])
    mut_strain_df.sort_index(axis=0, level=0, inplace=True)

    return mut_strain_df


def main(data_dir, results_dir):
    # load raw data into pandas DataFrame
    raw_singleton_summary = load_data(data_dir)
    # raw_singleton_summary.to_csv(results_dir + 'raw_singleton_summary')

    # remove indel mutations and mutations affecting one chromosome
    filtered_singletons = filter_raw_data(raw_singleton_summary)
    # filtered_singletons.to_csv(results_dir + 'filtered_singletons')

    # reformat kmer into 3 cols: '5-base', '3-base', 'snv'
    # formatted_singletons = format_kmers(filtered_singletons)
    # formatted_singletons.to_csv(results_dir + 'formatted_singletons')

    # df comparing mutation spectrum for each strain
    mutations_strains = mutations_by_strains_df(filtered_singletons)
    mutations_strains.to_csv(results_dir + 'mutation_strains')
