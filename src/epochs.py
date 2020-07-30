import re
import numpy as np
import pandas as pd


def load_meta_data(f_name, strains):
    df = pd.read_csv(f_name, header=0)
    df.dropna(axis=1, how='all', inplace=True)
    df.set_index('Expanded name', inplace=True)

    # filter strains
    df = df.loc[strains, :]

    return df


def extract_gen(s):
    s = s.split('+')

    if len(s) == 1:
        s = s[0]

        if s[0] != 'F':
            return s

        n = re.findall(r'\d+', s)
        n = list(map(int, n))
        return sum(n)
    else:
        s1 = s[0]
        s2 = s[1]
        n = re.findall(r'\d+', s1)
        n += re.findall(r'\d+', s2)
        n = list(map(int, n))
        return sum(n)


def gen_at_seq(df):
    gen = df['Generation at sequencing']
    gen.rename('generation', inplace=True)

    print(gen.loc[gen.isnull()])
    gen.dropna(inplace=True)

    gen = gen.map(extract_gen)

    return gen


def main(meta_data_csv, snv_df):
    meta_df = load_meta_data(meta_data_csv, snv_df.columns.get_level_values(0))
    generations = gen_at_seq(meta_df)

    return meta_df, generations
