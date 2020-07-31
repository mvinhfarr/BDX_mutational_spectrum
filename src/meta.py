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
    if pd.isnull(s):
        return None
    elif s[0] != 'F':
        return None
    else:
        n = re.findall(r'\d+', s)
        n = list(map(int, n))
        return sum(n)


def gen_at_seq(df):
    gen = df['Generation at sequencing']
    gen.rename('gen', inplace=True)
    gen = gen.map(extract_gen)
    print(gen.loc[gen.isnull()])
    return gen


def get_epoch(df):
    epoch = df['Epoch']
    epoch.rename('epoch', inplace=True)
    epoch.replace(to_replace='1.5', value='1', inplace=True)

    return epoch


def main(meta_data_csv, snv_df):
    meta_df = load_meta_data(meta_data_csv, snv_df.columns.get_level_values(0))
    generations = gen_at_seq(meta_df)
    epochs = get_epoch(meta_df)

    return meta_df, generations, epochs
