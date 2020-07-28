import numpy as np
import pandas as pd


def add_leading_zero(s):
    bxd = s[0:3]
    num = s[3:]

    if len(num) == 1:
        s = bxd + '00' + num
    elif len(num) == 2:
        s = bxd + '0' + num

    return s


def convert_strain_name(name):
    s = name.split('_')

    if len(s) == 2:
        s[0] = add_leading_zero(s[0])
        s = s[0]
        return s
    elif len(s) == 3:
        s[0] = add_leading_zero(s[0])
        s = s[0] + '/' + s[1]
        return s
    else:
        return name


def number_strains(name):
    s = name.split('_')[0]

    if s[0:3] == 'BXD':
        num = s[3:]
        print(name + ': ' + num)
    else:
        print(False)


def load_meta_data(meta_data_csv, snv_df):
    df = pd.read_csv(meta_data_csv, header=0)
    df.dropna(axis=1, how='all', inplace=True)
    df.set_index('Expanded name', inplace=True)

    names = snv_df.columns
    temp = names.map(number_strains)

    names = names.map(convert_strain_name)
    # names = names.dropna()
    print(names)

    print(df.index)

    good_strains = names.intersection(df.index)
    bad_strains_meta = df.drop(good_strains, axis=0)
    bad_strains = names.drop(good_strains)
    print(good_strains)
    print(bad_strains)
    # print(snv_df['bxd_strain'])
