import os
import sys
import datetime
import dill
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def load_data():
    per_chr_singleton_dir = 'data/per_chr_singleton'

    chr_singletons = []

    for fname in os.listdir(per_chr_singleton_dir):
        chrom_mut_spec = pd.read_csv(fname, header=0, index_col=0)
        chr_singletons.append(chrom_mut_spec)

    raw_mut_spec = pd.concat(chr_singletons, axis=0)
    return raw_mut_spec


if __name__ == '__main__':
    raw_singleton_summary = load_data()
