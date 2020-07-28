import os
import numpy as np
import pandas as pd


def load_ht_data(ht_dir):
    strain_ht_data = []

    for f in os.listdir(ht_dir):
        path = os.path.join(ht_dir, f)
        ht_temp = pd.read_csv(path, index_col=0, header=None)
        strain_ht_data.append(ht_temp)

    df = pd.concat(strain_ht_data, axis=0)

    print(df)
