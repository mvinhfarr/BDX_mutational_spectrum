import os
import pandas as pd
import matplotlib.pyplot as plt

import preprocess
import visualize

os.chdir('..')

data_dir = 'data/per_chr_singleton'

results_df = 'out/dfs/'
results_figs = 'out/figs/'

# preprocess.main(data_dir, results_df)

raw_df = pd.read_csv(results_df+'raw_singleton_summary', index_col=0)
filtered_df = pd.read_csv(results_df+'filtered_singletons', index_col=0)
formatted_df = pd.read_csv(results_df+'formatted_singletons', index_col=0)
muts_strains = pd.read_csv(results_df+'mutation_strains', index_col=[0, 1, 2, 3])

visualize.mutation_spectrum_barcharts(muts_strains, show=True, save=True, results_dir=results_figs)

# preprocess.visualize(formatted_singletons)

# preprocess.visualize_strain_distrb(formatted_singletons)

# ratio_props = preprocess.mutation_spectrum_heatmap(filtered_singletons)
# ratio_props.to_csv(results_dir + 'ratio_props')

