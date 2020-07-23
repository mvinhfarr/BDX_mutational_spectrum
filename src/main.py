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
mutation_strain_df = pd.read_csv(results_df+'mutation_strains', index_col=[0, 1, 2, 3])

# snv_tots = mutations_strains.sum(axis=0, level=[0, 1, 2]).sum(axis=1)

# snv_tots.unstack([0]).plot(kind='bar')
# plt.bar(snv_tots.index)
# plt.show()

# snv_df, bl_df, dba_df = visualize.mutation_spectrum_barcharts(mutations_strains, show=True, save=False, results_dir=results_figs)

# bl_dba_ratios = bl_df['fracs'] / dba_df['fracs']

# the number of mutations per strain
mutations_per_strain = mutation_strain_df.sum(axis=0)
# mutation spectrum as a fraction per strain
# i.e. num of a mutation type per strain divided by tot num mutations per strain
mutation_fraction_df = mutation_strain_df / mutations_per_strain
# find the average mutation spectrum and them sum rows to collapse into 6 snv types
# e.g. C>T = 0.472115 means on average C>T comprised 47.5% of the mutations in a strain
snv_strain_avg_frac = mutation_fraction_df.mean(axis=1).sum(axis=0, level=0)

# sum of mutation spectrum --> collapse into 6 snvs and sum for all strains
snv_tots = mutation_strain_df.sum(axis=0, level=0).sum(axis=1)
# total total number of mutations
tot = snv_tots.sum()
# mutation spectrum as fraction of total mutations
# e.g. C>T = 0.470441 means 47% of the total mutations recorded were C to T
snv_frac = snv_tots / tot

haplo_snv_avg = mutation_fraction_df.mean(axis=1).sum(axis=0, level=[0, 3])

bl_snv_tots = mutation_strain_df.xs('BL', level='ht').sum(axis=0, level=0).sum(axis=1)
bl_tot = bl_snv_tots.sum()
bl_snv_frac = bl_snv_tots / bl_tot
dba_snv_tots = mutation_strain_df.xs('DBA', level='ht').sum(axis=0, level=0).sum(axis=1)
dba_tot = dba_snv_tots.sum()
dba_snv_frac = dba_snv_tots / dba_tot

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))

# snv_frac.plot(kind='bar', ax=ax1)
snv_strain_avg_frac.plot(kind='bar', ax=ax1)

haplo_snv_avg.unstack(1).plot(kind='bar', ax=ax4)
ax2.bar(bl_snv_frac.index, bl_snv_frac, align='edge', width=0.35, label='C57BL/6J')
ax2.bar(dba_snv_frac.index, dba_snv_frac, align='edge', width=-0.35, label='DBA/2J')

temp = mutation_strain_df.sum(axis=1)
temp.unstack([1]).plot(kind='bar', ax=ax3, stacked=True)

print(mutation_fraction_df)
print(haplo_snv_avg)

fig.tight_layout()
plt.show()

'''
snv_df = mutations_strains.sum(axis=1)
snv_df = snv_df / snv_df.sum()

bl_df = mutations_strains.sum(axis=1)

print(snv_df)
'''
# print(bl_df)
# print(dba_df)
# print(bl_dba_ratios)
