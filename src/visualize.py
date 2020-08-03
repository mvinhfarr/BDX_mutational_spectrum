import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb

params = {'legend.fontsize': 'x-small',
          'legend.framealpha': 0.5,
          'legend.title_fontsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'small'}
plt.rcParams.update(params)


def mutation_spectrum_barcharts_simple(muts_strains, show=True, save=False, results_dir=None):
    sum_lvl = 0

    snv_tots = muts_strains.sum(axis=0, level=sum_lvl).sum(axis=1)
    tot = snv_tots.sum()
    snv_frac = snv_tots / tot

    bl_snv_tots = muts_strains.xs('BL', level='ht').sum(axis=0, level=sum_lvl).sum(axis=1)
    bl_tot = bl_snv_tots.sum()
    bl_snv_frac = bl_snv_tots / bl_tot
    dba_snv_tots = muts_strains.xs('DBA', level='ht').sum(axis=0, level=sum_lvl).sum(axis=1)
    dba_tot = dba_snv_tots.sum()
    dba_snv_frac = dba_snv_tots / dba_tot

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    ax1.bar(snv_tots.index, snv_tots)

    ax2.bar(bl_snv_tots.index, bl_snv_tots, align='edge', width=0.35, label='C57BL/6J')
    ax2.bar(dba_snv_tots.index, dba_snv_tots, align='edge', width=-0.35, label='DBA/2J')

    ax3.bar(snv_frac.index, snv_frac)

    ax4.bar(bl_snv_frac.index, bl_snv_frac, align='edge', width=0.35, label='C57BL/6J')
    ax4.bar(dba_snv_frac.index, dba_snv_frac, align='edge', width=-0.35, label='DBA/2J')

    ax1.set_ylabel('Counts')
    ax1.set_xlabel('Mutation Type')
    ax1.set_title('Mutation Spectrum Counts')

    ax2.set_ylabel('Counts')
    ax2.set_xlabel('Mutation Type')
    ax2.set_title('Counts by Haplotype')
    ax2.legend()

    ax3.set_ylabel('Fraction')
    ax3.set_xlabel('Mutation Type')
    ax3.set_title('Mutation Spectrum Fractions')

    ax4.set_ylabel('Fraction')
    ax4.set_xlabel('Mutation Type')
    ax4.set_title('Fractions by Haplotype')
    ax4.legend()

    fig.tight_layout(pad=2.0)

    if save:
        fig.savefig(results_dir + 'mutation_spectrum_distr.pdf')
    if show:
        fig.show()

    snv_df = pd.concat([snv_tots, snv_frac], axis=1, keys=['counts', 'fracs'])
    bl_df = pd.concat([bl_snv_tots, bl_snv_frac], axis=1, keys=['counts', 'fracs'])
    dba_df = pd.concat([dba_snv_tots, dba_snv_frac], axis=1, keys=['counts', 'fracs'])

    return snv_df, bl_df, dba_df


def mutation_spectrum_barcharts(mutation_strain_df, show=True, save=False, results_dir=None):
    # the number of mutations per strain
    mutations_per_strain = mutation_strain_df.sum(axis=0)
    # mutation spectrum as a fraction per strain
    # i.e. num of a mutation type per strain divided by tot num mutations per strain
    mutation_fraction_df = mutation_strain_df / mutations_per_strain

    # find the average mutation spectrum and them sum rows to collapse into 6 snv types
    # e.g. C>T = 0.472115 means on average C>T comprised 47.5% of the mutations in a strain
    snv_strain_avg_frac = mutation_fraction_df.sum(axis=0, level=0).mean(axis=1)
    # average mutation spectrum divided by haplotype
    haplo_snv_avg = mutation_fraction_df.sum(axis=0, level=[0, 3]).mean(axis=1)

    # sum of mutation spectrum --> collapse into 6 snvs and sum for all strains
    snv_tots = mutation_strain_df.sum(axis=0, level=0).sum(axis=1)
    # total total number of mutations
    tot = snv_tots.sum()
    # mutation spectrum as fraction of total mutations
    # e.g. C>T = 0.470441 means 47% of the total mutations recorded were C to T
    snv_frac = snv_tots / tot

    per_ht_snv_tots = mutation_strain_df.sum(axis=0, level=[0, 3]).sum(axis=1)
    per_ht_tot = per_ht_snv_tots.sum(level=1)
    per_ht_snv_fracs = per_ht_snv_tots / per_ht_tot

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    snv_strain_avg_frac.plot(kind='bar', ax=ax1, rot=0, color='green', edgecolor='black', alpha=0.5, label='tot')
    haplo_snv_avg.unstack(1).plot(kind='bar', ax=ax1, rot=0, edgecolor='black', stacked=True, width=0.25, position=0)

    haplo_snv_avg.unstack(1).plot(kind='bar', ax=ax2, rot=0, edgecolor='black')

    snv_frac.plot(kind='bar', ax=ax3, rot=0, edgecolor='black')

    per_ht_snv_fracs.unstack(1).plot(kind='bar', ax=ax4, rot=0, edgecolor='black')

    ax1.set_ylabel('Fraction')
    ax1.set_xlabel('Mutation Type')
    ax1.set_title('Average of Mutation Spectrum as Fraction\n of per Strain Total')
    ax1.legend(title='ht')

    ax2.set_ylabel('Fraction')
    ax2.set_xlabel('Mutation Type')
    ax2.set_title('Average by Haplotype')
    ax2.legend(title='ht')

    ax3.set_ylabel('Fraction')
    ax3.set_xlabel('Mutation Type')
    ax3.set_title('Mutation Spectrum as Fraction of Total')

    ax4.set_ylabel('Fraction')
    ax4.set_xlabel('Mutation Type')
    ax4.set_title('Mutation Spectrum as Fraction of\n per Haplotype Total')
    ax4.legend(title='ht')

    fig.tight_layout()

    if save:
        fig.savefig(results_dir + 'mutation_spectrum_summary.pdf')
    if show:
        fig.show()

    return mutation_fraction_df, snv_strain_avg_frac, haplo_snv_avg, snv_frac, per_ht_snv_fracs


def strain_distrb(muts, epochs, gens, show=True, save=False):
    muts_per_strain = muts.sum(axis=0)
    # strain_counts.sort_index(inplace=True)

    epoch_vals = epochs.epoch.unique()

    muts_per_strain_per_gen = muts_per_strain / gens.gen

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    for n in epoch_vals:
        df = muts_per_strain.loc[epochs.epoch == n]
        ax1.bar(df.index, df, label=n)

    ax2.hist(muts_per_strain, bins=20)

    for n in epoch_vals:
        df = muts_per_strain_per_gen.loc[epochs.epoch == n]
        ax3.bar(df.index, df, label=n)

    ax4.hist(muts_per_strain_per_gen.dropna(), bins=20)
    # for n in epoch_vals:
    #     df = muts_per_strain_per_gen.loc[epochs.epoch == n]
    #     ax4.hist(df, bins=10, alpha=0.5, edgecolor='black', label=n)

    ax1.set_xticks([])
    ax1.set_xlabel('Strains')
    ax1.set_ylabel('# of Mutations')
    ax1.set_title('Number of Mutations per Strain')
    ax1.legend(title='epoch')

    ax2.set_xlabel('# SNVs')
    ax2.set_ylabel('# of Strains')
    ax2.set_title('Distribution of # SNVs per Strain')

    ax3.set_xticks([])
    ax3.set_xlabel('Strains')
    ax3.set_ylabel('# muts / # gens')
    ax3.set_title('Number of Mutations per Strain per Generation')
    ax3.legend(title='epoch')

    ax4.set_xlabel('# SNVs / # Gens')
    ax4.set_ylabel('# of Strains')
    ax4.set_title('Distribution of # SNVs per Strain per Generation')
    # ax4.legend(title='epoch')

    fig.tight_layout()

    if save:
        fig.savefig('out/strain_distrb.pdf')
    if show:
        fig.show()

    print(muts_per_strain)
    print(muts_per_strain_per_gen)
    return muts_per_strain, muts_per_strain_per_gen


def mutation_spectrum_heatmap(mut_df):
    mut_df = mut_df.sum(axis=1)
    tot = mut_df.sum()
    mut_frac = mut_df/tot

    bl_frac = mut_frac.xs('BL', level='ht')
    dba_frac = mut_frac.xs('DBA', level='ht')

    ratio_props = bl_frac/dba_frac

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    sb.heatmap(ratio_props.unstack(level=2), ax=ax1, cmap='bwr', cbar=True, square=True, center=1)
    ax1.hlines(range(4, 96, 4), *ax1.get_xlim())

    ax1.set_title('Ratio of proportions of SNVs between BL6: DBA')
    # ax1.set_ylabel('3\'')

    plt.tight_layout()

    # plt.savefig('out/ratio_props_heatmap.pdf')
    plt.close(1)
    plt.show()

    return mut_frac, ratio_props


def other_bar_charts(mutation_strain_df):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    temp = mutation_strain_df.sum(axis=1)
    temp.unstack([1, 2]).plot(kind='bar', ax=ax3, stacked=True)

    plt.show()


# needs df with epoch
def epoch_bar_charts(mutation_strain_df):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    epoch_mutation_counts = mutation_strain_df.sum(axis=1, level=1).sum(axis=0, level=0)
    epoch_tot = epoch_mutation_counts.sum(axis=0)
    epoch_mut_fracs = epoch_mutation_counts / epoch_tot

    strain_counts = mutation_strain_df.sum(axis=0, level=0)
    strain_tots = strain_counts.sum(axis=0)
    strain_fracs = strain_counts / strain_tots

    temp = strain_fracs.unstack().reset_index()
    temp.rename(columns={0: 'fraction'}, inplace=True)
    print(temp)

    sb.boxplot(x='snv', y='fraction', data=temp, hue='epoch', ax=ax1)

    epoch_mut_fracs.T.boxplot(ax=ax2)

    ax1.set_title('Mutation Spectrum across Epochs')
    ax1.set_ylabel('fractions')
    ax1.legend(title='epoch')

    # ax1 = plt.boxplot(epoch_mutation_counts)
    # print(epoch_mutation_counts)
    # print(epoch_mut_fracs)

    # print(strain_counts)
    print(strain_fracs)
    # print(strain_fracs.index)

    plt.tight_layout()
    plt.show()
