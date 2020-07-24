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


def mutation_spectrum_barcharts_v1(muts_strains, show=True, save=False, results_dir=None):
    sum_lvl = [0, 1, 2]

    snv_tots = muts_strains.sum(axis=0, level=sum_lvl).sum(axis=1)
    tot = snv_tots.sum()
    snv_frac = snv_tots / tot

    bl_snv_tots = muts_strains.xs('BL', level='ht').sum(axis=0, level=sum_lvl).sum(axis=1)
    bl_tot = bl_snv_tots.sum()
    bl_snv_frac = bl_snv_tots / bl_tot
    dba_snv_tots = muts_strains.xs('DBA', level='ht').sum(axis=0, level=sum_lvl).sum(axis=1)
    dba_tot = dba_snv_tots.sum()
    dba_snv_frac = dba_snv_tots / dba_tot

    if show:
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
            plt.savefig(results_dir + 'mutation_spectrum_distr.pdf')

        plt.show()

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
        plt.savefig(results_dir + 'mutation_spectrum_summary.pdf')
    if show:
        plt.show()

    return mutation_fraction_df, snv_strain_avg_frac, haplo_snv_avg, snv_frac, per_ht_snv_fracs


def visualize_strain_distrb(df):
    strain_counts = df['bxd_strain'].value_counts()
    fig = plt.figure(1)
    plt.bar(strain_counts.index, strain_counts)
    plt.xticks([])
    plt.xlabel('Strains')
    plt.ylabel('# Mutations')
    plt.title('Number of Mutations per Strains')

    plt.tight_layout()
    plt.savefig('out/numSNV_per_strain.pdf')

    fig2 = plt.figure(2)
    plt.hist(strain_counts, bins=20)
    plt.xlabel('# SNVs')
    plt.ylabel('# of Strains')
    plt.title('Distribution of # SNVs per Strain')
    # plt.savefig('out/distrb_snv_strains.pdf')
    # plt.show()
    print(strain_counts)


def mutation_spectrum_heatmap(ratio_props):
    ax = sb.heatmap(ratio_props.unstack(level=-1), cmap='bwr', cbar=True, square=True, vmax=1.16, vmin=0.84)
    ax.hlines(range(4, 96, 4), *ax.get_xlim())
    plt.title('Ratio of proportions of SNVs between BL6: DBA')
    plt.ylabel('3\'')
    plt.tight_layout()

    # plt.savefig('out/ratio_props_heatmap.pdf')
    plt.show()

    return ratio_props


def other_bar_charts(mutation_strain_df):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    temp = mutation_strain_df.sum(axis=1)
    temp.unstack([1]).plot(kind='bar', ax=ax3, stacked=True)

    # plt.show()
