import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def mutation_spectrum_barcharts(muts_strains, show=True, save=False, results_dir=None):
    snv_tots = muts_strains.sum(axis=0, level=0).sum(axis=1)
    tot = snv_tots.sum()
    snv_frac = snv_tots / tot

    bl_snv_tots = muts_strains.xs('BL', level='ht').sum(axis=0, level=0).sum(axis=1)
    bl_tot = bl_snv_tots.sum()
    bl_snv_frac = bl_snv_tots / bl_tot
    dba_snv_tots = muts_strains.xs('DBA', level='ht').sum(axis=0, level=0).sum(axis=1)
    dba_tot = dba_snv_tots.sum()
    dba_snv_frac = dba_snv_tots / dba_tot

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))

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
    if show:
        plt.show()

'''
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





def mutation_spectrum_heatmap(df):
    b_df = df[df['haplotype'] == 0]
    d_df = df[df['haplotype'] == 1]

    b_props = convert_into_proportion(b_df)
    d_props = convert_into_proportion(d_df)

    ratio_props = b_props/d_props

    # tempdf = np.random.random((16,16))

    ax = sb.heatmap(ratio_props.unstack(level=-1), cmap='bwr', cbar=True, square=True, vmax=1.16, vmin=0.84)
    ax.hlines(range(4, 96, 4), *ax.get_xlim())
    plt.title('Ratio of proportions of SNVs between BL6: DBA')
    plt.ylabel('3\'')
    plt.tight_layout()

    # plt.savefig('out/ratio_props_heatmap.pdf')
    plt.show()

    return ratio_props



'''

