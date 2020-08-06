import matplotlib.pyplot as plt

import preprocess
import visualize


def per_chrom_mut_spectrum(df, plot=True):
    per_chrom_muts = {}

    chroms = df.chrom.unique()

    for chr in chroms:
        chr_df = df[df.chrom == chr]
        chr_mut_strain = preprocess.mutations_by_strains_df(chr_df)

        per_chrom_muts[chr] = chr_mut_strain

    return per_chrom_muts


def plot(per_chrom_muts_dict):
    fig, ax = plt.subplots(1, len(per_chrom_muts_dict))

    n = 0
    for chr, df in per_chrom_muts_dict.items():
        if chr == 'chr12': continue
        ax[n] = visualize.mutation_spectrum_heatmap(df, per_chrom=True)
        n += 1

    plt.tight_layout()
    plt.show()

