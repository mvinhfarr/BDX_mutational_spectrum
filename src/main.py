import os
import distr_by_haplotype


os.chdir('..')
data_dir = 'data/per_chr_singleton'
results_dir = 'out/'

raw_singleton_summary = distr_by_haplotype.load_data(data_dir)
raw_singleton_summary.to_csv(results_dir + 'raw_singleton_summary')

filtered_singletons = distr_by_haplotype.filter_raw_data(raw_singleton_summary)
filtered_singletons.to_csv(results_dir + 'filtered_singletons')

formatted_singletons = distr_by_haplotype.convert_kmers(filtered_singletons)
formatted_singletons.to_csv(results_dir + 'formatted_singletons')

distr_by_haplotype.visualize(formatted_singletons)
