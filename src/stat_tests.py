import math
import numpy as np
import pandas as pd
from scipy import stats


def snv_chi_sq_test(snv, tot):
    obs = np.array([[snv.BL, tot.BL-snv.BL],
                    [snv.DBA, tot.DBA-snv.DBA]])
    chi2, p, dof, exp = stats.chi2_contingency(obs, correction=False)
    return p


# df --> mut_strain_df
def mut_spec_chi_sq(df):
    df = df.sum(axis=1).unstack(level=-1)
    tot = df.sum(axis=0)

    pvals = df.apply(lambda x: snv_chi_sq_test(x, tot), axis=1)
    return pvals.unstack(level=-1)


def ztest(x, n, p=0.5):
    z = (x/n - p) / math.sqrt(p*(1-p) / n)
    p = stats.norm.sf(abs(z))*2

    return p


# het filtered
def ab_binomial_test(df, p_ab=0.5):
    ab_scores = df.ab
    reads = df.dp
    successes = reads * ab_scores

    p_vals = []
    for strain, x, n in zip(df.index, successes, reads):
        binom_p = stats.binom_test(x, n, p=p_ab)
        ztest_p = ztest(x, n, p=p_ab)
        p_vals.append([strain, x, n, binom_p, ztest_p])

    p_vals = pd.DataFrame(p_vals, columns=['strain', 'x', 'n', 'binom_p', 'ztest_p'])

    return p_vals


def chr_windows_chi_sq_test(chr_windows):
    b6_tot_muts = chr_windows.b6_muts.sum()
    d2_tot_muts = chr_windows.d2_muts.sum()
    b6_tot_bp = chr_windows.b6_bp.sum()
    d2_tot_bp = chr_windows.d2_bp.sum()

    chr_windows['chi2_pval'] = np.nan

    for index, row in chr_windows.iterrows():
        try:
            obs = np.array([[row.b6_bp / row.b6_muts, (b6_tot_bp - row.b6_bp) / (b6_tot_muts - row.b6_muts)],
                            [row.d2_bp / row.d2_muts, (d2_tot_bp - row.d2_bp) / (d2_tot_muts - row.d2_muts)]])
            chi2, p, dof, exp = stats.chi2_contingency(obs, correction=True)
        except ValueError:
            print('window {} missing mutations'.format(index))
            continue
        except ZeroDivisionError:
            print('window {} error'.format(index))
            continue
        finally:
            chr_windows.loc[index, 'chi2_pval'] = p

    return chr_windows
