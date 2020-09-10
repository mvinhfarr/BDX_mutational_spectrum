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


def chr_windows_chi2_constant_ht_adjust(df):
    chr_windows = df.copy(deep=True)

    b6_tot_muts = chr_windows.b6_muts.sum()
    d2_tot_muts = chr_windows.d2_muts.sum()
    b6_tot_bp = chr_windows.b6_bp.sum()
    d2_tot_bp = chr_windows.d2_bp.sum()

    chr_windows['chi2_pval'] = np.nan

    for index, window in chr_windows.iterrows():
        if (window.b6_muts | window.d2_muts) == 0:
            chr_windows.loc[index, 'chi2_pval'] = None
            continue

        win_bp_ratio = window.b6_bp / window.d2_bp
        tot_bp_ratio = (b6_tot_bp - window.b6_bp) / (d2_tot_bp - window.d2_bp)

        adjusted_win_b6_muts = window.b6_muts
        adjusted_win_d2_muts = window.d2_muts * win_bp_ratio
        adjusted_tot_b6_muts = b6_tot_muts - window.b6_muts
        adjusted_tot_d2_muts = (d2_tot_muts-window.d2_muts) * tot_bp_ratio

        obs = np.array([[adjusted_win_b6_muts, adjusted_tot_b6_muts],
                        [adjusted_win_d2_muts, adjusted_tot_d2_muts]])
        chi2, p, dof, exp = stats.chi2_contingency(obs, correction=True)

        chr_windows.loc[index, 'chi2_pval'] = p

    return chr_windows


def chr_windows_chi2_elevated_win_ht_adjust(df):
    chr_windows = df.copy(deep=True)

    b6_tot_muts = chr_windows.b6_muts.sum()
    d2_tot_muts = chr_windows.d2_muts.sum()
    b6_tot_bp = chr_windows.b6_bp.sum()
    d2_tot_bp = chr_windows.d2_bp.sum()

    chr_windows['chi2_pval'] = np.nan

    for index, window in chr_windows.iterrows():
        if (window.b6_muts | window.d2_muts) == 0:
            chr_windows.loc[index, 'chi2_pval'] = None
            continue

        if window.b6_bp > window.d2_bp:
            win_bp_ratio = window.d2_bp / window.b6_bp
            tot_bp_ratio = (d2_tot_bp - window.d2_bp) / (b6_tot_bp - window.b6_bp)

            adjusted_win_b6_muts = window.b6_muts * win_bp_ratio
            adjusted_win_d2_muts = window.d2_muts
            adjusted_tot_b6_muts = (b6_tot_muts - window.b6_muts) * tot_bp_ratio
            adjusted_tot_d2_muts = d2_tot_muts - window.d2_muts
        else:
            win_bp_ratio = window.b6_bp / window.d2_bp
            tot_bp_ratio = (b6_tot_bp - window.b6_bp) / (d2_tot_bp - window.d2_bp)

            adjusted_win_b6_muts = window.b6_muts
            adjusted_win_d2_muts = window.d2_muts * win_bp_ratio
            adjusted_tot_b6_muts = b6_tot_muts - window.b6_muts
            adjusted_tot_d2_muts = (d2_tot_muts - window.d2_muts) * tot_bp_ratio

        obs = np.array([[adjusted_win_b6_muts, adjusted_tot_b6_muts],
                        [adjusted_win_d2_muts, adjusted_tot_d2_muts]])
        chi2, p, dof, exp = stats.chi2_contingency(obs, correction=True)

        chr_windows.loc[index, 'chi2_pval'] = p

    return chr_windows


def chr_windows_chi2_pre_win_adjust(df):
    def adjust_window(win):
        if win.b6_bp > win.d2_bp:
            win_bp_ratio = win.d2_bp / win.b6_bp
            adjusted_b6_muts = win.b6_muts * win_bp_ratio
            adjusted_d2_muts = win.d2_muts
        else:
            win_bp_ratio = win.b6_bp / win.d2_bp
            adjusted_b6_muts = win.b6_muts
            adjusted_d2_muts = win.d2_muts * win_bp_ratio
        return adjusted_b6_muts, adjusted_d2_muts

    chr_windows = df.copy(deep=True)

    adj_mut_counts = chr_windows.loc[:, ['chrom', 'window']]
    adj_mut_counts['b6_muts'], adj_mut_counts['d2_muts'] = 0, 0

    for index, window in chr_windows.iterrows():
        adj_b6_muts, adj_d2_muts = adjust_window(window)
        adj_mut_counts.loc[index, 'b6_muts'] = adj_b6_muts
        adj_mut_counts.loc[index, 'd2_muts'] = adj_d2_muts

    b6_tot_muts = adj_mut_counts.b6_muts.sum()
    d2_tot_muts = adj_mut_counts.d2_muts.sum()

    for index, row in adj_mut_counts.iterrows():
        if (row.b6_muts == 0) | (row.d2_muts == 0):
            chr_windows.loc[index, 'chi2_pval'] = None
            continue

        obs = np.array([[row.b6_muts, b6_tot_muts - row.b6_muts],
                        [row.d2_muts, d2_tot_muts - row.d2_muts]])
        chi2, p, dof, exp = stats.chi2_contingency(obs, correction=True)

        chr_windows.loc[index, 'chi2_pval'] = p

    return chr_windows


def chr_windows_chi2_elevated_bp_ht_adjust(df):
    chr_windows = df.copy(deep=True)

    b6_tot_muts = chr_windows.b6_muts.sum()
    d2_tot_muts = chr_windows.d2_muts.sum()
    b6_tot_bp = chr_windows.b6_bp.sum()
    d2_tot_bp = chr_windows.d2_bp.sum()

    chr_windows['chi2_pval'] = np.nan

    for index, window in chr_windows.iterrows():
        if (window.b6_muts | window.d2_muts) == 0:
            chr_windows.loc[index, 'chi2_pval'] = None
            continue

        if window.b6_bp > window.d2_bp:
            win_bp_ratio = window.d2_bp / window.b6_bp
            adjusted_win_b6_muts = window.b6_muts * win_bp_ratio
            adjusted_win_d2_muts = window.d2_muts
        else:
            win_bp_ratio = window.b6_bp / window.d2_bp
            adjusted_win_b6_muts = window.b6_muts
            adjusted_win_d2_muts = window.d2_muts * win_bp_ratio

        if (d2_tot_bp - window.d2_bp) < (b6_tot_bp - window.b6_bp):
            tot_bp_ratio = (d2_tot_bp - window.d2_bp) / (b6_tot_bp - window.b6_bp)
            adjusted_tot_b6_muts = (b6_tot_muts - window.b6_muts) * tot_bp_ratio
            adjusted_tot_d2_muts = d2_tot_muts - window.d2_muts
        else:
            tot_bp_ratio = (b6_tot_bp - window.b6_bp) / (d2_tot_bp - window.d2_bp)
            adjusted_tot_b6_muts = b6_tot_muts - window.b6_muts
            adjusted_tot_d2_muts = (d2_tot_muts - window.d2_muts) * tot_bp_ratio

        obs = np.array([[adjusted_win_b6_muts, adjusted_tot_b6_muts],
                        [adjusted_win_d2_muts, adjusted_tot_d2_muts]])
        chi2, p, dof, exp = stats.chi2_contingency(obs, correction=True)

        chr_windows.loc[index, 'chi2_pval'] = p

    return chr_windows
