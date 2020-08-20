import numpy as np
import pandas as pd
from scipy import stats


def snv_chi_sq_test(snv, tot):
    obs = np.array([[snv.BL, tot.BL-snv.BL],
                    [snv.DBA, tot.DBA-snv.DBA]])
    chi2, p, dof, exp = stats.chi2_contingency(obs)
    return p


# df --> mut_strain_df
def mut_spec_chi_sq(df):
    df = df.sum(axis=1).unstack(level=-1)
    tot = df.sum(axis=0)

    pvals = df.apply(lambda x: snv_chi_sq_test(x, tot), axis=1)
    return pvals.unstack(level=-1)
