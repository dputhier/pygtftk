"""
Contains various utility functions relative to the negative binomial distribution.
"""

import scipy
import scipy.stats
import numpy as np


def check_negbin_adjustment(obs,mean,var):
    """
    Considers a Negative Binomial distribution of given mean and var, and
    performs an adjustement test of this distrubution with regards to obs.

    If mean > 100, the negative binomial distributions is approximated by a
    normal distribution and the Kolmogorov-Smirnov test is used.
    Cannot be used if mean < 100 (which should not happen in our case)

    Returns the result of a KS test (a low p-value means the distributions
    are different).
    """
    norm_equivalent = scipy.stats.norm(mean,np.sqrt(var))
    result = scipy.stats.kstest(obs,norm_equivalent.cdf)
    return result


def log_nb_pval(k,mean,var):
    """
    Log p-value for a negative binomial of those moments.
    """
    r = mean**2 / (var-mean) ; p = 1/(mean/r + 1) # Conversion
    return scipy.stats.nbinom(r,p).logcdf(k)


def empirical_p_val(x,data):
    """
    Quick wrapper : returns the proportion of elements greater than x in the data
    """
    arr = np.array(data)
    where = np.where(arr > x)[0]
    return len(where)/len(arr)
