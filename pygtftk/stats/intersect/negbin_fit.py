"""
Contains various utility functions relative to the negative binomial distribution.
"""

import scipy
import scipy.stats
import numpy as np

from pygtftk.utils import message


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

    This is the two-sided p-value : it will return the minimum of the left-sided
    and right-sided p-value
    """
    # To prevent division by zero or negative r, if the mean is higher than
    # or equal to the variance, set it to variance + epsilon and send a warning
    if mean >= var :
        mean = var + 1E-4
        message("Computing log(p-val) for a Neg Binom with mean >= var ; mean was set to var + 1E-4")

    r = mean**2 / (var-mean) ; p = 1/(mean/r + 1) # Conversion

    left_pval = scipy.stats.nbinom(r,p).logcdf(k)
    right_pval = scipy.stats.nbinom(r,p).logsf(k)

    twosided_pval = min(left_pval,right_pval)

    return


def empirical_p_val(x,data):
    """
    Quick wrapper : empirical two-sided p value.

    Returns the proportion of elements greater than x or smaller than x in the data, whichever proportion is smaller.

    This can be used with any dataset, not just a negative-binomial-compliant one.
    """
    arr = np.array(data)

    higher = len(np.where(arr >= x)[0])
    lower = len(np.where(arr < x)[0])
    signif = min(higher, lower)

    return signif/len(arr)
