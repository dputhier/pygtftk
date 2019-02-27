"""
Contains various utility functions relative to the negative binomial distribution.
"""

import numpy as np
import scipy
import scipy.stats

from pygtftk.utils import message


def check_negbin_adjustment(obs, mean, var, bins_number = 15):
    """
    Considers a Negative Binomial distribution of given mean and var, and
    performs an adjustement test of this distrubution with regards to obs.

    Returns the result of a Cramer's V adjustment test.

    A classical Chi-squared test cannot be used, since in a typical use case
    sample size will be only 250 (shuffles) and as such expected frequency for
    each value will be well under 5. As such, we perform the Chi-Squared test
    on an equivalent "histogram" : instead of computing the expected number
    of variates exactly equal to, say, 22, we compute how many variates fall in
    the 20-30 range.
    Target bin number is ajustable, final bin number may vary due to rounding.

    Futhermore, the Chi-Squares test can artificially return (H0 false) if n is
    too large. Hence, we use Cramer's V test instead.

    The function returns 1 - V, where V is the Cramer's V test of a crosstab
    made by concatenating the observed counts and expected counts vectors.

    The V test is here used to determine whether there is an association between
    counts distribution and the status as 'expected' or 'observed'.
    As per Cramer (1948) a good fit should have a fit quality above 1 - 0.25 = 0.75
    because V > 0.25 would mean association in most cases.
    """

    if mean >= var:
        var = mean + 1E-4
        message("Computing log(p-val) for a Neg Binom with mean >= var ; var was set to mean + 1E-4")

    # Calculate r and p based on mean and var
    r = mean**2 / (var-mean)
    p = 1/(mean/r + 1)

    rv = scipy.stats.nbinom(r, p)

    # Unitary frequencies are too low, we must use bins.
    obs_range = max(obs) - min(obs)
    step_size = np.around(obs_range/bins_number)
    bin_size = max(1,int(step_size)) # If obs_range < bins_number, step size will be set to 1
    bins = range(min(obs),max(obs),bin_size)

    # There can be a bug later in the count table generation if the range is only 1
    if obs_range < 2 : bins = range(min(obs),min(obs)+2,bin_size)
    
    # Turn this binned distribution into frequencies
    obs_binned = np.digitize(obs, bins)

    counts = []
    for binned_value in range(len(bins)):
        count = sum(obs_binned == binned_value)
        counts.append(count)
    f_obs = np.array(counts)

    # Compute the expected frequencies : for each bin, sum the pmf
    f_exp = []
    for i in range(len(bins)) :
        total_freq = sum([rv.pmf(x) for x in range(bins[i-1],bins[i])]) # Same formula as np.digitize ; special case i=0 will return empty range
        f_exp.append(total_freq)
    f_exp = np.array(f_exp) * len(obs)

    # Remove leading zero in f_exp and f_obs and cast to a np array of integers
    f_exp = np.array(f_exp[1:]).astype(int) ; f_obs = np.array(f_obs[1:]).astype(int)

    f_exp = f_exp + 1E-100 # The table of expected frequencies must have no zeros.

    # Fuse them in a crosstab for Cramer's V test
    crosstab = np.concatenate([f_obs[:,np.newaxis], f_exp[:,np.newaxis]], axis=1)

    # Compute the test
    def cramers_V(crosstab):
        chi2 = scipy.stats.chi2_contingency(crosstab)[0]
        n = crosstab.sum()
        return np.sqrt(chi2 / (n*(min(crosstab.shape)-1)))

    ## Return (1 - V_score)
    # This way, a good fit has a score around 1 and a bad fit a score around 0
    result = 1 - cramers_V(crosstab)

    return result



def log_nb_pval(k, mean, var):
    """
    Log p-value for a negative binomial of those moments.

    This is the two-sided p-value : it will return the minimum of the left-sided
    and right-sided p-value

    NOTE : To prevent division by zero or negative r, if the mean is higher than
    or equal to the variance, set the variance to mean + epsilon and send a warning
    """

    if mean >= var:
        var = mean + 1E-4
        message("Computing log(p-val) for a Neg Binom with mean >= var ; var was set to mean + 1E-4")

    # Calculate r and p based on mean and var
    r = mean**2 / (var-mean)
    p = 1/(mean/r + 1)

    rv = scipy.stats.nbinom(r, p)

    left_pval = rv.logcdf(k)
    right_pval = rv.logsf(k)

    twosided_pval = min(left_pval, right_pval)

    return twosided_pval



def empirical_p_val(x, data):
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
