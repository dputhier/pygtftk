"""
Contains various utility functions relative to the negative binomial distribution.
"""

import mpmath
import numpy as np
import scipy
import scipy.stats

from pygtftk.stats.beta import BetaCalculator
from pygtftk.utils import message


def check_negbin_adjustment(obs, mean, var, bins_number=16):
    r"""
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

    Futhermore, the Chi-Squared test can artificially return (H0 false) if n is
    too large. Hence, we use Cramer's V test instead.

    The function returns 1 - V, where V is the Cramer's V test of a crosstab
    made by concatenating the observed counts and expected counts vectors.

    The V test is here used to determine whether there is an association between
    counts distribution and the status as 'expected' or 'observed'.
    As per Cramer (1948) a good fit should have a fit quality above 1 - 0.25 = 0.75
    because V > 0.25 would mean association in most cases.
    Typically, you will observe good association when len(obs) is in the thousands.
    When it is in the hundreds only, random chance may result in a okay-ish
    association (around 0.8). Conversely, if len(obs) is in the thousands, only
    fits above 0.9 should be considered good.

    :param obs: Observed values.
    :param mean: The mean for the negative binomial model.
    :param var: The variance for the negative binomial model.
    :param bins_number: The number of bins in histogram.

    :Example:

    >>> import numpy as np
    >>> import scipy.stats
    >>> from pygtftk.stats.negbin_fit import check_negbin_adjustment
    >>> np.random.seed(42)
    >>> mean = 100 ; var = 200
    >>> r = mean**2 / (var-mean) ; p = 1/(mean/r + 1)
    >>> rv = scipy.stats.nbinom(r, p)
    >>> obs = rv.rvs(2000)
    >>> fit = check_negbin_adjustment(obs, mean, var)
    >>> fit = np.around(fit, 2)
    >>> assert fit > 0.95

    """

    # Force cast obs to integers, just in case.
    obs = np.array(obs, dtype = int)

    if mean == 0:
        mean = 1
        message("Computing adjustment to a Neg Binom with mean = 0 ; mean was set to 1")
        # This is necessary, since r must be above 0.

    if mean >= var:
        var = mean + 1
        message("Computing adjustment to a Neg Binom with mean >= var ; var was set to mean+1")

    # Calculate r and p based on mean and var
    r = mean ** 2 / (var - mean)
    p = 1 / (mean / r + 1)

    rv = scipy.stats.nbinom(r, p)

    # We cannot compare unitary frequencies (frequencies of each individual
    # number), as they are too low. We must use bins.
    obs_range = max(obs) - min(obs)
    step_size = np.around(obs_range / bins_number)
    bin_size = max(1, int(step_size))  # If obs_range < bins_number, step size will be set to 1
    bins = np.arange(min(obs), max(obs), bin_size)

    # There can be a bug later in the count table generation if the range is only 1 or 2
    if obs_range < 3:
        bins = np.arange(min(obs), min(obs) + 3, bin_size)

    # Turn this binned distribution into frequencies
    obs_binned = np.digitize(obs, bins)

    # Remark : the last bin (maximum) will disappear. This is an acceptable loss on this kind of distribution.
    counts = []
    for binned_value in range(len(bins)):
        count = np.sum(
            np.where(obs_binned == binned_value, True, False)
        )
        counts.append(count)
    f_obs = np.array(counts)

    ## Compute the expected frequencies : across each bin, "sum the pmf" (obviously done by the difference of two cdf)
    f_exp = []
    for i in range(len(bins)):

        # Same formula as np.digitize ; special case i=0 will return empty range
        if i == 0:
            bin_left = bin_right = 0
        else:
            bin_left = bins[i - 1]
            bin_right = bins[i]

        total_freq = rv.cdf(bin_right - 1) - rv.cdf(bin_left - 1)

        f_exp.append(total_freq)

    f_exp = np.array(f_exp) * len(obs)

    # Remove leading zero in f_exp and f_obs and cast to a np array of integers
    f_exp = np.array(f_exp[1:], dtype = int)
    f_obs = np.array(f_obs[1:], dtype = int)

    f_exp = f_exp + 1E-100  # The table of expected frequencies must have no zeros.

    # Fuse them in a crosstab for Cramer's V test
    crosstab = np.concatenate([f_obs[:, np.newaxis], f_exp[:, np.newaxis]], axis=1)

    # Compute the test
    def cramers_V(crosstab):
        chi2 = scipy.stats.chi2_contingency(crosstab)[0]
        n = crosstab.sum()
        # print(crosstab.astype(int)) # For debug
        return np.sqrt(chi2 / (n * (min(crosstab.shape) - 1)))

    ## Return (1 - V_score)
    # This way, a good fit has a score around 1 and a bad fit a score around 0
    result = 1 - cramers_V(crosstab)

    return result


def negbin_pval(k, mean, var, precision=320, ft_type="Unknown"):
    r"""
    P-value for a negative binomial distribution of the given moments (mean, var).

    This is the two-sided p-value : it will return the minimum of the left-sided
    and right-sided p-value

    NOTE : To prevent division by zero or negative r, if the mean is higher than
    or equal to the variance, set the variance to mean + epsilon and send a warning

    :param k: the critical value for which the pvalue is computed.
    :param mean: The mean for the negative binomial model.
    :param var: The variance for the negative binomial model.
    :param precision: Floating point precision of mpmath. Should be at least 1000
    :param ft_type: The name of the feature to be tested (just for meaningful messages).

    >>> from pygtftk.stats.negbin_fit import negbin_pval
    >>> mean = 18400
    >>> var = 630200
    >>> k = 65630
    >>> pval = negbin_pval(k, mean, var)
    >>> import math
    >>> assert(math.isclose(pval,1.1999432787236828e-307))

    """

    if mean < 1:
        mean = 1
        msg = "Computing log(p-val) for a Neg Binom with mean < 1 ; mean was set to 1 (" + ft_type + ")"
        message(msg, type='WARNING')
    # This is necessary, since r must be above 0.

    if mean >= var:
        var = mean + 1
        msg = "Computing log(p-val) for a Neg Binom with mean >= var ; var was set to mean+1 (" + ft_type + ")"
        message(msg, type='WARNING')

    # Floating point precision of mpmath. Should be at least 320.
    mpmath.mp.dps = precision

    # Calculate r and p based on mean and var
    r = mpmath.mpf(mean ** 2 / (var - mean))
    p = mpmath.mpf(1 / (mean / r + 1))

    # To circumvent scipy floating point precision issues, we implement a
    # custom p-value calculation (see 'beta.py' for details)
    mybetacalc = BetaCalculator(use_log=True, precision=precision, ft_type=ft_type)
    incomplete_beta = mybetacalc.betainc(a=r, b=k + 1, x=p)
    complete_beta = mybetacalc.beta(a=r, b=k + 1)

    # Take the minimum of CDF and SF
    pval = 1 - (incomplete_beta / complete_beta)
    twosided_pval = min(pval, 1 - pval)

    # Convert back to Python float and return
    return float(twosided_pval)


def empirical_p_val(x, data):
    """
    Quick wrapper : empirical two-sided p value.

    Returns the proportion of elements greater than x (or smaller than x) in
    the data (resp. whichever proportion is smaller).
    This can be used with any dataset, not just a negative-binomial-compliant one.
    """
    arr = np.array(data)

    higher = len(np.where(arr > x)[0])
    lower = len(np.where(arr <= x)[0])
    signif = min(higher, lower)

    return float(signif / len(arr))


'''
import numpy as np

import pandas as pd
#import pybedtools
import scipy
import scipy.stats
from scipy import optimize

#import pymc3 as pm
#import theano.tensor as tt
#import theano

#from multiprocessing import Pool
#import functools as ft
#import time, sys

# This code is completely useless as simply calculating r and p from the
# mean and variance of the observed data is just as efficient.

# I keep it because it may be useful in other problems and provide an example
# of MLE fitting and Bayesian fitting in Python : may be useful in
# multivariate estimations, for example.


# ----------------------------- MLE fitting ---------------------------------- #
# Standard MLE fitting with scipy.optimize

def fit_negbin_mle(obs):

    def negbinom_likelihood_numpy(params,obs):
        # Distribution parameters
        log_mean, log_var = params
        mean, var = np.exp(log_mean), np.exp(log_var)
        r = mean**2 / (var-mean)
        p = 1/(mean/r + 1)

        # Calculating likelihood = sum of PMF
        def log_of_fact(n):
            n = n + 1E-300
            term_one = (n * np.log(n) - n)
            term_two = (1/6) * (8*(n**3) + 4*(n**2) + (1/30))
            term_thee = 0
            return term_one + term_two + (0.5*np.log(np.pi))
        # NOTE : using a worse approximation of factorials resulted in errors
        def log_of_binom(n,p): return log_of_fact(n) - log_of_fact(p) - log_of_fact(n-p)
        def custom_nb_logpmf(k): return log_of_binom(k+r-1,k) + k*np.log(1-p) + r*np.log(p)

        # We use a sum and not a product because we are working with the log pmf
        result = np.sum([custom_nb_logpmf(x) for x in obs])
        return -1 * result # NOTE Don't forget we want to maximize log-likelihood, so minimize -1*log-likelihood !

    E = np.mean(obs)
    V = np.var(obs, ddof=1)

    initial_guess = [np.log(E),np.log(V)]
    # WARNING : when artificially deluding it (ie setting initial guesses at [1,1]), it returns a result where fitted mean == fitted var !

    fitted = optimize.minimize(negbinom_likelihood_numpy,
                        initial_guess, method='Nelder-Mead', # NOTE Must use Nelder-mead
                        args=(obs),
                        options = {'maxiter':1000})
    fitted_log_E, fitted_log_V = fitted.x

    return np.exp(fitted_log_E), np.exp(fitted_log_V)




# ----------------------------- Bayesian fitting ----------------------------- #

# Bayesian fitting may be superfluous; it just inserts the info from tries
# that mean estimator is quickly good. So maybe offer it as an option for more
# robustness.

# Actually a simple optimisation returns the same result as bayesian for a fraction of the time.
# So keep Bayesian in this repository, but implement classical fitting in
# the release version.

def fit_pymc_model(obs):
    """
    Will fit a negative binomial distribution to this list of observations.
    """
    obs = np.array(obs)
    mean_obs = np.mean(obs)
    var_obs = np.var(obs)

    theano.gof.cc.get_module_cache().clear() # Clear theano cache, just in case

    def likelihood(obs, log_of_mean,var_obs,var_mult):
        # Distribution parameters
        mean = tt.exp(log_of_mean)
        var = var_obs*var_mult

        r = mean**2 / (var-mean)
        p = 1/(mean/r + 1)

        # Calculating likelihood = sum of PMF
        def log_of_fact(n):
            return (n * tt.log(n) - n) + (1/6) * (8*(n**3) + 4*(n**2) + (1/30)) + (0.5*tt.log(np.pi))
        def log_of_binom(n,p): return log_of_fact(n) - log_of_fact(p) - log_of_fact(n-p)
        def custom_nb_logpmf(k): return log_of_binom(k+r-1,k) + k*tt.log(1-p)+ r*tt.log(p)
        result = np.sum([custom_nb_logpmf(x) for x in obs])
        return result

    with pm.Model() as model:

        ### Priors
        # We know true mean will likely be close to the observed one
        # Furthermore, pymc has troubles with values>1E10 (memory overflow), so we use a log
        log_of_mean = pm.Normal('Log_of_esperance', mu=np.log(mean_obs), sd=0.25)
        # We know the real variance is likely lower the empirical one (although it could be higher)
        var_mult = pm.Normal('Variance_multiplier', mu=0.75, sd=0.25)

        # Compute the likelihood
        like = pm.DensityDist('Likelihood',
                    likelihood,
                    observed = dict(obs=obs,log_of_mean=log_of_mean,var_obs=var_obs,var_mult=var_mult))
        trace = pm.sample(1000, njobs=4,
                            init='adapt_diag')
        # NOTE See how to reduce verbosity
        # NOTE We use 'adapt_diag' as the initializer because 'jitter'
        # could result in a 'NaN bad energy at init' error

    pm.traceplot(trace)
    pm.summary(trace)

    # Return the fitted esperance and variance
    esperance_fitted = np.exp(pm.summary(trace).loc['Log_of_esperance','mean'])
    variance_fitted = pm.summary(trace).loc['Variance_multiplier','mean'] * var_obs

    return (esperance_fitted,variance_fitted)


'''
