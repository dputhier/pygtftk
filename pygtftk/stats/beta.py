
from pygtftk.utils import message


import mpmath
import numpy as np
import scipy
import scipy.stats

class BetaCalculator:
    r"""
    Computing the beta distribution in Python using mpmath.

    Using routines from the GNU Scientific Library (`beta_inc.c`):
    <Copyright (C) 2007 Brian Gough and Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman>

    Available under the GNU GPL v3 license.

    Formulas :
        self.betainc(A, B, p) <-->
        mpmath.betainc(A, B, 0, p) <-->
        mpmath.quad(lambdat: t**(A-1)*(1-t)**(B-1), [0, p])


    :Example:

    >>> from pygtftk.stats.beta import BetaCalculator
    >>> import mpmath
    >>> mycalc = BetaCalculator()
    >>> res_A = float(mpmath.beta(5,2))
    >>> res_B = float(mycalc.beta(5,2))
    >>> assert(res_A == res_B)
    >>> res_A = float(mpmath.betainc(a=5,b=2,x2=0.5))
    >>> res_B = float(mycalc.betainc(a=5,b=2,x=0.5))
    >>> assert(res_A == res_B)

    """

    def __init__(self, precision=320,
                 use_log=True,
                 itermax=50000,
                 ft_type="Unknown"):
        """
        :param precision: The number of significant digits in mpmath.
        :param use_log: Relevant when calculating beta (using mpmath.gamma or mpmath.loggamma).
        :param itermax: Max number of iterations when evaluating the continued fraction form of the incomplete Beta.

        """
        self.precision = precision
        self.use_log = use_log
        self.ft_type = ft_type

        # Set precision
        mpmath.mp.dps = self.precision

        # Incomplete beta estimation
        self.epsilon = mpmath.mpf("10E-" + str(precision))
        self.itermax = itermax

    # ----------------------------- Complete Beta ---------------------------- #

    def beta(self, a, b):
        """
        Uses either mpmath's gamma or log-gamma function to compute values of the beta function.
        """

        # HOTFIX prevent the original a, b from being numpy elements
        a = float(a);
        b = float(b)

        a, b = mpmath.mpf(a), mpmath.mpf(b)

        if self.use_log:
            beta = mpmath.exp(mpmath.loggamma(a) + mpmath.loggamma(b) - mpmath.loggamma(a + b))
        else:
            beta = mpmath.gamma(a) * mpmath.gamma(b) / mpmath.gamma(a + b)

        return beta

    # ----------------------------- Incomplete beta -------------------------- #

    def contfractbeta(self, a, b, x):
        """
        Evaluates the continued fraction form of the incomplete Beta function.

        Code translated from: GNU Scientific Library

        Uses the modified Lentz's method.

        You can see a representation of this form in the Digial Library of
        Mathematical functions <https://dlmf.nist.gov/8.17#SS5.p1>.
        The goal of the method is to calculate the successive 'd' terms,
        separately for odd and even.
        """

        a, b, x = mpmath.mpf(a), mpmath.mpf(b), mpmath.mpf(x)

        num_term = 1.0
        den_term = 1.0 - (a + b) * x / (a + 1.0)
        den_term = 1.0 / den_term

        cf = den_term

        for i in range(self.itermax + 1):
            k = i + 1
            coeff = k * (b - k) * x / (((a - 1.0) + 2 * k) * (a + 2 * k))

            # First step of the recurrence
            den_term = 1.0 + coeff * den_term
            num_term = 1.0 + coeff / num_term
            den_term = 1.0 / den_term

            delta_frac = den_term * num_term
            cf *= delta_frac

            coeff = -(a + k) * (a + b + k) * x / ((a + 2 * k) * (a + 2 * k + 1.0))

            # Second step
            den_term = 1.0 + coeff * den_term
            num_term = 1.0 + coeff / num_term
            den_term = 1.0 / den_term

            delta_frac = den_term * num_term
            cf *= delta_frac

            # Are we done ?
            if (abs(delta_frac - 1.0) < 2.0 * self.epsilon):
                return cf

        # If failed to converge, return our best guess but send a warning
        msg = 'a or b too large or given itermax too small for computing incomplete'
        msg += ' beta function ; pval may be slightly erroneous for feature (' + self.ft_type + ').'
        message(msg, type='WARNING')
        return cf

    def betaincreg(self, a, b, x):
        """
        betaincreg(a,b,x) evaluates the incomplete beta function (regularized).
        It requires a, b > 0 and 0 <= x <= 1.

        Code translated from: GNU Scientific Library
        """

        # In terms of methods, this function requires contfractbeta(), defined above.

        # HOTFIX prevent the original a, b or x from being numpy elements
        a = float(a);
        b = float(b);
        x = float(x)

        # Transpose a, x, x into mpmath objects.
        a, b, x = mpmath.mpf(a), mpmath.mpf(b), mpmath.mpf(x)

        def isnegint(X):
            return (X < 0) & (X == mpmath.floor(X))

        # Trivial cases
        if (x < 0) | (x > 1):
            raise ValueError("Bad x in betainc(a,b,x) - x must be between 0 and 1")
        elif isnegint(a) | isnegint(b) | isnegint(a + b):
            raise ValueError("Bad a or b in betainc(a,b,x) -- neither a, b nor a+b can be negative integers")
        elif x == 0:
            return 0
        elif x == 1:
            return 1

        else:

            # Factors in front of the continued fraction
            lnbeta = mpmath.loggamma(a) + mpmath.loggamma(b) - mpmath.loggamma(a + b)
            prefactor = -lnbeta + a * mpmath.log(x) + b * mpmath.log(1 - x)

            # Use continued fraction directly ...
            if x < (a + 1) / (a + b + 2):
                fraction = self.contfractbeta(a, b, x)
                result = mpmath.exp(prefactor) * fraction / a

            # ... or make a symmetry transformation first
            else:
                fraction = self.contfractbeta(b, a, 1 - x)
                term = mpmath.exp(prefactor) * fraction / b
                result = 1 - term

            return result

    def betainc(self, a, b, x):
        """
        Non-regularized incomplete beta. See betaincreg().
        """
        return self.betaincreg(a, b, x) * self.beta(a, b)



def fit_beta(obs):
    r"""
    Fits a four-parameter Beta distribution to the given list of observations
    using method-of-moments fitting.

    >>> from scipy.stats import beta
    >>> import numpy.testing as npt
    >>> import numpy as np
    >>> from pygtftk.stats.beta import fit_beta
    >>> a, b = 1., 2.
    >>> np.random.seed(seed=42)
    >>> obs = beta.rvs(a, b, size=10000)  # You need at least 10K for a good estimate (!)
    >>> ahat, bhat, mhat, chat = fit_beta(obs)
    >>> npt.assert_allclose((ahat, bhat), (a,b), rtol = 0.05)

    """

    mean = np.mean(obs)
    var = np.var(obs, ddof = 1)
    skewness = scipy.stats.skew(obs)
    kurtosis = scipy.stats.kurtosis(obs) # Note that we use Disher's definition, so it is really 'excess kurtosis', meaning kurtosis-3

    # First, estimate alpha and beta
    i = kurtosis - skewness**2 +2
    j = 1.5*skewness**2 - kurtosis
    nu = 3*i/j

    if skewness == 0:
        alpha = beta = (1.5*kurtosis +3) / (-kurtosis)

    else:
        term = 1/np.sqrt(
            1 + (16*(nu+1))/((nu+2)**2 * skewness**2)
        )

        if skewness <0:
            alpha = 0.5*nu*(1+term)
            beta = 0.5*nu*(1-term)
        if skewness >0:
            alpha = 0.5*nu*(1-term)
            beta = 0.5*nu*(1+term)

    # To estimate the spread = c-a
    spread = 0.5*np.sqrt(var)*np.sqrt(
        ((2+nu)**2)*(skewness**2)+16*(1+nu)
    )

    a = mean - (alpha/nu)*spread
    c = spread + a

    return alpha, beta, a, c


def beta_pval(k, obs,
            precision=320):
    r"""
    P-value for the given critical value against a beta distribution fitted to
    the given list of observations.

    This is the two-sided p-value : it will return the minimum of the left-sided
    and right-sided p-value

    :param k: the critical value whose p-value will be calculated
    :param obs: the list of observations to which the beta distribution will be fitted
    :param precision: Floating point precision of mpmath. Should be at least 1000  

    >>> from pygtftk.stats.beta import beta_pval
    >>> from scipy.stats import beta
    >>> import numpy.testing as npt
    >>> a, b = 1., 2.
    >>> np.random.seed(seed=42)
    >>> obs = beta.rvs(1, 2, size=10000)
    >>> k = 0.6
    >>> p = 1 - beta.cdf(k,a,b)
    >>> phat = beta_pval(k, obs)                        # Test the combined package
    >>> x = 0.9999
    >>> cp = beta.cdf(x,a,b)
    >>> mybetacalc = BetaCalculator()
    >>> cphat = mybetacalc.betaincreg(a=a, b=b, x=x)    # Test just the p-value
    >>> npt.assert_allclose(p, phat, rtol=0.1)          # Beta approximation is too imprecise for extreme p-values, however.
    >>> npt.assert_allclose(float(cp), float(cphat), rtol=0.05)
    
    """

    # Floating point precision of mpmath. Should be at least 320.
    mpmath.mp.dps = precision

    # Fit the distribution parameters
    alpha, beta, a, c = fit_beta(obs)

    # x is k normalized to be between 0 and 1
    x = mpmath.mpf((k-a)/(c-a))

    # Sanity checks: x must be between 0 and 1 included, and alpha and beta must be strictly positive
    alpha = max(1E-320, alpha)
    beta = max(1E-320, beta)
    x = max(0,x)
    x = min(1,x)

    # Custom p-value calculation (see 'beta.py' for details)
    mybetacalc = BetaCalculator(use_log=True, precision=precision)
    # For Beta, this is simply the *regularized* incomplete beta function
    incomplete_beta = mybetacalc.betaincreg(a=alpha, b=beta, x=x)

    # Take the minimum of CDF and SF
    pval = 1 - (incomplete_beta)
    twosided_pval = min(pval, 1 - pval)

    # Convert back to Python float and return
    return float(twosided_pval)