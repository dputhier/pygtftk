import mpmath

from pygtftk.utils import message


class BetaCalculator:
    r"""
    Computing the beta distribution in Python using mpmath.

    Using routines from the GNU Scientific Library (`beta_inc.c`):
    <Copyright (C) 2007 Brian Gough and  Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman>

    Available under the GNU GPL v3 license.

    Formulas :
        self.betainc(A, B, p) <-->
        mpmath.betainc(A, B, 0, p) <-->
        mpmath.quad(lambdat: t**(A-1)*(1-t)**(B-1), [0, p])

    >>> from pygtftk.stats.intersect.beta import BetaCalculator
    >>> import mpmath
    >>> mycalc = BetaCalculator()
    >>> res_A = float(mpmath.beta(5,2))
    >>> res_B = float(mycalc.beta(5,2))
    >>> assert(res_A == res_B)
    >>> res_A = float(mpmath.betainc(a=5,b=2,x2=0.5))
    >>> res_B = float(mycalc.betainc(a=5,b=2,x=0.5))
    >>> assert(res_A == res_B)
    """

    def __init__(self, precision=1500,
                 use_log=True,
                 itermax=100000,
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

        a, b = mpmath.mpf(a), mpmath.mpf(b)

        if self.use_log:
            beta = mpmath.exp(mpmath.loggamma(a) + mpmath.loggamma(b) - mpmath.loggamma(a + b))
            return beta
        else:
            beta = mpmath.gamma(a) * mpmath.gamma(b) / mpmath.gamma(a + b)
            return beta

    # ----------------------------- Incomplete beta -------------------------- #

    def contfractbeta(self, a, b, x):
        """
        Evaluates the continued fraction form of the incomplete Beta function.

        Code translated from: GNU Scientific Library

        Uses the modified Lentz's method.
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
            den_term = 1.0 + coeff * den_term;
            num_term = 1.0 + coeff / num_term;
            den_term = 1.0 / den_term;

            delta_frac = den_term * num_term;
            cf *= delta_frac;

            coeff = -(a + k) * (a + b + k) * x / ((a + 2 * k) * (a + 2 * k + 1.0));

            # Second step
            den_term = 1.0 + coeff * den_term;
            num_term = 1.0 + coeff / num_term;
            den_term = 1.0 / den_term;

            delta_frac = den_term * num_term;
            cf *= delta_frac;

            # Are we done ?
            if (abs(delta_frac - 1.0) < 2.0 * self.epsilon):
                return cf

        # If failed to converge, return our best guess but send a warning
        msg = 'a or b too large or given itermax too small for computing incomplete'
        msg += ' beta function ; pval may be slightly erroneous (' + self.ft_type + ').'
        message(msg, type='WARNING')
        return cf

    def betaincreg(self, a, b, x):
        """
        betaincreg(a,b,x) evaluates the incomplete beta function (regularized).
        It requires a, b > 0 and 0 <= x <= 1.

        Code translated from: GNU Scientific Library
        """

        # In terms of routines, this function requires contfractbeta(), defined above.

        a, b, x = mpmath.mpf(a), mpmath.mpf(b), mpmath.mpf(x)

        def isnegint(X):
            return (X < 0) & (X == mpmath.floor(X))

        # Trivial cases
        if (x < 0) | (x > 1):
            raise ValueError("Bad x in betainc(a,b,x) - x must be between 0 and 1")
        elif isnegint(a) | isnegint(b) | isnegint(a + b):
            raise ValueError("Bad a or b in betainc(a,b,x) -- neither a, b nor a+b can be negative integers")
        elif (x == 0):
            return 0
        elif (x == 1):
            return 1

        else:

            # Factors in front of the continued fraction
            lnbeta = mpmath.loggamma(a) + mpmath.loggamma(b) - mpmath.loggamma(a + b)
            prefactor = -lnbeta + a * mpmath.log(x) + b * mpmath.log(1 - x)

            # Use continued fraction directly ...
            if (x < (a + 1) / (a + b + 2)):
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
