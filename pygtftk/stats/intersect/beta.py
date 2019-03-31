import mpmath

from pygtftk.utils import message

class BetaCalculator:
    """
    Computing the beta distribution in Python using mpmath.

    Adapted from https://malishoaib.wordpress.com/2014/04/15/the-beautiful-beta-functions-in-raw-python/
    Comments deduced from "Numerical Recipes: The Art of Scientific Computing", Third Edition (2007), Teukolsky et al.

    Formulas :
        self.betainc(A, B, p) <-->
        mpmath.betainc(A, B, 0, p) <-->
        mpmath.quad(lambdat: t**(A-1)*(1-t)**(B-1), [0, p])
    """

    def __init__(self, precision = 1000, use_log = True, itermax = 10000):
        """
        Parameters :
        - precision is the number of significant digits in mpmath.
        - use_log is relevant when calculating beta (using mpmath.gamma or mpmath.loggamma)
        - itermax is the max number of iterations when evaluating the continued fraction form of the incomplete Beta
        """
        self.precision = precision
        self.use_log = use_log

        # Set precision
        mpmath.mp.dps = self.precision

        # Incomplete beta estimation
        self.epsilon = mpmath.mpf("10E-"+str(precision))
        self.itermax = itermax

        # NOTE : default itermax (10000) is suggested in "Numerical Recipes", Third Edition



    # ----------------------------- Complete Beta ---------------------------- #

    def beta(self,a,b):
        """Uses either mpmath's gamma or log-gamma function to compute values of beta function."""

        a,b = mpmath.mpf(a),mpmath.mpf(b)

        if self.use_log :
            beta = mpmath.exp(mpmath.loggamma(a) + mpmath.loggamma(b) - mpmath.loggamma(a+b))
            return beta
        else:
            beta = mpmath.gamma(a)*mpmath.gamma(b)/mpmath.gamma(a+b)
            return beta

    # ----------------------------- Incomplete beta -------------------------- #

    def contfractbeta(self,a,b,x):
        """ Evaluates the continued fraction form of the incomplete Beta function.
        (Code translated from: Numerical Recipes in C.)

        Uses the modified Lentz's method."""

        a,b,x = mpmath.mpf(a),mpmath.mpf(b),mpmath.mpf(x)

        # NOTE : ITBR means 'In The Book Recipe'

        bm = az = am = mpmath.mpf(1.0)

        # These q's will be used in factors that occur in the
        # fracion's coeffiecients
        qab = a+b
        qap = a+1.0
        qam = a-1.0

        bz = 1.0-qab*x/qap # This is 'd' ITBR

        for i in range(self.itermax+1):
            em = mpmath.mpf(i+1) # Proper iteration number from ITBR
            tem = em + em

            # First step of the recurrene (the even one)
            d = em*(b-em)*x/((qam+tem)*(a+tem)) # aa ITBR
            ap = az + d*am # d ITBR - several d=1/d steps simplified
            bp = bz+d*bm # c ITBR

            # Second step of the recurrence (the odd one)
            d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem))
            app = ap+d*az
            bpp = bp+d*bz

            aold = az
            am = ap/bpp
            bm = bp/bpp
            az = app/bpp
            bz = mpmath.mpf(1.0)

            # Are we done ?
            if (abs(az-aold)<(self.epsilon*abs(az))):
                return az

        # If failed to converge, return our best guess but send a warning
        message('a or b too large or given itermax too small for computing incomplete beta function ; pval may be erroneous', type='WARNING')
        return az

    def betaincreg(self, a, b, x):
        ''' betaincreg(a,b,x) evaluates the incomplete beta function (regularized)
        It requires a, b > 0 and 0 <= x <= 1.

        This function requires contfractbeta defined above.
        (Code translated from: Numerical Recipes in C.)'''

        a,b,x = mpmath.mpf(a),mpmath.mpf(b),mpmath.mpf(x)

        # Trivial cases
        if (x<0) | (x>1) : raise ValueError("Bad x in betainc()")
        elif (x == 0): return 0
        elif (x == 1): return 1

        else:

            # Factors in front of the continued fraction
            lbeta = mpmath.loggamma(a+b) - mpmath.loggamma(a) - mpmath.loggamma(b) + a * mpmath.log(x) + b * mpmath.log(1-x)

            # Use continued fraction directly ...
            if (x < (a+1) / (a+b+2)):
                return mpmath.mpf(mpmath.exp(lbeta) * self.contfractbeta(a, b, x) / a)
            # ... or make a symetry transformation first
            else:
                return mpmath.mpf(1 - mpmath.exp(lbeta) * self.contfractbeta(b, a, 1-x) / b)


    def betainc(self, a, b, x):
        """
        Non-regularized incomplete beta.
        """
        return self.betaincreg(a, b, x) * self.beta(a, b)
