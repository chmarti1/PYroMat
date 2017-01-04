"""PYroMat solver module

This module supplies classes for defining custom inversion routines for
class properties.  Each class defines a function-like object intended to
provide the inverse of thermodynamic properties through an iterative 
solver.  For example, if a species has a property method for enthalpy,
h(T), the solver would supply T(h).

    SOLVE1N     a Newton-iteration-based solver
    SOLVE1B     a bisection-method solver
"""

import pyromat as pyro
import numpy as np



class _proto_solver_(object):
    """The prototype solver class

SOLVE1 and SOLVE2 are built from this prototype class.  _proto_solver_
defines attributes:

    epsilon     fractional precision of the inversion (default = 1e-6)

    small       a tiny number representative of "numerically zero"
                (default = 1e-10)

    max_iter    integer maximum number of iterations allowed 
                (default = 100)

    defaults    A dictionary containing keyword arguments to be passed
                to the bound property method in addition to param

>>> _proto_solver_.__init__(f=O2.h, param='T',  paramval=1000.)
sets param to 'T', paramval to 1000., and f to oxygen's enthalpy
"""
    def __init__(self, 
            epsilon=1e-6, small=1e-10, max_iter=100, limits=None,
            defaults={}):
        # Set up the defaults
        self.defaults = defaults
        self.epsilon = epsilon
        self.small = small
        self.max_iter = max_iter







class solve1n(_proto_solver_):
    """SOLVE1N  a function-like object for inverting property functions
    solver1n = solve1n( param, )

SOLVE1N objects are function-like objects that can be configured to 
iteratively invert PYroMat property methods.  

"""


    def __init__(self,param,f=None,df=None,fdf=None,param_init=None,
            param_lim=None,prop_name=None, **kwarg):
        #
        # Call the parent init
        #
        super(solve1n,self).__init__(**kwarg)
        #
        # Address the 1N specific parameters one-by-one
        #
        # First initialize the 1N parameters
        self._f = None
        self._df = None
        self._fdf = None
        self._prop_name = prop_name
        self.param = None
        self.param_init=None
        self.param_lim=None

        # (1) the param keyword
        # This needs to be a string that indicates which keyword parameter
        # will be iterated on
        if not isinstance(param,str):
            raise pyro.utility.PMAnalysisError(
                "SOLVE1N objects' ""param"" keyword requires a string " + 
                "indicating\nthe property method to invert.")
        # All good; assign it
        self.param = param

        # (2) Now address the function definition
        # There are three possible ways to specify the function:
        #  -> In the simplest approach, a property function is supplied. The 
        #     solver will define its own function for estimating its derivative.
        #  -> The solver can be configured with an external function for 
        #     calculating the property's derivative (df).
        #  -> The solver can be configured with a single function that 
        #     calculates the property and its derivative together (fdf).
        # Regardless, the solver will bundle these into a single method, _fdf
        # that is responsible for the function evaluations.
        if fdf is not None:
            # Verify that it is callable
            if not hasattr(fdf,'__call__'):
                raise pyro.utility.PMAnalysisError(
                    "A SOLVE1N object was initialized with a fdf " + 
                    "(function and derivative)\n that is not callable.")
            # go ahead and assign fdf
            self._fdf = fdf
        elif f is not None:
            if not hasattr(f,'__call__'):
                raise pyro.utility.PMAnalysisError(
                    "A SOLVE1N object was initialized with an f function " + 
                    "that is not callable.""")
            # assign f
            self._f = f
            if df is None:                    
                self._fdf = self._fdf_f
            else:
                if not hasattr(df,'__call__'):
                    raise pyro.utility.PMAnalysisError(
                    "A SOLVE1N object was initialized with a df function " + 
                    "that is not callable.""")
                self._df = df
                self._fdf = self._fdf_f_df
        else:
            raise pyro.utility.PMAnalysisError(
                "A SOLVE1N object was initialized with no property function.")

        # (3) param_init and param_lim cannot simultaneously be None
        if param_init is None and param_lim is None:
            raise pyro.utility.PMAnalysisError(
                "A SOLVE1N object requires either a param_init or a param_lim "+
                "value.")

        # (4) param_init can be a function or a value
        if param_init is not None:
            if hasattr(param_init,'__call__'):
                self.param_init = param_init
            else:
                try:
                    self.param_init = np.array(param_init,dtype=float).squeeze()
                except:
                    raise pyro.utility.PMAnalysisError(
                    "A SOLVE1N object received an illegal value for param_init.")
                if self.param_init.size != 1:
                    raise pyro.utility.PMAnalysisError(
                    "A SOLVE1N object expected a scalar value for param_init, "+
                    "but received\na vector.")

        # (5) param_lim can be a two-element iterable or a function 
        if param_lim is not None:
            if hasattr(param_lim,'__call__'):
                self.param_lim = param_lim
            else:
                try:
                    self.param_lim = np.array(param_lim,dtype=float).squeeze()
                except:
                    raise pyro.utility.PMAnalysisError(
                    "A SOLVE1N object received an illegal value for param_lim.")
                if self.param_lim.size != 2:
                    raise pyro.utility.PMAnalysisError(
                    "A SOLVE1N object expected a two-element iterable for " + 
                    "param_lim,\n but got something else.")
                elif self.param_lim[0]>=self.param_lim[1]:
                    raise pyro.utility.PMAnalysisError(
                    "The upper and lower limits supplied to a SOLVE1N object " +
                    "appear to be\nreversed.")



    def __repr__(self):
        return '<solver {:s}({:s}) for {:s}>'.format(
            self.param,
            self.prop_name(),
            self._f.__self__.data['id'])


    def __call__(self,y=None,**kwarg):

        # If the property value isn't specified it must be in kwarg; get it
        propname = self.prop_name()
        if y is None:
            if propname in kwarg:
                y = thisarg.pop(propname)
            else:
                raise pyro.utility.PMAnalysisError(
                "The property value '{:s}' was not specified.".format(propname))

        # if the argument is iterable, then treat it as an array
        if hasattr(y,'__iter__'):
            return np.array([self(thisy,**kwarg) for thisy in y])


        # We will let "thisarg" be a dictionary of arguments to pass to the
        # property funciton.  It will include the iteration parameter, 
        # self.param, the solver defaults, self.defaults, and any extras the
        # user passed in the function call.
        # Let x be the value of the parameter under iteration.
        # Let y be the target value for self._f

        # Start with the defaults and override with any new arguments
        thisarg = self.defaults.copy()
        thisarg.update(kwarg)

        # Establish the limits
        if self.param_lim is not None:
            if hasattr(self.param_lim,'__call__'):
                limits = self.param_lim(**thisarg)
            else:
                limits = self.param_lim.copy()
        else:
            limits = None

        # Establish the initial conditions
        if self.param_init is not None:
            if hasattr(self.param_init,'__call__'):
                x = self.param_init(**thisarg)
            else:
                x = self.param_init
        else:
            x = np.average(limits)

        # Establish the error threshold for convergence
        threshold = np.max(self.small,np.abs(self.epsilon*y))

        # Define fake old x and error values
        aerror_old = float('inf')
        x_old = x
        dx = 0.

        # Start the iteration
        for count in range(self.max_iter):
            # Make an initial call to the function
            thisarg[self.param] = x
            f,df = self._fdf(**thisarg)
            error = f - y
            aerror = np.abs(error)
            # Test for convergence
            if aerror < threshold:
                return x

            # Catch a bad guess - if the error has increased,
            # don't keep going; back-track
            if aerror > aerror_old:
                dx /= 2.
                x = x_old + dx
            else:
                # If the guess is good, project the next x value and shift "old"
                # values.
                dx = -error / df
                x_old = x
                aerror_old = aerror
                x += dx
                # enforce the limits
                if limits is not None:
                    x = max(min(x,limits[1]), limits[0])
        raise pyro.utility.PMAnalysisError(
            "The solution failed to converge " + repr(self))


    def _fdf_f(self, **kwarg):
        """Provide f() and df() from f only"""
        # Evaluate the property
        f0 = self._f(**kwarg)
        # Get the parameter and perturb it
        x0 = kwarg[self.param]
        dx = max(abs(x0*self.epsilon), self.small)
        x1 = x0+dx
        kwarg[self.param] = x1
        # Evaluate the function again
        f1 = self._f(**kwarg)
        # Restore the property to its original value
        kwarg[self.param] = x0
        # Calculate the derivative
        df0 = (f1-f0)/dx
        return f0,df0

    def _fdf_f_df(self, **kwarg):
        """Provide f() and df() from f and df"""
        # This is just a wrapper function
        return self._f(**kwarg), self._df(**kwarg)


    def prop_name(self):
        if self._prop_name is not None:
            return self._prop_name
        elif hasattr(self._f,'__name__'):
            return self._f.__name__
        elif hasattr(self._f,'__func__') and hasattr(self._f.__func__,'__name__'):
            return self._f.__func__.__name__
        else:
            raise pyro.utility.PMAnalysisError(
"""Could not determine the property method/function name.""")
