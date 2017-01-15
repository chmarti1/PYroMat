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


class _vectorize_args(object):
    """A wrapper iterator for solver target values and keyword arguments

>>> Y = [1,2,3]
>>> kwarg = {'a':[1,2,3], 'b':[4,5,6], 'c':7}
>>> varg = _vectorize_args(Y,kwarg)
>>> for thisy,thisarg in varg:
...     print(thisy,thisarg)
...
8, {'a':1, 'b':4, 'c':7}
9, {'a':2, 'b':5, 'c':7}
10, {'a':3, 'b':6, 'c':7}

The 'ignore' keyword is an optional list of string key values to ignore
when iterating.  They will still be passed, but without iteration.

>>> varg = _vectorize_args(Y, kwarg, ignore=['b'])
>>> for thisy,thisarg in varg:
...     print(thisy,thisarg)
...
8, {'a':1, 'b':[4,5,6], 'c':7}
9, {'a':2, 'b':[4,5,6], 'c':7}
10, {'a':3, 'b':[4,5,6], 'c':7}
"""
    def __init__(self,y,args,ignore=[], pop=[]):
        # The original dictionary
        self.args=args
        self.y = y
        self.copy={}
        # We're iterating; where are we in the lists
        self.index = -1
        # This is a list of keys over which to iterate
        self.iterkeys = set(args.keys())
        # Is Y an iterable too?
        self.itery = False
        ignore = set(ignore)
        # What is the length of the iteration?
        self.N = 1
        # Check the y length; is it iterable?
        if hasattr(y,'__len__') and hasattr(y,'__iter__'):
            keylen = len(y)
            if keylen==1:
                self.itery = False
            else:
                self.N = keylen
                self.itery = True
        else:
            self.N = 1
            self.itery = False
        # Check the keys one-by-one.  Are they iterable?
        for key in self.iterkeys:
            if hasattr(args[key],'__len__') and hasattr(args[key],'__iter__'):
                keylen = len(args[key])
                if keylen==1:
                    ignore.add(key)
                elif self.N == 1:
                    self.N = keylen
                elif self.N != len(args[key]):
                    raise pyro.utility.PMParamError(
                    'Iterable length for "' + key + '" does not match other '+
                    'iterables')
            else:
                ignore.add(key)

        self.iterkeys -= ignore
        # Initialize the iterator copy with the ignored keys we won't be 
        # touching
        for key in ignore:
            if key in args:
                self.copy[key] = args[key]


    def __len__(self):
        return self.N

    def __iter__(self):
        self.index=-1
        return self

    def __next__(self):
        self.index+=1
        if self.index>=self.N:
            raise StopIteration
        for key in self.iterkeys:
            self.copy[key] = self.args[key][self.index]
        if self.itery:
            return self.y[self.index], self.copy
        else:
            return self.y, self.copy

    next = __next__


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

"""
    def __init__(self, 
            epsilon=1e-6, small=1e-10, max_iter=30,
            defaults={}):
        # Set up the defaults
        self.defaults = defaults
        self.epsilon = epsilon
        self.small = small
        self.max_iter = max_iter



class solve1n(_proto_solver_):
    """SOLVE1N  a function-like object for inverting property functions
    solver1n = solve1n( param, f=property_function, ...)

SOLVE1N objects are function-like objects that can be configured to 
iteratively invert PYroMat property methods.  This is an example solver
to calculate temperature ('T') from enthalpy ('h') for nitrogen.

>>> N2 = pyro.get('N2')
>>> T_h = pyro.solve.solve1n('T', f=N2.h, param_init=1000.)
>>> h = N2.h(543.)      # Calculate enthalpy at T=543K
>>> T_h(h)              # Go backwards
543.00020136574233

The 'param' argument is the only strictly mandatory argument.  It 
specifies the string name of the property's parameter that is to be
calculated.  In the example above, the SOLVE1N object, 'T_h', will 
iteratively call N2.h(T=..) because 'f' specified the N2.h method, and
param was set to 'T'.

The 'f', 'df', and 'fdf' parameters accept callable objects (methods 
and functions) for evaluating a property (f), its derivative with 
respect to the parameter of interest (df), or both at once (fdf).  More 
information is below.

Either 'param_init' or 'param_lim' is required.  'param_init' either 
specifies an initial guess for the iteration or a callable (function or
method) that returns the initial guess.  'param_lim' is used to specify
the upper and lower bound on legal values of the parameter.  If it is a 
two-element iterable, it will be of the form ['min', 'max'].  If 
'param_init' is specified as a callable (function or method), then it 
will be expected to return an interable with the same format.

All keywords passed when a SOLVE1N object are passed directly to the 
property function along with the parameter being iterated.  In the above
example using nitrogen,
>>> T_h(h=500., p=2.)
>>> T_h(500., p=2.)
are both valid calls to 'T_h' that will result in iterative calls
    N2.h(T=param_value, p=2.)

Any calls to param_init() or param_lim() callables use the same keywords
but with the parameter removed.  In the above example,
    param_init(p=2.)
    param_lim(p2.)
would be used to generate the initial value and the limits on the 
iteration.

There are three modes of operation for the SOLVE1N iterators:
[f only] With the property function only.  
    This is the simplest to set up, but it tends to have poorer 
    computational performance.  In this mode, a property function is 
    specified, and its derivative is estimated (for the Newton 
    iteration algorithm) by a small perturbation method.  Below is an 
    example for configuring a temperature-from-enthalpy function.
    The above nitrogen example uses this arrangement.

[f and df] With the property function and its derivative
    In this mode, the property function or method is specified along 
    with a separate function or method for calculating its derivative
    with respect to the parameter of interest.  In the nitrogen example,
    this would take the form:

>>> T_h = pyro.solve.solve1n('T', f=N2.h, df=N2.cp, param_lim=[300., 3000.])

    'df' should be a callable (function or method) that accepts all the
    same keywords and arguments as 'f'

[fdf] A function that returns both the property and its derivative
    There are many cases where this will be the most efficient method.
    When 'fdf' is specified, neither 'f', nor 'df' should be specified.
    'fdf' should return a tuple whose first value is the property and 
    the second value should be the derivative with respect to the 
    property.  For the nitrogen example, it might appear:

>>> def n2fdf(**kwarg):
...     return N2.h(**kwarg), N2.cp(**kwarg)
...
>>> T_h = pyro.solve.solve1n('T', fdf=n2fdf, param_init=1000.)

There are a number of keyword options that can be passed to the 
initializer to configure the behavior of the SOLVE1N object.
    prop_name   is a string indicating the name of the property being
                inverted.  By default it is None, and the name will be
                generated by the prop_name() function from the string 
                name of the function passed to 'f' or 'fdf'.  This 
                value is stored in the '_prop_name' attribute.  Code
                should never access the '_prop_name' attribute directly,
                but through the prop_name() method.

A number of the behaviors of the solver can be manipulated by changing 
the member attributes manually or through keyword arguments to the
initializer.  For example:
    epsilon     is 1e-6 by default, and indicates the fractional 
                precision to which the property should be evaluated.
    small       is 1e-10 by default, and indicates a number that is 
                numerically small enough to be effectively zero.
    max_iter    indicates the maximum number of iteration steps allowed
                before raising a failure to converge exception. This is
                set to 30 by default.
    defaults    is a dictionary of keywords and values to be passed on
                to functions regardless of whether they are specified 
                when the SOLVE1N object is called.  In the nitrogen 
                example, "defaults={'p':1.01325}" passes atmospheric
                pressure when 'p' is not explicitly defined by a call to
                'T_h'.

There are also some "undocumented" attributes that might come in handy
for development, but that shouldn't be relied on for exported code, and
that are not exposed through the initializer.
    _verbose    is False by default, but when True, the solver will 
                print an iteration history table to stdout.
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
        self._verbose = False
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
                    self.param_init = float(param_init)
                except:
                    raise pyro.utility.PMAnalysisError(
                    "A SOLVE1N object received an illegal value for param_init.")

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

        # We will let "thisarg" be a dictionary of arguments to pass to the
        # property funciton.  It will include the iteration parameter, 
        # self.param, the solver defaults, self.defaults, and any extras the
        # user passed in the function call.
        # Let x be the value of the parameter under iteration.
        # Let y be the target value for self._f

        # Start with the defaults and override with any new arguments
        args = self.defaults.copy()
        args.update(kwarg)

        propname = self.prop_name()
        if propname in args:
            if y is not None:
                raise pyro.utility.PMAnalysisError(
                "The property value '{:s}' was specified twice.".format(propname))
            y = args.pop(propname)
        elif y is None:
            raise pyro.utility.PMAnalysisError(
            "The property value '{:s}' was not specified.".format(propname))

        # Iterate over any iterable arguments
        arg_iter = _vectorize_args(y,args)
        # initialize the result
        done = np.zeros((len(arg_iter),))
        
        for thisy,thisarg in arg_iter:
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
                    x = float(self.param_init)
            else:
                x = np.average(limits)

            # Establish the error threshold for convergence
            threshold = max(self.small,abs(self.epsilon*thisy))

            # Define fake old x and error values
            aerror_old = float('inf')
            x_old = float(x)
            dx = 0.

            # Start the iteration
            for count in range(self.max_iter):
                # Make an initial call to the function
                thisarg[self.param] = x
                f,df = self._fdf(**thisarg)
                error = f - thisy
                if self._verbose:
                    if count%10==0:
                        print("{:4s}{:>15s}{:>15s}{:>15s}{:>15s}".format(
                        'step','x','f','df','error'))
                    print("{:4d}{:15.4e}{:15.4e}{:15.4e}{:15.4e}".format(
                        int(count),float(x),float(f),float(df),float(error)))

                aerror = np.abs(error)
                # Test for convergence
                if aerror < threshold:
                    done[arg_iter.index] = x
                    break

                # Catch a bad guess - if the error has increased,
                # don't keep going; back-track
                # If aerror evaluates to None or NAN, the comparison will always
                # fail.  Using the logical not operation catches that exception too.
                if not (aerror < aerror_old):
                    dx /= 2.
                    x = x_old + dx
                    if self._verbose:
                        print("  Recovering from bad guess")
                else:
                    # If the guess is good, project the next x value and shift "old"
                    # values.
                    dx = -error / df
                    # Force floats in case numpy or integer types show up somewhere
                    x_old = float(x)
                    aerror_old = float(aerror)
                    x += dx
                    # enforce the limits
                    if limits is not None:
                        x = max(min(x,limits[1]), limits[0])
            if count>=self.max_iter:
                raise pyro.utility.PMAnalysisError(
                    "The solution failed to converge " + repr(self))

        if len(arg_iter) == 1:
            done = np.float(done)
        return done


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
