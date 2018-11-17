# MP1
#   PYroMat Multi-phase generalist class
#   Calculates physical properties from a fit for the helmholtz free 
#   energy in terms of density and temperature.

import numpy as np
import pyromat as pm



class mp1(pm.reg.__basedata__):
    """The PYroMat multi-phase generalist class 1

Provides property methods:
    cp()    Isobaric specific heat
    cv()    Isochoric specific heat
    gam()   Specific heat ratio
    e()     Internal energy
    h()     Enthalpy
    s()     Entropy
    hsd()   Enthalpy, entropy, and density

Provides equation-of-state methods:
    T()     Temperature
    p()     Pressure
    d()     Density
    
All of the above methods accept a standardized call signature, which 
accepts any of the following arguments:
    T       Temperature
    p       Pressure
    d       Density
    x       Quality
For example, enthalpy might be called
    h(T,p)
    h(T,d)
    h(T,x)
or any other combination of two.  See the method documentation for more
details.
    
There are also saturation property methods:
    es()    Saturation internal energy
    hs()    Saturation enthalpy
    ss()    Saturation entropy
    
And saturation equations of state methods:
    Ts()    Saturaiton temperature
    ds()    Saturation density
    ps()    Saturation pressure
    
All saturation properties accept T OR p as arguments.  See the method 
documentation for more details.

*** MORE DOCUMENTATION ***
MP1 models thermo-physical properties of a liquid-gas system using a 
general fit for helmholtz free energy.  These "Spand & Wagner" fits are 
evaluated in a polynomial form with exponential post factors.

The MP1 class is divided into three layers of methods (routines).  

--- USER ROUTINES ---
Accept data in any format (array or scalar) and in whatever units are
configured in the PYroMat configuration object.  These routines rely on
_argparse and _sat_argparse to standardize their call signatures, to
convert to the correct units, and to enforce that all inner routines
receive correctly broadcast ndarray objects.

Values from these methods are returned in appropriately broadcast arrays
in the correctly configured units.

--- INNER ROUTINES ---
These methods presume that all arguments are numpy arrays and that they
are in a common unit system.  This prevents redundant calls to the unit
conversion functions as MP1 methods call one another.
    Energy:     J
    Matter      kg
    Pressure:   Pa
    Temperature:K
    
VERY rarely, should these routines by called by the user.  They are 
faster than the user routines because they do not have the overhead of
unit conversions, array broadcasting, and call signature conversion, but
they have stringent requirements on the format of data.  Users should
beware.

1) All arguments must be a numpy NDARRAY object of dimension 1 or 
    greater.
2) Array broadcasting must be done BEFORE passing arguments to the inner
    routines.
3) The above units MUST be respected regardless of PYroMat's settings.
4) Many of these functions also return their derivatives to facilitate
    numerical inversion.  Check the documentation to verify the call
    signature of each inner routine BEFORE implementing it in your code.

--- PRIMATIVE ROUTINES ---
Methods that have been labeled as primative routines should UNDER NO
CIRCUMSTANCES be called by the user.  They accept non-dimensionalized
arguments and return non-dimensional parameters.  These are encapsulated
as independent methods either because they are complicated and need to 
be called by a number of other methods, or because separating them made
sense for numerical efficiency.  In summary: these aren't the methods
you're looking for.

--- DATA DICTIONARY ---
The MP1 data dictionary must have certain data "groups" to define the 
various empirical fits.  Each group is a dictionary (within the 
dictionary) that defines the various parameters necessary for at least
one of the inner methods.

PSgroup         Saturated pressure data group
    Tscale      Temperature scale for normalizing T in the fit
    pscale      Pressure scale for re-scaling the result
    coef        a coefficient group to be passed to _poly1()
If tt = T/Tscale
    tt*log(ps/pscale) = p(1-tt)
where p() is the polynomial defined by the coef list.

DSLgroup        Saturated liquid density data group
    Tscale      Temperature scale for normalizing T in the fit
    dscale      Density scale for re-scaling the result
    coef        a coefficient group to be passed to _poly1()
If tt = T/Tscale
    dsl / dscale = p(1-tt)
where p() is the polynomial defined by the coef list.

DSVgroup        Saturated vapor density data group; a dict containing:
    Tscale      Temperature scale for normalizing T in the fit
    dscale      Density scale for re-scaling the result
    coef        a coefficient group to be passed to _poly1()
If tt = T/Tscale
    log(dsv/dscale) = p(1-tt)
where p() is the polynomial defined by the coef list.

AOgroup        Helmholtz free energy ideal gas group; a dict containing:
    Tscale      Temperature scale for normalizing T
    dscale      density scale for normalizing d
    logt        a scalar coefficient of a log(tt) term
    coef        a coefficient list to be passed to _poly1()
If tt = Tscale/T    <=== INVERSE!
and dd = d/dscale
    ao = log(d) + LOGT*log(tt) + p(tt)
    Ao = ao * R * T
where LOGT is the coefficient defined by the 'logt' parameter, and p is
the polynomial defined by the coef list

ARgroup         Helmholtz free energy residual group; a dict containing:
    Tscale      Temperature scale for normalizing T
    dscale      density scale for normalizing d
    coef0       a nested list of coefficient lists
    coef1       an optional nested list of coefficients
    coef2       an optional nested list of coefficients

The Tscale and dscale are used to non-dimensionalize temperature and 
density.
tt = Tscale/T    <=== INVERSE!
dd = d/dscale

Each element of coef0 is, itself a coefficient list intended to be 
passed to _poly2().  After the first element, each individual polynomial
is multiplied by exp(-dd**k) where k is the index in the coef list.
    ar0 = p0(tt,dd) + exp(-dd)*p1(tt,dd) + exp(-dd**2)*p2(tt,dd) + ...
    Ar0 = ar0 * R * T
    
coef1 is an optional list of lists of coefficients forming a matrix
[...
    [ c, d, t, a, b, gam, ep ], ...
]
    ar1 = c * dd**d * tt**t * exp(-a*(dd-ep)**2 - b*(tt-gam)**2) + ...
    Ar1 - ar1 * R * T
If coef1 is defined it will be combined with the other coefficients
to form the residual.  If coef1 is not defined, it will be ignored.

coef2 is an optional list of lists of coefficients forming a matrix
[...
    [ c, a, b, m, A, B, C, D ], ...
]
    X = ((1-tt) + A*((dd-1)**2)**(0.5/m))**2 + B*((dd-1)**2)**a
    ar2 = c * X**b * d * exp(-C*(dd-1)**2 - D*(tt-1)**2) + ...
    Ar2 - ar2 * R * T

Additionally, there are a number of parameters that define static 
properties

Tlim            A two-element list of the upper and lower temperatures
                for which the data set is valid.
plim            A two-element list of the upper and lower pressures for
                which the data set is valid.
dlim            A two-element list the represent practical maximum and
                minimum densities over the entire data set.  These are 
                NOT guaranteed limits of validity.
Tc, pc, dc      Critical temperature, pressure, and density
Tt, pt          Triple-point temperature and pressure
R               Ideal gas constant 8.314 / mw
mw              Molecular weight
contents        A dictionary with a key for each atom and a value for 
                its count in the molecule.  For example, CO2 would 
                have content = {'C':1, 'O':2}
                
There are also the typical mandatory PYroMat meta data elements:
id              What substance is this?
doc             Where did it come from?
class           What class should be used to evaluate the data?
"""
        


    def _poly2(self,x,y,group,diff=2):    
        """Polynomial evaluation
(p, px, py, pxx, pxy, pyy) = _poly(x,y,coef,diff=2)

Evaluates a polynomial on x and y and its derivatives.
x       x value
y       y value
param   coefficient list
diff    the highest order derivative to evaluate (0,1, or 2)

Returns
p       polynomial value at p(x,y)
px      dp/dx
py      dp/dy
pxx     d2p/dx2
pxy     d2p/dxdy
pyy     d2p/dy2

The behavior of poly2 is very much the same as poly1, but for funcitons
of two variables.  The pre- and post- exponents for poly2 are expected
to be lists or tuples: [prex, prey], [postx, posty]
The innermost lists defining the terms must contain three elements:
[powxN, powyN, coefN].  The coefficients must be sorted by x-power and
then by y-power in descending order.

A coefficient list might appear

coef = [
    [
        [prex, prey], 
        [postx, posty],
        [powxN, powyN, coefN],
        ...
        [powx0, powy0, coef0]
    ],
    [
        [prex, prey], 
        [postx, posty],
        [powxN, powyN, coefN],
        ...
        [powx0, powy0, coef0]
    ]
]
    
The pre-exponents are applied to the arguments to the polynomial, and
the post-exponents are applied after the polynomial is evaluated, so 
that the value returned is
    x**postx * y**posty * p(x**prex, y**prey)

Starting with the third element (element 2) of the coefficient list,
each element of coef is a three-element list defining a term in the 
polynomial; the x-exponent, the y-exponent, and the corresponding
coefficient.  It must be sorted in descending order by the first column
and then the second column.

For example, the list,
[[[1,1], [0,0], [1, 1, 0.1], [0, 2, 0.2], [0, 1, 1.2], [0, 0, 0.5]]]
corresponds to the polynomial
p(x,y) = .5 + 1.2y + .2y**2 + 0.1xy

Efficient polynomial evaluation algorithms are normally restricted to
positive integer exponents, but many thermodynamic property models use 
much more interesting polynomials.  The pre- and post- exponents can be
used to acheive a much wider range of functions.

For example,
    p(x,y) = x**(-1.5) + x**(3.5)
might be expressed as a coefficient list
    [[[0.5,1], [-1.5, 0], [10, 0, 1], [0, 0, 1]]]
The pre- and post- exponents make this equivalent to
    xx = x**0.5
    p(x,y) = x**(-1.5) (xx**10 + 1)
which is equivalent to the original polynomial, except that the core of
the evaluation algorithm only operates on positive integers.
"""
    
        g = 0.  # total group
        gx = 0.
        gy = 0.
        gxx = 0.
        gxy = 0.
        gyy = 0.

        for coef in group:
            # initialize the final polynomial and its derivatives
            p = 0.  # total polynomial
            px = 0.
            py = 0.
            pxx = 0.
            pxy = 0.
            pyy = 0.
            
            # collect the pre- and post-exponentials
            prex,prey = coef[0]
            postx,posty = coef[1]
            
            # Apply the pre-exponentials
            if prex!=1.:
                x_0 = x**prex
                if diff>0:
                    x_1 = x_0*prex/x
                else:
                    x_1 = 0.
                if diff>1:
                    x_2 = x_1*(prex-1.)/x
                else:
                    x_2 = 0.
            else:
                x_0 = x
                x_1 = 1.
                x_2 = 0.
                
            if prey!=1.:
                y_0 = y**prey
                if diff>0:
                    y_1 = y_0*prey/y
                else:
                    y_1 = 0.
                if diff>1:
                    y_2 = y_1*(prey-1.)/y
                else:
                    y_2 = 0.
            else:
                y_0 = y
                y_1 = 1.
                y_2 = 0.

            # From here, we loop over terms of the form a*(x**ii)*(y**jj)
            # If a particular ii,jj combination is not found in the data, then
            # its coefficient is treated as zero.
            # What is the largest ii?
            II = coef[2][0]
            
            # On which coefficient are we currently operating?
            index = 2
            # This is a flag that indicates the active index was used in 
            # the last loop, so it needs to be incremented.
            
            for ii in range(II,-1,-1):
                # If the current x-exponent is the same one represented in
                # the active coefficient row, then calculate q.
                if index<len(coef) and coef[index][0] == ii:
                    # For this value of ii, what is the largest jj?
                    JJ = coef[index][1]
                    # q is a sub-polynomial on y that represents the 
                    # variation on y of all terms that share the same
                    # power in x.  This inner loop is much like the loop
                    # on ii, except that it looks at both the x and y 
                    # exponents.
                    q = 0
                    qy = 0
                    qyy = 0
                    for jj in range(JJ,-1,-1):
                        if diff > 1:
                            qyy = 2*qy + y_0*qyy
                        if diff > 0:
                            qy = q + y_0*qy
                        # If the current y-exponent is represented in the 
                        # active coefficient row, then fold it into the q
                        # expansion.
                        if index<len(coef) and coef[index][0] == ii and coef[index][1] == jj:
                            q = coef[index][2] + y_0*q
                            # increment the active index
                            index += 1
                        else:
                            q *= y_0
                        
                    # Fold the current q values into the p expansion
                    # Update the highest derivatives first since they depend
                    # on the historical values of the lower derivatives
                    if diff > 1:
                        pyy = qyy + x_0 * pyy
                        pxx = 2*px + x_0 * pxx
                        pxy = py + x_0 * pxy
                    if diff > 0:
                        px = p + x_0 * px
                        py = qy + x_0 * py
                    p = q + x_0 * p
                # If the current x exponent is not represented, execute a
                # p-expansion with zero q.
                else:
                    if diff > 0:
                        if diff > 1:
                            pyy = x_0 * pyy
                            pxx = 2*px + x_0 * pxx
                            pxy = py + x_0 * pxy
                        px = p + x_0 * px
                        py = x_0 * py
                    p = x_0 * p
                    
            # Modify the derivatives for the pre-exponnetials
            if prex!=1.:
                if diff>0:
                    if diff>1:
                        pxx = pxx*x_1*x_1 + px*x_2
                        pyy = pyy*y_1*y_1 + py*y_2
                        pxy = pxy*x_1*y_1
                    px *= x_1
                    py *= y_1
                
            # Apply the post-exponentials
            if postx!=0:
                f = x**postx
                if diff>0:
                    fx = postx*f/x
                    if diff>1:
                        fxx = fx*(postx-1)/x
                        pxx = pxx*f + 2.*px*fx + p*fxx
                        pyy = pyy*f
                        pxy = pxy*f + py*fx
                    px = px*f + p*fx
                    py = py*f
                p *= f
            if posty!=0:
                f = y**posty
                if diff>0:
                    fy = posty*f/y
                    if diff>1:
                        fyy = fy*(posty-1)/y
                        pyy = pyy*f + 2.*py*fy + p*fyy
                        pxx = pxx*f
                        pxy = pxy*f + px*fy
                    py = py*f + p*fy
                    px = px*f
                p *= f
            
            # If the group has only one coefficient set, just return
            if len(group) == 1:
                return p,px,py,pxx,pxy,pyy
                
            g += p
            if diff>0:
                gx += px
                gy += py
                if diff>1:
                    gxx += pxx
                    gxy += pxy
                    gyy += pyy
                    
        return g,gx,gy,gxx,gxy,gyy


    def _poly1(self,x,group,diff=2):    
        """Polynomial evaluation
(p, px, pxx) = _poly1(x,coef,diff=2)

Evaluates a polynomial on x and y and its derivatives.
x       x value
coef    coefficient list
diff    the highest order derivative to evaluate (0,1, or 2)

Returns
p       polynomial value at p(x)
px      dp/dx
pxx     d2p/dx2

When diff is less than 2, the corresponding values of px and pxx are
returned as 0.  The default for diff is 2 to protect against careless
treatment as if these values ARE zero, but reducing diff will make poly1
execute more efficiently.

The coefficient list represents a nested list structure defining a poly-
nomial series.  The inner-most lists contain two elements, specifying
an integer exponent and a coefficient: [power, coefficient] for each
term.  They are contained in a list that represents groups of terms.
The groups must be sorted by power from highest to lowest.

Each group is lead by two elements that define a pre- and post- 
exponents.
    [pre, post, [powN, coefN], ... , [pow0, coef0]]

This defines a polynomial of the form
    x**post * p(x**pre)
    
The powers must be integers, but no such restriciton exists on the pre-
and post- exponents.  This permits efficient evaluaiton of polynomials
with rational exponents.

The highest level list contains a list of these groups, so that separate
pre- and post- exponentials may be applied to certain terms of the 
polynomial.

coef = [
    [
        pre, 
        post,
        [
            [powN, coefN],
            ...,
            [pow0, coef0]
        ]
    ],
    [
        pre, 
        post,
        [
            [powN, coefN],
            ...,
            [pow0, coef0]
        ]
    ]
]

In a simple example, the polynomial,
    p(x) = 2*x**-1.5 - x**0.5
    
might be specified
[[  0.5, -1.5, [0, 2.], [4, -1.]]]
"""
        g = 0.
        gx = 0.
        gxx = 0.
    
        for coef in group:
            # initialize the final polynomial and its derivatives
            p = 0.  # total polynomial
            px = 0.
            pxx = 0.

            # collect the pre- and post-exponentials
            pre = coef[0]
            post = coef[1]
            
            # Apply the pre-exponentials
            if pre!=1.:
                x_0 = x**pre
                if diff>0:
                    x_1 = x_0*pre/x
                else:
                    x_1 = 0.
                if diff>1:
                    x_2 = x_1*(pre-1.)/x
                else:
                    x_2 = 0.
            else:
                x_0 = x
                x_1 = 1.
                x_2 = 0.

            # From here, we loop over terms of the form a*(x**ii)
            # If a particular ii,jj combination is not found in the data, then
            # its coefficient is treated as zero.
            # What is the largest ii?
            II = coef[2][0]
            
            # On which coefficient are we currently operating?
            index = 2
            # This is a flag that indicates the active index was used in 
            # the last loop, so it needs to be incremented.
            
            for ii in range(II,-1,-1):
                # If the current x-exponent is the same one represented in
                # the active coefficient row, then calculate q.
                if index<len(coef) and coef[index][0] == ii:
                    # Fold the current coefficient into the p expansion
                    # Update the highest derivatives first since they depend
                    # on the historical values of the lower derivatives
                    if diff>0:
                        if diff > 1:
                            pxx = 2*px + x_0 * pxx
                        px = p + x_0 * px
                    p = coef[index][1] + x_0 * p
                    index += 1
                # If the current x exponent is not represented, execute a
                # p-expansion with zero q.
                else:
                    if diff > 0:
                        if diff > 1:
                            pxx = 2*px + x_0 * pxx
                        px = p + x_0 * px
                    p = x_0 * p
                    
            # Modify the derivatives for the pre-exponnetials
            if pre!=1.:
                if diff>0:
                    if diff>1:
                        pxx = pxx*x_1*x_1 + px*x_2
                    px *= x_1
                
            # Apply the post-exponentials
            if post!=0:
                f = x**post
                if diff>0:
                    fx = post*f/x
                    if diff>1:
                        fxx = fx*(post-1)/x
                        pxx = pxx*f + 2.*px*fx + p*fxx
                    px = px*f + p*fx
                p *= f

            if len(group)==1:
                return p,px,pxx

            g += p
            if diff>0:
                gx += px
                if diff>1:
                    gxx += pxx
                    
        return g,gx,gxx
    
    
    def _iter1(self, fn, prop, y, x, Ids, xmin, xmax,
                ep=1e-6, Nmax=20, fx_index=1, verbose=False,
                param={}):
        """Invert an inner routine.
        
    _iter1(fn, prop, y, x, Ids, xmin, xmax,)
    
Iteration is performed in-place on the x array.

*** Required Parameters ***
fn          The inner routine (method) to be inverted.  It must have a 
            call signature 
                f, fx0, ... = fn(x0, x1, ..., diff)
            where f is the value of fn, and fx0 is the derivative of fn
            with respect to prop0. The fx_index keyword can be used to
            change where fx is found in the returned tuple.  By default
            it is 1.
prop        The string keyword index of the property to be calculated.
y           An array of N target values for fn().
x           The result array.  It should be an N-element floating point
            array that has already been initialized.
Ids         A down-select boolean index array; only x[Ids],y[Ids] will 
            be evaluated.  This allows iteration in-place on data sets 
            where only a portion of the data require iteration.  If y is
            a floating point array with N elements, Ids must be a bool
            array with N elements.  It will specify a down-selected 
            data set with M elements, where M<=N.
xmin, xmax  Upper and lower limit arrays for the x values.  These must
            broadcastable to match x and y.  Even values outside of the
            down-select region should have legal values.
*** Optional Parameters ***
ep          Epsilon; fractional error permitted in y (default 1e-6)
Nmax        Maximum number of iterations (default 20)
fx_index    The location of the property derivative in the call 
            signature (default 1)
param       A dicitonary of keyword arguments are passed directly to the 
            inner routine being inverted.

"""
        # As the iteration progresses, the number of True elements in 
        # Ids will decrease until they are all false
        # There are some important intermediate values that will also
        # require indexable arrays
        dx = np.zeros_like(y, dtype=float)
        error = np.zeros_like(dx, dtype=float)
        IooB = np.zeros_like(Ids, dtype=bool)

        if verbose:
            print('Iterating on "' + prop + '"')
            print('Target values:')
            print(y)
            print('Limits:')
            print(xmin,xmax)
            print("GO!")
            print('x', 'yvalue', 'dydx', 'dx', 'Ids')

        arg = param.copy()
        count = 0
        while Ids.any():
            # Build the new argument list
            for k,v in param.items():
                # For any array arguments, shrink them along with Ids
                if isinstance(v,np.ndarray):
                    arg[k] = v[Ids]
            # Shrink the primary property array
            arg[prop] = x[Ids]
            # Evaluate the funciton and isolate its derivative
            FF = fn( diff=1, **arg)
            yy = FF[0]
            yyx = FF[fx_index]
            # note that x[Ids], yy, yyx, and all the other floating 
            # intermediates are now in m-space; the sub-set of values
            # still under iteration.
            # Calculate the error, the linear change in x, and the new x
            error[Ids] = y[Ids] - yy
            dx[Ids] = error[Ids] / yyx
            if verbose:
                print(x, yy, yyx, dx, Ids)
            x[Ids] += dx[Ids]
            # An out-of-bounds index
            #IooB = np.logical_or( x < xmin, x > xmax)
            IooB[Ids] = np.logical_or( x[Ids] < xmin[Ids], x[Ids] > xmax[Ids])
            count_oob = 0
            while IooB.any():
                dx[IooB] /= 2.
                x[IooB] -= dx[IooB]
                IooB[Ids] = np.logical_or( x[Ids] < xmin[Ids], x[Ids] > xmax[Ids])
                # Prevent a while-loop-trap
                count_oob += 1
                if count_oob>Nmax:
                    raise pm.utility.PMAnalysisError(
                        'iter1_() failed to produce a guess that was in-bounds')
            
            # Check the iteration convergence
            Ids[Ids] = abs(error[Ids]) > abs(ep*y[Ids])
            # Prevent a while-loop-trap
            count += 1
            if count>Nmax:                
                pm.utility.print_warning(\
                    'iter1_() failed to converge for %d elements after %d attempts'%(\
                    Ids.sum(), Nmax))
                return

            
    def _tpiter1(self, T, p, fn, diff=1):
        """T,p iterator wrapper
    _tpiter1(T,p,fn,diff=1)
    
    This wrapper function evaluates density from temperature and 
pressure and returns the fn(T,d,diff) inner routine.  

    d = self._d(T,p)
    return fn(T,d)
    
When diff is 1 (as it needs to be for _iter1 to work properly), the 
property's partial derivatives need to be shifted from constant-density
into constant temperature space.
    
For example, a call to _iter1 to calculate temperature while specifying
entropy and pressure might appear

    self._iter1( self._tpiter,  # Don't use _s, use _tpiter
        'T',                    # We want to calculate T
        svalues,                # Here are the target entropy values
        T,                      # The pre-initialized T array
        Ids,                    # The pre-initialized down-select array
        Tmin, Tmax,             # T bounds
        param={'fn':self._s, 'p':pvalues})
"""
        # Convert T,p into T,d
        d = self._d(T,p)
        # Assume a standard inner property routine call signature
        y,yt,yd = fn(T,d,diff=diff)
        # If the derivative is requested, we need to shift from constant
        # density to constant pressure.
        if diff>0:
            _,pt,pd = self._p(T,d,diff=1)
            temp = pt/pd
            yyt = yt - yd*temp
            yyp = yd + yt/temp
        # Do not support higher derivatives than 1
        return y,yyt,yyp
        

    def _ao(self, tt, dd, diff=2):
        """Dimensionless ideal gas helmholtz free energy (primative routine)
Evaluates an ideal gas equation of the form
    a = log(dd) + logt*log(tt) + p(t) + ... c log(1-exp(-theta*tt)) + ...
    
where
    dd = d / dscale
    tt = Tscale / T

In the AOgroup dictionary defined by the mp1 data, the polynomial, p,
is defined by the 'coef0' list.  This list should should be readable
by the _poly1() method.  The 'logt' constant defines the coefficient
of the log(tt) term.  

The log/exp expansion is defined by the 'coef1' list.  Each element of
'coef1' should be a two-element list or tuple containing [theta, c]. 
    
This is a PRIMATIVE ROUTINE.  The arguments must already be 
nondimensionalized, and the returned values are non-dimensionalzied.
"""
        
        # Start with the logarithmic terms
        A = np.log(dd)
        At = 0.
        Ad = 0.
        Att = 0.
        Atd = 0.
        Add = 0.
        
        if diff>0:
            At = 0
            Ad = 1./dd
            if diff>1:
                Add = -Ad/dd
                Att = 0.
                Atd = 0.
        
        if 'logt' in self.data['AOgroup']:
            A += self.data['AOgroup']['logt'] * np.log(tt)
            if diff>0:
                pt = self.data['AOgroup']['logt']/tt
                At += pt
                if diff>1:
                    Att += -pt/tt
                    
        # Move on to the polynomial expansion
        if 'coef0' in self.data['AOgroup']:
            p,pt,ptt = self._poly1(tt,self.data['AOgroup']['coef0'],diff)
            A+=p
            if diff>0:
                At += pt
                if diff>1:
                    Att += ptt
        
        # Now the nested log/exp expansion
        if 'coef1' in self.data['AOgroup']:
            for theta,c in self.data['AOgroup']['coef1']:
                e = np.exp(-theta*tt)
                p = np.log(1-e)
                A += c*p
                if diff>0:
                    pt = theta*e/(1.-e)
                    At += c*pt
                    if diff>1:
                        ptt = -pt*(theta + pt)
                        Att += c*ptt
                        
        return A, At, Ad, Att, Atd, Add


    def _ar(self, tt, dd, diff=2):
        """Dimensionless residual helmhotz free energy (primative routine)
Each fit in the group is of the form
    a = exp(-dd**k) * pk(tt, dd)
    
when dd = d / dscale, tt = Tscale / T
    
    A,Ad,At,Add,Adt,Att = _ar(tt, dd, order=2)

Returns the Helmholtz free energy and its derivatives up to diff.

This is a PRIMATIVE ROUTINE.  The arguments must already be 
nondimensionalized, and the returned values are non-dimensionalzied.
"""
        # Sparse polynomial evaluation is roughly 2x as fast as dense
        # polynomial evaluation for the R134a polynomials.  The 
        # explicitly defined algorithm is also roughly 2x as fast.  the
        # p_d() algorithm for R134a evaluated in about 80us on a 4 core
        # AMD A10-9700B R7

        ARgroup = self.data['ARgroup']

        # First evaluate the polynomial without an exponential coefficient
        A,At,Ad,Att,Atd,Add = self._poly2(tt,dd,ARgroup['coef0'][0],diff)
        
        k=0
        ddk = 1.
        for term in ARgroup['coef0'][1:]:
            p,pt,pd,ptt,ptd,pdd = self._poly2(tt, dd, term, diff)
            
            k += 1
            ddk *= dd
            e = np.exp(-ddk)
            if diff>0:
                # Calculate the derivative of exp(-dd**k) without
                # the exponential; it will be multiplied through next.
                ed = -k*ddk/dd

                if diff>1:
                    # Calculate the second derivative of exp(-dd**k) without
                    # the exponential; it will be multiplied through next.
                    edd = ed*((k-1.)/dd + ed)
                    # Multiply the exponential into all the terms
                    pdd = e*(p*edd + 2.*pd*ed + pdd)
                    ptd = e*(ptd + pt*ed)
                    ptt *= e
                    
                    Att += ptt
                    Atd += ptd
                    Add += pdd
                    
                pd = e*(p*ed + pd)
                pt *= e
                
                At += pt
                Ad += pd
                
            p *= e
            
            A += p
    
        # Evaluate AR1: c * dd**d * tt**t * exp(-a*(dd-gam)**2 - b*(tt-ep)**2)
        if 'coef1' in ARgroup:
            #This is the original table order
            #for c,d,t,a,b,gam,ep in ARgroup['coef1']:
            for t,d,b,a,gam,ep,c in ARgroup['coef1']:
                ddm1 = dd-ep
                ttm1 = tt-gam
                e = np.exp(-a*ddm1**2 - b*ttm1**2)
                p = c * dd**d * tt**t

                if diff>0:
                    pt = t*p/tt
                    pd = d*p/dd
                    et = -2*b*ttm1*e
                    ed = -2*a*ddm1*e
                    if diff>1:
                        ptt = -(t-1)*pt/tt
                        pdd = -(d-1)*pd/dd
                        ptd = d*pt/dd
                        ett = -2*b*(e + ttm1*et)
                        edd = -2*a*(e + ddm1*ed)
                        etd = -2*b*ttm1*ed
                        
                        Att += e*ptt + 2*et*pt + ett*p
                        Atd += e*ptd + ed*pt + et*pd + etd*p
                        Add += e*pdd + 2*ed*pd + edd*p
                    At += e*pt + et*p
                    Ad += e*pd + ed*p
                A += e*p
        
        if 'coef2' in ARgroup:
            for a,b,m,AA,BB,CC,DD,c in ARgroup['coef2']:
                ddm1 = dd-1
                ttm1 = tt-1
                m = 1./m
                
                # Construct the distance function terms inside-out.
                # This method allows the derivatives to be efficiently
                # constructed along with the algebra; preventing 
                # redundant power operations.
                # Start with the inner-most term, and borrow p as the
                # temporary variable for construction
                
                # p = A(dd-1)**m
                p = AA*(ddm1*ddm1)**(m/2.)
                if diff>0:
                    pt = 0
                    pd = p*m/ddm1
                    if diff>1:
                        ptt = 0
                        ptd = 0
                        pdd = pd*(m-1)/ddm1
                        
                # p = (1-tt) + A(dd-1)**m
                p -= ttm1
                if diff>0:
                    pt -= 1
                    
                # p = [(1-tt) + A(dd-1)**m]**2
                if diff>0:
                    if diff>1:
                        ptt = 2*pt*pt # + 2*p*ptt (but ptt=0)
                        ptd = 2*pt*pd # + 2*p*ptd (but ptd=0)
                        xdd = 2*pd*pd + 2*p*pdd
                    pt = 2*p*pt
                    pd = 2*p*pd
                p = p*p
                
                # p = [(1-tt) + A(dd-1)**m]**2 + B(dd-1)**2a
                # borrow e for the new term
                e = BB*(ddm1*ddm1)**a
                p += e
                if diff>0:
                    ed = 2*a*e/ddm1
                    pd += ed
                    if diff>1:
                        pdd += (2*a-1)*ed/ddm1
                        
                # e = {[(1-tt) + A(dd-1)**m]**2 + B(dd-1)**2a}**b
                # borrow e for the new term
                e = p**b
                if diff>0:
                    et = b*e*pt/p
                    ed = b*e*pd/p
                    if diff>1:
                        ett = (b-1)*et*pt/p + et*ptt/pt
                        etd = (b-1)*et*pd/p + et*ptd/pt
                        edd = (b-1)*ed*pd/p + ed*pdd/pd
                
                # p = c * dd * {[(1-tt) + A(dd-1)**m]**2 + B(dd-1)**2a}**b
                # This moves the intermediate value in e back into p
                p = c * dd * e
                if diff>0:
                    pt = c*dd*et
                    pd = c*(e + dd*ed)
                    if diff>1:
                        ptt = c*dd*ett
                        ptd = c*(et + dd*etd)
                        pdd = c*(2*ed + dd*edd)
                
                # Finally, construct the exponential
                # e = exp(-C*(dd-1)**2 - D*(tt-1)**2)
                e = np.exp(-CC*ddm1**2 - DD*ttm1**2)
                if diff>0:
                    et = -2*DD*ttm1*e
                    ed = -2*CC*ddm1*e
                    if diff>1:
                        ett = -2*DD*(ttm1*et + e)
                        etd = (-2*CC*ddm1)*et
                        edd = -2*CC*(ddm1*ed + e)
                # Finally done; add the result to A
                A += p*e
                if diff>0:
                    At += pt*e + p*et
                    Ad += pd*e + p*ed
                    if diff>1:
                        Att += ptt*e + 2*pt*et + p*ett
                        Atd += ptd*e + pt*ed + pd*et + p*etd
                        Add += pdd*e + 2*pd*ed + p*edd
                
        
        return A,At,Ad,Att,Atd,Add

    def _satfit(self, tt, fn, coef, diff=0):
        """Generic saturated property fit
    s, st, stt = _satfit(tt, fn=0, diff=0)

fn is an integer indicating which property fit form to use
1   Basic polynomial on 1-tt
    coef is interpreted by poly1
2   exp(poly(1-tt))
    coef is interpreted by poly1, and the result is passed to np.exp()
3   exp(1/tt * poly(1-tt))
    coef is interpreted by poly1, the result is multiplied by 1/tt, and
    passed to np.exp()
"""
        
        if fn == 1:
            p,pt,ptt = self._poly1(1-tt, coef, diff=diff)
            if diff>0:
                pt = -pt
        elif fn == 2:
            p,pt,ptt = self._poly1(1-tt, coef, diff=diff)
            p = np.exp(p)
            if diff>0:
                pt = -pt
                if diff>1:
                    ptt = p*(ptt + pt*pt)
                pt *= p
        elif fn == 3:
            p,pt,ptt = self._poly1(1-tt, coef, diff=diff)
            invt = 1./tt
            p*=invt
            if diff>0:
                pt = invt*(-pt-p)
                if diff>1:
                    ptt = invt*(ptt-2*pt)
            p=np.exp(p)
            if diff>0:
                if diff>1:
                    ptt = p*(ptt + pt*pt)
                pt = p * pt
        return p,pt,ptt
        
        
    def _dsv(self,T,diff=0):
        """Saturated vapor density
"""
        Tscale = self.data['DSVgroup']['Tscale']
        dscale = self.data['DSVgroup']['dscale']
        
        d,dt,dtt = self._satfit( 
                T/Tscale,
                self.data['DSVgroup']['fn'],
                self.data['DSVgroup']['coef'],
                diff)
        # Rescale 
        d *= dscale
        if diff>0:
            dscale /= Tscale
            dt *= dscale
            if diff>1:
                dtt *= dscale/Tscale
        
        return d,dt,dtt
        
        
    def _dsl(self,T,diff=0):
        """Saturated liquid density
"""
        
        Tscale = self.data['DSLgroup']['Tscale']
        dscale = self.data['DSLgroup']['dscale']
        
        d,dt,dtt = self._satfit( 
                T/Tscale,
                self.data['DSLgroup']['fn'],
                self.data['DSLgroup']['coef'],
                diff)
        # Rescale 
        d *= dscale
        if diff>0:
            dscale /= Tscale
            dt *= dscale
            if diff>1:
                dtt *= dscale/Tscale
        
        return d,dt,dtt
        
        
    def _ps(self,T,diff=0):
        """Saturation pressure (inner routine)
    ps, ps_T, ps_TT = _ps(T, diff=0)
    
Presumes temperature is in Kelvin, reports pressure in Pa
"""
        Tscale = self.data['PSgroup']['Tscale']
        pscale = self.data['PSgroup']['pscale']
        
        p,pt,ptt = self._satfit( 
                T/Tscale,
                self.data['PSgroup']['fn'],
                self.data['PSgroup']['coef'],
                diff)
        # Rescale 
        p *= pscale
        if diff>0:
            pscale /= Tscale
            pt *= pscale
            if diff>1:
                ptt *= pscale/Tscale
        
        return p,pt,ptt
        
        

    def _Ts(self,p):
        
        # Initialize the result array
        T = np.ones_like(p, dtype=float) * \
                0.5*(self.data['Tt'] + self.data['Tc'])
        T,Tmin,Tmax = np.broadcast_arrays(T, self.data['Tt'], self.data['Tc'])
        
        # Create a down-select array
        Ids = np.logical_and(
                p >= self.data['pt'],
                p <= self.data['pc'])
        # Execute the iteration
        self._iter1(
                self._ps,           # Invert the saturation pressure
                'T',                # Solve for temperature
                p,                  # such that _ps(T) = p
                T,                  # The initial T values
                Ids,                # The down-select array
                Tmin,               # Minimum at the triple temp.
                Tmax)               # Maximum at the critical temp.
        return T

        
    def _p(self, T, d, diff=0):
        """Calculate pressure from (T,d) (inner routine)
    p, pt, pd = _p
    
_p() does NOT handle cases where d is "under the dome."  _p() expects
sub-critical densities to be either purely liquid or vapor.
"""
        p = 0.
        pt = 0.
        pd = 0.

        Tscale = self.data['ARgroup']['Tscale']
        dscale = self.data['ARgroup']['dscale']
        R = self.data['R']
        # Calculate dimensionless arrays
        tt = Tscale/T
        dd = d/dscale
        # Calculate the Helmholtz free energy
        _,_,ard,_,artd,ardd = self._ar(tt,dd,diff+1)
        p = T*d*R*(1. + dd*ard)
        if diff>0:
            pt = R*d*(1 + dd*ard - tt*dd*artd)
            pd = R*T*(1 + 2*dd*ard + dd*dd*ardd)

        return p,pt,pd
        
        
    def _d(self,T,p):
        """Density iterator - calculate density from T,p (inner routine)
T and p MUST be ndarrays
"""
        # Benchmarking shows that calls to _p() with fewer than 100
        # data points are all equivalently expensive; even when 
        # utilizing only a single thread.  As a result, iterations must
        # under no circumstances be conducted in series.  This bisection
        # algorithm acts on all valid data in parallel.
        
        # Create a down-select array
        I = np.ones_like(T, dtype=bool)
        # And initialize a solution array
        d = np.zeros_like(T, dtype=float)
        # Initialize upper and lower iteration densities
        da = np.zeros_like(T, dtype=float)
        db = np.zeros_like(T, dtype=float)
        
        # Separate out sub-critical and super-critical values for 
        # initial conditions.  For temperatures that are super-critical, 
        # use the extreme density limits of the data set.
        Itest = T>=self.data['Tc']
        da[Itest] = self.data['dlim'][0]
        db[Itest] = self.data['dlim'][1]
        d[Itest] = 0.5*(self.data['dlim'][0] + self.data['dlim'][1])
        # For temperatures that are sub-critical, detect whether the 
        # state is liquid or gaseous.  Set Itest to sub-critical.  
        Itest = np.logical_not(Itest)
        if Itest.any():
            # Now, isolate the vapor points; set the upper density to the
            # saturated vapor density FORCE Istate to be an ndarray
            Istate = np.zeros_like(T, dtype=bool)
            Istate[Itest] = p[Itest] < self._ps(T[Itest], 0)[0]
            da[Istate] = self.data['dlim'][0]
            db[Istate] = self._dsv(T[Istate], 0)[0]
            d[Istate] = db[Istate] - da[Istate]
            # Move the saturation bounds by 5%
            db[Istate] += .05 * d[Istate]
            d[Istate] = 0.5*d[Istate] + da[Istate]
            # Now, isolate the liquid points; set the lower density to the
            # saturated liquid density
            Istate[Itest] = np.logical_not(Istate[Itest])
            da[Istate] = self._dsl(T[Istate], 0)[0]
            db[Istate] = self.data['dlim'][1]
            d[Istate] = db[Istate] - da[Istate]
            da[Istate] -= .05*d[Istate]
            d[Istate] = db[Istate] - 0.5*d[Istate]
        
            # Release the memory from Istate
            del Istate

        # perform the iteration
        self._iter1(
                self._p,
                'd',
                p,
                d,
                I,
                da,
                db,
                fx_index = 2,
                param={'T':T})
                
        return d
        
        
    def _T(self,d,p,sat=False):
        """Temperature iterator - calculate temperature from d,p (inner routine)
d and p MUST be ndarrays

    T = _T(d,p,sat=False)

Unlike _p(), _T() DOES handle cases where d is "under the dome."  These
calculations are relatively expensive, but they are necessary to the _T
inversion process.  When sat is set to True, these intermediate 
calculations are returned to prevent redundent saturation property calls

    T,dsL,dsV,Isat = _T(d,p,sat=True)
    
dsL and dsV are the saturation densities at p
Isat is a boolean index array that is True at points where d is between
    dsL and dsV.

Calling _T() should be avoided when possible, since it is one of the
more expensive iterators.  It requires iterative steps to calculate
the saturation properties in terms of pressure AND the EOS has to be
inverted to calculate T
"""
        # Benchmarking shows that calls to _p() with fewer than 100
        # data points are all equivalently expensive; even when 
        # utilizing only a single thread.  As a result, iterations must
        # under no circumstances be conducted in series.  This bisection
        # algorithm acts on all valid data in parallel.
        
        # Initialize a down-select array
        I = np.ones_like(d, dtype=bool)
        # Initialize a saturation index array
        Isat = np.zeros_like(I, dtype=bool)
        # Initialize a result array
        T = np.zeros_like(d, dtype=float)
        # Initialize upper and lower iteration densities
        Ta = np.zeros_like(d, dtype=float)
        Tb = np.zeros_like(d, dtype=float)
        # Saturaiton density arrays
        dsL = np.zeros_like(d, dtype=float)
        dsV = np.zeros_like(d, dtype=float)
        
        # Separate out sub-critical and super-critical values for 
        # initial conditions.  For pressures that are super-critical, 
        # use the extreme temperature limits of the data set.
        Itest = np.logical_and(p>=self.data['pc'], I)
        Ta[Itest] = self.data['Tlim'][0]
        Tb[Itest] = self.data['Tlim'][1]
        T[Itest] = 0.5*(self.data['Tlim'][0] + self.data['Tlim'][1])
        
        # For pressures that are sub-critical, detect whether the 
        # state is liquid or gaseous.  Set Itest to sub-critical.  
        Itest = np.logical_not(Itest)
        if Itest.any():
            # Now, identify the points in liquid, vapor, and mixed states
            # First, we'll need the saturation temperatures... this is 
            # a numerically expensive process since Ts() is iterative.
            # Let Ta temporarily be the saturation temperature
            Ta[Itest] = self._Ts(p[Itest])
            dsL[Itest] = self._dsl(Ta[Itest], 0)[0]
            dsV[Itest] = self._dsv(Ta[Itest], 0)[0]
        
            # Now, identify the liquid points
            Isat[Itest] = d[Itest] > dsL[Itest]
            # Shift the saturation temperature to Tb
            Tb[Isat] = Ta[Isat]
            Ta[Isat] = self.data['Tlim'][0]
            T[Isat] = Tb[Isat] - Ta[Isat]
            # Grow the boundary by 5%
            Tb[Isat] += 0.05*T[Isat]
            T[Isat] = Ta[Isat] + 0.5*T[Isat]
            
            # Now, identify the vapor points
            Isat[Itest] = d[Itest] < dsV[Itest]
            # Leave Ta as the saturation temperature
            Tb[Isat] = self.data['Tlim'][1]
            T[Isat] = Tb[Isat] - Tb[Isat]
            # Grow the boundary by 5%
            Ta[Isat] -= 0.05*T[Isat]
            T[Isat] = Tb[Isat] - 0.5*T[Isat]
            
            # Now, get the saturated states
            Isat[Itest] = np.logical_and(
                    d[Itest] >= dsV[Itest],
                    d[Itest] <= dsL[Itest])
            # We now have the solution at these points.
            # Assign the value to T
            T[Isat] = Ta[Isat]
            # Put safe values in Ta and Tb
            Tb[Isat] = self.data['Tlim'][1]
            Ta[Isat] = self.data['Tlim'][0]
            # Eliminate these from the down-select array - no iteraiton required.
            I[Isat] = False
        
        self._iter1(
                self._p,
                'T',
                p,
                T,
                I,
                Ta,
                Tb,
                param={'d':d})
        
        if sat:
            return T, dsL, dsV, Isat
        return T
        
        
    def _sat_argparse(self, T=None, p=None):
        """A standard argument parsing scheme for all user-layer saturation properties
    T,dL,dV = _sat_argparse(T=None, p=None)
    
Enforces that all returned parameters are numpy arrays with at least one
dimension.  Accepts T and p as scalars or array-like objects in 
[unit_temperature] and [unit_pressure] respectively.
    
Returns
T   the temperature in K
dL and dV are the liquid and vapor densities in kg/m3
"""
        if p is None:
            if T is None:
                T = pm.config['def_T']
            T = pm.units.temperature_scale(
                    np.asarray(T, dtype=float), 
                    to_units='K')
            if T.ndim==0:
                T = np.reshape(T, (1,))
            if (T<self.data['Tt']).any() or (T>self.data['Tc']).any():
                raise pm.utility.PMParamError(
                        'Saturation properties are not available at ' +
                        'temperatures beyond the triple or critical points.')
        elif T is None:
            p = pm.units.pressure(
                    np.asarray(p, dtype=float), 
                    to_units='Pa')
            if p.ndim==0:
                p = np.reshape(p, (1,))
            if (p<self.data['pt']).any() or (p>self.data['pc']).any():
                raise pm.utility.PMParamError(
                        'Saturation properties are not available at ' +
                        'pressures beyond the triple or critical points.')
            T = self._Ts(p)
        else:
            raise pm.utility.PMParamError(
                'Saturation temperature and pressure cannot be simultaneously specified')

        dL = self._dsl(T,0)[0]
        dV = self._dsv(T,0)[0]
        return T, dL, dV
        
        
    def _argparse(self, T=None, p=None, d=None, x=None):
        """Present a standard argument scheme for all user-layer property methods
    T,d1,d2,x,I = _argparse(T=None, p=None, d=None, x=None)

T the temperature array in K.
d1 and d2 are densities in kg/m3.  If the conditions are under the dome,
    then d1 is the liquid, and d2 is the vapor density.
x is the quality.  If conditions are not under the dome, it can be 
    ignored.
I is a boolean array whose elements will be True for conditions that are
    under the dome.  Its values will be False otherwise.

returns temperature, density and quality at which the property is to be
evaluated as independent arrays.  When the conditions are under the dome
d1 and d2 represent the liquid and vapor densities respectively.  At all 
other conditions, x<0 and d1 == d2.
"""
        # Case out the possible combinations
        # There are 11 possible pairs: Unspecified, T, p, d, x
        # 
        # 1) Convert the inputs to an array of dimension 1 or greater
        # 2) Convert the inputs into the correct units
        # 3) Broadcast the arrays appropriately
        # 4) Calculate T,d1,d2,x, and I
        
        # First, assign defaults if necessary
        # count the assigned parameters
        nparam = ((0 if T is None else 1) + 
                (0 if p is None else 1) + 
                (0 if d is None else 1) + 
                (0 if x is None else 1))
        if nparam == 1:
            if T is None:
                T = pm.config['def_T']
            else:
                p = pm.config['def_p']
        elif nparam == 0:
            T = pm.config['def_T']
            p = pm.config['def_p']
        elif nparam > 2:
            raise utility.PMParameterError(
                    'Specifying more than two simultaneous parameters is illegal.')
        
        if T is not None:
            # Make a copy of the array only if necessary
            T = pm.units.temperature_scale(
                    np.asarray(T, dtype=float), 
                    to_units='K')
            # Force T to have at least one dimension
            if T.ndim == 0:
                T = np.reshape(T, (1,))
            # Ensure that T is a legal value
            if ((T<self.data['Tlim'][0]).any() or 
                    (T>self.data['Tlim'][1]).any()):
                raise pm.utility.PMParamError('MP1: Temperature is out-of-bounds.')
            
            # T,p
            # If p is defined, then none of the conditions are under the
            # dome.  Use _d() to invert into density units.
            if p is not None:
                # Convert pressure
                p = pm.units.pressure(
                        np.asarray(p, dtype=float), 
                        to_units='Pa')
                if p.ndim==0:
                    p = np.reshape(p,(1,))
                if ((p<self.data['plim'][0]).any() or 
                        (p>self.data['plim'][1]).any()):
                    raise pm.utility.PMParamError('MP1: Pressure is out-of-bounds.')
                    
                # Force compatible arrays
                T,p = np.broadcast_arrays(T,p)
                d1 = self._d(T,p)
                d2 = d1
                x = np.broadcast_to(-1, T.shape)
                I = np.broadcast_to(False, T.shape)
                return T, d1, d2, x, I
            
            # T,d
            # If T,d then points CAN be under the dome
            elif d is not None:
                # Convert density
                d1 = pm.units.matter(
                        np.asarray(d, dtype=float), 
                        self.data['mw'], to_units='kg')
                if d1.ndim == 0:
                    d1 = np.reshape(d1,(1,))
                d1 = pm.units.volume(d1, 
                        to_units='m3', exponent=-1)
                # broadcast the arrays
                T,d1 = np.broadcast_arrays(T,d1)
                # Isolate the sub-critical temperatures
                I = np.asarray(T<self.data['Tc'])
                # Calculate the saturation densities
                dsL = self._dsl(T[I])[0]
                dsV = self._dsv(T[I])[0]
                # Identify the densities that are under the dome
                Isat = np.logical_and(
                            d1[I] < dsL, d1[I] > dsV)
                # Modify I to include only the points that are saturated
                I[I] = Isat
                # If there are any densities under the dome
                if Isat.any():
                    # Broadcasting can cause elements of d1 to refer to
                    # common locations in memory.  If there are 
                    # saturated points, we need to modify d1, so this 
                    # will force d1 to be a fully populated array.
                    d1 = d1.copy()
                    # Calculate the quality
                    x = -np.ones_like(T, dtype=float)
                    x[I] = (dsL[Isat] / d1[I] - 1.) / (dsL[Isat] / dsV[Isat] - 1.)
                    # Update the densities
                    d2 = d1.copy()
                    d1[I] = dsL[Isat]
                    d2[I] = dsV[Isat]
                else:
                    d2 = d1
                    x = np.broadcast_to(-1,T.shape)
                return T, d1, d2, x, I
                
            # T,x
            # If quality is defined, then the points MUST be saturated
            elif x is not None:
                T,x,I = np.broadcast_arrays(T,x,True)
                if (T>self.data['Tc']).any():
                    raise pm.utility.PMParamError(
                        'Quality cannot be specified above the critical temperature.')
                d1 = self._dsl(T,0)[0]
                d2 = self._dsv(T,0)[0]
                return T,d1,d2,x,I
                
            # This should never happen
            else:
                raise pm.utility.PMAnalysisError(
                        'Unhandled case in _argparse()')
        # p
        # If p is the primary parameter
        elif p is not None:
            # Convert p to the correct units
            pm.units.pressure(
                    np.asarray(p, dtype=float), 
                    to_units='Pa')
            if p.ndim==0:
                p = np.reshape(p, (1,))
            # Force p to have dimension 1 or greater
            # Ensure that p is a legal value
            if ((p<self.data['plim'][0]).any() or 
                    (p>self.data['plim'][1]).any()):
                raise pm.utility.PMParamError('MP1: Pressure is out-of-bounds.')
                
            # p,d
            # Pressure, density is an expensive combination since it
            # involves iteration to determine the saturation properties
            # AND to recover temperature.
            if d is not None:
                d1 = pm.units.matter(
                        np.asarray(d, dtype=float), 
                        self.data['mw'], to_units='kg')
                if d1.ndim==0:
                    d1 = np.reshape(d1, (1,))
                d1 = pm.units.volume(d1, to_units='m3', exponent=-1)
                # Broadcast the arrays
                d1,p = np.broadcast_arrays(d1,p)
                # This one's an expensive funciton call
                # Get temperature and the saturation densities
                T,dsL,dsV,Isat = self._T(d1,p,sat=True)
                # If there are any saturated points
                if Isat.any():
                    # Broadcasting can cause elements of d1 to refer to
                    # common locations in memory.  If there are 
                    # saturated points, we need to modify d1, so this 
                    # will force d1 to be a fully populated array.
                    d1 = d1.copy()
                    # Calculate the quality
                    x = -np.ones_like(p, dtype=float)
                    x[Isat] = (dsL[Isat] / d1[Isat] - 1.) / (dsL[Isat] / dsV[Isat] - 1.)
                    # Update the densities
                    d2 = d1.copy()
                    d1[Isat] = dsL[Isat]
                    d2[Isat] = dsV[Isat]
                else:
                    d2 = d1
                    x = np.broadcast_to(-1,T.shape)
                return T,d1,d2,x,Isat
            
            # p,x
            # If quality is defined, we are saturated
            elif x is not None:
                # Ensure that p is sub-critical
                if (p>self.data['pc']).any():
                    raise pm.utility.PMParamError('Quality cannot be specified at pressures above the critical point.')
                    
                p,x,I = np.broadcast_arrays(p,x,True)
                # T is just the saturation temperature
                T = self._Ts(p)
                d1 = self._dsl(T,0)[0]
                d2 = self._dsv(T,0)[0]
                return T, d1, d2, x, I
                                
            # This should never happen
            else:
                raise pm.utility.PMAnalysisError(
                        'Unhandled case in _argparse()')
                
        # d
        elif d is not None:
            # d,x
            # This combination is not supported!
            # This represents an expensive inversion problem that is
            # not certain to have a solution, and it is highly unusual 
            # to specify quality AND density.
            if x is not None:
                raise pm.utility.PMParamError(
                    'Specifying properties by density and quality is not currently supported.')                        
            # This should never happen
            else:
                raise pm.utility.PMAnalysisError(
                        'Unhandled case in _argparse()')
        
        raise pm.utility.PMAnalysisError(
                    'Unhandled case in _argparse()')



    def _e(self,T,d,diff=0):
        """Internal energy (inner routine)
    e,eT,ed = _e(T,d,diff=0)
"""
        eT = 0.
        ed = 0.

        # The IG part        
        R = self.data['R']
        Tscale = self.data['AOgroup']['Tscale']
        dscale = self.data['AOgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,at,_,att,atd,_ = self._ao(tt,dd,diff+1)
        
        e = at
        if diff>0:
            eT = tt*tt*att
            ed = atd
        
        # The residual part
        Tscale = self.data['ARgroup']['Tscale']
        dscale = self.data['ARgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,at,ad,att,atd,add = self._ar(tt,dd,diff+1)
        e += at
        e *= R*Tscale
        if diff>0:
            eT += tt*tt*att
            eT *= -R
            ed += atd
            ed *= R*Tscale/dscale

        return e,eT,ed


    def _h(self,T,d,diff=0):
        """enthalpy (inner routine)
    h,hT,hd = _h(T,d,diff=0)
"""
        hT = 0.
        hd = 0.

        # The IG part        
        R = self.data['R']
        Tscale = self.data['AOgroup']['Tscale']
        dscale = self.data['AOgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,at,_,att,atd,_ = self._ao(tt,dd,diff+1)
        
        h = 1. + tt*at
        if diff>0:
            hT = 1. - tt*tt*att
            hd = tt*atd
        
        # The residual part
        Tscale = self.data['ARgroup']['Tscale']
        dscale = self.data['ARgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,at,ad,att,atd,add = self._ar(tt,dd,diff+1)
        h += dd*ad + tt*at
        h *= R*T
        if diff>0:
            hT += dd*ad - tt*(tt*att + dd*atd)
            hT *= R
            hd += ad + dd*add + tt*atd
            hd *= R*T/dscale

        return h,hT,hd


    def _s(self,T,d,diff=0):
        """entropy (inner routine)
    s,sT,sd = _s(T,d,diff=0)
"""
        sT = 0.
        sd = 0.

        # The IG part        
        R = self.data['R']
        Tscale = self.data['AOgroup']['Tscale']
        dscale = self.data['AOgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        a,at,ad,att,atd,_ = self._ao(tt,dd,diff+1)
        
        s = tt*at - a
        if diff>0:
            sT = tt*tt*att
            sd = tt*atd - ad
        
        # The residual part
        Tscale = self.data['ARgroup']['Tscale']
        dscale = self.data['ARgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        a,at,ad,att,atd,_ = self._ar(tt,dd,diff+1)
        s += tt*at - a
        s *= R
        if diff>0:
            sT += tt*tt*att
            sT *= -R/T
            sd += tt*atd - ad
            sd *= R/dscale

        return s,sT,sd
        
        
    def _cp(self,T,d):
        """Isobaric specific heat (inner routine)
    cp = _cp(T,d)
"""

        # The IG part        
        R = self.data['R']
        Tscale = self.data['AOgroup']['Tscale']
        dscale = self.data['AOgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,_,att,_,_ = self._ao(tt,dd,2)
        
        cp = -tt*tt*att
        
        # The residual part
        Tscale = self.data['ARgroup']['Tscale']
        dscale = self.data['ARgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,ad,att,atd,add = self._ar(tt,dd,2)

        temp = 1.+dd*(ad-tt*atd)
        cp += -tt*tt*att + temp*temp/(1.+dd*(2.*ad+dd*add))
        cp *= R
        return cp
        
        
    def _cv(self,T,d):
        """Isochoric specific heat (inner routine)
    cv = _cv(T,d)
"""

        # The IG part        
        R = self.data['R']
        Tscale = self.data['AOgroup']['Tscale']
        dscale = self.data['AOgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,_,att,_,_ = self._ao(tt,dd,2)
        
        cv = tt*tt*att
        
        # The residual part
        Tscale = self.data['ARgroup']['Tscale']
        dscale = self.data['ARgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,_,att,_,_ = self._ar(tt,dd,2)

        cv += tt*tt*att
        cv *= -R
        return cv


    #               #
    # USER ROUTINES #
    #               #
    
    #               #
    # Data limits   #
    #               #
    def Tlim(self, p=None):
        """Return the temperature limits for the data set
    Tmin, Tmax = Tlim(p=None)
    
Tlim accepts pressure as an argument for extensibility, but the MP1 
class has homogeneous temperature limits.

Returns the temperature limits in [unit_temperature].
"""
        return pm.units.temperature_scale(
            np.asarray(self.data['Tlim']),
            from_units='K')
        
        
    def plim(self, T=None):
        """Returns the pressure limits for the data set
    pmin, pmax = plim(T=None)
    
plim accepts temperature as an argument for extensibility, but the MP1 
class has homogeneous pressure limits.

Returns the pressure limits in [unit_pressure]
"""
        return pm.units.pressure(
            np.asarray(self.data['plim']),
            from_units='Pa')
        
    #                               #
    # Critical and triple points    #
    #                               #
        
    def critical(self, density=False):
        """Critical point
    Tc, pc = critical()
    
To also return the density, set the 'density' keyword to True

    Tc, pc, dc = critical(density=True)
    
Returns the critical temperature, pressure, and density in 
[unit_temperature], [unit_pressure], [unit_matter/unit_volume]
"""
        if density:
            return  pm.units.temperature_scale( \
                        self.data['Tc'], from_units='K'),\
                    pm.units.pressure( \
                        self.data['pc'], from_units='Pa'), \
                    pm.units.volume(\
                        pm.units.matter( \
                            self.data['dc'], \
                            self.data['mw'], \
                            from_units='kg'),\
                        from_units='m3', exponent=-1)
                    
        return  pm.units.temperature_scale( \
                    self.data['Tc'], from_units='K'),\
                pm.units.pressure( \
                    self.data['pc'], from_units='Pa')
        
        
    def triple(self):
        """Triple point
    Tt, pt = triple()
    
Returns the triple temperature and pressure in a tuple pair in
[unit_temperature], [unit_pressure]
"""
        return  pm.units.temperature_scale( \
                    self.data['Tt'], from_units='K'),\
                pm.units.pressure( \
                    self.data['pt'], from_units='Pa')
        
    #                       #
    # Saturaiton properties #
    #                       #
    
    def ps(self, T=None):
        """Saturation pressure
    psat = ps(T)
    
Returns the saturaiton pressure in [unit_pressure]

Calls to ps() are MUCH faster than calls to Ts(), so when given a choice
ps() should always be preferred.  The MP1 class exposes ps() as a 
fundamental empirical relationship, while Ts() has to perform iterative 
numerical inversion.
"""
        if T is None:
            T = pm.config['def_T']

        # Replace T with an array of the correct units
        T = pm.units.temperature_scale(
                np.asarray(T, dtype=float), 
                to_units='K')
        # Exclude points outside the triple-critical range
        if np.logical_or( T<self.data['Tt'], T>self.data['Tc'] ).any():
            raise pm.utility.PMParamError(
                'Saturation properties are not ' +
                'available above the critical point Tc=%f K or below the '%self.data['Tc'] +
                'triple point Tt=%f K.'%self.data['Tt'] )

        return pm.units.pressure(self._ps(T)[0], from_units='Pa')
        
        
    def Ts(self,p=None):
        """Saturation temperature
    Tsat = Ts(p)
    
Uses Newton iteration to calculate Ts from the _ps() inner method
"""
        if p is None:
            p = pm.config['def_p']

        # Replace p with an array of the correct units
        p = pm.units.pressure(
                np.asarray(p, dtype=float), 
                to_units='Pa')
        # Force p to have at least 1 dimension
        if p.ndim==0:
            p = np.reshape(p, (1,))
        
        # Exclude points outside the triple-critical range
        if np.logical_or( p<self.data['pt'], p>self.data['pc'] ).any():
            raise pm.utility.PMParamError(
                'Saturation properties are not ' +
                'available above the critical point pc=%f bar or below the '%(self.data['pc']/1e5) +
                'triple point pt=%f bar.'%(self.data['pt']/1e5) )
        
        return pm.units.temperature_scale( \
            self._Ts(p), from_units='K')
        
        
    def ds(self, *varg, **kwarg):
        """Saturation density
    dsL, dsV = ds(T)
    
If no keyword is specified, saturation properties interpret the argument
as temperature.  However, pressure can be specified as well

    dsL, dsV = ds(p=pvalue)
    
Returns the liquid (dsL) and vapor (dsV) saturation density in units
[unit_matter / unit_volume]
"""
        _,dL,dV = self._sat_argparse(*varg, **kwarg)
        # Get a conversion factor
        conv = pm.units.matter(1., self.data['mw'],
                from_units='kg')
        conv = pm.units.volume(conv, from_units='m3', exponent=-1)
        dL *= conv
        dV *= conv
        return dL, dV


    def es(self, *varg, **kwarg):
        """Saturation internal energy
    esL, esV = es(T)

If no keyword is specified, saturation properties interpret the argument
as temperature.  However, pressure can be specified as well

    esL, esV = es(p=pvalue)
    
Returns the liquid (esL) and vapor (esV) saturation internal energy in
units [unit_energy / unit_matter]
"""
        T,dL,dV = self._sat_argparse(*varg, **kwarg)
        esL = self._e(T,dL,0)[0]
        esV = self._e(T,dV,0)[0]
        
        # Get a conversion factor
        conv = pm.units.energy(1., from_units='J')
        conv = pm.units.matter(conv, self.data['mw'],
                from_units='kg', exponent=-1)
        esL *= conv
        esV *= conv
        return esL, esV


    def hs(self, *varg, **kwarg):
        """Saturation enthalpy
    hsL, hsV = hs(T)
    
If no keyword is specified, saturation properties interpret the argument
as temperature.  However, pressure can be specified as well

    hsL, hsV = hs(p=pvalue)
    
Returns the liquid (hsL) and vapor (hsV) saturation enthalpy in
units [unit_energy / unit_matter]
"""
        T,dL,dV = self._sat_argparse(*varg, **kwarg)
        hsL = self._h(T,dL,0)[0]
        hsV = self._h(T,dV,0)[0]
        
        # Get a conversion factor
        conv = pm.units.energy(1., from_units='J')
        conv = pm.units.matter(conv, self.data['mw'],
                from_units='kg', exponent=-1)
        hsL *= conv
        hsV *= conv
        return hsL, hsV
        
        
    def ss(self, *varg, **kwarg):
        """Saturation entropy
    ssL, ssV = ss(T,p)
    
If no keyword is specified, saturation properties interpret the argument
as temperature.  However, pressure can be specified as well

    ssL, ssV = ss(p=pvalue)
    
Returns the liquid (ssL) and vapor (ssV) saturation entropy in
units [unit_energy / unit_matter / unit_temperature]
"""
        T,dL,dV = self._sat_argparse(*varg, **kwarg)
        ssL = self._s(T,dL,0)[0]
        ssV = self._s(T,dV,0)[0]
        
        # Get a conversion factor
        conv = pm.units.energy(1., from_units='J')
        conv = pm.units.matter(conv, self.data['mw'],
                from_units='kg', exponent=-1)
        conv = pm.units.temperature(conv,
                from_units='K', exponent=-1)
        ssL *= conv
        ssV *= conv
        return ssL, ssV


    #                       #
    # EOS properties T,p,d  #
    #                       #
    
    def p(self, quality=False, *varg, **kwarg):
        """Pressure
    p = p(T=None, p=None, d=None, x=None, quality=False)

Calculates the pressure in [unit_pressure]

If the optional "quality" keyword is set to True, then the two-phase
mixture quality will also be returned.

    p,x = p(...,quality=True)
    
For points that are not "under the dome" quality will be computed to be
-1.
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        # Use d2.  In theory, p(d1) = p(d2), but the liquid is so stiff
        # that small numerical errors cause huge pressure errors
        # The problem is solved when the vapor density is used instead.
        # In all other conditions d1=d2
        p = self._p(T,d2,0)[0]
        
        pm.units.pressure(p, from_units='Pa', inplace=True)
        
        if quality:
            return p,x
        return p
        
        
    def d(self, *varg, **kwarg):
        """Density
    d = d(T=None, p=None, d=None, x=None)
    
Calculates density in [unit_matter / unit_volume]
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        if I.any():
            d1[I] = (1.-x[I])/d1[I]
            d1[I] += x[I]/d2[I]
            d1[I] = 1. / d1[I]
            
        pm.units.matter(d1, self.data['mw'], from_units='kg', inplace=True)
        pm.units.volume(d1, from_units='m3', inplace=True, exponent=-1)
        return d1
        
        
    def T(self, *varg, **kwarg):
        """Temperature
    T = T(T=None, p=None, d=None, x=None)
    
Calculates temperature in [unit_temperature]
"""
        T,_,_,_,_ = self._argparse(*varg, **kwarg)
        return T
        
        
    def e(self, *varg, **kwarg):
        """Internal energy  e(T=None, p=None, d=None, x=None)
From any two of the provided primary properties
    
e   Int. Energy [unit_energy / unit_matter]
T   Temperature [unit_temperature]
p   Pressure    [unit_pressure]
d   Density     [unit_matter / unit_volume]
x   Quality     [dimensionless]
"""
        quality=False
        if 'quality' in kwarg:
            quality = kwarg.pop('quality')
            
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        e = self._e(T,d1,0)[0]
        if I.any():
            e[I] *= (1.-x[I])
            e[I] += self._e(T[I],d2[I],0)[0] * x[I]
        # Convert the units back to user space
        pm.units.energy(e, from_units='J', inplace=True)
        pm.units.matter(e, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        if quality:
            return e,x
        return e
        
        
        
    def h(self, *varg, **kwarg):
        """Enthalpy  h(T=None, p=None, d=None, x=None)
From any two of the provided primary properties
    
h   Enthalpy    [unit_energy / unit_matter]
T   Temperature [unit_temperature]
p   Pressure    [unit_pressure]
d   Density     [unit_matter / unit_volume]
x   Quality     [dimensionless]
"""
        quality=False
        if 'quality' in kwarg:
            quality = kwarg.pop('quality')
            
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        h = self._h(T,d1,0)[0]
        if I.any():
            h[I] *= (1.-x[I])
            h[I] += self._h(T[I],d2[I],0)[0] * x[I]
        # Convert the units back to user space
        pm.units.energy(h, from_units='J', inplace=True)
        pm.units.matter(h, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        if quality:
            return h,x
        return h


    def s(self, *varg, **kwarg):
        """Entropy  s(T=None, p=None, d=None, x=None)
From any two of the provided primary properties
    
    s = mp1.s( ... )

If the optional keyword "quality" is set to True, then a quality array
will also be returned

    s,x = mp1.s( ..., quality=True)

    gamma,x = mp1.gam( ..., quality=True)
s   Entropy     [unit_energy / unit_matter / unit_temperature]
T   Temperature [unit_temperature]
p   Pressure    [unit_pressure]
d   Density     [unit_matter / unit_volume]
x   Quality     [dimensionless]
"""
        quality=False
        if 'quality' in kwarg:
            quality = kwarg.pop('quality')
            
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        s = self._s(T,d1,0)[0]
        if I.any():
            s[I] *= (1.-x[I])
            s[I] += self._s(T[I],d2[I],0)[0] * x[I]
        # Convert the units back to user space
        pm.units.energy(s, from_units='J', inplace=True)
        pm.units.matter(s, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        pm.units.temperature(s, from_units='K', 
                exponent=-1, inplace=True)
        if quality:
            return s,x
        return s


    def hsd(self, *varg, **kwarg):
        """Enthalpy, Entropy, Density
    h,s,d = hsd(T=None, p=None, d=None, x=None)

If the optional keyword "quality" is set to True, then a quality array
will also be returned

    h,s,d,x = mp1.hsd( ..., quality=True)

Calculates the three most commonly used parameters at once.  This 
method represents substantial computational savings over calling the
methods independently.
"""
        quality=False
        if 'quality' in kwarg:
            quality = kwarg.pop('quality')
            
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        
        # There is no inner hsd funciton.  
        # We have to do this the hard way.
        
        # The IG part        
        R = self.data['R']
        Tscale = self.data['AOgroup']['Tscale']
        dscale = self.data['AOgroup']['dscale']
        tt = Tscale / T
        dd = d1 / dscale
        a,at,_,_,_,_ = self._ao(tt,dd,1)
        
        h = 1. + tt*at
        s = tt*at - a
        
        # The residual part
        Tscale = self.data['ARgroup']['Tscale']
        dscale = self.data['ARgroup']['dscale']
        tt = Tscale / T
        dd = d1 / dscale
        a,at,ad,_,_,_ = self._ar(tt,dd,1)
        h += dd*ad + tt*at
        s += tt*at - a

        # If there are data under the dome
        if I.any():
            temp = 1-x[I]
            h[I] *= temp
            s[I] *= temp
            
            # The IG part        
            R = self.data['R']
            Tscale = self.data['AOgroup']['Tscale']
            dscale = self.data['AOgroup']['dscale']
            tt = Tscale / T[I]
            dd = d2[I] / dscale
            a,at,_,_,_,_ = self._ao(tt,dd,1)
            
            h[I] += (1. + tt*at)*x[I]
            s[I] += (tt*at - a)*x[I]
            
            # The residual part
            Tscale = self.data['ARgroup']['Tscale']
            dscale = self.data['ARgroup']['dscale']
            tt = Tscale / T
            dd = d2 / dscale
            a,at,ad,_,_,_ = self._ar(tt,dd,1)
            h[I] += (dd*ad + tt*at)*x[I]
            s[I] += (tt*at - a)*x[I]
            # Modify density
            d1[I] = temp/d1[I] 
            d1[I] += x[I]/d2[I]
            d1[I] = 1./d1[I]
            
        s *= R
        h *= R*T
        
        conv = pm.units.energy(1.,from_units='J')
        conv = pm.units.matter(conv, self.data['mw'], from_units='kg')
        h*=conv
        conv = pm.units.temperature(conv, from_units='K')
        s*=conv
        pm.units.matter(d1, self.data['mw'],from_units='kg',inplace=True)
        pm.units.volume(d1, from_units='m3', exponent=-1, inplace=True)
        
        if quality:
            return h,s,d1,x
        return h,s,d1
        

    def cp(self, *varg, **kwarg):
        """Isobaric Specific Heat  cp(T=None, p=None, d=None, x=None)
From any two of the provided primary properties

    cp = mp1.cp( ... )

If the optional keyword "quality" is set to True, then a quality array
will also be returned

    cp,x = mp1.cp( ..., quality=True)
    
cp  Sp. heat    [unit_energy / unit_matter / unit_temperature]
T   Temperature [unit_temperature]
p   Pressure    [unit_pressure]
d   Density     [unit_matter / unit_volume]
x   Quality     [dimensionless]
"""
        quality=False
        if 'quality' in kwarg:
            quality = kwarg.pop('quality')
            
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        cp = self._cp(T,d1)
        if I.any():
            cp[I] *= (1.-x[I])
            cp[I] += self._cp(T[I],d2[I]) * x[I]
        # Convert the units back to user space
        pm.units.energy(cp, from_units='J', inplace=True)
        pm.units.matter(cp, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        pm.units.temperature(cp, from_units='K', 
                exponent=-1, inplace=True)
        if quality:
            return cp, x
        return cp


    def cv(self, *varg, **kwarg):
        """Isochoric Specific Heat  cv(T=None, p=None, d=None, x=None)
From any two of the provided primary properties

    cv = mp1.cv( ... )

If the optional keyword "quality" is set to True, then a quality array
will also be returned

    cv,x = mp1.cv( ..., quality=True)
    
cv  Sp. heat    [unit_energy / unit_matter / unit_temperature]
T   Temperature [unit_temperature]
p   Pressure    [unit_pressure]
d   Density     [unit_matter / unit_volume]
x   Quality     [dimensionless]
"""
        quality=False
        if 'quality' in kwarg:
            quality = kwarg.pop('quality')
            
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        cv = self._cv(T,d1)
        if I.any():
            cv[I] *= (1.-x[I])
            cv[I] += self._cv(T[I],d2[I]) * x[I]
        # Convert the units back to user space
        pm.units.energy(cv, from_units='J', inplace=True)
        pm.units.matter(cv, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        pm.units.temperature(cv, from_units='K', 
                exponent=-1, inplace=True)
        if quality:
            return cv, x
        return cv
        
        
    def gam(self, quality=False, *varg, **kwarg):
        """Specific Heat Ratio gam(T=None, p=None, d=None, x=None)
From any two of the provided primary properties

    gamma = mp1.gam( ... )

If the optional keyword "quality" is set to True, then a quality array
will also be returned

    gamma,x = mp1.gam( ..., quality=True)
    
gam Sp. heat ratio [dless]
T   Temperature [unit_temperature]
p   Pressure    [unit_pressure]
d   Density     [unit_matter / unit_volume]
x   Quality     [dimensionless]
"""
        quality=False
        if 'quality' in kwarg:
            quality = kwarg.pop('quality')
            
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        cv = self._cv(T,d1)
        cp = self._cp(T,d1)
        if I.any():
            cv[I] *= (1.-x[I])
            cp[I] *= (1.-x[I])
            cv[I] += self._cv(T[I],d2[I]) * x[I]
            cp[I] += self._cp(T[I],d2[I]) * x[I]
        
        if quality:
            return cp/cv, x
        return cp/cv


    def T_s(self, s, p=None, quality=False, debug=False):
        """Temperature from entropy
    T = T_s(s, p=None, quality=False)

Pressure is optional, but entropy is mandatory.

The optional keyword flag, quality, will cause quality to be returned
along with temperature.

    T,x = T_s(s, p=None, quality=True)
"""
        # Prepare the s array
        s = pm.units.energy(
                np.asarray(s,dtype=float),
                to_units='J')
        pm.units.matter(s, self.data['mw'],
                to_units='kg', exponent=-1, inplace=True)
        pm.units.temperature(s,
                to_units='K', exponent=-1, inplace=True)
        if s.ndim == 0:
            s = np.reshape(s, (1,))

        # Set a default pressure?
        if p is None:
            p = pm.config['def_p']

        p = pm.units.pressure(
                np.asarray(p, dtype=float),
                to_units='Pa')
        # Enforce pressure limits
        if ((p<self.data['plim'][0]).any() or
                (p>self.data['plim'][1]).any()):
            raise pm.utility.PMParamError(
                'MP1: Pressure is out-of-bounds.')
        if p.ndim == 0:
            p = np.reshape(p, (1,))
            
        # broadcast
        s,p = np.broadcast_arrays(s,p)
        # Initialize results
        T = np.zeros_like(s, dtype=float)
        x = -np.ones_like(s,dtype=float)
        # Some important intermediates
        Isat = np.zeros_like(s, dtype=bool)
        Ta = np.zeros_like(s, dtype=float)
        Tb = np.zeros_like(s, dtype=float)
        
        
        # Start with super-critical points
        I = p >= self.data['pc']
        Ta[I] = self.data['Tlim'][0]
        Tb[I] = self.data['Tlim'][1]
        T[I] = 0.5*(Ta[I] + Tb[I])
        # Now work with the sub-critical points
        I = np.logical_not(I)
        if I.any():
            # Get the saturation temperatures
            Tsat = self._Ts(p[I])
            # And densities
            dsL = self._dsl(Tsat,0)[0]
            dsV = self._dsv(Tsat,0)[0]
            # finally, get the saturation entropies
            ssL = self._s(Tsat,dsL,0)[0]
            ssV = self._s(Tsat,dsV,0)[0]

            # Isolate points that are liquid
            Isat[I] = s[I] < ssL
            Ta[Isat] = self.data['Tlim'][0]
            Tb[Isat] = Tsat[Isat]
            T[Isat] = 0.5*(Ta[Isat] + Tb[Isat])
            # T[Isat] = Tb[Isat] - Ta[Isat]
            # Grow the boundary by 2%
            # Tb[Isat] += 0.05*T[Isat]
            #T[Isat] = Ta[Isat] + 0.5*T[Isat]

            # Isolate points that are vapor
            Isat[I] = s[I] > ssV
            Ta[Isat] = Tsat[Isat]
            Tb[Isat] = self.data['Tlim'][1]
            T[Isat] = 0.5*(Ta[Isat] + Tb[Isat])
            # T[Isat] = Tb[Isat] - Ta[Isat]
            # Grow the boundary by 2%
            # Ta[Isat] -= 0.05*T[Isat]
            # T[Isat] = Tb[Isat] - 0.5*T[Isat]

            # Finally, isolate points that are saturated
            Isat[I] = np.logical_and( s[I]<=ssV, s[I]>=ssL )
            T[Isat] = Tsat[Isat]
            if quality:
                x[Isat] = (s[Isat] - ssL)/(ssV-ssL)

        # Isat is now a down-select array
        Isat = np.logical_not(Isat)
        self._iter1(
                self._tpiter1,
                'T',
                s,
                T,
                Isat,
                Ta, Tb,
                param={'fn':self._s, 'p':p},
                verbose=debug)
                
        if quality:
            return T,x
        return T


    def T_h(self, h, p=None, quality=False, debug=False):
        """Temperature from enthalpy
    T = T_h(h, p=None, quality=False)

Pressure is optional, but enthalpy is mandatory.

The optional keyword flag, quality, will cause quality to be returned
along with temperature.

    T,x = T_h(s, p=None, quality=True)
"""
        # Prepare the h array
        h = pm.units.energy(
                np.asarray(h,dtype=float),
                to_units='J')
        pm.units.matter(h, self.data['mw'],
                to_units='kg', exponent=-1, inplace=True)
        if h.ndim == 0:
            h = np.reshape(h, (1,))

        # Set a default pressure?
        if p is None:
            p = pm.config['def_p']

        p = pm.units.pressure(
                np.asarray(p, dtype=float),
                to_units='Pa')
        # Enforce pressure limits
        if ((p<self.data['plim'][0]).any() or
                (p>self.data['plim'][1]).any()):
            raise pm.utility.PMParamError(
                'MP1: Pressure is out-of-bounds.')
        
        if p.ndim == 0:
            p = np.reshape(p, (1,))

        # broadcast
        h,p = np.broadcast_arrays(h,p)
        # Initialize results
        T = np.zeros_like(h, dtype=float)
        x = -np.ones_like(h,dtype=float)
        # Some important intermediates
        Isat = np.zeros_like(h, dtype=bool)
        Ta = np.zeros_like(h, dtype=float)
        Tb = np.zeros_like(h, dtype=float)
        
        
        # Start with super-critical points
        I = p >= self.data['pc']
        Ta[I] = self.data['Tlim'][0]
        Tb[I] = self.data['Tlim'][1]
        T[I] = 0.5*(Ta[I] + Tb[I])
        # Now work with the sub-critical points
        I = np.logical_not(I)
        if I.any():
            # Get the saturation temperatures
            Tsat = self._Ts(p[I])
            # And densities
            dsL = self._dsl(Tsat,0)[0]
            dsV = self._dsv(Tsat,0)[0]
            # finally, get the saturation entropies
            hsL = self._h(Tsat,dsL,0)[0]
            hsV = self._h(Tsat,dsV,0)[0]

            # Isolate points that are liquid
            Isat[I] = h[I] < hsL
            Ta[Isat] = self.data['Tlim'][0]
            Tb[Isat] = Tsat[Isat]
            T[Isat] = Tb[Isat] - Ta[Isat]
            # Grow the boundary by 5%
            Tb[Isat] += 0.05*T[Isat]
            T[Isat] = Ta[Isat] + 0.5*T[Isat]
            
            # Isolate points that are vapor
            Isat[I] = h[I] > hsV
            Ta[Isat] = Tsat[Isat]
            Tb[Isat] = self.data['Tlim'][1]
            T[Isat] = Tb[Isat] - Ta[Isat]
            # Grow the boundary by 5%
            Ta[Isat] -= 0.05*T[Isat]
            T[Isat] = Tb[Isat] - 0.5*T[Isat]

            # Finally, isolate points that are saturated
            Isat[I] = np.logical_and( h[I]<=hsV, h[I]>=hsL )
            T[Isat] = Tsat[Isat]
            if quality:
                x[Isat] = (h[Isat] - hsL)/(hsV-hsL)
                
        # Isat is now a down-select array
        Isat = np.logical_not(Isat)
        self._iter1(
                self._tpiter1,
                'T',
                h,
                T,
                Isat,
                Ta, Tb,
                param={'fn':self._h, 'p':p},
                verbose=debug)
                
        if quality:
            return T,x
        return T
