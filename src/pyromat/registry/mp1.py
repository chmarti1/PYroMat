# MP1
#   PYroMat Multi-phase generalist class
#   Calculates physical properties from a fit for the helmholtz free 
#   energy in terms of density and temperature.
#
#   VERSION 2 (7/2020)
#   - fixed bug in _tpiter(): bad partial derivative conversion
#   - added singularity detection to _ar() for improved performance
#       near the critical point
#   

import numpy as np
import pyromat as pm
import os,sys


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

Most property methods also accept the "quality" as an optional keyword.  When
it is set to True, the property will also return the vapor/liquid mixture
quality in a tuple with the property value.  For example,
    h,x = h(T,d,quality=True)

There are also saturation property methods:
    es()    Saturation internal energy
    hs()    Saturation enthalpy
    ss()    Saturation entropy
    
And saturation equations of state methods:
    Ts()    Saturaiton temperature
    ds()    Saturation density
    ps()    Saturation pressure

Saturation methods accept either temperature or pressure as an argument.  
The density saturation method returns both liquid and vapor densities in a 
tuple pair.  See their in-line documentaiton for more details.

What is the chemical formula and charge of the substance?
    atoms() Returns a dictionary specifying the chemical composition.

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
    
Inner routines begin with a "_" to emphasize that they are not part of
the standard interface.  Most property functions are wrappers for inner
routine property functions, so they may call each other when needed.  
Inner routine property functions (like _h, _s, _p, etc...) have standard
call signatures that require temperature and density, and return the 
property and its derivatives to temperature and density.
    
VERY rarely, these routines might be called by the user.  They are 
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
    coef        a coefficient group to be passed to _satfit()
    fn          An integer index identifying the fit form to use 
                (see the _satfit method for details)

DSLgroup        Saturated liquid density data group
    Tscale      Temperature scale for normalizing T in the fit
    dscale      Density scale for re-scaling the result
    coef        a coefficient group to be passed to _poly1()
    fn          An integer index identifying the fit form to use 
                (see the _satfit method for details)

DSVgroup        Saturated vapor density data group; a dict containing:
    Tscale      Temperature scale for normalizing T in the fit
    dscale      Density scale for re-scaling the result
    coef        a coefficient group to be passed to _poly1()
    fn          An integer index identifying the fit form to use 
                (see the _satfit method for details)

AOgroup        Helmholtz free energy ideal gas group; a dict containing:
    Tscale      Temperature scale for normalizing T
    dscale      density scale for normalizing d
    logt        a scalar coefficient of a log(tt) term
    coef0       a coefficient list to be passed to _poly1() to build p0 below
    coef1       a simple Nx2 coefficient list used to build q(tt) below
If tt = Tscale/T    <=== INVERSE!
and dd = d/dscale
    ao = log(d) + LOGT*log(tt) + p0(tt) + q(tt)
        q1(tt) = sum_k coef1[k,1] * log(1 - exp(-tt*coef[k,0]))
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
    [ t, d, b, a, gam, ep, c ], ...
]
    ar1 = c * dd**d * tt**t * exp(-a*(dd-ep)**2 - b*(tt-gam)**2) + ...
    Ar1 = ar1 * R * T
If coef1 is defined it will be combined with the other coefficients
to form the residual.  If coef1 is not defined, it will be ignored.

coef2 is an optional list of lists of coefficients forming a matrix
[...
    [ a, b, m, A, B, C, D, c ], ...
]
    X = ((1-tt) + A*((dd-1)**2)**(0.5/m))**2 + B*((dd-1)**2)**a
    ar2 = c * X**b * d * exp(-C*(dd-1)**2 - D*(tt-1)**2) + ...
    Ar2 = ar2 * R * T

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
atoms           A dictionary with a key for each atom and a value for 
                its count in the molecule.  For example, CO2 would 
                have content = {'C':1, 'O':2}
                
There are also the typical mandatory PYroMat meta data elements:
id              What substance is this?
doc             Where did it come from?
class           What class should be used to evaluate the data?
"""
    
    
    def _test(self, tab, sattab, report=None, basic=False):
        """Test the MP1 class model
    _test(tab, sattab)     # Prints to stdout
        OR
    _test(tab, sattab, report_file='/path/to/report')  # Prints to a file
        OR
    _test(tab, sattab, report_file=open_file_descriptor)   # Prints to a file
    
tab and sattab are nested lists or 2D numpy arrays forming tables of "truth"
data used for validation of the core data and property methods.

If the optional "basic" keyword is set to True, only the data integrity checks
are completed (see below).

** tab **
The TAB table is used to test the core properties, and should have columns 
and units

    T (K)   p (Pa)  d (kg/m3)   cp (kJ/kg/K)    s (kJ/kg/K)     h (kJ/kg)

The units are selected to match those typically used in the publication of
so-called Span and Wagner equations of state.

** sattab **
The SATTAB table is used to test the saturation property functions, and should
have columns and units

    T (K)   p (Pa)  dL (kg/m3)  dV (kg/m3)

where dL and dV are the liquid and vapor densities respectively.

Test criteria:
(0) Data Integrity
    0.1 AOgroup must contain positive scalars, Tscale, dscale
    0.2 AOgroup's logt parameter must be a scalar
    0.2 AOgroup's coef0 parameter must be a legal poly1 group
    0.4 AOgroup's coef1 parameter must be a table with two columns
    0.5 ARgroup must contain Tscale, dscale
    0.6 ARgroup's coef0 must be a list of legal poly2 groups
    0.7 ARgroup's coef1 must be a table with seven columns
    0.8 ARgroup's coef2 must be a table with eight columns
    0.9 PSgroup must contain positive scalars Tscale, pscale and integer, fn
    0.10 PSgroup's coef must be a legal poly1 group
    0.11 DSLgroup must contain positive scalars Tscale, dscale and integer, fn
    0.12 DSLgroup's coef must be a legal poly1 group
    0.13 DSVgroup must contain positive scalars Tscale, dscale and integer, fn
    0.14 DSVgroup's coef must be a legal poly1 group
(1) Saturation 
    1.1 Saturation pressure agrees to within 0.01%
    1.2 Liquid saturation density agrees to within 0.01%
    1.3 Vapor saturation density agrees to within 0.01%
(2) Inverse Saturation
    2.1 Saturation temperature agrees to within 0.01%
(3) Equation of State
    3.1 Density must agree to within 0.01%
    3.2 Pressure must agree to within 0.01%
    3.3 Temperature must agree to within 0.01%
(4) Core Properties
    4.1 Specific heat (cp) must agree to within 0.01%
    4.2 Entropy must agree to within 0.01%
    4.3 Enthalpy must agree to within 0.01%
(5) Inverse Properties
    5.1 Temperature from entropy and pressure must agree to within 0.01%
    5.2 Temperature from entropy and density must agree to within 0.01%
    5.3 Density from entropy and temperature must agree to within 0.01%
    5.4 Temperature from enthalpy and pressure must agree to within 0.01%
    5.5 Temperature from enthalpy and density must agree to within 0.01%

"""
        # Recurse with a fresh file descriptor if the file is a string
        if isinstance(report, str):
            with open(report, 'w') as ff:
                return self._test(report=ff)
        elif report is None:
            report = sys.stdout
            
        result = True
        
        report.write('PYroMat version: ' + pm.config['version'] + '\n')
        report.write('Species: ' + self.data['id'] + '\n')
        
        # CRITERION 0
        # Basic Data integrity
        def _check_poly1(coef):
            # Loop through the coefficient groups
            for cgi,cc in enumerate(coef):
                if not isinstance(cc, (list,tuple)):
                    return True, f'Coef. group {cgi} was neither a list nor a tuple.'
                # Check pre- and post-exponents
                pre = cc[0]
                post = cc[1]
                if not isinstance(pre, (int,float)):
                    return True, f'In coef. group {cgi} pre-exponent must be scalar. Found: {pre}'
                if not isinstance(post, (int,float)):
                    return True, f'In coef. group {cgi} post-exponent must be scalar. Found: {post}'
                for ti, term in enumerate(cc[2:]):
                    if not isinstance(term,(list,tuple)):
                        return True, f'In coef. group {cgi}, the {ti} term is neither a list nor a tuple: {term}'
                    if len(term) != 2:
                        return True, f'In coef. group {cgi}, the {ti} term should have 2 elements.  Found: {term}'
                    if not isinstance(term[0], int) or term[0] <0:
                        return True, f'In coef. group {cgi}, the {ti} term exponent was not a non-negative integer.  {term}'

            return False, ''
            
        def _check_poly2(coef):
            # Loop through the coefficient groups
            for cgi,cc in enumerate(coef):
                if not isinstance(cc, (list,tuple)):
                    return True, f'Coef. group {cgi} was neither a list nor a tuple.'
                # Check pre- and post-exponents
                pre = cc[0]
                post = cc[1]
                if not isinstance(pre, list) or len(pre) != 2:
                    return True, f'In coef. group {cgi} pre-exponent must be a two-element list. Found: {pre}'
                if not isinstance(post, list) or len(post) != 2:
                    return True, f'In coef. group {cgi} post-exponent must be a two-element list. Found: {post}'
                prex,prey = pre
                postx,posty = post
                if not isinstance(prex, (int,float)) or not isinstance(prey, (int,float)) or not isinstance(postx, (int,float)) or not isinstance(posty, (int,float)):
                    return True, f'In coef. group {cgi}, pre- and post-exponents must be scalars.  Pre: {pre}, Post: {post}'
                
                for ti, term in enumerate(cc[2:]):
                    if not isinstance(term,(list,tuple)):
                        return True, f'In coef. group {cgi}, the {ti} term is neither a list nor a tuple: {term}'
                    if len(term) != 3:
                        return True, f'In coef. group {cgi}, the {ti} term should have 3 elements.  Found: {term}'
                    if not isinstance(term[0], int) or term[0] <0:
                        return True, f'In coef. group {cgi}, the {ti} term exponent was not a non-negative integer.  {term}'

            return False, ''

        #0.1 AOgroup must contain positive scalars, Tscale, dscale
        error = False
        for test in ['Tscale', 'dscale']:
            if test in self.data['AOgroup']:
                pass
            elif self.data['AOgroup'][test] > 0:
                pass
            else:
                error = True
                break
        if error:
            report.write('[FAILED]    0.1 AOgroup must contain positive scalars, Tscale, dscale\n')
            report.write('            ' + test + ' = ' + repr(self.data['AOgroup'][test]) + '\n')
            result = False
        else:
            report.write('[passed]    0.1 AOgroup must contain positive scalars, Tscale, dscale\n')
            
        #0.2 AOgroup's logt parameter must be a scalar
        if 'logt' in self.data['AOgroup'] and not isinstance(self.data['AOgroup']['logt'], (float,int)):
            report.write('[FAILED]    0.2 AOgroup logt parameter must be a scalar\n')
            report.write('            logt = ' + repr(self.data['AOgroup']['logt']) + '\n')
            result = False
        else:
            report.write('[passed]    0.2 AOgroup logt parameter must be a scalar\n')
        
        #0.3 AOgroup's coef0 parameter must be a legal poly1 group
        error = False
        if 'coef0' in self.data['AOgroup']:
            error,message = _check_poly1(self.data['AOgroup']['coef0'])
        if error:
            report.write('[FAILED]    0.3 AOgroup coef0 parameter must be a legal poly1 group\n')
            report.write('            ' + message + '\n')
            result = False
        else:
            report.write('[passed]    0.3 AOgroup coef0 parameter must be a legal poly1 group\n')
        
        #0.4 AOgroup's coef1 parameter must be a table with two columns
        error = False
        if 'coef1' in self.data['AOgroup']:
            if not isinstance(self.data['AOgroup']['coef1'], (list,tuple)):
                error = True
                message = 'coef1 was not iterable.'
            else:
                for row in self.data['AOgroup']['coef1']:
                    if not isinstance(row, (list,tuple)) or len(row) != 2:
                        error = True
                        message = 'coef1 has at least one row without 2 columns'
        if error:
            report.write('[FAILED]    0.4 AOgroup coef1 parameter must be a table with two columns\n')
            report.write('            ' + message + '\n')
            result = False
        else:
            report.write('[passed]    0.4 AOgroup coef1 parameter must be a table with two columns\n')            
            
        #0.5 ARgroup must contain Tscale, dscale
        error = False
        for test in ['Tscale', 'dscale']:
            if test in self.data['ARgroup']:
                pass
            elif self.data['ARgroup'][test] > 0:
                pass
            else:
                error = True
                break
        if error:
            report.write('[FAILED]    0.5 ARgroup must contain positive scalars, Tscale, dscale\n')
            report.write('            ' + test + ' = ' + repr(self.data['ARgroup'][test]) + '\n')
            result = False
        else:
            report.write('[passed]    0.5 ARgroup must contain positive scalars, Tscale, dscale\n')
            
        #0.6 ARgroup's coef0 must be a list of legal poly2 groups
        error = False
        message = ''
        if 'coef0' in self.data['ARgroup']:
            if not isinstance(self.data['ARgroup']['coef0'], (list,tuple)):
                error = True
                message = 'coef0 was not iterable.'
            else:
                for cgi,cg in enumerate(self.data['ARgroup']['coef0']):
                    error,message = _check_poly2(cg)
                    if error:
                        message = 'Error in term number ' + str(cgi) + '\n            ' + message
                        break
        if error:
            result = False
            report.write('[FAILED]    ARgroup coef0 must be a list of legal poly2 groups\n')
            report.write('            ' + message + '\n')
        else:
            report.write('[passed]    ARgroup coef0 must be a list of legal poly2 groups\n')
            
        #0.7 ARgroup's coef1 must be a table with seven columns
        error = False
        message = ''
        if 'coef1' in self.data['ARgroup']:
            if not isinstance(self.data['ARgroup']['coef1'], (list,tuple)):
                error = True
                message = 'coef1 was not iterable.'
            else:
                for row in self.data['ARgroup']['coef1']:
                    if len(row) != 7:
                        error = True
                        message = 'coef1 has at least one row without 7 columns'
        if error:
            result = False
            report.write('[FAILED]    ARgroup coef1 must be a table with seven columns\n')
            report.write('            ' + message + '\n')
        else:
            report.write('[passed]    ARgroup coef1 must be a table with seven columns\n')
            
        #0.8 ARgroup's coef2 must be a table with eight columns
        error = False
        message = ''
        if 'coef2' in self.data['ARgroup']:
            if not isinstance(self.data['ARgroup']['coef2'], (list,tuple)):
                error = True
                message = 'coef2 was not iterable.'
            else:
                for row in self.data['ARgroup']['coef2']:
                    if len(row) != 8:
                        error = True
                        message = 'coef2 has at least one row without 8 columns'
        if error:
            result = False
            report.write('[FAILED]    ARgroup coef2 must be a table with eight columns\n')
            report.write('            ' + message + '\n')
        else:
            report.write('[passed]    ARgroup coef2 must be a table with eight columns\n')
        
        #0.9 PSgroup must contain positive scalars Tscale, pscale and integer, fn
        #0.10 PSgroup's coef must be a legal poly1 group
        #0.11 DSLgroup must contain positive scalars Tscale, dscale and integer, fn
        #0.12 DSLgroup's coef must be a legal poly1 group
        #0.13 DSVgroup must contain positive scalars Tscale, dscale and integer, fn
        #0.14 DSVgroup's coef must be a legal poly1 group
        
        if basic:
            return result
        
        # First, extract the sattab columns
        sattab = np.asarray(sattab, dtype=float)
        T = pm.units.temperature_scale(sattab[:,0], from_units='K')
        p = pm.units.pressure(sattab[:,1], from_units='Pa')
        dL = pm.units.matter(self.data['mw'], sattab[:,2], from_units='kg')
        pm.units.volume(dL, from_units='m3', inplace=True)
        dV = pm.units.matter(self.data['mw'], sattab[:,3], from_units='kg')
        pm.units.volume(dV, from_units='m3', inplace=True)
        
        # CRITERION 1
        # Saturation numerical integrity
        report.write('\n1. Saturation Properties\n')
        result = pm.utility.proptest(self.ps, {'T':T}, truth=p, ep=.0001, 
                text='1.1 Saturation pressure must agree to within .01%', 
                report=report) or result
        result = pm.utility.proptest(self.ds, {'T':T}, truth=dL, ep=.0001, 
                text='1.2 Liquid saturation density must agree to within .01%', 
                report=report, findex=0) or result
        result = pm.utility.proptest(self.ds, {'T':T}, truth=dV, ep=.0001, 
                text='1.3 Vapor saturation density must agree to within .01%', 
                report=report, findex=1) or result
        # CRITERION 2
        # Inverse saturation
        report.write('2. Inverse Saturation Properties\n')
        result = pm.utility.proptest(self.Ts, {'p':p}, truth=T, ep=.0001, 
                text='2.1 Saturation temperature must agree to within .01%', 
                report=report, findex=1) or result
        # Throw away the saturation table values; we're done with those
        # Switch to the core property table
        T = pm.units.temperature_scale(tab[:,0], from_units='K')
        p = pm.units.pressure(tab[:,1], from_units='Pa')
        d = pm.units.matter(self.data['mw'], tab[:,2], from_units='kg')
        pm.units.volume(d, from_units='m3', exponent=-1, inplace=True)
        cp = pm.units.energy(tab[:,3], from_units='kJ')
        pm.units.temperature(cp, from_units='K', exponent=-1, inplace=True)
        pm.units.matter(self.data['mw'], cp, from_units='kg', exponent=-1, inplace=True)
        s = pm.units.energy(tab[:,4], from_units='kJ')
        pm.units.temperature(s, from_units='K', exponent=-1, inplace=True)
        pm.units.matter(self.data['mw'], s, from_units='kg', exponent=-1, inplace=True)
        h = pm.units.energy(tab[:,5], from_units='kJ')
        pm.units.matter(self.data['mw'], h, from_units='kg', exponent=-1, inplace=True)
        
        # CRITERION 3
        # Equation of state
        
        return result
        

    def _poly2(self,x,y,group,diff=2):    
        """Polynomial evaluation (primative routine)
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

The behavior of poly2 is very much the same as poly1, but for functions
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
        """Polynomial evaluation (primative routine)
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
                ep=1e-6, Nmax=20, fx_index=1, 
                verbose=False, param={}):
        """Modified Newton iteration on a 1D inner routine. (primative routine)
        
    _iter1(fn, prop, y, x, Ids, xmin, xmax,)

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
            array.
Ids         A down-select boolean index array; only x[Ids],y[Ids] will 
            be evaluated.  This allows iteration in-place on data sets 
            where only a portion of the data require iteration.  If y is
            a floating point array with N elements, Ids must be a bool
            array with N elements.  It will specify a down-selected 
            data set with M elements, where M<=N.
xmin, xmax  Upper and lower limit arrays for the x values.  These must
            broadcastable to match x and y.  Even values outside of the
            down-select region should have legal values.  Note that 
            these arrays are volatile, and will be written to by the
            bisection process.
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
            print('x', 'yvalue', 'dydx', 'dx', 'Ids')


        # Create an argument dictionary
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
            while IooB[Ids].any():
                dx[IooB] /= 2.
                x[IooB] -= dx[IooB]
                IooB[Ids] = np.logical_or( x[Ids] < xmin[Ids], x[Ids] > xmax[Ids])
                # Prevent a while-loop-trap
                count_oob += 1
                if count_oob>Nmax:
                    raise pm.utility.PMAnalysisError(
                        '_iter1() failed to produce a guess that was in-bounds')
            
            # Check the iteration convergence
            Ids[Ids] = abs(error[Ids]) > abs(ep*y[Ids])
            # Prevent a while-loop-trap
            count += 1
            if count>Nmax:                
                pm.utility.print_warning(\
                    '_iter1() failed to converge for %d elements after %d attempts'%(\
                    Ids.sum(), Nmax))
                return

    
    def _hybrid1(self, fn, prop, y, x, Ids, xmin, xmax,
                ep=1e-6, Nmax=20, fx_index=1, 
                verbose=False, paranoid=True, param={}):
        """Hybrid numerical inversion of an inner routine (primative routine)
        
    _hybrid1(fn, prop, y, x, Ids, xmin, xmax,)

This hybrid iteration algorithm is named for being a hybrid of biseciton
and Newton iteration.  On "well behaved" functions it converges as 
quickly as the Newton algorithm, but on "badly behaved" functions, it 
is extremely stable.

Iteration is performed in-place on the x array until fn(x) == y.  The 
hybrid1 algorithm depends on the xmax and xmin values to bracket a 
solution.  The funciton, fn, and its derivative are evaluated at the 
maximum and minimum, and Newton's method is used to generate two 
candidate next guesses.  The point bisecting the maximum and minimum is
calculated, providing a third candidate guess.  Of the three candidates,
the one in the middle is selected for the next iteration step.  If the
middle point lies outside of xmax and xmin, then the bisection point is
selected instead.

Once a next guess is selected, the function and its derivative are 
evaluated there.  This guess is used to replace either xmin or xmax, 
just like would be done in a bisection algorithm, but then the three-
candidate voting algorithm is repeated.  Since the calculated next guess
of the boundaries is unchanged from the last iteration step, only one
function evaluation is required per step, making the computational cost
comparable with Newton's method.

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
            array.
Ids         A down-select boolean index array; only x[Ids],y[Ids] will 
            be evaluated.  This allows iteration in-place on data sets 
            where only a portion of the data require iteration.  If y is
            a floating point array with N elements, Ids must be a bool
            array with N elements.  It will specify a down-selected 
            data set with M elements, where M<=N.
xmin, xmax  Upper and lower limit arrays for the x values.  These must
            broadcastable to match x and y.  Even values outside of the
            down-select region should have legal values.  Note that 
            these arrays are volatile, and will be written to by the
            bisection process.
*** Optional Parameters ***
ep          Epsilon; fractional error permitted in y (default 1e-6)
Nmax        Maximum number of iterations (default 20)
fx_index    The location of the property derivative in the call 
            signature (default 1)
paranoid    If any of the candidate guesses is out of bounds, the revert to
            bisection.  Otherwise, only test the median guess.  Paranoid 
            operation can be essential in functions with +/- slope inflections
            in the domain.
param       A dicitonary of keyword arguments are passed directly to the 
            inner routine being inverted.

"""
        #================================#
        # Initialize intermediate arrays #
        #================================#
        # Produce arrays of candidate guesses xa and xb are produced by 
        # the Newton algorithm from xmin and xmax respectively.  
        # xc is produced by bisection.
        xa = np.zeros_like(x, dtype=float)
        xb = np.zeros_like(x, dtype=float)
        xc = np.zeros_like(x, dtype=float)
        
        # Make local copies of xmax and xmin
        xmax = np.array(xmax)
        xmin = np.array(xmin)
        
        # The Iab, Ibc, and Ica indices are used to store comparsion
        # truth values for sorting the candidate solutions, and Iwork
        # is used to assign the values
        Iab = np.zeros_like(Ids, dtype=bool)
        Ibc = np.zeros_like(Ids, dtype=bool)
        Ica = np.zeros_like(Ids, dtype=bool)
        Iwork = np.zeros_like(Ids, dtype=bool)
        Iaoob = np.zeros_like(Ids, dtype=bool) # Out-of-bounds arrays
        Iboob = np.zeros_like(Ids, dtype=bool) # 
        Iswap = np.zeros_like(Ids, dtype=bool) # which were swapped?
        
        if verbose:
            print("Fn: " + repr(fn.__name__))
            print("param: " + repr(param))
        
        # Initialize an argument dicitonary
        arg = param.copy()
        
        # Build the argument list
        for k,v in param.items():
            # For any array arguments, shrink them along with Ids
            if isinstance(v,np.ndarray):
                arg[k] = v[Ids]
        # Now, we'll evalaute the funciton at the limits
        # Start at the minimum
        arg[prop] = xmin[Ids]
        FF = fn(diff=1, **arg)
        yy = FF[0]
        yyx = FF[fx_index]
        # Calculate the first candidate solution
        xa[Ids] = xmin[Ids] + (y[Ids] - yy)/yyx
        # If f(xmin) > f(xmax) then the nominal slope of the curve is negative
        # That means that these boundaries will need to be updated in reverse
        # of the other boundaries.
        Iswap[Ids] = yy >= y[Ids]
        
        # Now, evaluate at the maximum 
        arg[prop] = xmax[Ids]
        FF = fn(diff=1, **arg)
        yy = FF[0]
        yyx = FF[fx_index]
        # Calculate the second candidate solution
        xb[Ids] = xmax[Ids] + (y[Ids] - yy)/yyx
        
        # Verify that the limits bracket a solution
        if not np.logical_xor(Iswap[Ids], yy >= y[Ids]).all():
            if verbose:
                print(" xmin  xmax  xa  xb  xc ")
                print(xmin, xmax, xa, xb, xc)
                print(" y")
                print(y)
            pm.utility.print_error('_HYBRID1: Failure to bracket a solution. Check function arguments to be sure they reference a valid state. This error usually occurs in inversion routines very close to the saturation line or if the properties are out-of-range.')
            raise pm.utility.PMParamError('_HYBRID1: At least one max/min value does not bracket a solution!')
        
        # Calculate the thrid candidate solution
        xc[Ids] = 0.5*(xmin[Ids] + xmax[Ids])
        
        if verbose:
            print(" xmin  xmax  xa  xb  xc ")
        
        count = 0
        while Ids.any():
            if count>Nmax:
                pm.utility.print_warning(f'_HYBRID1: Failed to converge for {Ids.sum()} elements in {Nmax} iterations.')
                return
            
            if verbose:
                print(xmin, xmax, xa, xb, xc)
            
            # Clean the worker indexes
            Iab[:] = False
            Ibc[:] = False
            Ica[:] = False
            Iwork[:] = False
            
            # The last step has established three candidate solutions
            # Which should we select?  First, compare the three candidate
            # solutions to determine which is in the middle
            Iab[Ids] = xa[Ids] < xb[Ids]
            Ibc[Ids] = xb[Ids] < xc[Ids]
            Ica[Ids] = xc[Ids] < xa[Ids]
            
            # Now, assign all values for which xa is the next guess
            Iwork[Ids] = Iab[Ids] == Ica[Ids]
            x[Iwork] = xa[Iwork]
            # Now, assign all values for which xb is the next guess
            Iwork[Ids] = Iab[Ids] == Ibc[Ids]
            x[Iwork] = xb[Iwork]
            # Now, assign all value for which xc is the next guess
            Iwork[Ids] = Ibc[Ids] == Ica[Ids]
            x[Iwork] = xc[Iwork]
            # finally, deal with the xa and xb out-of-bounds case
            if paranoid:
                # In paranoid mode, either xa and xb being out of bounds
                # forces xc to be selected
                Iwork[Ids] = np.logical_or(xa[Ids] < xmin[Ids], xa[Ids] > xmax[Ids])
                x[Iwork] = xc[Iwork]
                Iwork[Ids] = np.logical_or(xb[Ids] < xmin[Ids], xb[Ids] > xmax[Ids])
                x[Iwork] = xc[Iwork]
            else:
                Iwork[Ids] = np.logical_or(x[Ids] < xmin[Ids], xa[Ids] > xmax[Ids])
                x[Iwork] = xc[Iwork]
                
                        
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
            
            # use xc as a temporary variable
            # First, calculate the size of the change in x
            xc[Ids] = x[Ids] + (y[Ids] - yy) / yyx
            # x is now the next guess
            # xc is its next projected guess
            
            # Where is the guess?
            # Should it be stored in xmin?
            # xor with Iswap forces a swap when needed
            Iwork[Ids] = np.logical_xor(yy <= y[Ids], Iswap[Ids])
            xmin[Iwork] = x[Iwork]
            xa[Iwork] = xc[Iwork]
            # Test OOB
            Iaoob[Iwork] = np.logical_or(xc[Iwork] <= xmin[Iwork], xc[Iwork] >= xmax[Iwork])
            # or in xmax
            Iwork[Ids] = np.logical_not(Iwork[Ids])
            xmax[Iwork] = x[Iwork]
            xb[Iwork] = xc[Iwork]
            Iboob[Iwork] = np.logical_or(xc[Iwork] <= xmin[Iwork], xc[Iwork] >= xmax[Iwork])            
            # Calculate the new bisection point
            xc[Ids] = 0.5*(xmax[Ids] + xmin[Ids])
            
            # Check for convergence
            # if xmax-xmin is small OR
            # if the y error is small
            Ids[Ids] = np.logical_and(\
                    np.abs((xmax[Ids] - xmin[Ids])/x[Ids]) > ep,\
                    np.abs((yy - y[Ids])/y[Ids]) > ep)
                        
            # Prevent a while-loop-trap
            count += 1

        if verbose:
            print(f"Converged for all elements in {count} iterations.")


    def _tditer(self, T, d, fn, diff=1, debug=False):
        """T,d iterator wrapper (primative routine)
    _tditer(T,d,fn, diff=1, debug=False)
    
    This wrapper function evaluates a property inner routine from
temperature and density.  It is intended to be used to allow ID 
iteration on a property with respect to temperature with constant 
density (isochoric).  

_tditer accpets three required arguments:
    T   Temperature numpy array in Kelvin
    d   Density numpy array in kg/m3
    fn  Property inner routine to be evaluated (e.g. self._h or self._s)
    
For example, a call to _hybrid1 to calculate temperature while 
specifying entropy and pressure might appear

    self._hybrid1( self._tditer, # Don't use _s, use _tditer
        'T',                    # We want to calculate T
        svalues,                # Here are the target entropy values
        T,                      # The T array
        Ids,                    # The pre-initialized down-select array
        Tmin, Tmax,             # T bounds
        param={'fn':self._s, 'd':dvalues})

Optional parameters (and their defaults) are:
    debug (False)   
    Has no effect.  It is included only to provide the same call signature as
    _tpiter().
    
    diff (1)
    When 0, returns no derivatives.  When 1, returns the first-order 
    derivatives.
"""
        y = np.empty_like(T, dtype=float)
        yt = np.empty_like(T, dtype=float)
        yd = np.empty_like(T, dtype=float)
        dsL = np.empty_like(T, dtype=float)
        dsLt = np.empty_like(T, dtype=float)
        dsV = np.empty_like(T, dtype=float)
        dsVt = np.empty_like(T, dtype=float)

        # Only check sub-critical temperatures for saturation
        I = T < self.data['Tc']
        dsL[I],dsLt[I],_ = self._dsl(T[I],diff=diff)
        dsV[I],dsVt[I],_ = self._dsv(T[I],diff=diff)
        # Down-select for densities between the saturation properties
        I[I] = np.logical_and(d[I] <= dsL[I], d[I] >= dsV[I])
        if I.any():
            # Use xd as a temporary result
            # Start with the denominator
            xd = 1./(dsL[I]/dsV[I] - 1.)
            # Calculate quality
            x = (dsL[I]/d[I] - 1.) * xd
            if diff:
                # Construct quality's deriv. w.r.t. temperature in steps
                # First, the most complicated term: the denominator's derivative
                xt = x * xd * (dsL[I]/dsV[I]) * (dsVt[I]/dsV[I] - dsLt[I]/dsL[I])
                # Continue to construct xt
                xt += xd*dsLt[I]/d[I]
                # Finalize xd
                xd *= -dsL[I] / d[I] / d[I]
            # Evaluate the saturation properties
            y[I],yt[I],yd[I] = fn(T[I],dsV[I],diff=diff)
            yy,yyt,yyd = fn(T[I],dsL[I],diff=diff)
            # Calculate the mixture properties and derivatives
            if diff:
                yt[I] = x*yt[I] + xt*y[I] + (1-x)*yyt - xt*yy
                yd[I] = x*yd[I] + xd*y[I] + (1-x)*yyd - xd*yy
            y[I] = x*y[I] + (1-x)*yy
        # Now deal with points that are not under the dome.
        I = np.logical_not(I)
        if I.any():
            y[I],yt[I],yd[I] = fn(T[I],d[I],diff=diff)
            
        return y,yt,yd
        

    def _tpiter(self, T, p, fn, diff=1, debug=False):
        """T,p iterator wrapper (primative routine)
    _tpiter(T,p,fn, diff=1, debug=False)
    
    This wrapper function evaluates a property inner routine from 
temperature and pressure.  It is intended to be used to allow 1D 
iteration on a property with respect to temperature with constant 
pressure (isobaric).  Because the property functions require density
and temperature, _tpiter implements an inner hybrid1 iteration routine
to determine the density for each value of temperature.

When diff is 1 (as it needs to be for _hybrid1 to work properly), the 
property's partial derivatives need to be shifted from constant-density
into constant pressure space.

_tpiter accpets three required arguments:
    T   Temperature numpy array in Kelvin
    p   Pressure numpy array in Pa
    fn  Property inner routine to be evaluated (e.g. self._h or self._s)
    
For example, a call to _hybrid1 to calculate temperature while 
specifying entropy and pressure might appear

    self._hybrid1( self._tpiter, # Don't use _s, use _tpiter
        'T',                    # We want to calculate T
        svalues,                # Here are the target entropy values
        T,                      # The T array
        Ids,                    # The pre-initialized down-select array
        Tmin, Tmax,             # T bounds
        param={'fn':self._s, 'p':pvalues})
        
Optional parameters (and their defaults) are:
    debug (False)   
    Passed to the verbose parameter of the inner _hybrid1 routine.  This 
    produces quite a bit of text in most iterations.  It should not be used 
    by most users.
    
    diff (1)
    When 0, returns no derivatives.  When 1, returns the first-order 
    derivatives.
"""
        
        d = np.empty_like(p, dtype=float)
        I = np.zeros_like(p, dtype=bool)
        da = np.empty_like(p, dtype=float)
        db = np.empty_like(p, dtype=float)
        
        # We need to establish boundaries for the isothermal density iterations
        # Start with super-critical points
        Isat = T >= self.data['Tc']
        # use 10% of the IG density for the lower boundary
        da[Isat] = 0.1 * p[Isat] / (self.data['R'] * T[Isat])
        db[Isat] = self.data['dlim'][1]
        
        # Now, move on to the sub-critical points
        Isat = np.logical_not(Isat)
        if Isat.any():
            # Use da and db as temporaries for the saturation densities
            ps = self._ps(T[Isat])[0]
            db[Isat] = self._dsv(T[Isat])[0]    # vapor 
            da[Isat] = self._dsl(T[Isat])[0]    # liquid
            
            # Start with the vapor points
            # Use densities slightly inside the dome
            # It is important to use specific volumes at low temperatures
            # The vapor densities are so low, that even a small liquid percent
            # dominates the value.
            I[Isat] = p[Isat] < ps
            db[I] = 1./(0.98 / db[I] + 0.02 / da[I])
            da[I] = 0.1 * p[I] / (self.data['R'] * T[I])
            
            # Now, the liquid points
            I[Isat] = np.logical_not(I[Isat])
            da[I] = 0.98 * da[I] + 0.02 * db[I]
            db[I] = self.data['dlim'][1]
        
        # Re-initialize the down-select vector to True
        I[:] = True
        
        self._hybrid1(self._p, 'd', p, d, I, da, db,
                ep=1e-6, Nmax=30, fx_index=2, 
                verbose=debug, param={'T':T})
        
        # Assume a standard inner property routine call signature
        y,yt,yd = fn(T,d,diff=diff)
        # If the derivative is requested, we need to shift from constant
        # density to constant pressure.
        if diff>0:
            _,pt,pd = self._p(T,d,diff=1)
            # Correct the partial derivatives of the property to be 
            #  with respect to T,p instead of T,d.  Since d is used for
            #  density, let's use D for derivative and _T or _d for
            #  partial derivatives
            # The property, y,
            #   Dy(T,d) = y_T DT + y_d Dd    <== as evaluated by fn()
            # Pressure, p,
            #   Dp(T,d) = p_T DT + p_d Dd    <== as evaluated by _p()
            # So, differentials in density w.r.t. temperature while 
            # holding pressure constant, Dp = 0, and
            #   Dd/DT | p=const = -p_T / p_d
            # Differentials with density w.r.t. pressure while holding 
            # temperature constant, DT = 0, and
            #   Dd/Dp | T=const = 1/p_d
            # Therefore, 
            #   Dd = (-p_T / p_d) DT + (1 / p_d) Dp
            # Finally,
            #   Dy = (y_T - y_d p_T / p_d) DT + (y_d / p_d) Dp
            yt = yt - yd * pt / pd
            yp = yd / pd
        # Do not support higher derivatives than 1
        return y,yt,yp
        

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
                
                # So long as m > 0, this will be fine, even when ddm1==0
                # For all data examined so far, m > 0
                # p = A(dd-1)**m
                p = AA*(ddm1*ddm1)**(0.5*m)
                
                if diff>0:
                    # Detect indices where ddm1 is zero.  This will cause 
                    # singularities in differentiation
                    Izero = (np.abs(ddm1) < 1e-6)
                    Inzero = np.logical_not(Izero)

                    # Case out the temperature and density derivatives separately
                    # The density derivatives are complicated by singularity!
                    pt = 0
                    if diff>1:
                        ptt = 0
                        ptd = 0
                    pd = np.zeros_like(p)
                    pdd = np.zeros_like(p)
                    # Case out the nearly singular densities
                    # In virtually all cases, the function can be evaluated near the
                    # singularity, but it is more expensive.
                    if Izero.any():
                        # if 0.5*m < 1, there is an unreconcilable singularity!
                        if m < 2:
                            pm.utility.print_warning('_ar():: this model is very nearly singular near critical density.')
                        # borrow e for a temporary variable
                        e[Izero] = AA * m * (ddm1[Izero]*ddm1[Izero])**(0.5*m-1)
                        pd[Izero] = e[Izero] * ddm1[Izero]
                        if diff>1:
                            pdd[Izero] = e[Izero] * (m-1)
                    elif Inzero.any():
                        pd[Inzero] = m * p[Inzero] / ddm1[Inzero]
                        if diff>1:
                            pdd[Inzero] = pd[Inzero] * (m-1)/ddm1[Inzero]
                        
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
                    ed = np.zeros_like(e)
                    edd = np.zeros_like(e)
                    # This term can have a singularity from ddm1=>0 too
                    if Izero.any():
                        ed[Izero] = BB*2*a*(ddm1[Izero]*ddm1[Izero])**(a-1)*ddm1[Izero]
                        if diff>1:
                            edd[Izero] = BB*2*a*(2*a-1)*(ddm1[Izero]*ddm1[Izero])**(a-1)
                    if Inzero.any():
                        ed[Inzero] = 2*a*e[Inzero]/ddm1[Inzero]
                        if diff>1:
                            edd[Inzero] = (2*a-1)*ed[Inzero]/ddm1[Inzero]
                    pd += ed
                    if diff>1:
                        pdd += edd
                
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
        """Generic saturated property fit (primative routine)
    s, st, stt = _satfit(tt, fn=0, diff=0)

tt = T / Tc

returns saturation property normalized by its critical value

fn is an integer indicating which property fit form to use
0   poly(tt)
    coef is interpreted by poly1
1   poly(1-tt)
    coef is interpreted by poly1
2   exp(poly(1-tt))
    coef is interpreted by poly1, and the result is passed to np.exp()
3   exp(1/tt * poly(1-tt))
    coef is interpreted by poly1, the result is multiplied by 1/tt, and
    passed to np.exp()
"""
        
        if fn == 0:
            p,pt,ptt = self._poly1(tt, coef, diff=diff)
        elif fn == 1:
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
        """Saturated vapor density (inner routine)
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
        """Saturated liquid density (inner routine)
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
        """Saturated temperature from pressure (inner routine)"""
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
        
        
    def _d(self,T,p,debug=False):
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
        #da[Itest] = self.data['dlim'][0]
        # Produce a minimum density from one tenth the ideal gas relationship
        da[Itest] = 0.1 * p[Itest] / (self.data['R'] * T[Itest])
        db[Itest] = self.data['dlim'][1]
        #d[Itest] = 0.5*(self.data['dlim'][0] + self.data['dlim'][1])
        # For temperatures that are sub-critical, detect whether the 
        # state is liquid or gaseous.  Set Itest to sub-critical.  
        Itest = np.logical_not(Itest)
        if Itest.any():
            # Now, isolate the vapor points; set the upper density to the
            # saturated vapor density FORCE Istate to be an ndarray
            Istate = np.zeros_like(T, dtype=bool)
            Istate[Itest] = p[Itest] < self._ps(T[Itest], 0)[0]
            #da[Istate] = self.data['dlim'][0]
            # Produce a minimum density from half the ideal gas relationship
            da[Istate] = 0.5 * p[Istate] / (self.data['R'] * T[Istate])
            db[Istate] = self._dsv(T[Istate], 0)[0]
            #d[Istate] = db[Istate] - da[Istate]
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
        
        # perform the iteration
        #self._iter1(
        self._hybrid1(
                self._p,
                'd',
                p,
                d,
                I,
                da,
                db,
                fx_index = 2,
                param={'T':T},
                verbose=debug)
                
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
        
        self._hybrid1(
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
        # There is ONE special case when three parameters are allowed
        # T,p,x is a special case - NO D!
        elif nparam == 3 and d is not None:
            raise utility.PMPParameterError(
                    'Density may not be specified simultaneously with two other properties.')
        elif nparam > 3:
            raise utility.PMParameterError(
                    'Specifying more than two simultaneous parameters is illegal (except for T,p,x).')
        
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
                    
                # Deal with the special case that T,p,x are specified
                if x is not None:
                    x = np.asarray(x, dtype=float)
                else:
                    x = np.array(-1, dtype=float)
                if x.ndim==0:
                    x = np.reshape(x,(1,))
                    
                # Force compatible arrays
                T,p,x = np.broadcast_arrays(T,p,x)
                # Which points are not under the dome?
                I = x<0.
                
                # Initialize densities
                d1 = np.zeros_like(T)
                d2 = np.zeros_like(T)
                # deal with the non-saturated points
                d1[I] = self._d(T[I],p[I])
                d2[I] = d1[I]
                # Deal with densities under the dome
                I = np.logical_not(I)
                d1[I] = self._dsl(T[I])[0]
                d2[I] = self._dsv(T[I])[0]

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
            p = pm.units.pressure(
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
    
    def atoms(self):
        """Return a dictionary specifying the chemical composition of the substance.
    aa = atoms()
    
The dictionary keys are the symbols of atoms and their corresponding values 
are the integer quantities in the chemical formula.  For example
    aa = {'C':1, 'O':2}
would represent carbon dioxide.
"""
        aa = self.data.get('atoms')
        if aa is None:
            raise pm.utility.PMDataError('The substance does not have atomic composition data: ' + self.data['id'])
        return aa.copy()
    
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
    # General Properties            #
    #                               #
    
    def mw(self):
        """Molecular weight
    mw = mw()
    
Returns the molecular weight of the substance in
    [unit_mass / unit_molar]
"""
        mw = self.data['mw']
        mw = pm.units.mass(mw, from_units='kg')
        mw = pm.units.molar(mw, from_units='kmol', exponent=-1)
        return mw
    
    def R(self):
        """Ideal gas constant
    R = R()
    
Returns the ideal gas constant in
    [unit_energy / unit_matter / unit_temperature]
    
The mp1 data set includes a values for R lifted from the original data set.
The gas constant can be independently calculated from the universal gas 
constant or more precisely from the Boltzmann constant.  
    R = Ru / mw         # mw = molecular weight
        OR
    R = k * Na / mw     # Na = avagadro's number
    
The value returned by R is based on the value stored in the species data,
from which all other properties are constructed.
"""
        # R is stored in in J/kg/K
        R = pm.units.energy(self.data['R'], from_units = 'J')
        R = pm.units.matter(R, self.data['mw'], from_units='kg', exponent=-1)
        R = pm.units.temperature(R, from_units='K', exponent=-1)
        return R
        
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
        pm.units.temperature_scale(T, from_units='K', inplace=True)
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


    def T_s(self, s, p=None, d=None, quality=False, debug=False):
        """Temperature from entropy
    T = T_s(s, p=p)
        OR
    T = T_s(s, d=d)

If neither pressure nor density is specified, the default pressure will be 
used (config['def_p']).  

The optional keyword flag, quality, will cause quality to be returned
along with temperature.

    T,x = T_s(s, p=p, quality=True)
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
        if p is None and d is None:
            p = pm.config['def_p']

        ###############
        # If isobaric #
        ###############
        if p is not None:
            p = pm.units.pressure(
                    np.asarray(p, dtype=float),
                    to_units='Pa')
                    
            if p.ndim == 0:
                p = np.reshape(p, (1,))

            # Enforce pressure limits
            if ((p<self.data['plim'][0]).any() or
                    (p>self.data['plim'][1]).any()):
                raise pm.utility.PMParamError(
                    'MP1.T_s: Pressure is out-of-bounds.')

            # broadcast
            s,p = np.broadcast_arrays(s,p)
            # Initialize results
            T = np.empty_like(s, dtype=float)
            if quality:
                x = np.full_like(s, -1, dtype=float)
            # Some important intermediates
            Isat = np.zeros_like(s, dtype=bool)
            Ta = np.empty_like(s, dtype=float)
            Tb = np.empty_like(s, dtype=float)
            Tsat = np.empty_like(s, dtype=float)
            
            # Start with super-critical and sub-triple points
            I = np.logical_or(p >= self.data['pc'], p<=self.data['pt'])
            Ta[I] = self.data['Tlim'][0]
            Tb[I] = self.data['Tlim'][1]
            
            # Now work with the sub-critical points
            I = np.logical_not(I)
            if I.any():
                # Get the saturation temperatures
                Tsat[I] = self._Ts(p[I])
                # And densities (use ssL and ssV as temporaries)
                ssL = self._dsl(Tsat[I],0)[0]
                ssV = self._dsv(Tsat[I],0)[0]
                # finally, get the saturation entropies
                ssL = self._s(Tsat[I],ssL,0)[0]
                ssV = self._s(Tsat[I],ssV,0)[0]

                # Isolate points that are liquid
                Isat[I] = s[I] < ssL
                Ta[Isat] = self.data['Tlim'][0]
                Tb[Isat] = Tsat[Isat]*(1-1e-6)
                
                # Isolate points that are vapor
                Isat[I] = s[I] > ssV
                Ta[Isat] = Tsat[Isat]*(1+1e-6)
                Tb[Isat] = self.data['Tlim'][1]

                # Finally, isolate points that are saturated
                Isat[I] = np.logical_and( s[I]<=ssV, s[I]>=ssL )
                Ta[Isat] = Tsat[Isat]
                Tb[Isat] = Tsat[Isat]
                T[Isat] = Tsat[Isat]
                if quality:
                    x[Isat] = (s[Isat] - ssL)/(ssV-ssL)

            # Isat is now a down-select array
            Isat = np.logical_not(Isat)
    #        self._iter1(
            self._hybrid1(
                    self._tpiter,
                    'T',
                    s,
                    T,
                    Isat,
                    Ta, Tb,
                    param={'fn':self._s, 'p':p, 'debug':debug},
                    verbose=debug,
                    Nmax=50)
                    
        ################
        # If isochoric #
        ################
        elif d is not None:
            d = pm.units.matter(
                    np.asarray(d, dtype=float),
                    self.data['mw'],
                    to_units='kg')
            d = pm.units.volume(d, to_units='m3', exponent=-1)
                    
            if d.ndim == 0:
                d = np.reshape(d, (1,))

            # Enforce density limits
            if ((d<self.data['dlim'][0]).any() or
                    (d>self.data['dlim'][1]).any()):
                raise pm.utility.PMParamError(
                    'MP1.T_s: Density is out-of-bounds.')
            
            # set up the iteration
            # broadcast
            s,p = np.broadcast_arrays(s,p)
            # Initialize results
            T = np.empty_like(s, dtype=float)
            x = -np.ones_like(s,dtype=float)
            # Some important intermediates
            Isat = np.ones_like(s, dtype=bool)
            Ta = np.full_like(s, self.data['Tlim'][0], dtype=float)
            Tb = np.full_like(s, self.data['Tlim'][1], dtype=float)
            
    #        self._iter1(
            self._hybrid1(
                    self._tditer,
                    'T',
                    s,
                    T,
                    Isat,
                    Ta, Tb,
                    param={'fn':self._s, 'd':d, 'debug':debug},
                    verbose=debug,
                    Nmax=50)
            
            # Finally, reconstruct quality.  This is somewhat inefficient 
            # since saturation and quality were already calculated as part
            # of the iteration, but getting these intermediate parameters is
            # more irritating than it's worth for now.
            # Detect points with sub-critical temperatures.
            if quality:
                I = T < self.data['Tc']
                # Re-purpose Ta and Tb to hold the saturation densities
                Ta[I],_,_ = self._dsl(T[I])   # Liquid
                Tb[I],_,_ = self._dsv(T[I])   # Vapor
                I[I] = np.logical_and(d[I] >= Tb[I], d[I] <= Ta[I])
                if I.any():
                    x[I] = (Ta[I]/d[I] - 1.)/(Ta[I]/Tb[I] - 1.)
                
        # Convert the results
        pm.units.temperature_scale(T, from_units='K', inplace=True)
        
        if quality:
            return T,x
        return T


    def d_s(self, s, T=None, quality=False, debug=False):
        """Density from entropy
    d = d_s(s,T=T)
    
If temperature is not specified, the default temperature will be used 
(config['def_T']).

The optional keyword flag, quality, will cause quality to be returned along
with pressure.
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
        if T is None:
            T = pm.config['def_T']

        ###########################################
        # Isothermal (isobaric not yet supported) #
        ###########################################

        T = pm.units.temperature_scale(np.asarray(T, dtype=float), to_units='K')
        if T.ndim == 0:
            T = np.reshape(T, (1,))
            
        # Enforce temperature limits
        if ((T<self.data['Tlim'][0]).any() or
                (T>self.data['Tlim'][1]).any()):
            raise pm.utility.PMParamError(
                'MP1.p_s: Temperature is out-of-bounds.')

        # broadcast
        s,T = np.broadcast_arrays(s,T)
        # Initialize results
        d = np.empty_like(s, dtype=float)
        if quality:
            x = np.full_like(s, -1, dtype=float)
        # Some important intermediates
        da = np.empty_like(s, dtype=float)
        db = np.empty_like(s, dtype=float)
        ssL = np.empty_like(s, dtype=float)
        ssV = np.empty_like(s, dtype=float)
        p = np.empty_like(s, dtype=float)
        Isat = np.zeros_like(s, dtype=bool)
        
        # for points with sub-critical temperatures, we will test for
        # saturation
        I = T < self.data['Tc']
        if I.any():
            # Use upper and lower density limits as temporaries
            # da=liquid, db=vapor
            da[I] = self._dsl(T[I])[0]
            db[I] = self._dsv(T[I])[0]
            ssL[I] = self._s(T[I], da[I])[0]
            ssV[I] = self._s(T[I], db[I])[0]

            # Start with liquid points. 
            Isat[I] = s[I] < ssL
            #da[I] = da[I]  # the lower bound is already correct; sat.liq.
            db[Isat] = self.data['dlim'][1]
            
            # Now, move on to vapor points
            Isat[I] = s[I] > ssV
            da[Isat] = self.data['dlim'][0]
            #db[I] = db[I]  # the upper bound is already correct; sat.vap.
            
            # Finally, deal with points that are saturated
            Isat[I] = np.logical_and(ssL <= s[I], s[I] <= ssV)
            # Isat now indicates the saturated points
            # We already know the solution there; no iteraiton needed
            x[Isat] = (s[Isat]-ssL[Isat])/(ssV[Isat]-ssL[Isat])
            d[Isat] = 1./ ( 1./da[Isat] * (1-x[Isat]) + 1./db[Isat] * x[Isat])

        # Isat is now True for all non-saturated points
        Isat = np.logical_not(Isat)
        
        # For most data sets, the minim density is zero, but that is not
        # actually a legal value.  If necessary, iterate on the lower 
        # bound until the solution is bracketed
        I = (da == 0.)
        if I.any() and Isat.any():
            # Adjust the density a bit
            da[I] = np.full_like(s, .9999 * self.data['dlim'][0] + .0001 * self.data['dlim'][1], dtype=float)
            # calculate the entropy at the minimum and test to ensure inclusion
            s_ = self._s(T=T[I],d=da[I])[0]
            I[I] = s_ < s[I]
            count = 0
            while I.any():
                # dump out if this has been going on for too long
                count += 1
                if count > 20:
                    raise pm.utility.PMAnalysisError('Failed to find a lower-density bound to bracket a solution.')
                da[I] = 0.5*(self.data['dlim'][0] + da[I])
                s_[I] = self._s(T=T[I], d=da[I])[0]
                I[I] = s_[I] < s[I]
        
        self._hybrid1(
                self._s,
                'd',
                s,
                d,
                Isat,
                da, db,
                param={'T':T},
                verbose=debug,
                fx_index=2,
                Nmax=50)

        pm.units.matter(self.data['mw'], d, from_units='kg', inplace=True)
        pm.units.volume(d, from_units='m3', inplace=True, exponent=-1)
        if quality:
            return d,x
        return d
            



    def T_h(self, h, p=None, d=None, quality=False, debug=False):
        """Temperature from entropy
    T = T_s(s, p=p)
        OR
    T = T_s(s, d=d)

If neither pressure nor density is specified, the default pressure will be 
used (config['def_p']).  

The optional keyword flag, quality, will cause quality to be returned
along with temperature.

    T,x = T_s(s, p=p, quality=True)
"""
        # Prepare the s array
        h = pm.units.energy(
                np.asarray(h,dtype=float),
                to_units='J')
        pm.units.matter(h, self.data['mw'],
                to_units='kg', exponent=-1, inplace=True)

        if h.ndim == 0:
            h = np.reshape(h, (1,))

        # Set a default pressure?
        if p is None and d is None:
            p = pm.config['def_p']

        ###############
        # If isobaric #
        ###############
        if p is not None:
            p = pm.units.pressure(
                    np.asarray(p, dtype=float),
                    to_units='Pa')
                    
            if p.ndim == 0:
                p = np.reshape(p, (1,))

            # Enforce pressure limits
            if ((p<self.data['plim'][0]).any() or
                    (p>self.data['plim'][1]).any()):
                raise pm.utility.PMParamError(
                    'MP1.T_s: Pressure is out-of-bounds.')

            # broadcast
            h,p = np.broadcast_arrays(h,p)
            # Initialize results
            T = np.empty_like(h, dtype=float)
            if quality:
                x = np.full_like(h, -1, dtype=float)
            # Some important intermediates
            Isat = np.zeros_like(h, dtype=bool)
            Ta = np.empty_like(h, dtype=float)
            Tb = np.empty_like(h, dtype=float)
            Tsat = np.empty_like(h, dtype=float)
            
            # Start with super-critical and sub-triple points
            I = np.logical_or(p >= self.data['pc'], p<=self.data['pt'])
            Ta[I] = self.data['Tlim'][0]
            Tb[I] = self.data['Tlim'][1]
            
            # Now work with the sub-critical points
            I = np.logical_not(I)
            if I.any():
                # Get the saturation temperatures
                Tsat[I] = self._Ts(p[I])
                # And densities (use ssL and ssV as temporaries)
                hsL = self._dsl(Tsat[I],0)[0]
                hsV = self._dsv(Tsat[I],0)[0]
                # finally, get the saturation entropies
                hsL = self._h(Tsat[I],hsL,0)[0]
                hsV = self._h(Tsat[I],hsV,0)[0]

                # Isolate points that are liquid
                Isat[I] = h[I] < hsL
                Ta[Isat] = self.data['Tlim'][0]
                Tb[Isat] = Tsat[Isat]*(1-1e-6)
                
                # Isolate points that are vapor
                Isat[I] = h[I] > hsV
                Ta[Isat] = Tsat[Isat]*(1+1e-6)
                Tb[Isat] = self.data['Tlim'][1]

                # Finally, isolate points that are saturated
                Isat[I] = np.logical_and( h[I]<=hsV, h[I]>=hsL )
                Ta[Isat] = Tsat[Isat]
                Tb[Isat] = Tsat[Isat]
                T[Isat] = Tsat[Isat]
                if quality:
                    x[Isat] = (h[Isat] - hsL)/(hsV-hsL)

            # Isat is now a down-select array
            Isat = np.logical_not(Isat)
    #        self._iter1(
            self._hybrid1(
                    self._tpiter,
                    'T',
                    h,
                    T,
                    Isat,
                    Ta, Tb,
                    param={'fn':self._h, 'p':p, 'debug':debug},
                    verbose=debug,
                    Nmax=50)
                    
        ################
        # If isochoric #
        ################
        elif d is not None:
            d = pm.units.matter(
                    np.asarray(d, dtype=float),
                    self.data['mw'],
                    to_units='kg')
            d = pm.units.volume(d, to_units='m3', exponent=-1)
                    
            if d.ndim == 0:
                d = np.reshape(d, (1,))

            # Enforce density limits
            if ((d<self.data['dlim'][0]).any() or
                    (d>self.data['dlim'][1]).any()):
                raise pm.utility.PMParamError(
                    'MP1.T_s: Density is out-of-bounds.')
            
            # set up the iteration
            # broadcast
            h,p = np.broadcast_arrays(h,p)
            # Initialize results
            T = np.empty_like(h, dtype=float)
            x = -np.ones_like(h,dtype=float)
            # Some important intermediates
            Isat = np.ones_like(h, dtype=bool)
            Ta = np.full_like(h, self.data['Tlim'][0], dtype=float)
            Tb = np.full_like(h, self.data['Tlim'][1], dtype=float)
            
    #        self._iter1(
            self._hybrid1(
                    self._tditer,
                    'T',
                    h,
                    T,
                    Isat,
                    Ta, Tb,
                    param={'fn':self._h, 'd':d, 'debug':debug},
                    verbose=debug,
                    Nmax=50)
            
            # Finally, reconstruct quality.  This is somewhat inefficient 
            # since saturation and quality were already calculated as part
            # of the iteration, but getting these intermediate parameters is
            # more irritating than it's worth for now.
            # Detect points with sub-critical temperatures.
            if quality:
                I = T < self.data['Tc']
                # Re-purpose Ta and Tb to hold the saturation densities
                Ta[I],_,_ = self._dsl(T[I])   # Liquid
                Tb[I],_,_ = self._dsv(T[I])   # Vapor
                I[I] = np.logical_and(d[I] >= Tb[I], d[I] <= Ta[I])
                if I.any():
                    x[I] = (Ta[I]/d[I] - 1.)/(Ta[I]/Tb[I] - 1.)
                
        # Convert the results
        pm.units.temperature_scale(T, from_units='K', inplace=True)
        
        if quality:
            return T,x
        return T
