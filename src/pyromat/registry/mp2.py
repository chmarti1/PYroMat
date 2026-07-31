# MP1
#   PYroMat Multi-phase generalist class
#   Calculates physical properties from a fit for the helmholtz free 
#   energy in terms of density and temperature.


import numpy as np
from numpy import linalg as la
import pyromat as pm
import os,sys


#
# Helper Functions
#

def ndxgen(xc, N, r):
    """Generate a 1D array of N values in [0,1] with density r about xc
    x, ci = ndxgen(xc, N, r)
    
Constructs an array of N data increasing from 0 to 1 with a relative 
density, r, at 0 < xc < 1.  The nominal density of data is N, so the 
relative denstiy, r*N.  

x   The array of values
ci  The index corresponding precisely to xc.
"""

    # Construct a dimensionless piece-wise fit of two quadratics
    # joined at the critical point (d-less, xc)
    Nc = int(xc * N)
    N1 = N-1
    
    A = np.matrix([[ 0, 0, 1],
                   [ Nc*Nc, Nc, 1],
                   [ 2*Nc, 1, 0]], dtype=float)
    B = np.array([0, xc, 1./r/N1])
    c1 = np.linalg.solve(A,B)

    A = np.matrix([[N1*N1, N1, 1],
                   [ Nc*Nc, Nc, 1],
                   [ 2*Nc, 1, 0]], dtype=float)
    B = np.array([1., xc, 1./r/N1])
    c2 = np.linalg.solve(A,B)

    # Generate the array
    x = np.empty(N, dtype=float)
    x[:Nc] = np.polyval(c1, np.arange(0,Nc))
    x[Nc:] = np.polyval(c2, np.arange(Nc,N))
    return x, Nc

def interp_scalar(x, x0, x1, f0, f1):
    """Perform 1D linear interpolation between two explicitly provided points
    f = interp_scalar(x, x0, x1, f0, f1)
    
Performs the standard linear interpolation on the line segment between 
the point pairs, (x0,f0) to (x1,f1), using the equation
    f = (x - x0)/(x1 - x0) * (f1 - f0) + f0

Detects precise equality with x0 or x1 to return precisely f0 or f1
respectively.
"""
    # If x lies precisely on the nodes, return the f values precisely
    # Strangely, single-line if-else structures give about a 10% speedup
    # even if they are less readable
    return f0 if x == x0 else f1 if x == x1 else f0 + (f1-f0)*(x-x0)/(x1-x0)


def crossing2(I):
    """Find elements of a 2D boolean grid with dissimilar corners for _mapsearch2
    
    J = crossing2(I)

If I is an m by n 2D array of boolean values, representing a logical 
test values at the nodes of a 2D data map, crossing2() identifies the
elements with at least one False and at least one True node.

    False   
        +---+ True
        |   |
   True +---+
            True
            
When I.shape is (m,n), J.shape is (m-1, n-1).

This is relegated to a helper function because it needs to be performed
twice by _mapsearch2() -- once on fdata and once on gdata.
"""
    a = I[:-1,:-1]      # Lower-left
    b = I[:-1,1:]       # Lower-right
    c = I[1:,:-1]       # Upper-left
    d = I[1:,1:]        # Upper-right
    # not a*b*c*d is true iff at least one node is False
    # a+b+c+d is true iff at least one node is True
    return (~(a*b*c*d)) * (a+b+c+d)



class mp2(pm.reg.__basedata__):
    """The PYroMat multi-phase generalist class 1

** Available Property Methods **
MP2 provides property methods:
    a()     Speed of sound
    cp()    Isobaric specific heat
    cv()    Isochoric specific heat
    gam()   Specific heat ratio
    e()     Internal energy
    f()     Free (Helmholtz) energy
    g()     Gibbs energy
    h()     Enthalpy
    s()     Entropy
    T()     Temperature
    p()     Pressure
    d()     Density
    v()     Specific volume
    x()     Quality
    state() Calculates all properties
    
All of the above methods accept a standardized call signature, which 
accepts any of the following arguments: T, p, d, v, e, h, s, x

For example, enthalpy might be called
    h(T=300., p=1.01325)
    h(T=300., d=990.)
    h(T=300., x=0.5)
    h(s=6., p=2.5)

In the back end, all properties are calculated from temperature and density,
so providing this interface flexibility has a numerical cost.  Once T and d
are known, additional property evaluations should always be made in terms
of them.

Most property pairs are supported, but several are not.  For example, e,
s, and h must be specified with a "basic" property; T, d, p, v, or x.  
This limitation is to prevent the costly numerical iteration that occurs
when two "higher" properties need to be simultaneously inverted.  

Furthermore, since it is impossible to specify a saturated mixture with
temperature and pressure alone, there is a special case, which permits 
three properties: T, p, x.  When x is negative, it is ignored, but for all
points where it is between 0 and 1, pressure is ignored, and presumed to
be the saturation pressure at the specified temperature.  For performance
reasons, this condition is not tested, so if it is violated, an error will
not be raised.

Most property methods also accept the "quality" as an optional keyword.  When
it is set to True, the property will also return the vapor/liquid mixture
quality in a tuple with the property value.  For example,
    h,x = h(T,d,quality=True)

** Saturation Properties **
There are also saturation property methods:
    es()    Saturation internal energy
    hs()    Saturation enthalpy
    ss()    Saturation entropy
    
And saturation equations of state methods:
    Ts()    Saturaiton temperature
    ds()    Saturation densities
    vs()    Saturation specific volumes
    ps()    Saturation pressure

Saturation methods accept either temperature or pressure as an argument.  
The density saturation method returns both liquid and vapor densities in a 
tuple pair.  See their in-line documentaiton for more details.

** Other Properties **
There are other methods that return useful information, but that do 
not depend on the state.
    atoms() Returns a dictionary specifying the chemical composition.
    mw()    Returns the molecular weight/mass
    R()     Returns the ideal gas constant
    Tlim()  Returns [Tmin, Tmax] valid temperature range
    plim()  Returns [pmin, pmax] valid pressure range
    critical()  Returns the state at the critical point
    triple()    Returns the state at the triple point

*** MORE DOCUMENTATION ***
MP1 models thermo-physical properties of a liquid-gas system using a 
general fit for helmholtz free energy.  These "Span & Wagner" fits are 
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
The MP2 data dictionary must have certain data "groups" to define the 
various empirical fits.  Each group is a dictionary (within the 
dictionary) that defines the various parameters necessary for at least
one of the inner methods.

AOgroup        Helmholtz free energy ideal gas group; a dict containing:
    Tscale      Temperature scale for normalizing T
    dscale      density scale for normalizing d
    logt        a scalar coefficient of a log(tt) term
    coef0       a coefficient list to be passed to _poly1() to build p0 below
    coef1       a simple Nx2 coefficient list used to build q(tt) below
If tt = Tscale/T    <=== INVERSE!
and dd = d/dscale
    ao = log(d) + LOGT*log(tt) + TLOGT*tt*log(tt) + p0(tt) + q(tt)
        q(tt) = sum_k coef1[k,1] * log(1 - exp(-tt*coef[k,0]))
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

There are optional tabular data elements that allow designers to 
explicitly store tabular data.  If they are omitted, the data will be
automatically generated from AOgroup and ARgroup.  The properties of
the automatically generated table can also be specified.

Tdata           A 1D list with (m) elements specifying the temperatures 
                of the table entries.  Units must be Kelvin.
NT              An integer specifying the number of temperature data 
                to automatically generate if 'Tdata' is absent. 
                Defaults to 101 if absent.
rT              The relative density of temperature data near the
                critical point.  Defaults to 2 if absent.
ddata           A 1D list with (n) elements specifying the densities 
                for the table entries.  Units must be kg/m^3.
Nd              An integer specifying the number of temperature data 
                to automatically generate if 'Tdata' is absent. 
                Defaults to 101 if absent.
rd              The relative density of temperature data near the
                critical point.  Defaults to 2 if absent.
hdata           A 2D list with (m x n) entries of enthalpy evaluated at
                h(Tdata, ddata).  hdata may not be specified if either
                Tdata or ddata were not specified.  Units must be J/kg.
sdata           A 2D list with (m x n) entries of entropy evaluated at
                s(Tdata, ddata).  sdata may not be specified if either
                Tdata or ddata were not specified.  Units must be J/kg/K.
                
Tsdata, psdata, dsLdata, dsVdata
    1D lists specifying

Additionally, there are a number of parameters that define global 
properties (true at all states)

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
                have atoms = {'C':1, 'O':2}
                
There are also the typical mandatory PYroMat meta data elements:
id              What substance is this?
doc             Where did it come from?
class           What class should be used to evaluate the data?
"""

    def __init__(self, *arg, **kwarg):
        # Call the basedata class initializer
        pm.reg.__basedata__.__init__(self, *arg, **kwarg)

            
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
            if test in self.data['IGgroup']:
                pass
            elif self.data['IGgroup'][test] > 0:
                pass
            else:
                error = True
                break
        if error:
            report.write('[FAILED]    0.1 AOgroup must contain positive scalars, Tscale, dscale\n')
            report.write('            ' + test + ' = ' + repr(self.data['IGgroup'][test]) + '\n')
            result = False
        else:
            report.write('[passed]    0.1 AOgroup must contain positive scalars, Tscale, dscale\n')
            
        #0.2 AOgroup's logt parameter must be a scalar
        if 'logt' in self.data['IGgroup'] and not isinstance(self.data['IGgroup']['logt'], (float,int)):
            report.write('[FAILED]    0.2 AOgroup logt parameter must be a scalar\n')
            report.write('            logt = ' + repr(self.data['IGgroup']['logt']) + '\n')
            result = False
        else:
            report.write('[passed]    0.2 AOgroup logt parameter must be a scalar\n')
        
        #0.3 AOgroup's coef0 parameter must be a legal poly1 group
        error = False
        if 'coef0' in self.data['IGgroup']:
            error,message = _check_poly1(self.data['IGgroup']['coef0'])
        if error:
            report.write('[FAILED]    0.3 AOgroup coef0 parameter must be a legal poly1 group\n')
            report.write('            ' + message + '\n')
            result = False
        else:
            report.write('[passed]    0.3 AOgroup coef0 parameter must be a legal poly1 group\n')
        
        #0.4 AOgroup's coef1 parameter must be a table with two columns
        error = False
        if 'coef1' in self.data['IGgroup']:
            if not isinstance(self.data['IGgroup']['coef1'], (list,tuple)):
                error = True
                message = 'coef1 was not iterable.'
            else:
                for row in self.data['IGgroup']['coef1']:
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
            if test in self.data['Rgroup']:
                pass
            elif self.data['Rgroup'][test] > 0:
                pass
            else:
                error = True
                break
        if error:
            report.write('[FAILED]    0.5 ARgroup must contain positive scalars, Tscale, dscale\n')
            report.write('            ' + test + ' = ' + repr(self.data['Rgroup'][test]) + '\n')
            result = False
        else:
            report.write('[passed]    0.5 ARgroup must contain positive scalars, Tscale, dscale\n')
            
        #0.6 ARgroup's coef0 must be a list of legal poly2 groups
        error = False
        message = ''
        if 'coef0' in self.data['Rgroup']:
            if not isinstance(self.data['Rgroup']['coef0'], (list,tuple)):
                error = True
                message = 'coef0 was not iterable.'
            else:
                for cgi,cg in enumerate(self.data['Rgroup']['coef0']):
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
        if 'coef1' in self.data['Rgroup']:
            if not isinstance(self.data['Rgroup']['coef1'], (list,tuple)):
                error = True
                message = 'coef1 was not iterable.'
            else:
                for row in self.data['Rgroup']['coef1']:
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
        if 'coef2' in self.data['Rgroup']:
            if not isinstance(self.data['Rgroup']['coef2'], (list,tuple)):
                error = True
                message = 'coef2 was not iterable.'
            else:
                for row in self.data['Rgroup']['coef2']:
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
        


    def _poly2(self,x,y,pcoef,diff=2):    
        """Polynomial evaluation (primative routine)
(p, px, py, pxx, pxy, pyy) = _poly(x,y,pcoef,diff=2)

Evaluates a polynomial on x and y and its derivatives.
x       x value
y       y value
pcoef   coefficient dictionary/list
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
The coefficient lists defining the terms must contain three elements:
[<i>, <j>, <c>].  The coefficients must be sorted by x-power and
then by y-power in descending order.  The powers must be non-negative
integers.

A coefficient dictionary might appear

coef = {
    'pre': [<xpre>, <ypre>],
    'post': [<xpost>, <ypost>],
    'coef':[
        [<i>, <j>, <c>],
        ...
    ]
]

The pre-exponents are applied to the arguments to the polynomial, and
the post-exponents are applied after the polynomial is evaluated, so 
that the value returned is
    X = x**xpre
    Y = y**ypre
    p = c00 + ... cij * X**i * Y**j + ...
    output = x**xpost * y**post * p

For example, the polynomial,
    p(x,y) = .5 + 1.2y + .2y**2 + 0.1xy
might be represented by the dict:
pcoef = {
    'coef':[
        [1, 0, 0.1],
        [0, 2, 0.2],
        [0, 1, 1.2],
        [0, 0, 0.5]
    ]
}

In this example the 'pre' and 'post' values are omitted because they are
not necessary.

Efficient polynomial evaluation algorithms are normally restricted to
positive integer exponents, but many thermodynamic property models use 
much more interesting polynomials.  The pre- and post- exponents can be
used to acheive a much wider range of functions.

For example,
    p(x,y) = x**(-1.5) + x**(3.5)
might be expressed with the dictionary
pcoef = {
    'pre': [0.5, 1],
    'post':[-1.5, 0],
    'coef':[
        [10, 0, 1],
        [0, 0, 1]
    ]
}

which is equivalent to the original polynomial, except that the core of
the evaluation algorithm only operates on positive integers.
"""
        if isinstance(pcoef, list):
            g = 0.  # total group
            gx = 0.
            gy = 0.
            gxx = 0.
            gxy = 0.
            gyy = 0.
            
            for this in pcoef:
                p,px,py,pxx,pxy,pyy = self._poly2(x,y,this,diff)
                
                g += p
                gx += px
                gy += py
                gxx += pxx
                gxy += pxy
                gyy += pyy
                
            return g,gx,gy,gxx,gxy,gyy

        # initialize the polynomial and its derivatives
        p = 0.  # total polynomial
        px = 0.
        py = 0.
        pxx = 0.
        pxy = 0.
        pyy = 0.
        
        # collect the pre-exponentials
        if 'pre' in pcoef:
            prex,prey = pcoef['pre']
        else:
            prex = 1
            prey = 1
        # Apply the pre-exponentials
        if prex != 1:
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
            
        if prey!=1:
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
        coef = pcoef['coef']
        ncoef = len(coef)
        imax = coef[0][0]
        
        # On which coefficient are we currently operating?
        index = 0
        # This is a flag that indicates the active index was used in 
        # the last loop, so it needs to be incremented.
        
        for ii in range(imax,-1,-1):
            # If the current x-exponent is the same one represented in
            # the active coefficient row, then calculate q.
            if index<ncoef and coef[index][0] == ii:
                # For this value of ii, what is the largest jj?
                jmax = coef[index][1]
                # q is a sub-polynomial on y that represents the 
                # variation on y of all terms that share the same
                # power in x.  This inner loop is much like the loop
                # on ii, except that it looks at both the x and y 
                # exponents.
                q = 0
                qy = 0
                qyy = 0
                for jj in range(jmax,-1,-1):
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
        if prex!=1 or prey!=1:
            if diff>0:
                if diff>1:
                    pxx = pxx*x_1*x_1 + px*x_2
                    pyy = pyy*y_1*y_1 + py*y_2
                    pxy = pxy*x_1*y_1
                px *= x_1
                py *= y_1
        
        # Collect the post-exponentials
        if 'post' in pcoef:
            postx,posty = pcoef['post']
        else:
            postx = 0
            posty = 0
        
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
        
        return p,px,py,pxx,pxy,pyy


    def _poly1(self,x,pcoef,diff=2):    
        """Polynomial evaluation (primative routine)
(p, px, pxx) = _poly1(x,pcoef,diff=2)

Evaluates a polynomial on x and y and its derivatives.
x       x value
pcoef   coefficient list/dictionary
diff    the highest order derivative to evaluate (0,1, or 2)

Returns
p       polynomial value at p(x)
px      dp/dx
pxx     d2p/dx2

When diff is less than 2, the corresponding values of px and pxx are
returned as 0.  The default for diff is 2 to protect against careless
treatment as if these values ARE zero, but reducing diff will make poly1
execute more efficiently.

The pcoef parameter is either a dictionary or a list of dictionaries
representing a polynomial series.  Each dictionary defines a polynomial
with optional pre- and post-exponents, defining a polynomial of the form
    X = x**pre
    p = c0 + c1*X + c2*X**2 + ... cn*X**n + ...
    output = x**post * p

When the pcoef parameter is a list of dictionaries, each dictionary 
defines one of these polynomials, and the results are summed together.

Each dictionary is of the form:
{
    'pre':<pre>,
    'post':<post>,
    'coef':[
        [<n>, <c>],
        [<n>, <c>],
        ...
    ]
}
<n> and <c> represent the exponent and corresponding coefficient.  The
coefficient list is sparse, so zero-value coefficients are simply 
omitted.  THE LIST MUST BE IN DECENDING ORDER OF EXPONENTS.

If the pre- and post-terms are omitted, then they are ignored.  This is
equivalent to pre=1 and post=0.

In a simple example, the polynomial,
    p(x) = 2*x**-1.5 - x**0.5
    
might be specified
pcoef = {
    'pre': 0.5,
    'post': -1.5,
    'coef': [[4,-1], [0,2.]]
}
"""
        if isinstance(pcoef, list):
            g = 0.
            gx = 0.
            gxx = 0.
            for this in pcoef:
                p,px,pxx = self._poly1(x,this,diff=diff)
                g += p
                gx += px
                gxx += pxx
            return g,gx,gxx

    
        # initialize the final polynomial and its derivatives
        p = 0.  # total polynomial
        px = 0.
        pxx = 0.
       
        # Apply the pre-exponentials
        if 'pre' in pcoef and pcoef['pre'] != 1:
            pre = pcoef['pre']
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
            pre = 1.
            x_0 = x
            x_1 = 1.
            x_2 = 0.

        # From here, we loop over terms of the form a*(x**ii)
        # If a particular ii,jj combination is not found in the data, then
        # its coefficient is treated as zero.
        # What is the largest ii?
        coef = pcoef['coef']
        imax = coef[0][0]
        ncoef = len(coef)
        
        # On which coefficient are we currently operating?
        index = 0
        # Loop through all polynomial powers
        for ii in range(imax,-1,-1):
            # If the current x-exponent is the same one represented in
            # the active coefficient row, then calculate q.
            if index<ncoef and coef[index][0] == ii:
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
        if 'post' in pcoef and pcoef['post'] != 0:
            post = pcoef['post']
            f = x**post
            if diff>0:
                fx = post*f/x
                if diff>1:
                    fxx = fx*(post-1)/x
                    pxx = pxx*f + 2.*px*fx + p*fxx
                px = px*f + p*fx
            p *= f

        return p,px,pxx


    def _mapsearch1(self, xdata, fdata, fvalue=0):
        """Search 1D map for an inverse estimates (primative routine)
    x = mapsearch1(xdata, fdata)
        OR
    x = mapsearch1(xdata, fdata, fvalue)
    
Uses tabulated data to generate estimates for x in the 1D inversion 
problem
    f(x) = fvalue

The fdata is a 1D array of tabulated values of f(x) with an identically
sized array of corresponding xdata.  This is notably distinct from 1D 
interpolation because the fdata map does not need to be monotonically 
increasing.  Instead, the algorithm performs a global search by 
explicitly comparing all node values,
    fvalue < f_i
Because this is vectorized and performed in the compiled Numpy back-end,
exhaustive explicit search in tabulated data is significantly faster
than iterative calculation using the full Span and Wagner models, and
cost is roughly linear with data set size.

Elements containing a solution are identified when one node is below or
equal to the value, and the other is greater than the value.  The 
solution estimate is extrated by linear interpolation.

SEE ALSO:
    _mapsearch2()
"""
        I = fvalue < fdata
        xi = np.nonzero(I[:-1] ^ I[1:])[0]
        xi1 = xi+1
        x = xdata[xi] + (xdata[xi1] - xdata[xi]) * (fvalue - fdata[xi]) / (fdata[xi1] - fdata[xi])
        return x

    def _mapsearch2(self, xdata, ydata, fdata, gdata, fvalue, gvalue):
        r"""Search 2D map for inverse estimates (primative routine)
    x,y,xi,yi = mapsearch2(xdata, ydata, fdata, gdata, fvalue, gvalue)
    
Uses tabulated data to generate an estimates for x,y in the 2D inversion
problem
    f(x,y) = fvalue
    g(x,y) = gvalue

ARGUMENTS:
xdata, ydata
    One-dimensional array-like containing grid values for the x- and y-
    coordinates.  The sizes of the x- and y-data arrays must match the 
    fdata and gdata arrays (see below).
    
fdata, gdata
    Two-dimensional array-like containing tabulated values for f(x,y) 
    and g(x,y).  The indices should be arranged so that
        fdata[i,j] = f(xdata[i], ydata[j])
        gdata[i,j] = g(xdata[i], ydata[j])
        
fvalue, gvalue
    Numpy arrays with the same shape containing values for properties,
    fdata and gdata.

    
RETURNS: 
x,y
    Arrays of the same shape as fvalue and gvalue that approximate 
    solutions to the problem
        f(x,y) =approx= fvalue
        g(x,y) =approx= gvalue

xi,yi
    One-dimensional integer arrays identifying the indices of the 
    elements in which the estimated solution was identified.  Care must
    be taken, because the actual solution may lie in a neighboring 
    element -- especially when estimates are very near the element edge.

DESCRIPTION:

The fdata and gdata are 2D arrays of tabulated values of f(x,y) and 
g(x,y) in a rectangular grid of x and y values.  This is notably 
distinct from 2D interpolation because the maps, fdata and gdata, do not
need to be monotonically increasing.  The algorithm performs a global 
search by explicitly comparing all node values:
    fvalue < f_ij
    gvalue < g_ij

Grid elements containing potential solutions are identified as those 
with at least one node above and below the target values for both f() 
and g().  Then estimates are generated by finding the approximate 
intersections of the paths in x,y implied by the f() and g() 
constraints inside the element.  First, the element's edges are 
interpolated to find estimates for two points where f(x,y)==fvalue and
g(x,y)==gvalue.  The intersection (if one exists) of the two resulting
line segments is interpreted as the estimated solution.

    +--x----+           +---x---+
    |  |    |           |  /    |
    |  \_,.-o           | /   ,-o
    o-' |   |           x'   /  |
    +---x---+           +---o---+
    Intersection        No Intersection

SPEED AND STABILITY:

Because this is vectorized and performed in the compiled Numpy back-end,
exhaustive explicit search in tabulated data is significantly faster
than iterative calculation using the full Span and Wagner models. Tests
using a single core of an AMD Ryzen 9 7900 show that property evaluation
is roughly equivalent with 5e6 (five million) floating point comparisons.
On RISC systems, without machine-level vectorized comparison operations,
vector comparison may be significantly slower, but most maps only 
contain roughly 1e4 (ten thousand) elements.  In general, an iteration 
saved by a better initial guess is worth MANY bulk floating point 
comparisons.

Though speed is certainly a benefit, the real reason to use maps is to
produce initial guesses close enough to the actual solution so that the 
faster (and simpler) Newton-Rapson root polishing algorithm can be used 
without fear of numerical stability problems in higher dimensions. 
Beyond a minimum performance threshold, PYroMat values reliabile 
convergence and robust identification of all possible solutions more 
highly than speed.

The real limitation of _mapsearch2() is that its inputs are inherently 
scalar, allowing only one fvalue, gvalue pair at a time.  This means 
_mapsearch2() is vectorized by a for a loop to work on datasets, which 
always bodes poorly for performance.  Most users seem to use PYroMat on 
datasets smaller than the back-end maps, so it is better to vectorize 
the map search than to vectorize the value inputs. 

ABOUT SOLUTION SEGMENT INTERPOLATION:

Solution segment interpolation was selected over the usual bilinear 
interpolation because of its linearity.  Bilinear 2D element 
interpolation is obnoxious to invert because of its nonlinear xy term,
which can cause saddle points and other irritating issues.  However, 
solution segment interpolation still suffers from problems, which are 
mitigated in this algorithm:
(1) When the solution lies precisely on a node, one line segment 
    vanishes, leading to a singular problem.  This is mitigated by 
    explicitly testing for precise equality at the nodes.
(2) When solution estimates lie very close to the element edge, tiny 
    numerical errors can cause redundant estimates from neighboring
    elements or the estimate can be omitted altogether.  When estimates
    are a small distance from an element's edge (even if it is very 
    slightly outside) it is included.  If the corresponding neighboring
    element also appears as a candidate, it is deselected to prevent
    redundant reporting.
(3) When line segments are very nearly parallel, the intersection 
    problem  becomes singular.  The determinant of the 2x2 matrix is
    calculated in a separate step, and the process is halted if it is 
    too small to possibly generate a reasonable solution.  This approach
    also prevents wasting time calculating the intersections of lines 
    that clearly have no chance of intersecting inside the element.
(4) "Saddle" elements have diagonal nodes on one side of the target 
    value and anti-diagonal nodes on the other.  The interpolation of
    the two implied solution path segments is ambiguous, the existence
    of a solution is uncertain, and it is likely to be very nearly 
    singular.  For the purposes of PYroMat's numerical problems, these
    cases are detected and discarded.
    
A number of versions of _mapsearch2() were tested. This version simply 
returns the first solution discovered.  Other versions faithfully 
reported multiple candidate solutions if they were discovered.  Since
the top layer of PYroMat does not currently permit reporting multiple
solutions, this funcitonality was discarded.  It might be recovered in
later versions if it is needed.

SEE ALSO:
    _mapsearch1(), _mapsearch2(), _mapsearch2x(), _mapsearch2y()
"""
        # Define an increment for small values
        # For most systems, eps is about 2.2e-16, so small will be about
        # 2.2e-12.  This is the number we use to detect dimensionless
        # proximity to the element boundary.
        small = np.finfo(float).eps * 1e4
        # Initialize lists for the result values
        x = np.empty_like(fvalue, dtype=float)
        y = np.empty_like(fvalue, dtype=float)
        XI = np.empty_like(fvalue, dtype=int)
        YI = np.empty_like(fvalue, dtype=int)
        
        for index in range(fvalue.size):
            fv = fvalue.flat[index]
            gv = gvalue.flat[index]

            # Generate a boolean array indicating candidate elements with a solution
            # Bulk element comparison seems expensive, but it is not on a 
            # system with vectorized processing.  Bulk comparisons like this
            # are remarkably cheap. 
            fI = fv < fdata
            gI = gv < gdata
            I = crossing2(fI) * crossing2(gI)

            fail = True

            # For each element that contains a crossing in both f and g
            for xi,yi in zip(*np.nonzero(I)):
                # Only continue if this candidate is still flagged
                # Elements can be unflagged as the algorithm progresses if a neighbor
                # has claimed a point on the border or in the corner.
                # This conditional was removed when the code was modified to return the
                # first solution discovered.  Uncomment and indent if the code needs to
                # return multiple solutions.
                #<<==>>
                #if I[xi,yi]:
                # Indices for the other four nodes in this element
                xi1 = xi+1
                yi1 = yi+1
                # Identify the two f-edge crossings [(x,y), ...]
                fc = []
                # Track the indices of the neighboring elements in case the
                # solution is very near to the element's boundary.  Only 
                # the neighbors of the f-segment are tracked.
                # This code
                #<<==>>
                #neighbor = []
                # Test each of the edges for a crossing of f()
                # Bottom edge
                if fI[xi,yi] != fI[xi1,yi]:
                    xx = interp_scalar(fv, fdata[xi,yi], fdata[xi1,yi], xdata[xi], xdata[xi1])
                    fc.append(np.array((xx,ydata[yi])))
                    #<<==>>
                    #neighbor.append((xi, yi-1))
                # Left edge
                if fI[xi,yi] != fI[xi,yi1]:
                    yy = interp_scalar(fv, fdata[xi,yi], fdata[xi,yi1], ydata[yi], ydata[yi1])
                    fc.append(np.array((xdata[xi], yy)))
                    #<<==>>
                    #neighbor.append((xi-1, yi))
                # Top edge
                if fI[xi,yi1] != fI[xi1,yi1]:
                    xx = interp_scalar(fv, fdata[xi,yi1], fdata[xi1,yi1], xdata[xi], xdata[xi1])
                    fc.append(np.array((xx,ydata[yi1])))
                    #<<==>>
                    #neighbor.append((xi, yi+1))
                # Right edge
                if fI[xi1,yi] != fI[xi1,yi1]:
                    yy = interp_scalar(fv, fdata[xi1,yi], fdata[xi1,yi1], ydata[yi], ydata[yi1])
                    fc.append(np.array((xdata[xi1], yy)))
                    #<<==>>
                    #neighbor.append((xi+1, yi))
                # Identify the two g-edge crossings [(x,y), ...]
                gc = []
                # Test each of the edges for a crossing of g()
                # Bottom edge
                if gI[xi,yi] != gI[xi1,yi]:
                    xx = interp_scalar(gv, gdata[xi,yi], gdata[xi1,yi], xdata[xi], xdata[xi1])
                    gc.append(np.array((xx,ydata[yi])))
                # Left edge
                if gI[xi,yi] != gI[xi,yi1]:
                    yy = interp_scalar(gv, gdata[xi,yi], gdata[xi,yi1], ydata[yi], ydata[yi1])
                    gc.append(np.array((xdata[xi], yy)))
                # Top edge
                if gI[xi,yi1] != gI[xi1,yi1]:
                    xx = interp_scalar(gv, gdata[xi,yi1], gdata[xi1,yi1], xdata[xi], xdata[xi1])
                    gc.append(np.array((xx,ydata[yi1])))
                # Right edge
                if gI[xi1,yi] != gI[xi1,yi1]:
                    yy = interp_scalar(gv, gdata[xi1,yi], gdata[xi1,yi1], ydata[yi], ydata[yi1])
                    gc.append(np.array((xdata[xi1], yy)))
                # At this point, fc and gc list (x,y) coordinates for 
                # the points along the element edge where crossings occur
                # Meanwhile, neighbor lists the (xi,yi) indices of the
                # elements that share the edges where f() has a solution
                # We'll use neighbor to resolve conflict over solutions
                # very close to the edges.
                
                # Detect the saddle case
                if len(fc) != 2 or len(gc) != 2:
                    # For now, warn the user, and DO NOT append the case
                    pm.utility.print_warning('mp2._mapsearch2: Discarded a potential solution near a saddle point.  If you believe this was a legitimate solution, please report the code that generated this warning to the PYroMat GitHub issues page.')
                # Two edges have intersections for each function
                else:
                    fx0 = fc[0]
                    fdx = fc[1] - fc[0]
                    gx0 = gc[0]
                    gdx = gc[1] - gc[0]
                    #print('')
                    #print('xi,yi,x,y:', xi,yi,xdata[xi], ydata[yi])
                    #print('fdata values:', fdata[xi,yi], fdata[xi1,yi], fdata[xi,yi1], fdata[xi1,yi1])
                    #print('gdata values:', gdata[xi,yi], gdata[xi1,yi], gdata[xi,yi1], gdata[xi1,yi1])
                    
                    # Check for a solution precisely at the corner
                    if (fdx == 0).all():
                        if (gc[0] == fx0).all() or (gc[1] == fx0).all():
                            # When this code was modified to merely return the first solution discovered,
                            # these lines were commented out.  Return them if multiple solutions are 
                            # desired in the future.
                            #I[*neighbor[0]] = False
                            #I[*neighbor[1]] = False
                            x.flat[index] = fx0[0]
                            y.flat[index] = fx0[1]
                            XI.flat[index] = xi
                            YI.flat[index] = yi
                            break
                    # Ignore gdx == 0 cases - we'll catch corners with fdx == 0
                    elif not (gdx == 0).all():                        
                        # Solve for a dimensionless number, s
                        # The linear problem is
                        #   fdx * s + fx0 - gdx * r - gx0 = 0
                        # So, solving for scalars, r and s, leads to a matrix
                        #   A = [fdx  -gdx]
                        #   B = -fx0 + gx0
                        #   A * [s r]' = B
                        # Rather than use the solve algorithm, calculate the
                        # determinant explicitly to detect the very nearly
                        # parallel case
                        det = -fdx[0]*gdx[1] + fdx[1]*gdx[0]
                        B = -fx0 + gx0
                        # s without dividing by det yet
                        s = -B[0]*gdx[1] + B[1]*gdx[0]
                        # If the determinant is small, there is no need to keep going
                        if 2*abs(det) > abs(s):
                            s /= det
                            # If the solution lies in the element or very
                            # slightly outside of it, log the potential 
                            # solution.
                            if -small < s < 1+small:
                                # Store the solution
                                x.flat[index] = fx0[0] + s*fdx[0]
                                y.flat[index] = fx0[1] + s*fdx[1]
                                # If the solution is very near a boundary, remove
                                # the neighbor element as a candidate to prevent
                                # redundant solutions.
                                # This was removed when the code was modified to only
                                # return the first solution discovered.  Uncomment it
                                # if multiple solutions are desired in the future
                                #if -small < s < small:
                                #    I[*neighbor[0]] = False
                                #if 1-small < s < 1+small:
                                #    I[*neighbor[1]] = False
                                XI.flat[index] = xi
                                YI.flat[index] = yi
                                break
        return x,y,XI,YI


    def _mapsearch2x(self, xdata, ydata, fdata, yvalue, fvalue, indices=True):
        r"""Search 2D map for inverse estimates (primative routine)
    x, xi, yi = _mapsearch2x(xdata, ydata, fdata, yvalue, fvalue)
        OR
    x, xi, yi = _mapsearch2x(..., indices=False)
    
Uses tabulated data to generate an estimate for x in the 2D inversion
problem
    f(x,yvalue) = fvalue

ARGUMENTS:
xdata, ydata
    One-dimensional array-like containing grid values for the x- and y-
    coordinates.  The sizes of the x- and y-data arrays must match the 
    fdata and gdata arrays (see below).
    
fdata
    Two-dimensional array-like containing tabulated values for f(x,y).  
    The indices should be arranged so that
        fdata[i,j] = f(xdata[i], ydata[j])
        gdata[i,j] = g(xdata[i], ydata[j])
        
yvalue
    An array of y-values to interpolate from the table.  The dimensions
    must match the dimensions of fvalue.
    
fvalue
    An array of f-values to interpolate from the table.  The dimensions
    must match the dimensions of yvalue.
    
RETURNS: 
x
    An array with the same dimensions as fvalue and yvalue approximating
    the inversion solution based on interpolation of the data given.
        
xi, yi
    Scalar integer indices of the element where the solution was 
    discovered in the table.

DESCRIPTION:

Similarly to _mapsearch2, _mapsearch2x looks for intersections of the
curves implied by
    f(x, y) = fvalue
    y = yvalue
cross.  Inside of elements, the f(x,y)=fvalue curve is interpolated 
linearly between the points where it crosses along the element edges.

Unlike _mapsearch2, _mapsearch2x does not need to search the entire 
domain for solutions - it only performs operations on the row of 
elements implied by the y-value.  As a result, it is faster.

SEE ALSO:
    _mapsearch1(), _mapsearch2(), _mapsearch2x(), _mapsearch2y()
"""
        # Define an increment for small values
        # For most systems, eps is about 2.2e-16, so small will be about
        # 2.2e-12.  This is the number we use to detect dimensionless
        # proximity to the element boundary.
        small = np.finfo(float).eps * 1e4
        # Initialize result arrays
        x = np.empty_like(fvalue, dtype=float)
        XI = np.empty_like(fvalue, dtype=int)
        YI = np.searchsorted(ydata, yvalue, side='right')-1
        for index in range(fvalue.size):
            yv = yvalue.flat[index]
            fv = fvalue.flat[index]
            yi = YI.flat[index]
            yi1 = yi + 1
            # Compare the values of only the appropriate row
            fI = fv < fdata[:, yi:yi+2]
            # Detect elements with a crossing
            I = crossing2(fI)
            for xi in np.nonzero(I)[0]:
                xi1 = xi+1
                # Initialize some crossing parameters
                fc = []
                #neighbor = []
                # Proceed only if the element is still flagged
                #<<==>>
                #if I[xi,0]:
                # Detect the edges
                # Bottom Edge
                if fI[xi,0] != fI[xi1,0]:
                    xx = interp_scalar(fv, fdata[xi,yi], fdata[xi1,yi], xdata[xi], xdata[xi1])
                    fc.append(np.array([xx, ydata[yi]]))
                    #<<==>>
                    #neighbor.append(None)
                # Left Edge
                if fI[xi,0] != fI[xi,1]:
                    yy = interp_scalar(fv, fdata[xi,yi], fdata[xi,yi1], ydata[yi], ydata[yi1])
                    fc.append(np.array([xdata[xi], yy]))
                    #<<==>>
                    #neighbor.append((xi-1, 0))
                # Top Edge
                if fI[xi,1] != fI[xi1,1]:
                    xx = interp_scalar(fv, fdata[xi,yi1], fdata[xi1,yi1], xdata[xi], xdata[xi1])
                    fc.append(np.array([xx, ydata[yi1]]))
                    #<<==>>
                    #neighbor.append(None)
                # Right Edge
                if fI[xi1,0] != fI[xi1,1]:
                    yy = interp_scalar(fv, fdata[xi1,yi], fdata[xi1,yi1], ydata[yi], ydata[yi1])
                    fc.append(np.array([xdata[xi1], yy]))
                    #<<==>>
                    #neighbor.append((xi1, 0))
                # Detect the saddle case
                if len(fc) != 2:
                    # For now, warn the user, and DO NOT append the case
                    pm.utility.print_warning('mp2._mapsearch2x: Discarded a potential solution near a saddle point.  If you believe this was a legitimate solution, please report the code that generated this warning to the PYroMat GitHub issues page.')
                # Two edges have intersections for each function
                else:
                    fx0 = fc[0]
                    fdx = fc[1] - fc[0]
                    # Detect precise equality at a corner
                    if (fdx == 0).all() and fx0[1] == yv:
                        x.flat[index] = fx0[0]
                        XI.flat[index] = xi
                        break
                        #<<==>>
                        #if neighbor[0] is not None:
                        #    I[*neighbor[0]] = False
                        #if neighbor[1] is not None:
                        #    I[*neighbor[1]] = False
                    else:
                        # Calculate the distance along the f=0 curve to intersect 
                        # Perform the calculations in two steps - leave the division
                        # for last, so we can detect nearly singular problems
                        s = yv - fx0[1]
                        det = fdx[1]
                        
                        if 2*abs(det) > abs(s):
                            s /= det
                            if -small < s < 1+small:
                                x.flat[index] = fx0[0] + fdx[0] * s
                                XI.flat[index] = xi
                                break
                            # Clear the flag for a neighbor if the solution is very near an edge
                            #<<==>>
                            #if -small < s < small and neighbor[0] is not None:
                            #    I[*neighbor[0]] = False
                            #if 1-small < s < 1+small and neighbor[1] is not None:
                            #    I[*neighbor[1]] = False
        return x, XI, YI
        
    def _mapsearch2y(self, xdata, ydata, fdata, xvalue, fvalue):
        r"""Search 2D map for inverse estimates (primative routine)
    y, xi, yi = mapsearch2x(xdata, ydata, fdata, yvalue, fvalue)
    
Uses tabulated data to generate an estimate for y in the 2D inversion
problem
    f(xvalue,y) = fvalue

ARGUMENTS:
xdata, ydata
    One-dimensional array-like containing grid values for the x- and y-
    coordinates.  The sizes of the x- and y-data arrays must match the 
    fdata and gdata arrays (see below).
    
fdata
    Two-dimensional array-like containing tabulated values for f(x,y).  
    The indices should be arranged so that
        fdata[i,j] = f(xdata[i], ydata[j])
        gdata[i,j] = g(xdata[i], ydata[j])
        
xvalue
    The scalar value of x used to interpolate the table.
    
fvalue
    The scalar value of f() for which we are searching.
    
RETURNS: 
y
    One-dimensional array, such that each y value represents a distinct 
    estimated solution.  This implies that, for every entry in 
    the y array,
        f(xvalue,y) =approx= fvalue

xi
    Scalar integer index indicating the row in which xvalue was found.

yi
    One-dimensional array containing indices of the elements in which 
    a solution was found.  

DESCRIPTION:

Similarly to _mapsearch2, _mapsearch2y looks for intersections of the
curves implied by
    f(x, y) = fvalue
    x = xvalue
cross.  Inside of elements, the f(x,y)=fvalue curve is interpolated 
linearly between the points where it crosses along the element edges.

Unlike _mapsearch2, _mapsearch2y does not need to search the entire 
domain for solutions - it only performs operations on the column of 
elements implied by the x-value.  As a result, it is faster.

SEE ALSO:
    _mapsearch1(), _mapsearch2(), _mapsearch2x(), _mapsearch2y()
"""
        # Define an increment for small values
        # For most systems, eps is about 2.2e-16, so small will be about
        # 2.2e-12.  This is the number we use to detect dimensionless
        # proximity to the element boundary.
        small = np.finfo(float).eps * 1e4
        # Initialize result arrays
        y = np.empty_like(fvalue, dtype=float)
        YI = np.empty_like(fvalue, dtype=int)
        XI = np.searchsorted(xdata, xvalue, side='right')
        for index in range(fvalue.size):
            fv = fvalue.flat[index]
            xv = xvalue.flat[index]
            xi1 = XI.flat[index]
            xi = xi1 - 1
            # Compare the values of only the appropriate row
            fI = fv < fdata[xi:xi+2, :]
            # Detect elements with a crossing
            I = crossing2(fI)
            for yi in np.nonzero(I)[1]:
                yi1 = yi+1
                # Initialize some crossing parameters
                fc = []
                #<<==>>
                #neighbor = []
                # Proceed only if the element is still flagged
                if I[0,yi]:
                    # Detect the edges
                    # Bottom Edge
                    if fI[0,yi] != fI[1,yi]:
                        xx = interp_scalar(fv, fdata[xi,yi], fdata[xi1,yi], xdata[xi], xdata[xi1])
                        fc.append(np.array([xx, ydata[yi]]))
                        #<<==>>
                        #neighbor.append((0,yi-1))
                    # Left Edge
                    if fI[0,yi] != fI[0,yi1]:
                        yy = interp_scalar(fv, fdata[xi,yi], fdata[xi,yi1], ydata[yi], ydata[yi1])
                        fc.append(np.array([xdata[xi], yy]))
                        #<<==>>
                        #neighbor.append(None)
                    # Top Edge
                    if fI[0,yi1] != fI[1,yi1]:
                        xx = interp_scalar(fv, fdata[xi,yi1], fdata[xi1,yi1], xdata[xi], xdata[xi1])
                        fc.append(np.array([xx, ydata[yi1]]))
                        #<<==>>
                        #neighbor.append((0,yi1))
                    # Right Edge
                    if fI[1,yi] != fI[1,yi1]:
                        yy = interp_scalar(fv, fdata[xi1,yi], fdata[xi1,yi1], ydata[yi], ydata[yi1])
                        fc.append(np.array([xdata[xi1], yy]))
                        #<<==>>
                        #neighbor.append(None)
                    # Detect the saddle case
                    if len(fc) != 2:
                        # For now, warn the user, and DO NOT append the case
                        pm.utility.print_warning('mp2._mapsearch2y: Discarded a potential solution near a saddle point.  If you believe this was a legitimate solution, please report the code that generated this warning to the PYroMat GitHub issues page.')
                    # Two edges have intersections for each function
                    else:
                        fx0 = fc[0]
                        fdx = fc[1] - fc[0]
                        # Detect precise equality at a corner
                        if (fdx == 0).all() and fx0[0] == xv:
                            y.flat[index] = fx0[1]
                            YI.flat[index] = yi
                            break
                            #if neighbor[0] is not None:
                            #    I[*neighbor[0]] = False
                            #if neighbor[1] is not None:
                            #    I[*neighbor[1]] = False
                        else:
                            # Calculate the distance along the f=0 curve to intersect 
                            # Perform the calculations in two steps - leave the division
                            # for last, so we can detect nearly singular problems
                            s = xv - fx0[0]
                            det = fdx[0]
                            
                            if 2*abs(det) > abs(s):
                                s /= det
                                if -small < s < 1+small:
                                    y.flat[index] = fx0[1] + fdx[1] * s
                                    YI.flat[index] = yi
                                    break
                                # Clear the flag for a neighbor if the solution is very near an edge
                                #<<==>>
                                #if -small < s < small and neighbor[0] is not None:
                                #    I[*neighbor[0]] = False
                                #if 1-small < s < 1+small and neighbor[1] is not None:
                                #    I[*neighbor[1]] = False
        return y, XI, YI



    def _Tsatiter(self, T, p, dL, dV, Ids, Nmax=20, ep=1e-6):
        """Iterates on Maxwell's criteria while holding T constant (primative routine)
    _Tsatiter(T, p, dL, dV, Ids)

T       Saturation temperature used to specify the state.
p       Pressure.  These values are overwritten without being used.
dL      Liquid density.
dV      Vapor density.
Ids     Downselect array.  This is an array of booleans the same size 
        and shape as the property arrays.  Iteration is only performed
        on the corresponding elements set to True.  As states converge,
        the corresponding values are set to False.

Initial guesses for the state are taken from the values in dL, dV, and
T.  The values in p are overwritten.

Optional keywords are:
Nmax        Maximum number of iterations allowed. (def = 20)
ep          Fractional error allowed for convergence (def = 1e-6)
"""
        # Initialize an error vector and a jacobian matrix
        e = np.empty(T.shape + (2,1), dtype=float)
        J = np.empty(T.shape + (2,2), dtype=float)
        fail = True
        for count in range(Nmax):
            # Create down-selected views
            T_ = T[Ids]
            dL_ = dL[Ids]
            dV_ = dV[Ids]
            
            gL,gLt,gLd = self._g(T_, dL_, diff=1)
            gV,gVt,gVd = self._g(T_ ,dV_, diff=1)
            pL,pLt,pLd = self._p(T_, dL_, diff=1)
            p[Ids],pVt,pVd = self._p(T_, dV_, diff=1)
            
            # Error vector
            # The vapor pressure is stored in p
            e[Ids,0,0] = gL - gV
            e[Ids,1,0] = pL - p[Ids]
            # Jacobian
            J[Ids,0,0] = gLd
            J[Ids,0,1] = -gVd
            J[Ids,1,0] = pLd
            J[Ids,1,1] = -pVd
            # Overwrite error with the perturbation to the estimates
            e[Ids,:] = np.linalg.solve(J[Ids,:],e[Ids,:])
            
            # Update unknowns
            dL[Ids] -= e[Ids,0,0]
            dV[Ids] -= e[Ids,1,0]
            
            # Detect convergence
            Ids[Ids] = np.logical_or( np.abs(e[Ids,0,0]) > ep*dL[Ids],
                    np.abs(e[Ids,1,0]) > ep*dV[Ids] )
            
            if not Ids.any():
                fail = False
                break;

                                
        if fail:
            raise pm.utility.PMAnalysisError(f'_Tsatiter: Failed to converge in {Nmax} iterations.')
        

    def _dVsatiter(self, T, p, dL, dV, Ids, Nmax=20, ep=1e-6):
        """Iterates on Maxwell's criteria while holding dV constant (primative routine)
    _dVsatiter(T, p, dL, dV, Ids)

T       Saturation temperature used to specify the state.
p       Pressure.  These values are overwritten without being used.
dL      Liquid density.
dV      Vapor density.
Ids     Downselect array.  This is an array of booleans the same size 
        and shape as the property arrays.  Iteration is only performed
        on the corresponding elements set to True.  As states converge,
        the corresponding values are set to False.

Initial guesses for the state are taken from the values in dL, dV, and
T.  The values in p are overwritten.

Optional keywords are:
Nmax        Maximum number of iterations allowed. (def = 20)
ep          Fractional error allowed for convergence (def = 1e-6)
"""
        # Initialize arrays for the linear algebra
        e = np.empty((Ids.size,) + (2,1), dtype=float)
        J = np.empty((Ids.size,) + (2,2), dtype=float)
        fail = True
        for count in range(Nmax):
            # Create down-selected views
            T_ = T[Ids]
            dL_ = dL[Ids]
            dV_ = dV[Ids]
            
            # Evaluate the properties at the liquid and vapor lines
            gL,gLt,gLd = self._g(T_, dL_,1)
            gV,gVt,gVd = self._g(T_, dV_,1)
            pL,pLt,pLd = self._p(T_ ,dL_ ,1)
            p[Ids],pVt,pVd = self._p(T_ ,dV_ ,1)
            
            # Build the Jacobian on temperature and liquid density
            J[Ids,0,0] = gLt-gVt
            J[Ids,0,1] = gLd
            J[Ids,1,0] = pLt-pVt
            J[Ids,1,1] = pLd
            # Build the error vector
            e[Ids,0,0] = gL-gV        # Gibbs error
            e[Ids,1,0] = pL-p[Ids]    # Pressure error
            # Solve.  Ovewrite error with the estimate perturbation
            e[Ids,:] = np.linalg.solve(J[Ids,:],e[Ids,:])
            # Update temperature and density
            T[Ids] -= e[Ids,0,0]
            dL[Ids] -= e[Ids,1,0]
            
            # Test for convergence
            Ids[Ids] = np.logical_or(np.abs(e[Ids,0,0]) > ep*T_, np.abs(e[Ids,1,0]) > ep*dL_)
            
            # If all points have converged
            if not Ids.any():
                fail = False
                break
        if fail:
            raise pm.utility.PMAnalysisError(f'_dVsatiter: Failed to converge in {Nmax} iterations.')

    def _dLsatiter(self, T, p, dL, dV, Ids, Nmax=20, ep=1e-6):
        """Iterates on Maxwell's criteria while holding dL constant (primative routine)
    _dLsatiter(T, p, dL, dV, Ids)

T       Saturation temperature used to specify the state.
p       Pressure.  These values are overwritten without being used.
dL      Liquid density.
dV      Vapor density.
Ids     Downselect array.  This is an array of booleans the same size 
        and shape as the property arrays.  Iteration is only performed
        on the corresponding elements set to True.  As states converge,
        the corresponding values are set to False.

Initial guesses for the state are taken from the values in dL, dV, and
T.  The values in p are overwritten.

Optional keywords are:
Nmax        Maximum number of iterations allowed. (def = 20)
ep          Fractional error allowed for convergence (def = 1e-6)
"""
        # Initialize arrays for the linear algebra
        e = np.empty(T.shape + (2,1), dtype=float)
        J = np.empty(T.shape + (2,2), dtype=float)
        fail = True
        for count in range(Nmax):
            # Create down-selected views
            T_ = T[Ids]
            dL_ = dL[Ids]
            dV_ = dV[Ids]
            
            # Evaluate the properties at the liquid and vapor lines
            gL,gLt,gLd = self._g(T_, dL_,1)
            gV,gVt,gVd = self._g(T_, dV_,1)
            pL,pLt,pLd = self._p(T_ ,dL_ ,1)
            p[Ids],pVt,pVd = self._p(T_ ,dV_ ,1)
            
            # Build the Jacobian on temperature and liquid density
            J[Ids,0,0] = gLt-gVt
            J[Ids,0,1] = -gVd
            J[Ids,1,0] = pLt-pVt
            J[Ids,1,1] = -pVd
            # Build the error vector
            e[Ids,0,0] = gL-gV        # Gibbs error
            e[Ids,1,0] = pL-p[Ids]    # Pressure error
            # Solve.  Ovewrite error with the estimate perturbation
            e[Ids,:] = np.linalg.solve(J[Ids,:],e[Ids,:])
            # Update temperature and density
            T[Ids] -= e[Ids,0,0]
            dV[Ids] -= e[Ids,1,0]
            
            # Test for convergence
            Ids[Ids] = np.logical_or(np.abs(e[Ids,0,0]) > ep*T_, np.abs(e[Ids,1,0]) > ep*dV_)
            
            # If all points have converged
            if not Ids.any():
                fail = False
                break
        if fail:
            raise pm.utility.PMAnalysisError(f'_dLsatiter: Failed to converge in {Nmax} iterations.')


    def _psatiter(self, T, p, dL, dV, Ids, Nmax=20, ep=1e-6):
        """Iterates on Maxwell's criteria while holding p constant (primative routine)
    _psatiter(T, p, dL, dV, Ids)

T       Saturation temperature.
p       Pressure used to determine the saturation state.
dL      Liquid density.
dV      Vapor density.
Ids     Downselect array.  This is an array of booleans the same size 
        and shape as the property arrays.  Iteration is only performed
        on the corresponding elements set to True.  As states converge,
        the corresponding values are set to False.

T, dL, and dV hold initial guesses for the saturation properties, while
the values in p are treated as a constraint.  Values in p are not 
changed.

Optional keywords are:
Nmax        Maximum number of iterations allowed. (def = 20)
ep          Fractional error allowed for convergence (def = 1e-6)
"""


        # Initialize arrays for the linear algebra
        e = np.empty(T.shape + (3,1), dtype=float)
        J = np.empty(T.shape + (3,3), dtype=float)
        fail = True
        for count in range(Nmax):
            # Generate views of the updated down-selected variables
            T_ = T[Ids]
            dL_ = dL[Ids]
            dV_ = dV[Ids]
            p_ = p[Ids]
            
            gL,gLt,gLd = self._g(T[Ids],dL[Ids],diff=1)
            gV,gVt,gVd = self._g(T[Ids],dV[Ids],diff=1)
            pL,pLt,pLd = self._p(T[Ids],dL[Ids],diff=1)
            pV,pVt,pVd = self._p(T[Ids],dV[Ids],diff=1)
            
            # Error vector
            e[Ids,0,0] = gL - gV
            e[Ids,1,0] = pL - p_
            e[Ids,2,0] = pV - p_
            # Jacobian
            J[Ids,0,0] = gLt-gVt
            J[Ids,0,1] = gLd
            J[Ids,0,2] = -gVd
            
            J[Ids,1,0] = pLt
            J[Ids,1,1] = pLd
            J[Ids,1,2] = 0.
            
            J[Ids,2,0] = pVt
            J[Ids,2,1] = 0.
            J[Ids,2,2] = pVd
            # Overwrite error with the perturbation to the estimates
            e[Ids,:] = np.linalg.solve(J[Ids,:],e[Ids,:])
            
            # Update the variables
            T[Ids] -= e[Ids,0,0]
            dL[Ids] -= e[Ids,1,0]
            dV[Ids] -= e[Ids,2,0]
            
            # Detect convergence
            Ids[Ids] = np.logical_or( np.abs(e[Ids,0,0]) > ep*T[Ids],
                        np.logical_or( np.abs(e[Ids,1,0]) > ep*dL[Ids],
                        np.abs(e[Ids,2,0]) > ep*dV[Ids]))
            
            if not Ids.any():
                fail = False
                break
            
        if fail:
            raise pm.utility.PMAnalysisError(f'_psatiter: Failed to converge in {Nmax} iterations.')

    def _satiter2(self, T, p, dL, dV, x, fn0, fn1, f0value, f1value, Ids, Nmax=10, ep=1e-6):
        """Two-property saturation iteration (primative routine)
    _satiter2(self, T, dL, dV, x, fn0, fn1, f0value, f1value, Ids, Nmax=10, ep=1e-6)

Iteratively calculates the two-phase mixture conditions where a pair of
properties have the prescribed values.  

T           Temperature array
p           Pressure array
dL          Saturated liquid density array
dV          Saturated vapor density array
x           Quality array
fn0         Property method 0
fn1         Property method 1
f0value     Property method 0
f1value     Property method 1 used to calculate x
Ids         Boolean down-select array

This algorithm iteratively solves the problem
    p(T, dL) = p(T, dV)
    g(T, dL) = g(T, dV)
    (1-x) f0(T, dL) + x f0(T, dV) = f0value
    (1-x) f1(T, dL) + x f1(T, dV) = f1value
Because x can be explicitly calculated in each iteration, and because 
some properties (p and g) are constant across the dome, there are 
benefits to eliminating x during iteration, so
    p(T, dL) = p(T, dV)
    g(T, dL) = g(T, dV)
    (f0value - f0L) (f1V - f1L) = (f1value - f1L) (f0V - f0L)

T, dL, and dV contain initial guesses for the saturation conditions.  
p and x are calculated explicitly and will be overwritten.  Quality is
calculated as
    x = (f1value - f1L) / (f1V - f1L)
    
If pressure is one of the properties, it should never be passed as f1,
since pV == pL.  It is MUCH faster to use _psatiter instead.
"""

        E = np.empty(T.shape + (3,1), dtype=float)
        J = np.empty(T.shape + (3,3), dtype=float)

        count = 0
        while Ids.any():
            count += 1
            if count > Nmax:
                raise pm.utility.PMParamError(
                        f'mp2._satiter2: Failed to converge after {Nmax} iterations.')
            
            TT = T[Ids]
            DL = dL[Ids]
            DV = dV[Ids]
            # Evaluate the properties
            pL,pLt,pLd = self._p(TT,DL,diff=1)
            pV,pVt,pVd = self._p(TT,DV,diff=1)
            gL,gLt,gLd = self._g(TT,DL,diff=1)
            gV,gVt,gVd = self._g(TT,DV,diff=1)
            f0L,f0Lt,f0Ld = fn0(TT,DL,diff=1)
            f0V,f0Vt,f0Vd = fn0(TT,DV,diff=1)
            f1L,f1Lt,f1Ld = fn1(TT,DL,diff=1)
            f1V,f1Vt,f1Vd = fn1(TT,DV,diff=1)

            # Property deltas across the dome
            df0 = f0V - f0L
            vf0 = f0value[Ids] - f0L
            df1 = f1V - f1L
            vf1 = f1value[Ids] - f1L

            E[Ids, 0, 0] = pL - pV       # Maxwell, pressure
            E[Ids, 1, 0] = gL - gV       # Maxwell, gibbs energy
            E[Ids, 2, 0] = vf0*df1 - vf1*df0     # Quality constraint
            
            J[Ids, 0, 0] = pVt - pLt
            J[Ids, 0, 1] = -pLd
            J[Ids, 0, 2] = pVd
            
            J[Ids, 1, 0] = gVt - gLt
            J[Ids, 1, 1] = -gLd
            J[Ids, 1, 2] = gVd
            
            J[Ids, 2, 0] = -f1Lt*df0 + vf1*(f0Vt - f0Lt) + f0Lt*df1 - vf0*(f1Vt - f1Lt)
            J[Ids, 2, 1] = -f1Ld*df0 + vf1*f0Ld + f0Ld*df1 + vf0*f1Ld
            J[Ids, 2, 2] = vf1*f0Vd - vf0*f1Vd
            
            delta = np.linalg.solve(J, E)
            T[Ids] += delta[Ids,0,0]
            dL[Ids] += delta[Ids,1,0]
            dV[Ids] += delta[Ids,2,0]
            p[Ids] = pV
            x[Ids] = vf1/df1
            # Update convergence criteria
            Ids[Ids] = (delta[Ids,0,0] > TT*ep) + (delta[Ids,1,0] > DL*ep) + (delta[Ids,2,0] > DV*ep)
            

    def _Titer(self, T, d, fn, fvalue, Ids, Nmax=10, ep=1e-6):
        """Constant-temperature iteration (primative routine)
    _Titer(T, d, fn, fvalue, Ids)

While holding temperature constant, iterates on density to match a 
property with its target values.

Arguments are:
    T       Temperature array
    d       density array
    f       inner property function
    fvalue  array of target property values
    Ids     boolean down-select array
Optional keyword arguments are:
    Nmax    Maximum number of iterations before declaring failure
    ep      epsilon or fractional precision to require in density

The iteration is performed in-place, so the initival values in the T and
d arrays are used as the initial guesses for iteration.  As iteration 
progresses, the Ids boolean array is modified to reflect values that 
have converged.
"""
        count = 0
        while Ids.any():
            count += 1
            # Only permit Nmax iterations
            if count > Nmax:
                raise pm.utility.PMParamError(
                        f'mp2._Titer: Failed to converge after {Nmax} iterations.')
            
            DD = d[Ids]
            TT = T[Ids]
            
            f,ft,fd = fn(T=TT, d=DD, diff=1)
            dd = (fvalue[Ids] - f) / fd
            d[Ids] += dd
            Ids[Ids] = np.abs(dd) > ep * DD


    def _diter(self, T, d, fn, fvalue, Ids, Nmax=10, ep=1e-6):
        """Constant-density iteration (primative routine)
    _diter(T, d, fn, fvalue, Ids)

While holding density constant, iterates on temperature to match a 
property with its target values.

Arguments are:
    T       Temperature array
    d       density array
    fn      inner property function
    fvalue  array of target property values
    Ids     boolean down-select array
Optional keyword arguments are:
    Nmax    Maximum number of iterations before declaring failure
    ep      epsilon or fractional precision to require in density

The iteration is performed in-place, so the initival values in the T and
d arrays are used as the initial guesses for iteration.  As iteration 
progresses, the Ids boolean array is modified to reflect values that 
have converged.
"""
        count = 0
        while Ids.any():
            count += 1
            # Only permit Nmax iterations
            if count > Nmax:
                raise pm.utility.PMParamError(
                        f'mp2._diter: Failed to converge after {Nmax} iterations.')
            
            DD = d[Ids]
            TT = T[Ids]
            
            f,ft,fd = fn(T=TT, d=DD, diff=1)
            dT = (fvalue[Ids] - f) / ft
            T[Ids] += dT
            Ids[Ids] = np.abs(dT) > ep * TT


    def _iter2(self, T, d, f0, f1, f0value, f1value, Ids, Nmax=10, ep=1e-6):
        """Constant-density iteration (primative routine)
    _iter2(T, d, f0, f1, f0value, f1value Ids)

Iterate on both temperature and density to obtain a pair of property 
values.

Arguments are:
    T       Temperature array
    d       density array
    f0      inner property function
    f1      inner property function
    f0value array of target property values
    f1value array of target property values
    Ids     boolean down-select array
Optional keyword arguments are:
    Nmax    Maximum number of iterations before declaring failure
    ep      epsilon or fractional precision to require in density

The iteration is performed in-place, so the initival values in the T and
d arrays are used as the initial guesses for iteration.  As iteration 
progresses, the Ids boolean array is modified to reflect values that 
have converged.
"""
        count = 0
        while Ids.any():
            count += 1
            # Only permit Nmax iterations
            if count > Nmax:
                raise pm.utility.PMParamError(
                        f'mp2._iter2: Failed to converge after {Nmax} iterations.')
            
            DD = d[Ids]
            TT = T[Ids]
            
            # We'll use f and g as placeholder function values
            f,ft,fd = f0(T=TT, d=DD, diff=1)
            g,gt,gd = f1(T=TT, d=DD, diff=1)
            
            # Calculate error arrays
            ef = f0value[Ids] - f
            eg = f1value[Ids] - g
            # Determinants
            det = (ft*gd - fd*gt)
            # Calculate the changes in temperature and density
            dT = (ef*gd-eg*fd)/det
            dd = (-ef*gt+eg*ft)/det
            # Apply the changes
            T[Ids] += dT
            d[Ids] += dd
            # Update the convergence tests
            Ids[Ids] = np.logical_and(np.abs(dT) > ep * TT, np.abs(dd) > ep * DD)


        


    def _fo(self, tt, dd, diff=2):
        """Dimensionless ideal gas helmholtz free energy (primative routine)
Evaluates an ideal gas equation of the form
    a = log(dd) + logt*log(tt) + tlogt*tt*log(tt) + p(t) + ... 
            + c log(1-exp(-theta*tt)) + ...
    
where
    dd = d / dscale
    tt = Tscale / T

In the IGgroup dictionary defined by the mp1 data, the polynomial, p,
is defined by the 'coef0' list.  This list should should be readable
by the _poly1() method.  The 'logt' and 'tlogt' constants define the 
coefficients of the log(tt) and tt*log(tt) terms.  They are optional.

The log/exp expansion is defined by the 'coef1' list.  Each element of
'coef1' should be a two-element list or tuple containing [theta, c]. 
    
This is a PRIMATIVE ROUTINE.  The arguments must already be 
nondimensionalized, and the returned values are non-dimensionalzied.
"""
        
        # Start with the logarithmic terms
        # Log of density is easy - no coefficients needed
        F = np.log(dd)
        Ft = 0.
        Fd = 0.
        Ftt = 0.
        Ftd = 0.
        Fdd = 0.
        
        IGgroup = self.data['IGgroup']
        
        if diff>0:
            Ft = 0
            Fd = 1./dd
            if diff>1:
                Fdd = -Fd/dd
                Ftt = 0.
                Ftd = 0.
        # The logt term does not appear in all models, but it does in many
        logt = None
        coef = IGgroup.get('logt')
        if coef is not None:
            # We might need logt again - keep it for later
            logt = np.log(tt)
            F += coef * logt
            if diff>0:
                pt = coef/tt
                Ft += pt
                if diff>1:
                    Ftt += -pt/tt
        
        # In rare cases, the ideal gas model also includes a t * ln(t) term
        coef = IGgroup.get('tlogt')
        if coef is not None:
            # Don't repeat the logt call if it has already been made
            if logt is None:
                logt = np.log(tt)
            F += coef * tt * logt
            if diff>0:
                Ft += coef * (logt + 1)
                if diff>1:
                    Ftt += coef/tt
        
        
        # Move on to the polynomial expansion
        gr = IGgroup.get('group0')
        if gr is not None:
            p,pt,ptt = self._poly1(tt,gr,diff)
            F+=p
            if diff>0:
                Ft += pt
                if diff>1:
                    Ftt += ptt
        
        # Now the nested log/exp expansion
        gr = IGgroup.get('group1')
        if gr is not None: 
            for theta,c in gr:
                e = np.exp(-theta*tt)
                p = np.log(1-e)
                F += c*p
                if diff>0:
                    pt = theta*e/(1.-e)
                    Ft += c*pt
                    if diff>1:
                        ptt = -pt*(theta + pt)
                        Ftt += c*ptt
                        
        return F, Ft, Fd, Ftt, Ftd, Fdd


    def _fr(self, tt, dd, diff=2):
        """Dimensionless residual helmhotz free energy (primative routine)
Each fit in the group is of the form
    f = exp(-dd**k) * pk(tt, dd)
    
when dd = d / dscale, tt = Tscale / T
    
    F,Fd,Ft,Fdd,Fdt,Ftt = _fr(tt, dd, order=2)

Returns the Helmholtz free energy and its derivatives up to diff.

This is a PRIMATIVE ROUTINE.  The arguments must already be 
nondimensionalized, and the returned values are non-dimensionalzied.
"""
        # Sparse polynomial evaluation is roughly 2x as fast as dense
        # polynomial evaluation for the R134a polynomials.  The 
        # explicitly defined algorithm is also roughly 2x as fast.  the
        # p_d() algorithm for R134a evaluated in about 80us on a 4 core
        # AMD A10-9700B R7

        Rgroup = self.data['Rgroup']

        gr = Rgroup['group0'][0]
        # First evaluate the polynomial without an exponential coefficient
        F,Ft,Fd,Ftt,Ftd,Fdd = self._poly2(tt,dd,gr,diff)
        
        k=0
        ddk = 1.
        for gr in Rgroup['group0'][1:]:
            p,pt,pd,ptt,ptd,pdd = self._poly2(tt, dd, gr, diff)
            
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
                    
                    Ftt += ptt
                    Ftd += ptd
                    Fdd += pdd
                    
                pd = e*(p*ed + pd)
                pt *= e
                
                Ft += pt
                Fd += pd
                
            p *= e
            
            F += p
    
        # Evaluate group1: c * dd**d * tt**t * exp(-a*(dd-ep)**2 - b*(tt-gam)**2)
        gr = Rgroup.get('group1')
        if gr is not None:
            #This is the original table order
            #for c,d,t,a,b,gam,ep in ARgroup['coef1']:
            for t,d,b,a,gam,ep,c in gr:
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
                        ptt = (t-1)*pt/tt
                        pdd = (d-1)*pd/dd
                        ptd = d*pt/dd
                        ett = -2*b*(e + ttm1*et)
                        edd = -2*a*(e + ddm1*ed)
                        etd = -2*b*ttm1*ed
                        
                        Ftt += e*ptt + 2*et*pt + ett*p
                        Ftd += e*ptd + ed*pt + et*pd + etd*p
                        Fdd += e*pdd + 2*ed*pd + edd*p
                    Ft += e*pt + et*p
                    Fd += e*pd + ed*p
                F += e*p
        
        gr = Rgroup.get('group2')
        if gr is not None:
            for a,b,m,AA,BB,CC,DD,c in gr:
                ddm1 = dd-1
                ttm1 = tt-1
                # The model uses 1/m.  We will invert it only once.
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
                    # Detect indices where ddm1 is zero.  If any elements
                    # of ddm1 are small, then we will use a different 
                    # method for calculating derivatives of this term.
                    Ismall = (np.abs(ddm1) < 1e-6).any()
                    
                    # Case out the nearly singular densities
                    # When ddm1 is very small, we can't use the trick of
                    # dividing by it to prevent repeated calls to **
                    # c'est la vie
                    if Ismall:
                        # if m < 1, there is an irreconcilable singularity!
                        if m < 1:
                            pm.utility.print_warning('_ar():: m<1 in the ar2 term. This causes signularities in derivatives near critical density.')
                        pd = AA * m * (ddm1*ddm1)**(0.5*m-0.5)
                        if diff>1:
                            if m < 2:
                                pm.utility.print_warning('_ar():: m<2 in the ar2 term. This causes signularities in second derivatives near critical density.')
                            pdd = AA * m * (m-1) * (ddm1*ddm1)**(0.5*m-1)
                    else:
                        pd = m * p / ddm1
                        if diff>1:
                            pdd = pd * (m-1)/ddm1
                        
                # p = (1-tt) + A(dd-1)**m
                p -= ttm1
                if diff>0:
                    # It's OK to use a scalar here.  We're going to 
                    # multiply it by p in a moment, and that will broad-
                    # cast it appropriately.
                    pt = -1.
                    # Forcing ptt and ptd to zero is not necessary
                    # we already know to ignore them in the next term.
                    #if diff>1:
                        #ptt = np.zeros_like(p)
                        #ptd = np.zeros_like(p)
                
                # Now, square the whole term
                # p = [(1-tt) + A(dd-1)**m]**2
                if diff>0:
                    if diff>1:
                        ptt = 2*pt*pt # + 2*p*ptt (but ptt=0)
                        ptd = 2*pt*pd # + 2*p*ptd (but ptd=0)
                        pdd = 2*pd*pd + 2*p*pdd
                    pt = 2*p*pt
                    pd = 2*p*pd
                p = p*p
                
                # p = [(1-tt) + A(dd-1)**m]**2 + B(dd-1)**2a
                # borrow e for the new term
                e = BB*(ddm1*ddm1)**a
                p += e
                if diff>0:
                    # This term can have a singularity from ddm1 near 0 too
                    if Ismall:
                        # If a is small, then this can be singular!
                        if a < 0.5:
                            pm.utility.print_warning('_ar():: a<0.5 in the ar2 term. This causes singularities in derivatives near critical density.')
                        
                        pd += BB*2*a*(ddm1*ddm1)**(a-0.5)
                        if diff>1:
                            if a < 1.:
                                pm.utility.print_warning('_ar():: a<1 in the ar2 term. This causes singularities in second derivatives near critical density.')
                            pdd += BB*2*a*(2*a-1)*(ddm1*ddm1)**(a-1)
                    else:
                        ed = 2*a*e/ddm1
                        pd += ed
                        if diff>1:
                            edd = (2*a-1)*ed/ddm1
                            pdd += edd

                # e = {[(1-tt) + A(dd-1)**m]**2 + B(dd-1)**2a}**b
                # borrow e for the new term
                # p is now Delta, raise it to the b power.
                e = p**b
                if diff>0:
                    # We're done with Ismall, so repurpose it to detect
                    # small values of delta.
                    Ismall = (np.abs(p)<1e-6).any()
                    if Ismall:
                        if b<1:
                            pm.utility.print_warning('_ar():: b<1 in the ar2 term. This causes singularities in derivatives very close to the critical point.')
                        # Use a temporary to hold this expensive intermediate
                        temp = b*p**(b-1)
                        et = temp*pt
                        ed = temp*pd
                    else:
                        et = b*e*pt/p
                        ed = b*e*pd/p
                        
                    if diff>1:
                        if Ismall:
                            if b<2:
                                pm.utility.print_warning('_ar():: b<2 in the ar2 term. This causes singularities in second derivatives very close to the critical point.')
                            # First, borrow temp, which is currently holding b*p**(b-1)
                            ett = temp*ptt
                            etd = temp*ptd 
                            edd = temp*pdd
                            
                            # Now, calculate the second terms
                            temp = b*(b-1)*p**(b-2)
                            ett += temp*pt*pt
                            etd += temp*pt*pd
                            edd += temp*pd*pd 
                        else:
                            ett = ((b-1)*et*pt + b*e*ptt)/p
                            etd = ((b-1)*et*pd + b*e*ptt)/p
                            edd = ((b-1)*ed*pd + b*e*pdd)/p
                
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
                # Finally done; add the result to F
                F += p*e
                if diff>0:
                    Ft += pt*e + p*et
                    Fd += pd*e + p*ed
                    if diff>1:
                        Ftt += ptt*e + 2*pt*et + p*ett
                        Ftd += ptd*e + pt*ed + pd*et + p*etd
                        Fdd += pdd*e + 2*pd*ed + p*edd
                
        
        return F,Ft,Fd,Ftt,Ftd,Fdd


    def _build_sattab(self, step=0.02, ep=1e-6, verbose=False):
        """Generate saturation table values (primative routine)
    sattab = _buil_sattab(step=0.02, epsilon=1e-6, verbose=False, aslist=False)
    
Constructs a table of values for the temperature, liquid density, vapor
density, and pressure saturation curves.  This works by beginning at the
critical point and following the Maxwell criteria until the triple point
temperature is reached.

Returns a dictionary with numpy array entries:
    T       Saturation temperature array
    p       Saturation pressure array
    g       Gibbs free energy
    dL      Density, liquid
    dV      Density, vapor
    hL      Enthalpy, liquid
    hV      Enthalpy, vapor
    sL      Entropy, liquid
    sV      Entropy, vapor


Optional parameters:
step    The step is the approximate length of the step taken between 
        table entries.  The step size is measured as the magnitude of 
        the vector d [T/Tc, dL/dc, dV/dc].  The temperature and density
        elements are scaled by their critical values.  So, the default
        step size 0.02 corresponds to a step size around 2% of the 
        critical value.
        
ep      Epsilon or fractional error allowed.  Default is 1e-6
        
verbose Print progress to stdout? Default is False

** DESCRIPTION **
_build_sattab() is the first method called to discover the properties
of a mp2 class dataset.  It begins at the critical point and cautiously
explores the saturation curve in intervals defined by the step value.

Once a point on the saturation curve is known, the next (more distant
from the critical point) is approximated by perturbing the vector, 
    x = [T, dL, dV]^T
tangent to the saturation line.  The tangent is determined from the 
Jacobian of the Maxwell criteria, and its magnitude is adjusted to obey
the step magnitude.
    g(T, dL) = g(T, dV)
    p(T, dL) = p(T, dV)

For states near the critical point, the vapor density is held constant
and Newton-Rhapson iteration is used to polish temperature and liquid
density.  This is done because the saturation curve is _very_ sensitive
to small errors in temperature close to the critical point.  

After the saturation curve's tangent line inclines to the point where
dimensionless changes in vapor density are slower than dimensionless 
changes in temperature, the iteration transitions to be constant-density
    ddV / dc < dT / Tc
Far from the critical point, the saturation curves transition to be very
sensitive to small errors in density, but are much more tolerant to 
temperature.

Iterating in series like this produces an algorithm that is quite robust
but relatively slow.  Because it is only called when a substance is
initially imported, it is treated as a tolerable cost.
"""

        Tt = self.data['Tt']
        pc = self.data['pc']
        Tc = self.data['Tc']
        dc = self.data['dc']
        
        # Initialize the outputs
        Ts_array = [Tc]
        dsL_array = [dc]
        dsV_array = [dc]
        ps_array = [pc]
        
        if verbose:
            print('T pc dL dV')
            print('Critical Point:')
            print(f'{Tc:8.2f} {pc:12.4e} {dc:8.2f} {dc:12.4e}')
        
        # Perform the iteration in two steps.  Very close to the critical
        # point, the temperature is nearly constant, so we'll perturb in
        # density steps.  When the vapor density is less than 20% of the
        # critical density, we'll transition to constant temperature.
        
        # We'll need some linear algebra
        J = np.empty((2,2), dtype=float)
        B = np.empty((2,), dtype=float)
        
        # Initialize scalar saturation state
        T = np.array([Tc])
        dL = np.array([dc])
        dV = np.array([dc])
        p = np.array([pc])
        Ids = np.array([1],dtype=bool)
        # Create an initial perturbation of the densities
        # Do not perturb temperature
        dL += step * dc / 1.414
        dV -= step * dc / 1.414
        fail = True
        for count in range(200):
            # Iterate with constant dV
            Ids[0] = True
            self._dVsatiter(T, p, dL, dV, Ids)
            
            Ts_array.insert(0, T[0])
            dsL_array.insert(0, dL[0])
            dsV_array.insert(0, dV[0])
            ps_array.insert(0, p[0])
            
            if verbose:
                print(f'{T[0]:8.2f} {p[0]:12.4e} {dL[0]:8.2f} {dV[0]:12.4e}')
            
            # Perturb the solution to the next interval
            # Assume a unity change in dV, calculate other changes
            ddV = -1.
            # Use the Maxwell criteria and its derivatives to construct
            # a Jacobian and a perturbation vector assuming a unity 
            # change in vapor density.
            gL,gLt,gLd = self._g(T=T,d=dL,diff=1)
            gV,gVt,gVd = self._g(T=T,d=dV,diff=1)
            pL,pLt,pLd = self._p(T=T,d=dL,diff=1)
            pV,pVt,pVd = self._p(T=T,d=dV,diff=1)
            J[0,0] = gLt[0] - gVt[0]
            J[0,1] = gLd[0]
            J[1,0] = pLt[0] - pVt[0]
            J[1,1] = pLd[0]
            B[0] = gVd[0]*ddV
            B[1] = pVd[0]*ddV
            # Solve for the corresponding changes in T and dL
            x = np.linalg.solve(J,B)
            dT = x[0]
            ddL = x[1]
            # Rescale the steps so that the metric T/Tc, d/dc is equal to step
            scale = step / np.sqrt(dT*dT/Tc/Tc + (ddL*ddL + ddV*ddV)/dc/dc)
            dT *= scale
            ddL *= scale
            ddV *= scale
            
            T += dT
            dL += ddL
            dV += ddV
            
            # Detect the exit condition
            # When the fractional change in temperature is larger than 
            # the fractional change in vapor density, transition to
            # constant-temperature iteration.
            if abs(dT/Tc) > abs(ddV/dc):
                fail=False
                break
            
        if fail:
            pm.utility.print_error('This error should never appear in a release - please report this on the PYroMat Github Issues page.')
            raise pm.utility.PMDataError('_build_sattab: Iteration froze near the critical point.' )
        
        if verbose:
            print('Transitioning to constant-temperature:')
        
        fail=True
        for count in range(200):
            # Iterate with constant T
            Ids[0] = True
            self._Tsatiter(T, p, dL, dV, Ids)
            
            Ts_array.insert(0, T[0])
            dsL_array.insert(0, dL[0])
            dsV_array.insert(0, dV[0])
            ps_array.insert(0, p[0])
            
            if verbose:
                print(f'{T[0]:8.2f} {p[0]:12.4e} {dL[0]:8.2f} {dV[0]:12.4e}')
            
            # Perturb the solution to the next interval
            # Assume a unity change in dV, calculate other changes
            dT = -1.
            # Use the Maxwell criteria and its derivatives to construct
            # a Jacobian and a perturbation vector assuming a unity 
            # change in vapor density.
            gL,gLt,gLd = self._g(T=T,d=dL,diff=1)
            gV,gVt,gVd = self._g(T=T,d=dV,diff=1)
            pL,pLt,pLd = self._p(T=T,d=dL,diff=1)
            pV,pVt,pVd = self._p(T=T,d=dV,diff=1)
            J[0,0] = gLd[0]
            J[0,1] = -gVd[0]
            J[1,0] = pLd[0]
            J[1,1] = -pVd[0]
            B[0] = (gVt[0] - gLt[0])*dT
            B[1] = (pVt[0] - pLt[0])*dT
            # Solve for the corresponding changes in T and dL
            x = np.linalg.solve(J,B)
            ddL = x[0]
            ddV = x[1]
            # Rescale the steps so that the metric T/Tc, d/dc is equal to step
            scale = step / np.sqrt(dT*dT/Tc/Tc + (ddL*ddL + ddV*ddV)/dc/dc)
            dT *= scale
            ddL *= scale
            ddV *= scale

            # Detect the exit condition
            # If the next guess would be beyond the triple point, halt
            if T[0] + dT < Tt:
                fail=False
                break

            T += dT
            dL += ddL
            dV += ddV
            
        if fail:
            pm.utility.print_error('This error should never appear in a release - please report this on the PYroMat Github Issues page.')
            raise pm.utility.PMDataError('_build_sattab: Iteration froze near the triple point.' )
        
        scale = (Tt - T[0]) / dT
        ddL *= scale
        ddV *= scale
        
        T[0] = Tt
        dL += ddL
        dV += ddV
        Ids[0] = True
        
        self._Tsatiter(T, p, dL, dV, Ids)
    
        if verbose:
            print('Triple Point:')
            print(f'{T[0]:8.2f} {p[0]:12.4e} {dL[0]:8.2f} {dV[0]:12.4e}')
        
        Ts_array.insert(0, T[0])
        dsL_array.insert(0, dL[0])
        dsV_array.insert(0, dV[0])
        ps_array.insert(0, p[0])
        
        if verbose:
            print(f'Used {len(Ts_array)} points.')
            print('Populating property lists...')
        
        # Convert to Numpy arrays
        Ts_array = np.array(Ts_array)
        ps_array = np.array(ps_array)
        dsL_array = np.array(dsL_array)
        dsV_array = np.array(dsV_array)
        
        self._sattable = {
            'T':Ts_array, 
            'p':ps_array, 
            'dL':dsL_array, 
            'dV':dsV_array,
            'eL':self._e(T=Ts_array,d=dsL_array)[0],
            'eV':self._e(T=Ts_array,d=dsV_array)[0],
            'hL':self._h(T=Ts_array,d=dsL_array)[0],
            'hV':self._h(T=Ts_array,d=dsV_array)[0],
            'sL':self._s(T=Ts_array,d=dsL_array)[0],
            'sV':self._s(T=Ts_array,d=dsV_array)[0],
            'g':self._g(T=Ts_array,d=dsV_array)[0],
            'fL':self._f(T=Ts_array,d=dsL_array)[0],
            'fV':self._f(T=Ts_array,d=dsV_array)[0]
        }
        
        if verbose:
            print('Done')


    def _build_tab(self, NT=100, Nd=100, verbose=False):
        """Generate lookup tables (primative routine)
    _build_tab(NT=100, Nd=100)

Accepts arguments, NT and Nd, which specify a nominal number of 
temperature and density points in the grid.  The actual number may be
significantly more (but not less) in order to obey certain spacing 
conditions (See GRID GENERATION below).

Automatically generates a dictionary attribute ``self._table'' 
containing keyword member arrays:
    'T'     1D temperature array in ascending order
    'd'     1D density array in ascending order
    'p'     2D pressure array
    'e'     2D internal energy array
    'h'     2D enthalpy array
    's'     2D entropy array
    'g'     2D Gibbs energy array
    'f'     2D Helmholts energy array
    'cI'    A 2-tuple of integers indicating the indices of T and d
            containing precisely the critical point.

The properties are evaluated in a rectangular grid using nodes defined
by the 1D T and d arrays.  2D arrays are indexed such that, for a 
property, f(T,d),
    f[i,j] = f(T[i], d[j])
This can be accomplished using Numpy broadcasting rules 
    TT,dd = np.meshgrid(T,d,indexing='ij')
        OR
    TT,dd = np.ix_(T,d)

The critical coordinates, ``cI'' can be used
    Tc == T[cI[0]]
    dc == d[cI[1]]
    pc == p[*cI]
    ec == e[*cI] 
        ... and so on ...

** GRID GENERATION **
Property surfaces have severe curvature with repeated local maxima and 
minima ``under the dome,'' which can pose severe problems to naive 
numerical routines.  To ensure good characterization of the surface, the
T,d grid is chosen with the following rules:

1) Density and temperature arrays must be in ascending order.
2) No step between any two temperatyre or density values may be larger 
    than (Tmax-Tmin)/NT or (dmax-dmin)/Nd respectively.
3) The critical point must be represented precisely as a node.
4) For each temperature below the critical temperature, there must be a
    density value corresponding precisely to the liquid and vapor 
    saturation densities.
5) For each density between the triple point (liquid and vapor) 
    densities, there must be a temperature value corresponding precisely
    to the corresponding saturation temperature.

These rules naturally result in dense grid groupings near the vertical
and horizontal portions of the saturation line, where initial estimates
are most important for numerical convergence.

They also have the useful side effect that segments of the T and d 
arrays precisely represent the saturation curve.  Saturation temperature
density triples can be formed by incrementing from the cI indices.
    Ts = T[cI[0] - k]
    dsV = d[cI[1] - k]
    dsL = d[cI[1] + k]
If the property surfaces were plotted with against indices instead of 
temperature and density values, the saturation points would lie on an 
isosceles triangle descending on either side of the critical point.

** ZERO DENSITY **
There is no minimum density value.  Asymptotically, the Span and Wagner
model converges to ideal gas properties at zero density.  Pressure is
simply zero.  Though energy and enthalpy converge, entropy diverges.  
For the purposes of interpolation to produce a useful initial guess for 
numerical convergence, s, g, and f are set to their respective values at
the triple point vapor density (very low density) instead of zero 
density.
"""

        Tc = self.data['Tc']
        Tmin,Tmax = self.data['Tlim']
        dc = self.data['dc']
        dmin,dmax = self.data['dlim']
        pc = self.data['pc']
        # The nominal temperature step - use to determine density values
        Tstep = (Tmax - Tmin)/NT
        # Generate a nominal density step
        dstep = (dmax - dmin)/Nd

        if verbose:
            print(f'Tmin={Tmin}, Tc={Tc}, Tmax={Tmax}, Tstep={Tstep}')
            print(f'dmin={dmin}, dc={dc}, dmax={dmax}, dstep={dstep}')

        # Generate an array of sub-critical points
        Ts = np.linspace(Tmin, Tc, 1+int(np.ceil((Tc-Tmin)/Tstep)))
        if verbose:
            print(f'Preliminary sub-critical temperature array: {len(Ts)} temperatures.')

        dsV = np.interp(Ts, self._sattable['T'], self._sattable['dV'])
        # Repeatedly bisect the density steps until they are all smaller than dstep
        # Very near the critical point, a high density of temperatures
        # will be needed.
        for count in range(10):
            I = np.nonzero((dsV[1:] - dsV[:-1]) > dstep)[0]+1
            if len(I) == 0:
                break
            if verbose:
                print(f'Refinement step {count}/10: Bisecting {len(I)} intervals.')
            dnew = (dsV[I] + dsV[I-1])*0.5
            Tnew = np.interp(dnew, self._sattable['dV'], self._sattable['T'])
            dsV = np.insert(dsV, I, dnew)
            Ts = np.insert(Ts, I, Tnew)
        
        if verbose:
            print(f'Using {len(Ts)} sub-critical temperatures.')
            print('Constant-temperature polishing far from the critical point...')
        # Generate liquid density values
        dsL = np.interp(Ts, self._sattable['T'], self._sattable['dL'])
        # Create an empty psat array
        ps = np.empty_like(Ts)
        # Polish with constant-temperature far from the critical point
        I = dsV < 0.5 * dc
        Ids = np.array(I)
        self._Tsatiter(Ts, ps, dsL, dsV, Ids)
        # Polish with constant-density near the critical point
        if verbose:
            print('Constant-vapor-density polishing near the critical point...')
        Ids = np.logical_not(I)
        Ids[-1] = False     # Do not polish the critical point
        self._dVsatiter(Ts, ps, dsL, dsV, Ids)


        # Build temperature and density arrays to flesh out the remainder
        # of the parameter space
        T = np.concatenate((
                Ts,
                np.linspace(Ts[-1], Tmax, 1+int(np.ceil((Tmax-Ts[-1])/Tstep)))[1:]
            ))
        # Record the index for the critical point
        Tci = len(Ts)-1
        
        d = np.concatenate((
                [dsV[0]],       # d=0 is the lowest valid density, but will generate singularities  We'll override it later.
                dsV,            # saturated vapor points
                np.flip(dsL[:-1]), # Reverse the liquid densities to be in ascending order and leave off the critical point
                np.linspace(dsL[0], dmax, 1+int(np.ceil((dmax-dsL[0])/dstep)))[1:]    # Use even space throughout the remainder
            ))
        # Record the index for the critical point
        dci = len(dsV)

        if verbose:
            print(f'Appended the rest of the domain.  NT={len(T)}, Nd={len(d)}.')
            print(f'Generating property data...')

        TT,dd = np.broadcast_arrays(*np.ix_(T,d))

        # Generate state data
        p = self._p(T=TT, d=dd)[0]
        e = self._e(T=TT, d=dd)[0]
        h = self._h(T=TT, d=dd)[0]
        s = self._s(T=TT, d=dd)[0]
        g = self._g(T=TT, d=dd)[0]
        f = self._f(T=TT, d=dd)[0]
        
        if verbose:
            print('Evaluating the zero-density limit...')
        # Restore the minimum density to zero
        d[0] = 0.
        # Override the minimum density values
        R = self.data['R']
        Tscale = self.data['IGgroup']['Tscale']
        dscale = self.data['IGgroup']['dscale']
        tt = Tscale / T
        one = np.broadcast_to(np.array(1.), T.shape)
        _,ft,_,_,_,_ = self._fo(tt, one, 1)
        p[:,0] = 0.
        e[:,0] = ft*(R*Tscale)
        h[:,0] = (ft*tt + 1.)*R*T
        # s, g, and f diverge in reality.
        s[:,0] = float('inf')
        g[:,0] = float('-inf')
        f[:,0] = float('-inf')
        
        if verbose:
            print('Interpolating two-phase mixture data...')
        for k in range(1, Tci+1):
            # Indices for the saturation temperature and density
            iT = Tci-k
            iL = dci-k
            iV = dci+k
            
            dmix = d[iL+1:iV]
            xV = ((1./dmix)-(1./d[iL]))/((1./d[iV])-(1./d[iL]))
            xL = 1-xV
            # First, broadcast the constant properties, p and g
            p[iT, iL+1:iV] = p[iT,iV]
            g[iT, iL+1:iV] = g[iT,iV]
            # Next, use quality to calculate the mixture properties
            e[iT, iL+1:iV] = e[iT, iL]*xL + e[iT, iV]*xV
            h[iT, iL+1:iV] = h[iT, iL]*xL + h[iT, iV]*xV
            s[iT, iL+1:iV] = s[iT, iL]*xL + s[iT, iV]*xV
            f[iT, iL+1:iV] = f[iT, iL]*xL + f[iT, iV]*xV
        
        # Build the table dictionary
        self._table = {'T':T, 'd':d, 'cI':(Tci, dci), 'p':p, 'e':e, 'h':h, 's':s, 'g':g, 'f':f}
        if verbose:
            print('Done.')


    def _ds(self, T, diff=0):
        """Calculate saturated liquid and vapor density (inner routine)
    dL,dV,dLT,dVT = _ds(T, diff=0)
    
Iteratively determines the saturation 
"""
        # Create an iteration downselect array
        # and an out-of-bounds array
        I = np.logical_and(T < self.data['Tt'], T > self.data['Tc'])

        dL = np.empty_like(T, dtype=float)
        dV = np.empty_like(T, dtype=float)
        dLt = None
        dVt = None
        if diff:
            dLt = np.empty_like(T, dtype=float)
            dVt = np.empty_like(T, dtype=float)
            dLt[I] = pm.config['def_oob']
            dVt[I] = pm.config['def_oob']
        
        dL[I] = pm.config['def_oob']
        dV[I] = pm.config['def_oob']
        
        I = np.logical_not(I)
        
        Ts = self.sattab['T']
        dsL = self.sattab['dL']
        dsV = self.sattab['dV']
        
        temp = T[I]
        ii = np.searchsorted(Ts, temp)
        temp -= Ts[ii-1]
        temp /= (T[ii] - T[ii-1])
        temp1 = 1 - temp
        dL[I] = temp * dsL[ii-1] + temp1 * dsL[ii]
        dV[I] = temp * dsV[ii-1] + temp1 * dsV[ii]

        p = np.empty_like(T)
        self._Tsatiter(T, p, dL, dV, I)
        return dL, dV
        
        
    def _ps(self,T,diff=0):
        """Saturation pressure (inner routine)
    ps, ps_T, ps_TT = _ps(T, diff=0)
    
Presumes temperature is in Kelvin, reports pressure in Pa
"""
        group = self.data['PSgroup']
        Tscale = group['Tscale']
        pscale = group['pscale']
        
        p,pt,ptt = self._satfit( 
                T/Tscale,
                group['fn'],
                group['poly'],
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
        T,Tmin,Tmax = np.broadcast_arrays(T, self.data['Tt']*.99, self.data['Tc'])
        
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

        Tscale = self.data['Rgroup']['Tscale']
        dscale = self.data['Rgroup']['dscale']
        R = self.data['R']
        # Calculate dimensionless arrays
        tt = Tscale/T
        dd = d/dscale
        # Calculate the Helmholtz free energy
        _,_,ard,_,artd,ardd = self._fr(tt,dd,diff+1)
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
            # Move the saturation bounds by 1%
            db[Istate] *= 1.01
            # Now, isolate the liquid points; set the lower density to the
            # saturated liquid density
            Istate[Itest] = np.logical_not(Istate[Itest])
            da[Istate] = self._dsl(T[Istate], 0)[0]
            db[Istate] = self.data['dlim'][1]
            # Reduce the lower density by 1%
            da[Istate] *= 0.99
        
        # Iteratively reduce da until all points are bracketed
        Itest = self._p(T,da,0)[0] > p
        while Itest.any():
            da[Itest]/=2.
            Itest[Itest] = self._p(T[Itest], da[Itest],0)[0] > p[Itest]
        
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
                Nmax=50,
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
        Itest = np.asarray(p>=self.data['pc'], dtype=bool)
        Ta[Itest] = self.data['Tlim'][0]
        Tb[Itest] = self.data['Tlim'][1]
        
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
            # Grow the boundary by 1%
            Tb[Isat] *= 1.01
            
            # Now, identify the vapor points
            Isat[Itest] = d[Itest] < dsV[Itest]
            # Leave Ta as the saturation temperature
            Tb[Isat] = self.data['Tlim'][1]
            # Grow the boundary by 1%
            Ta[Isat] = np.maximum(0.99*Ta[Isat], self.data['Tlim'][0])
            
            # Now, get the saturated states
            Isat[Itest] = np.logical_and(
                    d[Itest] >= dsV[Itest],
                    d[Itest] <= dsL[Itest])
            # We now have the solution at these points.
            # Assign the value to T
            T[Isat] = Ta[Isat]
            # Put safe values in Ta and Tb... just in case
            Tb[Isat] = self.data['Tlim'][1]
            Ta[Isat] = self.data['Tlim'][0]
            # Eliminate these from the down-select array - no iteraiton required.
            I[Isat] = False
        
        # Note from v2.2.0... It is necessary to use _tditer instead of
        # using _p directly. Even when p is super-critical, when d is 
        # under the dome, the lower temeprature guess reverts to a sub-
        # critical state, and the _p() values diverge wildly there.  The
        # ideal future fix would be to invert the dsL or dsV lines to 
        # find the actual minimum T at the specified density, but for 
        # v2.2.1, we will revert to _tditer().
        self._hybrid1(
                self._tditer,
                'T',
                p,
                T,
                I,
                Ta,
                Tb,
                param={'d':d, 'fn':self._p})
        
        if sat:
            return T, dsL, dsV, Isat
        return T
        
        
    def _sat_argparse(self, T=None, p=None, Nmax=20, ep=1e-6):
        """A standard argument parsing scheme for all user-layer saturation properties
    T,p,dL,dV = _sat_argparse(T=None, p=None)
    
Enforces that all returned parameters are numpy arrays with at least one
dimension.  Accepts T and p as scalars or array-like objects in 
[unit_temperature] and [unit_pressure] respectively.
    
Returns
T   the temperature in K
dL and dV are the liquid and vapor densities in kg/m3
"""
        Ts = self.data['sattab']['T']
        ps = self.data['sattab']['p']
        dsL = self.data['sattab']['dL']
        dsV = self.data['sattab']['dV']
        if p is None:
            if T is None:
                T = pm.config.def_T()
            T = pm.units.temperature_scale(
                    np.asarray(T, dtype=float), 
                    to_units='K')
            if T.ndim==0:
                T = np.reshape(T, (1,))
            
            # Initialize results
            p = np.full_like(T, pm.config['def_oob'])
            dL = np.full_like(T, pm.config['def_oob'])
            dV = np.full_like(T, pm.config['def_oob'])
            
            # Detect points that are precisely equal to the critical point
            Ids = (T == self.data['Tc'])
            p[Ids] = self.data['pc']
            dL[Ids] = self.data['dc']
            dV[Ids] = self.data['dc']
            
            # Detect points that are in-bounds            
            Ids = np.logical_and(T >= self.data['Tt'], T < self.data['Tc'])
            # Happy calculating!
            self._Tsatiter(T,p,dL,dV,Ids)
                
        elif T is None:
            p = pm.units.pressure(
                    np.asarray(p, dtype=float), 
                    to_units='Pa')
            if p.ndim==0:
                p = np.reshape(p, (1,))

            # Initialize results
            T = np.full_like(p, pm.config['def_oob'])
            dL = np.full_like(p, pm.config['def_oob'])
            dV = np.full_like(p, pm.config['def_oob'])
            
            # Detect points that are precisely equal to the critical point
            Ids = (p == self.data['pc'])
            T[Ids] = self.data['Tc']
            dL[Ids] = self.data['dc']
            dV[Ids] = self.data['dc']
            
            # Detect points that are in-bounds
            Ids = np.logical_and(p >= self.data['pt'], p < self.data['pc'])
            # Happy calculating!
            self._psatiter(T,p,dL,dV,Ids)

        else:
            raise pm.utility.PMParamError(
                '_sat_argparse: Saturation temperature and pressure cannot be simultaneously specified')

        return T, p, dL, dV
        
        
    def _argparse(self, *varg, **kwarg):
        """Present a standard argument scheme for all user-layer property methods
    T,d1,d2,x,I = _argparse( .. keyword arguments ..)

Accepts keyword arguments:
    e   internal energy
    f   free energy
    g   Gibbs energy
    h   enthalpy
    s   entropy
    T   temperature
    p   pressure
    d   density
    v   specific volume
    x   quality

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
        
        # 1) Handle varg and kward and their defaults
        # 2) Apply the argument rules...
        #   2.1: All arguments must be legal
        #   2.2: There are only two arguments unless one is x
        #   2.3: x may only be specified with T, g, or p
        #   2.4: Energy properties, T, e, h, f, and g may not be specified together
        #   2.5: d and v may not be specified together 
        #   
        # 3) Convert the arguments to arrays with dim 1 or greater
        # 4) Convert to standard units
        # 5) Check for out-of-bounds on basic arguments
        # 6) Replace specific volume with density if it appears
        # 7) Case out the possible combinations
        #   7.1: x is specified
        #       7.1.1: x,T,p
        #       7.1.2: x,T
        #       7.1.3: x,p
        #   7.2: Two inverse properties
        #   7.3: One inverse property
        #       7.3.1: T,?
        #       7.3.2: d,?
        #   7.4: No inverse properties
        #       7.4.1: T,d
        #       7.4.2: Unhandled Exception
        # 8) Broadcast the arrays appropriately
        # 9) Calculate T,d1,d2,x, and I

        # 1) Handle varg and kwarg and apply defaults

        # If varg is specified, assign its values to T,p,x
        Nargs = len(varg)
        if Nargs > 0:
            if 'T' in kwarg:
                raise pm.utility.PMParamError('T was specified both positionally and with a keyword.')
            kwarg['T'] = varg[0]
            if Nargs > 1:
                if 'p' in kwarg:
                    raise pm.utility.PMParamError('p was specified both positionally and with a keyword.')
                kwarg['p'] = varg[1]
                if Nargs > 2:
                    if 'x' in kwarg:
                        raise pm.utility.PMParamError('x was specified both positionally and with a keyword.')
                    kwarg['x'] = varg[2]
                    if Nargs > 3 :
                        raise pm.utility.PMParamError('Property calls with more than two arguments require keywords.')

        # Re-count the number of arguments -- in kwarg this time
        # We'll assign default properties to T,p to ensure there are 
        # at least two arguments.
        Nargs = len(kwarg)
        if Nargs == 1:
            if 'T' not in kwarg:
                kwarg['T'] = pm.config.def_T()
            else:
                kwarg['p'] = pm.config.def_p()
        elif Nargs == 0:
            kwarg['T'] = pm.config.def_T()
            kwarg['p'] = pm.config.def_p()
        
        # 2) Apply the argument rules
        # Re-measure the number of arguments and use sets to enforce
        # the remaining rules
        Nargs = len(kwarg)
        args = set(kwarg.keys())
        # inverse_methods is a map between the property names that require
        # iteration and the inner method that calculates it.  Inverse 
        # args is a set of their names that will be used for argument 
        # parsing
        inverse_methods = {'p':self._p, 'e':self._e, 'h':self._h, 's':self._s, 'f':self._f, 'g':self._g}
        inverse_args = set(inverse_methods.keys())
        # basic_args are the remaining legal arguments that do not need
        # iteration (OK, p does, but it's special). 
        # legal_args are all arguments that can be legally accepted.
        basic_args = set(['T','d','v','x'])
        legal_args = inverse_args.union(basic_args)
        # Find the number of inverse arguments
        inverse_args &= args
        Ninv = len(inverse_args)
        
        
        # 2.1: All arguments must be "legal" recognized arguments
        these_args = args - legal_args
        if these_args:
            message = 'Unrecognized propert(y/ies):'
            prefix = '  '
            for name in these_args:
                message += prefix + name
                prefix = ', '
            raise pm.utility.PMParamError(message)
        
        # Special rules applying to quality
        if 'x' in args:
            # 2.2: There may only be 2 arguments excluding 'x'
            if Nargs > 3:
                raise pm.utility.PMParamError(
                    'Specifying more than two simultaneous parameters is illegal (except x with T and p).')
            # 2.3: 'x' may only be specified with T or p
            if args - {'x', 'T', 'p'}:
                raise pm.utility.PMParamError(
                        'Quality may only be specified with temperature and/or pressure.')
        # 2.2: There may only be 2 arguments excluding 'x'
        elif Nargs > 2:
            raise pm.utility.PMParamError(
                    'Specifying more than two simultaneous parameters is illegal (except x with T,p or g,p).')
        # 2.4: T, e, h, f, and g may not be specified together
        if len(args.intersection({'T', 'e', 'h', 'f', 'g'})) > 1:
            raise pm.utility.PMParamError(
                    'Energy parameters, T, e, h, f, or g, may not be specified as a pair.')
        # 2.5: Density and specific volume cannot be specified together
        if 'v' in args and 'd' in args:
            raise pm.utility.PMParamError('Density (d) and specific volume (v) cannot be specified together.')

        
        # 3) Convert all arguments to numpy arrays
        #    The asarray function does NOT copy the array if it is already
        #    a numpy array.
        for name,value in kwarg.items():
            value = np.asarray(value, dtype=float)
            if value.ndim == 0:
                value = np.reshape(value, (1,))
            kwarg[name] = value
        
        # 4) Convert the units appropriately
        #   This step will only make a copy of the array if the units need
        #   to be converted.  Otherwise, the array is passed through verbatim
        #   As a result, the input array will ONLY be copied if it needs to
        #   be reshaped, converted, or retyped.
        # 5) Check for out-of-bounds on the converted values
        #   Checking before arrays are broadcast minimizes the number of
        #   elements that need to be inspected
        # 6) Replace v with d if it appears
        if 'T' in kwarg:
            kwarg['T'] = pm.units.temperature_scale(kwarg['T'], to_units='K')
            # Test for out-of-bounds
            Ioob = np.logical_or(kwarg['T'] < self.data['Tlim'][0], 
                    kwarg['T'] > self.data['Tlim'][1])
            if Ioob.all():
                pm.utility.print_warning('All of the temperature values are out-of-bounds for this substance.'
                        'Legal values are between {:f} and {:f} Kelvin.'.format(*self.data['Tlim']))
                raise pm.utility.PMParamError('_ARGPARSE: Temperature values were all out of range.')
            elif Ioob.any():
                kwarg['T'][Ioob] = pm.config['def_oob']
                pm.utility.print_warning('Some temperature values were out-of-bounds for this substance.')
        if 'p' in kwarg:
            kwarg['p'] = pm.units.pressure(kwarg['p'], to_units='Pa')
            # Test for out-of-bounds
            Ioob = np.logical_or(kwarg['p'] < self.data['plim'][0], 
                    kwarg['p'] > self.data['plim'][1])
            if Ioob.all():
                pm.utility.print_warning('All of the pressure values are out-of-bounds for this substance.'
                        'Legal values are between {:f} and {:f} Pascals.'.format(*self.data['plim']))
                raise pm.utility.PMParamError('_ARGPARSE: Pressure values were all out of range.')
            elif Ioob.any():
                kwarg['p'][Ioob] = pm.config['def_oob']
                pm.utility.print_warning('Some pressure values were out-of-bounds for this substance.')
        if 'd' in kwarg:
            value = pm.units.volume(kwarg['d'], to_units='m3', exponent=-1)
            kwarg['d'] = pm.units.matter(value, self.data['mw'], to_units='kg')
        if 'v' in kwarg:
            # Convert and replace with d at the same time
            value = pm.units.volume(kwarg['v'], to_units='m3')
            kwarg['d'] = 1./pm.units.matter(value, self.data['mw'], to_units='kg', exponent=-1)
            # Update the keywords and argument sets to reflect the
            # substitution.
            args.add('d')
            basic_args.add('d')
            # Keep v - it is sometimes useful
        for this in ['h', 'e', 'f', 'g']:
            if this in kwarg:
                value = kwarg[this]
                value = pm.units.energy(value, to_units='J')
                value = pm.units.matter(value, self.data['mw'], to_units='kg', exponent=-1)
                kwarg[this] = value
        if 's' in kwarg:
            value = kwarg['s']
            value = pm.units.energy(value, to_units='J')
            value = pm.units.matter(value, self.data['mw'], to_units='kg', exponent=-1)
            value = pm.units.temperature(value, to_units='K', exponent=-1)
            kwarg['s'] = value
        # x is dimensionless - no need to convert anything
        if 'x' in kwarg:
            if (kwarg['x'] > 1).any() or (kwarg['x'] < -1).any():
                raise pm.utility.PMParamError('Quality was found to be outside of the range -1,1.')

        
        # 7: Case out the different property combinations
        # 7.1: x is specified
        if 'x' in kwarg:
            if 'T' in kwarg:
                # 7.1.1: T,p,x
                # This is the special case that lets T,p specify any state
                if 'p' in kwarg:
                    T,p,x = np.broadcast_arrays(kwarg['T'], kwarg['p'], kwarg['x'])
                    I = (x >= 0)
                    d1 = np.empty_like(T)
                    d2 = np.empty_like(T)
                    # Calculate densities for saturated states
                    if I.any():
                        d1[I] = np.interp(T[I], self._sattable['T'], self._sattable['dL'], left=config['def_oob'], right=config['def_oob'])
                        d2[I] = np.interp(T[I], self._sattable['T'], self._sattable['dV'], left=config['def_oob'], right=config['def_oob'])
                        Ids = I.copy()
                        self._Tsatiter(T,p,d1,d2,Ids)
                    # Calculate densities for non-saturated states
                    Ids = np.logical_not(I)
                    if Ids.any():
                        # Calculate densities for non-saturated points
                        d2[Ids],Ti,di = self._mapsearch2y(self._table['T'], self._table['d'], self._table['p'], T[Ids], p[Ids])
                        self._Titer(T, d2, self._p, p, Ids.copy())
                        d1[Ids] = d2[Ids]
                    return T, d1, d2, x, I
                # 7.1.2: T,x
                else:
                    T,x = np.broadcast_arrays(kwarg['T'], kwarg['x'])
                    p = np.empty_like(T)                    
                    I = (x >= 0)
                    if not I.all():
                        raise pm.utility.PMParamError(
                            'Found x < 0.  Only two-phase mixtures can be specified with T,x.  All values of x must be [0,1].')
                    dL = np.interp(T, self._sattable['T'], self._sattable['dL'], left=pm.config['def_oob'], right=pm.config['def_oob'])
                    dV = np.interp(T, self._sattable['T'], self._sattable['dV'], left=pm.config['def_oob'], right=pm.config['def_oob'])
                    self._Tsatiter(T,p,dL,dV,I.copy())
                    return T, dL, dV, x, I
            # 7.1.3: p,x
            else:
                p,x = np.broadcast_arrays(kwarg['p'], kwarg['x'])
                I = (x >= 0)
                if not I.all():
                    raise pm.utility.PMParamError(
                        'Found x < 0.  Only two-phase mixtures can be specified with p,x.  All values of x must be [0,1].')
                dL = np.interp(p, self._sattable['p'], self._sattable['dL'], left=pm.config['def_oob'], right=pm.config['def_oob'])
                dV = np.interp(p, self._sattable['p'], self._sattable['dV'], left=pm.config['def_oob'], right=pm.config['def_oob'])
                T = np.interp(p, self._sattable['p'], self._sattable['T'], left=pm.config['def_oob'], right=pm.config['def_oob'])
                self._psatiter(T,p,dL,dV,I.copy())
                return T, dL, dV, x, I
            
        # 7.2: Two inverse properties
        elif Ninv > 1:
            # Isolate the property strings, their methods, and their value arrays
            f0str = args.pop()
            f1str = args.pop()
            fn0 = inverse_methods[f0str]
            fn1 = inverse_methods[f1str]
            f0value, f1value = np.broadcast_arrays(kwarg[f0str], kwarg[f1str])
            # Look up estimates for T and d in the property tables
            T,d2,Ti,di = self._mapsearch2(self._table['T'], self._table['d'], self._table[f0str], self._table[f1str], f0value, f1value)
            x = np.full_like(T, -1.)
            d1 = np.empty_like(d2)
            # Test for entries under the dome
            k = self._table['cI'][0] - Ti
            diL = self._table['cI'][1] + k
            diV = self._table['cI'][1] - k
            I = (k>0) * (diV <= di) * (di < diL)
            if I.any():
                # g,p iteration will fail under the dome
                if args == {'g', 'p'}:
                    raise pm.utility.PMParamError(
                            'mp2._argparse: Received g and p in or very close to a two-phase mixture: numerically singular.')
                # Obtain estimates for saturation densities
                d1[I] = np.interp(T[I], self._sattable['T'], self._sattable['dL'])
                d2[I] = np.interp(T[I], self._sattable['T'], self._sattable['dV'])
                # Constant-pressure iteration under the dome is a special case
                # Pressure gives us temperature and densities explicitly,
                # Then x can be calculated from f1value
                if f0str == 'p':
                    self._psatiter(T, f0value, d1, d2, I.copy())
                    f1L,_,_ = fn1(T[I], d1[I], diff=0)
                    f1V,_,_ = fn1(T[I], d2[I], diff=0)
                    x[I] = (f1value[I] - f1L)/(f1V - f1L)
                elif f1str == 'p':
                    self._psatiter(T, f1value, d1, d2, I.copy())
                    f0L,_,_ = fn0(T[I], d1[I], diff=0)
                    f0V,_,_ = fn0(T[I], d2[I], diff=0)
                    x[I] = (f0value[I] - f0L)/(f0V - f0L)
                # For other property combinations, it will be necessary to 
                # iterate.
                else:
                    self._satiter2(T, np.empty_like(T), d1, d2, x, fn0, fn1, f0value, f1value, I.copy())
            # Detect states that are not quite under the dome, but very
            # close.  These will have converged to out-of-bounds values
            # for x.
            Ids = np.zeros_like(I, dtype=bool)
            Ids[I] = (x[I] < 0)
            if Ids.any():
                d2[Ids] = d1[Ids]
                I[Ids] = False
                x[Ids] = -1
            Ids[I] = (x[I] > 1)
            if Ids.any():
                d1[Ids] = d2[Ids]
                I[Ids] = False
                x[Ids] = -1
            
            # All other states
            Ids = np.logical_not(I)
            if Ids.any():
                self._iter2(T, d2, fn0, fn1, f0value, f1value, Ids.copy())
                d1[Ids] = d2[Ids]
            return T,d1,d2,x,I
        # 7.3: One inverse property
        elif Ninv > 0:
            # 7.3.1: T,?
            if 'T' in kwarg:
                pass
            # 7.3.2: d,?
            elif 'd' in kwarg:
                pass
            # UNHANDLED CASE
            else:
                pass
        # 7.4: T,d
        else:
            pass
        
        message = 'Please report a bug: Unhandled event [MASTER] in mp2._argparse with args:'
        prefix = ' '
        for name in args:
            message += prefix + name
            prefix = ', '
        raise pm.utility.PMParamError(message)



    def _e(self,T,d,diff=0):
        """Internal energy (inner routine)
    e,eT,ed = _e(T,d,diff=0)
"""
        eT = None
        ed = None

        # The IG part        
        R = self.data['R']
        Tscale = self.data['IGgroup']['Tscale']
        dscale = self.data['IGgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,at,_,att,atd,_ = self._fo(tt,dd,diff+1)
        
        e = at
        if diff>0:
            eT = tt*tt*att
            ed = atd/dscale
        
        # The residual part
        Tscale = self.data['Rgroup']['Tscale']
        dscale = self.data['Rgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,at,ad,att,atd,add = self._fr(tt,dd,diff+1)
        e += at
        e *= R*Tscale
        if diff>0:
            eT += tt*tt*att
            eT *= -R
            ed += atd/dscale
            ed *= R*Tscale

        return e,eT,ed


    def _h(self,T,d,diff=0):
        """enthalpy (inner routine)
    h,hT,hd = _h(T,d,diff=0)
"""
        hT = None
        hd = None

        # The IG part        
        R = self.data['R']
        Tscale = self.data['IGgroup']['Tscale']
        dscale = self.data['IGgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,at,_,att,atd,_ = self._fo(tt,dd,diff+1)
        
        h = 1. + tt*at
        if diff>0:
            hT = 1. - tt*tt*att
            hd = tt*atd/dscale
        
        # The residual part
        Tscale = self.data['Rgroup']['Tscale']
        dscale = self.data['Rgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,at,ad,att,atd,add = self._fr(tt,dd,diff+1)
        h += dd*ad + tt*at
        h *= R*T
        if diff>0:
            hT += dd*ad - tt*(tt*att + dd*atd)
            hT *= R
            hd += (ad + dd*add + tt*atd)/dscale
            hd *= R*T

        return h,hT,hd

    def _s(self,T,d,diff=0):
        """entropy (inner routine)
    s,sT,sd = _s(T,d,diff=0)
"""
        sT = None
        sd = None

        # The IG part        
        R = self.data['R']
        Tscale = self.data['IGgroup']['Tscale']
        dscale = self.data['IGgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        a,at,ad,att,atd,_ = self._fo(tt,dd,diff+1)
        
        s = tt*at - a
        if diff>0:
            sT = tt*tt*att
            sd = (tt*atd - ad)/dscale
        
        # The residual part
        Tscale = self.data['Rgroup']['Tscale']
        dscale = self.data['Rgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        a,at,ad,att,atd,_ = self._fr(tt,dd,diff+1)
        s += tt*at - a
        s *= R
        if diff>0:
            sT += tt*tt*att
            sT *= -R/T
            sd += (tt*atd - ad)/dscale
            sd *= R

        return s,sT,sd

    def _f(self, T, d, diff=0):
        """Free energy
    f,ft,fd = _f(T,d,diff=0)
    
"""
        # The IG part        
        R = self.data['R']
        Tscale = self.data['IGgroup']['Tscale']
        dscale = self.data['IGgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        a,at,ad,_,_,_ = self._fo(tt,dd,diff)

        f = a
        ft = None
        fd = None
        if diff:
            ft = a - tt*at
            fd = ad/dscale
        
        # The residual part
        Tscale = self.data['Rgroup']['Tscale']
        dscale = self.data['Rgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        a,at,ad,_,_,_ = self._fr(tt,dd,diff)

        f += a
        f *= R*T
        if diff:
            ft += a - tt*at
            ft *= R
            fd += ad/dscale
            fd *= R*T
            
        return f,ft,fd
        


    def _g(self, T, d, diff=0):
        """Gibbs energy
    g,gt,gd = _g(T,d,diff=0)
    
"""
        # The IG part        
        R = self.data['R']
        Tscale = self.data['IGgroup']['Tscale']
        dscale = self.data['IGgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        a,at,ad,_,atd,add = self._fo(tt,dd,diff+1)
        
        g = a + 1.
        gt = None
        gd = None
        if diff:
            gt = a + 1. - tt*at
            gd = ad/dscale
        
        # The residual part
        Tscale = self.data['Rgroup']['Tscale']
        dscale = self.data['Rgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        a,at,ad,_,atd,add = self._fr(tt,dd,diff+1)

        g += a + dd*ad
        g *= R*T
        if diff:
            gt += a + dd*ad - tt*(at + dd*atd)
            gt *= R
            gd += (2*ad + dd*add)/dscale
            gd *= R*T
            
        return g,gt,gd
        
    def _a(self,T,d):
        """Speed of sound (inner routine)
    a = _a(T,d)
"""
        R = self.data['R']
        Tscale = self.data['IGgroup']['Tscale']
        dscale = self.data['IGgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,_,att,_,_ = self._fo(tt,dd,2)

        # We'll build this in three terms
        # b - c*c/d
        B = 1       # The IG portion of b and c are simple
        C = 1
        D = tt * tt * att
        
        # The residual part
        Tscale = self.data['Rgroup']['Tscale']
        dscale = self.data['Rgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,ad,att,atd,add = self._fr(tt,dd,2)
        B += dd*(2*ad + dd*add)
        C += dd*(ad - tt*atd)
        D += tt * tt * att

        return np.sqrt(R * T * (B - C*C/D))

        
    def _cp(self,T,d):
        """Isobaric specific heat (inner routine)
    cp = _cp(T,d)
"""

        # The IG part        
        R = self.data['R']
        Tscale = self.data['IGgroup']['Tscale']
        dscale = self.data['IGgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,_,att,_,_ = self._fo(tt,dd,2)
        
        cp = -tt*tt*att
        
        # The residual part
        Tscale = self.data['Rgroup']['Tscale']
        dscale = self.data['Rgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,ad,att,atd,add = self._fr(tt,dd,2)

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
        Tscale = self.data['IGgroup']['Tscale']
        dscale = self.data['IGgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,_,att,_,_ = self._fo(tt,dd,2)
        
        cv = tt*tt*att
        
        # The residual part
        Tscale = self.data['Rgroup']['Tscale']
        dscale = self.data['Rgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,_,att,_,_ = self._fr(tt,dd,2)

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

Calls to ps() are MUCH faster than calls to Ts(), so when given a choice,
specifying saturation states with temperature should always be preferred.
The MP1 class exposes ps() as an empirical relationship, while Ts() has 
to perform iterative numerical inversion.

Unlike the other saturation properties, ps() and Ts() only accept one
argument and only return one value - each calculates the one in terms
of the other.
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
        
        
    def Ts(self, p=None):
        """Saturation temperature
    Tsat = Ts(p)
    
Calculates the saturation temperature in terms of the pressure.  

Unlike the other saturation properties, ps() and Ts() only accept one
argument and only return one value - each calculates the one in terms
of the other.
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
        
    def vs(self, *varg, **kwarg):
        """Saturation specific volume
    vsL, vsV = vs(T)
    
If no keyword is specified, saturation properties interpret the argument
as temperature.  However, pressure can be specified as well

    vsL, vsV = vs(p=pvalue)
    
Returns the liquid (vsL) and vapor (vsV) saturation density in units
[unit_volume / unit_matter]
"""
        dL,dV = self.ds(*varg, **kwarg)
        return 1./dL, 1./dV


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
    
    def p(self, *varg, quality=False, **kwarg):
        """Pressure
    p(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Additionally, if the optional keyword, "quality" is set to True, the 
quality of the liquid/vapor mixture is also returned
    e,x = e(..., quality=True)

Returns pressure in unit_pressure
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        # Use d2.  In theory, p(d1) = p(d2), but the liquid is so stiff
        # that small numerical errors cause huge pressure errors
        # The problem is solved when the vapor density is used instead.
        # In all other conditions d1=d2
        p = self._p(T,d2,0)[0]
        
        p = pm.units.pressure(p, from_units='Pa')
        
        if quality:
            return p,x
        return p
        
        
    def d(self, *varg, quality=False, **kwarg):
        """Density
    d(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Additionally, if the optional keyword, "quality" is set to True, the 
quality of the liquid/vapor mixture is also returned
    e,x = e(..., quality=True)

Returns density in unit_matter / unit_volume
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        if I.any():
            d1[I] = (1.-x[I])/d1[I]
            d1[I] += x[I]/d2[I]
            d1[I] = 1. / d1[I]
            
        d1 = pm.units.matter(d1, self.data['mw'], from_units='kg')
        d1 = pm.units.volume(d1, from_units='m3', exponent=-1)
        if quality:
            return d1,x
        return d1
        
        
    def v(self, *varg, quality=False, **kwarg):
        """specific volume
    v(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Additionally, if the optional keyword, "quality" is set to True, the 
quality of the liquid/vapor mixture is also returned
    v,x = v(..., quality=True)

Returns volume in unit_volume / unit_matter
"""
        d,x = self.d(*varg, quality=True, **kwarg)
        if quality:
            return 1./d, x
        return 1./d
        
    def T(self, *varg, quality=False, **kwarg):
        """Temperature
    T(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns temperature in unit_temperature

In many applications, it is also necessary to calculate quality to 
completely specify the state, and since it is an intermediate for any
property calculation, it can be returned as well.  If the optional 
"quality" keyword argument is set to True, x is appended in a tuple to 
save an unnecessary redundant call to x().

    T,x = T(..., quality=True)
"""
        T,_,_,x,_ = self._argparse(*varg, **kwarg)
        T = pm.units.temperature_scale(T, from_units='K')
        if quality:
            return T,x
        return T
        
    def x(self, *varg, **kwarg):
        """Quality
    x(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns quality, which is a dimensionless number between 0 and 1 for 
saturated mixtures and -1 for all other states.

In many applications quality is one of a few important properties.  To
avoid redundant function calls, consider using the "quality" keyword in
another property method or the state() method.
"""
        _,_,_,x,_ = self._argparse(*varg, **kwarg)
        return x
        
    #                    #
    # Property functions #
    #                    #
    
    def state(self, *varg, **kwarg):
        """The state method calculates all available properties at once.
        
    sd = state(...)
    
The properties are returned in a dictionary with keys:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    f   free energy         unit_energy / unit_matter
    g   gibbs energy        unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature
    cp  const. p sp. ht.    unit_energy / unit_matter / unit_temperature
    cv  const. v sp. ht.    unit_energy / unit_matter / unit_temperature
    
Like all of the other property functions, arguments may be any two of
T, p, d, v, e, h, s, and x.  

Because calculating cv for saturation conditions is more computationally
expensive, and because users rarely need this property, state() will
return NaN for cv at saturated conditions.  This is a deliberate design
decision to preserve the speed and simplicitly of the state() method.  
For users who do want true constant-volume specific heat is still 
available by calling the cv() method directly.
"""
        
        # Parse the arguments
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        R = self.data['R']
        
        # Initialize the output
        out = {}
        
        # Start with the vapor (d2) half of the calculation
        # In saturated cases, d2 should always be used to caluclate 
        # pressure
        # The IG part        
        Tscale = self.data['IGgroup']['Tscale']
        dscale = self.data['IGgroup']['dscale']
        tt = Tscale / T
        dd = d2 / dscale
        a,at,ad,att,atd,add = self._fo(tt,dd,2)
        
        p = 1.
        e = at
        h = 1. + tt*at
        s = tt*at - a
        cp = -tt*tt*att
        cv = tt*tt*att
        
        # The residual part
        Tscale = self.data['Rgroup']['Tscale']
        dscale = self.data['Rgroup']['dscale']
        tt = Tscale / T
        dd = d2 / dscale
        a,at,ad,att,atd,add = self._fr(tt,dd,2)

        p += dd*ad
        p *= T*d2*R
        e += at
        e *= R*Tscale
        h += dd*ad + tt*at
        h *= R*T
        s += tt*at - a
        s *= R
        temp = 1.+dd*(ad-tt*atd)
        cp += -tt*tt*att + temp*temp/(1.+dd*(2.*ad+dd*add))
        cp *= R
        cv += tt*tt*att
        cv *= -R
        
        # Before we go back and calculate the liquid properties,
        # go ahead and store the vapor calculations
        out['p'] = p
        out['T'] = T
        out['d'] = d1
        out['x'] = x
        out['e'] = e
        out['f'] = e - T*s
        out['g'] = h - T*s
        out['h'] = h
        out['s'] = s
        out['cp'] = cp
        out['cv'] = cv
        
        # Finish with the liquid (d1) half of the calculation
        # The IG part        
        Tscale = self.data['IGgroup']['Tscale']
        dscale = self.data['IGgroup']['dscale']
        tt = Tscale / T[I]
        dd = d1[I] / dscale
        a,at,ad,att,atd,add = self._fo(tt,dd,2)
        
        e = at
        h = 1. + tt*at
        s = tt*at - a
        cp = -tt*tt*att
        cv = tt*tt*att
        
        # The residual part
        Tscale = self.data['Rgroup']['Tscale']
        dscale = self.data['Rgroup']['dscale']
        tt = Tscale / T[I]
        dd = d1[I] / dscale
        a,at,ad,att,atd,add = self._fr(tt,dd,2)

        e += at
        e *= R*Tscale
        h += dd*ad + tt*at
        h *= R*T[I]
        s += tt*at - a
        s *= R
        temp = 1.+dd*(ad-tt*atd)
        cp += -tt*tt*att + temp*temp/(1.+dd*(2.*ad+dd*add))
        cp *= R
        cv += tt*tt*att
        cv *= -R
        
        # Finally, calculate the mixture properties with the appropriate
        # quality.
        out['cp'][I] = np.inf
        out['cv'][I] = np.nan
        out['e'][I] = out['e'][I]*(x[I]) + e*(1-x[I])
        out['h'][I] = out['h'][I]*(x[I]) + h*(1-x[I])
        out['s'][I] = out['s'][I]*(x[I]) + s*(1-x[I])
        # Overwrite the helmholtz function with the mixture values
        out['f'][I] = out['e'][I] - out['T'][I]*out['s'][I]
        # Gibbs energy is constant across an equilibrium phase transition
        # d is not weighted by x - v is.
        out['d'][I] = 1./((1-x[I])/d1[I] + x[I]/d2[I])
        
        # Apply unit conversions
        c1 = pm.units.energy(1., from_units='J')
        c1 = pm.units.matter(c1, self.data['mw'], from_units='kg', exponent=-1)
        out['e'] *= c1
        out['h'] *= c1
        out['f'] *= c1
        out['g'] *= c1
        c1 = pm.units.temperature(c1, from_units='K',exponent=-1)
        out['s'] *= c1
        out['cp'] *= c1
        out['cv'] *= c1
        out['gam'] = out['cp'] / out['cv']
        out['gam'][I] = np.inf
        out['p'] = pm.units.pressure(out['p'], from_units='Pa')
        out['T'] = pm.units.temperature_scale(out['T'], from_units='K')
        c1 = pm.units.volume(1., from_units='m3', exponent=-1)
        c1 = pm.units.matter(c1, self.data['mw'], from_units='kg')
        out['d'] *= c1
        out['v'] = 1./out['d']
        return out
        
        
    def e(self, *varg, quality=False, **kwarg):
        """Internal energy
    e(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Additionally, if the optional keyword, "quality" is set to True, the 
quality of the liquid/vapor mixture is also returned
    e,x = e(..., quality=True)

Returns energy in unit_energy / unit_matter
"""
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
        
    def f(self, *varg, quality=False, **kwarg):
        """Free (Helmholtz) energy
    f(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Additionally, if the optional keyword, "quality" is set to True, the 
quality of the liquid/vapor mixture is also returned
    f,x = f(..., quality=True)

Returns free energy in unit_energy / unit_matter
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        f = self._f(T,d1,0)[0]
        if I.any():
            f[I] *= (1.-x[I])
            f[I] += self._f(T[I],d2[I],0)[0] * x[I]
        # Convert the units back to user space
        pm.units.energy(f, from_units='J', inplace=True)
        pm.units.matter(f, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        if quality:
            return f,x
        return f

    def g(self, *varg, quality=False, **kwarg):
        """Gibbs energy
    g(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Additionally, if the optional keyword, "quality" is set to True, the 
quality of the liquid/vapor mixture is also returned
    g,x = g(..., quality=True)

Returns free energy in unit_energy / unit_matter
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        g = self._g(T,d1,0)[0]
        if I.any():
            g[I] *= (1.-x[I])
            g[I] += self._g(T[I],d2[I],0)[0] * x[I]
        # Convert the units back to user space
        pm.units.energy(g, from_units='J', inplace=True)
        pm.units.matter(g, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        if quality:
            return g,x
        return g    
        
    def h(self, *varg, quality=False, **kwarg):
        """Enthalpy
    h(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Additionally, if the optional keyword, "quality" is set to True, the 
quality of the liquid/vapor mixture is also returned
    h,x = h(..., quality=True)

Returns enthalpy as unit_energy / unit_matter
"""
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


    def s(self, *varg, quality=False, **kwarg):
        """Entropy
    s(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Additionally, if the optional keyword, "quality" is set to True, the 
quality of the liquid/vapor mixture is also returned
    s,x = s(..., quality=True)

Returns entropy in unit_energy / unit_matter / unit_temperature
"""
            
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


    def a(self, *varg, quality=False, **kwarg):
        """Speed of sound
    a(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Additionally, if the optional keyword, "quality" is set to True, the 
quality of the liquid/vapor mixture is also returned
    s,x = s(..., quality=True)

Returns speed of sound in unit_length / unit_time

The speed of sound in a two-phase mixture is not currently defined.  
Normally, the saturated state forms two separate regions of vapor and 
liquid, each with its own speed of sound, which should be calculated at
the saturation line.
"""
        
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        a = self._a(T,d1)
        if I.any():
            a[I] = pm.config['def_oob']
        # Convert the units back to user space
        pm.units.length(a, from_units='m', inplace=True)
        pm.units.time(a, from_units='s', inplace=True, exponent=-1)
        if quality:
            return a,x
        return a


    def hsd(self, *varg, quality = False, **kwarg):
        """Enthalpy, Entropy, Density
** Deprecated - Use state() **
        
    h,s,d = hsd(...)
        OR
    h,s,d,x = hsd(..., quality=True)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Additionally, if the optional keyword, "quality" is set to True, the 
quality of the liquid/vapor mixture is also returned
    e,x = e(..., quality=True)

"""
            
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        
        # There is no inner hsd funciton.  
        # We have to do this the hard way.
        
        # The IG part        
        R = self.data['R']
        Tscale = self.data['IGgroup']['Tscale']
        dscale = self.data['IGgroup']['dscale']
        tt = Tscale / T
        dd = d1 / dscale
        a,at,_,_,_,_ = self._fo(tt,dd,1)
        
        h = 1. + tt*at
        s = tt*at - a
        
        # The residual part
        Tscale = self.data['Rgroup']['Tscale']
        dscale = self.data['Rgroup']['dscale']
        tt = Tscale / T
        dd = d1 / dscale
        a,at,ad,_,_,_ = self._fr(tt,dd,1)
        h += dd*ad + tt*at
        s += tt*at - a

        # If there are data under the dome
        if I.any():
            temp = 1-x[I]
            h[I] *= temp
            s[I] *= temp
            
            # The IG part        
            R = self.data['R']
            Tscale = self.data['IGgroup']['Tscale']
            dscale = self.data['IGgroup']['dscale']
            tt = Tscale / T[I]
            dd = d2[I] / dscale
            a,at,_,_,_,_ = self._fo(tt,dd,1)
            
            h[I] += (1. + tt*at)*x[I]
            s[I] += (tt*at - a)*x[I]
            
            # The residual part
            Tscale = self.data['Rgroup']['Tscale']
            dscale = self.data['Rgroup']['dscale']
            tt = Tscale / T[I]
            dd = d2[I] / dscale
            a,at,ad,_,_,_ = self._fr(tt,dd,1)
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
        

    def cp(self, *varg, quality=False, **kwarg):
        """Constant-pressure specific heat
    cp(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Additionally, if the optional keyword, "quality" is set to True, the 
quality of the liquid/vapor mixture is also returned
    cp,x = cp(..., quality=True)

Note that constant-pressure specific heat is theoretically infinite for
saturated liquid-vapor mixtures.  cp() returns +Inf for any states that
are under the dome.

Returns specific heat in unit_energy / unit_matter / unit_temperature
"""
            
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        cp = self._cp(T,d1)
        if I.any():
            cp[I] = np.inf
        # Convert the units back to user space
        pm.units.energy(cp, from_units='J', inplace=True)
        pm.units.matter(cp, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        pm.units.temperature(cp, from_units='K', 
                exponent=-1, inplace=True)
        if quality:
            return cp, x
        return cp


    def cv(self, *varg, quality=False, **kwarg):
        """Constant-volume specific heat
    cv(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Additionally, if the optional keyword, "quality" is set to True, the 
quality of the liquid/vapor mixture is also returned
    cv,x = cv(..., quality=True)
    
The cv() method is unique in that it provides slightly different 
behaviors from its corresponding value returned by the state() method.
The state() method does not calculate specific heats of any kind for 
saturated conditions.  Meanwhile, cv() uses the derivatives saturation
density and internal energy to calculate the total mixture specific 
heat.  Applications that require this behavior should use cv() 
explicitly instead of depending on state().

Returns specific heat in unit_energy / unit_matter / unit_temperature
"""
        
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        cv = self._cv(T,d1)
        if I.any():
            # How do the saturation densities change with temperature?
            _,dVT,_ = self._dsv(T[I], diff=1)
            _,dLT,_ = self._dsl(T[I], diff=1)
            # How does x change with temperature
            temp = d1[I]/d2[I]
            xT = (dLT/d1[I]*(1-x[I]) + temp*dVT/d2[I]*x) / (temp-1)
            # Grab the saturation sensitivities
            eL,eLT,eLd = self._e(T[I],d1[I],diff=1)
            eV,eVT,eVd = self._e(T[I],d2[I],diff=1)
            # Calculate the true isochoric specific heat for the
            # two-phase mixture
            cv[I] = (eLT+eLd*dLT)*(1-x) + (eVT+eVd*dVT)*x + (eV-eL)*xT
            
        # Convert the units back to user space
        pm.units.energy(cv, from_units='J', inplace=True)
        pm.units.matter(cv, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        pm.units.temperature(cv, from_units='K', 
                exponent=-1, inplace=True)
        if quality:
            return cv, x
        return cv
        
        
    def gam(self, *varg, quality=False, **kwarg):
        """Specific heat ratio
    gam(...)

All properties accept two other properties as flexible inputs.
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    x   quality             dimensionless
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Additionally, if the optional keyword, "quality" is set to True, the 
quality of the liquid/vapor mixture is also returned
    gam,x = gam(..., quality=True)

Returns specific heat ratio, which is dimensionless
"""
            
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        cv = self._cv(T,d1)
        cp = self._cp(T,d1)
        if I.any():
            cp[I] = np.inf
        
        if quality:
            return cp/cv, x
        return cp/cv


    def T_s(self, s, p=None, d=None, quality=False, debug=False):
        """Temperature from entropy
** Deprecated - use T() **

    T = T_s(s, p=p)
        OR
    T = T_s(s, d=d)

If neither pressure nor density is specified, the default pressure will be 
used (config['def_p']).  

The optional keyword flag, quality, will cause quality to be returned
along with temperature.

    T,x = T_s(s, p=p, quality=True)
"""
        if p is not None:
            return self.T(s=s,p=p,quality=quality)
        elif d is not None:
            return self.T(s=s,d=d,quality=quality)
        p = pm.config['def_p']
        return self.T(s=s, p=p)


    def d_s(self, s, T=None, quality=False, debug=False):
        """Density from entropy
** Deprecated - use d() **

    d = d_s(s,T=T)
    
If temperature is not specified, the default temperature will be used 
(config['def_T']).

The optional keyword flag, quality, will cause quality to be returned along
with pressure.
"""
        if T is not None:
            return self.d(s=s, T=T, quality=quality)
        return self.d(s=s, quality=quality)
            



    def T_h(self, h, p=None, d=None, quality=False, debug=False):
        """Temperature from entropy
** Deprecated - use T() **

    T = T_s(s, p=p)
        OR
    T = T_s(s, d=d)

If neither pressure nor density is specified, the default pressure will be 
used (config['def_p']).  

The optional keyword flag, quality, will cause quality to be returned
along with temperature.

    T,x = T_s(s, p=p, quality=True)
"""
        if p is not None:
            return self.T(h=h,p=p,quality=quality)
        elif d is not None:
            return self.T(h=h,d=d,quality=quality)
        p = pm.config['def_p']
        return self.T(h=h, p=p)
