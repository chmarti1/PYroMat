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

def interp_multiple(x, xdata, *varg):
    """Perform 1D linear interpolation with multiple simultaneous ydata sets
    y0, y1, y2, ... = interp_multiple(x, xdata, y0data, y1data, ...)

This is equivalent to 
    y0 = interp(x, xdata, y0data)
    y1 = interp(x, xdata, y1data)
    ...
However, because redundant calls to interp cause the xdata array to be
searched repeatedly, this algorithm is far more efficient for two or
more interpolations.  

In tests with data arrays with 10,000 data elements gave these results
with arbitrary time units:
    interp()            1.00
    searchsorted()      0.94
    interp_multiple()   1.28    (with one data set)
    interp_multiple()   1.37    (with two data sets)
    interp_multiple()   1.46    (with three data sets)

The majority of the algorithm's time is spent on searching xdata, so
the efficiency lost by performing the interpolation in uncompiled code
is more than made up by stashing the search result.  Each additional
dataset only costs about 9% of one call to interp().
"""
    i1 = np.minimum(np.searchsorted(xdata, x), len(xdata)-1)
    i0 = i1-1
    # Advanced indexing makes an array copy.  We need x0 more than once
    # so stash the copy for efficiency.  We'll overwrite it as soon as 
    # we're done with it.
    t0 = xdata[i0]
    # Use dimensionless parameters, t0 and t1 = 1 - t0
    t1 = (x-t0)/(xdata[i1]-t0)
    t0[:] = 1-t1
    # Initialize the output
    output = []
    for ydata in varg:
        output.append(t1*ydata[i1] + t0*ydata[i0])
    return tuple(output)

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



class mp2(pm.reg.PYroMatModel):
    """The PYroMat multi-phase generalist class 2

** PROPERTY METHODS **
MP2 provides property methods:
    a()     Speed of sound
    cp()    Isobaric specific heat
    cv()    Isochoric specific heat
    d()     Density
    e()     Internal energy
    f()     Free (Helmholtz) energy
    g()     Gibbs energy
    gam()   Specific heat ratio
    h()     Enthalpy
    s()     Entropy
    T()     Temperature
    p()     Pressure
    d()     Density
    v()     Specific volume
    x()     Quality
    state() Calculates most properties
    
All of the above methods accept a standardized call signature.  See the
_argparse() method documentation for a detailed description:
    import pyromat as pm
    S = pm.get('AN_MP2_SUBATANCE')
    help(S._argparse)

For example, enthalpy might be called
    h(T=300., p=1.01325)
    h(T=300., d=990.)
    h(T=300., x=0.5)
    h(s=6., p=2.5)

In the back end, all properties are calculated from temperature and density,
so providing this interface flexibility has a numerical cost.  Once T and d
are known, additional property evaluations should always be made in terms
of them.

Most property pairs are supported, but several are not.  

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

** SATURATION PROPERTY METHODS **
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

It is also possible to obtain saturation properties by from the non-
saturated property methods (for example by passing the liquid and vapor
densities), but this can be less numerically precise than calling the 
saturation methods.  Standard property methods do not ``understand'' that
the state is constrained to be _precisely_ on the saturation line unless, 
quality is given explicitly.  If the algorithm required iteration, it will 
merely be ``close'' to the line.

** OTHER PROPERTY METHODS **
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
MP2 models thermo-physical properties of a liquid-gas system using a 
general fit for helmholtz free energy.  These "Span & Wagner" fits are 
evaluated in a polynomial form with exponential post factors.

The MP2 class is divided into three layers of methods (routines).  

--- USER ROUTINES ---
Accept data in any format (array or scalar) and in whatever units are
configured in the PYroMat configuration object.  These routines rely on
_argparse and _sat_argparse to standardize their call signatures, to
convert to the correct units, and to enforce that all inner routines
receive correctly broadcast ndarray objects.

Values from these methods are returned in appropriately broadcast arrays
in the correctly configured units.

--- INNER ROUTINES ---
Configured for speed and efficiency, these methods presume that all 
arguments are properly broadcast numpy arrays and that they are in a 
common unit system.  This prevents redundant calls to the _argparse() 
layey.  Units used by the MP2 back-end are:
    Energy:     J
    Matter      kg
    Pressure:   Pa
    Temperature:K
    
Inner routines begin with a "_" to emphasize that they are not part of
the standard interface, and their first line contains the text 
``(inner-routine)''.  Most property functions are wrappers for inner
routine property functions, so they may call each other when needed.  
Inner routine property functions (e.g. _h, _s, _p, etc...) have standard
call signatures that require temperature and density, and return the 
property and its derivatives to temperature and density.  For example:
    h, ht, hd = _h(T=T, d=d, diff=1)
h is enthalpy, ht is the derivative of enthallpy with respect to 
temeprature while holding density constant, and hd is the derivative of
enthalpy with respect to density while holding temperature constant.

VERY rarely, these routines might be called by the user to achieve 
supperior numerical performance.  Inner routines are significantly 
faster than the user routines because they do not have the overhead of
unit conversions, array broadcasting, and casing out the property 
combinations.  However, they have stringent requirements on their 
arguments:

1) All arguments must be a numpy NDARRAY object of dimension 1 or 
    greater.
2) Array broadcasting must be done BEFORE passing arguments to the inner
    routines.
3) The above units MUST be respected regardless of PYroMat's settings.
4) Many of these functions also return their derivatives to facilitate
    numerical inversion.  Check the documentation to verify the call
    signature of each inner routine BEFORE implementing it in your code.
5) Inner routines make no error checking for out-of-bounds or 
    saturation.  This allows access to metastable states, but it also
    allows users to naively query states that return utter nonsense.

--- PRIMITIVE ROUTINES ---
Methods that have been labeled as primitive routines should UNDER NO
CIRCUMSTANCES be called by the user.  They accept non-dimensionalized
arguments and return non-dimensional parameters.  These are encapsulated
as independent methods either because they are complicated and need to 
be called by a number of other methods, or because separating them made
sense for numerical efficiency.  They are also subject to change without 
warning in future upgrades.  In summary: these aren't the methods you're
looking for.

--- DATA DICTIONARY ---
The MP2 data dictionary must have certain data "groups" to define the 
various empirical fits.  Each group is a dictionary (within the 
dictionary) that defines the various parameters necessary for at least
one of the inner methods.

The Helmholtz free energy is nondimensionalized by RT, and calculated 
from groups briefly summarized below in terms of dimensionless 
temperature and density,
    tt = Tc / T     <== INVERSE!
    dd = d / dc

Data dictionaries have sub-dictionaries with members listed below:
-- Ideal Gas Group --
IGgroup         Helmholtz free energy ideal gas group; a dict containing:
    logt        a scalar coefficient of a log(tt) term
    coef0       a coefficient list to be passed to _poly1() to build p0 
    coef1       a simple Nx2 coefficient list used to build q(tt) below

The formula for the ideal gas portion of free energy is:
    fo = log(d) + LOGT*log(tt) + TLOGT*tt*log(tt) + p0(tt) + q(tt)
        q(tt) = sum_k coef1[k,1] * log(1 - exp(-tt*coef[k,0]))
    Fo = fo * R * T
where LOGT is the coefficient defined by the 'logt' parameter, and p is
the polynomial defined by the coef list

-- Residual Group --
Rgroup          Helmholtz free energy residual group; a dict containing:
    coef0       a nested list of coefficient lists
    coef1       an optional nested list of coefficients
    coef2       an optional nested list of coefficients

Each element of coef0 is, itself a coefficient list intended to be 
passed to _poly2().  After the first element, each individual polynomial
is multiplied by exp(-dd**k) where k is the index in the coef list.
    
    fr0 
   ----- = p0(tt,dd) + exp(-dd)*p1(tt,dd) + exp(-dd**2)*p2(tt,dd) + ...
    R T
    
coef1 is an optional list of lists of coefficients forming a matrix
[...
    [ t, d, b, a, gam, ep, c ], ...
]
    fr1 
   ----- = c * dd**d * tt**t * exp(-a*(dd-ep)**2 - b*(tt-gam)**2) + ...
    R T
    
If coef1 is defined it will be combined with the other coefficients
to form the residual.  If coef1 is not defined, it will be ignored.

coef2 is an optional list of lists of coefficients forming a 2D array
[...
    [ a, b, m, A, B, C, D, c ], ...
]
In the evaluation of coef2, there is an intermediate term, X

    X = ((1-tt) + A*((dd-1)**2)**(0.5/m))**2 + B*((dd-1)**2)**a
    
used to calculate the dimensionless free energy

    fr2 
   ----- = c * X**b * d * exp(-C*(dd-1)**2 - D*(tt-1)**2) + ...
    R T

Additionally, there are numerous parameters that define global 
properties -- parameters that do not vary with state.

Tlim            A two-element list of the upper and lower temperatures
                for which the data set is valid.
plim            A two-element list of the upper and lower pressures for
                which the data set is valid.
dlim            A two-element list the represent practical maximum and
                minimum densities over the entire data set.  These are 
                NOT guaranteed limits of validity.
Tc, dc          Critical temperature and density
Tt              Triple-point temperature
R               (optional) Ideal gas constant in J/kg/K
mw              Molecular weight
atoms           A dictionary with a key for each atom and a value for 
                its count in the molecule.  For example, CO2 would 
                have atoms = {'C':1, 'O':2}
                
If the gas constant, R, is not provided, it will be calculated from 
molecular weight and pm.units.const_Ru.  Providing it as a data value
allows the model to be evaluated using precisely the same parameters 
used by the authors of the original models.
                
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
        

    ########################
    #                      #
    #  Numerical Routines  #
    #                      #
    ########################

    def _poly2(self,x,y,pcoef,diff=2):    
        """Polynomial evaluation (primitive routine)
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
        """Polynomial evaluation (primitive routine)
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
        """Search 1D map for an inverse estimates (primitive routine)
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

    def _mapsearch2(self, f0data, f1data, f0value, f1value, zde0=0, zde1=0):
        r"""Search 2D map for inverse estimates (primitive routine)
    T, d, Isat, Ioob = mapsearch2(f0data, f1data, f0value, f1value)
    
Uses tabulated data to generate an estimates for x,y in the 2D inversion
problem
    f0(T,d) = f0value
    f1(T,d) = f1value

ARGUMENTS:    
f0data, f1data
    Two-dimensional array containing tabulated values for f0(T,d) and 
    f1(T,d) from the substance's _table dict.  The indices should be 
    arranged so that
        f0[i,j] = f0(T[i], d[j])
        f1[i,j] = f1(T[i], d[j])
    where the T and d arrays are the temperature and density values 
    itemized in the substance _table dict.
    
f0value, f1value
    Numpy arrays with the same shape containing values for properties,
    f0data and f1data.
    
zde0, zde1  (Default 0)
    Zero-density extrapolation method -- an integer specifying how 
    values found to line between density index 0 and 1 should be 
    treated.  Enthalpy and internal energy converge to their ideal gas
    values, but entropy and any property derived from it diverges like 
    ln(d).  The following values are accepted:
    0 - Use standard linear interpolation (default)
            f(d) = f(d[1])-f(d[0]) * (d-d[0]) / (d[1]-d[0])
    1 - Use entropy extrapolation: 
            f(d) = f(d=d[1]) - R*ln(d/d[1])
    2 - Use free energy extrapolation:
            f(d) = f(d=d[1]) + T*R*ln(d/d[1])
    
RETURNS: 
T,d
    Arrays of the same shape as fvalue and gvalue that approximate 
    solutions to the problem
        f0(T,d) =approx= f0value
        f1(T,d) =approx= f1value

Isat
    A boolean array of the same shape as the xvalue and yvalue arrays,
    indicating states at which the estimated solution is either 
    saturated or very nearly saturated.  If Isat is False, the state is
    definiately NOT saturated.
    
Ioob
    A boolean array of the same shape as the xvalue and yvalue arrays,
    indicating states that are out-of-bounds of the substance data map.

DESCRIPTION:

The f0data and f1data are 2D arrays of tabulated values of f0(x,y) and 
f1(x,y) in a rectangular grid of x and y values.  This is notably 
distinct from 2D interpolation because the maps, fdata and gdata, do not
need to be monotonically increasing.  The algorithm performs a global 
search by explicitly comparing all node values:
    f0value < f0_ij
    f1value < f1_ij

Grid elements containing potential solutions are identified as those 
with at least one node above and below the target values for both f0() 
and f1().  Then estimates are generated by finding the approximate 
intersections of the paths in x,y implied by the f0() and f1() 
constraints inside the element.  First, the element's edges are 
interpolated to find estimates for two points where f0(x,y)==f0value and
f1(x,y)==f1value.  The intersection (if one exists) of the two resulting
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

ABOUT LINE SEGMENT INTERPOLATION:

Line segment interpolation was selected over the usual bilinear 
interpolation because of its linearity.  Bilinear 2D element 
interpolation is obnoxious to invert because of its nonlinear xy term,
which can cause saddle points and other irritating issues.  However, 
line segment interpolation still suffers from problems, which are 
mitigated in this algorithm:
(1) When the solution lies precisely on a node, one line segment 
    vanishes, leading to a singular problem.  This is mitigated by 
    explicitly testing for precise equality at the nodes.
(2) When solution estimates lie very close to the element edge, tiny 
    numerical errors can cause redundant estimates from neighboring
    elements or the estimate can be omitted altogether.  When estimates
    are a small distance from an element's edge (even if it is very 
    slightly outside) it is included.  
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
    cases are detected and discarded with a warning.
    
A number of versions of _mapsearch2() were tested. This version simply 
returns the first solution discovered.  Other versions faithfully 
reported multiple candidate solutions if they were discovered.  The 
design intent is for _argparse() to weed out cases that might have 
multiple solutions.  Still, some special cases (especially h,s) have
strage edge cases where multiple solutions creep in.

SEE ALSO:
    _mapsearch1(), _mapsearch2(), _mapsearch2x(), _mapsearch2y()
"""
        # Define an increment for small values
        # For most systems, eps is about 2.2e-16, so small will be about
        # 2.2e-12.  This is the number we use to detect dimensionless
        # proximity to the element boundary.
        small = np.finfo(float).eps * 1e4
        # Initialize lists for the result values
        T = np.full_like(f0value, pm.config['def_oob'], dtype=float)
        d = np.full_like(f0value, pm.config['def_oob'], dtype=float)
        TI = np.empty_like(f0value, dtype=int)
        DI = np.empty_like(f0value, dtype=int)
        Isat = np.zeros_like(f0value, dtype=bool)
        Ioob = np.ones_like(f0value, dtype=bool)
        
        # Retrieve the T and d tabular arrays
        Tdata = self._table['T']
        ddata = self._table['d']
        Tci, dci = self._table['cI']
        
        # Keep a flag to indicate whether the user should be warned about
        # out-of-bounds values
        warn = pm.config['warning_verbose']
        for index in range(f0value.size):
            f0v = f0value.flat[index]
            f1v = f1value.flat[index]

            # Generate a boolean array indicating candidate elements with a solution
            # Bulk element comparison seems expensive, but it is not on a 
            # system with vectorized processing.  Bulk comparisons like this
            # are remarkably cheap. 
            f0I = f0v < f0data
            f1I = f1v < f1data
            I = np.logical_and(crossing2(f0I), crossing2(f1I))

            # For each element that contains a crossing in both f and g
            for Ti,di in zip(*np.nonzero(I)):
                Ti1 = Ti+1
                di1 = di+1

                # Identify the two f-edge crossings [(x,y), ...]
                f0cross = self._intersect(f0data, f0I, f0v, Ti, di, zde=zde0)
                # Identify the two g-edge crossings [(x,y), ...]
                f1cross = self._intersect(f1data, f1I, f1v, Ti, di, zde=zde1)
                # At this point, f0cross and f1cross list (x,y) coordinates for 
                # the points along the element edge where crossings occur
                # Meanwhile, neighbor lists the (xi,yi) indices of the
                # elements that share the edges where f() has a solution
                # We'll use neighbor to resolve conflict over solutions
                # very close to the edges.
                
                # Detect the saddle case
                if len(f0cross) != 2 or len(f1cross) != 2:
                    # For now, warn the user, and DO NOT append the case
                    pm.utility.print_warning('mp2._mapsearch2: Discarded a potential solution near a saddle point.  If you believe this was a legitimate solution, please report the code that generated this warning to the PYroMat GitHub issues page.')
                # Two edges have intersections for each function
                else:
                    fx0 = f0cross[0]
                    fdx = f0cross[1] - f0cross[0]
                    gx0 = f1cross[0]
                    gdx = f1cross[1] - f1cross[0]
                    #print('')
                    #print('xi,yi,x,y:', xi,yi,xdata[xi], ydata[yi])
                    #print('fdata values:', fdata[xi,yi], fdata[xi1,yi], fdata[xi,yi1], fdata[xi1,yi1])
                    #print('gdata values:', gdata[xi,yi], gdata[xi1,yi], gdata[xi,yi1], gdata[xi1,yi1])
                    
                    # Check for a solution precisely at the corner
                    if (fdx == 0).all():
                        if (f1cross[0] == fx0).all() or (f1cross[1] == fx0).all():
                            # When this code was modified to merely return the first solution discovered,
                            # these lines were commented out.  Return them if multiple solutions are 
                            # desired in the future.
                            #I[*neighbor[0]] = False
                            #I[*neighbor[1]] = False
                            T.flat[index] = fx0[0]
                            d.flat[index] = fx0[1]
                            TI.flat[index] = Ti
                            DI.flat[index] = di
                            Ioob.flat[index] = False
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
                                T.flat[index] = fx0[0] + s*fdx[0]
                                d.flat[index] = fx0[1] + s*fdx[1]
                                TI.flat[index] = Ti
                                DI.flat[index] = di
                                Ioob.flat[index] = False
                                break

        # If operating verbosely, warn the user about out-of-bounds elements
        if Ioob.any() and pm.config['warning_verbose']:
            pm.utility.print_warning('mp2._mapsearch2: Failed to find value(s) in the table. Result is out-of-bounds.')
        

        # Identify any element indices under the dome
        k = self._table['cI'][0] - TI
        dLi = self._table['cI'][1] + k
        dVi = self._table['cI'][1] - k
        Isat = (TI>=0) * (k>0) * (dVi <= DI) * (DI < dLi)

        return T,d,Isat,Ioob


    def _intersect(self, fdata, fI, fvalue, Ti, di, zde=0):
        r"""Helper method for the _mapsearch2() method (primitive routine)
    [(T0,d0), (T1,d1)] = _intersect(fdata, fvalue, fI, Ti, di, zde=0)
    
Calculates the (T,d) coordinates of points on an element's edges where
the specified property data interpolates to equal the scalar fvalue.

          Ti,di+1
            +--------+ Ti+1,di+1
            |        o (T1,d1)
    (T0,d0) o        |
            |        |
      Ti,di +--------+
                     Ti+1,di

_intersect() accepts arguments:

fdata
    The 2D property data array taken from the _table dictionary

fI
    A 2D array of boolean values indicating the result of the comparison
    fI = fvalue < fdata.  _mapsearch2() has already performed this 
    operation, so repeating it to determine which edges are crossed is
    redundant.

fvalue
    The desired scalar value of the property.
    
Ti, di 
    Temperature and density indices for the element being searched.  As
    in the figure above, the element index corresponds to the indices
    of the lower-left node in the rectangular element.
    
zde     (Default 0)
    Zero-density extrapolation algorithm to use
    0 - Use standard linear interpolation
    1 - Use entropy extrapolation
    2 - Use free-energy extrapolation

Returns a list of two-element tuples.  If no intersections are found,
the list is empty.  Two intersections are expected, but four are 
possible in saddle node cases.
"""
        ddata = self._table['d']
        Tdata = self._table['T']
        R = self.data['R']
        # Stash the upper indices
        Ti1 = Ti + 1
        di1 = di + 1
        
        # Initialize the result
        fcross = []
        # Test each of the edges for a crossing of f()
        # Bottom edge
        if fI[Ti,di] != fI[Ti1,di]:
            # Bottom-edge intersection is not possible with zde enabled
            # Never extrapolate.
            TT = interp_scalar(fvalue, fdata[Ti,di], fdata[Ti1,di], Tdata[Ti], Tdata[Ti1])
            fcross.append(np.array((TT,ddata[di])))
        # Left edge
        if fI[Ti,di] != fI[Ti,di1]:
            # Entropy extrapolation
            if di == 0 and zde == 1:
                dd = ddata[1] * np.exp((fdata[Ti,1] - fvalue)/R)
            # Free-energy extrapolation
            elif di == 0 and zde == 2:
                dd = ddata[1] * np.exp((fvalue - fdata[Ti,1])/R/Tdata[Ti])
            else:
                dd = interp_scalar(fvalue, fdata[Ti,di], fdata[Ti,di1], ddata[di], ddata[di1])
            fcross.append(np.array((Tdata[Ti], dd)))
        # Top edge
        if fI[Ti,di1] != fI[Ti1,di1]:
            # There is no need to perform extrapolation on the top edge under any circumstances
            TT = interp_scalar(fvalue, fdata[Ti,di1], fdata[Ti1,di1], Tdata[Ti], Tdata[Ti1])
            fcross.append(np.array((TT,ddata[di1])))
        # Right edge
        if fI[Ti1,di] != fI[Ti1,di1]:
            # Entropy extrapolation
            if di == 0 and zde == 1:
                dd = ddata[1] * np.exp((fdata[Ti1,1] - fvalue)/R)
            # Free-energy extrapolation
            elif di == 0 and zde == 2:
                dd = ddata[1] * np.exp((fvalue - fdata[Ti1,1])/R/Tdata[Ti1])
            else:
                dd = interp_scalar(fvalue, fdata[Ti1,di], fdata[Ti1,di1], ddata[di], ddata[di1])
            
            fcross.append(np.array((Tdata[Ti1], dd)))
        return fcross

    def _dmapsearch2(self, fdata, dvalue, fvalue, zde=0):
        r"""Constant-density 2D map search (primitive routine)
    T, Isat, Ioob = _dmapsearch2(fdata, dvalue, fvalue, zde=0)
    
Uses tabulated data to generate an estimate for T in the 2D inversion
problem
    f(T,dvalue) = fvalue

ARGUMENTS:
fdata
    Two-dimensional array-like containing tabulated values for f(T,d).  
    The indices should be arranged so that
        fdata[i,j] = f(T[i], d[j])
    where T and d are the tabulated temperature and density values in 
    the substance _table dict.
        
dvalue
    An array of density values to use when scanning the table.
    
fvalue
    An array of f-values to interpolate from the table.  The dimensions
    must match the dimensions of dvalue.
    
zde     (0)
    Zero-density extrapolation method -- an integer specifying how 
    values found to line between density index 0 and 1 should be 
    treated.  Enthalpy and internal energy converge to their ideal gas
    values, but entropy and any property derived from it diverges like 
    ln(d).  The following values are accepted:
    0 - Use standard linear interpolation (default)
            f(d) = f(d[1])-f(d[0]) * (d-d[0]) / (d[1]-d[0])
    1 - Use entropy extrapolation: 
            f(d) = f(d=d[1]) - R*ln(d/d[1])
    2 - Use free energy extrapolation:
            f(d) = f(d=d[1]) + T*R*ln(d/d[1])
        
RETURNS: 
T
    An array of temperatures that approximately solve the problem.
        

Isat
    A boolean array of the same shape as the xvalue and yvalue arrays,
    indicating states at which the estimated solution is either 
    saturated or very nearly saturated.  If Isat is False, the state is
    definiately NOT saturated.
    
Ioob
    A boolean array of the same shape as the xvalue and yvalue arrays,
    indicating states that are out-of-bounds of the substance data map.
    
DESCRIPTION:

Similarly to _mapsearch2, _dmapsearch2 looks for intersections of the
curves implied by
    f(T, d) = fvalue
    d = dvalue
cross.  Inside of elements, the f(T,d)=fvalue curve is interpolated 
linearly between the points where it crosses along the element edges.

Unlike _mapsearch2, _dmapsearch2 does not need to search the entire 
domain for solutions - it only performs operations on the row of 
elements implied by the d-value.  As a result, it is faster.

SEE ALSO:
    _mapsearch1(), _mapsearch2(), _dmapsearch2(), _Tmapsearch2()
"""
        # Define an increment for small values
        # For most systems, eps is about 2.2e-16, so small will be about
        # 2.2e-12.  This is the number we use to detect dimensionless
        # proximity to the element boundary.
        small = np.finfo(float).eps * 1e4

        # Get Tdata and ddata
        Tdata = self._table['T']
        ddata = self._table['d']
        # Initialize result arrays
        T = np.full_like(fvalue, pm.config['def_oob'], dtype=float)
        TI = np.full_like(fvalue, -1, dtype=int)
        # Find indices for the density locations in the array
        DI = np.searchsorted(ddata, dvalue, side='right')-1
        Ioob = np.ones_like(fvalue, dtype=bool)
        Isat = np.zeros_like(fvalue, dtype=bool)

        for index in range(fvalue.size):
            # Scalar density and property values
            dv = dvalue.flat[index]
            fv = fvalue.flat[index]
            # Scalar density index
            di = DI.flat[index]
            di1 = di + 1
            
            # Case out the density location
            # If it is out-of-bounds, do nothing
            if dv < ddata[0] or dv > ddata[-1]:
                pass
            # If entropy zero-density extrapolation is selected
            elif zde == 1 and di == 0:
                # Extrapolate to form a function of temperature along
                # the constant-density line
                fex = fdata[:, 1] - self.data['R'] * np.log(dv / ddata[1])
                # Test for crossings with the property value
                fI = fv < fex
                I = fI[:-1] != fI[1:]
                Ti = np.nonzero(I)[0]
                # If at least one crossing is identified, take the lowest
                if len(Ti) > 0:
                    Ti = Ti[0]
                    Ti1 = Ti+1
                    # Interpolate to identify the temperature estimate
                    T.flat[index] = interp_scalar(fv, fex[Ti], fex[Ti1], Tdata[Ti], Tdata[Ti1])
                    TI.flat[index] = Ti
                    Ioob.flat[index] = False
                # If there are no crossings, do nothing -- this is oob
            # If free-energy zero-density-extrapolation is selected
            elif zde == 2 and di == 0:
                # Extrapolate to form a function of temperature along
                # the constant-density line
                fex = fdata[:, 1] + Tdata * self.data['R'] * np.log(dv / ddata[1])
                # Test for crossings with the property value
                I = np.diff(fv < fex)
                Ti = np.nonzero(I)[0]
                # If at least one crossing is identified, take the lowest
                if len(Ti) > 0:
                    Ti = Ti[0]
                    Ti1 = Ti+1
                    # Interpolate to identify the temperature estimate
                    T.flat[index] = interp_scalar(fv, fex[Ti], fex[Ti1], Tdata[Ti], Tdata[Ti1])
                    TI.flat[index] = Ti
                    Ioob.flat[index] = False
            # The standard linear interpolation algorithm
            else:
                # Compare the values of only the appropriate row
                fI = fv < fdata[:, di:di+2]
                # Detect elements with a crossing
                I = crossing2(fI)
                for Ti in np.nonzero(I)[0]:
                    Ti1 = Ti+1
                    # Initialize some crossing parameters
                    fcross = []
                    # Detect the edges
                    # Bottom Edge
                    if fI[Ti,0] != fI[Ti1,0]:
                        TT = interp_scalar(fv, fdata[Ti,di], fdata[Ti1,di], Tdata[Ti], Tdata[Ti1])
                        fcross.append(np.array([TT, ddata[di]]))
                    # Left Edge
                    if fI[Ti,0] != fI[Ti,1]:
                        dd = interp_scalar(fv, fdata[Ti,di], fdata[Ti,di1], ddata[di], ddata[di1])
                        fcross.append(np.array([Tdata[Ti], dd]))
                    # Top Edge
                    if fI[Ti,1] != fI[Ti1,1]:
                        TT = interp_scalar(fv, fdata[Ti,di1], fdata[Ti1,di1], Tdata[Ti], Tdata[Ti1])
                        fcross.append(np.array([TT, ddata[di1]]))
                    # Right Edge
                    if fI[Ti1,0] != fI[Ti1,1]:
                        dd = interp_scalar(fv, fdata[Ti1,di], fdata[Ti1,di1], ddata[di], ddata[di1])
                        fcross.append(np.array([Tdata[Ti1], dd]))
                    # Detect the saddle case
                    if len(fcross) != 2:
                        # For now, warn the user, and DO NOT append the case
                        pm.utility.print_warning('mp2._dmapsearch2: Discarded a potential solution near a saddle point.  If you believe this was a legitimate solution, please report the code that generated this warning to the PYroMat GitHub issues page.')
                    # Two edges have intersections for each function
                    else:
                        fx0 = fcross[0]
                        fdx = fcross[1] - fcross[0]
                        # Detect precise equality at a corner
                        if (fdx == 0).all() and fx0[1] == dv:
                            T.flat[index] = fx0[0]
                            TI.flat[index] = Ti
                            Ioob.flat[index] = False
                            break
                        else:
                            # Calculate the distance along the f=0 curve to intersect 
                            # Perform the calculations in two steps - leave the division
                            # for last, so we can detect nearly singular problems
                            s = dv - fx0[1]
                            det = fdx[1]
                            
                            if 2*abs(det) > abs(s):
                                s /= det
                                if -small < s < 1+small:
                                    T.flat[index] = fx0[0] + fdx[0] * s
                                    TI.flat[index] = Ti
                                    Ioob.flat[index] = False
                                    break
        
        if pm.config['warning_verbose'] and Ioob.any():
            pm.utility.print_warning('mp2._dmapsearch2: Property value(s) were out-of-bounds.')
                    
        # Identify any element indices under the dome
        k = self._table['cI'][0] - TI
        dLi = self._table['cI'][1] + k
        dVi = self._table['cI'][1] - k
        Isat = (TI>=0) * (k>0) * (dVi <= DI) * (DI < dLi)
            
        return T, Isat, Ioob
        
    def _Tmapsearch2(self, fdata, Tvalue, fvalue):
        r"""Search 2D map for inverse estimates (primitive routine)
    d, Isat, Ioob = Tmapsearch2(fdata, Tvalue, fvalue)
    
Uses tabulated data to generate an estimate for y in the 2D inversion
problem
    f(T, d) = fvalue
    T = Tvalue

ARGUMENTS:
fdata
    Two-dimensional array-like containing tabulated values for f(T,d).  
    The indices should be arranged so that
        fdata[i,j] = f(T[i], d[j])
    where T and d are the tabulated temperature and density values in 
    the substance _table dict.
        
Tvalue
    An array of temperature values to use when scanning the table.
    
fvalue
    An array of f-values to interpolate from the table.  The dimensions
    must match the dimensions of Tvalue.
    
RETURNS: 
T
    An array of temperatures that approximately solve the problem.
        

Isat
    A boolean array of the same shape as the xvalue and yvalue arrays,
    indicating states at which the estimated solution is either 
    saturated or very nearly saturated.  If Isat is False, the state is
    definiately NOT saturated.
    
Ioob
    A boolean array of the same shape as the xvalue and yvalue arrays,
    indicating states that are out-of-bounds of the substance data map.
    
DESCRIPTION:

Similarly to _mapsearch2, _dmapsearch2 looks for intersections of the
curves implied by
    f(T, d) = fvalue
    T = Tvalue
cross.  Inside of elements, the f(T,d)=fvalue curve is interpolated 
linearly between the points where it crosses along the element edges.

Unlike _mapsearch2, _Tmapsearch2 does not need to search the entire 
domain for solutions - it only performs operations on the row of 
elements implied by the d-value.  As a result, it is faster.

SEE ALSO:
    _mapsearch1(), _mapsearch2(), _dmapsearch2(), _Tmapsearch2()
"""
        # Define an increment for small values
        # For most systems, eps is about 2.2e-16, so small will be about
        # 2.2e-12.  This is the number we use to detect dimensionless
        # proximity to the element boundary.
        small = np.finfo(float).eps * 1e4
        
        # Get the temperature and density data
        Tdata = self._table['T']
        ddata = self._table['d']
        
        # Initialize result arrays
        d = np.empty_like(fvalue, dtype=float)
        DI = np.empty_like(fvalue, dtype=int)
        TI = np.searchsorted(Tdata, Tvalue, side='right')
        Isat = np.zeros_like(fvalue, dtype=bool)
        Ioob = np.ones_like(fvalue, dtype=bool)
        
        for index in range(fvalue.size):
            fv = fvalue.flat[index]
            Tv = Tvalue.flat[index]
            # Halt if the temeprature value is out-of-bounds
            if Tdata[0] <= Tv <= Tdata[-1]:
                Ti1 = TI.flat[index]
                Ti = Ti1 - 1
                # Compare the values of only the appropriate row
                fI = fv < fdata[Ti:Ti+2, :]
                # Detect elements with a crossing
                I = crossing2(fI)
                for di in np.nonzero(I)[1]:
                    di1 = di+1
                    # Initialize some crossing parameters
                    fcross = []
                    # Detect the edges
                    # Bottom Edge
                    if fI[0,di] != fI[1,di]:
                        TT = interp_scalar(fv, fdata[Ti,di], fdata[Ti1,di], Tdata[Ti], Tdata[Ti1])
                        fcross.append(np.array([TT, ddata[di]]))
                    # Left Edge
                    if fI[0,di] != fI[0,di1]:
                        dd = interp_scalar(fv, fdata[Ti,di], fdata[Ti,di1], ddata[di], ddata[di1])
                        fcross.append(np.array([Tdata[Ti], dd]))
                    # Top Edge
                    if fI[0,di1] != fI[1,di1]:
                        TT = interp_scalar(fv, fdata[Ti,di1], fdata[Ti1,di1], Tdata[Ti], Tdata[Ti1])
                        fcross.append(np.array([TT, ddata[di1]]))
                    # Right Edge
                    if fI[1,di] != fI[1,di1]:
                        dd = interp_scalar(fv, fdata[Ti1,di], fdata[Ti1,di1], ddata[di], ddata[di1])
                        fcross.append(np.array([Tdata[Ti1], dd]))
                    # Detect the saddle case
                    if len(fcross) != 2:
                        # For now, warn the user, and DO NOT append the case
                        pm.utility.print_warning('mp2._Tmapsearch2: Discarded a potential solution near a saddle point.  If you believe this was a legitimate solution, please report the code that generated this warning to the PYroMat GitHub issues page.')
                    # Two edges have intersections for each function
                    else:
                        fx0 = fcross[0]
                        fdx = fcross[1] - fcross[0]
                        # Detect precise equality at a corner
                        if (fdx == 0).all() and fx0[0] == Tv:
                            d.flat[index] = fx0[1]
                            DI.flat[index] = di
                            Ioob.flat[index] = False
                            break
                        else:
                            # Calculate the distance along the f=0 curve to intersect 
                            # Perform the calculations in two steps - leave the division
                            # for last, so we can detect nearly singular problems
                            s = Tv - fx0[0]
                            det = fdx[0]
                            
                            if 2*abs(det) > abs(s):
                                s /= det
                                if -small < s < 1+small:
                                    d.flat[index] = fx0[1] + fdx[1] * s
                                    DI.flat[index] = di
                                    Ioob.flat[index] = False
                                    break
        if pm.config['warning_verbose'] and Ioob.any():
            pm.utility.print_warning('mp2._Tmapsearch2: Property value(s) were out-of-bounds.')\
                    
        # Identify any element indices under the dome
        k = self._table['cI'][0] - TI
        dLi = self._table['cI'][1] + k
        dVi = self._table['cI'][1] - k
        Isat = (TI>=0) * (k>0) * (dVi <= DI) * (DI < dLi)
        
        return d, Isat, Ioob



    def _Tsatiter(self, T, dL, dV, Ids, Nmax=20, ep=1e-6, debug=False):
        """Iterates on Maxwell's criteria while holding T constant (primitive routine)
    _Tsatiter(T, dL, dV, Ids)

T       Saturation temperature used to specify the state.
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
        for count in range(Nmax):
            # Create down-selected views
            T_ = T[Ids]
            dL_ = dL[Ids]
            dV_ = dV[Ids]
            # Get the dimensionless 
            argL = self._ff(T_, dL_, diff=2)
            argV = self._ff(T_, dV_, diff=2)
            
            gL,gLt,gLd = self._g(*argL, diff=1)
            gV,gVt,gVd = self._g(*argV, diff=1)
            pL,pLt,pLd = self._p(*argL, diff=1)
            pV,pVt,pVd = self._p(*argV, diff=1)

            # Initialize an error vector and a jacobian matrix
            e = np.empty(T_.shape + (2,1), dtype=float)
            J = np.empty(T_.shape + (2,2), dtype=float)            
            # Error vector
            # The vapor pressure is stored in p
            e[:,0,0] = gL - gV
            e[:,1,0] = pL - pV
            # Jacobian
            J[:,0,0] = gLd
            J[:,0,1] = -gVd
            J[:,1,0] = pLd
            J[:,1,1] = -pVd
            # Overwrite error with the perturbation to the estimates
            delta = np.linalg.solve(J,e)
            
            # Update unknowns
            dL_ -= delta[:,0,0]
            dV_ -= delta[:,1,0]
            
            # Test for densities that have overshot the critical point
            Ioob = (dV_ > self.data['dc']) + (dL_ < self.data['dc'])
            inner_count = 0
            while Ioob.any():
                inner_count += 1
                if inner_count > Nmax:
                    raise pm.utility.PMAnalysisError(f'mp2._Tsatiter: Crossed the critical point, and failed to produce a valid estimate after {Nmax} divisions!')
                if debug:
                    print(f'  Overstep correction {inner_count}')
                e[Ioob,...] /= 2
                dL_[Ioob] += delta[Ioob,0,0]
                dV_[Ioob] += e[Ioob,1,0]
                Ioob = (dV_ > self.data['dc']) + (dL_ < self.data['dc'])

            if debug:
                print(f'**{count}**')
                print('dL:', dL_)
                print('dV:', dV_)
                print('delta:,', delta)
            
            # Update results
            dV[Ids] = dV_
            dL[Ids] = dL_
            
            # Detect convergence
            Ids[Ids] = np.logical_or( np.abs(delta[:,0,0]) > ep*dL_,
                    np.abs(delta[:,1,0]) > ep*dV_ )
            
            if not Ids.any():
                return

        raise pm.utility.PMAnalysisError(f'_Tsatiter: Failed to converge in {Nmax} iterations.')
        

    def _dVsatiter(self, T, dL, dV, Ids, Nmax=20, ep=1e-6, debug=False):
        """Iterates on Maxwell's criteria while holding dV constant (primitive routine)
    _dVsatiter(T, dL, dV, Ids)

T       Saturation temperature used to specify the state.
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
        for count in range(Nmax):
            # Create down-selected views
            T_ = T[Ids]
            dL_ = dL[Ids]
            dV_ = dV[Ids]
            
            # Evaluate the properties at the liquid and vapor lines
            argL = self._ff(T_, dL_, diff=2)
            argV = self._ff(T_, dV_, diff=2)
            
            gL,gLt,gLd = self._g(*argL, diff=1)
            gV,gVt,gVd = self._g(*argV, diff=1)
            pL,pLt,pLd = self._p(*argL ,diff=1)
            pV,pVt,pVd = self._p(*argV ,diff=1)
            
            # Initialize arrays for the linear algebra
            e = np.empty(T_.shape + (2,1), dtype=float)
            J = np.empty(T_.shape + (2,2), dtype=float)            
            # Build the Jacobian on temperature and liquid density
            J[:,0,0] = gLt-gVt
            J[:,0,1] = gLd
            J[:,1,0] = pLt-pVt
            J[:,1,1] = pLd
            # Build the error vector
            e[:,0,0] = gL-gV        # Gibbs error
            e[:,1,0] = pL-pV        # Pressure error
            # Solve.  Ovewrite error with the estimate perturbation
            delta = np.linalg.solve(J,e)
            # Update temperature and density
            T_ -= delta[:,0,0]
            dL_ -= delta[:,1,0]
            # Test for densities that have overshot the critical point
            Ioob = (dL_ < self.data['dc']) + (T_ > self.data['Tc'])
            inner_count = 0
            while Ioob.any():
                inner_count += 1
                if inner_count > Nmax:
                    raise pm.utility.PMAnalysisError(f'mp2._dVsatiter: Crossed the critical point, and failed to produce a valid estimate after {Nmax} divisions!')
                if debug:
                    print(f'  Overstep correction {inner_count}')
                delta[Ioob,...] /= 2
                T_[Ioob] += delta[Ioob,0,0]
                dL[Ioob] += delta[Ioob,1,0]
                Ioob = (dL_ < self.data['dc']) + (T_ > self.data['Tc'])
            
            if debug:
                print(f'**{count}**')
                print('T:', T_)
                print('dL:', dL_)
                print('delta:,', delta)
            
            # Update the results
            T[Ids] = T_
            dL[Ids] = dL_
            
            # Test for convergence
            Ids[Ids] = np.logical_or(np.abs(delta[:,0,0]) > ep*T_, np.abs(delta[:,1,0]) > ep*dL_)
            
            # If all points have converged
            if not Ids.any():
                return
        
        raise pm.utility.PMAnalysisError(f'_dVsatiter: Failed to converge in {Nmax} iterations.')

    def _dLsatiter(self, T, p, dL, dV, Ids, Nmax=20, ep=1e-6, debug=False):
        """Iterates on Maxwell's criteria while holding dL constant (primitive routine)
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
        for count in range(Nmax):
            # Create down-selected views
            T_ = T[Ids]
            dL_ = dL[Ids]
            dV_ = dV[Ids]
            
            # Evaluate the properties at the liquid and vapor lines
            argL = self._ff(T_, dL_, diff=2)
            argV = self._ff(T_, dV_, diff=2)
            
            gL,gLt,gLd = self._g(*argL, diff=1)
            gV,gVt,gVd = self._g(*argV, diff=1)
            pL,pLt,pLd = self._p(*argL ,diff=1)
            pV,pVt,pVd = self._p(*argV ,diff=1)
            
            # Initialize arrays for the linear algebra
            e = np.empty(T.shape + (2,1), dtype=float)
            J = np.empty(T.shape + (2,2), dtype=float)
            
            # Build the Jacobian on temperature and liquid density
            J[:,0,0] = gLt-gVt
            J[:,0,1] = -gVd
            J[:,1,0] = pLt-pVt
            J[:,1,1] = -pVd
            # Build the error vector
            e[:,0,0] = gL-gV        # Gibbs error
            e[:,1,0] = pL-p[Ids]    # Pressure error
            # Solve.  Ovewrite error with the estimate perturbation
            delta = np.linalg.solve(J,e)
            # Update temperature and density
            T_ -= delta[:,0,0]
            dV_ -= delta[:,1,0]
            # Test for densities that have overshot the critical point
            Ioob = (dV_ > self.data['dc']) + (T_ > self.data['Tc'])
            inner_count = 0
            while Ioob.any():
                inner_count += 1
                if inner_count > Nmax:
                    raise pm.utility.PMAnalysisError(f'mp2._dLsatiter: Crossed the critical point, and failed to produce a valid estimate after {Nmax} divisions!')
                if debug:
                    print(f'  Overstep correction {inner_count}')
                e[Ioob,...] /= 2
                T_[Ioob] += delta[Ioob,0,0]
                dV_[Ioob] += delta[Ioob,1,0]
                Ioob = (dV_ > self.data['dc']) + (T_ > self.data['Tc'])
            
            if debug:
                print(f'**{count}**')
                print('T:', T_)
                print('dV:', dV_)
                print('delta:,', e[Ids])
            
            # Test for convergence
            Ids[Ids] = np.logical_or(np.abs(delta[:,0,0]) > ep*T_, np.abs(delta[:,1,0]) > ep*dV_)
            
            # If all points have converged
            if not Ids.any():
                return
        
        raise pm.utility.PMAnalysisError(f'_dLsatiter: Failed to converge in {Nmax} iterations.')


    def _psatiter(self, T, dL, dV, p, Ids, Nmax=20, ep=1e-6, debug=False):
        """Iterates on Maxwell's criteria while holding p constant (primitive routine)
    _psatiter(T, dL, dV, p, Ids)

T       Saturation temperature.
dL      Liquid density.
dV      Vapor density.
p       Pressure used to determine the saturation state.
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
        fail = True
        for count in range(Nmax):
            # Generate views of the updated down-selected variables
            T_ = T[Ids]
            dL_ = dL[Ids]
            dV_ = dV[Ids]
            p_ = p[Ids]
            
            # Evaluate the properties at the liquid and vapor lines
            argL = self._ff(T_, dL_, diff=2)
            argV = self._ff(T_, dV_, diff=2)
            
            gL,gLt,gLd = self._g(*argL, diff=1)
            gV,gVt,gVd = self._g(*argV, diff=1)
            pL,pLt,pLd = self._p(*argL ,diff=1)
            pV,pVt,pVd = self._p(*argV ,diff=1)

            # Initialize arrays for the linear algebra
            e = np.empty(T_.shape + (3,1), dtype=float)
            J = np.empty(T_.shape + (3,3), dtype=float)
            
            # Error vector
            e[:,0,0] = gL - gV
            e[:,1,0] = pL - p_
            e[:,2,0] = pV - p_
            # Jacobian
            J[:,0,0] = gLt-gVt
            J[:,0,1] = gLd
            J[:,0,2] = -gVd
            
            J[:,1,0] = pLt
            J[:,1,1] = pLd
            J[:,1,2] = 0.
            
            J[:,2,0] = pVt
            J[:,2,1] = 0.
            J[:,2,2] = pVd
            # Calculate change in the variables
            delta = np.linalg.solve(J,e)
            
            # Update the variables
            T_ -= delta[:,0,0]
            dL_ -= delta[:,1,0]
            dV_ -= delta[:,2,0]
            
            # Test for densities that have overshot the critical point
            Ioob = (dV_ > self.data['dc']) + (dL_ < self.data['dc']) + (T_ > self.data['Tc'])
            inner_count = 0
            while Ioob.any():
                inner_count += 1
                if inner_count > Nmax:
                    raise pm.utility.PMAnalysisError(f'mp2._psatiter: Crossed the critical point, and failed to produce a valid estimate after {Nmax} divisions!')
                if debug:
                    print(f'  Overstep correction {inner_count}')
                delta[Ioob,...] /= 2
                T_[Ioob] += delta[Ioob,0,0]
                dL_[Ioob] += delta[Ioob,1,0]
                dV_[Ioob] += delta[Ioob,2,0]
                Ioob = (dV_ > self.data['dc']) + (dL_ < self.data['dc']) + (T_ > self.data['Tc'])
            
            if debug:
                print(f'**{count}**')
                print('T:', T_)
                print('dL:', dL_)
                print('dV:', dV_)
                print('delta:,', e)
            
            # Update the results
            T[Ids] = T_
            dL[Ids] = dL_
            dV[Ids] = dV_
            
            # Detect convergence
            Ids[Ids] = np.logical_or( np.abs(e[:,0,0]) > ep*T_,
                        np.logical_or( np.abs(e[:,1,0]) > ep*dL_,
                        np.abs(e[:,2,0]) > ep*dV_))
            if not Ids.any():
                return
            
        raise pm.utility.PMAnalysisError(f'_psatiter: Failed to converge in {Nmax} iterations.')


    def _gsatiter(self, T, dL, dV, g, Ids, Nmax=20, ep=1e-6, debug=False):
        """Iterates on Maxwell's criteria while holding p constant (primitive routine)
    _gsatiter(T, dL, dV, g, Ids)

T       Saturation temperature.
dL      Liquid density.
dV      Vapor density.
g       Gibbs energy used to determine the state
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
        # Make a copy of the 
        Istash = Ids.copy()
        
        fail = True
        for count in range(Nmax):
            # Generate views of the updated down-selected variables
            T_ = T[Ids]
            dL_ = dL[Ids]
            dV_ = dV[Ids]
            g_ = g[Ids]
            
            # Evaluate the properties at the liquid and vapor lines
            argL = self._ff(T_, dL_, diff=2)
            argV = self._ff(T_, dV_, diff=2)
            
            gL,gLt,gLd = self._g(*argL, diff=1)
            gV,gVt,gVd = self._g(*argV, diff=1)
            pL,pLt,pLd = self._p(*argL ,diff=1)
            pV,pVt,pVd = self._p(*argV ,diff=1)

            # Initialize arrays for the linear algebra
            e = np.empty(T_.shape + (3,1), dtype=float)
            J = np.empty(T_.shape + (3,3), dtype=float)
            
            # Error vector
            e[:,0,0] = pL - pV
            e[:,1,0] = gL - g_
            e[:,2,0] = gV - g_
            # Jacobian
            J[:,0,0] = pLt-pVt
            J[:,0,1] = pLd
            J[:,0,2] = -pVd
            
            J[:,1,0] = gLt
            J[:,1,1] = gLd
            J[:,1,2] = 0.
            
            J[:,2,0] = gVt
            J[:,2,1] = 0.
            J[:,2,2] = gVd
            # Calculate change in the variables
            delta = np.linalg.solve(J,e)
            
            # Update the variables
            T_ -= delta[:,0,0]
            dL_ -= delta[:,1,0]
            dV_ -= delta[:,2,0]
            
            # Test for densities that have overshot the critical point
            Ioob = (dV_ > self.data['dc']) + (dL_ < self.data['dc']) + (T_ > self.data['Tc'])
            inner_count = 0
            while Ioob.any():
                inner_count += 1
                if inner_count > Nmax:
                    raise pm.utility.PMAnalysisError(f'mp2._gsatiter: Crossed the critical point, and failed to produce a valid estimate after {Nmax} divisions!')
                if debug:
                    print(f'  Overstep correction {inner_count}')
                delta[Ioob,...] /= 2
                T_[Ioob] += delta[Ioob,0,0]
                dL_[Ioob] += delta[Ioob,1,0]
                dV_[Ioob] += delta[Ioob,2,0]
                Ioob = (dV_ > self.data['dc']) + (dL_ < self.data['dc']) + (T_ > self.data['Tc'])
            
            if debug:
                print(f'**{count}**')
                print('T:', T_)
                print('dL:', dL_)
                print('dV:', dV_)
                print('delta:,', e)
            
            # Update the results
            # NOTE: This stores the pressure value PRIOR to applying the last delta
            T[Ids] = T_
            dL[Ids] = dL_
            dV[Ids] = dV_
            
            # Detect convergence
            Ids[Ids] = np.logical_or( np.abs(e[:,0,0]) > ep*T_,
                        np.logical_or( np.abs(e[:,1,0]) > ep*dL_,
                        np.abs(e[:,2,0]) > ep*dV_))
            
            if not Ids.any():
                return
            
        raise pm.utility.PMAnalysisError(f'_gsatiter: Failed to converge in {Nmax} iterations at {np.sum(Ids)} values.')


    def _satiter2(self, T, dL, dV, fn0, fn1, f0value, f1value, Ids, Nmax=50, ep=1e-6, debug=False):
        """Two-property saturation iteration (primitive routine)
    _satiter2(self, T, dL, dV, fn0, fn1, f0value, f1value, Ids, Nmax=20, ep=1e-6)

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

        count = 0
        while Ids.any():
            count += 1
            if count > Nmax:
                raise pm.utility.PMParamError(
                        f'mp2._satiter2: Failed to converge after {Nmax} iterations at {np.sum(Ids)} value(s).')
            
            TT = T[Ids]
            DL = dL[Ids]
            DV = dV[Ids]
            
            # Evaluate the properties
            argL = self._ff(TT, DL, diff=2)
            argV = self._ff(TT, DV, diff=2)
            
            pL,pLt,pLd = self._p(*argL, diff=1)
            pV,pVt,pVd = self._p(*argV, diff=1)
            gL,gLt,gLd = self._g(*argL, diff=1)
            gV,gVt,gVd = self._g(*argV, diff=1)
            f0L,f0Lt,f0Ld = fn0(*argL, diff=1)
            f0V,f0Vt,f0Vd = fn0(*argV, diff=1)
            f1L,f1Lt,f1Ld = fn1(*argL, diff=1)
            f1V,f1Vt,f1Vd = fn1(*argV, diff=1)

            # Property deltas across the dome
            df0 = f0V - f0L
            vf0 = f0value[Ids] - f0L
            df1 = f1V - f1L
            vf1 = f1value[Ids] - f1L

            E = np.empty(TT.shape + (3,1), dtype=float)
            J = np.empty(TT.shape + (3,3), dtype=float)

            E[:, 0, 0] = pV - pL       # Maxwell, pressure
            E[:, 1, 0] = gV - gL       # Maxwell, gibbs energy
            E[:, 2, 0] = vf0*df1 - vf1*df0     # Quality constraint
            
            J[:, 0, 0] = pVt - pLt
            J[:, 0, 1] = -pLd
            J[:, 0, 2] = pVd
            
            J[:, 1, 0] = gVt - gLt
            J[:, 1, 1] = -gLd
            J[:, 1, 2] = gVd
            
            J[:, 2, 0] = -f0Lt*df1 + vf0*(f1Vt - f1Lt) + f1Lt*df0 - vf1*(f0Vt - f0Lt)
            J[:, 2, 1] = -f0Ld*df1 - vf0*f1Ld + f1Ld*df0 + vf1*f0Ld
            J[:, 2, 2] = vf0*f1Vd - vf1*f0Vd
            
            delta = np.linalg.solve(J, E)
            TT -= delta[:,0,0]
            DL -= delta[:,1,0]
            DV -= delta[:,2,0]

            # Test for densities that have overshot the critical point
            Ioob = (DV > self.data['dc']) + (DL < self.data['dc']) + (TT > self.data['Tc'])
            inner_count = 0
            while Ioob.any():
                inner_count += 1
                if inner_count > Nmax:
                    raise pm.utility.PMAnalysisError(f'mp2._satiter2: Crossed the critical point, and failed to produce a valid estimate after {Nmax} divisions!')
                if debug:
                    print(f'  Overstep correction {inner_count}')
                delta[Ioob,...] /= 2
                TT[Ioob] += delta[Ioob,0,0]
                DL[Ioob] += delta[Ioob,1,0]
                DV[Ioob] += delta[Ioob,2,0]
                Ioob = (DV > self.data['dc']) + (DL < self.data['dc']) + (TT > self.data['Tc'])
            
            if debug:
                print(f'**{count}**')
                print('T:', TT)
                print('dL:', DL)
                print('dV:', DV)
                print('delta:,', delta)
            
            T[Ids] = TT
            dL[Ids] = DL
            dV[Ids] = DV
            
            # Update convergence criteria
            Ids[Ids] = (delta[:,0,0] > TT*ep) + (delta[:,1,0] > DL*ep) + (delta[:,2,0] > DV*ep)
            
            
    def _dsatiter2(self, T, dL, dV, d, fn, fvalue, Ids, ep=1e-6, Nmax=20, debug=False):
        """Iterate on saturation properties to achieve mix density and one inverse (primitive routine)
    _dsatiter2(T, dL, dV, x, d, fn, fvalue, Ids, ep=1e-6, Nmax=20)
    
T       Temperature array used as an initial guess
dL      Saturated liquid array used as an initial guess
dV      Saturated vapor array used as an initial guess
d       Target density mixture array - not written to
fn      The inverse property's method
fvalue  The inverse property value array
Ids     Boolean down-select array

**DESCRIPTION**
Solves the problem 
    p(T,dV) = p(T,dL)
    g(T,dV) = g(T,dV)
    1/d = x/dV + (1-x)/dL
    fvalue = x*fn(T,dV) + (1-x)*fn(T,dL)
    
Given a guess for T, dL, and dV, the quality required to respect the 
density constraint can be calculated explicitly, leaving three nonlinear
constraints.
"""

        count = 0
        while Ids.any():
            count += 1
            if count > Nmax:
                raise pm.utility.PMParamError(
                        f'mp2._dsatiter2: Failed to converge after {Nmax} iterations.')
            
            TT = T[Ids]
            DL = dL[Ids]
            DV = dV[Ids]
            # Evaluate the properties
            # Evaluate the properties
            argL = self._ff(TT, DL, diff=2)
            argV = self._ff(TT, DV, diff=2)
            
            pL,pLt,pLd = self._p(*argL, diff=1)
            pV,pVt,pVd = self._p(*argV, diff=1)
            gL,gLt,gLd = self._g(*argL, diff=1)
            gV,gVt,gVd = self._g(*argV, diff=1)
            fL,fLt,fLd = fn(*argL, diff=1)
            fV,fVt,fVd = fn(*argV, diff=1)

            # Calculate quality and its derivatives from density
            dd = d[Ids]
            den = DL/DV-1
            xV = (DL/dd - 1)/den
            xVL = 1./dd/den - xV/DV/den
            xVV = xV*DL/DV/DV/den
            xL = 1 - xV
            xLL = -xVL
            xLV = -xVV
            
            E = np.empty(TT.shape + (3,1), dtype=float)
            J = np.empty(TT.shape + (3,3), dtype=float)
            
            E[:, 0, 0] = pL - pV       # Maxwell, pressure
            E[:, 1, 0] = gL - gV       # Maxwell, gibbs energy
            E[:, 2, 0] = fvalue[Ids] - xV*fV - xL*fL
            
            J[:, 0, 0] = pVt - pLt
            J[:, 0, 1] = -pLd
            J[:, 0, 2] = pVd
            
            J[:, 1, 0] = gVt - gLt
            J[:, 1, 1] = -gLd
            J[:, 1, 2] = gVd
            
            J[:, 2, 0] = xV*fVt + xL*fLt
            J[:, 2, 1] = xVL*fV + xLL*fL + xL*fLd
            J[:, 2, 2] = xVV*fV + xV*fVd + xLV*fL
            
            delta = np.linalg.solve(J, E)
            T[Ids] += delta[:,0,0]
            dL[Ids] += delta[:,1,0]
            dV[Ids] += delta[:,2,0]
            
            if debug:
                print(f'**{count}**')
                print('T:', TT)
                print('dL:', DL)
                print('dV:', DV)
                print('delta:,', delta[Ids])
            
            # Update convergence criteria
            Ids[Ids] = (delta[:,0,0] > TT*ep) + (delta[:,1,0] > DL*ep) + (delta[:,2,0] > DV*ep)

    def _Titer(self, T, d, fn, fvalue, Ids, Nmax=20, ep=1e-6, debug=False):
        """Constant-temperature iteration (primitive routine)
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
            
            arg = self._ff(TT,DD,diff=2)
            f,ft,fd = fn(*arg, diff=1)
            
            dd = (fvalue[Ids] - f) / fd
            d[Ids] += dd
            if debug:
                print(f'**{count}**')
                print('d:', DD)
                print('delta:', dd)
            
            Ids[Ids] = np.abs(dd) > ep * DD


    def _diter(self, T, d, fn, fvalue, Ids, Nmax=20, ep=1e-6, debug=False):
        """Constant-density iteration (primitive routine)
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
            
            arg = self._ff(TT,DD,diff=2)
            f,ft,fd = fn(*arg, diff=1)
            
            dT = (fvalue[Ids] - f) / ft
            T[Ids] += dT
            if debug:
                print(f'**{count}**')
                print('T:', TT)
                print('delta:', dT)
            
            Ids[Ids] = np.abs(dT) > ep * TT


    def _iter2(self, T, d, f0, f1, f0value, f1value, Ids, Nmax=20, ep=1e-6, debug=False):
        """Constant-density iteration (primitive routine)
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
                        f'mp2._iter2: Failed to converge after {Nmax} iterations at {np.sum(Ids)} value(s).')
            
            DD = d[Ids]
            TT = T[Ids]
            
            # We'll use f and g as placeholder function values
            arg = self._ff(TT,DD, diff=2)
            f,ft,fd = f0(*arg, diff=1)
            g,gt,gd = f1(*arg, diff=1)
            
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
            if debug:
                print(f'**{count}**')
                print('T:', T)
                print('d:', d)
                print('deltas:', dT, dd)
            
            # Update the convergence tests
            Ids[Ids] = np.logical_and(np.abs(dT) > ep * TT, np.abs(dd) > ep * DD)

    ##############################
    #                            #
    # Fundamental Property Model #
    #                            #
    ##############################

    def _fo(self, tt, dd, diff=2):
        """Dimensionless ideal gas helmholtz free energy (primitive routine)
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
        """Dimensionless residual helmhotz free energy (primitive routine)
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


    def _ff(self, T, d, diff=2):
        """Wrapper function for the dimensionless free energy methods
    tt,dd,a,at,ad,att,atd,add = _ff(T,d,diff=2)
    
Sums the free energy and its derivatives from the ideal gas and residual
components.  See _fo() and _fr() for more information.

Arguments
T       Temperature array in Kelvin
d       Density array in kg/m3

Returns 
tt, dd  -   Dimensionless temperature and density
a,at,ad,att,atd,add - Dimensionless free energy and its derivatives

The return signature is such that any property can be efficiently 
evaluated by passing the returned tuple as its ordered-value argumnet.

For example, this sequence calculates both enthalpy and internal energy
with only one call to the dimensionless back-end.

arg = _ff(T,d)
h = _h(*arg)
e = _e(*arg)
"""
        tt = self.data['Tc'] / T
        dd = d / self.data['dc']
        a,at,ad,att,atd,add = self._fo(tt,dd,diff=diff)
        b,bt,bd,btt,btd,bdd = self._fr(tt,dd,diff=diff)
        return tt, dd, a+b, at+bt, ad+bd, att+btt, atd+btd, add+bdd


    ###################
    #                 #
    #  Build Methods  #
    #                 #
    ###################

    def _build_sattab(self, step=0.02, ep=1e-6, verbose=False, debug=False):
        """Generate saturation table values (primitive routine)
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

        # We'll identify the lower temperature limit based on the triple
        # point or the data limits -- whichever is higher
        Tt = self.data['Tlim'][0]
        if 'Tt' in self.data:
            Tt = max(self.data['Tt'], Tt)
        Tc = self.data['Tc']
        dc = self.data['dc']
        
        # Initialize the outputs
        Ts_array = [Tc]
        dsL_array = [dc]
        dsV_array = [dc]
        
        if verbose:
            print('T pc dL dV')
            print('Critical Point:')
            print(f'{Tc:8.2f} {dc:8.2f} {dc:12.4e}')
        
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
        Ids = np.array([True],dtype=bool)
        # Create an initial perturbation of the densities
        # Do not perturb temperature
        dL += step * dc / 1.414
        dV -= step * dc / 1.414
        fail = True
        for count in range(200):
            # Iterate with constant dV
            Ids[0] = True
            self._dVsatiter(T, dL, dV, Ids, debug=debug)
            
            Ts_array.insert(0, T[0])
            dsL_array.insert(0, dL[0])
            dsV_array.insert(0, dV[0])
            
            if verbose:
                print(f'{T[0]:8.2f} {dL[0]:8.2f} {dV[0]:12.4e}')
            
            # Perturb the solution to the next interval
            # Assume a unity change in dV, calculate other changes
            ddV = -1.
            # Use the Maxwell criteria and its derivatives to construct
            # a Jacobian and a perturbation vector assuming a unity 
            # change in vapor density.
            argL = self._ff(T=T, d=dL, diff=2)
            argV = self._ff(T=T, d=dV, diff=2)
            gL,gLt,gLd = self._g(*argL,diff=1)
            gV,gVt,gVd = self._g(*argV,diff=1)
            pL,pLt,pLd = self._p(*argL,diff=1)
            pV,pVt,pVd = self._p(*argV,diff=1)
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
            self._Tsatiter(T, dL, dV, Ids, debug=debug)
            
            Ts_array.insert(0, T[0])
            dsL_array.insert(0, dL[0])
            dsV_array.insert(0, dV[0])
            
            if verbose:
                print(f'{T[0]:8.2f} {dL[0]:8.2f} {dV[0]:12.4e}')
            
            # Perturb the solution to the next interval
            # Assume a unity change in dV, calculate other changes
            dT = -1.
            # Use the Maxwell criteria and its derivatives to construct
            # a Jacobian and a perturbation vector assuming a unity 
            # change in vapor density.
            argL = self._ff(T=T, d=dL, diff=2)
            argV = self._ff(T=T, d=dV, diff=2)
            gL,gLt,gLd = self._g(*argL,diff=1)
            gV,gVt,gVd = self._g(*argV,diff=1)
            pL,pLt,pLd = self._p(*argL,diff=1)
            pV,pVt,pVd = self._p(*argV,diff=1)
            J[0,0] = gLd[0]
            J[0,1] = -gVd[0]
            J[1,0] = pLd[0]
            J[1,1] = -pVd[0]
            B[0] = (gVt[0] - gLt[0])*dT
            B[1] = (pVt[0] - pLt[0])*dT
            # Solve for the corresponding changes in T and dL
            x = np.linalg.solve(J,B)
            ddL = x[0]
            dVV = -x[1]/dV/dV     # Near the triple point, we'll perterb vapor volume instead of density
            # Rescale the steps so that the metric T/Tc, d/dc is equal to step
            scale = step / np.sqrt(dT*dT/Tc/Tc + ddL*ddL/dc/dc)
            dT *= scale
            ddL *= scale
            dVV *= scale

            # Detect the exit condition
            # If the next guess would be beyond the triple point, halt
            if T[0] + dT < Tt:
                fail=False
                break

            T += dT
            dL += ddL
            dV = 1./(1./dV + dVV)   # Perterb volume rather than density
            
        if fail:
            pm.utility.print_error('This error should never appear in a release - please report this on the PYroMat Github Issues page.')
            raise pm.utility.PMDataError('_build_sattab: Iteration froze near the triple point.' )
        
        scale = (Tt - T[0]) / dT
        ddL *= scale
        dVV *= scale
        
        T += dT
        dL += ddL
        dV = 1./(1./dV + dVV)
        Ids[0] = True
        
        self._Tsatiter(T, dL, dV, Ids)
    
        if verbose:
            print('Triple Point:')
            print(f'{T[0]:8.2f} {dL[0]:8.2f} {dV[0]:12.4e}')
        
        Ts_array.insert(0, T[0])
        dsL_array.insert(0, dL[0])
        dsV_array.insert(0, dV[0])
        
        if verbose:
            print(f'Used {len(Ts_array)} points.')
            print('Populating property lists...')
        
        # Convert to Numpy arrays
        Ts_array = np.array(Ts_array)
        dsL_array = np.array(dsL_array)
        dsV_array = np.array(dsV_array)
        
        argL = self._ff(Ts_array, dsL_array, diff=1)
        argV = self._ff(Ts_array, dsV_array, diff=1)
        
        self._sattable = {
            'T':Ts_array, 
            'p':self._p(*argV, diff=0)[0], 
            'dL':dsL_array, 
            'dV':dsV_array,
            'eL':self._e(*argL, diff=0)[0],
            'eV':self._e(*argV, diff=0)[0],
            'hL':self._h(*argL, diff=0)[0],
            'hV':self._h(*argV, diff=0)[0],
            'sL':self._s(*argL, diff=0)[0],
            'sV':self._s(*argV, diff=0)[0],
            'g':self._g(*argV, diff=0)[0],
            'fL':self._f(*argL, diff=0)[0],
            'fV':self._f(*argV, diff=0)[0]
        }
        
        if verbose:
            print('Done')


    def _build_tab(self, NT=100, Nd=100, verbose=False):
        """Generate lookup tables (primitive routine)
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
2) No step between any two temperature or density values may be larger 
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
        # Polish with constant-temperature far from the critical point
        I = dsV < 0.5 * dc
        Ids = np.array(I)
        self._Tsatiter(Ts, dsL, dsV, Ids)
        # Polish with constant-density near the critical point
        if verbose:
            print('Constant-vapor-density polishing near the critical point...')
        Ids = np.logical_not(I)
        Ids[-1] = False     # Do not polish the critical point
        self._dVsatiter(Ts, dsL, dsV, Ids)


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
        arg = self._ff(T=TT, d=dd, diff=1)
        p = self._p(*arg,diff=0)[0]
        e = self._e(*arg,diff=0)[0]
        h = self._h(*arg,diff=0)[0]
        s = self._s(*arg,diff=0)[0]
        
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
            # First, broadcast the constant properties, p
            p[iT, iL+1:iV] = p[iT,iV]
            # Next, use quality to calculate the mixture properties
            e[iT, iL+1:iV] = e[iT, iL]*xL + e[iT, iV]*xV
            h[iT, iL+1:iV] = h[iT, iL]*xL + h[iT, iV]*xV
            s[iT, iL+1:iV] = s[iT, iL]*xL + s[iT, iV]*xV
        
        # Build the table dictionary
        self._table = {'T':T, 'd':d, 'cI':(Tci, dci), 'p':p, 'e':e, 'h':h, 's':s}
        if verbose:
            print('Done.')

    def _build(self, force=False):
        """Build the back-end tables (inner routine)
    _build(force=False)

Creates the _table and _sattable attributes and populates them with data
required for the proper function of the MP2 class.  If they are already
created, they will not be re-created unless the force keyword is set to
True.
"""
        if force or ('_table' not in self.__dict__) or ('_sattable' not in self.__dict__):
            if pm.config['dat_verbose']:
                print(f'{self.sid()}: Populating saturation table...')
            self._build_sattab()
            if pm.config['dat_verbose']:
                print(f'{self.sid()}: Populating table...')
            self._build_tab()
            if pm.config['dat_verbose']:
                print(f'{self.sid()}: Done.')            

    ##############################
    #                            #
    #  Inner Saturation Methods  #
    #                            #
    ##############################


    def _Tsat(self, T, debug=False):
        """Calculate saturation state from temperature (inner routine)
    T, dL, dV = _Tsat(T)
    
Calculates saturation densities and their derivatives from temperature.  
    T       A numpy array of temperatures. All values must be between
            the triple point and critical point.  For speed, this is 
            not verified, so unexpected behaviors or failures will 
            result.
Returns:
    T       The original T array
    dL      Saturated liquid density
    dV      Saturated vapor density

**DESCRIPTION**
Interpolates the _sattable data for initial guesses of the saturation 
properties, then polishes with _Tsatiter().  The return values and their
order are to preserve a standard call signature for all saturation 
routines.

**SEE ALSO**
    _Tsat(), _psat(), _gsat(), _dLsat(), _dVsat()
"""
        dL, dV = interp_multiple(T, self._sattable['T'], self._sattable['dL'], self._sattable['dV'])

        I = np.ones_like(T, dtype=bool)
        self._Tsatiter(T, dL, dV, I, debug=debug)
        return T, dL, dV
        
        
    def _psat(self, p, debug=False):
        """Saturation state from pressure (inner routine)
    T, dL, dV = _psat(p)
    
Calculates the saturation state from pressure
    p       A Numpy array of pressure values.  All values must be 
            between the triple point and critical point.  For speed, 
            this is not verified, so unexpected behaviors or failures 
            will result.
Returns:
    T       The saturation temperature array
    dL      Saturated liquid density array
    dV      Saturated vapor density array

**DESCRIPTION**
Interpolates the _sattable data for initial guesses of the saturation 
properties, then polishes with _psatiter().  The return values and their
order are to preserve a standard call signature for all saturation 
routines.

**SEE ALSO**
    _Tsat(), _psat(), _gsat(), _dLsat(), _dVsat()
"""
        T, dL, dV = interp_multiple(p, self._sattable['p'], 
                self._sattable['T'], self._sattable['dL'], self._sattable['dV'])
        I = np.ones_like(p, dtype=bool)
        self._psatiter(T, dL, dV, p, I, debug=debug)
        return T, dL, dV
        
    def _gsat(self, g, debug=False):
        """Saturation state from Gibbs energy (inner routine)
    T, dL, dV = _gsat(g)
    
Calculates the saturation state from Gibbs energy
    g       A Numpy array of Gibbs energy values.  All values must be 
            between the triple point and critical point.  For speed, 
            this is not verified, so unexpected behaviors or failures 
            will result.
Returns:
    T       The saturation temperature array
    p       The same array passed to _psat()
    dL      Saturated liquid density array
    dV      Saturated vapor density array

**DESCRIPTION**
Interpolates the _sattable data for initial guesses of the saturation 
properties, then polishes with _gsatiter().  The return values and their
order are to preserve a standard call signature for all saturation 
routines.

**SEE ALSO**
    _Tsat(), _psat(), _gsat(), _dLsat(), _dVsat()
"""
        T, dL, dV = interp_multiple(g, self._sattable['g'], 
                self._sattable['T'], self._sattable['dL'], self._sattable['dV'])
        I = np.ones_like(g, dtype=bool)
        self._gsatiter(T, dL, dV, g, I, debug=debug)
        return T, dL, dV
        
    def _dLsat(self, dL, debug=False):
        """Saturation state from liquid density (inner routine)
    T, dL, dV = _dLsat(dL)
    
Calculates the saturation state from pressure
    dL      A Numpy array of saturated liquid density values.  All 
            values must be between the triple point and critical point 
            liquid densities.  For speed, this is not verified, so 
            unexpected behaviors or failures will result.
Returns:
    T       The saturation temperature array
    dL      The same saturated liquid array passed to _dLsat()
    dV      Saturated vapor density array

**DESCRIPTION**
Interpolates the _sattable data for initial guesses of the saturation 
properties, then polishes with _dLsatiter().  The return values and their
order are to preserve a standard call signature for all saturation 
routines.

**SEE ALSO**
    _Tsat(), _psat(), _gsat(), _dLsat(), _dVsat()
"""
        T, dV = interp_multiple(dL, np.flip(self._sattable['dL']), 
                np.flip(self._sattable['T']), np.flip(self._sattable['dV']))

        I = np.ones_like(dL, dtype=bool)
        self._dLsatiter(T, dL, dV, I, debug=debug)
        return T, dL, dV


    def _dVsat(self, dV, debug=False):
        """Saturation state from liquid density (inner routine)
    T, dL, dV = _dVsat(dV)
    
Calculates the saturation state from pressure
    dV      A Numpy array of saturated vapor density values.  All values
            must be between the triple point and critical point vapor 
            densities.  For speed, this is not verified, so unexpected 
            behaviors or failures will result.
Returns:
    T       The saturation temperature array
    dL      The same saturated liquid array passed to _dLsat()
    dV      Saturated vapor density array

**DESCRIPTION**
Interpolates the _sattable data for initial guesses of the saturation 
properties, then polishes with _dVsatiter().  The return values and their
order are to preserve a standard call signature for all saturation 
routines.

**SEE ALSO**
    _Tsat(), _psat(), _gsat(), _dLsat(), _dVsat()
"""
        T, dL = interp_multiple(dV, self._sattable['dV'], self._sattable['T'], self._sattable['dL'])

        I = np.ones_like(dV, dtype=bool)
        self._dVsatiter(T, dL, dV, I, debug=debug)
        return T, dL, dV


    ############################
    #                          #
    #  Inner Property Methods  #
    #                          #
    ############################


    def _R(self):
        """Obtain the ideal gas constant
    R = _R()
    
Returns the ideal gas constant in units J/kg/K.  

Most published models explicitly provide a value for R that was used in 
the model's development.  When specified along with the substance's 
molecular weight, these can be in small numerical contradiction with the
precise definition of the universal ideal gas constant (defined 
precisely by Boltzmann's constant and Avagadro's number).  

The theoretical relationship is always

    R = Ru / mw

When R is in J/kg/K, Ru is in J/kmol/K, and mw is in kg/kmol.

If the data dictionary includes a value for R, it is returned verbatim.
If not, R is calculated from the value in pm.units.const_Ru.
"""
        R = self.data.get('R')
        if R is None:
            R = 1000 * pm.units.const_Ru / self.data['mw']
        return R
        
        
    def _p(self,tt,dd,a,at,ad,att,atd,add, diff=0):
        """Calculate pressure from (T,d) (inner routine)
    p, pt, pd = _p(T, d, diff=0)
    
Accepts arguments:
T   Temperature array in Kelvin
d   Density array in kg/m3

Returns:
p   Pressure
pt  Derivative with respect to temperature
pd  Derivative with respect to density
"""
        p = 0.
        pt = 0.
        pd = 0.

        R = self._R()
        Tc = self.data['Tc']
        dc = self.data['dc']

        dd2 = dd*dd
        p = (dd2*ad/tt) * (dc*R*Tc)
        if diff>0:
            pt = dd2*(ad - tt*atd) * (dc*R)
            pd = dd*(2*ad + dd*add)/tt * (R*Tc)

        return p,pt,pd
        

    def _e(self,tt,dd,a,at,ad,att,atd,add, diff=0):
        """Internal Energy (inner routine)
    e,eT,ed = _e(tt,dd,a,at,ad,att,atd,add, diff=0)

Calculates internal energy and its derivatives.  This inner routine 
accepts dimensionless temperature, density, Helmholtz free energy, and
its derivatives, but returns values that are rescaled to units J, K, kg.

tt, dd
    Dimensionless temperature and dimensionless density: Tc/T, d/dc.
    
a, at, ad, att, atd, add
    Dimensionless Helmholtz free energy and its derivatives with respect
    to dimensionless temperature and density.  These are the values as
    returned by _fo() and _fr().  Written as "alpha" in the handbook.
    
diff    (default=0)
    The highest derivative to calculate.  Inner routines must be either
    0 or 1 -- 2 or higher is not accepted.  The value passed to _fo()
    and _fr() must be one greater than the value passed here.
    
Returns
e   [J/kg]      Internal energy
eT  [J/kg/K]    Differential with respect to temperature
ed  [J.m3/kg2]  Differential with respect to density
"""
        eT = None
        ed = None
        
        R = self._R()
        Tc = self.data['Tc']
        dc = self.data['dc']
        
        e = at*R*Tc
        if diff>0:
            eT = -tt*tt*att*R
            ed = atd*R*Tc/dc
        
        return e,eT,ed


    def _h(self,tt,dd,a,at,ad,att,atd,add, diff=0):
        """Enthalpy (inner routine)
    h,hT,hd = _e(tt,dd,a,at,ad,att,atd,add, diff=0)

Calculates enthalpy and its derivatives.  This inner routine 
accepts dimensionless temperature, density, Helmholtz free energy, and
its derivatives, but returns values that are rescaled to units J, K, kg.

tt, dd
    Dimensionless temperature and dimensionless density: Tc/T, d/dc.
    
a, at, ad, att, atd, add
    Dimensionless Helmholtz free energy and its derivatives with respect
    to dimensionless temperature and density.  These are the values as
    returned by _fo() and _fr().  Written as "alpha" in the handbook.
    
diff    (default=0)
    The highest derivative to calculate.  Inner routines must be either
    0 or 1 -- 2 or higher is not accepted.  The value passed to _fo()
    and _fr() must be one greater than the value passed here.
    
Returns
h   [J/kg]      Enthalpy
hT  [J/kg/K]    Differential with respect to temperature
hd  [J.m3/kg2]  Differential with respect to density
"""
        hT = None
        hd = None
        
        R = self._R()
        Tc = self.data['Tc']
        dc = self.data['dc']
        
        h = at + dd*ad/tt
        h *= R*Tc
        if diff>0:
            hT = dd*ad - tt*(tt*att + dd*atd)
            hT *= R
            hd = (ad + dd*add)/tt + atd
            hd *= R*Tc/dc
        
        return h,hT,hd

    def _s(self,tt,dd,f,ft,fd,ftt,ftd,fdd, diff=0):
        """Entropy (inner routine)
    s,sT,sd = _s(tt,dd,f,ft,fd,ftt,ftd,fdd, diff=0)

Calculates entropy and its derivatives.  This inner routine 
accepts dimensionless temperature, density, Helmholtz free energy, and
its derivatives, but returns values that are rescaled to units J, K, kg.

tt, dd
    Dimensionless temperature and dimensionless density: Tc/T, d/dc.
    
f, ft, fd, ftt, ftd, fdd
    Dimensionless Helmholtz free energy and its derivatives with respect
    to dimensionless temperature and density.  These are the values as
    returned by _fo() and _fr().
    
diff    (default=0)
    The highest derivative to calculate.  Inner routines must be either
    0 or 1 -- 2 or higher is not accepted.  The value passed to _fo()
    and _fr() must be one greater than the value passed here.
    
Returns
s   [J/kg/K]    Entropy
sT  [J/kg/K2]   Differential with respect to temperature
sd  [J.m3/K/kg2] Differential with respect to density
"""
        sT = None
        sd = None
        
        R = self._R()
        Tc = self.data['Tc']
        dc = self.data['dc']
        
        s = R*(tt*ft - f)
        if diff>0:
            sT = -tt*tt*tt*ftt*(R/Tc)
            sd = (tt*ftd - fd)*(R/dc)
        return s,sT,sd

    def _f(self,tt,dd,a,at,ad,att,atd,add, diff=0):
        """Free (Helmholtz) Energy (inner routine)
    f,fT,fd = _f(tt,dd,a,at,ad,att,atd,add, diff=0)

Calculates free energy and its derivatives.  This inner routine 
accepts dimensionless temperature, density, Helmholtz free energy, and
its derivatives, but returns values that are rescaled to units J, K, kg.

tt, dd
    Dimensionless temperature and dimensionless density: Tc/T, d/dc.
    
a, at, ad, att, atd, add
    Dimensionless Helmholtz free energy and its derivatives with respect
    to dimensionless temperature and density.  These are the values as
    returned by _fo() and _fr().  Written as "alpha" in the handbook.
    
diff    (default=0)
    The highest derivative to calculate.  Inner routines must be either
    0 or 1 -- 2 or higher is not accepted.  The value passed to _fo()
    and _fr() must be one greater than the value passed here.
    
Returns
f   [J/kg]      Helmholtz free energy
fT  [J/kg/K]    Differential with respect to temperature
fd  [J.m3/kg2]  Differential with respect to density
"""
        fT = None
        fd = None
        
        R = self._R()
        Tc = self.data['Tc']
        dc = self.data['dc']
        
        f = (a/tt) * (R*Tc)
        if diff:
            fT = (a - tt*at)*R
            fd = (ad/tt)*(R*Tc/dc)
        return f, fT, fd

    def _g(self,tt,dd,a,at,ad,att,atd,add, diff=0):
        """Free (Gibbs) Energy (inner routine)
    g,gT,gd = _g(tt,dd,a,at,ad,att,atd,add, diff=0)

Calculates Gibbs free energy and its derivatives.  This inner routine 
accepts dimensionless temperature, density, Helmholtz free energy, and
its derivatives, but returns values that are rescaled to units J, K, kg.

tt, dd
    Dimensionless temperature and dimensionless density: Tc/T, d/dc.
    
a, at, ad, att, atd, add
    Dimensionless Helmholtz free energy and its derivatives with respect
    to dimensionless temperature and density.  These are the values as
    returned by _fo() and _fr().  Written as "alpha" in the handbook.
    
diff    (default=0)
    The highest derivative to calculate.  Inner routines must be either
    0 or 1 -- 2 or higher is not accepted.  The value passed to _fo()
    and _fr() must be one greater than the value passed here.
    
Returns
g   [J/kg]      Gibbs free energy
gT  [J/kg/K]    Differential with respect to temperature
gd  [J.m3/kg2]  Differential with respect to density
"""
        gT = None
        gd = None
        
        R = self._R()
        Tc = self.data['Tc']
        dc = self.data['dc']

        g = (a + dd*ad)/tt
        g *= R*Tc
        if diff:
            gT = a + dd*ad - tt*at - tt*dd*atd
            gT *= R
            gd = (2*ad + dd*add)/tt
            gd *= (R*Tc/dc)
            
        return g,gT,gd
        
    def _a(self,tt,dd,a,at,ad,att,atd,add):
        """Speed of sound (inner routine)
    A = _a(tt,dd,a,at,ad,att,atd,add)

Calculates speed of sound.  This inner routine accepts dimensionless 
temperature, density, Helmholtz free energy, and its derivatives, but 
returns values that are rescaled to units J, K, kg.

tt, dd
    Dimensionless temperature and dimensionless density: Tc/T, d/dc.
    
a, at, ad, att, atd, add
    Dimensionless Helmholtz free energy and its derivatives with respect
    to dimensionless temperature and density.  These are the values as
    returned by _fo() and _fr().  Written as "alpha" in the handbook.
    
*NOTE*
This routine does not return the derivatives of speed of sound.  Always
pass diff=2 to _fo() and _fr() when using this property.

Returns
A   [m/s]      Wave speed
"""
        R = self._R()
        Tc = self.data['Tc']
        dc = self.data['dc']

        C = dd*ad
        B = C - tt*dd*atd
        A = 2*C + dd*dd*add - B*B/(tt*tt*att)
        return np.sqrt(R * Tc * A / tt)

        
    def _cp(self,tt,dd,a,at,ad,att,atd,add):
        """Constant-pressure specific heat (inner routine)
    cp = _cp(tt,dd,a,at,ad,att,atd,add)

Calculates constant-pressure specific heat.  This inner routine 
accepts dimensionless temperature, density, Helmholtz free energy, and
its derivatives, but returns values that are rescaled to units J, K, kg.

tt, dd
    Dimensionless temperature and dimensionless density: Tc/T, d/dc.
    
a, at, ad, att, atd, add
    Dimensionless Helmholtz free energy and its derivatives with respect
    to dimensionless temperature and density.  These are the values as
    returned by _fo() and _fr().  Written as "alpha" in the handbook.
    
*NOTE*
This routine does not return the derivatives of speed of sound.  Always
pass diff=2 to _fo() and _fr() when using this property.

Returns
cp  [J/kg/K]    Constant-pressure specific heat
"""
        R = self._R()
        Tc = self.data['Tc']
        dc = self.data['dc']

        C = dd*ad
        B = C - tt*atd
        cp = -tt*tt*att + B*B/(2*C + dd*dd*add)
        return R*cp
        
        
    def _cv(self,tt,dd,a,at,ad,att,atd,add):
        """Constant-volume specific heat (inner routine)
    cv = _cv(tt,dd,a,at,ad,att,atd,add)

Calculates constant-volume specific heat.  This inner routine 
accepts dimensionless temperature, density, Helmholtz free energy, and
its derivatives, but returns values that are rescaled to units J, K, kg.

tt, dd
    Dimensionless temperature and dimensionless density: Tc/T, d/dc.
    
a, at, ad, att, atd, add
    Dimensionless Helmholtz free energy and its derivatives with respect
    to dimensionless temperature and density.  These are the values as
    returned by _fo() and _fr().  Written as "alpha" in the handbook.
    
*NOTE*
This routine does not return the derivatives of specific heat.  Always
pass diff=2 to _fo() and _fr() when using this property.

Returns
cv  [J/kg/K]    Constant-volume specific heat
"""
        R = self._R()
        Tc = self.data['Tc']
        dc = self.data['dc']

        return -R * tt*tt*att

    #########################
    #                       #
    #  Argparse Algorithms  #
    #                       #
    #########################
        
    def _sat_argparse(self, T=None, p=None):
        """A standard argument parsing scheme for all user-layer saturation properties
    T,dL,dV = _sat_argparse(T=None, p=None)
    
Enforces that all returned parameters are numpy arrays with at least one
dimension.  Accepts T or p as scalars or array-like objects in 
[unit_temperature] and [unit_pressure] respectively.
    
Returns
T   the temperature in K
dL and dV are the liquid and vapor densities in kg/m3

** DESCRIPTION **
Calls the _Tsat() or _psat() method based on the information provided, 
while also asserting the correct units, array dimension, and out-of-
bounds checking.

The _argparse algorithm does not return pressure, because it is not 
needed in all property calculations, and it is not always needed as a
part of specifying the state, so there are many cases in which 
calculating it is wasted effort.  However, pressure is always a 
necessary calculation when iterating on the saturation state, so it is
returned to prevent a potentially redundant calculation.
"""
        if p is None:
            if T is None:
                T = pm.config.def_T()
            T = pm.units.temperature_scale(
                    np.asarray(T, dtype=float), 
                    to_units='K')
            if T.ndim==0:
                T = np.reshape(T, (1,))
            
            # Check for values out-of-bounds
            # For now, we'll also exclude the critical point, because 
            # the sat algorithms can't accept it.
            Ioob = np.logical_or(T < self.data['Tt'], T >= self.data['Tc'])
            if Ioob.any():
                # Initialize results
                dL = np.empty_like(T, dtype=float)
                dV = np.empty_like(T, dtype=float)
                # Check for precise equality with the critical point
                I = (T == self.data['Tc'])
                if I.any():
                    dL[I] = self.data['dc']
                    dV[I] = self.data['dc']
                # Check for points out of bounds and not precisely critical
                I = Ioob ^ I
                if I.any():
                    if pm.config['warning_verbose']:
                        pm.utility.print_warning('_mp2._sat_argparse: Saturation properties are not available beyond the triple or critical points.')
                    dL[I] = pm.config['def_oob']
                    dV[I] = pm.config['def_oob']
                # Then, work on the in-bounds values
                Ioob = np.logical_not(Ioob)
                _,dL[Ioob],dV[Ioob] = self._Tsat(T[Ioob])
            else:
                _,dL,dV = self._Tsat(T)
                
        elif T is None:
            p = pm.units.pressure(
                    np.asarray(p, dtype=float), 
                    to_units='Pa')
            if p.ndim==0:
                p = np.reshape(p, (1,))

            # Check for values out-of-bounds
            # For now, we'll also exclude the critical point, because 
            # the sat algorithms can't accept it.
            Ioob = np.logical_or(p < self._sattable['p'][0], p >= self._sattable['p'][-1])
            if Ioob.any():
                # Initialize results
                T = np.empty_like(p, dtype=float)
                dL = np.empty_like(p, dtype=float)
                dV = np.empty_like(p, dtype=float)
                # Check for precise equality with the critical point
                I = (p == self._sattable['p'][-1])
                if I.any():
                    T[I] = self.data['Tc']
                    dL[I] = self.data['dc']
                    dV[I] = self.data['dc']
                # Check for points out of bounds and not precisely critical
                I = Ioob ^ I
                if I.any():
                    if pm.config['warning_verbose']:
                        pm.utility.print_warning('_mp2._sat_argparse: Saturation properties are not available beyond the triple or critical points.')
                    T[I] = pm.config['def_oob']
                    dL[I] = pm.config['def_oob']
                    dV[I] = pm.config['def_oob']
                # Then, work on the in-bounds values
                Ioob = np.logical_not(Ioob)
                T[Ioob],dL[Ioob],dV[Ioob] = self._psat(p[Ioob])
            else:
                T,dL,dV = self._psat(p)
        else:
            raise pm.utility.PMParamError(
                '_sat_argparse: Saturation temperature and pressure cannot be simultaneously specified')

        return T, dL, dV
        
        
    def _argparse(self, *varg, debug=False, **kwarg):
        """Present a standard argument scheme for all user-layer property methods
    T,d1,d2,x,I = _argparse( .. property arguments ..)

This method processes the arguments passed to property methods, allowing
users to specify the thermodynamic state using flexible combinations of 
properties.

Below are the keyword arguments accepted and the corresponding :
    e   internal energy
    h   enthalpy
    s   entropy
    T   temperature
    p   pressure
    d   density
    v   specific volume
    x   quality

**SPECIFYING A STATE**
Users may specify the state in three ways:
(1) Positional arguments (with no keywords) are always interpreted as
    temperature and pressure.  For example:
        f(304.1, 1.9)
    is interpreted T=304.1 [unit_temperature] and p=1.9 [unit_pressure]
    If positional arguments are omitted, the missing argument is set to
    its default in pm.config:
        def_T
        def_T_unit
        def_p
        def_p_unit
    For example:
        f(304.1)
    is interpreted T=304.1 [unit_temperature] and p=[def_p] [def_p_unit]

(2) Keyword arguments are used to identify the properties being passed.
    For example:
        f(h=192.1, p=14.1)
    is interpreted as enthalpy and pressure in their respective units.
    Only two arguments may be specified this way.

(3) As a special case, (T,p,x) may be specified as a keyword argument 
    triple.  While other properties (except g) are capable of specifying
    a two-phase mixture, (T,p) cannot.  For example,
        f(T=[284., 285., 286.], p=1., x=[-1, 0.45, -1])
    specifies constant-pressure states on either side of saturation 
    with one two-phase mixture in the middle.
    **NOTE** 
    Technically, the state is over-specified if T,p,x are all specified 
    together.  Temperature is used to specify the saturation state, and 
    pressure IS NOT CHECKED for consistency.  User beware.

**ARRAYS**
Like in example 3 above, users may pass some or all properties as array-
like objects.  These are automatically built into Numpy arrays of the
appropriate shape, and all properties are automatically made compatible 
using Numpy's broadcasting rules.  If broadcasting fails, Numpy's 
back-end will raise an error.  

Example 1: An array of states
    T = [200., 300., 400.]
    s = [1.9, 2.0, 2.1]
    f(T=T, s=s)
Here, the two properties (temperature and entropy) are passed as lists
with the same dimensions.  No broadcasting is needed - this is 
interpreted as three states (T=200., s=1.9), (T=300., s=2.0), and so on.

Example 2: Broadcast arrays
    T = [[200.], [300.], [400.]]
    s = [1.9, 2.0, 2.1]    
    f(T=T, s=s)
Here, temperature has been modified to appear as a column vector.  
Because its values proceed in a different dimension (Numpy's "axis"),
this is interpreted as the nine states that result from the different 
combinations of T and s values.  See Numpy's ix_() or meshgrid() 
functions for automatically generating these kinds of arrays.

Example 3: Constant parameter
    T = 300
    s = [1.9, 2.0, 2.1]
    f(T=T, s=s)
This is interpreted as three states, all of which have temperature 300.
Note that temperature is lazily expressed as an integer - it will be
automatically promoted to a floating point.

**DISALLOWED COMBINATIONS**

Most property combinations are allowed, but some combinations are either
numerically unstable, or they do not theoretically define a unique 
state.
(1) Density and specific volume may not be specified together - they are
    redundant expressions of the same property.
(2) No two properties from the "energy" set may be specified together:
    {'T', 'e', 'h'}.  These either fail to describe a unique state 
    (meaning that there are multiple states that can be found with the 
    same values of these properties) or the state is very poorly
    defined (meaning that the resulting numerical inversion problem
    is very nearly singular).
(3) Quality may only be specified with temperature and/or pressure.  
    Specifying a entropy or an energy property (like enthalpy) with 
    quality does not define a unique state.  For example, there are 
    multiple states at which the same quality and enthalpy can be found.
(4) All combinations with Gibbs energy (g) and free energy (f) are 
    disallowed for specifying a state.  There is no combination with
    other properties that forms a well defined state over the entire 
    domain.

**BACK END**
The mp2 class back-end calculates properties exclusively in units
    unit_matter         kg
    unit_mass           kg
    unit_molar          kmol
    unit_temperature    K
    unit_pressure       Pa
    unit_length         m
    unit_volume         m3
    unit_energy         J
    unit_time           s

All inputs are coerced into floating point Numpy arrays of at least 
dimension 1, so broadcasting has the highest chance of success.  Unit
conversion of the outputs automatically in mp2's standard units (see
below).

Returns:
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
        # Always make sure the tables have been built
        self._build(force=False)
        
        # 1) Handle varg and kward and their defaults
        # 2) Apply the argument rules...
        #   2.1: All arguments must be legal
        #   2.2: There are only two arguments unless one is x
        #   2.3: x may only be specified with T or p
        #   2.4: Energy properties, T, e, h, f, and g may not be specified together
        #   2.5: d and v may not be specified together 
        #   
        # 3) Convert the arguments to arrays with dim 1 or greater
        # 4) Convert to standard units
        # 5) Check for out-of-bounds on basic arguments
        # 6) Replace specific volume with density if it appears
        # 7) Case out the possible combinations
        # Even though p and g are "standard" inverse properties, they
        # they are treated as special cases, because their values are
        # constant under the dome -- they behave differently in inverse
        # calculations.
        #   7.1: x is specified
        #       7.1.1: x,T,p
        #       7.1.2: x,T
        #       7.1.3: x,p
        #   7.2: T,?
        #       7.2.1: T,d      <== Easiest case -- already primary properties
        #       7.2.2: T,p      <== Special case because T,p is constant under the dome
        #       7.2.3: T + inverse
        #   7.3: d,?
        #       d + inverse
        #   7.4: p,?
        #       p + inverse     <== Special case 
        #   7.5: ?,?
        #       Any two remaining inverse
        # 
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
        inverse_methods = {'p':self._p, 'e':self._e, 'h':self._h, 's':self._s}
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
        # 2.4: T, e, and h may not be specified together
        if len(args.intersection({'T', 'e', 'h'})) > 1:
            raise pm.utility.PMParamError(
                    'Energy parameters, T, e, h, f, or g, may not be specified as a pair.')
        # 2.5: Density and specific volume cannot be specified together
        if 'v' in args and 'd' in args:
            raise pm.utility.PMParamError('Density (d) and specific volume (v) cannot be specified together.')
        # 2.6: p may not be specified with f or g
        if 'p' in args:
            if 'f' in args or 'g' in args:
                raise pm.utility.PMParamError('Pressure (p) may not be specified with Gibbs energy (g) nor with free energy (f).')
        
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
        if 'p' in kwarg:
            kwarg['p'] = pm.units.pressure(kwarg['p'], to_units='Pa')
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
            # Remove v --  it will be as if the user passed d instead
            args.remove('v')
            basic_args.remove('v')
            del kwarg['v']
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
                    d1 = np.empty_like(T, dtype=float)
                    d2 = np.empty_like(T, dtype=float)
                    # Check for out-of-bounds
                    TT = T[I]
                    Ioob = (TT < self.data['Tt']) + (self.data['Tc'] <= TT)
                    if Ioob.any():
                        if pm.config['warning_verbose']:
                            pm.utility.print_warning('mp2._argparse: Specified non-negative quality and temperatures out of [Tt,Tc].')
                        d1[I][Ioob] = pm.config['def_oob']
                        d2[I][Ioob] = pm.config['def_oob']
                        I[I][Ioob] = False
                        TT = T[I]
                    # Calculate densities for saturated states
                    if I.any():
                        _, d1[I], d2[I] = self._Tsat(TT)
                    # Calculate densities for non-saturated states
                    Ids = np.logical_not(I)
                    if Ids.any():
                        # Calculate densities for non-saturated points
                        d2[Ids],Isat,Ioob = self._Tmapsearch2(self._table['p'], T[Ids], p[Ids])
                        self._Titer(T, d2, self._p, p, Ids.copy())
                        d1[Ids] = d2[Ids]
                    return T, d1, d2, x, I
                # 7.1.2: T,x
                else:
                    T,x = np.broadcast_arrays(kwarg['T'], kwarg['x'])
                    d1 = np.empty_like(T, dtype=float)
                    d2 = np.empty_like(T, dtype=float)
                    if (x<0).any():
                        raise pm.utility.PMParamError(
                            'mp2._argparse(): Found x<0.  Only two-phase mixtures can be specified with T,x.  All values of x must be [0,1].')
                    # Detect out-of-bounds
                    I = (T < self.data['Tt']) + (self.data['Tc'] <= T)
                    if I.any():
                        if pm.config['warning_verbose']:
                            pm.utility.print_warning('mp2._argparse(): With (T,x) found temperatures below Tt or above Tc.')
                        d1[I] = pm.config['def_oob']
                        d2[I] = pm.config['def_oob']
                    I = np.logical_not(I)
                    _, d1[I], d2[I] = self._Tsat(T[I])
                    return T, d1, d2, x, I
            # 7.1.3: p,x
            else:
                p,x = np.broadcast_arrays(kwarg['p'], kwarg['x'])
                d1 = np.empty_like(p, dtype=float)
                d2 = np.empty_like(p, dtype=float)
                T = np.empty_like(p, dtype=float)
                if (x<0).any():
                    raise pm.utility.PMParamError(
                        'Found x<0.  Only two-phase mixtures can be specified with p,x.  All values of x must be [0,1].')
                # Detect out-of-bounds
                I = (p < self.data['pt']) + (self._sattable['p'][-1] <= p)
                if I.any():
                    pm.utility.print_warning('mp2._argparse: With (p,x) found pressures beyond the triple or critical points.')
                    d1[I] = pm.config['def_oob']
                    d2[I] = pm.config['def_oob']
                    T[I] = pm.config['def_oob']
                I = np.logical_not(I)
                T[I], d1[I], d2[I] = self._psat(p[I])
                return T, d1, d2, x, I
        # 7.2: T,?
        elif 'T' in kwarg:
            # 7.2.1: T,d
            if 'd' in kwarg:
                # broadcast the arrays
                T,d = np.broadcast_arrays(kwarg['T'],kwarg['d'])
                x = np.full_like(T, -1.)
                # By default, d1 and d2 are merely pointers to d
                # If there are 2-phase points, this behavior will be overridden
                d1 = d
                d2 = d
                # Identify sub-critical temperatures
                I = (T < self.data['Tc'])
                if I.any():
                    _,dL,dV = self._Tsat(T[I])
                    # Down-select to the densities that are under the dome
                    dd = d[I]
                    # Of the down-selected states, which are actually 2-phase?
                    Isat = np.logical_and(dV < dd, dd < dL)
                    I[I] = Isat
                    # If there are any two-phase mixture points
                    if Isat.any():
                        # Make copies of the d array
                        d1 = d.copy()
                        d2 = d.copy()
                        # Down-select the vapor, liquid, and mixture densities
                        dV = dV[Isat]
                        dL = dL[Isat]
                        dd = dd[Isat]
                        d1[I] = dL
                        d2[I] = dV
                        # Calculate liquid volume
                        dL = 1./dL
                        # Calculate quality
                        x[I] = (1./dd - dL) / (1./dV - dL)
                return T,d1,d2,x,I
            # 7.2.2: T,p
            elif 'p' in kwarg:
                T,p = np.broadcast_arrays(kwarg['T'],kwarg['p'])
                x = np.full_like(T, -1.)
                # T,p cannot be used to specify a saturated state - no need to check
                d2, _, I = self._Tmapsearch2(self._table['p'], T, p)
                # Iterate only on states that are in-bounds
                I = np.logical_not(I)
                self._Titer(T, d2, self._p, p, I)
                # All I values should now be False
                d1 = np.copy(d2)
                return T,d1,d2,x,I
            # 7.2.3: T + inverse
            else:
                # Isolate the inverse property and its method
                args.remove('T')
                fstr = args.pop()
                fn = inverse_methods[fstr]
                # Broadcast the arrays
                T,fvalue = np.broadcast_arrays(kwarg['T'], kwarg[fstr])
                # Search the table for a density to match
                d2, Isat, Ioob = self._Tmapsearch2(self._table[fstr], T, fvalue)
                # Initialize quality and d1
                x = np.full_like(T, -1.)
                d1 = np.empty_like(d2, dtype=float)
                # Deal with states that are saturated or nearly saturated
                if Isat.any():
                    TT = T[Isat]
                    # Calculate the saturation densities
                    _, dL, dV = self._Tsat(TT)
                    d1[Isat] = dL
                    d2[Isat] = dV
                    # Calculate the inverse property's saturation properties
                    fL,_,_ = fn(*self._ff(T=TT,d=dV,diff=1))
                    fV,_,_ = fn(*self._ff(T=TT,d=dL,diff=1))
                    # Deduce quality from fvalue
                    x[Isat] = (fvalue[Isat] - fL)/(fV - fL)
                    # Some of these will be points that are merely near
                    # the dome and not actually under it.  
                    Ids = np.zeros_like(T, dtype=bool)
                    # If out on the liquid side, use liquid density
                    Ids[Isat] = x[Isat]<0
                    if Ids.any():
                        x[Ids] = -1
                        Isat[Ids] = False
                        d2[Ids] = d1[Ids]
                    # If out on the vapor side, use vapor density
                    Ids[Isat] = x[Isat]>1
                    if Ids.any():
                        x[Ids] = -1
                        Isat[Ids] = False
                        d1[Ids] = d2[Ids]
                # On all non-saturated points, iterate
                Ids = np.logical_not(Isat)
                # Deal with any out-of-bounds points
                if Ioob.any():
                    if pm.config['warning_verbose']:
                        pm.utility.print_warning('mp2._argparse(): Found T,? property combinations that were out-of-bounds.')
                    Ids[Ioob] = False
                    # d2 will already be set by the mapsearch algorithm
                    d1[Ioob] = pm.config['def_oob']
                    # Leave temperature as specified
                self._Titer(T, d2, fn, fvalue, Ids.copy())
                d1[Ids] = d2[Ids]
                return T,d1,d2,x,Isat
        # 7.3: d + inverse
        elif 'd' in kwarg:
            args.remove('d')
            fstr = args.pop()
            fn = inverse_methods[fstr]
            # Broadcast to the appropriate dimensions
            d,fvalue = np.broadcast_arrays(kwarg['d'], kwarg[fstr])
            # Initialize d2, d1, and x
            # For now, d, d2, and d1 are separate, because some d values
            # can represent two-phase mixtures.  We'll dole out the d
            # values appropriately once we know which are under the dome
            d2 = np.empty_like(d, dtype=float)
            d1 = np.empty_like(d, dtype=float)
            x = np.full_like(d, -1.)
            # Identify estimates for T
            # Use entropy extrapolation if property is s
            if fstr == 's':
                T, Isat, Ioob = self._dmapsearch2(self._table[fstr], d, fvalue, zde=1)
            # Otherwise, keep normal extrapolation
            else:
                T,Isat,Ioob = self._dmapsearch2(self._table[fstr], d, fvalue)
            # Investigate states that may be saturated
            if Isat.any():
                # Calculate saturated densities at our best guess for T
                _, d1[Isat], d2[Isat] = self._Tsat(T[Isat])
                self._dsatiter2(T, d1, d2, d, fn, fvalue, Isat.copy())
                # Calculate quality
                xx = (d1[Isat]/d[Isat] - 1)/(d1[Isat]/d2[Isat] - 1)
                # Points that were merely very close to saturated will
                # converge with quality out of bounds
                Ids = np.logical_or(xx<0, xx>1)
                if Ids.any():
                    xx[Ids] = -1
                    x[Isat] = xx
                    Isat[Isat] = np.logical_not(Ids)
                else:
                    x[Isat] = xx
                
            # Down-select only points that are not saturated
            Ids = np.logical_not(Isat)
            # Check for out-of-bounds states
            if Ioob.any():
                # Remove out-of-bounds points from iteration
                Ids[Ioob] = False
                # Temperature will already be set by the mapsearch
                # Leave density as-specified.
                d1[Ioob] = d[Ioob]
                d2[Ioob] = d[Ioob]
            
            self._diter(T, d, fn, fvalue, Ids.copy())
            d1[Ids] = d[Ids]
            d2[Ids] = d[Ids]
            return T, d1, d2, x, Isat
        # At this stage, there are two inverse properties
        # 7.4: p + inverse
        elif 'p' in kwarg:
            # Isolate the other inverse property
            args.remove('p')
            fstr = args.pop()
            fn = inverse_methods[fstr]
            # Broadcast to the appropriate dimensions
            p,fvalue = np.broadcast_arrays(kwarg['p'], kwarg[fstr])
            # Initialize results
            d1 = np.empty_like(p, dtype=float)
            x = np.full_like(p, -1.)
            # Find an initial guess for the state
            zde1 = 1 if fstr == 's' else 0
            T, d2, Isat, Ioob = self._mapsearch2(self._table['p'], self._table[fstr], p, fvalue, zde1=zde1)
            if Isat.any():
                # Establish the saturation states
                TT,dL,dV = self._psat(p[Isat])
                # Get the saturation values
                argL = self._ff(T=TT, d=dL, diff=1)
                argV = self._ff(T=TT, d=dV, diff=1)
                fL = fn(*argL, diff=0)[0]
                fV = fn(*argV, diff=0)[0]
                # Deduce quality from fvalue
                x[Isat] = (fvalue[Isat] - fL)/(fV - fL)
                # Some of these will be points that are merely near
                # the dome and not actually under it.  
                Ids = np.zeros_like(T, dtype=bool)
                # If out on the liquid side, use liquid density for the iteration
                Ids[Isat] = x[Isat]<0
                if Ids.any():
                    x[Ids] = -1
                    Isat[Ids] = False
                    d2[Ids] = d1[Ids]
                # If out on the vapor side, use vapor density for the iteration
                Ids[Isat] = x[Isat]>1
                if Ids.any():
                    x[Ids] = -1
                    Isat[Ids] = False
                    d1[Ids] = d2[Ids]
            # On all non-saturated points, iterate
            Ids = np.logical_not(Isat)
            # Deal with any out-of-bounds points
            if Ioob.any():
                if pm.config['warning_verbose']:
                    pm.utility.print_warning('mp2._argparse(): Found p,? property combinations that were out-of-bounds.')
                Ids[Ioob] = False
                # d2 will already be set by the mapsearch algorithm
                d1[Ioob] = pm.config['def_oob']
                # Leave temperature as specified
            self._iter2(T, d2, self._p, fn, p, fvalue, Ids.copy(), debug=True)
            d1[Ids] = d2[Ids]
            return T,d1,d2,x,Isat
        # 7.5: Two inverse properties
        else:
            # Isolate the property strings, their methods, and their value arrays
            f0str = args.pop()
            f1str = args.pop()
            fn0 = inverse_methods[f0str]
            fn1 = inverse_methods[f1str]
            f0value, f1value = np.broadcast_arrays(kwarg[f0str], kwarg[f1str])
            # Look up estimates for T and d in the property tables
            # Detect whether zero-density extrapolation is needed
            zde0 = 1 if f0str == 's' else 0
            zde1 = 1 if f1str == 's' else 0
            T,d2,Isat,Ioob = self._mapsearch2(self._table[f0str], self._table[f1str], f0value, f1value, zde0=zde0, zde1=zde1)
            print('mapsearch2:',T,d2,Isat,Ioob)
            x = np.full_like(T, -1.)
            d1 = np.empty_like(d2, dtype=float)
            if Isat.any():
                # Stash a copy of the original T-values so we can recover from failed iteration
                # Obtain estimates for saturation densities using our best guess for T
                _, d1[Isat], d2[Isat] = self._Tsat(T[Isat])
                self._satiter2(T, d1, d2, fn0, fn1, f0value, f1value, Isat.copy(), debug=True)
                # Calculate quality from the converged values
                argL = self._ff(T=T[Isat], d=d1[Isat], diff=1)
                argV = self._ff(T=T[Isat], d=d2[Isat], diff=1)
                f0L = fn0(*argL, diff=0)[0]
                f0V = fn0(*argV, diff=0)[0]
                x[Isat] = (f0value[Isat] - f0L)/(f0V - f0L)
                # Detect states that are not quite under the dome, but very
                # close.  These will have converged to out-of-bounds values
                # for x.
                Ids = np.zeros_like(Isat, dtype=bool)
                # On the vapor edge
                Ids[Isat] = x[Isat] > 1
                if Ids.any():
                    # Use the vapor density
                    d1[Ids] = d2[Ids]
                    x[Ids] = -1
                    Isat[Ids] = False
                # On the liquid edge
                Ids[Isat] = x[Isat] < 0
                if Ids.any():
                    # Use the liquid density
                    d2[Ids] = d1[Ids]
                    x[Ids] = -1
                    Isat[Ids] = False
            # Transition to working on non-saturated states
            Ids = np.logical_not(Isat)
            
            # Deal with out-of-bounds
            if Ioob.any():
                T[Ioob] = pm.config['def_oob']
                d1[Ioob] = pm.config['def_oob']
                d2[Ioob] = pm.config['def_oob']
                x[Ioob] = -1
                Ids[Ioob] = False
            # Finally, iterate on any non-saturated points
            if Ids.any():
                self._iter2(T, d2, fn0, fn1, f0value, f1value, Ids.copy(), debug=True)
                print('iter2:', T, d2)
                d1[Ids] = d2[Ids]
                
            return T,d1,d2,x,Isat
            
        message = 'Please report a bug: Unhandled event [MASTER] in mp2._argparse with args:'
        prefix = ' '
        for name in args:
            message += prefix + name
            prefix = ', '
        raise pm.utility.PMParamError(message)

    ########################
    #                      #
    #  User-Layer Methods  #
    #                      #
    ########################
    
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
        R = pm.units.energy(self._R(), from_units = 'J')
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
        pc = self._sattable['p'][-1]
        if density:
            return  pm.units.temperature_scale( \
                        self.data['Tc'], from_units='K'),\
                    pm.units.pressure( \
                        pc, from_units='Pa'), \
                    pm.units.volume(\
                        pm.units.matter( \
                            self.data['dc'], \
                            self.data['mw'], \
                            from_units='kg'),\
                        from_units='m3', exponent=-1)
                    
        return  pm.units.temperature_scale( \
                    self.data['Tc'], from_units='K'),\
                pm.units.pressure( \
                    pc, from_units='Pa')
        
        
    def triple(self):
        """Triple point
    Tt, pt = triple()
    
Returns the triple temperature and pressure in a tuple pair in
[unit_temperature], [unit_pressure]
"""
        Tt = self.data.get('Tt')
        if Tt is None or Tt < self.data['Tlim'][0]:
            raise pm.utility.PMParamError('mp2.triple: This dataset does not include the triple point.')
        # Always make sure the tables have been built
        self._build(force=False)
        pt = self._sattable['p'][0]
        return  pm.units.temperature_scale( \
                    Tt, from_units='K'),\
                pm.units.pressure( \
                    pt, from_units='Pa')
        
    #                       #
    # Saturaiton properties #
    #                       #
    
    def ps(self, *varg, **kwarg):
        """Saturation pressure
    psat = ps(T)
        OR
    psat = ps(p=p)
    
Saturation line properties accept either T or p as keyword arguments.  

The optional diff keyword argument is 0 by default.  When set to 1 or
True, the temperature derivative of the saturation 

Returns:
psat    The saturaiton pressure in [unit_pressure]
pT      The derivative with respect to temperature in units
        [unit_pressure / unit_temperature]
"""
        T,_,dV = self._sat_argparse(*varg, **kwarg)
        arg = self._ff(T,dV,diff=1)
        p,_,_ = self._p(*arg, diff=0)
        pm.units.pressure(p, from_units='Pa', inplace=True)
        return p
        
        
    def Ts(self, *varg, **kwarg):
        """Saturation temperature
    Tsat = Ts(p)
    
Calculates the saturation temperature in terms of the pressure.  

Unlike the other saturation properties, ps() and Ts() only accept one
argument and only return one value - each calculates the one in terms
of the other.
"""
        T,_,_ = self._sat_argparse(*varg, **kwarg)
        pm.units.temperature_scale(T, from_units='K', inplace=True)
        return T
        
        
    def ds(self, *varg, **kwarg):
        """Saturation density
    dsL, dsV = ds(...)
    
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
    vsL, vsV = vs(...)
    
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
    esL, esV = es(...)

If no keyword is specified, saturation properties interpret the argument
as temperature.  However, pressure can be specified as well

    esL, esV = es(p=pvalue)
    
Returns the liquid (esL) and vapor (esV) saturation internal energy in
units [unit_energy / unit_matter]
"""
        T,dL,dV = self._sat_argparse(*varg, **kwarg)
        arg = self._ff(T,dL,diff=1)
        esL = self._e(*arg,diff=0)[0]
        
        arg = self._ff(T,dV,diff=1)
        esV = self._e(*arg,diff=0)[0]
        
        # Get a conversion factor
        conv = pm.units.energy(1., from_units='J')
        conv = pm.units.matter(conv, self.data['mw'],
                from_units='kg', exponent=-1)
        esL *= conv
        esV *= conv
        return esL, esV


    def hs(self, *varg, **kwarg):
        """Saturation enthalpy
    hsL, hsV = hs(...)
    
If no keyword is specified, saturation properties interpret the argument
as temperature.  However, pressure can be specified as well

    hsL, hsV = hs(p=pvalue)
    
Returns the liquid (hsL) and vapor (hsV) saturation enthalpy in
units [unit_energy / unit_matter]
"""
        T,dL,dV = self._sat_argparse(*varg, **kwarg)
        arg = self._ff(T,dL,diff=1)
        hsL = self._h(*arg,diff=0)[0]
        
        arg = self._ff(T,dV,diff=1)
        hsV = self._h(*arg,diff=0)[0]
        
        # Get a conversion factor
        conv = pm.units.energy(1., from_units='J')
        conv = pm.units.matter(conv, self.data['mw'],
                from_units='kg', exponent=-1)
        hsL *= conv
        hsV *= conv
        return hsL, hsV
        
        
    def ss(self, *varg, **kwarg):
        """Saturation entropy
    ssL, ssV = ss(...)
    
If no keyword is specified, saturation properties interpret the argument
as temperature.  However, pressure can be specified as well

    ssL, ssV = ss(p=pvalue)
    
Returns the liquid (ssL) and vapor (ssV) saturation entropy in
units [unit_energy / unit_matter / unit_temperature]
"""
        T,dL,dV = self._sat_argparse(*varg, **kwarg)
        arg = self._ff(T,dL,diff=1)
        ssL = self._s(*arg,diff=0)[0]
        
        arg = self._ff(T,dV,diff=1)
        ssV = self._s(*arg,diff=0)[0]
        
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
    p = p(...)
        OR
    p,x = p(..., quality=True)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns pressure in [unit_pressure] found in pm.config

If the optional ``quality'' keyword is set to True, the quality is also
returned to save a redundant call to x().

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        # Use d2.  In theory, p(d1) = p(d2), but the liquid is so stiff
        # that small numerical errors cause huge pressure errors
        # The problem is solved when the vapor density is used instead.
        # In all other conditions d1=d2
        arg = self._ff(T,d2,diff=1)
        p = self._p(*arg,diff=0)[0]
        
        p = pm.units.pressure(p, from_units='Pa')
        
        if quality:
            return p,x
        return p
        
        
    def d(self, *varg, quality=False, **kwarg):
        """Density
    d = d(...)
        OR
    d,x = d(..., quality=True)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns density in [unit_matter / unit_volume] found in pm.config

If the optional ``quality'' keyword is set to True, the quality is also
returned to save a redundant call to x().

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        if I.any():
            xx = x[I]
            dL = d1[I]
            dV = d2[I]
            d1[I] = (1.-xx)/dL
            d1[I] += xx/dV
            d1[I] = 1. / d1[I]
            
        d1 = pm.units.matter(d1, self.data['mw'], from_units='kg')
        d1 = pm.units.volume(d1, from_units='m3', exponent=-1)
        if quality:
            return d1,x
        return d1
        
        
    def v(self, *varg, quality=False, **kwarg):
        """Specific volume
    v = v(...)
        OR
    v,x = v(..., quality=True)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns specific volume in [unit_volume / unit_matter] found in 
pm.config

If the optional ``quality'' keyword is set to True, the quality is also
returned to save a redundant call to x().

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
        d,x = self.d(*varg, quality=True, **kwarg)
        if quality:
            return 1./d, x
        return 1./d
        
    def T(self, *varg, quality=False, **kwarg):
        """Temperature
    T = T(...)
        OR
    T,x = T(..., quality=True)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns temperature in [unit_temperature] found in pm.config

If the optional ``quality'' keyword is set to True, the quality is also
returned to save a redundant call to x().

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
        T,_,_,x,_ = self._argparse(*varg, **kwarg)
        T = pm.units.temperature_scale(T, from_units='K')
        if quality:
            return T,x
        return T
        
    def x(self, *varg, **kwarg):
        """Quality
    x = x(...)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns temperature in [unit_temperature] found in pm.config

Quality is usually needed in conjunction with other properties, and it 
is calculated as an intermediate parameter in the back-end anyway.  See
the optional ``quality'' keyword of the other property methods to save
a redundant call.

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
        _,_,_,x,_ = self._argparse(*varg, **kwarg)
        return x
        
    #                    #
    # Property functions #
    #                    #
    
    def state(self, *varg, **kwarg):
        """The state method calculates most available properties at once.
        
    sd = state(...)
    
Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.
    
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
    e = e(...)
        OR
    e,x = e(..., quality=True)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns internal energy in [unit_energy / unit_matter] found in 
pm.config.

If the optional ``quality'' keyword is set to True, the quality is also
returned to save a redundant call to x().

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        arg = self._ff(T,d2,diff=1)
        e = self._e(*arg,diff=0)[0]
        # If there are points under the dome
        if I.any():
            xx = x[I]
            e[I] *= xx
            arg = self._ff(T[I], d1[I], diff=1)
            e[I] += self._e(*arg, diff=0)[0] * (1-xx)
        # Convert the units back to user space
        pm.units.energy(e, from_units='J', inplace=True)
        pm.units.matter(e, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        if quality:
            return e,x
        return e
        
    def f(self, *varg, quality=False, **kwarg):
        """Free (Helmholtz) energy
    f = f(...)
        OR
    f,x = f(..., quality=True)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns free energy in [unit_energy / unit_matter] found in pm.config.

If the optional ``quality'' keyword is set to True, the quality is also
returned to save a redundant call to x().

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        arg = self._ff(T,d2,diff=1)
        f = self._f(*arg,diff=0)[0]
        # If there are points under the dome
        if I.any():
            xx = x[I]
            f[I] *= xx
            arg = self._ff(T[I], d1[I], diff=1)
            f[I] += self._f(*arg, diff=0)[0] * (1-xx)
        # Convert the units back to user space
        pm.units.energy(f, from_units='J', inplace=True)
        pm.units.matter(f, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        if quality:
            return f,x
        return f

    def g(self, *varg, quality=False, **kwarg):
        """Gibbs energy
    g = g(...)
        OR
    g,x = g(..., quality=True)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns Gibbs energy in [unit_energy / unit_matter] found in pm.config

If the optional ``quality'' keyword is set to True, the quality is also
returned to save a redundant call to x().

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        arg = self._ff(T,d2,diff=1)
        g = self._g(*arg,diff=0)[0]
        # Ignore points under the dome -- Gibbs energy is constant
        # Convert the units back to user space
        pm.units.energy(g, from_units='J', inplace=True)
        pm.units.matter(g, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        if quality:
            return g,x
        return g    
        
    def h(self, *varg, quality=False, **kwarg):
        """Temperature
    h = h(...)
        OR
    h,x = h(..., quality=True)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns enthalpy in [unit_energy / unit_matter] found in pm.config

If the optional ``quality'' keyword is set to True, the quality is also
returned to save a redundant call to x().

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        arg = self._ff(T,d2,diff=1)
        h = self._h(*arg,diff=0)[0]
        # If there are points under the dome
        if I.any():
            xx = x[I]
            h[I] *= xx
            arg = self._ff(T[I], d1[I], diff=1)
            h[I] += self._h(*arg, diff=0)[0] * (1-xx)
        # Convert the units back to user space
        pm.units.energy(h, from_units='J', inplace=True)
        pm.units.matter(h, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        if quality:
            return h,x
        return h


    def s(self, *varg, quality=False, **kwarg):
        """Entropy
    s = s(...)
        OR
    s,x = s(..., quality=True)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns entropy in [unit_energy / unit_matter / unit_temperature] found 
in pm.config.

**NOTE**
The entropy is calculated assuming the two-phase mixture is stratified -
the vapor is fully separated above the liquid.  The entropy of an finely
mixed liquid-vapor (e.g. a mist or cloud) is higher.

If the optional ``quality'' keyword is set to True, the quality is also
returned to save a redundant call to x().

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
            
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        arg = self._ff(T,d2,diff=1)
        s = self._s(*arg,diff=0)[0]
        # If there are points under the dome
        if I.any():
            xx = x[I]
            s[I] *= xx
            arg = self._ff(T[I], d1[I], diff=1)
            s[I] += self._s(*arg, diff=0)[0] * (1-xx)
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
    a = a(...)
        OR
    a,x = a(..., quality=True)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns speed of sound in [unit_length / unit_time] found in pm.config

**NOTE**
The speed of sound in a two-phase mixture is not well defined.  If the
mixture is stratified (fully separated vapor over liquid) each volume
has its own speed of sound determined by the saturation properties.
If the mixture is a finely-mixed mist, bubble field, or something in
the middle, wave propagation is a far more complex business.

As a result, a() returns the out-of-bounds value config['def_oob'] in
saturated states.  If the vapor and liquid speeds of sound are needed,
use the saturation densities to obtain their values there.

If the optional ``quality'' keyword is set to True, the quality is also
returned to save a redundant call to x().

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
        
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        arg = self._ff(T,d2,diff=2)
        a = self._a(*arg)
        # Speed of sound is not well defined under the dome.
        if I.any():
            a[I] = pm.config['def_oob']
        # Convert the units back to user space
        pm.units.length(a, from_units='m', inplace=True)
        pm.units.time(a, from_units='s', inplace=True, exponent=-1)
        if quality:
            return a,x
        return a
        

    def cp(self, *varg, quality=False, **kwarg):
        """Constant-pressure specific heat
    cp = cp(...)
        OR
    cp,x = cp(..., quality=True)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns specific heat in [unit_energy / unit_matter / unit_temperature] 
found in pm.config.

**NOTE**
Constant-pressure specific heat is infinite in saturated mixtures of any
kind.  This can cause some unexpected numerical problems if codes are 
not expecting these values.

If the optional ``quality'' keyword is set to True, the quality is also
returned to save a redundant call to x().

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
            
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        arg = self._ff(T,d2,diff=2)
        cp = self._cp(*arg)
        # Constant-pressure specific heat is infinite under the dome
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
    cv = cv(...)
        OR
    cv,x = cv(..., quality=True)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns specific heat in [unit_energy / unit_matter / unit_temperature] 
found in pm.config.

**NOTE**
The state() method returns out-of-bounds for constant-volume specific
heats of two-phase mixtures, but the cv() method calculates them 
correctly.  This is a deliberate design decision for speed in the 
state() method.  If cv() of two-phase mixtures is needed, users must use
cv() explicitly.

If the optional ``quality'' keyword is set to True, the quality is also
returned to save a redundant call to x().

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
        
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        arg = self._ff(T,d2,diff=2)
        cv = self._cv(*arg)
        # Constant-volume specific heat is a bit complicated under the dome
        if I.any():
            argL = self._ff(T[I], dL[I], diff=2)
            argV = self._ff(T[I], dV[I], diff=2)
            xx = x[I]
            # We'll need to calculate the derivative of quality with
            # respect to temperature.  To do that, we'll differentiate
            # the Maxwell criteria
            #   g(T,dL) = g(T,dV)
            #   p(T,dL) = p(T,dV)
            # Leads to
            #   (gLt-gVt)*dT = gLd*ddV - gVd*ddL
            #   (pLt-pVt)*dT = pLd*ddV - pVd*ddL
            # Matrix inversion gives ddV/dT and ddL/dT
            _,gLt,gLd = self._g(*argL,1)
            _,gVt,gVd = self._g(*argV,1)
            _,pLt,pLd = self._p(*argL,1)
            _,pVt,pVd = self._p(*argV,1)
            # This is only a 2x2, so we can do it "manually"
            temp = (gLd*pVd - pLd*gVd)
            gt = gLt - gVt
            pt = pLt - pVt
            dLT = (-pVd*gt + gVd*pt)/temp
            dVT = (-pLd*gt + gLd*pt)/temp

            # How does x change with temperature?  The process is 
            # constant volume, so the density is also constant.  Only
            # the saturation densities change.
            #     (dL/d ) - 1
            # x = -----------
            #     (dL/dV) - 1
            temp = dL/dV
            xT = (dLT * (1-xx)/dL + dVT * xx*temp/dV) / (temp-1)
            # Grab the saturation sensitivities
            eL,eLT,eLd = self._e(*argL,diff=1)
            eV,eVT,eVd = self._e(*argV,diff=1)
            # Calculate the true isochoric specific heat for the
            # two-phase mixture
            cv[I] = (eLT+eLd*dLT)*(1-xx) + (eVT+eVd*dVT)*xx + (eV-eL)*xT
            
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
    gam = gam(...)
        OR
    gam,x = gam(..., quality=True)

Query the _argparse() method's documentation for a detailed description
of the standard interface for specifying state.

Returns the dimensionless specific heat ratio.

**NOTE**
Constant-pressure specific heat is infinite in two-phase mixtures, and
so is gamma.  This can cause problems in codes that do not expect this
result.

If the optional ``quality'' keyword is set to True, the quality is also
returned to save a redundant call to x().

See also:
    a, cp, cv, d, e, f, g, gam, h, mw, p, R, s, T, v, x, state
"""
            
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        arg = self._ff(T,d2,diff=2)
        cp = self._cp(*arg)
        cv = self._cv(*arg)
        if I.any():
            cp[I] = np.inf
        
        if quality:
            return cp/cv, x
        return cp/cv

