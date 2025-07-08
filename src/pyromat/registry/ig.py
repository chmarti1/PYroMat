import pyromat as pm
import numpy as np
import os,sys
import re

#
# ig works in internal units
#   Energy: kJ
#   Matter: kmol
#   Pressure: Pa
#


class ig(pm.reg.__basedata__):
    """Ideal gas class using the Shomate equation of state.

** Available Properties **
IG has property methods:
  T()  temperature      (unit_temperature)
  p()  pressure         (unit_pressure)
  d()  density          (unit_matter / unit_volume)
  v()  specific volume  (unit_volume / unit_matter)
  cp() spec. heat       (unit_energy / unit_temperature / unit_matter)
  cv() spec. heat       (unit_energy / unit_temperature / unit_matter)
  gam()  spec. heat ratio (dless)
  e()  internal energy  (unit_energy / unit_matter)
  f()  free energy      (unit_energy / unit_matter)
  g()  Gibbs energy     (unit_energy / unit_matter)
  h()  enthalpy         (unit_energy / unit_matter)
  s()  entropy          (unit_energy / unit_temperature / unit_matter)
  state()       Calculates many properties!

These accept any of the following keyword arguments: T, p, d, v, e, h, s
  h(T=452.)
  h(d=1.3, p=3.4)
  h(s=6.2, p=3.4)

In the back end, most properties require only temperature except entropy,
which requires temperature and pressure.  When other properties are given,
temperature is needed, so additional calculations are necessary.  In the 
worst case, specifying entropy, enthalpy, or internal energy even requires
iteration.  For performance, once temperature and pressure become available,
they should always be used.  For enthalpy, pressure is not needed, so it 
may be omitted with no ill effect.  

** Other Properties **
There are some properties that do not depend on the state; they only need
to be converted to the relevant unit system.
  mw()      molecular weight (unit_mass / unit_molar)
  R()       gas constant
  Tlim()    a two-element array with the min,max temperatures 
            supported by the data set.
  atoms()   returns a dictionary with a key entry for each atom in
            the chemical formula and the corresponding integer 
            quantity of each.
              
For more information on any of these methods, access the in-line 
documentation using Python's built-in "help()" function.
"""

    def __init__(self,*arg,**kwarg):
        super(self.__class__,self).__init__(*arg,**kwarg)

        # Important constants
        self._pref_pa = 1e5
        # Initialize the species contents dictionary
        self._contents = None


    def _argparse(self, *varg, **kwarg):
        """Parse the arguments supplied to an IG property method
    T,p,d = _argparse(*varg, **kwarg)

_ARGPARSE automatically applies the default temperature and pressure,
def_T or def_p, from the pyromat.config system to deal with unspecified
parameters.  All inputs are re-cast as numpy arrays of at least one 
dimension and inputs are automatically converted from the configured 
user units into kJ, kmol, m^3.

The returned variables are arrays of temperature, T, pressure, p, and 
the density, d.  Temperature will always be returned, and at least one
of p and d will be populated as well, but one of them may be None.  
_argparse decides which to populate based on what is most efficient.
"""
        # 1) Handle varg and kward and their defaults
        # 2) Apply the argument rules...
        #   2.1: All arguments must be legal
        #   2.2: Only 2 arguments
        #   2.3: T, e, and h may not be specified together
        #   2.4: d and v may not be specified together 
        # 3) Convert the arguments to arrays with dim 1 or greater
        # 4) Convert to standard units
        # 5) Process equivalent arguments
        #   5.1: replace v with d
        #   5.2: replace e or h with T
        # 6) Check for out-of-bounds on basic arguments (This was shifted to last)
        # 7) Case out the possible combinations
        # 8) Broadcast the arrays appropriately
        # 9) Calculate T,p,d
        

        # Fancy tool for tracking iteration issues
        debug = False

        # 1) Handle varg and kwarg and apply defaults

        # If varg is specified, assign its values to T,p
        if len(varg) > 0:
            if 'T' in kwarg:
                raise pm.utility.PMParamError('T was specified both positionally and with a keyword.')
            kwarg['T'] = varg[0]
        if len(varg) > 1:
            if 'p' in kwarg:
                raise pm.utility.PMParamError('p was specified both positionally and with a keyword.')
            kwarg['p'] = varg[1]
        if len(varg) > 2:
            raise pm.utility.PMParamError('There are only two positional arguments: T, p.')

        # Count the number of arguments
        nargs = len(kwarg)
        if nargs == 1:
            if 'T' not in kwarg:
                kwarg['T'] = pm.config.def_T()
            else:
                kwarg['p'] = pm.config.def_p()
        elif nargs == 0:
            kwarg['T'] = pm.config.def_T()
            kwarg['p'] = pm.config.def_p()
        
        # 2) Apply the argument rules
        # Re-measure the number of arguments and use sets to enforce
        # the remaining rules
        nargs = len(kwarg)
        args = set(kwarg.keys())
        # A set of the legal arguments
        legal_args = set(['T','p','d','v','e','h','s'])
        
        # 2.1: There may only be 2 arguments
        if nargs>2:
            raise pm.utility.PMParamError(
                    'Specifying more than two simultaneous parameters is illegal.')
        
        # 2.2: All arguments must be "legal" recognized arguments
        these_args = args - legal_args
        if these_args:
            message = 'Unrecognized propert(y/ies):'
            prefix = '  '
            for name in these_args:
                message += prefix + name
                prefix = ', '
            raise pm.utility.PMParamError(message)
        
        # 2.3: T, e, and h may not be specified together
        these_args = set(['T', 'e', 'h']).intersection(args)
        if len(these_args) > 1:
            message = 'Properties may not be specified together:'
            prefix = ' '
            for name in these_args:
                message += prefix + name
                prefix = ', '
            raise pm.utility.PMParamError(message)
        
        # 2.4: Density and specific volume cannot be specified together
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
        # 5) Replace v with d and e/h with T
        if 'T' in kwarg:
            kwarg['T'] = pm.units.temperature_scale(kwarg['T'], to_units='K')
        if 'p' in kwarg:
            kwarg['p'] = pm.units.pressure(kwarg['p'], to_units='Pa')
        if 'd' in kwarg:
            value = pm.units.volume(kwarg['d'], to_units='m3', exponent=-1)
            kwarg['d'] = pm.units.matter(value, self.data['mw'], to_units='kmol')
        if 'v' in kwarg:
            # Convert and replace with d at the same time
            value = pm.units.volume(kwarg['v'], to_units='m3')
            kwarg['d'] = 1./pm.units.matter(value, self.data['mw'], to_units='kmol', exponent=-1)
        if 'h' in kwarg:
            value = kwarg['h']
            value = pm.units.energy(value, to_units='kJ')
            value = pm.units.matter(value, self.data['mw'], to_units='kmol', exponent=-1)
            Tlow = self.data['Tlim'][0]
            Thigh = self.data['Tlim'][-1]
            T = np.full(value.shape, 0.5*(Tlow+Thigh))
            I = np.ones(value.shape, dtype=bool)
            self._iter1(self._h, 'T', value, T, I, Tlow, Thigh)
            kwarg['T'] = T
        if 'e'  in kwarg:
            value = kwarg['e']
            value = pm.units.energy(value, to_units='kJ')
            value = pm.units.matter(value, self.data['mw'], to_units='kmol', exponent=-1)
            Tlow = self.data['Tlim'][0]
            Thigh = self.data['Tlim'][-1]
            T = np.full(value.shape, 0.5*(Tlow+Thigh))
            I = np.ones(value.shape, dtype=bool)
            self._iter1(self._e, 'T', value, T, I, Tlow, Thigh)
            kwarg['T'] = T
        if 's' in kwarg:
            value = kwarg['s']
            value = pm.units.energy(value, to_units='kJ')
            value = pm.units.matter(value, self.data['mw'], to_units='kmol', exponent=-1)
            value = pm.units.temperature(value, to_units='K', exponent=-1)
            kwarg['s'] = value

        # Convert R into J/kmol/K - use this for p = dRT
        # Do NOT use this for s, e, and h relationships
        R = 1000 * pm.units.const_Ru

        T = p = d = None
        # Entropy requires special iteration
        if 's' in kwarg:
            # If density is specified
            if 'd' in kwarg:
                s,d = np.broadcast_arrays(kwarg['s'], kwarg['d'])
                Tlow = self.data['Tlim'][0]
                Thigh = self.data['Tlim'][-1]
                T = np.full_like(s, 0.5*(Tlow+Thigh))
                I = np.ones_like(s,dtype=bool)
                self._iter1(self._sditer, 'T', s, T, I, Tlow, Thigh, param={'d':d})

            # If pressure is specified
            elif 'p' in kwarg:
                s,p = np.broadcast_arrays(kwarg['s'], kwarg['p'])
                # adjust entropy to the reference pressure
                s += pm.units.const_Ru * np.log(p / self._pref_pa)
                Tlow = self.data['Tlim'][0]
                Thigh = self.data['Tlim'][-1]
                T = np.full_like(s, 0.5*(Tlow+Thigh))
                I = np.ones_like(s,dtype=bool)
                self._iter1(self._s, 'T', s, T, I, Tlow, Thigh)
            # If temperature is specified
            elif 'T' in kwarg:
                s,T = np.broadcast_arrays(kwarg['s'], kwarg['T'])
                # Calculate the reference entropy at the specified temperature
                s0 = self._s(T)[0]
                p = self._pref_pa * np.exp((s0-s)/pm.units.const_Ru)
                # Otherwise, this is an illegal combination!
            else:
                raise pm.utility.PMParamError(f'Cannot simultaneously specify parameters: {args}')
        # If temperature is specified
        elif 'T' in kwarg:
            # There isn't much work to do
            if 'p' in kwarg:
                T,p = np.broadcast_arrays(kwarg['T'],kwarg['p'])
            elif 'd' in kwarg:
                T,d = np.broadcast_arrays(kwarg['T'],kwarg['d'])
            else:
                message = 'Please report a bug: Unhandled event [T] in ig2._argparse with args:'
                prefix = ' '
                for name in kwarg:
                    message += prefix + name
                    prefix = ', '
                raise pm.utility.PMParamError(message)
        # If pressure is specified
        elif 'p' in kwarg:
            if 'd' in kwarg:
                p,d = np.broadcast_arrays(kwarg['p'], kwarg['d'])
                T = p / (R * d)
            else:
                message = 'Please report a bug: Unhandled event [p] in ig2._argparse with args:'
                prefix = ' '
                for name in kwarg:
                    message += prefix + name
                    prefix = ', '
                raise pm.utility.PMParamError(message)
        else:
            message = 'Please report a bug: Unhandled event [MASTER] in ig2._argparse with args:'
            prefix = ' '
            for name in kwarg:
                message += prefix + name
                prefix = ', '
            raise pm.utility.PMParamError(message)
            
        # Test the temperatures for out-of-bounds
        I = np.logical_or(T < self.data['Tlim'][0], T > self.data['Tlim'][-1])
        if I.all():
            raise pm.utility.PMParamError('All of the specified states were out-of-bounds.  '  
                    'Legal temperatures for {} are between {} and {} Kelvin.'.format(self.data['id'], self.data['Tlim'][0], self.data['Tlim'][-1]))
        elif I.any():
            T[I] = pm.config['def_oob']
            pm.utility.print_warning('Some of the states were out of bounds - setting to config[\'def_oob\'].  '
                    'Legal temperatures for {} are between {} and {} Kelvin.'.format(self.data['id'], self.data['Tlim'][0], self.data['Tlim'][-1]))
        return T,p,d

    def _crange(self, T, index):
        """Return the down-select boolean array for temperature values in range "index".
    I = _crange(T,index)

I is a boolean array used for boolean indexing of the state properties.    
I will be returned such that Tlim[index] <= T[I] < Tlim[index]
"""
        # If the index corresponds to the last range segment, then make
        # the upper boundary inclusive.
        if index == len(self.data['Tlim'])-2:
            return np.logical_and(\
                    T >= self.data['Tlim'][index],
                    T <= self.data['Tlim'][index+1])
                    
        return np.logical_and(\
                    T >= self.data['Tlim'][index],
                    T < self.data['Tlim'][index+1])


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

        arg = param.copy()
        count = 0

        if verbose:
            print('x, yy, yyx, dx, Ids')
        while Ids.any():
            # Build the new argument list
            for k,v in param.items():
                # For any array arguments, shrink them along with Ids
                if isinstance(v,np.ndarray):
                    arg[k] = v[Ids]
            # Shrink the primary property array
            arg[prop] = x[Ids]
            # Evaluate the funciton and isolate its derivative
            FF = fn( diff=True, **arg)
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
            IooB = np.logical_or( x < xmin, x > xmax)
            count_oob = 0
            while IooB.any():
                dx[IooB] /= 2.
                x[IooB] -= dx[IooB]
                IooB = np.logical_or( x < xmin, x > xmax)
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
                
    def _sditer(self, T, d, diff=1):
        s,sT = self._s(T, diff)
        R = 1000 * pm.units.const_Ru
        s -= pm.units.const_Ru * np.log(d * R * T / self._pref_pa)
        if diff:
            sT -= pm.units.const_Ru/T
        else:
            sT = None
        return s,sT
        
        

    def _cp(self,T):
        """Constant pressure specific heat (inner routine)
    _cp(T)

Expects temperature in Kelvin and returns cp in kJ/kmol/K"""

        out = np.full_like(T,pm.config['def_oob'],dtype=float)
        # Loop through the available piece-wise temperature ranges
        # Use the _crange() method to identify the elements appropriate
        # for each.
        for ii in range(len(self.data['Tlim'])-1):
            C = self.data['C'][ii]
            I = self._crange(T, ii)
            t = T[I] / 1000.
            out[I] = C[0] + t*(C[1] + t*(C[2] + t*C[3])) + C[4] / (t*t)
        return out
        
    def _h(self, T, diff=False):
        """Enthalpy (inner routine)
    h = _h(T)
        OR
    h,hT = _h(T, diff=True)
    
Expects T in Kelvin and returns enthalpy in kJ/kmol
When diff=True, the partial derivative of enthalpy with respect to 
temperature with constant pressure is also returned.
"""
        out = np.full_like(T,pm.config['def_oob'],dtype=float)
        hT = None
        if diff:
            hT = np.full_like(T,pm.config['def_oob'],dtype=float)
        # Loop through the available piece-wise temperature ranges
        # Use the _crange() method to identify the elements appropriate
        # for each.
        for ii in range(len(self.data['Tlim'])-1):
            C = self.data['C'][ii]
            I = self._crange(T, ii)
            t = T[I] / 1000.
            out[I] = C[5] + t*(C[0] + t*(C[1]/2. + t*(C[2]/3. + t*C[3]/4.))) - C[4]/t
            if diff:
                hT[I] = C[0] + t*(C[1] + t*(C[2] + t*C[3])) + C[4]/(t*t)
        # Rescale for temperature integral
        out *= 1000.
        return out, hT
        
    def _e(self, T, diff=False):
        """Internal Energy (inner routine)
    e = _e(T)
        OR
    e,eT = _e(T, diff=True)
    
Expects T in Kelvin and returns energy in kJ/kmol
When diff=True, the partial derivative of energy with respect to 
temperature with constant pressure is also returned.
"""
        h,hT = self._h(T,diff)
        if diff:
            return h-pm.units.const_Ru*T, hT-pm.units.const_Ru
        return h-pm.units.const_Ru*T, None
        
    def _s(self, T, diff=False):
        """Entropy at reference pressure (inner routine)
    
    s, sT = _s(T,diff=True)
    
Expects T in Kelvin, and returns entropy in kJ/kmol.
When diff=True, the partial derivatives of entropy with respect to 
temperature is also returned.  Otherwise it is returned as None.
"""
        out = np.full_like(T,pm.config['def_oob'],dtype=float)
        sT = None
        if diff:
            sT = np.full_like(T,pm.config['def_oob'],dtype=float)
        # Loop through the available piece-wise temperature ranges
        # Use the _crange() method to identify the elements appropriate
        # for each.
        for ii in range(len(self.data['Tlim'])-1):
            C = self.data['C'][ii]
            I = self._crange(T, ii)
            t = T[I] / 1000.
            out[I] = C[6] + C[0]*np.log(t) + t*(C[1] + t*(C[2]/2. + t*C[3]/3.)) - C[4] / (2*t*t)
            if diff:
                sT[I] = C[0]/t + C[1] + t*(C[2] + t*C[3]) + C[4]/(t*t*t)
                # Rescale for temperature derivative
                sT[I] /= 1000
        return out, sT


    def _g(self, T, diff=False):
        """Gibbs energy at reference pressure (inner routine)
    
    g, gT = _g(T,diff=True)
    
Expects T in Kelvin, and returns Gibbs energy in kJ/kmol.
When diff=True, the partial derivatives of Gibbs energy with respect to 
temperature is also returned.  Otherwise it is returned as None.
"""
        out = np.full_like(T,pm.config['def_oob'],dtype=float)
        gT = None
        if diff:
            gT = np.full_like(T,pm.config['def_oob'],dtype=float)
        # Loop through the available piece-wise temperature ranges
        # Use the _crange() method to identify the elements appropriate
        # for each.
        for ii in range(len(self.data['Tlim'])-1):
            C = self.data['C'][ii]
            I = self._crange(T, ii)
            t = T[I] / 1000.
            logt = np.log(t)
            out[I] = C[5] - C[0]*t*logt + t*(C[0]-C[6] + t*(-0.5*C[1] + t*(-1./6.*C[2] + t*(-1./12.*C[3])))) - 0.5*C[4]/t
            if diff:
                gT[I] = -C[0]*(logt+1) + C[0]-C[6] + t*(-C[1] + t*(-0.5*C[2] + t*(-1./3.*C[3]))) + 0.5*C[4]/t/t
        out *= 1000.
        return out, gT

    def _f(self, T, diff=False):
        """Free (Helmholtz) energy at reference pressure (inner routine)
    
    f, fT = _f(T,diff=True)
    
Expects T in Kelvin, and returns free energy in kJ/kmol.
When diff=True, the partial derivatives of free energy with respect to 
temperature is also returned.  Otherwise it is returned as None.
"""
        g,gT = self._g(T, diff=diff)
        
        if diff:
            return g-pm.units.const_Ru*T, gT-pm.units.const_Ru
        return g-pm.units.const_Ru*T, None

    def _test(self, report=None):
        """Test the ig data model against a series of criteria
        
    _test()     # Prints results to stdout
        OR
    _test(report='/path/to/file')   # Writes a report file
        OR
    _test(report=open_file_descriptor)  # Appends to an open report file
    
Returns True when all criteria are satisfied and False otherwise.  Unlike ig2,
the ig class uses embedded tabulated data for numerical tests.  No external 
validation table is needed.

The criteria are:
1. Data integrity
    1.1 Tlim must be monotonic
    1.2 Tlim dimensions must match coefficient dimensions
            len(Tlim) == len(C)+1
    1.3 There must be 8 coefficients
            len(C[ii]) == 8 for all ii in range(len(C))
            
DATA INTEGRITY FAILURES ARE FATAL
FAILURES HERE WILL HALT THE TEST

2. Model continuity
    2.1 cp should not have discontinuities at the Tlim boundaries
    2.2 h should not have discontinuities
    2.3 s should not have discontinuities
    
3. Tabulated reference data should agree with the model in the currently
        configured unit system.
    3.1 cp should agree with tabulated reference data to within 0.01%
    3.2 s should agree with reference data to within 0.1 J/mol/K
    3.3 h should agree with reference data to wihtin 0.01 kJ/mol/K
    3.4 density at 1000K and 10bar should agree to within .01%.
    3.5 R should agree to within .01%.
    3.6 cv should match cp-R at all tabulated conditions to within .01%.
    3.7 e should match h-RT at all tabulated conditions to within .01%.
    3.8 gam should match cp/cv at all tabulated conditions to .0001
    3.9 T_h() should match T to within 0.01%
    3.10 T_s() should match T to within 0.01%
    
Optional keywords that configure the test:
keyword (default)   
Description

"""
        # Recurse with a fresh file descriptor if the file is a string
        if isinstance(report, str):
            with open(report, 'w') as ff:
                return self._test(report=ff)
        elif report is None:
            report = sys.stdout
        
        # Perform the checks
        result = True
                
        Tlim = self.data['Tlim']
        C = self.data['C']
        
        report.write(repr(self) + '\n1. Data integrity\n')
        
        # 1.1: Tlim must be monotonic
        test = True
        for t0,t1 in zip(Tlim[:-1], Tlim[1:]):
            test = test and (t0<t1)
        if test:
            report.write('[passed]')
        else:
            report.write('[FAILED]')
        report.write('    1.1: Tlim must increase monotonically\n')
        result = result and test
        
        # 1.2: Tlim and C dims must match
        test = (len(Tlim) == len(C)+1)
        result = result and test
        if test:
            report.write('[passed]')
        else:
            report.write('[FAILED]')
        report.write('    1.2: Tlim and C dimensions must be compatible\n')
        
        # 1.3: C must have 8 coefficients
        test = True
        for cc in C:
            test = test and (len(cc)==8)
        result = result and test
        if test:
            report.write('[passed]')
        else:
            report.write('[FAILED]')
        report.write('    1.3: All temperature regions must have 8 coefficients\n')
        
        if not result:
            report.write('[FATAL] Data integrity test failed. Aborting further checks.\n')
            return False
        
        # 2. tests for discontinuities
        report.write('2. Model continuity\n')
        cptest = []
        htest = []
        stest = []
        for T in Tlim[1:-1]:
            T = pm.units.temperature_scale(T, from_units='K')
            # Perturb the temperature by +/- 0.01%
            TT = T * np.array([0.9999, 1.0001])
            # Check specific heat continuity to within 0.01%
            cp = self.cp(TT)
            if np.abs(cp[0] - cp[1]) / cp[0] > .0001:
                cptest.append(T)
            # Check enthalpy continuity to within 0.01%
            h = self.h(TT)
            # Adjust for the known change in temperature
            h[1] -= cp[1] * T * .0001
            h[0] += cp[0] * T * .0001
            if np.abs(h[0] - h[1]) / h[0] > .0001:
                htest.append(T)
            s = self.s(TT)
            # Adjust for the known change in temperature
            s[1] -= cp[1] * .0001
            s[0] += cp[0] * .0001
            if np.abs(s[0] - s[1]) / s[0] > .0001:
                stest.append(T)
        
        if cptest:
            report.write('[FAILED]')
            result = False
        else:
            report.write('[passed]')
        report.write('    2.1 cp() must be continuous at piecewise boundaries.\n')
        if cptest:
            report.write('                Failure at T=')
            for T in cptest:
                report.write('%.2f,'%T)
            report.write('\n')
        
        if htest:
            report.write('[FAILED]')
            result = False
        else:
            report.write('[passed]')
        report.write('    2.2 h() must be continuous at piecewise boundaries.\n')
        if htest:
            report.write('                Failure at T=')
            for T in htest:
                report.write('%.2f,'%T)
            report.write('\n')

        if stest:
            report.write('[FAILED]')
            result = False
        else:
            report.write('[passed]')
        report.write('    2.3 s() must be continuous at piecewise boundaries.\n')
        if stest:
            report.write('                Failure at T=')
            for T in stest:
                report.write('%.2f,'%T)
            report.write('\n')
            
        # 3. Test for model consistency with tabulated values
        # Start by converting the table into the currently configured units
        report.write('3. Numerical consistency checks\n')
        TAB = np.array(self.data['TAB'])
        T = pm.units.temperature_scale(TAB[:,0], from_units='K')
        pref = pm.units.pressure(self._pref_pa, from_units='pa')
        cp_test = self.cp(T)
        cv_test = self.cv(T)
        h_test = self.h(T)
        e_test = self.e(T)
        gam_test = self.gam(T)
        s_test = self.s(T=T, p=pref)
        
        
        cp = pm.units.energy(TAB[:,1], from_units='J')
        pm.units.matter(cp, self.data['mw'], from_units='mol', inplace=True, exponent=-1)
        pm.units.temperature(cp, from_units='K', inplace=True, exponent=-1)
        
        s = pm.units.energy(TAB[:,2], from_units='J')
        pm.units.matter(s, self.data['mw'], from_units='mol', inplace=True, exponent=-1)
        pm.units.temperature(s, from_units='K', inplace=True, exponent=-1)
        # Convert the entropy error threshold too
        serr = pm.units.energy(0.1, from_units='J')
        serr = pm.units.matter(serr, self.data['mw'], from_units='mol', exponent=-1)
        serr = pm.units.temperature(serr, from_units='K', exponent=-1)
        
        h0 = self.h(pm.units.temperature_scale(298.15, from_units='K'))
        h = pm.units.energy(TAB[:,4], from_units='kJ')
        pm.units.matter(h, self.data['mw'], from_units='mol', inplace=True, exponent=-1)
        
        # Convert the enthalpy error threshold too
        # herr is also used below on energy checks
        herr = pm.units.energy(0.01, from_units='kJ')
        herr = pm.units.matter(herr, self.data['mw'], from_units='mol', exponent=-1)

        I = (np.abs(cp_test - cp)/cp > .001)
        if np.any(I):
            result = False
            report.write('[FAILED]')
        else:
            report.write('[passed]')
        report.write('    3.1 cp must agree with tabulated data to within 0.1%\n')
        if np.all(I):
            report.write('            Failed at all temperatures.\n')
        elif np.any(I):
            report.write('            Failed at T=')
            for tt in T[I]:
                report.write('%.2f,'%tt)
            report.write('\n')
            
        I = np.abs(s_test - s) > serr
        if np.any(I):
            report.write('[FAILED]')
            result = False
        else:
            report.write('[passed]')
        report.write('    3.2 s must agree with tabulated data to within 0.1 J/mol/K\n')
        if np.all(I):
            report.write('            Failed at all temperatures.\n')
        elif np.any(I):
            report.write('            Failed at T=')
            for tt in T[I]:
                report.write('%.2f,'%tt)
            report.write('\n')
        
        # 3.3 h must agree with tabulated data
        I = np.abs(h_test - h0 - h) > np.maximum(herr, np.abs(h) * .001)
        if np.any(I):
            report.write('[FAILED]')
            result = False
        else:
            report.write('[passed]')
        report.write('    3.3 h must agree with tabulated data to within 0.01 kJ/mol or 0.1%\n')
        if np.all(I):
            report.write('            Failed at all temperatures.\n')
        elif np.any(I):
            report.write('            Failed at T=')
            for tt in T[I]:
                report.write('%.2f,'%tt)
            report.write('\n')
        
        # 3.4 density at 1000K and 10bar should agree to within .01%.
        TT = 1000.
        pp = 10.
        dd = pp*1e5 / (pm.units.const_Ru * 1000.)
        dd = pm.units.matter(dd, self.data['mw'], from_units='mol')
        dd = pm.units.volume(dd, from_units='m3', exponent=-1)
        TT = pm.units.temperature_scale(TT, from_units = 'K')
        pp = pm.units.pressure(pp, from_units='bar')
        
        if np.abs(self.d(T=TT, p=pp) - dd)/dd > .0001:
            report.write('[FAILED]')
            result = False
        else:
            report.write('[passed]')
        report.write('    3.4 Density at 1000K and 10bar must match IG law to within .01%\n')
        
        # 3.5 R should agree to within .01%.
        R = pm.units.energy(pm.units.const_Ru, from_units='J')
        R = pm.units.matter(R, self.data['mw'], from_units='mol', exponent=-1)
        R = pm.units.temperature(R, from_units='K', exponent=-1)
        if np.abs(self.R() - R)/R > .0001:
            report.write('[FAILED]')
            result = False
        else:
            report.write('[passed]')
        report.write('    3.5 Ideal gas constant must match to within 0.01%\n')
        
        # 3.6 cv should match cp-R at all tabulated conditions to within .01%.
        I = np.abs(cv_test + self.R() - cp_test)/cv_test > .0001
        if np.any(I):
            report.write('[FAILED]')
            result = False
        else:
            report.write('[passed]')
        report.write('    3.6: cv == cp - R to within .01% at all tabulated values\n')
        if np.all(I):
            report.write('            Failed at all temperatures.\n')
        elif np.any(I):
            report.write('            Failed at T=')
            for tt in T[I]:
                report.write('%.2f,'%tt)
            report.write('\n')
        
        # 3.7 e should match h-RT at all tabulated conditions to within .01 kJ/mol.
        # Borrow the same herr as from the enthalpy checks
        I = np.abs(e_test + self.R()*T - h_test)/e_test > herr
        if np.any(I):
            report.write('[FAILED]')
            result = False
        else:
            report.write('[passed]')
        report.write('    3.7: e == h - R*T to within .01 kJ/mol at all tabulated values\n')
        if np.all(I):
            report.write('            Failed at all temperatures.\n')
        elif np.any(I):
            report.write('            Failed at T=')
            for tt in T[I]:
                report.write('%.2f,'%tt)
            report.write('\n')
    
        # 3.8 gam should match cp/cv at all tabulated conditions to .0001
        I = np.abs(gam_test - cp_test/cv_test) > .0001
        if np.any(I):
            report.write('[FAILED]')
            result = False
        else:
            report.write('[passed]')
        report.write('    3.8: gam = cp/cv to within .001 at all tabulated values\n')
        if np.all(I):
            report.write('            Failed at all temperatures.\n')
        elif np.any(I):
            report.write('            Failed at T=')
            for tt in T[I]:
                report.write('%.2f,'%tt)
            report.write('\n')
            
        # Trim the temperature, entropy and enthalpy arrays of the limiting boundaries
        T = T[1:-1]
        h = h_test[1:-1]
        s = s_test[1:-1]
        
        # 3.9 T_h should match T at all tabulated conditions to within .01%
        I = np.abs(T - self.T_h(h))/T > .0001
        if np.any(I):
            report.write('[FAILED]')
            result = False
        else:
            report.write('[passed]')
        report.write('    3.9: T == T_h(h) to within .01% at all tabulated values.\n')
        if np.all(I):
            report.write('            Failed at all temperatures.\n')
        elif np.any(I):
            report.write('            Failed at T=')
            for tt in T[I]:
                report.write('%.2f,'%tt)
            report.write('\n')
            
        # 3.10 T_s should match T at all tabulated conditions to within .01%
        I = np.abs(T - self.T_s(s,p=pref))/T > .0001
        if np.any(I):
            report.write('[FAILED]')
            result = False
        else:
            report.write('[passed]')
        report.write('    3.10: T == T_s(s) to within .01% at all tabulated values.\n')
        if np.all(I):
            report.write('            Failed at all temperatures.\n')
        elif np.any(I):
            report.write('            Failed at T=')
            for tt in T[I]:
                report.write('%.2f,'%tt)
            report.write('\n')
            
            
        return result
        

    def __test(self, report_file=None, report_level=2, basic=False):
        """Test the data and algorithm against tabulated data
    _test(report_file=None, report_level=2, basic=False)

returns True on success and False on failure.

report_file     The file name or file stream to which a test report 
                should be streamed.  (string or an open file object)
                When report_file is unspecified, it defaults to stdout.
report_level    Accepts an integer indicating the level of detail to
                print to the report file.  
                0 - Nothing, just return the boolean
                1 - Minimal, summarize the failures
                2 - Vital, detail the failures
                3 - All; report everything

basic           When True, only criterion 0 is assessed.

The following test criteria are used:
0)  The data should be correctly organized
    0a) Tlim must be monotonic
    0b) Tlim should have one more element than C has rows.
    0c) C should have 8 columns
    0d) Specific heat should have no discontinuities

1)  Tabulated specific heats should agree with cp()
    1a) Using molar units
    1b) Using mass units
    1c) Using a non-absolute temperature scale
    1d) Using an energy unit other than J

2)  Tabulated values should agree with h()
    2a) Using molar units
    2b) Using mass units
    2c) Using a non-absolute temperature scale
    2d) Using an energy unit other than J

3)  Tabulated values should agree with s()
    3a) Using molar units
    3b) Using mass units
    3c) Using a non-absolute temperature scale
    3d) Using an energy unit other than J

4)  R() should agree with explicit calculations for R
    4a) Using molar units
    4b) Using mass units

5)  The density should agree with explicit ideal gas density
    5a) Using mass units
    5b) Using a non-absolute temperature scale

6)  Constant-volume specific heat should equal constant pressure minus R

7)  Internal energy should equal enthalpy minus RT

8)  Specific heat ratio should equal an explicit evaluation of cp and cv

"""
        # Initialize a results dicitonary
        result = True
        # Assign some constants
        REP_ALL = 3
        REP_VITAL = 2
        REP_MINIMAL = 1
        REP_NONE = 0
        # Deal with the file 
        if report_file is None:
            ff = os.sys.stdout
        elif isinstance(report_file,str):
            ff = open(report_file)
        elif hasattr(report_file,'write'):
            ff = report_file
        else:
            raise pm.utility.PMParamError('Unrecognized file type')

        # Start with some basic stuff
        

        # Set up a numpy array with the tabulated values
        # Remove the first and last rows to prevent floating point
        # precision violations at the temperature limits
        TAB = np.array(self.data['TAB'])[1:-1,:]
        Tref = TAB[:,0]         # K
        cpref = TAB[:,1]        # J/mol/K
        sref = TAB[:,2]         # J/mol/K
        TAB[:,4] *= 1000.       # Convert h to J/mol in place
        href = TAB[:,4]         # J/mol

        ##
        pref = self._pref_pa / 1e5

        # Create a generic property test function
        # crit  criterion string for the report (1a, 2b, etc...)
        # prop  property method to call (cp, h, s)
        # args  the argv tuple to pass to the property function
        # ref   reference array
        # units = [e, m, p, t] units for energy, matter, pressure, and temperature
        # epsilon  The acceptable fractional error 
        # small    An acceptable total error (for values that can be close to 0)
        def _prop_test(crit, prop, args, ref, units, epsilon=.001, small=0.):
            pm.config['unit_energy'] = units[0]
            pm.config['unit_matter'] = units[1]
            pm.config['unit_pressure'] = units[2]
            pm.config['unit_temperature'] = units[3]
            
            if report_level >= REP_VITAL:
                ff.write(crit + '...')

            try:
                error = abs(prop(*args) - ref)
            except:
                if report_level>=REP_VITAL:
                    ff.write('[ERROR]\n')
                    if report_level>=REP_ALL:
                      ff.write('    ' + repr(os.sys.exc_info()[1].args))
                return False

            fail = (error > abs(epsilon*ref)) * (error > small)
            # Write to the report
            if fail.any():
                if report_level>=REP_VITAL:
                    ff.write('%10f...[FAIL]\n'%error.max())
                if report_level>=REP_ALL:
                    if fail.all():
                        ff.write('    Failed at all temperatures.\n')
                    else:
                        ff.write('    Failed at T(K)=')
                        for tt in Tref[fail]:
                            ff.write('%d,'%(int(tt)))
                        ff.write('\n')
            elif report_level>=REP_VITAL:
                ff.write('%10f...[pass]\n'%error.max())

            if fail.any():
                return False
            else:
                return True

        if report_level >= REP_MINIMAL:
            ff.write('Testing IG species %s\n'%self.data['id'])

        #==========================================================#
        # Criterion 0 - Data integrity checks                      #
        #==========================================================#
        subresult = True
        ##
        ## Criterion 0a - Tlim is monotonic
        ##
        if report_level >= REP_VITAL:
            ff.write('  0a...')
        if (np.diff(self.data['Tlim']) < 0.).any():
            subresult = False
            if report_level >= REP_VITAL:
                ff.write('[FAIL]\n')
        elif report_level >= REP_VITAL:
            ff.write('[pass]\n')
            
        ##
        ## Criterion 0b - Tlim and C are compatible lengths
        ##
        NLim = len(self.data['Tlim'])
        if report_level >= REP_VITAL:
            ff.write('  0b...')

        fail = len(self.data['C']) != NLim - 1
        if fail:
            subresult = False
            if report_level >= REP_VITAL:
                ff.write('[FAIL]\n')
        elif report_level >= REP_VITAL:
            ff.write('[pass]\n')

        ##
        ## Criterion 0c - 8-element coefficients
        ##
        if report_level >= REP_VITAL:
            ff.write('  0c...')

        local = True
        for C in self.data['C']:
            if len(C) != 8:
                local = False

        if not local:
            if report_level >= REP_VITAL:
                ff.write('[FAIL]\n')
        elif report_level >= REP_VITAL:
            ff.write('[pass]\n')

        subresult &= local


        ##
        ## Criterion 0d - No discontinuities in cp()
        ##
        if len(self.data['Tlim']) > 2:
            if report_level >= REP_VITAL:
                ff.write('  0d...')

            pm.config['unit_temperature']='K'
            pm.config['unit_energy']='J'
            pm.config['unit_matter']='mol'
            T = np.array(self.data['Tlim'][1:-1])
            clow = self.cp(T-.01)
            chigh = self.cp(T+.01)
            fail = (abs(chigh-clow) > .001*chigh)
            if fail.any():
                subresult = False
                if report_level >= REP_VITAL:
                    ff.write('[FAIL]\n')
                    if report_level >= REP_ALL:
                        ff.write('    Discontinuity at T(K)=')
                        for tt in T[fail]:
                            ff.write('%d,'%tt)
                        ff.write('\n')
            elif report_level >= REP_VITAL:
                ff.write('[pass]\n')


        if not subresult:
            if report_level >= REP_MINIMAL:
                ff.write('Criterion 0: Fundamental format [FAIL]\n')                
        elif report_level >= REP_MINIMAL:
            ff.write('Criterion 0: Fundamental format [pass]\n')

        if basic:
            return subresult

        result &= subresult
        #==========================================================#
        # Criterion 1 - cp should agree with tabulated data        #
        #==========================================================#
        subresult = True
        ##
        ## Criterion 1a - molar units
        ##
        units = ['J','mol','bar','K']
        tag = '  1a %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag,self.cp,(Tref,),cpref,units)
        ##
        ## Criterion 1b - mass units
        ##
        units = ['J','kg','bar','K']
        tag = '  1b %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag, self.cp, (Tref,),
            cpref*1000./self.data['mw'],units)
        ##
        ## Criterion 1c - A non-absolute temperature scale
        ##
        units = ['J','mol','bar','F']
        tag = '  1c %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag,self.cp,(Tref*1.8-459.67,),
            cpref/1.8,units)
        ##
        ## Criterion 1d - Energy units other than J
        ##
        units = ['BTU','mol','bar','K']
        tag = '  1d %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag,self.cp,(Tref,),
            cpref/1054.35,units)

        result &= subresult

        if report_level >= REP_MINIMAL:
            ff.write('Criterion 1: specific heat ')
            if subresult:
                ff.write('[pass]\n')
            else:
                ff.write('[FAIL]\n')
        #==========================================================#
        # Criterion 2 - h should agree with tabulated data         #
        #==========================================================#
        subresult = True
        ##
        ## Criterion 2a - molar units
        ##
        small = 5.
        units = ['J','mol','bar','K']
        tag = '  2a %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag,self.h,
            (Tref,None,False), href, units, small=small)
        ##
        ## Criterion 2b - mass units
        ##
        scale = 1000. / self.data['mw']
        small = 5. * scale
        units = ['J','kg','bar','K']
        tag = '  2b %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag,self.h, (Tref,None,False), 
            href*scale, units, small=small)
        ##
        ## Criterion 2c - A non-absolute temperature scale
        ##
        small = 5.
        units = ['J','mol','bar','F']
        tag = '  2c %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag,self.h,(Tref*1.8-459.67,None,False),
            href , units, small=small)
        ##
        ## Criterion 2d - Energy units other than J
        ##
        scale = 1./1054.35026
        small = 5. * scale
        units = ['BTU','mol','bar','K']
        tag = '  2d %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag,self.h, (Tref,None,False), 
            href*scale, units, small=small)

        result &= subresult

        if report_level >= REP_MINIMAL:
            ff.write('Criterion 2: ethalpy ')
            if subresult:
                ff.write('[pass]\n')
            else:
                ff.write('[FAIL]\n')

        #==========================================================#
        # Criterion 3 - s should agree with tabulated data         #
        #==========================================================#
        subresult = True
        p0 = self._pref_pa / 1e5
        ##
        ## Criterion 3a - molar units
        ##
        units = ['J','mol','bar','K']
        tag = '  3a %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag,self.s,
            (Tref,p0), sref, units)
        ##
        ## Criterion 3b - mass units
        ##
        scale = 1000. / self.data['mw']
        units = ['J','kg','bar','K']
        tag = '  3b %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag,self.s, (Tref,p0), 
            sref*scale, units)
        ##
        ## Criterion 3c - A non-absolute temperature scale
        ##
        units = ['J','mol','bar','F']
        tag = '  2c %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag,self.s,(Tref*1.8-459.67,p0),
            sref/1.8 , units)
        ##
        ## Criterion 3d - Energy units other than J
        ##
        scale = 1./1054.35026
        units = ['BTU','mol','bar','K']
        tag = '  3d %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag,self.s, (Tref,p0), 
            sref*scale, units, small=small)

        result &= subresult

        if report_level >= REP_MINIMAL:
            ff.write('Criterion 3: entropy ')
            if subresult:
                ff.write('[pass]\n')
            else:
                ff.write('[FAIL]\n')

        #==========================================================#
        # Criterion 4 - R should agree with explicit calculation   #
        #==========================================================#
        subresult = True
        ## 
        ## Criterion 4a - Using molar units
        ##
        pm.config['unit_energy'] = 'J'
        pm.config['unit_matter'] = 'mol'
        pm.config['unit_temperature'] = 'K'
        error = abs(self.R() - pm.units.const_Ru)
        fail = (error > .0001 * pm.units.const_Ru)
        subresult &= not fail
        if report_level >= REP_VITAL:
            ff.write('  4a   J mol bar   K...%10f...'%error)
            if fail:
                ff.write('[FAIL]\n')
            else:
                ff.write('[pass]\n')
        ## 
        ## Criterion 4b - Using mass units
        ##
        pm.config['unit_energy'] = 'J'
        pm.config['unit_matter'] = 'kg'
        pm.config['unit_temperature'] = 'K'
        R = pm.units.const_Ru * 1000. / self.data['mw']
        error = abs(self.R() - R)
        fail = (error > .0001 * R)
        subresult &= not fail

        if report_level >= REP_VITAL:
            ff.write('  4b   J  kg bar   K...%10f...'%error)
            if fail:
                ff.write('[FAIL]\n')
            else:
                ff.write('[pass]\n')

        if report_level >= REP_MINIMAL:
            ff.write('Criterion 4: ideal gas constant ')
            if subresult:
                ff.write('[pass]\n')
            else:
                ff.write('[FAIL]\n')

        result &= subresult        

        #==========================================================#
        # Criterion 5 - density should match IG density            #
        #==========================================================#
        subresult = True
        pm.config['unit_volume'] = 'm3'
        units = ['J','kg','bar','K']

        ##
        ## 5a - mass units
        ##
        p = 5.0
        R = pm.units.const_Ru * 1000. / self.data['mw']
        d = p * 1e5 / R / Tref
        subresult &= _prop_test('  5a          kg/m3', self.d, 
            (Tref,p), d, units)

        ##
        ## 5b - non-asolute temperature scale
        ## 
        units = ['J','kg','bar','F']
        subresult &= _prop_test('  5b          kg/m3', self.d, 
            (Tref*1.8-459.67,p), d, units)

        if report_level >= REP_MINIMAL:
            ff.write('Criterion 5: density ')
            if subresult:
                ff.write('[pass]\n')
            else:
                ff.write('[FAIL]\n')

        result &= subresult

        #==========================================================#
        # Criterion 6 - constant-volume specific heat              #
        #==========================================================#
        subresult = True
        units = ['J','mol','bar','K']
        R = pm.units.const_Ru
        tag = '   6 %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag, self.cv, (Tref,), cpref-R, units)

        if report_level >= REP_MINIMAL:
            ff.write('Criterion 6: constant volume specific heat ')
            if subresult:
                ff.write('[pass]\n')
            else:
                ff.write('[FAIL]\n')

        result &= subresult
        #==========================================================#
        # Criterion 7 - internal energy                            #
        #==========================================================#
        subresult = True
        units = ['J','mol','bar','K']
        R = pm.units.const_Ru
        tag = '   7 %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag, self.e, (Tref,None,False), href-R*Tref, units, small=10.)

        if report_level >= REP_MINIMAL:
            ff.write('Criterion 7: internal energy ')
            if subresult:
                ff.write('[pass]\n')
            else:
                ff.write('[FAIL]\n')

        result &= subresult
        #==========================================================#
        # Criterion 8 - specific heat ratio                        #
        #==========================================================#
        subresult = True
        units = ['J','mol','bar','K']
        R = pm.units.const_Ru
        tag = '   8 %3s %3s %3s %3s'%tuple(units)
        subresult &= _prop_test(tag, self.gam, (Tref,), cpref/(cpref-R), units)

        if report_level >= REP_MINIMAL:
            ff.write('Criterion 8: specific heat ratio')
            if subresult:
                ff.write('[pass]\n\n')
            else:
                ff.write('[FAIL]\n\n')

        result &= subresult

        if isinstance(report_file,str):
            ff.close()
        return result



    def Tlim(self):
        """Temperature limits
    (Tmin, Tmax) = Tlim()
Returns the temperature limits on the ig data set.

Accepts None
Returns unit_temperature
"""
        Tmin = pm.units.temperature_scale(self.data['Tlim'][0], from_units='K')
        Tmax = pm.units.temperature_scale(self.data['Tlim'][-1], from_units='K')
        return (Tmin,Tmax)

    ######################################
    # Pressure, density, and temperature #
    ######################################
    
    def d(self,*varg, **kwarg):
        """Density
    d(...)

All ideal gas properties accept two other properties as flexible inputs
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns density in unit_matter / unit_volume
"""
        T,p,d = self._argparse(*varg, **kwarg)
        if d is None:
            d = p / (1000*pm.units.const_Ru * T)
        scale = pm.units.matter(1., self.data['mw'], from_units='kmol')
        scale = pm.units.volume(scale, from_units='m3', exponent=-1)
        np.multiply(d, scale, out=d)
        return d
        
    def v(self,*varg, **kwarg):
        """Specific volume
    v(...)

All ideal gas properties accept two other properties as flexible inputs
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns volume in unit_volume / unit_matter
"""
        return 1./self.d(*varg, **kwarg)
        
    def T(self, *varg, **kwarg):
        """Temperature
    T(...)

All ideal gas properties accept two other properties as flexible inputs
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns temperature in unit_temperature
"""
        # T is a special case.  If there is only one parameter given,
        # the default should not be T=def_T, so we need to override
        # the behavior of _argparse.
        if len(varg) + len(kwarg) == 1:
            kwarg['p'] = pm.config['def_p']
        
        T,p,d = self._argparse(*varg, **kwarg)
        pm.units.temperature_scale(T, from_units='K', inplace=True)
        return T
        
    def p(self, *varg, **kwarg):
        """Pressure
    p(...)

All ideal gas properties accept two other properties as flexible inputs
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns pressure in unit_pressure
"""
        T,p,d = self._argparse(*varg, **kwarg)
        if p is None:
            p = d * 1000*pm.units.const_Ru * T
        pm.units.pressure(p, from_units='Pa', inplace=True)
        return p

    def mw(self, *varg, **kwarg):
        """Molecular weight
    mw(...)

Ignores the arguments are returns molecular weight as 
unit_mass / unit_molar
"""
        mw = pm.units.mass(self.data['mw'],from_units='g')
        mw = pm.units.molar(mw,from_units='mol',exponent=-1)
        return mw

    def R(self,*varg, **kwarg):
        """Ideal gas constant
    R(...)

Ignores the arguments are returns the gas constant as
unit_energy / unit_matter / unit_temperature
"""
        R = pm.units.energy(pm.units.const_Ru, from_units='J')
        R = pm.units.temperature(R, from_units='K', exponent=-1)
        R = pm.units.matter(R, self.data['mw'], from_units='mol', exponent=-1)
        return R

    def gam(self,*varg, **kwarg):
        """Specific heat ratio
    gam(...)

All ideal gas properties accept two other properties as flexible inputs
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns ideal gas ratio, which is dimensionless.
"""
        T,p,d = self._argparse(*varg, **kwarg)
        cp = self._cp(T)
        return cp/(cp-pm.units.const_Ru)


    def cp(self, *varg, **kwarg):
        """Constant-pressure specific heat
    cp(...)

All ideal gas properties accept two other properties as flexible inputs
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns specific heat in unit_energy / unit_matter / unit_temperature
"""
        T,p,d = self._argparse(*varg, **kwarg)
        # Apply the model
        out = self._cp(T)
        # calculate a conversion factor
        scale = pm.units.energy(1, from_units='kJ')
        scale = pm.units.temperature(scale, from_units='K', exponent=-1)
        scale = pm.units.matter(scale, self.data['mw'], from_units='kmol', exponent=-1)
        # Apply the conversion factor in-place and return
        np.multiply(out, scale, out=out)
        return out
        
    def cv(self,*varg, **kwarg):
        """Constant-volume specific heat
    p(...)

All ideal gas properties accept two other properties as flexible inputs
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns specific heat in unit_energy / unit_matter / unit_temperature
"""
        T,p,d = self._argparse(*varg, **kwarg)
        # Apply the model
        out = self._cp(T) - pm.units.const_Ru
        # calculate a conversion factor
        scale = pm.units.energy(1, from_units='kJ')
        scale = pm.units.temperature(scale, from_units='K', exponent=-1)
        scale = pm.units.matter(scale, self.data['mw'], from_units='kmol', exponent=-1)
        # Apply the conversion factor in-place and return
        np.multiply(out, scale, out=out)
        return out
        
    def h(self,*varg, **kwarg):
        """Enthalpy
    h(...)

All ideal gas properties accept two other properties as flexible inputs
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns enthalpy in unit_energy / unit_matter
"""
        T,p,d = self._argparse(*varg, **kwarg)
        # Apply the model
        out = self._h(T)[0]
        # calculate a conversion factor
        scale = pm.units.energy(1, from_units='kJ')
        scale = pm.units.matter(scale, self.data['mw'], from_units='kmol', exponent=-1)
        # Apply the conversion factor in-place and return
        np.multiply(out, scale, out=out)
        return out
        
    def g(self,*varg, **kwarg):
        """Gibbs energy
    g(...)
    
    g =def= h - T*s

All ideal gas properties accept two other properties as flexible inputs
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns enthalpy in unit_energy / unit_matter
"""
        T,p,d = self._argparse(*varg, **kwarg)
        if p is None:
            p = d * (1000*pm.units.const_Ru * T)
        # Apply the model
        out = self._g(T)[0] + T*pm.units.const_Ru * np.log(p/self._pref_pa)
        # calculate a conversion factor
        scale = pm.units.energy(1, from_units='kJ')
        scale = pm.units.matter(scale, self.data['mw'], from_units='kmol', exponent=-1)
        # Apply the conversion factor in-place and return
        np.multiply(out, scale, out=out)
        return out
            
            
    def f(self,*varg, **kwarg):
        """Free (Helmholtz) energy
    f(...)
    
    f =def= e - T*s

All ideal gas properties accept two other properties as flexible inputs
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns enthalpy in unit_energy / unit_matter
"""
        T,p,d = self._argparse(*varg, **kwarg)
        if p is None:
            p = d * (1000*pm.units.const_Ru * T)
        # Apply the model
        out = self._f(T)[0] + T*pm.units.const_Ru * np.log(p/self._pref_pa)
        # calculate a conversion factor
        scale = pm.units.energy(1, from_units='kJ')
        scale = pm.units.matter(scale, self.data['mw'], from_units='kmol', exponent=-1)
        # Apply the conversion factor in-place and return
        np.multiply(out, scale, out=out)
        return out
            

    def s(self,*varg,**kwarg):
        """Entropy
    s(...)

All ideal gas properties accept two other properties as flexible inputs
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns entropy in unit_energy / unit_matter / unit_temperature
"""
        T,p,d = self._argparse(*varg, **kwarg)
        if p is None:
            p = d * (1000*pm.units.const_Ru * T)
        # Apply the model
        out = self._s(T)[0] - pm.units.const_Ru * np.log(p/self._pref_pa)
        # calculate a conversion factor
        scale = pm.units.energy(1, from_units='kJ')
        scale = pm.units.temperature(scale, from_units='K', exponent=-1)
        scale = pm.units.matter(scale, self.data['mw'], from_units='kmol', exponent=-1)
        # Apply the conversion factor in-place and return
        np.multiply(out, scale, out=out)
        return out

    def e(self,*varg, **kwarg):
        """Internal energy
    e(...)

All ideal gas properties accept two other properties as flexible inputs
Below are the recognized keywords, their meaning, and the config entries
that determine their units.
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature

If no keywords are specified, the positional arguments are interpreted
as (T,p).  To configure their defaults, use the def_T and def_p config
entries.

Returns internal energy in unit_energy / unit_matter
"""
        T,p,d = self._argparse(*varg, **kwarg)
        # Apply the model
        out = self._h(T)[0] - pm.units.const_Ru*T
        # calculate a conversion factor
        scale = pm.units.energy(1, from_units='J')
        scale = pm.units.matter(scale, self.data['mw'], from_units='mol', exponent=-1)
        # Apply the conversion factor in-place and return
        np.multiply(out, scale, out=out)
        return out

    def state(self, *varg, **kwarg):
        """Calculate all properties at a state
    sd = state(...)

The properties are returned in a dictionary with keys:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume / unit_matter
    e   internal energy     unit_energy / unit_matter
    h   enthalpy            unit_energy / unit_matter
    s   entropy             unit_energy / unit_matter / unit_temperature
    cp  const. p sp. ht.    unit_energy / unit_matter / unit_temperature
    cv  const. v sp. ht.    unit_energy / unit_matter / unit_temperature
    
Like all of the other property functions, arguments may be any two of
T, p, d, v, e, h, and s.  
"""
        Ru = pm.units.const_Ru
        T,p,d = self._argparse(*varg, **kwarg)
        # Make sure we have both pressure and density
        if d is None:
            d = p / (1000 * Ru * T)
        elif p is None:
            p = 1000 * d * Ru * T

        # Now entropy
        s = self._s(T,False)[0] - Ru * np.log(p / self._pref_pa)

        # Enthalpy and specific heat at once
        h,cp = self._h(T,True)
        cv = cp - Ru
        gam = cp/cv
        e = h - Ru*T
        
        # Finally build the output
        out = {}
        out['T'] = pm.units.temperature_scale(T, from_units='K')
        out['p'] = pm.units.pressure(p, from_units='Pa')
        scale = pm.units.matter(1., self.data['mw'], from_units='kmol')
        scale = pm.units.volume(scale, from_units='m3', exponent=-1)
        out['d'] = scale * d
        out['v'] = 1./out['d']
        scale = pm.units.energy(1., from_units='kJ')
        scale = pm.units.matter(scale, self.data['mw'], from_units='kmol', exponent=-1)
        out['h'] = h * scale
        out['e'] = e * scale
        out['gam'] = gam
        scale = pm.units.temperature(scale, from_units='K', exponent=-1)
        out['s'] = s * scale
        out['cp'] = cp * scale
        out['cv'] = cv * scale
        return out



    def T_s(self,s,p=None, d=None):
        """Temperature as a function of entropy
** Deprecated - use T() **
        
    T = T_s(s)
        or
    T = T_s(s,p)

Accepts unit_energy / unit_matter / unit_temperature
        unit_pressure
        unit_matter / unit_volume
Returns unit_temperature
"""
        if p is not None:
            return self.T(s=s, p=p)
        elif d is not None:
            return self.T(s=s, d=d)
        return self.T(s=s)


    def T_h(self,h, p=None, d=None):
        """Temperature as a function of enthalpy
** Deprecated - use T() **

    T = T_h(h)
        or
    T = T_h(h,...)

Returns the temperature as a function of enthalpy and pressure.  Ideal 
gas enthalpy is not a function of pressure, so the p term is merely a
placeholder.

Accepts unit_energy / unit_matter / unit_temperature
        unit_pressure
Returns unit_temperature
"""
        if p is not None:
            return self.T(h=h, p=p)
        elif d is not None:
            return self.T(h=h, d=d)
        return self.T(h=h)


    def p_s(self,s,T=None):
        """Pressure as a function of entropy
** Deprecated - use p() **
        
    p = ig_instance.p_s(s)
        or
    p = ig_instance.p_s(s,...)

Returns the pressure as a function of entropy and temperature.

Accepts unit_energy / unit_matter / unit_temperature
        unit_temperature
Returns unit_pressure
"""
        if T is not None:
            return self.p(s=s,T=T)
        return self.p(s=s)
