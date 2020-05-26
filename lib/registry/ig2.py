import pyromat as pm
import numpy as np
import os



class ig2(pm.reg.__basedata__):
    """Ideal gas class using the NASA polynomial equation of state.
This class exposes properties through member methods.  All property 
functions accept temperature in PYroMat's unit_temperature and 
pressure in PYroMat's unit_pressure.  The following are the member 
methods and their unit conventions:
  cp() spec. heat       (unit_energy / unit_temperature / unit_matter)
  cv() spec. heat       (unit_energy / unit_temperature / unit_matter)
  d()  density          (unit_matter / unit_volume)
  e()  internal energy  (unit_energy / unit_matter)
  h()  enthalpy         (unit_energy / unit_matter)
  gam()  spec. heat ratio (dless)
  mw() molecular weight (unit_mass / unit_molar)
  R()  gas constant     (unit_energy / unit_temperature / unit_matter)
  s()  entropy          (unit_energy / unit_temperature / unit_matter)

There are also routines to invert properties; e.g. calculating 
temperature from enthalpy or from entropy and pressure.
  T_h()  temperature from enthalpy
  T_d()  temperature from density and pressure
  T_s()  temperature from entropy and pressure
  p_s()  pressure from entropy and temperature
  p_d()  pressure from density and temperature

Some meta-data on the species can be obtained using methods
  Tlim() a two-element array with the min,max temperatures supported by
         the data set.
  contents()  returns a dictionary with a key entry for each atom in
              the chemical formula and the corresponding integer 
              quantity of each.
"""

    def _argparse(self, T=None, p=None, d=None,\
        temperature=False, pressure=False, density=False):
        """Parse the arguments supplied to an IG2 property method
    T = _argparse(*varg, **kwarg)
        OR
    T,p,d = _argparse(*varg, **kwarg, temperature=True, pressure=True, 
                        density=True)
    
_ARGPARSE automatically applies the default temperature and pressure,
def_T or def_p, from the pyromat.config system to deal with unspecified
parameters.  All inputs are re-cast as numpy arrays of at least one 
dimension and inputs are automatically converted from the configured 
user units into kJ, kmol, m^3.

The returned variables are arrays of temperature, T, pressure, p, and 
the density, d.  The optional keywords TEMPERATURE, PRESSURE, and 
DENSITY are used to indicate which state variables should be returned.
They are always returned in the order T, p, d.  
"""
        nparam = ((0 if T is None else 1) + 
                (0 if p is None else 1) + 
                (0 if d is None else 1))
                
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

        # Perform the unit conversions, and format the arrays
        if T is not None:
            T = pm.units.temperature_scale(np.asarray(T,dtype=float), to_units='K')
            if T.ndim==0:
                T = np.reshape(T,(1,))
                
        if p is not None:
            p = pm.units.pressure(np.asarray(p,dtype=float), to_units='Pa')
            if p.ndim==0:
                p = np.reshape(p,(1,))
            
        if d is not None:
            d = pm.units.matter(np.asarray(d,dtype=float), to_units='kmol')
            d = pm.units.volume(d, to_units='m3', exponent=-1)
            if d.ndim==0:
                d = np.reshape(d, (1,))
        
        # Convert the IG constant to J/kmol/K
        R = 1000 * pm.units.const_Ru
        
        # Case out the specified state variables
        # There are three possible combinations
        if T is not None:
            # T,p
            if p is not None:
                # Broadcast the arrays
                T,p = np.broadcast_arrays(T,p)
                # Do we need density?
                if density:
                    d = p / (R*T)
            # T,d
            else:
                # Broadcast the arrays
                T,d = np.broadcast_arrays(T,d)
                # Do we need pressure?
                if pressure:
                    p = d*R*T
        # p,d
        else:
            # Broadcast the arrays
            p,d = np.broadcast_arrays(p,d)
            # Do we need temperature?
            if temperature:
                T = p / (R*d)
        
        out = []
        if temperature:
            out.append(T)
        if pressure:
            out.append(p)
        if density:
            out.append(d)
            
        if len(out)>1:
            return tuple(out)
        elif len(out)==1:
            return out[0]
        return
            

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



    def _cp(self, T):
        """Constant pressure specific heat
    _cp(T)

Expects temperature in Kelvin and returns cp in kJ/kmol/K
"""
        out = np.zeros_like(T,dtype=float)
        # Loop through the piece-wise temperature ranges
        for index in range(len(self.data['Tlim'])-1):
            # Which elements are in-range?
            I = self._crange(T,index)
            for c in self.data['C'][index][4::-1]:
                out[I] = c + out[I]*T[I]
        return pm.units.const_Ru * out
            

    def _h(self, T, diff=False):
        """Enthalpy
    h,hT = _h(T)

Expects temperature in Kelvin and returns h in kJ/kmol

If the optional keyword, diff, is True, then the first derivative of 
enthalpy is also returned; otherwise it is None.  This is more 
efficient than calculating specific heat separately.
"""
        out = np.zeros_like(T,dtype=float)
        dh = None
        if diff:
            dh = np.zeros_like(T,dtype=float)
        # Loop through the piece-wise temperature ranges
        for index in range(len(self.data['Tlim'])-1):
            # Which elements are in-range?
            I = self._crange(T,index)
            term = 5.
            for c in self.data['C'][index][4::-1]:
                temp = (c/term + out[I])
                if diff:
                    dh[I] = temp + T[I]*dh[I]
                out[I] = T[I]*temp
                term -= 1.
            out[I] += self.data['C'][index][5]
        if diff:
            dh *= pm.units.const_Ru
        return pm.units.const_Ru * out, dh


    def _s(self, T, diff=False):
        """Standard entropy (at p=pref)
    s,sT = _s(T)
    
Expects temperature in Kelvin, p in bar, and returns s in kJ/kmol/K 

If the optional keyword, diff, is True, then the first derivative of 
entropy is also returned; otherwise it is None.  This is more 
efficient than calculating ds/dT from specific heat separately.
"""
        out = np.zeros_like(T,dtype=float)
        ds = None
        if diff:
            ds = np.zeros_like(T,dtype=float)
        # Loop through the piece-wise temperature ranges
        for index in range(len(self.data['Tlim'])-1):
            # Which elements are in-range?
            I = self._crange(T,index)
            term = 4.
            for c in self.data['C'][index][4:0:-1]:
                if diff:
                    ds[I] = out[I] + T[I]*ds[I]
                out[I] = c/term + T[I]*out[I]
                term -= 1.
            if diff:
                ds[I] = out[I] + T[I]*ds[I]
                ds[I] += self.data['C'][index][0] / T[I]
                ds[I] *= pm.units.const_Ru
            out[I] = T[I]*out[I]\
                    + self.data['C'][index][6]\
                    + self.data['C'][index][0] * np.log(T[I])
        return pm.units.const_Ru * out, ds


    def contents(self):
        """Return the atomic contents dictionary
    C = contents()

Returns a dictionary, C, with keywords that are elements and integer
values representing the number of each present in the species ID.  There
is promise that PYroMat has a data entry for each of the elements named.

For example, the species ID ig.CO2 would return the dictionary
{'C':1, 'O':2}

This is entirely dissimilar from queries to the igmix class, which list
their constituents by their species ID.  These include the collection 
and and the chemical formula, and all constituents MUST have a valid 
species ID.

After the first call to contents(), the dictionary is stored in
the _contents member so that subsequent calls will not result in
redundant string parsing.  
"""
        return self.data['contents'].copy()


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

    #
    # EOS methods
    #
    def p(self,*varg, **kwarg):
        """Pressure from density and temperature
    p(T,d)

Returns the pressure as a function of density and temperature.  Omitted
temperature or pressure will be set to the PYroMat defaults 
config['def_T'] and config['def_p'].

Temperature in      [unit_temperature]
Density in          [unit_matter / unit_volume]
Returns pressure in [unit_pressure]
"""
        p = self._argparse(*varg, pressure=True, **kwarg)
        pm.units.pressure(p, from_units='Pa', inplace=True)
        return p

    
    def d(self,*varg, **kwarg):
        """Density
    d(T,p)
    
Returns density as a function of temperature and/or pressure.  Omitted
temperature or pressure will be set to the PYroMat defaults 
config['def_T'] and config['def_p'].

Temperature in      [unit_temperature]
Pressure in         [unit_pressure]
Returns density in  [unit_matter / unit_volume]
"""
        d = self._argparse(*varg, density=True, **kwarg)
        pm.units.volume(d, from_units='m3', exponent=-1, inplace=True)
        pm.units.matter(d, self.data['mw'], from_units='kmol', inplace=True)
        return d


    def T(self, *varg, **kwarg):
        """Temperature
    T(p,d)
    
Returns temperature as a function of pressure and density.  Omitted
temperature or pressure will be set to the PYroMat defaults 
config['def_T'] and config['def_p'].

Pressure in            [unit_pressure]
Density in             [unit_matter / unit_volume]
Returns temperature in [unit_temperature]
"""
        T = self._argparse(*varg, temperature=True, **kwarg)
        pm.units.temperature_scale(T, from_units='K', inplace=True)
        return T


    def cp(self, *varg, **kwarg):
        """Constant-pressure specific heat
    cp(T)   OR  cp(p=p, d=d)

Accepts any combination of state parameters that permit the calculation
of temperature.  Returns the constant-pressure specific heat.  Missing 
parameters will calculated from PYroMat's default temperature and 
pressure values in config['def_T'] and config['def_p'].

Temperature in      [unit_temperature]
pressure in         [unit_pressure]
density in          [unit_matter / unit_volume]
Specific heat in    [unit_energy / unit_matter / unit_temperature]
"""
        # Prep temperature and the result arrays
        T = self._argparse(*varg, temperature=True, **kwarg)
        out = self._cp(T)
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self.data['mw'], exponent=-1, from_units='kmol', inplace=True)
        pm.units.temperature(out, exponent=-1, from_units='K', inplace=True)
        return out
        
        
        
    def h(self,*varg, **kwarg):
        """Enthalpy
    h(T)   OR  h(p=p, d=d)

Accepts any combination of state parameters that permit the calculation
of temperature.  Returns the enthalpy.  

Missing parameters will calculated from PYroMat's default temperature
and pressure values in config['def_T'] and config['def_p'].

Temperature in      [unit_temperature]
pressure in         [unit_pressure]
density in          [unit_matter / unit_volume]
Returns enthalpy in [unit_energy / unit_matter]
"""

        # Prep temperature and the result arrays
        T = self._argparse(*varg, temperature=True, **kwarg)
        out = self._h(T)[0]
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self.data['mw'], exponent=-1, from_units='kmol', inplace=True)
        return out

    def s(self,*varg, **kwarg):
        """Entropy
    h(T)   OR  h(p=p, d=d)

Accepts any combination of state parameters that permit the calculation
of temperature.  Returns the constant-pressure specific heat.  Missing 
parameters will calculated from PYroMat's default temperature and 
pressure values in config['def_T'] and config['def_p'].

Temperature in  [unit_temperature]
pressure in     [unit_pressure]
density in      [unit_matter / unit_volume]
Returns in      [unit_energy / unit_matter / unit_temperature]
"""

        # Prep temperature and the result arrays
        T,p = self._argparse(*varg, temperature=True, pressure=True, **kwarg)
        out = self._s(T)[0] - pm.units.const_Ru * np.log(p/self.data['pref'])
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self.data['mw'], exponent=-1, from_units='kmol', inplace=True)
        pm.units.temperature(out, from_units='K', exponent=-1, inplace=True)
        return out
        

    def e(self, *varg, **kwarg):
        """Internal energy
    e(T)   OR  e(p=p, d=d)

Accepts any combination of state parameters that permit the calculation
of temperature.  Returns the internal energy.  

Missing parameters will calculated from PYroMat's default temperature 
and pressure values in config['def_T'] and config['def_p'].

Temperature in      [unit_temperature]
pressure in         [unit_pressure]
density in          [unit_matter / unit_volume]
Returns energy in   [unit_energy / unit_matter]
"""

        # Prep temperature and the result arrays
        T = self._argparse(*varg, temperature=True, **kwarg)
        out = self._h(T)[0] - pm.units.const_Ru*T
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self.data['mw'], exponent=-1, from_units='kmol', inplace=True)
        return out


    def cv(self,*varg, **kwarg):
        """Constant-volume specific heat
    cv(T)   OR  cv(p=p, d=d)

Accepts any combination of state parameters that permit the calculation
of temperature.  Returns the constant-volume specific heat.  Missing 
parameters will calculated from PYroMat's default temperature and 
pressure values in config['def_T'] and config['def_p'].

Temperature should be in  [unit_temperature]
pressure should be in     [unit_pressure]
density should be in      [unit_matter / unit_volume]
Returns specific heat in  [unit_energy / unit_matter / unit_temperature]
"""
        # Prep temperature and the result arrays
        T = self._argparse(*varg, temperature=True, **kwarg)
        out = self._cp(T) - pm.units.const_Ru
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self.data['mw'], exponent=-1, from_units='kmol', inplace=True)
        pm.units.temperature(out, exponent=-1, from_units='K', inplace=True)
        return out


    def mw(self,*varg, **kwarg):
        """Molecular weight
    mw(...)
Arguments to molecular mass/weight are ignored.  Molecular mass/weight
is returned as a scalar.

Returns molecular weight in [unit_mass / unit_molar]
"""
        mw = pm.units.mass(self.data['mw'],from_units='kg')
        mw = pm.units.molar(mw,from_units='kmol',exponent=-1)
        return mw

    def R(self,*varg, **kwarg):
        """Ideal gas constant
    R(...)
    
Arguments to ideal gas constant are ignored.  Ideal gas constant is 
returned as a scalar.

Returns R in [unit_energy / unit_temperature / unit_matter]
"""
        R = pm.units.energy(pm.units.const_Ru, from_units='kJ')
        R = pm.units.temperature(R, from_units='K', exponent=-1)
        R = pm.units.matter(R, self.data['mw'], from_units='kmol', exponent=-1)
        return R

    def gam(self, *varg, **kwarg):
        """Specific heat ratio
    gam(T)  OR  gam(p,d)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined.  Ideal gas specific 
heat ratio is not actually a function of p, but it is permitted as an 
argument for cross-compatibility between species' function calls.

Accepts unit_temperature
        unit_pressure
Returns dimensionless
"""
        # Prep temperature and the result arrays
        T = self._argparse(*varg, temperature=True, **kwarg)
        out = self._cp(T)[0]
        out = out / (out - pm.units.const_Ru)
        return out


    def p_s(self,s,T=None):
        """Pressure as a function of entropy
    p = p_s(s)
        or
    p = p_s(s,T)

Returns the pressure as a function of entropy and temperature.

Entropy in          [unit_energy / unit_matter / unit_temperature]
Temperature in      [unit_temperature]
Returns pressure in [unit_pressure]
"""
        s = pm.units.energy(np.asarray(s,dtype=float), to_units='kJ')
        s = pm.units.temperature(s, to_units='K', exponent=-1)
        s = pm.units.matter(s, self.data['mw'], to_units='kmol', exponent=-1)
        if s.ndim==0:
            s = np.reshape(s, (1,))
            
        if T is None:
            T = pm.config['def_T']
        T = pm.units.temperature(np.asarray(T, dtype=float), to_units='K')
        if T.ndim==0:
            T = np.reshape(T, (1,))
            
        s,T = np.broadcast_arrays(s,T)
        
        p0 = self.data['pref']
        s0 = self._s(T)[0]
        p = p0 * np.exp((s0 - s)/pm.units.const_Ru)
        pm.units.pressure(p, from_units='Pa', inplace=True)
        return p


    def T_s(self,s,p=None):
        """Temperature as a function of entropy
    T = T_s(s)
        or
    T = T_s(s,p)

Entropy in              [unit_energy / unit_matter / unit_temperature]
Returns temperature in  [unit_temperature]
"""
        if p is None:
            p = pm.config['def_p']
        p = pm.units.pressure(np.asarray(p, dtype=float), to_units='Pa')
        if p.ndim==0:
            p = np.reshape(p, (1,))
        
        s = pm.units.energy(np.asarray(s, dtype=float), to_units='kJ')
        s = pm.units.matter(s, self.data['mw'], to_units='kmol', exponent=-1)
        s = pm.units.temperature(s, to_units='K', exponent=-1)
        
        s,p = np.broadcast_arrays(s,p)
        # Return the entropy value to standard pressure
        s += pm.units.const_Ru * np.log(p/self.data['pref'])
        
        Ids = np.ones_like(s, dtype=bool)
        T = np.ones_like(s, dtype=float) * (self.data['Tlim'][0] + self.data['Tlim'][1]) * 0.5
        
        self._iter1(self._s, 'T', s, T, Ids, self.data['Tlim'][0], self.data['Tlim'][-1])
        pm.units.temperature_scale(T, from_units='K', inplace=True)
        return T


    def T_h(self,h):
        """Temperature as a function of enthalpy
    T = T_h(h)
        or
    T = T_h(h)

Returns the temperature as a function of enthalpy and pressure

Enthalpy is            [unit_energy / unit_matter]
Returns temperature as [unit_temperature]
"""
        h = pm.units.energy(np.asarray(h, dtype=float), to_units='kJ')
        h = pm.units.matter(h, self.data['mw'], to_units='kmol', exponent=-1)
        if h.ndim==0:
            h = np.reshape(h, (1,))
        
        Ids = np.ones_like(h, dtype=bool)
        T = np.ones_like(h, dtype=float) * (self.data['Tlim'][0] + self.data['Tlim'][1]) * 0.5
        
        self._iter1(self._h, 'T', h, T, Ids, self.data['Tlim'][0], self.data['Tlim'][-1])
        pm.units.temperature_scale(T, from_units='K', inplace=True)
        return T
        
