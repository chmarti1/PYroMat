import pyromat as pm
import numpy as np
######################
##                  ##
##  Mixture class   ##
##                  ##
######################

class igmix(pm.reg.__basedata__):
    """IGMIX  Ideal gas mixture class

The ideal gas mixture is comprised of components that are ideal gases.  
The properties are calculated by calling the property functions of the
constituents and performing the appropriate weighted averages (by mass
or by volume).  

IGMIX objects include functions for converting between the three basic
state variables.  Each can be calculated in terms of the other two.
  T()  temperature      (unit_temperature)
  p()  pressure         (unit_pressure)
  d()  density          (unit_matter / unit_volume)

IGMIX objects offer the following property functions:
  cp() spec. heat       (unit_energy / unit_temperature / unit_matter)
  cv() spec. heat       (unit_energy / unit_temperature / unit_matter)
  e()  internal energy  (unit_energy / unit_matter)
  h()  enthalpy         (unit_energy / unit_matter)
  gam()  spec. heat ratio (dless)
  mw() molecular weight (unit_mass / unit_molar)
  R()  gas constant     (unit_energy / unit_temperature / unit_matter)
  s()  entropy          (unit_energy / unit_temperature / unit_matter)
  X()  mole ratios      (dless)
  Y()  mass ratios      (dless)

There are also routines to invert properties; e.g. calculating 
temperature from enthalpy or from entropy and pressure.
  T_h()  temperature from enthalpy
  T_s()  temperature from entropy and pressure
  p_s()  pressure from entropy and temperature

The Tlim() method returns the intersection of all the supported 
temperature intervals of the constituents.
  Tlim()  temperature limits  (unit_temperature)

Like other species, the atoms() method returns a dictionary of the atomic 
constituents of the mixture.  However, igmix objects can have fractional 
values that indicate the average number of atoms per molecule.
  atoms()  atomic constituents (count/molecule)

For more information on any of these methods, access the in-line 
documentation using Python's built-in "help()" function.
"""

    def __init__(self,*arg,**kwarg):
        # Call the basedata class
        super(igmix,self).__init__(*arg,**kwarg)

        # Initialize the bootstrap flag
        # Initialization has to be split into two phases.  Calculation
        # of mean mixture properties cannot be completed until all other
        # members of the collection have been loaded.  The _bootsrap()
        # method is responsible for completing the process, and the _bs
        # flag indicates whether it has already been completed.
        self._bs = False
        # Initialize the static molar and mass fractions
        self._x = {}
        self._y = {}
        # Initialize a total mass and molar tally
        total_x = 0.
        total_y = 0.
        # Initialize the mean molecular weight
        self._mw = 0.
        # Initialize the mean reference pressure for entropy
        self._pref_bar = 0.
        # Initialize temperature limits
        self._Tlim = [float('-inf'), float('inf')]
        
        
    def _bootstrap(self):
        """Calculates internal parameters that are essential for the property functions
This operation cannot be completed by __init__ at load time because 
there is no no way to ensure that all of the constituent species have
already been loaded.  Instead, _boostrap() is called by the property
methods to be certain the relevant parameters have been calculated.  If
the _bs member flag has already been set, this method returns 
immediately.

Attribute   Description
_x          Dictionary of mole fractions
_y          Dictionary of mass fractions
_mw         Effective mean molecular mass (weight) in kg/kmol
_pref_bar   Effective log-mean reference pressure in bar
_Tlim       Lower and upper temperature limits of the most restrictive
            constintuent gas data in Kelvin
"""
        if self._bs:
            return
        
        self._bs = True
        
        self._x = {}
        self._y = {}
        # Initialize a total mass and molar tally
        total_x = 0.
        total_y = 0.
        # Initialize the mean molecular weight
        self._mw = 0.
        # Initialize the mean reference pressure for entropy
        self._pref_bar = 0.
        # Initialize temperature limits
        self._Tlim = [float('-inf'), float('inf')]

        for ss,qty in self.data['contents'].items():
            spec = pm.dat.data.get(ss)
            
            # Test for basic data integrity
            if spec is None:
                raise pm.utility.PMDataError('IGMIX: %s contains a species not in the collection: %s'%(self.data['id'], ss))
            # IG stores reference pressure as a member and in bar
            elif isinstance(spec, pm.reg.registry['ig']):
                spec_pref = spec._pref_bar
            # IG2 stores reference pressure in Pa and in the data dictionary
            elif isinstance(spec, pm.reg.registry['ig2']):
                spec_pref = spec.data['pref']/1e5
            else:
                raise pm.utility.PMDataError('IGMIX: %s contains a species class that is not supported: %s'%(self.data['id'], repr(spec)))
            
            spec_mw = spec.data['mw']
            spec_Tmin = spec.data['Tlim'][0]
            spec_Tmax = spec.data['Tlim'][-1]
            
            if self.data['bymass']:
                spec_y = qty
                spec_x = qty / spec_mw
            else:
                spec_y = qty * spec_mw
                spec_x = qty
                
            total_y += spec_y
            total_x += spec_x
            
            # Update parameters
            self._y[ss] = spec_y
            self._x[ss] = spec_x
            self._mw += spec_mw * spec_x
            self._pref_bar += np.log(spec_pref) * spec_x
            self._Tlim[0] = max(spec_Tmin, self._Tlim[0])
            self._Tlim[1] = min(spec_Tmax, self._Tlim[1])
            
        # Normalize weight by total molar contents
        self._mw /= total_x
        # Normalize and rescale the reference pressure by the total molar contents
        self._pref_bar = np.exp(self._pref_bar / total_x)
            
        # Loop through one more time to normalize by the mass and molar
        # totals.
        for ss in self._x:
            self._x[ss] /= total_x
            self._y[ss] /= total_y
            
        


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
            d = pm.units.matter(np.asarray(d,dtype=float), self.data['mw'], to_units='kmol')
            pm.units.volume(d, to_units='m3', exponent=-1, inplace=True)
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
        """Calculates the specific heat of the mixture from contents
    cp = _cp(T)
    
T must be a numpy array in Kelvin
"""
        # Initialize a result array
        out = np.zeros_like(T,dtype=float)
        
        # Calculate the property directly from the subordinate inner 
        # routines.  They use molar units by default.
        for ss,x in self._x.items():
            # Retrieve the species
            out += x * pm.dat.data[ss]._cp(T)[0]
            
        return out
        
        
    def _s(self, T, diff=False):
        """Calculates entropy at the reference pressure
    s0,sT = _s(T, diff=True)

s0 is the temperature contribution to entropy in kJ/kmol/K
sT is the derivative of s0 with respect to temperature.
T must be a numpy array in Kelvin.
"""
        # Initialize a result array
        s0 = np.zeros_like(T,dtype=float)
        sT = None
        if diff:
            sT = np.zeros_like(T,dtype=float)
        
        # Calculate the property directly from the subordinate inner 
        # routines.  They use molar units by default.
        for ss,x in self._x.items():
            # Retrieve the species data
            spec_s0,spec_sT = pm.dat.data[ss]._s(T, diff=diff)
            
            s0 += x * spec_s0
            if diff:
                sT += x * spec_sT
            
        return s0,sT
        
        
    def _h(self, T, diff=False):
        """Calculates enthalpy and its derivative
    h,hT = _h(T, diff=True)
    
T must be a numpy array in Kelvin
h is in kJ/kmol
hT is in kJ/kmol/K
"""
        # Initialize a result array
        h = np.zeros_like(T,dtype=float)
        hT = None
        if diff:
            hT = np.zeros_like(T,dtype=float)
        
        
        # Calculate the property directly from the subordinate inner 
        # routines.  They use molar units by default.
        for ss,x in self._x.items():
            # Retrieve the species data
            spec_h,spec_hT = pm.dat.data[ss]._h(T, diff=diff)
            
            h += x * spec_h
            if diff:
                hT += x * spec_hT
            
        return h,hT
        
        
    def atoms(self):
        """Return a dictionary specifying the chemical composition of the substance.
    aa = atoms()
    
This dictionary will be a sum of the atoms of the constituent gases weighted
by the mixture composition by volume.  This results in atom quantities that
are floating point instead of integer.  The number indicates the moles of
each atom contained per mole of mixed gas.
"""
        self._bootstrap()
        atoms = {}
        for ss,x in self._x.items():
            spec = pm.dat.data.get(ss)
            ss_atoms = spec.data.get('atoms')
            if ss_atoms is None:
                raise pm.utility.PMDataError('Mixture contains a species with no atomic data: ' + self.data['id'] + '-->' + spec.data['id'])
            for atom,qty in ss_atoms.items():
                if atom not in atoms:
                    atoms[atom] = 0
                atoms[atom] += qty * x
        return atoms
        
        
       
    def Tlim(self):
        """Temperature limits
    (Tmin, Tmax) = Tlim()
Returns the temperature limits on the ig data set.

Accepts None
Returns unit_temperature
"""
        self._bootstrap()
        T = np.array(self._Tlim, dtype=float)
        np.temperature_scale(T, from_units='K', inplace=True)
        return tuple(T)


    #
    #   Pressure, temperature, and density state functions
    #
    def d(self,*varg, **kwarg):
        """Density
    d(T,p)
        OR
    d(d=d,...)
    
Calculates density from any two of the T,p,d triple.  If any arguments 
are missing, they are replaced with their default values are replaced
by the defaults config['def_T'] or config['def_p'].

Accepts     Temperature [unit_temperature]
            Pressure    [unit_pressure]
Returns     Density     [unit_matter / unit_volume]
"""
        self._bootstrap()
        d = self._argparse(*varg, density=True, **kwarg)
        pm.units.matter(d, self._mw, from_units='kmol', inplace=True)
        pm.units.volume(d, from_units='m3', inplace=True, exponent=-1)
        return d
        
        
    def p(self, *varg, **kwarg):
        """Pressure
    d(T,p)
        OR
    d(p=p,...)
    
Calculates pressure from any two of the T,p,d triple.  If any arguments 
are missing, they are replaced with their default values are replaced
by the defaults config['def_T'] or config['def_p'].

Accepts     Temperature [unit_temperature]
            Density      [unit_matter / unit_volume]
Returns     Pressure     [unit_pressure]
"""
        self._bootstrap()
        p = self._argparse(*varg, pressure=True, **kwarg)
        pm.units.pressure(p, from_units='Pa', inplace=True)
        return p
        
        
    def T(self, *varg, **kwarg):
        """Temperature
    T(p=p,d=d)
        OR
    T(T=T,...)
    
Calculates temperature from any two of the T,p,d triple.  If any 
arguments are missing, they are replaced with their default values are 
replaced by the defaults config['def_T'] or config['def_p'].

Accepts     Pressure     [unit_pressure]
            Density      [unit_matter / unit_volume]
Returns     Temperature  [unit_temperature]
"""
        self._bootstrap()
        p = self._argparse(*varg, pressure=True, **kwarg)
        pm.units.pressure(p, from_units='Pa', inplace=True)
        return p

    #
    # Class property functions
    #
    def cp(self,*varg, **kwarg):
        """Constant-pressure specific heat
    cp(T)
        OR
    cp(T,p)
        OR
    cp(p=p, d=d)
    
Calculates the specific heat of the mixture from the specific heats of
the components.

Accepts:    Temperature [unit_temperature]
            Pressure    [unit_pressure]
            Density     [unit_matter / unit_volume]
Returns:    Spec. Heat  [unit_energy / unit_temperature / unit_matter]
"""
        self._boostrap()
        # Parse the arguments to isolate temperature in K
        T = self._argparse(*varg, temperature=True, **kwarg)

        out = self._cp(T)
        
        # Convert output
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self._mw, from_units='kmol', inplace=True, exponent=-1)
        pm.units.temperature(out, from_units='K', inplace=True, exponent=-1)
        return out


    def cv(self,T=None,p=None):
        """Constant-volume specific heat
    cp(T)
        OR
    cp(T,p)
        OR
    cp(p=p, d=d)
    
Calculates the specific heat of the mixture from the specific heats of
the components.  

Accepts:    Temperature [unit_temperature]
            Pressure    [unit_pressure]
            Density     [unit_matter / unit_volume]
Returns:    Spec. Heat  [unit_energy / unit_temperature / unit_matter]
"""
        self._boostrap()
        # Parse the arguments to isolate temperature in K
        T = self._argparse(*varg, temperature=True, **kwarg)
        
        out = self._cp(T) - pm.units.const_Ru
        
        # Convert output
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self._mw, from_units='kmol', inplace=True, exponent=-1)
        pm.units.temperature(out, from_units='K', inplace=True, exponent=-1)


    def h(self,*varg, **kwarg):
        """Enthalpy
    h(T)
        OR
    h(T,p)
        OR
    h(p=p, d=d)
    
Calculates the enthalpy at the specified state from the mixture 
components' entropies.  Missing state parameters will be assigned their
defaults from either config['def_T'] or config['def_p'].

Accepts:    Temperature [unit_temperature]
            Pressure    [unit_pressure]
            Density     [unit_matter / unit_volume]
Returns:    Enthalpy    [unit_energy / unit_matter]
"""
        self._bootstrap()
        T = self._argparse(*varg, temperature=True, **kwarg)
        out = np.zeros_like(T,dtype=float)
        for ss,x in self._x.items():
            out += x*self._h(T)[0]
        
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self._mw, from_units='kmol', inplace=True, exponent=-1)
        return out
        

    def e(self, *varg, **kwarg):
        """Internal energy
    e(T)
        OR
    e(T,p)
        OR
    e(p=p, d=d)
    
Calculates the internal energy at the specified state from the mixture 
components' entropies.  Missing state parameters will be assigned their
defaults from either config['def_T'] or config['def_p'].

Accepts:    Temperature [unit_temperature]
            Pressure    [unit_pressure]
            Density     [unit_matter / unit_volume]
Returns:    Int. Energy [unit_energy / unit_matter]
"""
        self._bootstrap()
        T = self._argparse(*varg, temperature=True, **kwarg)
        out = np.zeros_like(T,dtype=float)
        for ss,x in self._x.items():
            out += x*self._h(T)[0]
        out -= T * pm.units.const_Ru
        
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self._mw, from_units='kmol', inplace=True, exponent=-1)
        return out


    def mw(self,*varg, **kwarg):
        """Molecular weight (more correctly mass)
    mw(...)

Arguments are ignored.  The molecular mass of the mixture is calculated
as a weighted average of the molecular masses of the constituents as is
appropriate for use in the ideal gas law and for rescaling between molar
and mass-based intensive properties.

Accepts:    -
Returns:    Molecular mass [unit_mass / unit_molar]
"""
        self._bootstrap()
        out = pm.units.mass(self._mw, from_units='kg')
        out = pm.units.molar(out, from_units='kmol')
        return out


    def s(self, *varg, **kwarg):
        """Entropy
    s(T)
        OR
    s(T,p)
        OR
    s(p=p, d=d)
    
Calculates the entropy at the specified state from the mixture 
components' entropies.  Missing state parameters will be assigned their
defaults from either config['def_T'] or config['def_p'].

Accepts:    Temperature [unit_temperature]
            Pressure    [unit_pressure]
            Density     [unit_matter / unit_volume]
Returns:    Entropy     [unit_energy / unit_matter / unit_temperature]
"""
        self._bootstrap()
        T,p = self._argparse(*varg, temperature=True, pressure=True, **kwarg)
        s = self._s(T)[0] - pm.units.const_Ru * np.log(p/(1e5*self._pref_bar))
        pm.units.energy(s, from_units='kJ', inplace=True)
        pm.units.matter(s, self._mw, from_units='kmol', inplace=True, exponent=-1)
        pm.units.temperature(s, from_units='K', inplace=True, exponent=-1)
        return s


    def gam(self,*varg, **kwarg):
        """Specific heat ratio (gamma)
    gam(T)
        OR
    gam(T,p)
        OR
    gam(p=p, d=d)
    
Calculates the specific heat ratio (cp / cv) for the mixture.  Missing 
state parameters will be assigned their defaults from either 
config['def_T'] or config['def_p'].

Accepts:    Temperature [unit_temperature]
            Pressure    [unit_pressure]
            Density     [unit_matter / unit_volume]
Returns:    Gamma       [d-less]
"""
        self._bootsrap()
        T,p = self._argparse(*varg, temperature=True, **kwarg)
        out = self._cp(T)
        return out / (out - pm.units.const_Ru)


    def X(self):
        """Return a dictionary expressing the mixture composition by volume
    x = mix.X()

Where x is a dictionary with keys corresponding to the species in
the mixture and values indicating their respective mole fraction
in the mixture."""
        self._bootstrap()
        # Return the dictionary
        return self._x.copy()



    def Y(self):
        """Return a dictionary expressing the mixture composition by mass
    y = mix.Y()

Where y is a dictionary with keys corresponding to the species in
the mixture and values indicating their respective mass fraction
in the mixture."""
        self._bootstrap()
        # Return the dictionary
        return self._y.copy()


    def T_s(self,s,p=None, d=None, debug=False):
        """Temperature as a function of entropy
    T = T_s(s)
        or
    T = T_s(s,p)
        or
    T = T_s(s,d)

Accepts unit_energy / unit_matter / unit_temperature
        unit_pressure
        unit_matter / unit_volume
Returns unit_temperature
"""
        self._bootstrap()
        if p is None and d is None:
            p = pm.config['def_p']
                    
        s = pm.units.energy(np.asarray(s, dtype=float), to_units='kJ')
        s = pm.units.matter(s, self._mw, to_units='kmol', exponent=-1)
        s = pm.units.temperature(s, to_units='K', exponent=-1)
        if s.ndim == 0:
            s = np.reshape(s, (1,))
            
        # If isobaric
        if p is not None:
            p = pm.units.pressure(np.asarray(p, dtype=float), to_units='bar')
            if p.ndim==0:
                p = np.reshape(p, (1,))
            
            s,p = np.broadcast_arrays(s,p)
            # Adjust s by the pressure term
            s += pm.units.const_Ru * np.log(p/self._pref_bar)
            
            I = np.ones_like(s, dtype=bool)
            T = np.full_like(s, 0.5*(self._Tlim[0]+self._Tlim[-1]))
            self._iter1(self._s, 'T', s, T, I, self._Tlim[0], self._Tlim[-1], verbose=debug)
        # If isochoric
        else:
            d = pm.units.matter(np.asarray(d, dtype=float),
                    self._mw, to_units='kmol')
            d = pm.units.volume(d, to_units='m3', exponent=-1)
            
            s,d = np.broadcast_arrays(s,d)
            
            R = pm.units.const_Ru
            # Define a custom iterator function
            def fn(T,d,diff):
                sd = 0.
                s,sT = self._s(T,diff)
                
                s -= R*np.log(d * R * T / (self._pref_bar * 1e2))
                if diff:
                    sT -= R / T
                    sd = -R / d

                return s,sT,sd
            
            I = np.ones_like(s, dtype=bool)
            T = np.full_like(s, 0.5*(self._Tlim[0]+self._Tlim[-1]))
            self._iter1(fn, 'T', s, T, I, self._Tlim[0], self._Tlim[-1], param={'d':d}, verbose=debug)
            
        pm.units.temperature_scale(T, from_units='K')
        return T

    def T_s_(self,s,p=None):
        """Temperature as a function of entropy
    T = T_s(s)
        or
    T = T_s(s,p)

Accepts unit_energy / unit_matter / unit_temperature
        unit_pressure
Returns unit_temperature
"""
        self._bootstrap()
        s = pm.units.energy(np.asarray(s,dtype=float), to_units='kJ')
        s = pm.units.matter(s, self._mw, to_units='kmol', exponent=-1)
        s = pm.units.temperature(s, to_units='K', exponent=-1)
        if s.ndim==0:
            s = np.reshape(s, (1,))
            
        if p is None:
            p = pm.config['def_p']
        p = pm.units.pressure(np.asarray(p,dtype=float), to_units='Pa')
        
        s,p = np.broadcast_arrays(s,p)
        
        # Adjust entropy to the reference pressure
        s += pm.units.const_Ru * np.log(p / (1e5*self._pref_bar))
        
        T = np.full_like(s, 0.5*np.sum(self._Tlim))
        I = np.ones_like(s, dtype=bool)
        self._iter1(self._s, 'T', s, T, I, self._Tlim[0], self._Tlim[1])
        pm.units.temperature_scale(T, from_units='K', inplace=True)
        return T


    def T_h(self,h,p=None):
        """Temperature as a function of enthalpy
    T = T_h(h)
        or
    T = T_h(h,p)

Returns the temperature as a function of enthalpy and pressure

Accepts unit_energy / unit_matter
        unit_pressure
Returns unit_temperature
"""
        self._bootstrap()
        h = pm.units.energy(np.asarray(h,dtype=float), to_units='kJ')
        h = pm.units.matter(h, self._mw, to_units='kmol', exponent=-1)
        if h.ndim==0:
            h = np.reshape(h, (1,))

        # p will be ignored, so don't bother converting its units
        if p is None:
            p = pm.config['def_p']
        h,p = np.broadcast_arrays(h,p)
        
        T = np.full_like(h, 0.5*np.sum(self._Tlim))
        I = np.ones_like(h, dtype=bool)
        self._iter1(self._h, 'T', h, T, I, self._Tlim[0], self._Tlim[1])
        pm.units.temperature_scale(T, from_units='K', inplace=True)
        return T
        
        
    def p_s(self, s, T=None):
        """Pressure as a function of entropy
    p = p_s(s)
        OR
    p = p_s(s,T)

Returns the pressure as a function of entropy and temperature.

Accepts unit_energy / unit_matter / unit_temperature
        unit_temperature
Returns unit_pressure
"""
        self._bootstrap()
        s = pm.units.energy(np.asarray(s, dtype=float), to_units='kJ')
        s = pm.units.matter(s, self._mw, to_units='kmol', exponent=-1)
        s = pm.units.temperature(s, to_units='K', exponent=-1)
        if s.ndim==0:
            s = np.reshape(s,(1,))
        
        if T is None:
            T = pm.config['def_T']
        T = pm.units.temperature_scale(np.asarray(T, dtype=float), to_units='K')
        if T.ndim==0:
            T = np.reshape(T, (1,))
        
        s,T = np.broadcast_arrays(s,T)
        
        p = self._pref_bar * np.exp((self._s(T)[0]-s)/pm.units.const_Ru)
        
        pm.units.pressure(p, from_units='bar', inplace=True)
        return p
