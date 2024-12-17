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

** Available Property Methods **
IGMIX includes methods for calculating the same properties as ig and ig2
with a few extra unique to mixtures:
  T()  temperature      (unit_temperature)
  p()  pressure         (unit_pressure)
  d()  density          (unit_matter / unit_volume)
  cp() spec. heat       (unit_energy / unit_temperature / unit_matter)
  cv() spec. heat       (unit_energy / unit_temperature / unit_matter)
  e()  internal energy  (unit_energy / unit_matter)
  h()  enthalpy         (unit_energy / unit_matter)
  gam()  spec. heat ratio (dless)
  s()  entropy          (unit_energy / unit_temperature / unit_matter)
  X()  mole ratios      (dless)
  Y()  mass ratios      (dless)


** Other Properties **
The Tlim() method returns the intersection of all the supported 
temperature intervals of the constituents.
  Tlim()  temperature limits  (unit_temperature)

Molecular weight and ideal gas constant for a mixture are calculated as
an appropriately weighted average.  The molecular weight of a mixture is
the average mass per mole of molecules.
  mw() molecular weight (unit_mass / unit_molar)
  R()  gas constant     (unit_energy / unit_temperature / unit_matter)

Like other species, the atoms() method returns a dictionary of the atomic 
constituents of the mixture.  However, igmix objects can have fractional 
values that indicate the average number of atoms per molecule.
  atoms()  atomic constituents (count/molecule)

** Depreciated Methods **
Since version 2.2.0, the flexible interface has allowed passing enthalpy
and entropy directly to the T() and p() methods, so these methods are
no longer needed.  They will be removed with the next major revision 
change, so new code should no longer use them.
  T_h()  temperature from enthalpy
  T_s()  temperature from entropy and pressure
  p_s()  pressure from entropy and temperature

For more information on any of these methods, access the in-line 
documentation using Python's built-in "help()" function.
"""

    def __init__(self,*arg,**kwarg):
        # Call the basedata class
        super(igmix,self).__init__(*arg,**kwarg)

        # Initialize the bootstrap flag
        # Initialization has to be split into two phases.  Calculation
        # of mean mixture properties cannot be completed until all other
        # members of the collection have been loaded.  The _bootstrap()
        # method is responsible for completing the process, and the _bs
        # flag indicates whether it has already been completed.
        self._bs = False
        # Initialize the static molar and mass fractions
        self._x = {}
        self._y = {}
        # Initialize the mean molecular weight
        self._mw = 0.
        # Initialize the mean reference pressure for entropy
        self._pref_pa = 0.
        # Initialize temperature limits
        self._Tlim = [float('-inf'), float('inf')]
        
        
    def _bootstrap(self):
        """Calculates internal parameters that are essential for the property functions
This operation cannot be completed by __init__ at load time because 
there is no no way to ensure that all of the constituent species have
already been loaded.  Instead, _bootstrap() is called by the property
methods to be certain the relevant parameters have been calculated.  If
the _bs member flag has already been set, this method returns 
immediately.

Attribute   Description
_x          Dictionary of mole fractions
_y          Dictionary of mass fractions
_mw         Effective mean molecular mass (weight) in kg/kmol
_pref_bar   Effective log-mean reference pressure in bar used for entropy
_smix       The enthalpy of mixing 
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
        self._pref_pa = 0.
        # Initialize temperature limits
        self._Tlim = [float('-inf'), float('inf')]

        for ss,qty in self.data['contents'].items():
            spec = pm.dat.data.get(ss)
            
            # Test for basic data integrity
            if spec is None:
                raise pm.utility.PMDataError('IGMIX: %s contains a species not in the collection: %s'%(self.data['id'], ss))
            # IG stores reference pressure as a member and in bar
            elif isinstance(spec, pm.reg.registry['ig']):
                spec_pref = spec._pref_pa
            # IG2 stores reference pressure in Pa and in the data dictionary
            elif isinstance(spec, pm.reg.registry['ig2']):
                spec_pref = spec.data['pref']
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
            self._pref_pa += np.log(spec_pref) * spec_x
            self._Tlim[0] = max(spec_Tmin, self._Tlim[0])
            self._Tlim[1] = min(spec_Tmax, self._Tlim[1])
            
        # Normalize weight by total molar contents
        self._mw /= total_x
        # Normalize and rescale the reference pressure by the total molar contents
        self._pref_pa = np.exp(self._pref_pa / total_x)
        
        # Loop through one more time to normalize by the mass and molar
        # totals and calculate the entropy of mixing
        self._smix = 0.
        for ss in self._x:
            self._x[ss] /= total_x
            self._y[ss] /= total_y
            # Enthalpy of mixing
            self._smix -= self._x[ss] * np.log(self._x[ss])
        self._smix *= pm.units.const_Ru
            
        


    def _argparse(self, *varg, **kwarg):
        """Parse the arguments supplied to an IGMIX property method
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
            kwarg['d'] = pm.units.matter(value, self._mw, to_units='kmol')
        if 'v' in kwarg:
            # Convert and replace with d at the same time
            value = pm.units.volume(kwarg['v'], to_units='m3')
            kwarg['d'] = 1./pm.units.matter(value, self._mw, to_units='kmol', exponent=-1)
        if 'h' in kwarg:
            value = kwarg['h']
            value = pm.units.energy(value, to_units='kJ')
            value = pm.units.matter(value, self._mw, to_units='kmol', exponent=-1)
            T = np.full(value.shape, 0.5*(self._Tlim[0] + self._Tlim[-1]))
            I = np.ones(value.shape, dtype=bool)
            self._iter1(self._h, 'T', value, T, I, self._Tlim[0], self._Tlim[-1])
            kwarg['T'] = T
        if 'e'  in kwarg:
            value = kwarg['e']
            value = pm.units.energy(value, to_units='kJ')
            value = pm.units.matter(value, self._mw, to_units='kmol', exponent=-1)
            Tlow = self._Tlim[0]
            Thigh = self._Tlim[-1]
            T = np.full(value.shape, 0.5*(self._Tlim[0] + self._Tlim[-1]))
            I = np.ones(value.shape, dtype=bool)
            self._iter1(self._e, 'T', value, T, I, self._Tlim[0], self._Tlim[-1])
            kwarg['T'] = T
        if 's' in kwarg:
            value = kwarg['s']
            value = pm.units.energy(value, to_units='kJ')
            value = pm.units.matter(value, self._mw, to_units='kmol', exponent=-1)
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
                T = np.full_like(s, 0.5*(self._Tlim[0] + self._Tlim[-1]))
                I = np.ones_like(s,dtype=bool)
                self._iter1(self._sditer, 'T', s, T, I, self._Tlim[0], self._Tlim[-1], param={'d':d})
            # If pressure is specified
            elif 'p' in kwarg:
                s,p = np.broadcast_arrays(kwarg['s'], kwarg['p'])
                # adjust entropy to the reference pressure
                s += pm.units.const_Ru * np.log(p / self._pref_pa)
                T = np.full_like(s, 0.5*(self._Tlim[0] + self._Tlim[-1]))
                I = np.ones_like(s,dtype=bool)
                self._iter1(self._s, 'T', s, T, I, self._Tlim[0], self._Tlim[-1])
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
        I = np.logical_or(T < self._Tlim[0], T > self._Tlim[-1])
        if I.all():
            raise pm.utility.PMParamError('All of the specified states were out-of-bounds.  '  
                    'Legal temperatures for {} are between {} and {} Kelvin.'.format(self.data['id'], self._Tlim[0], self._Tlim[-1]))
        elif I.any():
            T[I] = pm.config['def_oob']
            pm.utility.print_warning('Some of the states were out of bounds - setting to config[\'def_oob\'].  '
                    'Legal temperatures for {} are between {} and {} Kelvin.'.format(self.data['id'], self._Tlim[0], self._Tlim[-1]))
        return T,p,d

    
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

    def _sditer(self, T, d, diff=1):
        s,sT = self._s(T, diff)
        R = 1000 * pm.units.const_Ru
        s -= pm.units.const_Ru * np.log(d * R * T / self._pref_pa)
        if diff:
            sT -= pm.units.const_Ru/T
        else:
            sT = None
        return s,sT

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
            out += x * pm.dat.data[ss]._cp(T)
            
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
        # Add the entropy of mixing
        s0 += self._smix
        return s0,sT
        
        
    def _g(self, T, diff=False):
        """Calculates Gibbs energy at the reference pressure
    g0,gT = _g(T, diff=True)

g0 is the Gibbs energy in kJ/kmol at the reference pressure
gT is the derivative of g0 with respect to temperature.
T must be a numpy array in Kelvin.
"""
        # Initialize a result array
        g0 = np.zeros_like(T,dtype=float)
        gT = None
        if diff:
            gT = np.zeros_like(T,dtype=float)
        
        # Calculate the property directly from the subordinate inner 
        # routines.  They use molar units by default.
        for ss,x in self._x.items():
            # Retrieve the species data
            spec_g0,spec_gT = pm.dat.data[ss]._g(T, diff=diff)
            
            g0 += x * spec_g0
            if diff:
                gT += x * spec_gT
        # Deduct the entropy of mixing
        g0 -= T * self._smix
        return g0,gT
        
        
    def _f(self, T, diff=False):
        """Calculates free energy at the reference pressure
    f0,fT = _f(T, diff=True)

f0 is the free (Helmholtz) energy in kJ/kmol at the reference pressure
fT is the derivative of g0 with respect to temperature.
T must be a numpy array in Kelvin.
"""
        # Initialize a result array
        f0 = np.zeros_like(T,dtype=float)
        fT = None
        if diff:
            fT = np.zeros_like(T,dtype=float)
        
        # Calculate the property directly from the subordinate inner 
        # routines.  They use molar units by default.
        for ss,x in self._x.items():
            # Retrieve the species data
            spec_f0,spec_fT = pm.dat.data[ss]._f(T, diff=diff)
            
            f0 += x * spec_f0
            if diff:
                fT += x * spec_fT
        # Deduct the entropy of mixing
        f0 -= T * self._smix
        return f0,fT
        
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
        
    def _e(self, T, diff=False):
        """Calculates internal energy and its derivative
    e,eT = _e(T, diff=True)
    
T must be a numpy array in Kelvin
e is in kJ/kmol
hT is in kJ/kmol/K
"""
        h,hT = self._h(T, diff)
        R = pm.units.const_Ru
        h -= R*T
        if diff:
            hT -= R
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
        pm.units.temperature_scale(T, from_units='K', inplace=True)
        return tuple(T)


    #
    #   Pressure, temperature, and density state functions
    #
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
        self._bootstrap()
        Ru = pm.units.const_Ru
        T,p,d = self._argparse(*varg, **kwarg)
        # Make sure we have both pressure and density
        if d is None:
            d = p / (1000 * Ru * T)
        pm.units.matter(d, self._mw, from_units='kmol', inplace=True)
        pm.units.volume(d, from_units='m3', inplace=True, exponent=-1)
        return d

    def v(self, *varg, **kwarg):
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
        return 1. / self.d(*varg, **kwarg)


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
        self._bootstrap()
        Ru = pm.units.const_Ru
        T,p,d = self._argparse(*varg, **kwarg)
        if p is None:
            p = 1000 * d * Ru * T
        pm.units.pressure(p, from_units='Pa', inplace=True)
        return p
        
        
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
        self._bootstrap()
        # T is a special case.  If there is only one parameter given, 
        # the default should not be T=def_T, so we need to override 
        # the behavior of _argparse. 
        if len(varg) + len(kwarg) == 1: 
            kwarg['p'] = pm.config['def_p'] 
        T,_,_ = self._argparse(*varg, **kwarg)
        pm.units.temperature_scale(T, from_units='K', inplace=True)
        return T

    #
    # Class property functions
    #
    
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
        self._bootstrap()
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
        scale = pm.units.matter(1., self._mw, from_units='kmol')
        scale = pm.units.volume(scale, from_units='m3', exponent=-1)
        out['d'] = scale * d
        out['v'] = 1./out['d']
        scale = pm.units.energy(1., from_units='kJ')
        scale = pm.units.matter(scale, self._mw, from_units='kmol', exponent=-1)
        out['h'] = h * scale
        out['e'] = e * scale
        out['gam'] = gam
        scale = pm.units.temperature(scale, from_units='K', exponent=-1)
        out['s'] = s * scale
        out['cp'] = cp * scale
        out['cv'] = cv * scale
        return out

    
    def cp(self,*varg, **kwarg):
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
        self._bootstrap()
        # Parse the arguments to isolate temperature in K
        T,_,_ = self._argparse(*varg, **kwarg)

        out = self._cp(T)
        
        # Convert output
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self._mw, from_units='kmol', inplace=True, exponent=-1)
        pm.units.temperature(out, from_units='K', inplace=True, exponent=-1)
        return out


    def cv(self,*varg,**kwarg):
        """Constant-volume specific heat
    cv(...)

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
        self._bootstrap()
        # Parse the arguments to isolate temperature in K
        T,_,_ = self._argparse(*varg, **kwarg)
        
        out = self._cp(T) - pm.units.const_Ru
        
        # Convert output
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self._mw, from_units='kmol', inplace=True, exponent=-1)
        pm.units.temperature(out, from_units='K', inplace=True, exponent=-1)
        return out


    def h(self, *varg, **kwarg):
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
        self._bootstrap()
        T,_,_ = self._argparse(*varg, **kwarg)
        out = self._h(T)[0]
        
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self._mw, from_units='kmol', inplace=True, exponent=-1)
        return out
        
        
    def g(self, *varg, **kwarg):
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

Returns Gibbs energy in unit_energy / unit_matter
"""
        self._bootstrap()
        T,p,d = self._argparse(*varg, **kwarg)
        if p is None:
            p = 1000 * pm.units.const_Ru * d * T
        out = self._g(T)[0] + T * pm.units.const_Ru * np.log(p/self._pref_pa)
        
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self._mw, from_units='kmol', inplace=True, exponent=-1)
        return out

    def e(self, *varg, **kwarg):
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
        self._bootstrap()
        T,_,_ = self._argparse(*varg, **kwarg)
        out = self._e(T)[0]
        
        pm.units.energy(out, from_units='kJ', inplace=True)
        pm.units.matter(out, self._mw, from_units='kmol', inplace=True, exponent=-1)
        return out

    def f(self, *varg, **kwarg):
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

Returns free energy in unit_energy / unit_matter
"""
        self._bootstrap()
        T,p,d = self._argparse(*varg, **kwarg)
        if p is None:
            p = 1000 * pm.units.const_Ru * d * T
        out = self._f(T)[0] + T * pm.units.const_Ru * np.log(p/self._pref_pa)
        
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
        out = pm.units.molar(out, from_units='kmol', exponent=-1)
        return out


    def R(self, *varg, **kwarg):
        """Ideal gas constant
    R(...)

Ignores the arguments are returns the gas constant as
unit_energy / unit_matter / unit_temperature
"""
        self._bootstrap()
        R = pm.units.energy(pm.units.const_Ru, from_units='J')
        R = pm.units.temperature(R, from_units='K', exponent=-1)
        R = pm.units.matter(R, self._mw, from_units='mol',
                                exponent=-1)
        return R


    def s(self, *varg, **kwarg):
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
        self._bootstrap()
        T,p,d = self._argparse(*varg, **kwarg)
        if p is None:
            p = 1000 * pm.units.const_Ru * d * T
        s = self._s(T)[0] - pm.units.const_Ru * np.log(p/self._pref_pa)
        pm.units.energy(s, from_units='kJ', inplace=True)
        pm.units.matter(s, self._mw, from_units='kmol', inplace=True, exponent=-1)
        pm.units.temperature(s, from_units='K', inplace=True, exponent=-1)
        return s


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
        self._bootstrap()
        T,_,_ = self._argparse(*varg, **kwarg)
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
