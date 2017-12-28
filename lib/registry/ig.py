import pyromat as pyro
import numpy as np
import os
import re



class ig(pyro.reg.__basedata__):
    """Ideal gas class using the Shomate equation of state.
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

    def __init__(self,*arg,**kwarg):
        super(self.__class__,self).__init__(*arg,**kwarg)

        # Important constants
        self._pref_bar = 1.0
        # Initialize the species contents dictionary
        self._contents = None


    def _crange(self, T):
        """Return the temperature range index and raise a meaningful exception
if it is out of range."""
        Tlim = self.data['Tlim']
        if T<Tlim[0]:
            raise pyro.utility.PMParamError(
                'Temperature, %f, is out of range for %s (%f to %f)'%(
                T, self.data['id'], Tlim[0], Tlim[-1]))
        for index in range(1,len(self.data['Tlim'])):
            if T<=Tlim[index]:
                return index-1
        raise pyro.utility.PMParamError(
            'Temperature, %f, is out of range for %s (%f to %f)'%(
            T, self.data['id'], Tlim[0], Tlim[-1]))


    def _invT(self, value, prop, p=None):
        """Inverse property function
    T = _invT(value, prop, dprop, p=None)

This is an internal service method used to provide T_s and T_h.
Iterate on values of temperature so that the property method, prop,
returns value.

The algorithm uses a modified Newton iteration.  The property derivative
is approximated numerically instead of relying on explicit derivatives.
The numerical cost of explicitly evaluating a function's derivatives is 
approximately equal to numerical approximation in 1D, but this 
drastically simplifies the implementation.

value   the property value to obtain
prop    the property method to be inverted
p       the pressure
"""
        if p is None:
            p = pyro.config['def_p']

        # Generic iteration parameters
        N = 100 # Maximum iterations
        small = 1e-8    # A "small" number
        epsilon = 1e-6  # Iteration precision

        it = np.nditer((None,value,p),op_flags=[['readwrite','allocate'],['readonly','copy'],['readonly','copy']],op_dtypes='float')
        for T_,y_,p_ in it:
            # Use Tk as the iteration parameter.  We will write to T_ later.
            # Initialize it to be in the center of the species' legal range.
            Tk = 0.5*(self.data['Tlim'][0] + self.data['Tlim'][-1])
            Tk1 = Tk
            # Initialize dT - the change in T
            dT = 0.
            # Calculate an error threshold
            thresh = max(small, abs(epsilon * y_))
            # Initialize absek1 - the absolute error from the last iteration
            #    Using +infty will force the error reduction criterion to be met
            abs_ek1 = float('inf')
            fail = True
            for count in range(N):
                ## CALL THE PROPERTY FUNCTION ##
                yk = prop(T=Tk,p=p_)
                # Calculate error
                ek = yk-y_
                abs_ek = abs(ek)
                # Test for convergence
                if abs_ek < thresh:
                    T_[...] = Tk
                    fail = False
                    break
                # If the error did not reduce from the last iteration
                elif abs_ek > abs_ek1:
                    dT /= 2.
                    Tk = Tk1 + dT
                # Continue normal iteration
                else:
                    # Shift out the old values
                    abs_ek1 = abs_ek
                    Tk1 = Tk
                    ## ESTIMATE THE DERIVATIVE ##
                    dT = max(small, epsilon*Tk)    # Make a small perturbation
                    dydx = (prop(T=Tk+dT,p=p_) - yk)/dT
                    # Calculate the next step size
                    dT = -ek / dydx
                    # Produce a tentative next value
                    Tk = Tk1 + dT
                    # Test the guess for containment in the temperature limits
                    # Shrink the increment until Tk is contained
                    while Tk<self.data['Tlim'][0] or Tk>self.data['Tlim'][-1]:
                        dT /= 2.
                        Tk = Tk1 + dT
            if fail:
                raise pyro.utility.PMAnalysisError('_invT() failed to converge!')
        return it.operands[0]

    def _test(self, report_file=None, report_level=2, basic=False):
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
            raise pyro.utility.PMParamError('Unrecognized file type')

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
        pref = self._pref_bar

        # Create a generic property test function
        # crit  criterion string for the report (1a, 2b, etc...)
        # prop  property method to call (cp, h, s)
        # args  the argv tuple to pass to the property function
        # ref   reference array
        # units = [e, m, p, t] units for energy, matter, pressure, and temperature
        # epsilon  The acceptable fractional error 
        # small    An acceptable total error (for values that can be close to 0)
        def _prop_test(crit, prop, args, ref, units, epsilon=.001, small=0.):
            pyro.config['unit_energy'] = units[0]
            pyro.config['unit_matter'] = units[1]
            pyro.config['unit_pressure'] = units[2]
            pyro.config['unit_temperature'] = units[3]
            
            if report_level >= REP_VITAL:
                ff.write(crit + '...')

            try:
                error = abs(prop(*args) - ref)
            except:
                if report_level>=REP_VITAL:
                    message = os.sys.exc_info()[1].message
                    ff.write('[ERROR]\n')
                    if report_level>=REP_ALL:
                      ff.write('    %s\n'%message)
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

            pyro.config['unit_temperature']='K'
            pyro.config['unit_energy']='J'
            pyro.config['unit_matter']='mol'
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
        p0 = self._pref_bar
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
        pyro.config['unit_energy'] = 'J'
        pyro.config['unit_matter'] = 'mol'
        pyro.config['unit_temperature'] = 'K'
        error = abs(self.R() - pyro.units.const_Ru)
        fail = (error > .0001 * pyro.units.const_Ru)
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
        pyro.config['unit_energy'] = 'J'
        pyro.config['unit_matter'] = 'kg'
        pyro.config['unit_temperature'] = 'K'
        R = pyro.units.const_Ru * 1000. / self.data['mw']
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
        pyro.config['unit_volume'] = 'm3'
        units = ['J','kg','bar','K']

        ##
        ## 5a - mass units
        ##
        p = 5.0
        R = pyro.units.const_Ru * 1000. / self.data['mw']
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
        R = pyro.units.const_Ru
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
        R = pyro.units.const_Ru
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
        R = pyro.units.const_Ru
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


    def contents(self):
        """Construct an atomic contents dictionary
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
        if self._contents is None:
            self._contents = {}
            # Get the chemical formula portion of the ID
            ID = self.data['id'].split('.')[-1]
            for key,value in re.findall('([A-Z][a-z]*)([0-9]*)', ID):
                if value:
                    self._contents[str(key)] = int(value)
                else:
                    self._contents[str(key)] = 1
                
        return self._contents


    def Tlim(self):
        """Temperature limits
    (Tmin, Tmax) = Tlim()
Returns the temperature limits on the ig data set.

Accepts None
Returns unit_temperature
"""
        Tmin = pyro.units.temperature_scale(self.data['Tlim'][0], from_units='K')
        Tmax = pyro.units.temperature_scale(self.data['Tlim'][-1], from_units='K')
        return (Tmin,Tmax)


    def cp(self,T=None,p=None):
        """Constant-pressure specific heat
    cp(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined.  Ideal gas specific 
heat is not actually a function of p, but it is permitted as an argument 
for cross-compatibility between species' function calls.

Accepts unit_temperature
        unit_pressure
Returns unit_energy / unit_matter
"""
        # Check for default values
        if T is None:
            T = pyro.config['def_T']
        # Don't bother checking for the p default value
        # It's value isn't used in the calculation, and None will 
        # still be broadcast correclty as a scalar.

        # Perform temperature conversion
        T = pyro.units.temperature_scale(T,to_units='K')

        # Calculate a scaling factor for the output
        scale = pyro.units.energy(from_units='J')
        scale = pyro.units.matter(scale,self.data['mw'],from_units='mol',exponent=-1)
        scale = pyro.units.temperature(scale,from_units='K',exponent=-1)

        # Create an iterator over T and out
        it = np.nditer((T,None),op_flags=[['readonly','copy'],['readwrite','allocate']],op_dtypes='float')

        for TT,oo in it:
            C = self.data['C'][self._crange(TT)]
            t = TT/1000.
            oo[...] = C[0] + t*(C[1] + t*(C[2] + t*C[3]))
            oo[...] += C[4]/t/t
            oo[...] *= scale
        # Broadcast the result to match the dims of p
        return it.operands[1]
        
    def h(self,T=None,p=None,hf=True):
        """Enthalpy
    h(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined.  Ideal gas enthalpy
is not actually a function of p, but it is permitted as an argument 
for cross-compatibility between species' function calls.

There is an optional parameter, hf, that is a boolean indicating whether
the enthalpy of formation should be included in the enthalpy.  If hf is
false, h() will calculate the enthalpy to be zero a T=298.15K (25C)

Accepts unit_temperature
        unit_pressure
Returns unit_energy / unit_matter
"""
        # Check for default values
        if T is None:
            T = pyro.config['def_T']
        # Don't bother checking for the p default value
        # It's value isn't used in the calculation, and None will 
        # still be broadcast correclty as a scalar.

        # Perform temperature conversion
        T = pyro.units.temperature_scale(T,to_units='K')

        # Calculate a scaling factor for the output
        scale = pyro.units.energy(from_units='kJ')
        scale = pyro.units.matter(scale,self.data['mw'],from_units='mol',exponent=-1)

        # Create an iterator over T and out
        it = np.nditer((T,None),op_flags=[['readonly','copy'],['readwrite','allocate']],op_dtypes='float')
        if hf:
            C7 = 0.
        else:
            C7 = self.data['C'][0][7]
        for TT,oo in it:
            C = self.data['C'][self._crange(TT)]
            t = TT/1000.
            oo[...] = C[5] + t*(C[0] + t*(C[1]/2. + t*(C[2]/3. + t*C[3]/4.)))
            oo[...] -= C[4]/t + C7
            oo[...] *= scale
        # Broadcast the result to match the dims of p
        return it.operands[1]

    def s(self,T=None,p=None):
        """Entropy
    s(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined. 

Accepts unit_temperature
        unit_pressure
Returns unit_energy / unit_matter / unit_temperature
"""
        # Check for default values
        if T is None:
            T = pyro.config['def_T']
        if not isinstance(T,np.ndarray):
            T = np.array(T)

        if p is None:
            p = pyro.config['def_p']
        if not isinstance(p,np.ndarray):
            p = np.array(p)

        # Perform temperature conversion
        T = pyro.units.temperature_scale(T,to_units='K')
        # and pressure conversion
        p = pyro.units.pressure(p,to_units='bar')

        # Calculate a scaling factor for the output
        scale = pyro.units.energy(from_units='J')
        scale = pyro.units.matter(scale,self.data['mw'],from_units='mol',exponent=-1)
        scale = pyro.units.temperature(scale,from_units='K',exponent=-1)

        # Create an iterator over T, p, and out
        it = np.nditer((T,p,None),op_flags=[['readonly','copy'],['readonly','copy'],['readwrite','allocate']],op_dtypes='float')

        for TT,pp,oo in it:
            C = self.data['C'][self._crange(TT)]
            t = TT/1000.
            oo[...] = C[6] + C[0]*np.log(t)
            oo[...] += t*(C[1] + t*(C[2]/2. + t*C[3]/3.))
            oo[...] -= C[4]/t/t/2.
            oo[...] -= pyro.units.const_Ru * np.log(pp/self._pref_bar)
            oo[...] *= scale
        # Broadcast the result to match the dims of p
        return it.operands[2]

    def e(self,T=None,p=None,hf=True):
        """Energy
    e(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined.  Ideal gas enthalpy
is not actually a function of p, but it is permitted as an argument 
for cross-compatibility between species' function calls.

There is an optional parameter, hf, that is a boolean indicating whether
the enthalpy of formation should be included in the enthalpy.  If hf is
False, h() will calculate the enthalpy to be zero a T=273.15K, and the
internal energy will also be adjusted.

Accepts unit_temperature
        unit_pressure
Returns unit_energy / unit_matter
"""
        # Check for default values
        if T is None:
            T = pyro.config['def_T']
        # Don't bother checking for the p default value
        # It's value isn't used in the calculation, and None will 
        # still be broadcast correclty as a scalar.

        # Perform temperature conversion
        T = pyro.units.temperature_scale(T,to_units='K')

        # Calculate a scaling factor for the output
        scale = pyro.units.energy(from_units='kJ')
        scale = pyro.units.matter(scale,self.data['mw'],from_units='mol',exponent=-1)

        # Create an iterator over T and out
        it = np.nditer((T,None),op_flags=[['readonly','copy'],['readwrite','allocate']],op_dtypes='float')

        R = 1e-3 * pyro.units.const_Ru
        if hf:
            C7 = 0.
        else:
            C7 = self.data['C'][0][7]
        for TT,oo in it:
            C = self.data['C'][self._crange(TT)]
            t = TT/1000.
            oo[...] = C[5] + t*(C[0] + t*(C[1]/2. + t*(C[2]/3. + t*C[3]/4.)))
            oo[...] -= C[4]/t + C7
            oo[...] -= TT * R
            oo[...] *= scale
        # Broadcast the result to match the dims of p
        return it.operands[1]


    def d(self,T=None,p=None):
        """Density
    d(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined.  It is important to 
note that density can be either a molar or mass density depending on the
"unit_matter" configuration directive.

Accepts unit_temperature
        unit_pressure
Returns unit_matter / unit_volume
"""
        # Check for default values
        if T is None:
            T = pyro.config['def_T']
        elif hasattr(T,'__iter__') and not isinstance(T,np.ndarray):
            T = np.array(T)
        if p is None:
            p = pyro.config['def_p']
        elif hasattr(p,'__iter__') and not isinstance(p,np.ndarray):
            p = np.array(p)
        
        p = pyro.units.pressure(p, to_units='Pa')
        T = pyro.units.temperature_scale(T, to_units='K')
        R = pyro.units.matter(pyro.units.const_Ru, self.data['mw'], from_units='mol', exponent=-1)

        return pyro.units.volume(p / R / T, from_units='m3', exponent=-1)


    def cv(self,T=None,p=None):
        """Constant-volume specific heat
    cv(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined.  Ideal gas specific 
heat is not actually a function of p, but it is permitted as an argument 
for cross-compatibility between species' function calls.

Accepts unit_temperature
        unit_pressure
Returns unit_energy/unit_matter/unit_temperature
"""
        return self.cp(T,p) - self.R()

    def mw(self,T=None,p=None):
        """Molecular weight
    mw(T,p)
Accepts temperature and pressure to conform with the property method 
prototype, but ignores their values.  This method returns a scalar value
in all cases.

Accepts unit_temperature
        unit_pressure
Returns unit_mass/unit_molar
"""
        mw = pyro.units.mass(self.data['mw'],from_units='g')
        mw = pyro.units.molar(mw,from_units='mol',exponent=-1)
        return mw

    def R(self,T=None,p=None):
        """Ideal gas constant
    R(T,p)
Accepts temperature and pressure to conform with the property method 
prototype, but ignores their values.  This method returns a scalar value
in all cases.

Accepts unit_temperature
        unit_pressure
Returns unit_energy/unit_temperature/unit_matter
"""
        R = pyro.units.energy(pyro.units.const_Ru, from_units='J')
        R = pyro.units.temperature(R, from_units='K', exponent=-1)
        R = pyro.units.matter(R, self.data['mw'], from_units='mol', exponent=-1)
        return R

    def gam(self,T=None,p=None):
        """Specific heat ratio
    gam(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined.  Ideal gas specific 
heat ratio is not actually a function of p, but it is permitted as an 
argument for cross-compatibility between species' function calls.

Accepts unit_temperature
        unit_pressure
Returns dimensionless
"""
        cp = self.cp(T,p)
        return cp/(cp-self.R())

    def p_d(self,d,T=None):
        """Pressure from density and temperature
    p = ig_instance.p_d(d)
        or
    p = ig_instance.p_d(d,T)

Returns the pressure as a function of density and temperature

Accepts unit_matter / unit_volume
        unit_temperature
Returns unit_pressure
"""
        p0 = pyro.config['def_p']
        d0 = self.d(T=T,p=def_p)
        return p0 * d / d0

    def p_s(self,s,T=None):
        """Pressure as a function of entropy
    p = ig_instance.p_s(s)
        or
    p = ig_instance.p_s(s,T)

Returns the pressure as a function of entropy and temperature.

Accepts unit_energy / unit_matter / unit_temperature
        unit_temperature
Returns unit_pressure
"""
        def_p = pyro.config['def_p']
        s0 = self.s(T=T,p=def_p)
        return def_p * np.exp((s0 - s)/self.R())


    def T_s(self,s,p=None):
        """Temperature as a function of entropy
    T = T_s(s)
        or
    T = T_s(s,p)

Accepts unit_energy / unit_matter / unit_temperature
        unit_pressure
Returns unit_temperature
"""
        return self._invT(s, self.s, p)


    def T_h(self,h,p=None):
        """Temperature as a function of enthalpy
    T = T_h(h)
        or
    T = T_h(h,p)

Returns the temperature as a function of enthalpy and pressure

Accepts unit_energy / unit_matter / unit_temperature
        unit_pressure
Returns unit_temperature
"""
        return self._invT(h, self.h, p)

    def T_d(self,d,p=None):
        """Temperature from density and pressure
    T = ig_instance.T_d(d)
        or
    T = ig_instance.T_d(d,p)

Returns the temperature as a function of density and pressure

Accepts unit_matter / unit_volume
        unit_pressure
Returns unit_temperature
"""
        T0 = pyro.config['def_T']
        d0 = self.d()
        # Convert to an absolute scale
        T0 = pyro.units.temperature_scale(T0,to_units='abs')
        return pyro.units.temperature_scale(T0 * d0 / d, from_units='abs')
