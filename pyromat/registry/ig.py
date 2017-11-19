import pyromat as pyro
import numpy as np



class ig(pyro.reg.__basedata__):
    """Ideal gas class using the Shomate equation of state"""

    def __init__(self,*arg,**kwarg):
        super(self.__class__,self).__init__(*arg,**kwarg)

        # Important constants
        self._pref_bar = 1.0
        self._Ru = 8314.

#        # Define inverstion routines
#        self.T_h = pyro.solve.solve1n('T',
#            f=self.h, df=self.cp,
#            param_lim = (self.data['Tmin'], self.data['Tmax']))
#
#        def ds(T,p=None):
#            return self.cp(T,p)/T 
#
#        self.T_s = pyro.solve.solve1n('T',
#            f=self.s, df=ds,
#            param_lim = (self.data['Tmin'], self.data['Tmax']))
        

    def _crange(self, T):
        """Return the temperature range index and raise a meaningful exception
if it is out of range."""
        Tlim = self.data['Tlim']
        if T<Tlim[0]:
            raise pyro.utility.PMParamError(
                'Temperature, %f, is out of range for %s (%f to %f)'%(
                T, self.data['id'], Tlim[0], Tlim[-1]))
        for index in range(1,len(self.data['Tlim'])):
            if T<Tlim[index]:
                return index-1
        raise pyro.utility.PMParamError(
            'Temperature, %f, is out of range for %s (%f to %f)'%(
            T, self.data['id'], Tlim[0], Tlim[-1]))


    def cp(self,T=None,p=None):
        """Constant-pressure specific heat
    cp(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined.  Ideal gas specific 
heat is not actually a function of p, but it is permitted as an argument 
for cross-compatibility between species' function calls.
"""
        # Check for default values
        if T is None:
            T = pyro.utility.get_config('def_T')
        # Don't bother checking for the p default value
        # It's value isn't used in the calculation, and None will 
        # still be broadcast correclty as a scalar.

        # Perform temperature conversion
        T = pyro.units.temperature_scale(T,to_units='K')

        # Calculate a scaling factor for the output
        scale = pyro.units.energy(from_units='J')
        scale = pyro.units.matter(scale,self.data['mw'],from_units='mol',exponent=-1)
        scale = pyro.units.temperature(scale,from_units='K',exponent=-1)

        # Initialize the result based on T
        # Since p doesn't play a role, the result can just be broadcast
        # to match p's dimensions at the end
        out = np.zeros_like(T)
        # Create an iterator over T and out
        it = np.nditer((T,out),op_flags=[['readonly'],['readwrite']])

        for TT,oo in it:
            C = self.data['C'][self._crange(TT)]
            t = TT/1000.
            oo[...] = C[0] + t*(C[1] + t*(C[2] + t*C[3]))
            oo[...] += C[4]/t/t
            oo[...] *= scale
        # Broadcast the result to match the dims of p
        return np.broadcast_to(out, np.broadcast(T,p).shape)
        
    def h(self,T=None,p=None):
        """Enthalpy
    h(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined.  Ideal gas enthalpy
is not actually a function of p, but it is permitted as an argument 
for cross-compatibility between species' function calls.
"""
        # Check for default values
        if T is None:
            T = pyro.utility.get_config('def_T')
        # Don't bother checking for the p default value
        # It's value isn't used in the calculation, and None will 
        # still be broadcast correclty as a scalar.

        # Perform temperature conversion
        T = pyro.units.temperature_scale(T,to_units='K')

        # Calculate a scaling factor for the output
        scale = pyro.units.energy(from_units='kJ')
        scale = pyro.units.matter(scale,self.data['mw'],from_units='mol',exponent=-1)

        # Initialize the result based on T
        # Since p doesn't play a role, the result can just be broadcast
        # to match p's dimensions at the end
        out = np.zeros_like(T)
        # Create an iterator over T and out
        it = np.nditer((T,out),op_flags=[['readonly'],['readwrite']])

        for TT,oo in it:
            C = self.data['C'][self._crange(TT)]
            t = TT/1000.
            oo[...] = C[5] + t*(C[0] + t*(C[1]/2. + t*(C[2]/3. + t*C[3]/4.)))
            oo[...] -= C[4]/t
            oo[...] *= scale
        # Broadcast the result to match the dims of p
        return np.broadcast_to(out, np.broadcast(T,p).shape)

    def s(self,T=None,p=None):
        """Entropy
    s(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined. 
"""
        # Check for default values
        if T is None:
            T = pyro.utility.get_config('def_T')
        if p is None:
            p = pyro.utility.get_config('def_p')
        # Don't bother checking for the p default value
        # It's value isn't used in the calculation, and None will 
        # still be broadcast correclty as a scalar.

        # Perform temperature conversion
        T = pyro.units.temperature_scale(T,to_units='K')
        # and pressure conversion
        p = pyro.units.pressure(p,to_units='bar')

        # Calculate a scaling factor for the output
        scale = pyro.units.energy(from_units='J')
        scale = pyro.units.matter(scale,self.data['mw'],from_units='mol',exponent=-1)
        scale = pyro.units.temperature(scale,from_units='K',exponent=-1)

        # Initialize the result based on T
        # Since p doesn't play a role, the result can just be broadcast
        # to match p's dimensions at the end
        out = np.zeros(np.broadcast(T,p).shape)
        # Calculate the gas constant
        R = self._Ru / self.data['mw']
        # Create an iterator over T, p, and out
        it = np.nditer((T,p,out),op_flags=[['readonly'],['readonly'],['readwrite']])

        for TT,pp,oo in it:
            C = self.data['C'][self._crange(TT)]
            t = TT/1000.
            oo[...] = C[6] + C[0]*np.log(t)
            oo[...] += t*(C[1] + t*(C[2]/2. + t*C[3]/3.))
            oo[...] -= C[4]/t/t/2.
            oo[...] -= R * np.log(pp/self._pref_bar)
            oo[...] *= scale
        # Broadcast the result to match the dims of p
        return out

    def e(self,T=None,p=None):
        """Internal Energy
    e(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined.  Ideal gas internal
energy is not actually a function of p, but it is permitted as an argument 
for cross-compatibility between species' function calls.
"""
        if T is None:
            T = pyro.utility.get_config('def_T')
        elif hasattr(T,'__iter__') and not isinstance(T,np.ndarray):
            T = np.array(T)
        return self.h(T,p) - T * self._Ru / self.data['mw'] / 1e3

    def d(self,T=None,p=None):
        """Density
    d(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined.
"""
        # Check for default values
        if T is None:
            T = pyro.utility.get_config('def_T')
        elif hasattr(T,'__iter__') and not isinstance(T,np.ndarray):
            T = np.array(T)
        if p is None:
            p = pyro.utility.get_config('def_p')
        elif hasattr(p,'__iter__') and not isinstance(p,np.ndarray):
            p = np.array(p)
        
        return p * 1e5 * self.data['mw'] / self._Ru / T

    def cv(self,T=None,p=None):
        """Constant-volume specific heat
    cv(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined.  Ideal gas specific 
heat is not actually a function of p, but it is permitted as an argument 
for cross-compatibility between species' function calls.
"""
        return self.cp(T,p) - self._Ru / self.data['mw']/1e3

    def mw(self,T=None,p=None):
        """Molecular weight
    mw(T,p)
Accepts temperature and pressure to conform with the property method 
prototype.  Returns a broadcasted array that conforms to the shape of
T and p.  This leverages Numpy's references for efficiency; only one
scalar value is returned, but the wrapper class makes it appear like
an array of any size.
"""
        return np.broadcast_to(self.data['mw'], np.broadcast(T,p).shape)

    def R(self,T=None,p=None):
        """Ideal gas constant
    R(T,p)
Accepts temperature and pressure to conform with the property method 
prototype.  Returns a broadcasted array that conforms to the shape of
T and p.  This leverages Numpy's references for efficiency; only one
scalar value is returned, but the wrapper class makes it appear like
an array of any size.
"""
        return self._Ru / self.mw(T,p)


    def k(self,T=None,p=None):
        """Specific heat ratio
    k(T,p)
Both arguments are optional, and will default to 'def_T' and 'def_p'
configuration parameters if they are left undefined.  Ideal gas specific 
heat ratio is not actually a function of p, but it is permitted as an 
argument for cross-compatibility between species' function calls.
"""
        cp = self.cp(T,p)
        return cp/(cp-self._Ru/self.data['mw']/1e3)
