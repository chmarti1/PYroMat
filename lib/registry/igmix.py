import pyromat as pyro
import numpy as np
######################
##                  ##
##  Mixture class   ##
##                  ##
######################

class igmix(pyro.reg.__basedata__):
    """IGMIX  Ideal gas mixture class

The ideal gas mixture is comprised of components that are ideal gases.  
The properties are calculated by calling the property functions of the
constituents and performing the appropriate weighted averages (by mass
or by volume).  

IGMIX objects offer the following property functions:
  cp() spec. heat       (unit_energy / unit_temperature / unit_matter)
  cv() spec. heat       (unit_energy / unit_temperature / unit_matter)
  d()  density          (unit_matter / unit_volume)
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

The Tlim() method returns the intersection of all the supported 
temperature intervals of the constituents.
  Tlim() temperature limits  (unit_temperature)
"""

    def __init__(self,*arg,**kwarg):
        # Call the basedata class
        super(igmix,self).__init__(*arg,**kwarg)

        # Initialize the static molar and mass fractions
        self._x = None
        self._y = None
      


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

        Tmin,Tmax = self.Tlim()

        it = np.nditer((None,value,p),op_flags=[['readwrite','allocate'],['readonly','copy'],['readonly','copy']],op_dtypes='float')
        for T_,y_,p_ in it:
            # Use Tk as the iteration parameter.  We will write to T_ later.
            # Initialize it to be in the center of the species' legal range.
            Tk = 0.5*(Tmin + Tmax)
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
                    while Tk<Tmin or Tk>Tmax:
                        dT /= 2.
                        Tk = Tk1 + dT
            if fail:
                raise pyro.utility.PMAnalysisError('_invT() failed to converge!')
        return it.operands[0]


    def Tlim(self):
        """Temperature limits
    (Tmin, Tmax) = Tlim()
Returns the temperature limits on the ig data set.

Accepts None
Returns unit_temperature
"""
        Tmin = float('-inf')
        Tmax = float('+inf')
        for ss in self.data['contents']:
            tmin,tmax = pyro.dat.data[ss].Tlim()
            Tmin = max(tmin,Tmin)
            Tmax = min(tmax,Tmax)
        return Tmin,Tmax

    #
    # Class property functions
    #
    def cp(self,T=None,p=None):
        """A function for calculating constant-pressure specific heat."""
        out = 0.
        # If matter is configured to mass, weight by mass
        if pyro.config['unit_matter'] in pyro.units.mass:
            W = self.Y()
        # Otherwise, weight by mole fraction
        else:
            W = self.X()
        for ss,f in W.items():
            out += f*pyro.dat.data[ss].cp(T,p)
        return out

    def cv(self,T=None,p=None):
        """A function for calculating constant-volume specific heat."""
        out = 0.
        # If matter is configured to mass, weight by mass
        if pyro.config['unit_matter'] in pyro.units.mass:
            W = self.Y()
        # Otherwise, weight by mole fraction
        else:
            W = self.X()
        for ss,f in W.items():
            out += f*pyro.dat.data[ss].cv(T,p)
        return out

    def d(self,T=None,p=None):
        """A function for calculating density."""
        out = 0.
        X = self.X()
        for ss in self.data['contents']:
            out += X[ss]*pyro.dat.data[ss].d(T,p)
        return out

    def h(self,T=None,p=None):
        """A function for calculating enthalpy."""
        out = 0.
        # If matter is configured to mass, weight by mass
        if pyro.config['unit_matter'] in pyro.units.mass:
            W = self.Y()
        # Otherwise, weight by mole fraction
        else:
            W = self.X()
        for ss,f in W.items():
            out += f*pyro.dat.data[ss].h(T,p)
        return out

    def e(self,T=None,p=None):
        """A function for calculating internal energy."""
        out = 0.
        # If matter is configured to mass, weight by mass
        if pyro.config['unit_matter'] in pyro.units.mass:
            W = self.Y()
        # Otherwise, weight by mole fraction
        else:
            W = self.X()
        for ss,f in W.items():
            out += f*pyro.dat.data[ss].e(T,p)
        return out

    def mw(self,T=None,p=None):
        """A function for calculating molecular weight."""
        out = 0.
        X = self.X()
        for ss in self.data['contents']:
            out += X[ss]*pyro.dat.data[ss].mw(T,p)
        return out

    def s(self,T=None,p=None):
        """A function for calculating entropy."""
        out = 0.
        # If matter is configured to mass, weight by mass
        if pyro.config['unit_matter'] in pyro.units.mass:
            W = self.Y()
        # Otherwise, weight by mole fraction
        else:
            W = self.X()
        for ss,f in W.items():
            out += f*pyro.dat.data[ss].s(T,p)
        return out

    def gam(self,T=None,p=None):
        """Specific heat ratio"""
        return self.cp(T,p) / self.cv(T,p)



    def X(self):
        """Return a dictionary expressing the mixture composition by volume
    x = mix.X()

Where x is a dictionary with keys corresponding to the species in
the mixture and values indicating their respective mole fraction
in the mixture."""
        # If this is the first call to X
        if self._x is None:
            self._x = {}
            # If the species is defined by mass (and not by mole)
            if self.data['bymass']:
                N = 0.  # Mole count
                for key,value in self.data['contents'].items():
                    temp = value / pyro.dat.data[key].data['mw']
                    # Add str() to convert out of unicode
                    self._x[str(key)] = temp
                    N += temp
            # If the species is defined by mole
            else:
                N = 0.  # Mole count
                for key,value in self.data['contents'].items():
                    # Add str() to convert out of unicode
                    self._x[str(key)] = value
                    N += value
            # Normalize the result
            for key in self._x:
                self._x[key] /= N

        # Return the dictionary
        return self._x



    def Y(self):
        """Return a dictionary expressing the mixture composition by mass
    y = mix.Y()

Where y is a dictionary with keys corresponding to the species in
the mixture and values indicating their respective mass fraction
in the mixture."""
        # If this is the first call to Y
        if self._y is None:
            self._y = {}
            # If the species is defined by mass (and not by mole)
            if self.data['bymass']:
                M = 0.  # Mass count
                for key,value in self.data['contents'].items():
                    # Add str() to convert out of unicode
                    self._y[str(key)] = value
                    M += value
            # If the species is defined by mole
            else:
                M = 0.  # Mole count
                for key,value in self.data['contents'].items():
                    temp = value * pyro.dat.data[key].data['mw']
                    # Add str() to convert out of unicode
                    self._y[str(key)] = temp
                    M += temp
            # Normalize the result
            for key in self._y:
                self._y[key] /= M

        # Return the dictionary
        return self._y


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

