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
  k()  spec. heat ratio (dless)
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

        # Define inverstion routines
        self.T_h = pyro.solve.solve1n('T',
            f=self.h, df=self.cp, param_init=1000.)

        def ds(T,p=None):
            return self.cp(T,p)/T 

        self.T_s = pyro.solve.solve1n('T',
            f=self.s, df=ds, param_init=1000.)        

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
        for ss,f in W.iteritems():
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
        for ss,f in W.iteritems():
            out += f*pyro.dat.data[ss].cv(T,p)
        return out

    def d(self,T=None,p=None):
        """A function for calculating density."""
        out = 0.
        X = self.X()
        for ss in self.data['contents']:
            out += X[ss]*pyro.dat.data[ss].d(T,p)
        return out

    def h(self,T=None,p=None,hf=True):
        """A function for calculating enthalpy."""
        out = 0.
        # If matter is configured to mass, weight by mass
        if pyro.config['unit_matter'] in pyro.units.mass:
            W = self.Y()
        # Otherwise, weight by mole fraction
        else:
            W = self.X()
        for ss,f in W.iteritems():
            out += f*pyro.dat.data[ss].h(T,p,hf)
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
        for ss,f in W.iteritems():
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
        for ss,f in W.iteritems():
            out += f*pyro.dat.data[ss].s(T,p)
        return out

    def k(self,T=None,p=None):
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
                for key,value in self.data['contents'].iteritems():
                    temp = value / pyro.dat.data[key].data['mw']
                    # Add str() to convert out of unicode
                    self._x[str(key)] = temp
                    N += temp
            # If the species is defined by mole
            else:
                N = 0.  # Mole count
                for key,value in self.data['contents'].iteritems():
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
                for key,value in self.data['contents'].iteritems():
                    # Add str() to convert out of unicode
                    self._y[str(key)] = value
                    M += value
            # If the species is defined by mole
            else:
                M = 0.  # Mole count
                for key,value in self.data['contents'].iteritems():
                    temp = value * pyro.dat.data[key].data['mw']
                    # Add str() to convert out of unicode
                    self._y[str(key)] = temp
                    M += temp
            # Normalize the result
            for key in self._y:
                self._y[key] /= M

        # Return the dictionary
        return self._y


