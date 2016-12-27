##############################################
##                                          ##
##  Ideal Gas Curve Fit Evaluation Class    ##
##                                          ##
##############################################
class igfit(__basedata__):
    """Ideal gas specific heat curve fit data class
"""

    # redefine the mandatory list
    mandatory = [
        'id',       # pyro-mandatory species identifier string
        'doc',      # pyro-mandatory documentation string
        'class',    # pyro-mandatory evaluation tag
        'C1',       # specific-heat coefficients (low-temperature)
        'C2',       # specific-heat coefficients (high-temperature)
        'h1',       # enthalpy integration coefficient (low-temperature)
        'h2',       # enthalpy integration coefficient (high-temperature)
        's1',       # entropy integration coefficient (low-temperature)
        's2',       # entropy integration coefficient (high-temperature)
        'Tmax',     # maximum valid temperature
        'Tmin',     # minimum valid temperature
        'Tsplit',   # threshold temperature between high and low coefficients
        'Pref',     # reference pressure (usually 1.01325bar)
        'mw',       # molecular weight
        ]
    
    
    #
    # IGFIT evaluation functions
    #
    def cp(self,T=None,p=None):
        """Constant-pressure specific heat"""

        (T,p) = self._vectorize(T,p)

        Tsplit = self.data['Tsplit']
        Tmax = self.data['Tmax']
        Tmin = self.data['Tmin']
        C1 = self.data['C1']
        C2 = self.data['C2']

        N = T.size

        if N==1:
            out = pyro.utility.np.zeros(())

            if T<=Tsplit:
                for k in range(4,-1,-1):
                    out = out*T + C1[k]
            else:
                for k in range(4,-1,-1):
                    out = out*T + C2[k]

        else:
            # Develop boolean arrays indicating which 
            I1 = T<=Tsplit
            I2 = T>Tsplit
            Inan = (T>Tmax) + (T<Tmin)
            
            # number of elements
            N = T.size
            # initialize the output
            out = pyro.utility.np.zeros(T.shape)

            T1 = T[I1]
            T2 = T[I2]

            for c in C1[-1::-1]:
                out[I1] = out[I1]*T1 + c
            for c in C2[-1::-1]:
                out[I2] = out[I2]*T2 + c

            out[Inan] = float('nan')

        return out



    def cv(self,T=None,p=None):
        """Constant-volume specific heat"""
        return self.cp(T,p) - self.R()



    def d(self,T=None,p=None):
        """Density."""
        (T,p) = self._vectorize(T,p)
        # Convert R into J from kJ
        # Convert p into Pa from bar
        # net result is 100
        return p*100/(self.R() * T)
        
        
    def h(self,T=None,p=None):
        """Enthalpy"""

        (T,p) = self._vectorize(T,p)

        Tsplit = self.data['Tsplit']
        Tmax = self.data['Tmax']
        Tmin = self.data['Tmin']
        C1 = self.data['C1']
        C2 = self.data['C2']
        h1 = self.data['h1']
        h2 = self.data['h2']

        N = T.size

        if N==1:
            out = pyro.utility.np.zeros(())

            if T<Tmin or T>Tmax:
                out = float('nan')
            elif T<=Tsplit:
                for k in range(4,-1,-1):
                    out = out*T + C1[k]/(k+1)
                out *= T
                out += h1
            else:
                for k in range(4,-1,-1):
                    out = out*T + C2[k]/(k+1)
                out *= T
                out += h2
        else:
            # Develop boolean arrays indicating which 
            I1 = T<=Tsplit
            I2 = T>Tsplit
            Inan = (T>Tmax) + (T<Tmin)
            
            # number of elements
            N = T.size
            # initialize the output
            out = pyro.utility.np.zeros(T.shape)

            T1 = T[I1]
            T2 = T[I2]
            for k in range(4,-1,-1):
                out[I1] = out[I1]*T1 + C1[k]/(k+1)
            for k in range(4,-1,-1):
                out[I2] = out[I2]*T2 + C2[k]/(k+1)
            out *= T
            out[I1] += h1
            out[I2] += h2
            out[Inan] = float('nan')

        return out

        

    def e(self,T=None,p=None):
        """Internal energy"""
        (T,p) = self._vectorize(T,p)
        return self.h(T,p) - self.R() * T

        
    def k(self,T=None,p=None):
        """Specific heat ratio"""
        (T,p) = self._vectorize(T,p)
        return 1.0 / ( 1.0 - self.R() / self.cp(T,p) )
        
        
    def mw(self,T=None,p=None):
        """Molecular weight"""
        # always return a scalar.
        return self.data['mw']
        
        
        
    def R(self,T=None,p=None):
        """Ideal gas constant"""
        Ru = 8.3144621
        return Ru / self.data['mw']
        
        
        
    def s(self,T=None,p=None):
        """Entropy"""

        (T,p) = self._vectorize(T,p)

        Tsplit = self.data['Tsplit']
        Tmax = self.data['Tmax']
        Tmin = self.data['Tmin']
        C1 = self.data['C1']
        C2 = self.data['C2']
        h1 = self.data['h1']
        h2 = self.data['h2']
        s1 = self.data['s1']
        s2 = self.data['s2']
        Pref = self.data['Pref']

        np = pyro.utility.np

        N = T.size

        if N==1:
            out = np.zeros(())

            if T<Tmin or T>Tmax:
                out = float('nan')
            elif T<=Tsplit:
                for k in range(4,0,-1):
                    out = out*T + C1[k]/k
                out *= T
                out += C1[0] * np.log(T)
                # add integration constants
                out += s1
            else:
                for k in range(4,0,-1):
                    out = out*T + C2[k]/k
                out *= T
                out += C2[0] * np.log(T)
                # add integration constants
                out += s2
        else:
            # Develop boolean arrays indicating which 
            I1 = T<=Tsplit
            I2 = T>Tsplit
            Inan = (T>Tmax) + (T<Tmin)
            
            # number of elements
            N = T.size
            # initialize the output
            out = np.zeros(T.shape)

            T1 = T[I1]
            T2 = T[I2]
            for k in range(4,0,-1):
                out[I1] = out[I1]*T1 + C1[k]/k
            for k in range(4,0,-1):
                out[I2] = out[I2]*T2 + C2[k]/k
            out *= T
            out[I1] += C1[0] * np.log(T1)
            out[I2] += C2[0] * np.log(T2)
            # add integration constants
            out[I1] += s1
            out[I2] += s2
            out[Inan] = float('nan')

        # pressure term
        out -= self.R() * np.log( p/Pref )

        return out

