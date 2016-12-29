##############################################
##                                          ##
##  Ideal Gas Tabular Data Class            ##
##                                          ##
##############################################

class igtab(__basedata__):
    """Ideal gas tabular data class
"""

    # redefine the mandatory list
    mandatory = [
        'id',       # pyro-mandatory species identifier string
        'doc',      # pyro-mandatory documentation string
        'class',    # pyro-mandatory evaluation tag
        'T',        # temperature array
        'h',        # enthalpy array
        's',        # entropy array
        'Pref',     # reference pressure (usually 1.01325bar)
        'mw',       # molecular weight
        ]
    


    def _lookup(self,T,Y):
        """Perform an array of table lookups
    igt._lookup(T,Y)

For points inside the data table, performs cubic interpolation.  For 
points at the extrema of the table, performs linear interpolation.  
Values outside the table's range will be returned as nan.

T must be a numpy array.  If it is a single-element array, then the 
lookup is performed right away.  If it is a multi-element array, then
the algorithm loops through the individual elements and recurses to 
perform the lookups individually.  

There are alternative approaches where T is sorted and lookups are 
performed in batches based on how they fall between table entries.  
The overhead for this kind of approach is relatively high, making the 
much simpler element-by-element approach preferable for most 
applications.
"""
        # if T is an array, loop through the elements
        # and recurse into each to perform the lookup
        if T.size > 1:
            out = pyro.utility.np.zeros( T.shape )
            for k in range(T.size):
                out[k] = self._lookup(T[k],Y)
            return out
        else:
            T.squeeze()
            X = self.data['T']
            N = len(X)
            # check that the element is in the table's range
            if T<X[0] or T>X[-1]:
                return float('nan')

            # search for the correct element with bisection
            a = 0
            b = int(N-1)
            while b-1>a:
                # v1.3 added int() to force integer in Python3
                c = int((a+b)/2)
                if T < X[c]:
                    b = c
                else:
                    a = c
            # if the element is at the extrema
            if (a==0) or (b==N-1):
                # use linear interpolation/extrapolation
                return Y[a] + (T-X[a])*(Y[b]-Y[a])/(X[b]-X[a])
            else:
                # use cubic interpolation
                X = X[a-1:a+3]
                Y = Y[a-1:a+3]
                TX = T - X
                out = 0.
                for ii in range(4):
                    term = Y[ii]
                    for jj in range(4):
                        if ii!=jj:
                            term *= (T - X[jj])/(X[ii] - X[jj])
                    out += term
                return out



    #
    # IGTAB evaluation functions
    #
    def cp(self,T=None,p=None):
        """Constant-pressure specific heat"""
        (T,p) = self._vectorize(T,p)
        return self._lookup(T,self.data['cp'])



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
        return self._lookup(T,self.data['h'])

        
    def e(self,T=None,p=None):
        """Internal energy"""
        return self.h(T,p) - self.R() * T
        

    def mw(self,T=None,p=None):
        """Molecular weight"""
        return self.data['mw']
        

    def k(self,T=None,p=None):
        """Specific heat ratio"""
        return 1.0 / ( 1.0 - self.R() / self.cp(T,p) )
        
        
    def R(self,T=None,p=None):
        """Ideal gas constant"""
        Ru = 8.3144621
        return Ru / self.data['mw']
        
        
        
    def s(self,T=None,p=None):
        """Entropy"""
        (T,p) = self._vectorize(T,p)
        Pref = self.data['Pref']
        return (self._lookup(T,self.data['s']) 
            - self.R()*pyro.utility.np.log(p/Pref))
