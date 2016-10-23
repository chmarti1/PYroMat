######################
##                  ##
##  Mixture class   ##
##                  ##
######################

class mixture(__basedata__):
    
    mandatory = [
        'id',       # pyro-mandatory species identifier string
        'doc',      # pyro-mandatory documentation string
        'class',    # pyro-mandatory evaluation tag
        'contents', # a dictionary indicating the makeup of the mixture
        'bymass',   # a boolean indicating whether by mass or volume
        ]




    #
    # Class-specific tests at init
    #
    def __test__(self):
        """Perform class-specific data checks
The last step of the __init__() function for the base class
is to execute the __test__ function, so every class needs one.
By default, it is empty, but it is intended to be redefined in
each class to make specific checks on data types and formats.
"""
        
        # First, make sure that 'contents' is a dictionary
        cont = self.data['contents']
        if not isinstance(cont,dict):
            pyro.utility.print_warning('The data for the mixture, "' + self.data['id'] + '", is corrupt or incomplete.  The "contents" element should be a dictionary whose keys name the species present, and whose values are the corresponding amounts.')
        else:
            for ss in cont:
                try:
                    cont[ss] = float(cont[ss])
                except:
                    pyro.utility.print_warning('The quantity of "' + ss + '" in the mixture "' + self.data['id'] + '" is invalid.')


        




    #
    # Class property functions
    #
    def cp(self,T=None,p=None):
        """A function for calculating constant-pressure specific heat."""
        out = 0.
        Y = self.Y()
        for ss in self.data['contents']:
            out += Y[ss]*pyro.dat.data[ss].cp(T,p)
        return out

    def cv(self,T=None,p=None):
        """A function for calculating constant-volume specific heat."""
        out = 0.
        Y = self.Y()
        for ss in self.data['contents']:
            out += Y[ss]*pyro.dat.data[ss].cv(T,p)
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
        Y = self.Y()
        for ss in self.data['contents']:
            out += Y[ss]*pyro.dat.data[ss].h(T,p)
        return out

    def e(self,T=None,p=None):
        """A function for calculating internal energy."""
        out = 0.
        Y = self.Y()
        for ss in self.data['contents']:
            out += Y[ss]*pyro.dat.data[ss].e(T,p)
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
        Y = self.Y()
        for ss in self.data['contents']:
            out += Y[ss]*pyro.dat.data[ss].s(T,p)
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
        contents = self.data['contents']
        N = 0   # temporary var for the total moles
        W = 0   # temporary var for each molecular weight

        if self.data['bymass']:
            # if the data in contents is mass-based
            X = {}
            for ss in contents:
                W = pyro.dat.data[ss].mw()
                X[ss] = contents[ss]/W
                N += X[ss]
            for ss in contents:
                X[ss] /= N
        else:
            # if the data in contents is mole based
            X = contents.copy()
            for ss in contents:
                N += contents[ss]
            for ss in contents:
                X[ss] /= N
        return X





    def Y(self):
        """Return a dictionary expressing the mixture composition by mass
    y = mix.Y()

Where y is a dictionary with keys corresponding to the species in
the mixture and values indicating their respective mass fraction
in the mixture."""
        contents = self.data['contents']
        M = 0   # temporary var for the total mass
        W = 0   # temporary var for each molecular weight

        if self.data['bymass']:
            # if the data in contents is mass-based
            Y = contents.copy()
            for ss in contents:
                M += contents[ss]
            for ss in contents:
                Y[ss] /= M
        else:
            # if the data in contents is mole based
            Y = {}
            for ss in contents:
                W = pyro.dat.data[ss].mw()
                Y[ss] = contents[ss]*W
                M += Y[ss]
            for ss in contents:
                Y[ss] /= M
        return Y





    def N(self):
        """Return a dictionary expressing the mixture composition by mole count
    n = mix.N()

Where n is a dictionary with keys corresponding to the species in
the mixture and values indicating their respective total mole counts 
in the mixture."""
        N = self.data['contents'].copy()
        if self.data['bymass']:
            for ss in N:
                N[ss] /= pyro.dat.data[ss].mw()
        return N





    def M(self):
        """Return a dictionary expressing the mixture composition by mass
    m = mix.M()

Where m is a dictionary with keys corresponding to the species in
the mixture and values indicating their respective total mass in the
mixture."""
        M = self.data['contents'].copy()
        if not self.data['bymass']:
            for ss in M:
                M[ss] *= pyro.dat.data[ss].mw()
        return M


