##############################################
##                                          ##
##  Example PYro data class definition     ##
##                                          ##
##############################################
class classname(__basedata__):
    """PYro data classes are used to store and interpret thermodynamic 
data.  To create a new data class, copy this template and edit 
to suit. 

__basedata__ is intended to be the parent class for all PYro 
data classes.  While it is not intended to be used directly, it 
establishes the basic structure and methods that are common to 
all data the data formats PYro can use.  

When the PYro package is loaded, files in the reg directory are
executed, and objects they create are loaded into the registry.  
The registry is a dictionary of data classes that interpret the data 
loaded from files.  Each data file is expected to contain a 'class'
element, which indicates which class in the registry should be
used to interpret the data.

There are a few requirements on these data classes:
(1)  Firstly, they must have a 'data' attribute, into which the raw 
data file contents will be loaded as a dictionary. 
(2)  They must have a 'mandatory' attribute, which is a list of string 
keys that that a proper 'data' dictionary must contain to be useable by 
that class.  Among other things, the __basetest__() function uses 
this list to check for bad data.
(3)  Data classes should expose the following methods for computing 
basic thermodynamic properties as a function of temperature (in Kelvin) 
and pressure (in bar)
    .cp()   -   constant pressure specific heat (in kJ/kg/K)
    .cv()   -   constant volume specific heat ( same )
    .d()    -   density (in kg/m**3)
    .h()    -   enthalpy (in kJ/kg)
    .e()    -   internal energy (in kJ/kg)
    .mw()   -   molecular weight (in kg/kmol)
    .s()    -   entropy (in kJ/kg/K)

The methods above must be functions of the form X(self,T,p), and must 
accept T and p as numpy arrays (as produced by the utility.vectorize()
function.  The igfit class should be a good example.

These files are executed in a context so that the pyro package root is
exposed to them verbatim.  They can access exception definitions, the 
error/warning messaging functions, and the various modules that are 
loaded in utility (like os,sys,numpy,etc...).  Simply evoke the utility
module in the code (e.g. pyro.utility.os.path.isdir('nope'))
"""

    # the 
    mandatory = [
        'id',       # pyro-mandatory species identifier string
        'doc',      # pyro-mandatory documentation string
        'class',    # pyro-mandatory evaluation tag
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
Basically, this is the bit of code that ensures that each data
element can be evaluated by the class methods.
"""
        pass

    #
    # Class property functions
    #
    def cp(self,T=None,p=None):
        """A function for calculating constant-pressure specific heat."""
        (T,p) = self._vectorize(T,p)
        pass
    def cv(self,T=None,p=None):
        """A function for calculating constant-volume specific heat."""
        pass
    def d(self,T=None,p=None):
        """A function for calculating density."""
        pass
    def h(self,T=None,p=None):
        """A function for calculating enthalpy."""
        pass
    def e(self,T=None,p=None):
        """A function for calculating internal energy."""
        pass
    def mw(self,T=None,p=None):
        """A function for calculating molecular weight."""
        pass
    def s(self,T=None,p=None):
        """A function for calculating entropy."""
        pass
