"""PYroMat data class registry

The registry is a dicitonary that resides within the 'reg' module.
It holds the classes responsible for interpreting thermodynamic
data.  When data is loaded, each file is required to contain a 
'class' string, which will be used as a key to look up the class
object intended to be used with the data. 

When the PYroMat package is loaded, all *.py files in the 'reg' 
directory are run, and the definitions in them are incorporated
into the registry dictionary.

The __basedata__ class is the only truely 'built-in' class.  In 
addition to defining the constructor responsible for incorporating
and checking data, it includes default definitions for all mandatory
class members.  User definitions should always point back to the
__basedata__ as a parent.  For ease of editing, '_example.py' shows
an example of a user-defined data class.

Chris Martin (c) 2015,2017
"""

# bring in the root package
import pyromat as pyro


# initialize the registry dicitonary
registry = {}



#####################################
#                                   #
#   Built in prototype data class   #
#                                   #
#####################################
class __basedata__(object):
    """This is the base Pyro data class.
This class is intended to be the basic building block for all Pyro 
data classes.  While it is not intended to be used directly, it 
establishes basic structure and methods that are common to all data 
the data formats Pyro can use.  This class is helpful because of the 
way Pyro is structured.

When a property function (such as pyro.h()) is called, the property 
funciton has no concept of the different data formats available 
(tabular lookup, ideal gas fits, etc...).  Instead, the data file for 
each species contains two essential pieces of information; 'id', which 
is the string by which the species will be identified, and 'class', 
which is a string indicating the dictionary entry in the pyro.registry 
dictionary.  

The registry is a dictionary of data classes (such as this one) that 
should be used to interpret the data loaded from files.  The files 
basically tell the pyro.load() function which of these classes should 
be used to interpret them.  In this way, the system is broadly 
expandable, and easily user-modified.

There are a few requirements on these data classes:
(1)  Firstly, they must have a 'data' attribute, into which the raw 
data file contents will be loaded as a dictionary. 
(2)  They must have a 'mandatory' attribute, which is a list of string 
keys that that a proper 'data' dictionary must contain to be useable by 
that class.  Among other things the basedata._test() function uses 
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
*   other   -   
The functions preceeded by an asterisk are optional.  Either they are 
not meaningful to all data classes (such as R), or, when absent, they 
can be automatically computed from other parameters (such as k = cp/cv).

The methods above must be functions of the form X(self,T,P), and must 
accept T and P as numpy arrays (as produced by the utility.vectorize()
function.  The igfit class should be a good example.

Other functions can be added to the base classes, but unless they are 
implemented in Pyro, they will only be evaluated for pure mixtures.
For example, the classes might expose speed of sound diffusion, or
viscosity functions, but they cannot be combined using typical mass-
based averages like the extensive properties.
"""


    mandatory = ['id','class','doc']


    def __init__(self,data):
        """Generic initializer for the base Pyro data class
  >>> bhd = pyro.utility.basePyrodata( data)
"""
        self.data = data.copy()
        self.__doc__ = self.data['doc']
        self.__basetest__()
        self.__test__()
        
        
        
        
    def __repr__(self):
        return self.data['class'] + ', ' + self.data['id']
        
        
        
        
    def __basetest__(self):
        """Test the data struct for basic Pyro requirements
Raises errors and prints meaningful messages if something important
is missing from the data.  It also checks that the essential methods
are present (see the basedata documentation for more 
information).  This is intended to be a fundamental test applied to 
all Pyro data classes.  For tests unique to each data class, 
define a custom __test__() function.  See the __test__() 
documentation for more details.
"""

        raise_error = False
    
        # Check that the mandatory list contains the essential Pyro entries
        if not 'id' in self.mandatory:
            self.mandatory.append('id')
            pyro.utility.print_warning(
'The ' + repr(self.__class__) + ' mandatory attribute list did not contain the ''id'' tag, so it was added durring test.  This is a requirement of all Pyro objects.  See the __basedata__ class documentation for more information.')
        if not 'class' in self.mandatory:
            self.mandatory.append('class')
            pyro.utility.print_warning(
'The ' + repr(self.__class__) + ' mandatory attribute list did not contain the ''id'' tag, so it was added durring test.  This is a requirement of all Pyro objects.  See the __basedata__ class documentation for more information.')

        # Make sure all the mandatory entries are present.
        out = 'The ' + repr(self.__class__) + ' data does not contain the following mandatory entries: ' 
        contents_error = False
        for mh in self.mandatory:
            if not mh in self.data:
                out += mh + ', '
                contents_error=True
        if contents_error:
            raise_error=True
            out = out[:-2] + '. The data file may be corrupt.'
            pyro.utility.print_error(out)
        
        out = 'The ' + repr(self.__class__) + ' data structure does not contain the mandatory methods for Pyro compatibility.  The class is missing the following methods: '
        if raise_error:
            raise pyro.utility.PyroDataError()



    def _vectorize(self,T,p,out_init=False,allow_scalar=True,def_T=None,def_p=None):
        """'vectorize' T and p inputs
    (T,p) = mydataclass._vectorize(T,p)

The vectorize function accepts numerical or indexable inputs for T and p
and converts them to numpy arrays.  They will be forced to be compatible
for numerical operations.  If T or p are None, they will be replaced by
their respective default values in pyro.config.

_vectorize() is configurable with the a set if optional keyword 
parameters.  They are listed below with their defaults in parentheses.

"out_init"  (False)
If True, _vectorize() will return an additional third parameter
containing a numpy array of zeros appropriately sized to serve as an 
output to a property calculation.

"allow_scalar"  (True)
If true, numerical or scalar values of T and p will be left as 
0-dimensional numpy arrays regardless of the shape of the other 
parameter.  When False, scalars will be repeated in an array matching
the shape of the other parameter.  When "allow_scalar" is False and both
T and p are scalars, they will be forced to be single-element 
one-dimensional arrays.

"def_T" and "def_p"  (None)
These are configurable default values for T and p.  If T or p are passed
as None, then their default values will be used instead.  If no default
is specified, then the 'def_T' and 'def_p' configuration parameters will
be used instead.

**A note about applicaiton**
In applications where T and p need to be used for math (like evaluating
a polynomial), it will likely be unnecessary to initialize an output,
and numpy's array arithmetic will gracefully handle a scalar interacting
with an array.  For those reasons, out_init defaults to False and 
allow_scalar defaults to True.

However, in applications where T and p need to be selectively indexed
and corresponding values assigned to the output, these behaviors aren't
helpful.  In these cases, allow_scalar should be set to False, and 
developers may want to add an output initialization.
"""

        # assign default values if none are provided
        if T is None:
            if def_T is None:
                T = pyro.utility.get_config('def_T')
            else:
                T = def_T
        if p is None:
            if def_p is None:
                p = pyro.utility.get_config('def_p')
            else:
                p = def_p

        # force T to be an array, but do not make a copy if it already is
        if not isinstance(T, pyro.utility.np.ndarray):
            T = pyro.utility.np.array(T)
        # force P to be an array
        if not isinstance(p, pyro.utility.np.ndarray):
            p = pyro.utility.np.array(p)

        # if scalars are disallowed
        if not allow_scalar:
            # force scalars into arrays
            if T.size==1 and p.size>1:
                temp = pyro.utility.np.zeros(p.shape)
                temp.fill(T)
                T = temp
            elif T.ndim==0:
                T = T.reshape((1,))
            
            if p.size==1 and T.size>1:
                temp = pyro.utility.np.zeros(T.shape)
                temp.fill(p)
                p = temp
            elif p.ndim==0:
                p = p.reshape((1,))

        # raise an error if the arrays are incompatible
        if p.size>1 and T.size>1 and p.shape!=T.shape:
            raise pyro.utility.PyroInputError(
            'T and p vectors are incompatible shapes')

        if out_init:
            if T.size>=p.size:
                N = T.shape
            else:
                N = p.shape
            return (T,p,pyro.utility.np.zeros(N))
        return (T,p)




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



    def psolve(self, Tinit=None, pinit=None, **kwargs):
        """Solve for pressure and temperature given property values
    (T,p) = subst.psolve( prop1=value1, prop2=value2 )
    
The properties specified in the keyword arguments to psolve() must 
correspond to property methods the class posesses (such as cp, h, or s) 
or to T or p themselves.  If only one property is encountered, then
pressure is implicitly assumed to be the value specified by the 'def_P'
configuration parameter (1.01325 by default).

The behavior of psolve() is modifiable through some optional keywords:
 keyword    Description (default value)
 Tinit      Initial temperature guess   (None, uses def_T)
 Tprec      Temperature fractional precision (0.1K)
 pinit      Initial pressure guess (None, uses def_P)
 pprec      Pressure precision (.001bar)

*** For Developers ***
To do its job, the __basedata__ implementation of psolve uses a modified 
two-dimensional newton-raphson method.  That works fine in most data 
classes, but it will struggle if there are discontinuities (like with 
phase changes.)  If you're developing for multi-phase data, you might 
need to write your own implementation of psolve.
"""
        if len(kwargs)==1:
            kwargs['p'] = pyro.utility.get_config('def_p')
        elif len(kwargs)==0:
            raise pyro.utility.PyroInputError('Constraints are required.')
        elif len(kwargs)>2:
            raise pyro.utility.PyroInputError(
                'Only two constraints may be specified.')
        
        # constrained properties
        constr = set(kwargs.keys())
        fun = []
        
        # construct a list with the property functions
        # and their respective target values
        if 'T' in constr:
            fun.append((
                (lambda T,p: T),    # a function that returns the temperature
                kwargs['T'],         # the temperature constraint
                0.1
                ))
            constr.remove('T')
        if 'p' in constr:
            fun.append((
                (lambda T,p: p),    # a function that returns the pressure
                kwargs['p'],         # the pressure constraint
                .0001
                ))
            constr.remove('p')
        # If there are no constraints left, then T and p were constrained
        # Just return the values directly
        if len(constr)==0:
            return (fun[0][1],fun[1][1])
        # Now add the constraints
        for this in constr:
            fun.append((
                self.__getattribute__(this),
                kwargs[this],
                .0001
                ))

        Tinit = pyro.utility.get_config('def_T')
        pinit = pyro.utility.get_config('def_p')
        return pyro.utility.newtoniter(
            f=fun[0][0], fval=fun[0][1], fprec=fun[0][2],
            g=fun[1][0], gval=fun[1][1], gprec=fun[0][2],
            T=Tinit, p=pinit )
#
#   Go load the contents of the reg directory
#


def regload(verbose = None):
    """regload - reloads the data class registry

regload() is automatically executed when Pyro is first imported,
but it can be re-run at any time to update the contents of the reg
directory, or to bring the registry up to date with changes to the
registry search path.

pyro.config parameters that affect the behavior of regload() are
'reg_dir'
    directories in which to search
'reg_verbose'
    operate verbosely?
'reg_overwrite'
    overwrite existing classes with redundant ones?
'reg_exist_fatal'
    exit with an error when a redundant class is discovered?
"""

    # initialize the registry
    global registry 
    registry = {}
    lead = 'regload->'

    # fetch the configuration parameters
    if verbose == None:
        verbose = pyro.utility.get_config('reg_verbose',dtype=bool)
    exist_fatal = pyro.utility.get_config('reg_exist_fatal',dtype=bool)
    exist_overwrite = pyro.utility.get_config('reg_overwrite',dtype=bool)

    # search each directory in the registry search path
    for loc in pyro.config['reg_dir']:
        cont = pyro.utility.os.listdir(loc)
        # modules to load should not begin with an underscore or period
        # modules to load should end with .py
        # modules should contain a single class matching the name of the file
        for fil in cont:
            f_go = fil[0]!='_' and fil[0]!='.'
            f_go = f_go & (len(fil)>3) & (fil[-3:]=='.py')
            # if the filename qualifies.
            if f_go:
                temp = {}
                thisfile = pyro.utility.os.path.join(loc,fil)
                if verbose:
                    pyro.utility.print_line('Examining file "' + thisfile + '"', lead)
                try:
                    with open(thisfile,'r') as ff:
                        exec(compile(ff.read(),thisfile,'exec'),globals(),temp)
                    # In Python 2.7, this used to work
                    # Use exec() for 3.4 compatibility
                    # execfile(thisfile,globals(),temp)
                except:
                    pyro.utility.print_warning(
'Failed to execute file: ' + pyro.utility.os.path.join(loc,fil) +
'.  Encountered exception: ' + repr(pyro.utility.sys.exc_info()[1]))

                # loop through all variables created in the file
                valid = False
                for new in temp:
                    if isinstance(temp[new],type) and issubclass(temp[new],__basedata__):
                        valid = True
                        # if the class is already registered, either raise 
                        # an exception, or throw a warning
                        if new in registry:
                            fullfil = pyro.utility.os.path.abspath( thisfile )
                            if exist_fatal:
                                pyro.utility.print_error(
'Encountered a redundant definition for data class "' + new + '" in file "' + 
fullfil + '"')
                                raise pyro.utility.PyroFileError()
                            elif exist_overwrite:
                                pyro.utility.print_warning(
'Overwriting a redundant definition for data class "' + new + 
'" with the definition in file "' + fullfil + '"')
                                registry[new] = temp[new]
                            else:
                                pyro.utility.print_warning(
'Ignoring a redundant definition for data class "' + new + '" in file "' + 
fullfil + '"')

                        # if everything is fine, add the class to the registry
                        else:
                            registry[new] = temp[new]
                            if verbose:
                                pyro.utility.print_line(
'Found class "' + new + '"', lead)
                if not valid:
                    pyro.utility.print_warning(
'File "' + pyro.utility.os.path.abspath(thisfile) + 
'" was found in a registry directory, but contained no data class definition.')



regload()
