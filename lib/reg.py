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

When a property function (such as spec.h()) is called, the property 
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
(1)  They must be children of the __basedata__ class
(2)  Their initializers must accept as their sole argument, the 
    dictionary result of the json.load() operation loading the species'
    data file.
(3)  The data dictionary must contain an "id" element, which is the 
    string species ID.
(4)  The data dictionary must contain a "class" element, which is the
    string name of the class from the PYroMat registry to use.
(5)  The data dictionary must contain a "doc" element, which describes
    the origins of the data.

While PYroMat does not explicitly impose rules on species ID or the 
call signatures of the property methods, there is a convention that 
all species IDs will be of the form "collection.formula"  Collection 
is a group of similar species who share basic assumptions or source 
data (like "ig" for ideal gas).  The formula is usually (but not always)
the chemical makeup of the species.  They element names should be in 
one- or two-character groups with the first character always upper-case
and the second (if present) always lower-case.  If the element is 
followed by an integer, that integer indicates the quantity.  If the
integer is omitted, 1 is implied.  For example, "CO2" indicates carbon
dioxide and "Co2" represents diatomic cobalt.
"""


    mandatory = []


    def __init__(self,data):
        """Generic initializer for the base Pyro data class
  >>> bhd = pyro.utility.basePyrodata( data)
"""
        #self.data = data.copy()
        self.data = data
        self.__basetest__()
        self.__doc__ = self.data['doc']
        
        
        
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
        pyro_mandatory = ['id','class','doc','fromfile']

        # Make sure all the mandatory entries are present.
        missing = ''
        for mh in pyro_mandatory:
            if not mh in self.data:
                missing += mh + ' '

        for mh in self.mandatory:
            if not mh in self.data:
                missing += mh + ' '
        
        if missing:
            message = \
'Mandatory entries are missing from a ' + repr(self.__class__) + ' file: '
            if 'fromfile' in self.data:
                message += self.data['fromfile']
            pyro.utility.print_error(message)
            pyro.utility.print_error(missing)
            raise pyro.utility.PyroDataError()


    def _vectorize(self,T,p,out_init=False,allow_scalar=True,def_T=None,def_p=None):
        """'vectorize' T and p inputs
** OBSOLETE **

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




    def psolve(self, Tinit=None, pinit=None, **kwargs):
        """*OBSOLETE* simply returns an error.
PSOLVE() was obsoleted in version 1.4 to improve numerical performance
and stability.  Use SOLVE1(), SOLVE2(), or a class specific solver.  If
it is vital to have psolve() support, install version 1.3.
"""
        raise pyro.utility.PMAnalysisError(
            "PSOLVE() was obsoleted in version 1.4\n" +
            "Use solve1() solve2() or a class-specific solver\n")





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
        verbose = pyro.config['reg_verbose']
    exist_fatal = pyro.config['reg_exist_fatal']
    exist_overwrite = pyro.config['reg_overwrite']

    # search each directory in the registry search path
    for loc in pyro.config['reg_dir']:
        # Expand references to the users' home directories
        # and environment variables
        loc = pyro.utility.os.path.expanduser(loc)
        loc = pyro.utility.os.path.expandvars(loc)
        loc = pyro.utility.os.path.abspath(loc)

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
                        exec(compile(ff.read(),thisfile,'exec'),temp)
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
                            fullfile = pyro.utility.os.path.abspath( thisfile )
                            if exist_fatal:
                                pyro.utility.print_error(
'Encountered a redundant definition for data class "' + new + '" in file "' + 
fullfile + '"')
                                raise pyro.utility.PyroFileError()
                            elif exist_overwrite:
                                pyro.utility.print_warning(
'Overwriting a redundant definition for data class "' + new + 
'" with the definition in file "' + fullfile + '"')
                                registry[new] = temp[new]
                            else:
                                pyro.utility.print_warning(
'Ignoring a redundant definition for data class "' + new + '" in file "' + 
fullfile + '"')

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




