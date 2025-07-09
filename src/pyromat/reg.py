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

Chris Martin (c) 2015,2017,2021,2022
"""

# bring in the root package
import pyromat as pm


# initialize the registry dicitonary
registry = {}



#####################################
#                                   #
#   Built in prototype data class   #
#                                   #
#####################################
class __basedata__(object):
    """This is the base PYroMat data class.
This class is intended to be the basic building block for all PM 
data classes.  While it is not intended to be used directly, it 
establishes basic structure and methods that are common to all data 
the data formats PYroMat can use.  

When a property function (such as spec.h()) is called, the property 
funciton has no concept of the different data formats available 
(tabular lookup, ideal gas fits, etc...).  Instead, the data file for 
each species contains two essential pieces of information; 'id', which 
is the string by which the species will be identified, and 'class', 
which is a string indicating the dictionary entry in the pm.registry 
dictionary.  

The registry is a dictionary of data classes (such as this one) that 
should be used to interpret the data loaded from files.  The files 
basically tell the load() function which of these classes should 
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
        """Generic initializer for the base PYroMat data class
  >>> bhd = __basedata__(data)
"""
        #self.data = data.copy()
        self.data = data
        self.__basetest__()
        self.__doc__ = self.data['doc']
        
        
    def __repr__(self):
        return '<' + self.data['class'] + ', ' + self.data['id'] + '>'

        
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
        mandatory = ['id','class','doc','fromfile']

        # Make sure all the mandatory entries are present.
        missing = ''
        for mh in mandatory:
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
            pm.utility.print_error(message)
            pm.utility.print_error(missing)
            raise pm.utility.PMDataError()

    def sid(self):
        """Returns the substance identifier string
    sid = subst.sid()

The substance identifier string uniquely identifies this substance in 
the PYroMat database.  It is a string in two parts identifying the 
collection and the substance separated by a period.  For example,
    ig.H2O      Ideal gas water
    mp.H2O      Multi-phase model water
    ig.CN2      CNN radical (CNN)
    ig.CN2_1    methanetetraylbis-Amidogen (NCN)
Note that substances with the same composition but dissimilar atomic 
arrangements are rare in PYroMat, but they do occur.  They are 
differentiated by a trailing underscore and index.
"""
        return self.data['id']

    def pmclass(self):
        """Return the PYroMat class string
    pmclassstr = subst.pmclass()
   
The pmclass string is used to identify which of the PYroMat data classes
to use for the substance data model.  This is a mandatory field, so this
is just a wrapper for subst.data['class'].
"""
        return self.data['class']
        

    def collection(self):
        """Return the name of the collection to which this substance belongs
    collectionstr = subst.collection()
   
The collection is just the portion of the substance identifier string
preceeding the '.' character.  If a substance is encountered without a
'.' character, then collection() will return an empty string.
"""
        parts = self.data['id'].split('.',1)
        if len(parts)==1:
            return ''
        return parts[0]
        
    def names(self):
        """Return a list of the substance names.
    namelist = subst.names()

If the "name" key is found in the dataset, a copy will be returned 
explicitly.  If not, an empty list will be retruned.
"""
        if 'names' in self.data:
            return self.data['names'].copy()
        return []
        
    def inchi(self):
        """Return the IUPAC InChI identifier string
    inchistr = subst.inchi()

If the "inchi" key is found in the dataset, it will be returned 
explicitly.  If not, an empty string will be returned.
"""
        if 'inchi' in self.data:
            return self.data['inchi']
        return ''
        
    def casid(self):
        """Return the CAS registry ID string
    casidstr = subst.casid()

If the "cas" key is found in the dataset, it will be returned 
explicitly.  If not, an empty string will be returned.
"""
        if 'cas' in self.data:
            return self.data['cas']
        return ''

    def atoms(self):
        """Return a dictionary specifying the chemical composition of the substance.
    aa = atoms()
    
The dictionary keys are the symbols of atoms and their corresponding values 
are the integer quantities in the chemical formula.  For example
    aa = {'C':1, 'O':2}
would represent carbon dioxide.

Ionization is represented by an entry for the free electron, 'e'.  For
example, Ne+ (ionized neon) has an atoms dictionary
    aa = {'Ne':1, 'e':-1}
indicating that the neutral neon atom has lost a single electron.

If the substance data does not include atomic data, an empty dictionary
is returned instead.
"""
        if 'atoms' in self.data:
            return self.data['atoms'].copy()
        return {}
        
    def hill(self):
        """Return a string with the Hill notation of the chemical compound
    hillstr = hill()
    
In Hill notation, chemical contents are listed in order [C][H][Others]
where "others" are listed in alphabetical order.
    
If the atoms() method returns a valid composition dictionary, it is used
to build the string.  Otherwise, hill() uses the substance ID string.
"""
        out = ''
        # If the atoms dictionary is available, use that.  This is 
        # preferred over atoms() because igmix instances will have 
        # confusing decimal quantities
        if 'atoms' in self.data:
            aa = self.data['atoms']
            contents = list(aa.keys())
            # The free electron is a special case
            if contents == ['e']:
                return 'e-'
            # Deal with carbon and hydrogen explicitly
            for this in ['C', 'H']:
                if this in contents:
                    out += this
                    if aa[this] != 1:
                        out += f"{aa[this]}"
                    contents.remove(this)
            # If there is charge, stash it for later
            charge = 0
            if 'e' in contents:
                charge = -aa['e']
                contents.remove('e')
            # all others in alphabetical order
            contents.sort()
            for this in contents:
                out += this
                if aa[this] != 1:
                    out += f'{aa[this]}'
            # Finally, add on any ionization
            if charge < 0:
                out += '-'
                charge = -charge
            elif charge > 0:
                out += '+'
            if charge > 1:
                out += f'{qty}'
        else:
            out = self.data['id'].split('.')[1]
            out = out.split('_')[0]
        return out
        
#
#   Go load the contents of the reg directory
#
def regload(verbose = None):
    """regload - reloads the data class registry

regload() is automatically executed when PYroMat is first imported,
but it can be re-run at any time to update the contents of the reg
directory, or to bring the registry up to date with changes to the
registry search path.

pm.config parameters that affect the behavior of regload() are
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
        verbose = pm.config['reg_verbose']
    exist_fatal = pm.config['reg_exist_fatal']
    exist_overwrite = pm.config['reg_overwrite']

    # search each directory in the registry search path
    for loc in pm.config['reg_dir']:
        # Expand references to the users' home directories
        # and environment variables
        loc = pm.utility.os.path.expanduser(loc)
        loc = pm.utility.os.path.expandvars(loc)
        loc = pm.utility.os.path.abspath(loc)

        cont = pm.utility.os.listdir(loc)
        # modules to load should not begin with an underscore or period
        # modules to load should end with .py
        # modules should contain a single class matching the name of the file
        for fil in cont:
            f_go = fil[0]!='_' and fil[0]!='.'
            f_go = f_go & (len(fil)>3) & (fil[-3:]=='.py')
            # if the filename qualifies.
            if f_go:
                temp = {}
                thisfile = pm.utility.os.path.join(loc,fil)
                if verbose:
                    pm.utility.print_line('Examining file "' + thisfile + '"', lead)
                try:
                    with open(thisfile,'r') as ff:
                        exec(compile(ff.read(),thisfile,'exec'),temp)
                    # In Python 2.7, this used to work
                    # Use exec() for 3.4 compatibility
                    # execfile(thisfile,globals(),temp)
                except:
                    pm.utility.print_warning(
'Failed to execute file: ' + pm.utility.os.path.join(loc,fil) +
'.  Encountered exception: ' + repr(pm.utility.sys.exc_info()[1]))

                # loop through all variables created in the file
                valid = False
                for new in temp:
                    if isinstance(temp[new],type) and issubclass(temp[new],__basedata__):
                        valid = True
                        # if the class is already registered, either raise 
                        # an exception, or throw a warning
                        if new in registry:
                            fullfile = pm.utility.os.path.abspath( thisfile )
                            if exist_fatal:
                                pm.utility.print_error(
'Encountered a redundant definition for data class "' + new + '" in file "' + 
fullfile + '"')
                                raise pm.utility.PMFileError()
                            elif exist_overwrite:
                                pm.utility.print_warning(
'Overwriting a redundant definition for data class "' + new + 
'" with the definition in file "' + fullfile + '"')
                                registry[new] = temp[new]
                            else:
                                pm.utility.print_warning(
'Ignoring a redundant definition for data class "' + new + '" in file "' + 
fullfile + '"')

                        # if everything is fine, add the class to the registry
                        else:
                            registry[new] = temp[new]
                            if verbose:
                                pm.utility.print_line(
'Found class "' + new + '"', lead)
                if not valid:
                    pm.utility.print_warning(
'File "' + pm.utility.os.path.abspath(thisfile) + 
'" was found in a registry directory, but contained no data class definition.')




