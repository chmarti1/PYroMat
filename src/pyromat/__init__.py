"""PYroMat    Thermodynamic property calculator for Python

PYroMat is an open-source Python-based software platform for retrieving
the physical properties of substances.  For complete documentation, 
visit "pyromat.org"

Chris Martin (c) 2015-2021
Released under the GNU General Publice License v3.0
  http://www.gnu.org/licenses/gpl-3.0.en.html

*** Retrieving data ***

To get started, retrieve a substance that is of interest
and commit it to a variable that you will use later. The
command below retrieves an object for Argon.
  >>> import pyromat as pm
  >>> Ar = pm.get('ig.Ar')

Once created, these objects can be called on to recover
thermodynamic properties given a temperature in K and
a pressure in bar.  The command below returns the 
enthalpy of argon at 325K and atmospheric pressure.
  >>> h = Ar.h(T=325, p=1.013)

To see what properties are supported, check the in-line
help for each species.
  >>> help(Ar)

*** What's available? ***

For a complete list of all available species, use the info()
function, 
  >>> pm.info()

To search for species, see the info() funciton's other features
  >>> help(pm.info)
"""

# This is the athoritative version number.
# utility.load_config() checks this value to establish the read-only version
# setup.py looks for this line to establish the version
# MUST be unindented
__version__ = "2.1.2"


# loading the PYroMat utility functions
from . import utility
# load the configuration
config = utility.PMConfig()

# import the dataclass registry
from . import reg
# import the module for handling data
from . import dat
# import the units module
from . import units
# By default, do not import the module for handling special applications
# This module has requirements beyond the base pyromat installation
#from . import aps

reg.regload()
dat.load()








def get( name ):
    """Return an object for a substance
    get( 'name' )

Returns a substance data class for the substance named.
"""
    if name in dat.data:
        return dat.data[name]
    
    utility.print_error('No substance named "' + str(name) + '" was found in the loaded data.')
    raise utility.PMParamError('Invaid substance name.')







def info( name=None, contains=None, collection=None, pmclass=None, verbose=True):
    """Print information on substance data
    info()
        or
    info(name=None, contains=None, collection=None, pmclass=None, verbose=True)

When the _info_ funciton is called without arguments, all substances currently
loaded are listed in a table along with their class id string and their 
properties.

The optional keyword arguments are used to filter the list.  Those substances 
shown must match all of the criteria given.

name (string)
    Search for substances with id strings that contain the string given.  If
    one of the substances matches exactly, only it will be returned.

contains (string or list of strings)
    Search for substances that contain the elements given.  For example, 
    contains = 'H' returns substances that contain hydrogen, while 
    contains = ['C', 'H', 'O'] returns substances that contain carbon, 
    hydrogen, and oxygen.  Specifying 'e' as an element indicates the free
    electron.  If it is present in either positive or negative quantities, it 
    indicates that the species is an ion.

collection (string)
    Searches for substances that belong to a particular colleciton.  This 
    string must appear before the '.' in their substance id string.  For 
    example, collection = 'ig' will only return ideal gas collection members.

pmclass (string)
    Searches for substances that are instances of the class specified by the
    class identifier string specified.  For example, 'ig', 'ig2', and 'mp1'
    are all valid values.

verbose (bool)
    When True (True by default) the info() function prints a table to stdout
    for a user to read.  When False, the info() funciton returns a list of the
    substance id strings that matched the search criteria.
    
Chris Martin (c) 2015,2017,2021
"""

    month = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

    # First, check for an exact name match
    if name is not None:
        candidate = dat.data.get(name)
        if candidate is not None:
            members = [name]
        else:
            members = list(dat.data.keys())
    # Otherwise, start with a list of all available substances
    else:
        # Create a list of the available data
        members = list(dat.data.keys())
    # Apply the relevant filters
    # Start with filters most likely to narrow things down...
    if name is not None:
        ii = 0
        while ii < len(members):
            candidate = members[ii]
            if name not in candidate:
                del members[ii]
            else:
                ii += 1
    if pmclass is not None:
        ii = 0
        while ii < len(members):
            candidate = members[ii]
            if dat.data[candidate].data['class'] != pmclass:
                del members[ii]
            else:
                ii += 1
                    
    if collection is not None:
        ii = 0
        while ii < len(members):
            candidate = members[ii]
            if not candidate.startswith(collection):
                del members[ii]
            else:
                ii += 1
    if contains is not None:
        if isinstance(contains, str):
            contains = [contains]
        ii = 0
        while ii < len(members):
            atoms = dat.data[members[ii]].atoms()
            for species in contains:
                if species not in atoms:
                    del members[ii]
                    ii-=1
                    break
            ii += 1

    if not verbose:
        return members

    # If operating verbosely, we'll need to print to stdout
    target = utility.sys.stdout

    # If the filters eliminated everything
    if not members:
        target.write('No substances matched the critera.\n')
        if name:
            target.write('  name: ' + name + '\n')
        if collection:
            target.write('  collection: ' + collection + '\n')
        if contains:
            target.write('  contains: ' + str(contains) + '\n')
    # If there's only one member left
    elif len(members) == 1:
        name = members.pop()

        target.write( '***\nInformation summary for substance: "' + name + '"\n***\n' )

        out = 'Uses class:   ' + dat.data[name].data['class'] + '\n'
        target.write( out )

        # identify the file
        fil = dat.data[name].data['fromfile']
        target.write('\nLoaded from:  ' + fil + '\n')

        # get the time last modified
        tt = utility.time.localtime( utility.os.path.getmtime(fil))
        out = '\nLast updated: {0.tm_hour}:{0.tm_min:02d} {1:s} {0.tm_mday}, {0.tm_year}\n\n'.format( tt, month[tt.tm_mon-1] )
        target.write( out )

        utility.print_line(dat.data[name].__doc__,'')

    else:

        # Print a version summary
        target.write('  PYroMat\nThermodynamic computational tools for Python\n')
        target.write('version: ' + str(config['version']) + '\n')

        # Generate a table of existing data
        # Print the ID string, the date modified, and the file path
        
        # first, obtain a sorted list of the loaded data
        members.sort()
        # find the longest id string
        idlen = 2
        classlen = 5
        for ss in members:
            idlen = max(idlen, len(ss))
        # A list of properties for wich to search
        proplist = ['T', 'p', 'd', 'cp', 'cv', 'gam', 'e', 'h', 's', 'mw', 'R', 's', 'T_h', 'T_s', 'p_s', 'd_s', 'X', 'Y']
        proplen = 0
        for prop in proplist:
            proplen += len(prop) + 1

        fmt = ' {:<' + str(idlen) + 's} : {:^5s} :'
        head = '-'*(proplen + 11 + idlen) + '\n' + \
            fmt.format('ID','class') + ' properties\n' + \
            '-'*(proplen + 11 + idlen) + '\n'

        index = 0
        for ss in members:
            if index%15 == 0:
                target.write(head)
            index+=1

            target.write( fmt.format(ss, dat.data[ss].data['class']) )
            for prop in proplist:
                if hasattr(dat.data[ss],prop):
                    target.write( ' ' + prop )
                else:
                    target.write( ' ' * (len(prop)+1) )
            target.write('\n')
