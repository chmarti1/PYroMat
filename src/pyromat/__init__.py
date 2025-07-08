"""PYroMat    Thermodynamic property calculator for Python

PYroMat is an open-source Python-based software platform for retrieving
the physical properties of substances.  For complete documentation, 
visit "pyromat.org"

Chris Martin (c) 2015-2022
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
__version__ = "2.2.6"


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








def get( idstr ):
    """Return an object for a substance
    get( 'name' )

Returns a substance data class for the substance named.
"""
    # If accessing by the id string
    out = dat.data.get(idstr)
    if out is None:
        utility.print_error('No substanced named "' + str(idstr) + '" was found in the loaded data.')
        raise utility.PMParamError('Invalid substance ID string')
    return out
    


def search(name=None, contains=None, collection=None, pmclass=None, cas=None, inchi=None, members=None):
    """Returns a set of substance instances that match a set of search criteria
    members = search( ... )

The result is a Python set, which is intended to allow for convenient
union, intersect, and other set operations between sequential search
results.

The optional "members" keyword allows the results of a previous search
to be further refined.

    members = search(members=old_members, ...)
    
If no members argument is supplied, then the search is performed against
all substances loaded into PYroMat's data dictionary.  Substances are
excluded unless they match all of the criteria given in the arguments.

The parameters of the search are specified through keyword arguments:

** name **
The name keyword expects a string. To match, at least one of the 
following must be true:
    1) The supplied name string may appear inside the substance ID 
        string. For example, name='H2O' would match id strings 'ig.H2O'
        and 'ig.B3FH2O3' and 'ig.H2O+'
    2) The supplied name string may appear inside any one of the 'names'
        list of the substance's data.  If the substance has no 'names'
        list, then this criterion is ignored.  This comparison is NOT
        case sensitive, but this is not a Google search - check your
        spelling.

** contains **
The contains keyword is used to specify the atomic contents of the 
substance.  There are three ways for the contains keyword to be used:
    1) If the value is a string, it is interpreted as an atom that must
        be found in the substance.  For example, contains='H' matches
        all substances that contain hydrogen and none that don't.
    2) If the value is a list of strings, each element of the list must
        be an atom found in the substance.  For example,
        contains=['H','O'] would match water but neither hydrogen nor
        oxygen.
    3) If the value is a dictionary, its keys are interpreted as atom
        name strings and its values are the amount of each that must be
        in the substance.  For example, contains={'H':1, 'O':2} would
        match both H2O and C2H2O, but not H2O2.  To allow the quantity
        of one of the atoms to be variable, it may be specified as None.
        For example, contains={'H':1, 'O':None} would match H2O, C2H2O
        and H2O2.

** collection **
The collection keyword is a string that must match the prefix of the
substance id string.  For example, collection='ig' matches 'ig.H2O' but
not 'mp.H2O'.

** pmclass **
The PYroMat class (pmclass) keyword specifies the data class used by
the members.  For example, pmclass='ig2' only matches substances that 
use the ig2 data class.

** inchi **
and
** cas **
The InChI and CAS identifiers are strings that uniquely identify a 
substance, but the same substance may be modeled in multiple collections.
For example, water inchi='InChI=1S/H2O/h1H2' or cas='7732-18-5' is 
listed in both multi-phase and ideal gas collections.
"""
    # If the members are not specified, build a working set
    if members is None:
        members = dat.data.values()
    
    newmembers = set()
    
    # Also define a case-adjusted name for all-lower-case
    name_lower = None
    if name is not None:
        name_lower = name.lower()
    
    for candidate in members:
        match = True
        # Start by reducing the set size with the simplest comparisons
        # The class is easy
        if match and pmclass is not None and pmclass != candidate.data['class']:
            match = False
            
        # The collection is also easy
        if match and collection is not None and not candidate.data['id'].startswith(collection):
            match = False
            
        # Inchi and cas identifiers
        if match and cas is not None and ('cas' not in candidate.data or cas != candidate.data['cas']):
            match = False
            
        if match and inchi is not None and ('inchi' not in candidate.data or inchi != candidate.data['inchi']):
            match = False
            
        # Next, move on to the name.  It is more complicated, but it is 
        # likely to whittle down the search quickly
        if match and name is not None:
            # If the name is contained in the id string, match right away
            if name in candidate.data['id']:
                pass
            # If the candidate has names, check there too
            elif 'names' in candidate.data:
                # We'll assume there was no match until we find one
                match = False
                for cname in candidate.data['names']:
                    if name_lower in cname.lower():
                        match = True
                        break
            # This is not the substance you're looking for
            else:
                match = False
        
        # Condition the contents argument into an atoms dictionary
        if isinstance(contains, str):
            atoms = {contains:None}
        elif isinstance(contains, dict):
            atoms = contains
        elif hasattr(contains,'__iter__'):
            atoms = {}
            for aa in contains:
                atoms[aa] = None
        else:
            atoms = None
        
        # Finally, deal with the contents
        if match and atoms is not None:
            if hasattr(candidate, 'atoms'):
                catoms = candidate.atoms()
                # Check the atomic contents one-by-one
                # and verify that they are in the correct quantity
                for atm,qty in atoms.items():
                    cqty = catoms.get(atm)
                    if cqty is None or (qty is not None and cqty != qty):
                        match = False
                        break
            # If atoms were specified, but the class doesn't support it,
            # fail the search
            else:
                match = False
            
        # If the candidate was matched, add it to the set
        if match:
            newmembers.add(candidate)
                
    return newmembers
                    


def info(members=None, target=None, **kwarg):
    """Print information on substance data
    info()
        or
    info(idstr)
        or
    info(members)
        or
    info( .. search arguments .. )

The info() function prints information on individual PYroMat data instances
and formatted tables summarizing groups of them.  When it is called without
any arguments, it prints a table summarizing ALL substances currently loaded.

When it is called with an iterable of PYroMat data instances like the sets
generated by a call to search(), info() generates a table summarizing only
those.  

Alternately, if a single data instance or a single substance ID string is 
passed to info(), a display of detailed information on that substance is 
generated.

If the "members" argument is omitted, then the remaining keyword arguments
are used for a call to search().  The results are then used to generate
the table.

See search() for more information.

If the target keyword is set to a file stream, the output will be sent 
there instead of standard out.
"""

    month = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

    if members is None:
        members = search(**kwarg)
    elif isinstance(members, reg.__basedata__):
        members = [members]
    elif isinstance(members, str):
        members = [get(members)]

    # If operating verbosely, we'll need to print to stdout
    if target is None:
        target = utility.sys.stdout

    # If the filters eliminated everything
    if not members:
        target.write('No substances matched the critera.\n' + repr(kwarg) + '\n')
    # If there's only one member left
    elif len(members) == 1:
        ss = members.pop()

        target.write( '***\nInformation summary for substance: "' + ss.data['id'] + '"\n***\n' )

        # If the instance provides its atomic contents
        if hasattr(ss, 'atoms'):
            aline = '    '
            nline = '    '
            charge = ''
            for aa,qty in ss.atoms().items():
                # Special case for charge (e)
                if aa == 'e':
                    if qty>0:
                        charge = '-'
                    else:
                        charge = '+'
                    qty = abs(qty)
                    if qty!=1:
                        charge += str(qty)
                        
                else:
                    aline += aa
                    nline += ' '*len(aa)
                    if qty == 1:
                        qty = ''
                    elif isinstance(qty, int):
                        qty = str(qty)
                    else:
                        qty = '{:.3f}'.format(qty)
                    aline += ' '*len(qty)
                    nline += qty
            target.write('\n' + aline + charge + '\n' + nline + '\n\n')

        fmt = '{:>18s} : '

        if 'names' in ss.data:
            lead = fmt.format('Names')
            for name in ss.data['names']:
                target.write(lead + name + '\n')
                lead = ' '*len(lead)

        target.write(fmt.format('Molecular Weight') + str(ss.mw()) + '\n')
        
        # Check for CAS and InChI strings
        if 'cas' in ss.data:
            target.write(fmt.format('CAS number') + ss.data['cas'] + '\n')
        if 'inchi' in ss.data:
            target.write(fmt.format('InChI string') + ss.data['inchi'] + '\n')
        
        target.write(fmt.format('Data class') + ss.data['class'] + '\n')

        # identify the file
        fil = ss.data['fromfile']
        target.write(fmt.format('Loaded from') + fil + '\n')

        # get the time last modified
        tt = utility.time.localtime( utility.os.path.getmtime(fil))
        target.write( fmt.format('Last updated') + '{0.tm_hour}:{0.tm_min:02d} {1:s} {0.tm_mday}, {0.tm_year}\n\n'.format( tt, month[tt.tm_mon-1] ))
        
        utility.print_line(ss.data['doc'],'')

    else:

        # Print a version summary
        target.write('  PYroMat\nThermodynamic computational tools for Python\n')
        target.write('version: ' + str(config['version']) + '\n')

        # Generate a table of existing data
        # Print the ID string, the date modified, and the file path
        
        # first, obtain a sorted list of the loaded data
        members = list(members)
        members.sort(key=lambda ss: ss.data['id'])
        # find the longest id string
        idlen = 2
        classlen = 5
        namelen = 4
        for ss in members:
            idlen = max(idlen, len(ss.data['id']))
            if 'names' in ss.data and len(ss.data['names']) > 0:
                namelen = max(namelen, len(ss.data['names'][0]))
                
        # A list of properties for wich to search
        proplist = ['T', 'p', 'd', 'v', 'cp', 'cv', 'gam', 'e', 'h', 's', 'mw', 'R', 'X', 'Y']
        proplen = 0
        for prop in proplist:
            proplen += len(prop) + 1

        fmt = ' {:<' + str(idlen) + 's} : {:^' + str(classlen) + 's} : {:' + str(namelen) + 's} :'
        head = '-'*(idlen + classlen + namelen + proplen + 11) + '\n' + \
            fmt.format('ID','class','name') + ' properties\n' + \
            '-'*(idlen + classlen + namelen + proplen + 11) + '\n'

        index = 0
        for ss in members:
            if index%15 == 0:
                target.write(head)
            index+=1

            if 'names' in ss.data and len(ss.data['names'])>0:
                target.write( fmt.format(ss.data['id'], ss.data['class'], ss.data['names'][0] ))
            else:
                target.write( fmt.format(ss.data['id'], ss.data['class'], ''))
                
            for prop in proplist:
                if hasattr(ss,prop):
                    target.write( ' ' + prop )
                else:
                    target.write( ' ' * (len(prop)+1) )
            target.write('\n')
