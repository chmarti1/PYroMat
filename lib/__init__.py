"""PYroMat    Thermodynamic property calculator for Python

PYroMat is an open-source Python-based software platform for retrieving
the physical properties of substances.  For complete documentation, 
visit "pyromat.org"

Chris Martin (c) 2015, 2017, 2018
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
"""

# This is the athoritative version number.
# utility.load_config() checks this value to establish the read-only version
# setup.py looks for this line to establish the version at install
# MUST be unindented
__version__ = "2.0.9"


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







def info( name=None ):
    """Print information on substance data
    info()
        or
    info( 'name' )

Lists all substances currently loaded in a table.
If called with a string indicating one of the 
substances, info() prints relevant information 
for the substance named.
"""

    month = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

    if name:
        if not (name in dat.data):    
            utility.print_error('No substance named "' + str(name) + '" was found in the loaded dat.')
            raise utility.PMParamError('Invaid substance name.')

        utility.sys.stdout.write( '***\nInformation summary for substance: "' + name + '"\n***\n' )

        out = 'Uses class:   ' + dat.data[name].data['class'] + '\n'
        utility.sys.stdout.write( out )

        # identify the file
        fil = dat.data[name].data['fromfile']
        utility.sys.stdout.write('\nLoaded from:  ' + fil + '\n')

        # get the time last modified
        tt = utility.time.localtime( utility.os.path.getmtime(fil))
        out = '\nLast updated: {0.tm_hour}:{0.tm_min:02d} {1:s} {0.tm_mday}, {0.tm_year}\n\n'.format( tt, month[tt.tm_mon-1] )
        utility.sys.stdout.write( out )

        utility.print_line(dat.data[name].__doc__,'')

    else:

        # Print a version summary
        utility.sys.stdout.write('  PYroMat\nThermodynamic computational tools for Python\n')
        utility.sys.stdout.write('version: ' + str(config['version']) + '\n')

        # Generate a table of existing data
        # Print the ID string, the date modified, and the file path
        
        # first, obtain a sorted list of the loaded data
        ids = list(dat.data.keys())
        ids.sort()
        # find the longest id string
        idlen = 2
        for ss in ids:
            idlen = max(idlen, len(ss))
        # A list of properties for wich to search
        proplist = ['cp', 'cv', 'd', 'e', 'gam', 'h', 'k', 'mw', 'p_d', 'p_s', 'R', 's', 'T_d', 'T_h', 'T_s', 'X', 'Y']
        proplen = 0
        for prop in proplist:
            proplen += len(prop) + 1

        fmt = ' {:>' + str(idlen) + 's} :'
        head = '-'*(proplen + 3 + idlen) + '\n' + \
            fmt.format('ID') + ' properties\n' + \
            '-'*(proplen + 3 + idlen) + '\n'

        index = 0
        for ss in ids:
            if index%15 == 0:
                utility.sys.stdout.write(head)
            index+=1

            utility.sys.stdout.write( fmt.format(ss) )
            for prop in proplist:
                if hasattr(dat.data[ss],prop):
                    utility.sys.stdout.write( ' ' + prop )
                else:
                    utility.sys.stdout.write( ' ' * (len(prop)+1) )
            utility.sys.stdout.write('\n')
