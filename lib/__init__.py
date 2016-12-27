"""PYroMat    Thermodynamic property calculator for Python

PYroMat includes a user-expandable library of thermodynamic
data on species and a suite of functions for working with 
them.  That includes calculating properties of pure sub-
stances and mixtures.

Chris Martin (c) 2015, 2017
Released under the GNU General Publice License v3.0
  http://www.gnu.org/licenses/gpl-3.0.en.html

*** Available data ***

Upon import, PYroMat automatically sets out to find its
constituent data files, and loads them into memory. Each
file contains data on a single substance or mixture, and
a summary of files discovered is printed.

For information on a single species or to print a summary
of available data, see the info() function.

*** Retrieving data ***

To get started, retrieve a substance that is of interest
and commit it to a variable that you will use later. The
command below retrieves an object for Argon.
  >>> import pyromat as pyro
  >>> Ar = pyro.get('Ar')

Once created, these objects can be called on to recover
thermodynamic properties given a temperature in K and
a pressure in bar.  The command below returns the 
enthalpy of argon at 325K and atmospheric pressure.
  >>> h = Ar.h(T=325, P=1.013)

For all data, methods exist for 
    cp  constant-pressure specific heat (kJ/kg/K)
    cv  constant-volume specific heat (kJ/kg/K)
    d   density (kg/m^3)
    e   internal energy (kJ/kg)
    h   enthalpy (kJ/kg)
    k   specific heat ratio (d-less)
    mw  molecular weight (kg/kmol)
    s   entropy (kJ/kg/K)

Functions accept arguments
    T   temperature (K)
    P   pressure (bar)

For some data, methods also exist to calculate
    R   ideal gas constant (kJ/kg/K)

The methods all support vector and array arguments.
  >>> T = numpy.arange(300,1000,100)
  >>> P = numpy.arange(.5,1.2,.1)
  >>> h = Ar.h( T )
  >>> d = Ar.d( T, P )
  >>> s = Ar.s( T, 1.013 )
  >>> mw = Ar.mw()

Note that the pressure argument is sometimes absent.  
Arguments that do not affect the output of a function
are left as optional.

The shape of the T and P arrays must be identical unless
one of them is a scalar.

*** Structure ***

All functions to which users are expected to need direct
access are exposed at the root level of the PYroMat package,
so that >>> dir(pyro) should expose all the functions 
typically needed.

There are four subordinate modules that obfuscate the more
basic functionality of the package; constants, dat, reg, 
and utility.  Each performs a specific service.

pyro.dat is a module for handling the substance data loaded
by the package.  Therein is the dictionary containing objects
for each of the available substances.  These are what is
returned by the pyro.get() function.  Functions to help
users customize their data directories are located here.
For more, see the pyro.dat documentation.

pyro.reg is a module for handling the PYroMat object class 
definitions in a modular but transparent way.  PYroMat is 
designed so that users can create their own data files and
even their own data classes without needing to modify the 
PYroMat code.  The class registry located in pyro.reg is
where these classes reside.  They are called on to create
the class objects found in the pyro.dat data dictionary.

pyro.utility is a collection of low-level functions to 
which the vast majority of users will never need access.
They are housed in the utility module to prevent them from
cluttering out more commonly used functions.
"""

# This is the athoritative version number.
# utility.load_config() checks this value to establish the read-only version
# setup.py looks for this line to establish the version at install
# MUST be unindented
__version__ = "1.3"


# loading the PYroMat utility functions
from . import utility

config = {}
utility.load_config()

# import the dataclass registry
from . import reg
# import the module for handling data
from . import dat










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
        ids = dat.data.keys()
        ids.sort()
        # find the longest id string
        idlen = 2
        for ss in ids:
            idlen = max(idlen, len(ss))

        # build a format string
        fline = '{:<' + str(idlen+2) + 's}{:<12s}{:s}\n'

        head = '-'*(idlen+25) + '\n'
        head = head + fline.format('ID','Modified', 'Type') + head
        index = 0

        for ss in ids:
            if not (index % 15):
                utility.sys.stdout.write( head )
            index += 1
            fil = dat.data[ss].data['fromfile']
            tt = utility.time.localtime( utility.os.path.getmtime(fil) )
            dstr = '{:d}/{:d}/{:d}'.format(tt.tm_mon, tt.tm_mday, tt.tm_year)
            utility.sys.stdout.write( fline.format(ss, dstr, dat.data[ss].data['class']))


