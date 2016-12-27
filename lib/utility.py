"""PYroMat utility module

This is a collection of miscelaneous objects that are used by the 
package, but that users don't need to access explicitly.  These modules
and functions are hidden in a subordinate "utility" module to keep them
from crowding out the functions and modules that users will need most.

Developers interested in adding their own data or classes may need to 
learn what lies herein, but most users will never need to care.

Chris Martin (c) 2015,2017
"""

import json
import sys
import numpy as np
import os
import time
# point back to the root package
import pyromat as pyro






####################################
# Error handling
#   These are classes and functions
# for PYroMat error handling
####################################
class Error(Exception):
    def __init__(self, value=''):
        self.value = value
    def __str__(self):
        return repr(self.value)


# if the data dictionary seems to be corrupt
# This error is raised in the data do not include the required parameters or
# if the values in those parameters do not conform to PYroMat's requirements.
class PMDataError(Error):
    pass

# if there is an error loading files
# This error class is raised if file loading is inhibited by permission errors,
# file-not-found, or any other OS problems
class PMFileError(Error):
    pass

# if there is an illegal combination of parameters
# This error class is intended to be used if methods/functions are called with 
# parameters or values that don't make sense.
class PMParamError(Error):
    pass

# An anylitical algorithm has failed with an error
# This error class is reserved for high level algorithms that may experience
# failures; e.g. an iterative algorithm that fails to converge.
class PMAnalysisError(Error):
    pass



def load_config(filename = None,verbose=None):
    """Load the configuration files
    load_config()

Bootstraps through the default configuration file in the
installation directory, the global configuration files 
indicated by the 'config_file' directive.

Checks for unrecognized parameters and illegal values.
"""
    lead = 'load_config-> '

    # if load_config is called with a filename, this is a recursive call
    if filename:
        # if verbose is specified, override the parameter
        if verbose==None:
            verbose = get_config('config_verbose',dtype=bool)
            
        if not os.path.isfile( filename ):
            if verbose:
                print_line('Could not find config file: ' + os.path.abspath(filename), lead)
            return

        if verbose:
            print_line('Reading config file: ' + os.path.abspath(filename), lead)

        lead += '   '

        temp_config = {}
        with open(filename,'r') as ff:
            exec(compile(ff.read(),filename,'exec'),globals(),temp_config)
        # In Python 2.7, this used to work
        # Use exec() for 3.4 compatibility
        # execfile(filename,globals(),temp_config)

        # check each parameter
        for param in temp_config:
            # if this is a new parameter, create an entry for it
            if not param in pyro.config:
                pyro.config[param] = []
                if verbose:
                    print_line('New parameter: "' + param + '"',lead)

            # if the parameter is writeable
            if isinstance(pyro.config[param],list):
                good = True
                if not isinstance(temp_config[param],list):
                    temp_config[param] = [ temp_config[param] ]
                for testme in temp_config[param]:
                    if not isinstance(testme, (str,int,float,bool)):
                        good = False
                        print_warning(
'Bad parameter value for "' + param + '" found in config file "' + filename + 
'."  Values must be strings, integers, or floats, or lists of strings, integers, or floats.  Ignoring this assignment.')
                        break
                # if there were no problems with the value
                if good:
                    pyro.config[param] += temp_config[param]
                    if verbose:
                        print_line('Set parameter: "' + param + '"',lead)
                        print_line('to value: ' + repr(temp_config[param]),lead)
            else:
                # If the parameter is not appendable
                print_warning('An illegal assignment to the "' + param + 
                    '" parameter was made in config file "' + filename + 
                    '.  This parameter is read only.')


    else:
        # if loadconfig is called without a filename, then this is the root call

        # assign default configuration parameters
        install_dir = os.path.dirname( pyro.__file__ )
        # Below was the installation directory deteciton method for version 1.2
        # The path attribute is controlled by the import system, but the master
        # __init__ file location is under the pyromat package's control.
        # install_dir = os.path.abspath(pyro.__path__[0])
        default_config = os.path.join( install_dir, 'defaults.py')

        pyro.config = {
            'install_dir': install_dir,
            'version' : pyro.__version__,
            'config_file' : [default_config],
            'config_verbose' : [False],
            'dat_dir' : [os.path.join( install_dir, 'data')],
            'dat_verbose' : [True],
            'dat_overwrite' : [True],
            'dat_exist_fatal' : [False],
            'dat_recursive' : [True],
            'reg_dir' : [os.path.join( install_dir, 'registry')],
            'reg_verbose' : [True],
            'reg_overwrite' : [True],
            'reg_exist_fatal' : [False],
            'def_T' : [300.],
            'def_p' : [1.01325],
            }

        if not os.path.exists(default_config):
            try:
                df = open(default_config,'a')
                df.close()
            except:
                print_error(
'The default config file was missing.  Unable to create one at "' + 
default_config + '".  Check for a permissions error.')
                raise PMFileError(default_config)

        try:
            load_config( default_config, verbose=verbose )
        except:
            # if the default config file fails for any reason
            print_error(
'There was an error loading the default configuration file, "' + 
default_config + '". Check its permissions and syntax.')
            raise PMFileError(default_config)

        
        k = 1
        found = [default_config]
        while k<len(pyro.config['config_file']):
            thisfile = os.path.abspath(pyro.config['config_file'][k])

            if thisfile in found:
                print_warning(
'Ignoring a repeated reference to config file: "' + thisfile + '"')
            else:
                load_config( pyro.config['config_file'][k] )
                found.append( thisfile )
            k+=1






def get_config( param, dtype=None, verbose=True ):
    """All config parameters that are writeable are contained
in lists.  The load_config() function simply appends
new values rather than overwriting old values, which
leaves a complete record of all values written since
the package was loaded.

For parameters that are only permitted one value, the
last one is the only one honored, but routines need to
be protected against the case that a user may overwrite
the list with a single value after load.

The get_config() function does its best to interpret 
the data stored in a config parameter and return a single
value.

If the dtype keyword is specified, it is treated as a
class to which to force the answer dtype(value)
"""
    if not param in pyro.config:
        if verbose:
            print_warning(
'Parameter "' + param + '" does not appear in the PYroMat configuration file.')
        return None

    if isinstance(pyro.config[param], list):
        done = pyro.config[param][-1]
    elif isinstance(pyro.config[param], (int,str,float,bool)):
        done = pyro.config[param]
    else:
        if verbose:
            print_warning(
'Could not retrieve a valid value for the PYroMat configuration parameter "' + 
param + '."')
        return None

    if dtype:
        return dtype(done)
    return done





def split_lines( text, lead='', tail='', width=74):
    """Split a string across multiple lines without 
dividing up words.  All combinations of whitespace 
characters are interpreted as a single space except
repeated newlines, which represent a paragraph break.

Optional keywords are:

'lead'  (def: lead='')
Indicates a string that will be inserted at the
beginning of each line.

'tail'  (def: tail='')
A string that will be inserted at the end of each
line.

'width' (def: width=74)
The maximum line width in characters.  If the lead 
and the tail add to more than 'width' characters,
then 'print_lines' returns -1.
"""

    NL = len(lead)
    NT = len(tail)
    tail += '\n'
    # adjust the width by the lead and tail
    Nmax = int(width - NL - NT)
    # make sure Nmax is positive
    if Nmax<=0:
        return -1

    out = lead

    pars = text.split('\n\n')
    for par in pars:
        words = par.split()
        Nline = 0
        
        for word in words:
            NW = len(word)
    
            # if the word is longer than a line, split it up
            while NW>Nmax:
                # if there is already text on the line
                # start a new line
                if Nline>0:
                    out+=tail
                out += lead + word[:Nmax] + tail + lead
                word = word[Nmax:]
                NW -= Nmax
                Nline = 0
            # if the word will not fit on the current line
            if Nline+1+NW>Nmax:
                out += tail + lead + word
                Nline = NW
            # if the current line is empty, omit the leading space
            elif Nline==0:
                out+=word
                Nline=NW
            # otherwise, add a space and the word
            else:
                out += (' '+word)
                Nline += 1+NW
        out += tail

    return out


def print_error(text):
    sys.stdout.write(split_lines(text,lead='PYroMat ERR:: '))

def print_warning(text):
    sys.stdout.write(split_lines(text,lead='PYroMat WARN:: '))

def print_line(text, lead):
    sys.stdout.write(split_lines(text,lead))








def load_file(filename,test=True,verbose=True):
    """Load a single file and return the dictionary
    loadfile(filename)
Suppress the default data check by setting the 'test'
keyword to False.
Suppress printing errors to stdout by setting the 
'verbose' keyword to False.
"""

    # open the file
    try:
        fil = open(filename,'r')
    except:
        print_error(
'Failed to open file ' + repr(filename) + 
'. The file does not exist, or there may be a permissions problem.')
        raise PMFileError(filename)

    # parse the file
    try:
        readin = json.load(fil)
    except:
        print_error('Could not parse file ' + repr(filename) + 
            '. Try replacing the file from a fresh installation.')
        fil.close()
        raise PMFileError(filename)
    fil.close()


    if test:
        failure=False
        musthave = ['id', 'doc', 'class']

        # check for essential contents
        for mh in musthave:
            if not mh in readin:
                failure=True
                if verbose:
                    print_error(
'File ' + repr(filename) + ' is missing the essential entry, ' + 
repr(mh) + '.')
        if failure:
            raise PMFileError(filename)

    return readin












def suppress_file( path, verbose=True ):
    """Rename a *.hpd file to suppress it
    suppress_file( 'path/to/file.hpd' )

Renames a file by appending a '~' after the extension.
The modification will cause load() to disregard the 
file, and most OSs will also hide the file.  Once a 
file is suppressed, it can be found by running 
load(check=True), which returns a list of suppressed
files in the search path as part of its output.
"""
    # Check that the file is a data file, and that it is
    # not already suppressed
    if (len(path)<=4) or (path[-4:]!='.hpd'):
        print_error('File ' + path + ' is not a *.hpd file.')
        raise PMFileError('Suppression failed.')
    try:
        os.rename(path,path+'~')
        if verbose:
            sys.stdout.write('Suppressed file ' + path + '\n')
    except:
        print_error('Failed to rename file ' + path + '.  Check permissions.')
        raise PMFileError('Suppression failed.')









def revive_file( path, verbose=True ):
    """Rename a *.hpd~ file to un-suppress it
    suppress_file( 'path/to/file.hpd~' )

Renames a file by removing the '~' after the extension.
The modification will allow load() to see the file.
"""
    # Check that the file is a data file, and that it is
    # already suppressed
    if (len(path)<=5) or (path[-5:]!='.hpd~'):
        print_error('File ' + path + ' is not a suppressed *.hpd file.')
        raise PMFileError('Revival failed.')
    try:
        os.rename(path,path[:-1])
        if verbose:
            sys.stdout.write('Revived file ' + path + '\n')
    except:
        print_error('Failed to revive file ' + path + '.  Check permissions.')
        raise PMFileError('Revival failed.')









def red_repair(RED=None, verbose=True, action='s', select=0):
    """Repair reduntant file conflict
    red_repair()
        or
    red_repair(RED, verbose=True, action='s', select=0)

Repairs data file redundancy conflicts by either
deleting or suppressing all but one of the files
listed in conflict.  RED is supposed to be the 
'redundancy' entry of the output of the load()
function run with the check=True option set. If
RED is left None, then red_repair() will call
load() itself and continue.

If run verbosely, the action and select options
will be overridden by prompted user inputs.  If
run quietly, the ACTION option accepts any of
three case-insensitive characters; 's'uppress, 
'd'elete, or 'i'gnore.

Data files are suppressed by appending a '~' 
character after their extension so that the load
function will disregard them.  Ignoring files
renders the red_repair() function useless.

The select index is an integer indicating which 
file in the list will be left alone.
"""

    # If RED is unspecified, call load()
    if RED == None:
        status = pyro.dat.load(check=True,verbose=False)
        RED = status['redundant']


    # loop through the reported problems    
    for ss in RED:

        ######
        #
        # construct a menu to prompt the user for options
        #
        ######
        if verbose:
            # Print a statement of the conflict including a numbered list of files
            N = len(RED[ss])
            sys.stdout.write('\nREDUNDANCY found for entry "{:s}"\n'.format(ss))
            for index in range(N):
                sys.stdout.write( '({:d})   {:s}\n'.format(index,RED[ss][index]) )

            #sys.stdout.write('='*35 + '\n')
            # print acceptable options for addressing the conflict
            # The menus are a set of nested loops
            # They execute until the user supplies acceptable inputs
            f_go = True
            while f_go:
                #
                # prompt the user for an action
                sys.stdout.write(
"""  ++
(I)gnore:       The conflict will remain unresolved
(D)elete:       Delete redundant files
(S)uppress:     Modify filenames to prevent loading
""")
                while f_go:
                    action = raw_input('(I/D/S):').lower()
                    # accept the input when it is i, d, or s
                    f_go = not (action in 'ids')
                    # if the input is unacceptable, scold the user
                    if f_go:
                        sys.stdout.write(
'   Please enter ''I'', ''D'', or ''S''\n')
                #
                # prompt the user for a file selection
                #
                f_go = (action in 'ds')
                while f_go:
                    sys.stdout.write('  ++\n')
                    out = 'Select a file to KEEP by its index above.  '
                    if action == 'd':
                        out+= 'All other files will be DELETED.  '
                    else:
                        out+= 'All other files'' will be RENAMED to suppress them during load.  '
                    out+= 'Your selection will be executed immediately upon entry.'
                    print_line(out,'')
                    sys.stdout.write('(B)ack:         Returns to the previous menu.\n')
                    select = raw_input('(0-{:d}/B):'.format(N-1)).lower()

                    if select == 'b':
                        f_go = False
                    else:
                        try:
                            select = int(select)
                        except:
                            select = -1
                            sys.stdout.write('! Please enter an integer.\n')
                        
                        # make sure the selection is in range
                        f_go = (select < 0) or (select >= N)
                        if f_go:
                            sys.stdout.write('! Please make a selection in the range specified.\n')
                    
                # keep executing the menu loop if the selection was 'back'
                f_go = select == 'b'
                if f_go:
                    select = 0

        # if the action is 'delete'
        if action == 'd':
            # execute the repair
            for index in range(N):
                if index != select:
                    try:
                        fil = RED[ss][index]
                        os.remove(fil)
                        if verbose:
                            sys.stdout.write(
'Removed file "' + fil + '"\n')
                    except:
                        print_warning(
'Failed to remove file "' + fil + 
'". Check permissions and check that the file is not open.')
        elif action == 's':
            # execute the repair
            for index in range(N):
                if index != select:
                    try:
                        fil = RED[ss][index]
                        suppress_file(fil,verbose=verbose)
                    except:
                        print_warning(
'Failed to suppress file "' + fil + 
'". Check permissions and check that the file is not open.')




    
def newtoniter(f,g,T,p,fval=0.,gval=0.,fprec=1e-3,gprec=1e-3,N=100):
    """Perform newton iteration on f(T,p) and g(T,p)
    (T,p) = newtoniter(f,g,T,p,fval=0.,gval=0.)
    
Solves the system of equations
  fval = f(T,p)
  gval = g(T,p)
given an initial guess, T0 and p0.  The desired precision of the solution
is configurable through optional keyword parameters, 'fprec' and 'gprec',
which specify property precision in their respective units.  These 
default to 1e-3, which may not be appropriate to all properties.

The maximum number of iterations is limited by the parameter 'N', and is
100 by default.
"""

    Fval = np.array([fval,gval])
    Tpert = 0.1
    ppert = .001    
    
    def F(T,p):
        return (f(T,p)-fval, g(T,p)-gval)
    
    (f0,g0) = F(T,p)
    # test the new approximation
    # The old error
    e0 = f0*f0 + g0*g0

    go = True
    count = 0
    while go and count<N:
        count += 1
        T1 = T+Tpert
        p2 = p+ppert
        (f1,g1) = F(T1,p)
        (f2,g2) = F(T,p2)
        
        J = np.array([[ (f1-f0)/Tpert, (f2-f0)/ppert], [(g1-g0)/Tpert, (g2-g0)/ppert]])
        try:
            # calculate the motion of the approximation
            dX = np.linalg.solve(J,-np.array((f0,g0)))
        except:
            raise PMParamError('The constraints appear to be incompatible.')
        
        
        Ttest = T + dX[0]
        ptest = p + dX[1]
        (f0,g0) = F(Ttest,ptest)
        # the new error
        etest = f0*f0 + g0*g0

        # if the iteration does not pass the error test
        # Using "not <" instead of ">=" is deliberate
        # If etest is nan or inf or -inf, it will fail any comparison
        # Using ">=" would cause the iteration to halt incorrectly
        # Using "not <" will continue iteration until a valid value is found
        while ( not (etest < e0) and count<N):
            count += 1
            dX /= 2.
            Ttest = T + dX[0]
            ptest = p + dX[1]
            (f0,g0) = F(Ttest,ptest)
            etest = f0*f0 + g0*g0

        # update the status variables        
        T = Ttest
        p = ptest
        e0 = f0*f0 + g0*g0

        
        go = abs(f0)>=fprec or abs(g0)>=gprec
        
    if count==N:
        raise PMAnalysisError('newtoniter() failed to converge.')
    return (T,p)
