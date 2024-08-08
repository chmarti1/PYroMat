"""PYroMat utility module

This is a collection of miscelaneous objects that are used by the 
package, but that users don't need to access explicitly.  These modules
and functions are hidden in a subordinate "utility" module to keep them
from crowding out the functions and modules that users will need most.

Developers interested in adding their own data or classes may need to 
learn what lies herein, but most users will never need to care.

Chris Martin (c) 2015,2017,2021
"""

import json
import sys
import numpy as np
import os
import traceback as tb
import time
# point back to the root package
import pyromat as pm






####################################
# Error handling
#   These are classes and functions
# for PYroMat error handling
####################################
class Error(Exception):
    pass

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


class PMConfigEntry:
    """PYroMat Configuration Entry

This class stores various meta information about configuration entries
such as read/write permissions, entry types, and default values.

The behavior of a configuration entry is defined by four parameters:

default     The value assigned to this parameter will be used to 
            initialize the entry.

append      When False, writing to the entry will overwrite the previous
            value.  When True, the entry's value will be a list, and all
            write operations will append new values to the list.

write       When False, an error will be raised if the entry's write()
            method is called.

etype       When this parameter is not None, it will be treated as a 
            type cast to convert values
            self.value = etype(value)
            For example, str, int, and float are common etype entries.
"""
    def __init__(self, default=None, append=False, write=True, etype=None):
        # Temporarily allow writing regardless of the requested mode
        self.write_allowed = True
        self.etype = etype
        self.append = append
        # The type requirement is enforced when the default is stored
        # The set_default() method also respects this behavior.
        self.set_default(default)
        self.restore_default()
        self.write_allowed = write

    def __repr__(self):
        return 'PMConfigEntry(' + repr(self.value) + ')'

    def write(self, value):
        """Write a value to the entry
    pmentry.write(value)
"""
        # Enforce write permissions
        if not self.write_allowed:
            raise PMParamError('Entry is read-only.')

        # Deal with appended iterables
        if self.append and isinstance(value, (list,tuple)):
            for vv in value:
                self.write(vv)
            return

        # Enforce type specifier
        if self.etype:
            try:
                tvalue = self.etype(value)
            except:
                raise PMParamError('Expected %s, but got %s'%(repr(self.etype), repr(value)))
        else:
            tvalue = value

        # Append or overwrite?
        if self.append:
            self.value.append(tvalue)
        else:
            self.value = tvalue
        

    def read(self):
        """Read a value from the entry
    pmentry.read()

This is equivalent to
    pmentry.value
"""
        return self.value


    def restore_default(self):
        """Set the value to its original default"""
        # Enforce write permissions
        if not self.write_allowed:
            return
        # Append?
        if self.append:
            self.value = [self.default]
        else:
            self.value = self.default


    def set_default(self,default):
        """Set the entry's default value
    pmentry.set_default( new_default )

This operation is only allowed on writable entries.
"""
        # Enforce write permissions
        if not self.write_allowed:
            raise PMParamError('Default cannot be set. Entry is read-only.')
        # Enforce type specifier
        if self.etype:
            try:
                self.default = self.etype(default)
            except:
                raise PMParamError('Expected %s, but got %s'%(repr(self.etype), repr(default)))
        else:
            self.default = value

    
class PMConfig:
    """PYroMat Configuration Class

This behaves much like a dictionary that enforces the PYroMat 
configuration rules.  To read or modify configuration parameters, access
them 
    config['unit_temperature'] = 'F'

if config['version'].split('.')[0] < 2:
    raise Exception('PYroMat is old.')

As of version 2.2.0, the PYroMat config instance also supports 
iteration.  For example, this prints all entry names and their current
values:

    for key in config:
        print(key, config[key])
        
"""
    def __init__(self, load=True):
        # detect the package installation directory
        install_dir = os.path.dirname( pm.__file__ )
        install_dir = os.path.abspath(install_dir)
        # Below was the installation directory deteciton method for version 1.2
        # The path attribute is controlled by the import system, but the master
        # __init__ file location is under the pyromat package's control.
        # install_dir = os.path.abspath(pm.__path__[0])

        # Point to the default configuration file
        default_config = os.path.join( install_dir, 'config.py')
        # directories
        data_dir = os.path.join( install_dir, 'data')
        reg_dir = os.path.join( install_dir, 'registry')

        self.entries = {
            'install_dir': PMConfigEntry(default=install_dir, write=False, etype=str),
            'version' : PMConfigEntry(default=pm.__version__, write=False, etype=str),
            'config_file' : PMConfigEntry(default=default_config, append=True, etype=str),
            'config_verbose' : PMConfigEntry(default=False, etype=bool),
            'warning_verbose' : PMConfigEntry(default=True, etype=bool),
            'error_verbose' : PMConfigEntry(default=True, etype=bool),
            'dat_dir' : PMConfigEntry(default=data_dir, append=True, etype=str),
            'dat_verbose' : PMConfigEntry(default=True, etype=bool),
            'dat_overwrite' : PMConfigEntry(default=True, etype=bool),
            'dat_exist_fatal' : PMConfigEntry(default=False, etype=bool),
            'dat_recursive' : PMConfigEntry(default=True, etype=bool),
            'reg_dir' : PMConfigEntry(default=reg_dir, append=True, etype=str),
            'reg_verbose' : PMConfigEntry(default=True, etype=bool),
            'reg_overwrite' : PMConfigEntry(default=True, etype=bool),
            'reg_exist_fatal' : PMConfigEntry(default=False, etype=bool),
            'def_T' : PMConfigEntry(default=298.15, etype=float),
            'def_T_unit' : PMConfigEntry(default='K', etype=str),
            'def_p' : PMConfigEntry(default=1.01325, etype=float),
            'def_p_unit' : PMConfigEntry(default='bar', etype=str),
            'def_oob' : PMConfigEntry(default=np.nan, etype=float),
            'unit_force' : PMConfigEntry(default='N', etype=str),
            'unit_energy' : PMConfigEntry(default='kJ', etype=str),
            'unit_temperature' : PMConfigEntry(default='K', etype=str),
            'unit_pressure' : PMConfigEntry(default='bar', etype=str),
            'unit_molar' : PMConfigEntry(default='kmol', etype=str),
            'unit_volume' : PMConfigEntry(default='m3', etype=str),
            'unit_length' : PMConfigEntry(default='m', etype=str),
            'unit_mass' : PMConfigEntry(default='kg', etype=str),
            'unit_time' : PMConfigEntry(default='s', etype=str),
            'unit_matter' : PMConfigEntry(default='kg', etype=str),
        }
        if load:
            self.load()

    def __repr__(self):
        out = ''
        # Detect the longest parameter name
        justify = 0
        for k in self.entries:
            justify = max(justify,len(k))
        # Build a format string that will right-align the parameter names
        fmt = '%' + str(justify) + 's : %s\n'
        # And print each configuration value on a line, sorted by parameter
        parameters = list(self.entries.keys())
        parameters.sort()
        for k in parameters:
            # Construct a representation of the value
            temp = repr(self.entries[k].value)
            # If it is too long to fit on a line, cut it off
            if len(temp)+justify > 72:
                temp = temp[:66-justify]+'...'
            out += fmt%(k,temp)
        return out


    def load(self,filename = None,verbose=None):
        """Load the configuration files
    pmconfig.load()

Bootstraps through the default configuration file in the
installation directory, the global configuration files 
indicated by the 'config_file' directive.

Checks for unrecognized parameters and illegal values.
"""
        lead = 'PMConfig-> '

        # if load is called with a filename, this is a recursive call
        if filename:
            # if verbose is specified, override the configuration parameter
            if verbose is None:
                verbose = self['config_verbose']
                
            if not os.path.isfile( filename ):
                if verbose:
                    print_line('Could not find config file: ' + os.path.abspath(filename), lead)
                return

            if verbose:
                print_line('Reading config file: ' + os.path.abspath(filename), lead)

            lead += '   '

            temp_config = {}
            with open(filename,'r') as ff:
                exec(ff.read(),{},temp_config)
            # In Python 2.7, the code below used to work
            # Use exec() for 3.4 compatibility
            # Also, allowing access to globals() is a security problem.
            # execfile(filename,globals(),temp_config)

            self.update(temp_config)
            if verbose:
                message = 'Found entries: '
                for item in temp_config:
                    message += item + ' '
                print_line(message, lead)

        else:
            # if loadconfig is called without a filename, then this is the root call
            k = 0
            found = []
            # As the files are loaded in, the config_file list will grow,
            # so the length needs to be detected after each iteration.
            # This method allows configuration files to be read in tiers;
            # Those that are identified in the earlier config files will be
            # read first and will have low precedence.  Those that are 
            # identified deeper in a series of config files will be read last
            # and will have high precedence.
            cfiles = self.entries['config_file'].value
            while k<len(cfiles):
                # Expand references to the users' home directories
                # and environment variables
                thisfile = os.path.expanduser(cfiles[k])
                thisfile = os.path.expandvars(thisfile)
                thisfile = os.path.abspath(thisfile)

                # If the config file is already in the found list, do not
                # load it.  This could result in an endless loop.
                if thisfile in found:
                    print_warning(
    'Ignoring a repeated reference to config file: "' + thisfile + '"')
                else:
                    self.load( thisfile )
                    found.append( thisfile )
                k+=1


    def update(self, new):
        """Read in parameters from a dictionary
    pmconfig.update( new_dictionary )

This function is equivalent to
for item in new_dictionary:
    pmconfig[item] = new_dictionary[item]
"""
        for item in new:
            self[item] = new[item]

    def restore_default(self, item=None):
        """Return a parameter to its default value
    pmconfig.restore_default(item)
        Or
    pmconfig.restore_default()

When item is a configuration parameter string, that configuration 
parameter is restored to its default value.  When item is omitted,
all parameters are restored to their default.

This function is equivalent to
pmconfig.entries[item].apply_default()
"""
        if item is None:
            for item, entry in self.entries.iteritems():
                entry.restore_default()
        else:
            self.entries[item].restore_default()


    def __getitem__(self, item):
        """Return the configuration value for an item
    value = config[item]
"""
        if item not in self.entries:
            raise PMParamError('%s is not a PYroMat configuration parameter'%repr(item))
        return self.entries[item].value


    def __setitem__(self,item,value):
        """Set the value of a configuration parameter
    config[item] = value
"""
        if item not in self.entries:
            raise PMParamError('%s is not a PYroMat configuration parameter'%repr(item))
        try:
            self.entries[item].write(value)
        except:
            print_error('Failed to write to configuration parameter, %s'%repr(item))
            tb.print_exception(*sys.exc_info())
        
    def __contains__(self,item):
        return self.entries.__contains__(item)

    def __iter__(self):
        return self.entries.__iter__()

    def def_T(self):
        """Return the default temperature in the current units
    T = def_T()
    
Uses the def_T_unit and def_T configuration fields to convert the default
temperature into the currently configured temperature units.
"""
        return pm.units.temperature_scale(self['def_T'], from_units = self['def_T_unit'])
        
    def def_p(self):
        """Return the default pressure in the current units
    p = def_p()
    
Uses the def_p_unit and def_p configuration fields to convert the default
temperature into the currently configured temperature units.
"""
        return pm.units.pressure(self['def_p'], from_units = self['def_p_unit'])
        
        


def get_config( param, dtype=None, verbose=True ):
    """**DEPRECIATED**
This is only left for reverse compatibility.  The present implementation
permits direct access to pyro.config.

Equivalent to 
pm.config[param]

dtype and verbose keywords are now ignored.
"""
    return pm.config[param]





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
    if pm.config['error_verbose']:
        sys.stdout.write(split_lines(text,lead='PM ERR: '))

def print_warning(text):
    if pm.config['warning_verbose']:
        sys.stdout.write(split_lines(text,lead='PM WARN: '))

def print_line(text, lead):
    sys.stdout.write(split_lines(text,lead))








def load_file(filename,test=True,verbose=True):
    """Load a single data file and return the dictionary
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
        status = pm.dat.load(check=True,verbose=False)
        RED = status['redundant']

    # Choose an input function based on the Python version
    if sys.version_info.major == 2:
        minput = raw_input
    else:
        minput = input

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
                    action = minput('(I/D/S):').lower()
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
                    select = minput('(0-{:d}/B):'.format(N-1)).lower()

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
            N = len(RED[ss])
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




def proptest(fn, kwarg, truth, ep, text, report, percent=True, findex=None):
    """Test a property method against truth values
    result = proptest(fn, kwarg, truth, ep, text report, percent=True)
    
PROPTEST is a utility function for performing numerical integrity checks on
property methods.  Returns True if the test is successful and False if not.

    fn
The method under test.  For functions that return tuples of array values 
(like saturation properties), the optional findex index should be set.  It 
is treated as the index of the element in the tuple that will be tested.

    kwarg
The argument tuple and keyword argument dictionary to pass to fn

    truth
The array of truth values expected from fn(**kwarg).
    
    ep
The fractional permissible error.  By setting the optional percent keyword to
False, ep is interpreted as absolute error instead of fractional.
    
    text
A string (usually a single-line) to insert in the report detailing the test 
being run.  A successful test appears in the reort as:
    "[passed]    __string__inserted__here__"
A failed test appears in the report as:
    "[FAILED]    __string__inserted__here__"
    " ... postmortem table ... "
A newline will be inserted at the end of the text string, so do not add one 
unless you want there to be a double-line-break.

    report
An open file descriptor to which report summary text should be printed.
"""
    if findex is None:
        test = fn(**kwarg)
    else:
        test = fn(**kwarg)[findex]
        
    test,truth = np.broadcast_arrays(test,truth)
    error = test - truth
    if percent:
        error /= truth
    I = np.nonzero(np.abs(error) > ep)[0]
    
    # If the test failed, do post-mortem
    if I.size:
        # We'll log the result shape to use for broadcasting
        shape = test.shape
        # Build an ordered list of the values to print in the table
        labels = ['test', 'truth', 'error']
        table = [test, truth, error]
        for index,value in kwarg.items():
            labels.append(index)
            table.append(np.broadcast_to(value, shape))
        # Now build the output
        # Print that the test failed
        report.write('[FAILED]    ' + text + '\n')
        # Print a table of the values that failed.
        for index,value in zip(labels, table):
            report.write('{: >15s}: '.format(index))
            for vv in value[I]:
                report.write('{: >15.6e}'.format(vv))
            report.write('\n')
        return False
        
    report.write('[passed]    ' + text + '\n')
    return True
