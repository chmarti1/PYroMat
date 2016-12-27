"""PYROMAT.DAT

Responsible for handling all data operations, the data
module supplies a dictionary where all loaded data are
housed, and a number of methods for data manipulation
and diagnostics.

Chris Martin (c) 2015,2017
"""


# load the root of the module
import pyromat as pyro
utility = pyro.utility
reg = pyro.reg






######################################
##
##  The master data dictionary
##
######################################
data = {}






#############################
##
##  Data manipulation methods
##
#############################
def load(datasource=None, check=None, verbose=None):
    """Import all *.hpd files in a directory
    load()
        or
    load('/path/to/files/')
        or
    load('/path/to/file.hpd')
        or
    info = load(check=True)

By default, load() will use the values in pyro.config{'dat_dir'} to 
find the files to load.  PYroMat automatically adds the 'data' directory
found in the installation directory.  Alternatively, if the first 
argument can be an explicit path to a file or a directory in which to 
find files.  All *.hpd files will be opened.

The pyro.config parameters that affect load() are:
'dat_verbose'
    Print messages to stdout? Override by setting
    the 'verbose' keyword argument.
'dat_exist_fatal'
    Raise an error if a data set already exists 
    in memory?  (default=False)
'dat_overwrite'
    Overwrite existing data? (default=True)
'dat_recursive'
    Recurse into subdirectories? (default=True)

The separate keyword argument, 'check' prompts
load() to run a data test instead of actually
loading data.

If load() is run in 'check' mode, it returns information about the data
files and data currently in memory.  The returned parameters are 
contained in a dictionary with the following keys:

changed
    A list of identifiers whose data in memory
    does not match the data currently in its 
    file.
added
    A list of identifiers appearing in memory but
    not appearing in the files.  These have been
    added since load.
removed
    A list of identifiers appearing in the 
    files but not in memory.  These have been
    removed since load.
redundant
    A list of identifiers for which there are
    multiple files.  These are redundant.
suppressed
    A list of files with a .hpd~ extension. 
    These files are ignored at load and may have
    been renamed by the utility.red_repair()
    function to correct a redundant definition.
data
    The loaded data set.
"""

    lead = 'load-> '
    root = False    # is this the root recursion in a check?

    # fetch the configuration parameters
    if verbose == None:
        verbose = utility.get_config('dat_verbose',dtype=bool)
    exist_fatal = utility.get_config('dat_exist_fatal',dtype=bool)
    exist_overwrite = utility.get_config('dat_overwrite',dtype=bool)
    recursive = utility.get_config('dat_recursive',dtype=bool)

    # If the load function is called with check=True, then it's time to 
    # make a few changes to the typical operation.  All recursive calls
    # will be called with check equal to the dictionary containing the
    # check results. 
    if check:
        if not isinstance(check, dict):
            check = {'changed':[], 'added':pyro.dat.data.keys(), 'removed':[], 'redundant':{}, 'suppressed':[], 'bad':[], 'data':{}}
            root = True
        loadto = check['data']
        CH = check['changed']
        ADD = check['added']
        REM = check['removed']
        RED = check['redundant']
        SUP = check['suppressed']
        BAD = check['bad']
    else:
        # if this is real, load the data into the hotpy data dictionary
        loadto = data


    # load is recursive. Unless load is called explicitly
    # with a path to a file to load, it will expand directories and
    # call itself with specific file names.  This arrangement prevents 
    # multiple copies of the same algorithm and saves an extra helper 
    # function.
    if datasource:
        datasource=utility.os.path.abspath(datasource)
        # if the data source is a directory
        if utility.os.path.isdir(datasource):
            # list the contents of the directory
            contents = utility.os.listdir(datasource)
            contents.sort()
            out='In directory ' + repr(datasource) + ' found files: '
            for this in contents:
                this_long = utility.os.path.join(datasource,this)
                # if recursion is enabled, and we come across a directory
                if recursive and utility.os.path.isdir(this_long):
                    #
                    # recurse into sub-directories
                    load(this_long,check=check,verbose=verbose)
                # if this is a file and it has the .hpd extension
                elif len(this)>4 and this[-4:]=='.hpd':
                    #
                    # recurse with the actual file name
                    load(this_long,check=check,verbose=verbose)

                    # assemble an output string
                    if verbose:
                        out += (this+', ')

                elif check and len(this)>5 and this[-5:]=='.hpd~':
                    #
                    # note if there are suppressed files
                    SUP.append(this_long)

            if verbose:
                utility.print_line('',lead)
                utility.print_line(out,lead)

        # if the data source is a file
        elif utility.os.path.isfile(datasource):


            # 
            # Load the data!
            # All recursive calls end here.
            #
            
            # If running in check mode, test to see if the load fails
            if check:
                try:
                    temp = utility.load_file(datasource)
                except:
                    BAD.append(datasource)
                    return
            else:
                try:
                    temp = utility.load_file(datasource)
                except:
                    return

            # Log the file source in the loaded data dictionary
            temp['fromfile'] = datasource

            # test for existance
            if temp['id'] in loadto:
                # if this identifier already exists in the loaded data

                if check:
                    # if run in check mode, we need to note the redundancy
                    if temp['id'] in RED:
                        # if the ID is already in the redundancy dictionary,
                        # add this file to the respective list
                        RED[temp['id']].append(datasource)
                    else:
                        # if this is the first time this file appears in the
                        # redundancy dictionary, construct a new list
                        RED[temp['id']] = [ loadto[temp['id']].data['fromfile'], datasource ]

                if exist_fatal:
                    # panic!
                    utility.print_error('Found an existing entry for ' + repr(temp['id']))
                    raise utility.HotPyDataError()

                elif not exist_overwrite:
                    # stop this recursion branch instead of overwriting the data
                    if verbose:
                        utility.print_warning('Found an existing entry for ' + repr(temp['id']) 
                            + '. Ignoring.')
                    return
                elif verbose:
                    # keep going, but warn the user that data is being overwritten
                    utility.print_warning('Found an existing entry for ' + repr(temp['id']) 
                        + '. Overwriting.')

            # only execute the rest of the data checks if the file isn't redundant
            elif check:
                # if run in check mode, it's time to do some extra work
                # if the identifier is in the add list, remove it
                if temp['id'] in ADD:
                    ADD.remove(temp['id'])
                # next, check to see if this identifier is in the hotpy data
                if not (temp['id'] in pyro.dat.data):
                    # this identifier seems to have been removed
                    REM.append(temp['id'])
                # finally check the data for identity
                elif data[temp['id']].data != temp:
                    CH.append(temp['id'])


            # look for the data class in the registry
            if temp['class'] in reg.registry:
                dataclass = reg.registry[temp['class']]
            else:
                utility.print_error('Species ' + repr(temp['id']) + ' called for data class ' + repr(temp['class']) + '.  That class does not exist in the registry.  The data file is corrupt or out of date.')
                raise utility.PMDataError()
            # write the data
            loadto[temp['id']] = dataclass(temp)

        else:
            # does not exist
            utility.print_warning('Data file or directory does not exist: ' + repr(datasource))

    # if called without an argument, use the config directories
    else:        

        for dd in pyro.config['dat_dir']:
            load(dd, check=check, verbose=verbose)


    if check and root:
        # if in check mode and this is the root of the recursion tree...
        if verbose:
            # if running verbosely, print a summary of the findings
            utility.print_line( '' , lead)
            # started with data that have changed
            utility.print_line( 'Found ' + str(len(check['changed'])) + ' elements that may have been CHANGED since load:', lead)
            out = ''
            for ss in check['changed']:
                out += ss + ', '
            utility.print_line( out[:-2], lead + '  ')

            # now data that have been added
            utility.print_line( 'Found ' + str(len(check['added'])) + ' elements that may have been ADDED since load:', lead)
            out = ''
            for ss in check['added']:
                out += ss
            utility.print_line( out, lead + '  ')

            # now data that have been added
            utility.print_line( 'Found ' + str(len(check['removed'])) + ' elements that may have been REMOVED since load:', lead)
            out = ''
            for ss in check['removed']:
                out += ss
            utility.print_line( out, lead + '  ')

            # now data with redundant definitions
            utility.print_line( 'Found ' + str(len(check['redundant'])) + ' files with REDUNDANT identifiers that are causing conflicts:', lead)
            for ss in check['redundant']:
                utility.print_line( ss+':', lead )
                for fil in check['redundant'][ss]:
                    utility.print_line( fil, lead + '  ')
            utility.print_line('',lead)

            # now suppressed files
            utility.print_line( 'Found ' + str(len(check['suppressed'])) + ' SUPPRESSED files in the search path:', lead)
            for fil in check['suppressed']:
                utility.print_line( fil, lead + '  ')
            utility.print_line('',lead)

            # now bad files
            utility.print_line( 'Found ' + str(len(check['bad'])) + ' BAD files in the search path:', lead)
            for fil in check['bad']:
                utility.print_line( fil, lead + '  ')
        return check











def clear():
    """Empty the data dictionary."""
    pyro.dat.data = {}










def new(newdata):
    """Create a new data entry
    new(newdata)

Creates a new entry in the HotPy data dictionary for a substance defined
by the 'newdata' dictionary.

The dictionary should have the data elements necessary for defining its
type.  The minimum elements are 'id', 'class', and 'doc', establishing 
the name, type, and description for the entry respectively.  
Additionally, the data dictionary should have whatever data is required 
by the class.

The new() function will load the data element into the data dictionary 
and add the 'fromfile' descriptor as if it had been created by the load()
function. 
"""
    lead = 'new-> '
    if ('id' in newdata) and ('class' in newdata):
        if newdata['class'] in reg.registry:
            data[newdata['id']] = reg.registry[newdata['class']]( newdata )
        else:
            utility.print_error(
'Could not find the class "' + newdata['class'] + '" in the registry.', lead)
            raise utility.HotPyDataError('Class not found')

    else:
        utility.print_error(
'The class data does not contain either the ''id'' or the ''class'' key.')
        raise utility.HotPyDataError('Bad data dictionary')











def updatefiles(dest=None, verbose=True, deletefiles=False):
    """Bring the files into agreement with data changes
    updatefiles()
        or
    updatefiles('/path/to/data/dir')

Changes to the data dictionary can be recorded 
permanently by running the updatefiles() function.
Changes are detected with the load(check=True)
function. Changed files are overwritten, added entries
are used to generate new .hpd files, and removed 
entries will have their corresponding files deleted.
Any entries that show redundancy errors will result
in a prompt for how to repair the conflict.  
Files in conflicts that are not resolved will be 
disqualified from other operations.

Existing files that qualify for an update will be 
replaced based on their entry in the 'fromfile' 
keyword.  Data entries with no record of their parent 
file or that were created from scratch will be placed 
in the first directory in the constants.DATADIR list
unless an alternate location is explicitly specified 
in the save() function call.

It is strongly recommended that the operation be 
performed verbosely, where the user is prompted for 
permission to write each file.  If the application
demands it, updatefiles() can be run quietly
    update(verbose=False)
In this mode, it will use the flag found in the
'deletefiles' keyword to decide whether to delete
or suppress all files without prompting the user
for a decision.
"""
    
    lead = 'updatefiles-> '
    
    if dest==None:
        dest = pyro.config['dat_dir'][0]
    
    status = load(check=True,verbose=False)

    # check for redundancy problems
    if verbose and status['redundant']:
        pyro.utility.print_line('Resolving redundancy issues in the files.',lead)
        pyro.utility.print_line(' ',lead)
    if status['redundant']:
        # First, repair redundancies
        pyro.utility.red_repair(status['redundant'], verbose=verbose)



    # check for bad files
    if verbose and status['bad']:
        pyro.utility.print_line(
'Resolving issues with files that failed to load.  You will be prompted to delete, suppress, or ignore each of them.  Ignoring a file is safest, but will not correct the problem.  Selecting "delete" will permanently remove the file.  Suppression adds a "~" to the file name so load() will ignore it.',lead)
        pyro.utility.print_line(' ',lead)
        for ss in status['bad']:
            # print the file name under consideration
            # and prompt the user for a choice of what to do
            pyro.utility.sys.stdout.write(ss)
            charin = 'z'
            while not charin in 'dsi':
                charin = raw_input('(D)elete/(S)uppress/(I)gnore:').lower()
            if charin == 'd':
                try:
                    pyro.utility.os.remove(ss)
                    pyro.utility.print_line('Removed.',lead)
                except:
                    pyro.utility.print_warning(
'Deletion failed. Ignoring. Check permissions.')

            elif charin == 's':
                try:
                    pyro.utility.suppress_file(ss)
                except:
                    pass

            elif charin == 'i':
                pyro.utility.print_line(
'Ignoring.  This problem will persist until corrected.',lead)

    elif status['bad']:
        for ss in status['bad']:
            if deletefiles:
                try:
                    pyro.utility.os.remove(ss)
                except:
                    pyro.utility.print_warning(
'Failed to remove file ' + ss + '.  Ignoring.  Check permissions.')
            else:
                try:
                    pyro.utility.suppress_files(ss)
                except:
                    pass


    if status['bad'] or status['redundant']:
        # Re run the check
        status = load(check=True, verbose=False)


    # next, save any new data
    if verbose and status['added']:
        pyro.utility.print_line('Saving data that were added since load.',lead)
        pyro.utility.print_line(' ',lead)
    for ss in status['added']:
        this = pyro.dat.data[ss]
        if ('fromfile' in this.data) and this.data['fromfile']:
            fil = this.data['fromfile']
        else:
            fil = dest
            if dest[-1]!=pyro.utility.os.path.sep:
                fil += pyro.utility.os.path.sep
            fil +=  this.data['id'] + '.hpd'
        fil = pyro.utility.os.path.abspath(fil)

        FIL = None
        try:
            FIL = open(fil,'w')
            pyro.utility.json.dump(this.data,FIL,sort_keys=True,indent=4)
            FIL.close()
        except:
            pyro.utility.print_warning(
'Failed to create file: ' + fil + '.  Ignoring.  Check permissions and re-run to correct.')
            if FIL:
                FIL.close()
        if verbose:
            pyro.utility.print_line('Wrote new file: ' + fil, lead)
            pyro.utility.print_line('',lead)



            
    # now save changes
    if verbose and status['changed']:
        pyro.utility.print_line('Updating files that were changed since load.',lead)
        pyro.utility.print_line(' ',lead)
    for ss in status['changed']:
        this = pyro.dat.data[ss]
        # get the file target from the loaded file
        fil = status['data'][ss].data['fromfile']
        # get the fromfile field from the data in memory
        dfil = '* none *'
        if 'fromfile' in this.data:
            dfil = this.data['fromfile']
        # If the fromfile fields don't agree, it could be a sign that these
        # aren't really from the same file.  Prompt the user for an option.
        if (fil != dfil) and verbose:
            pyro.utility.print_line(
'The data loaded for ' + ss + 
' seems to have originated from a source other than the file found by load().',lead)
            pyro.utility.print_line('Load() found: ' + fil, lead)
            pyro.utility.print_line('In memory:    ' + dfil, lead)
            
            charin = 'z'
            while not charin in 'yn':
                charin = raw_input('Overwrite? (y/n):').lower()

        try:
            FIL = open(fil,'w')
            pyro.utility.json.dump(this.data,FIL,sort_keys=True,indent=4)
            FIL.close()
            if verbose:
                pyro.utility.print_line('Updated file: ' + fil, lead)
                pyro.utility.print_line('',lead)
        except:
            pyro.utility.print_warning(
'Failed to update file: ' + fil + 
'.  Ignoring.  Check permissions and re-run to correct.')


    
    # remove files
    if verbose and status['removed']:
        pyro.utility.print_line(
'Correcting files that were not represented in the data.',lead)
        pyro.utility.print_line('',lead)
    # deletion or suppression?
    if deletefiles:
        for ss in status['removed']:
            fil = status['data'][ss].data['fromfile']
            charin='y'
            if verbose:
                pyro.utility.sys.stdout.write(fil + '\n')
                charin = 'z'
                while not charin in 'yn':
                    charin = raw_input('Delete? (y/n):').lower()
            if charin=='y':
                try:
                    pyro.utility.os.remove(fil)
                except:
                    pyro.utility.print_warning(
'Failed to delete file: ' + fil + '.  Ignoring.  Check permissions and re-run to correct.')
    else:
        for ss in status['removed']:
            fil = status['data'][ss].data['fromfile']
            try:
                utility.suppress_file(fil,verbose=verbose)
            except:
                pyro.utility.print_warning(
'Failed to suppress file: ' + fil + '.  Ignoring.  Check permissions and re-run to correct.')



#
# Run the load function
#
load()


