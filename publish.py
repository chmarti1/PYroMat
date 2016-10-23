# Not part of the PYro distribution
# publish.py moves the contents of the directory to an appropriately named
# directory and compresses it into a .zip and .tar.bz2 distribution.
# ***
# This file is intended to be run from the pyro directory.
# ***
#
# This is a record of the SVN revisions corresponding to each public release
# r31 v1.1  First public release
# r64 v1.2 Added steam, unified to lower-case p

import shutil
import os
import time

# We'll figure out what the working version is from the __version__ member
# string in the code itself.  In what file should we look for the __version__ definition?
vfile = '__init__.py'
# The variable name containing the version string in install_init
vlookfor = '__version__'

# The distribution will be a directory named PYroMat-[version].  The root 
# directory contents are determined by the "toroot" list.
# Inside the root directory will be an installation directory.  Its contents
# are determined by the "toinstall" list.
install_dir = 'lib'
# Inventory of things to send to the distribution root
toroot = [
    'setup.py',     # Python distribution/install file
    'README.md',    # Basic documentation
    'CHANGELOG.md',    # Well, it's a changelog, silly
]
# Inventory of things to send to the lib directory
tolib = [
    '__init__.py',  # The root module
    'dat.py',       # Data module
    'reg.py',       # Registry module
    'utility.py',   # Utility module
    'defaults.py',  # Default configuration file
    'registry',     # Registry directory
    'data',         # Data directory
]

# for debugging - how much are we doing here?
# 0 - dry run; generate a report and halt
# 1 - dry run; report files to copy and halt
# 2 - clean and copy; clean old builds and make the new one
# 3 - make the zip and tarball
dolevel = 3


# Start a log file
with open('publish.log','w+') as logfile:
    # Make a standard output format
    def logline(addline):
        output = '[' + time.strftime('%H:%M:%S %m-%d-%Y') + '] '+ addline + '\n'
        sys.stdout.write( output )
        logfile.write( output )

    logline('Starting a build with DOLEVEL == {:d}'.format(dolevel))
    # Within the destination directory, there will be a directory containing
    # the files to be installed on the users' systems.  What should it be called?

 
    # scans __init__.py for the __version__ string
    NS = {vlookfor:'0.0'}   # Create a dummy locals dictionary with __version__
    with open(vfile,'r') as init_file:
        for line in init_file:
            # If the line starts with __version__
            if line[:len(vlookfor)] == vlookfor:
                # Execute it - it will assign the version string
                exec(line,globals(),NS)
                break
    # copy the result        
    version = str(NS[vlookfor])
    logline( "Working version: " + version )

    # Where is the destination directory?
    base_name = 'pyro-' + version
    dest_dir = os.path.abspath(base_name)
    logline( "Build root: " + dest_dir )
    install_dir = os.path.join(dest_dir,install_dir)
    logline( "Library: " + install_dir )
    linux_file = os.path.abspath(base_name + '.tar.bz2')
    logline( "Linux output: " + linux_file )
    win_file = os.path.abspath( base_name + '.zip' )
    logline( "Windows output: " + win_file)
    
    if dolevel<1:
        logline( "Exited cleanly with DOLEVEL {:d}".format(dolevel) )
        raise Exception()
    ###
    # End DOLEVEL 0
    ###
    try:
        # Test the files and directories
        logline("Testing TOROOT files")
        rooterrors = 0
        for this in toroot:
            if os.path.lexists(this):
                if os.path.isfile(this):
                    logline( "   File: " + os.path.abspath(this))
                elif os.path.isdir(this):
                    logline( "   Directory: " + os.path.abspath(this))
                else:
                    logline( " ? Unknown: " + os.path.abspath(this))
                    rooterrors+=1
            else:
                logline( " X NOT FOUND: " + os.path.abspath(this))
                rooterrors += 1

        logline("Testing TOLIB files")
        liberrors = 0
        for this in tolib:
            if os.path.lexists(this):
                if os.path.isfile(this):
                    logline( "   File: " + os.path.abspath(this))
                elif os.path.isdir(this):
                    logline( "   Directory: " + os.path.abspath(this))
                else:
                    logline( " ? Unknown: " + os.path.abspath(this))
                    liberrors+=1
            else:
                logline( " X NOT FOUND: " + os.path.abspath(this))
                liberrors += 1

        logline("Found {:d} root errors and {:d} library errors".format(
                rooterrors,liberrors))
    except:
        logline("UNHANDLED EXCEPTION in DOLEVEL {:d}".format(dolevel))
        output = sys.exc_info()
        logline(output[1])
        raise output[0](output[1])
    
    if rooterrors or liberrors:
        logline("FATAL ERROR.  Halting.")
        raise Exception()
    elif dolevel<2:
        logline( "Exited cleanly with DOLEVEL {:d}".format(dolevel) )
        raise Exception()
    ###
    # End DOLEVEL 1
    ###
    try:
        # If the directory exists, remove it
        if shutil.os.path.isdir(dest_dir):
            logline("Removing existing destination " + dest_dir)
            shutil.rmtree(dest_dir)

        # Make the parallel directory
        logline("Creating " + dest_dir)
        os.mkdir(dest_dir)
        # make the install dir
        logline("Creating " + install_dir)
        os.mkdir(install_dir)
            
        # Now COPY!
        for this in toroot:
            dest = shutil.os.path.join(dest_dir,this)
            if os.path.isdir(this):
                logline( "Copying directory: " + this + " >> " + dest )
                shutil.copytree(this,dest)
            else:
                logline( "Copying file: " + this + " >> " + dest )
                shutil.copy(this, dest)

        for this in tolib:
            dest = shutil.os.path.join(install_dir,this)
            if os.path.isdir(this):
                logline( "Copying directory: " + this + " >> " + dest )
                shutil.copytree(this,dest)
            else:
                logline( "Copying file: " + this + " >> " + dest )
                shutil.copy(this, dest)

    except:
        logline("UNHANDLED EXCEPTION in DOLEVEL {:d}".format(dolevel))
        output = sys.exc_info()
        logline(output[1])
        raise output[0](output[1])
    
    if dolevel<3:
        logline( "Exited cleanly with DOLEVEL {:d}".format(dolevel) )
        raise Exception()

    ###
    # End DOLEVEL 2
    ###   

    try:
        # Now compress
        if os.path.isfile(linux_file):
            logline( "Removing old Linux install: " + linux_file )
            os.remove(linux_file)
        logline( "Generating Linux install: " + linux_file )
        os.system('tar -cvjf {0} {1}'.format(linux_file,base_name))

        if os.path.isfile(win_file):
            logline( "Removing old Windows install: " + win_file )
            os.remove(win_file)
        logline( "Generating Windows install: " + win_file )
        os.system('zip -r {0} {1}'.format(win_file,base_name))

    except:
        logline("UNHANDLED EXCEPTION in DOLEVEL {:d}".format(dolevel))
        output = sys.exc_info()
        logline(output[1])
        raise output[0](output[1])
    

    logline( "Done." )
