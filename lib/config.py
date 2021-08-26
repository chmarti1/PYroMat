#
# PYroMat configuration file
#   Chris Martin (c) 2015,2017,2021
#
# When PYroMat loads config files, it automatically looks for 'defaults.py' 
# in the installation directory.  Here, administrators can configure 
# advanced options like allowing users to keep their own config files 
# that override the defaults or maintaining their own registries and data
# in their home directories.
#
# Configuration files are ordinary python scripts that are executed when 
# PYroMat decides to load its configuration options.  PYroMat inspects the 
# locals() dictionary to see what (if any) variables were created.  If 
# the parameters are recognized, they are loaded into the pyromat.config 
# dictionary.  This means that script files are just ordinary python code.
#
# If PYroMat finds variable names that aren't recognized configuration 
# parameters, it throws an error.  That means that while the configuration
# files are ordinary code, there are strict rules for what variables can be
# defined and what their values can be.
#
# There are three types of parameters:
# 1 - Normal parameters require a specific data type, but each time they are
#     written, the prior value is lost.
# e.g.
#>  dat_verbose = True  # will overwrite the default of False
#
# 2 - Appended parameters never lose data, but contain a list of all prior 
#     values.  If values are assigned to these parameters in multiple files, 
#     they are simply appended to the list.  Within a single configuration file, 
#     the normal variable assignment rules apply, so writing to a variable twice 
#     will destroy the original value.  To write multiple values at once, assign
#     a list to the parameter
# e.g.
#>  # This will add two new configuration files to the load sequence
#>  config_file = ["/etc/pyromat/config.py", "/home/$USER/.pyromat/config.py"]
#
# 3 - Read-only parameters supply a value for information purposes only.
# e.g.
#>  install_dir = "/home/malicious_user/malicious_project/"   # Won't work
#>  version = "596.0.1" # Will throw an error - with good reason



#** Configuration files **
# This tells PYroMat where to find other configuration files. Keep in mind 
# that the last file loaded will be given precedence.  PYroMat resolves Windows
# and Unix references to home and environment variables, so references to ~, 
# $USER, etc. have the desired result.
#
#> config_file = '/path/to/my/file'
#
#   or, perhaps
#
#> config_file = '/etc/pyromat.py'
#
#   also popular is
#
#> config_file = '~/.pyromat/config.py'
#
#   but you can always get explicit
#
#> config_file = ['/path/to/config1', '/path/to/config2']


# Should load_config() print its activity to stdout?  If you are setting up a
# system where users are allowed their own configuration options in their own
# home directories, this might be helpful for debugging.
#
# Configuration options only take effect AFTER they are loaded, so this file's
# activity will not be affected by this directive.
config_verbose = False


#** Data directories **
# The data directory list contains "install_dir/data" by default.  Users
# can add data from their own directories as well.  It is not possible
# to overwrite previous entries - you can only add additional 
# directories.
#
#> dat_dir = '~/.pyromat/data'
#
#   or
#
#> dat_dir = ['/usr/share/pyromat/data', '/home/bestuserever/.pyromat/data']

# Should the PYroMat dat module print its load activity to stdout?
#
dat_verbose = False

# If two data files are found with the same 'id', what do we do; overwrite or 
# ignore?  This determines the order of precedence.  In overwrite mode, data 
# files entries in data directories specified later in the dat_dir list (like 
# a user's home directory) will overwrite those higher up in the tree.  This 
# allows users to experiment with their own data classes that replace the built-
# in files.  If you want the files distributed with PYroMat to take precedence, 
# you should change this to False.
dat_overwrite = True

# If two data files are found with the same name, should we exit with an
# error?  If you are really want to put a stop to users overriding the existing
# data files, set this directive to True.
dat_exist_fatal = False

# Should we descend into sub-directories of the directories listed in
# the dat_dir parameter?  Unfortunately, there is not currently a way to apply 
# this setting locally.  PYroMat's recursion is an all-or-none.
dat_recursive = True


#** Registry behavior **
# By default, the registry will consist of class definitions found in 
# "install_dir/registry".  However, users may want to define their own data
# classes from their own registry directories.  It is not possible to 
# overwrite previous entries - you can only add additional registry
# directories.
#
# System administrators should be VERY careful of the latter 
# implementation.  bestuserever really needs to be an extra special
# user, because she will be able to write scripts that will be executed
# as whatever user loads PYroMat.
#
#> reg_dir = '~/.pyromat/reg'
#
#   or
#
#> reg_dir = ['/usr/share/pyromat/reg', '/home/bestuserever/.pyromat/reg']

# Should the pyromat.reg.regload() function print its activity to stdout?
reg_verbose = False

# If redundant class definitions are discovered, should we overwrite the
# old one, or should we ignore the new one?  This serves the same function as
# the dat_overwrite directive.
reg_overwrite = True

# Should a redundant class definition cause PYroMat to exit with an error?  This
# serves the same function as the dat_exist_fatal directive
reg_exist_fatal = False

# What is the default temperature and pressure that property functions
# should use when entries are omitted?  These must be in the same units
# specified by unit_pressure and unit_temperature
#> def_T = 298.15
#> def_p = 1.01325


# In what units should functions accept and return values?  A list of 
# the accepted units appears below each directive
#> unit_force = 'N'
#> unit_energy = 'J'
#> unit_temperature = 'K'
#> unit_pressure = 'bar'
#> unit_molar = 'kmol'
#> unit_volume = 'm3'
#> unit_length = 'm'
#> unit_mass = 'kg'
#> unit_time = 's'
#> unit_matter = 'kg'
#
# Matter can be either molar OR mass units.  This specifies the unit to 
# be used by intensive properties.
#
# For a list of available units, see the pyromat.units documentation or
# type pyromat.units.show()

