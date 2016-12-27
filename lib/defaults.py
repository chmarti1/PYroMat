#
# PYroMat defaults configuration file
#   Chris Martin (c) 2015,2017
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
# PYroMat uses an unusual configuration loading scheme where no configuration 
# parameter is ever truly overwritten.  Instead, values are appended to a
# list of all values ever assigned to the parameter.  In this way, 
# parameters that are intended to be configurable for multiple values 
# need no special treatment, and it is easy to see a parameter's history.
#
# The line "myparam = 'value'" in a config file will result in
# >>> import pyromat as pyro
# >>> pyro.config{'myparam'}
# ['initial_value', 'value']
#
# where "initial_value" indicates the hard-coded default given to the 
# parameter when the pyro.utility.load_config() function is called.
#
# See the pyro.utility.get_config() function for conveniently accessing
# the configuration parameters.  Thanks to the get_config() function, a
# well-intentioned user who does this:
#
# >>> pyro.config{'myparam'} = 'value'
#
# will not cause any problems.  If load_config() is called again, this
# parameter will now be treated as read-only, and config files will not 
# be allowed to modify to myparam.



#** Configuration files **
# This tells PYroMat where to find other configuration files. Keep in mind 
# that the last file loaded will be given precedence.
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
#   buy you can always get explicit
#
#> config_file = ['/path/to/config1', '/path/to/config2']

# Should load_config() print its activity to stdout?  
#
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

# If two data files are found with the same 'id', what do we do; 
# overwrite or ignore?
dat_overwrite = True

# If two data files are found with the same name, should we exit with an
# error?
dat_exist_fatal = False

# Should we descend into sub-directories of the directories listed in
# the dat_dir parameter?
dat_recursive = True


#** Registry behavior **
# By default, the registry will consist of class definitions found in 
# "install_dir/reg".  However, users may want to define their own data
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
# old one, or should we ignore the new one?
reg_overwrite = True

# Should a redundant class definition cause PYroMat to exit with an error?
reg_exist_fatal = False

# What is the default temperature and pressure that property functions
# should use when entries are omitted?
#> def_T = 300.
#> def_p = 1.01325
