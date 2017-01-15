# PYro changelog

## Version 1.1: 
Original Release, including ideal gas data

## Version 1.2: 
- Add the 'steam' substance and its class, if97.
- Added new functionality to reg.__baseclass__._vectorize()
- Corrected inconsistent capitalization p/P usage across classes.  Lower case p 
    is now used for pressure in all classes.  This may cause reverse 
    compatibility issues between 1.1 and 1.2 if your code calls out P in its 
    arguments.
- Corrected a bug in psolve() caused by the same inconsistency
- Changed the "def_T" configuration parameter to 300. Some species' data are not
    defined below 300, so the default of 298.15 caused an error.
- Corrected all IG pressure references to be 1 bar instead of 1 atm to be 
    consistent with the source JANAF reference data.
- Replaced def_P with def_p to be consistent with the lower case pressure 
    definition.

## Version 1.3:
- Ported to GIT
- Changed package name from pyro to pyromat to avoid collision with existing 
    (https://pythonhosted.org/Pyro4/)
- Edited inline documentation to reflect changes
- Cleaned up some code in the data module (functionality should not change)
- Eliminated the input error class and reverted to the parameter error
- Changed the method for detecting the installation directory 
    utility.load_config()
- Eliminated a bug that caused the IGTAB class _lookup() method to fail in 
    Python3

## Version 1.4:
- Added the solve module
- Obsoleted the psolve() function
- Added inverse relations T_h() and T_s() to all classes
- Added inverse relations p_h() and p_s() to the ideal gas classes
- Changed the namespace handling so that registry class definitions behave like
    normal files.
- Added methods for calculating the inverse polynomials to the IF-97 class


