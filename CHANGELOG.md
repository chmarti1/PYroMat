# PYroMat changelog

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

## Version 2.0.1
This is the first version that deliberately breaks reverse compatibility.  Every time reverse compatibility is not preserved, the major version number will increase.

- Provides a uniform array handling behavior for all classes
    - All methods now accept arrays or array-like objects of any dimension
    - All methods attempt reduce their returns to scalars if possible
    - All methods return values that obey Numpy broadcasting rules
    - Some methods will return a Python float and others return a Numpy array scalar depending on implementation
- Provides a new configuration class, for in-line configuration
    - config[] is now a dictionary-like object and NOT a dictionary.
    - config[] enforces its own rules about type etc.
    - config scripts still behave identically.
- Provides a new unit conversion module, "units"
    - config["unit_XXX"] now configures global unit defaults for unit XXX.
    - Added unit conversion routines
    - in-line documentation for all property methods reports the config unit parameters on which they depend.
- Modifies the info() function to print a table of supported properties
- Adds Tlim() and plim() to IF-97 (steam)
- Migrates specific heat ratio method k() to gam() in preparation for k() to become thermal conductivity

## Version 2.0.4
The version increments between 2.0.1 and 2.0.4 were primarily spent correcting issues with the python package index and documentation.
- Corrected a bug that prevented `info()` from displaying in Python 3.

## Version 2.0.5
- Corrected a bug preventing T_s() in the if97 (steam) class from returning quality correctly.

## Version 2.0.6
- Corrected a bug that caused hsd() in the if97 (steam) class to crash
- Corrected a bug in the configuration object that prevented a parameter summary from displaying correctly in Python 3.

## Version 2.0.7
- Added the mp1 multi-phase class
- Migrated steam away from if97 to the mp1 class
- Added CO2 and R134a (C2H2F4) to the multiphase data
- Added ig2 class incorporating the NASA polynomial
- Added hundreds of new ideal gas data including a better ig.H2O data set

## Version 2.0.8
- Corrected unit conversion in mp1 density output

## Version 2.0.9
- Skipped due to an error uploading to the Python package index

## Version 2.0.10
- Temperature units in T_s, T_h, and T in the mp1 class were not being converted: fixed.

## Version 2.0.11
- Corrected a bug in units.matter(); inplace directives were not being honored in certain cases, causing odd errors in T_s and T_h()

## Version 2.0.12
- Corrected data error in CO2 that caused errors near the saturation curve
- Added support for specifying T,p, and x simultaneously; p is ignored unless x is less than 0
- Added the aps module with support for calculating the performance of thermodynamic cycles

## Version 2.0.13
- Transitioned to the `_hybrid1()` inversion algorithm for mp1's `T_s` and `T_h` algorithms

## Version 2.1.0
- Completed total transition to and validation of `_hybrid1()` in the `mp1` class  
- Added `d_s()` to the `mp1` class  
- Standardized `igmix`, `ig`, and `ig2` classes to use the same argument parsing rules  
- Standardized `igmix`, `ig`, and `ig2` inverse property routines `T_s`, `T_h`, `p_s`.  
- Added filtering capabilities to `info()`  
- Added `atoms()` to all classes to retrieve atomic composition data.  
- Updated the `dat.updatefiles()` function to correct bugs in Python 3.  
- Changed `PMConfigEntry.apply_default()` to `restore_default()`  
- Changed `PMConfig.write()` algorithm to deal more gracefully with appendable entries.  
- Added function type 0 to `mp1._satfit` for future melting line support.  
- Added `mp.N2` and `mp.CH4` multiphase models.  
- Added `R()` and `mw()` to the `mp1` class.  
- Adjusted all multi-phase R-values to agree with Ru / mw to numerical (double) precision.  For most properties, the change will not be noticeable, but it brings some into better agreement with published values.  In all cases, the change is on the order .01% or smaller.

## Version 2.1.1
- Corrected a bug in `setup.py` that prevented the ideal gas mixture data from being installed correctly.  

## Version 2.1.2
- Corrected a bug in `mp1.T_s` and `mp1.T_h` that caused a crash when calculating quality of arrays in isobaric mode.

## Version 2.1.3
- Corrected a bug in `ig2.s` that caused a crash with pressure arrays

## Version 2.1.4
- Corrected keystroke errors in `igmix` methods
- Corrected a bug in `igmix` and `ig2` that caused incorrect values when working with arrays in `cp` and `gam`

## Version 2.1.5
- Corrected a bug in `mp1.d_s()` that caused crashes above the critical point

## Version 2.1.6
- Moved `mp1.T()` to use `_tditer()` to address bracketing errors at super-critical states.

## Version 2.1.7
- Suppressed accidental verbose output in `mp1._T()`

## Version 2.1.8
- Improved convergence performance of `mp1.d()` very close to the phase transition.  Opened the bracketing interval by 1% to compensate for small numerical errors.

## Version 2.1.9
- Corrected an array handling bug that caused intermittent crashes in `mp1.d_s()`.

## Version 2.1.10
- Corrected a bug in `igmix.T()` that caused it to return pressure
- Corrected a bug: `igmix._argparse()` tried to use `data['mw']` for density conversion.

## Version 2.2.0
- Corrected small errors in mp1 model properties near the critical point.
  This bug resulted in strange things like negative cp values, but it was only found VERY close to the critical point.  Otherwise, errors were so small that the properties still passed validation checks against reference data.
- Added iteration to the pyromat configuration class
- Added the specific volume property to all classes
- Added the warning_verbose configuration entry to allow users to mute warning messages
- Transitioned to a fully flexible property argument format that accepts combinations of h,s,e,T,p,d,v,x.
- Obsoleted "inverse" property methods (like T_s and T_h) - they are still available for reverse compatibility.
- Changed the out-of-bounds detection behavior: now returns config["def_oob"] (default np.nan) on those array elements that are out-of-bounds.
- Added O2, R1234ze to the multiphase collection
- Corrected a number of minor bugs reported since the last release.  These are documented on the PYroMat github issues page.

## Version 2.2.1
- Issued bugfixes for github issues 41, 42, 43
  - Corrected a boneheaded typo in igmix that should have been caught in testing.
  - Reverted to `_tditer()` in `mp1._T()`
  - Eliminated the upper 1% grace range in `Ts()` to prevent imaginary values past Tc.

## Version 2.2.2
- Changed the call signature for `T_s`, `T_h` functions to address issue 52
  *NOTE* inverse methods like `T_s` and `T_h` are now deprecated.
- Reassessed all multiphase `dlim` values in the core data. This addresses issues 44, 45, and 46.

## Version 2.2.3
- Updated the README to adopt recommendations made by the JOSS community - specifically to include recommendations for community involvement.
- Added the optional `pip install pyromat[dev]` option, which requires the `pytest` package.

## Version 2.2.4
- Corrected a bug reported in issue 64 where inverse routines were not returning the correct units.

## Version 2.2.5
- Modified `ig`, `ig2`, and `igmix` classes to accept enthalpy and entropy simultaneously (github issue 83)
- Added entropy of mixing to the `igmix` s() calculation (github issue 92)
- Added the `igtools` module with dynamic ideal gas mixture support
- Added the `ismass()` function in the units module
- Added `sid()` to the `__basedata__` class (all data instances)
- Added `hill()` to the `__basedata__` class
- Added `def_p_unit` and `def_T_unit` to the configuration system (github issue 68)
- Added `def_p()` and `def_T()` as methods to the `PMConfig` class to automatically handle the default units
- Added `a()` speed of sound method to all classes
- Corrected the error in the R1234ze saturation line data (github issue 86)
- Added Helmholtz (`f`) and Gibbs (`g`) energies to all classes
- Wrote the _ds() precision inner saturation property method; fully functional, but not yet used.

## Version 2.2.6
- Changed `mp1.state()` to adapt to small errors in the saturation properties (github issue 104).
- Corrected crashes in some arithmetic combinations with `igtools.IgtMix` (github issue 103).
- Added unhandled cases to `ig.f()`, `ig.g()`, `ig2.f()`, and `ig2.g()` when `_argparse` does not return `p` (github issue 100).
- Forced `IgtMix.T()`, `IgtMix.p()`, and `IgtMix.d()` to convert views into copies (github issue 101).
- Created `IgtMix._mw()` to handle internal calcualtions without unit conversions (github issue 108).
- Added scalar division to the `IgtMix` class.
- Corrected bugs and made minor improvements to `__basedata__.hill()` suggested by @allrob23.
- Corrected bugs in `_mp1._sat_argparse()`, `Ts()`, and `ps()` that ignored the default state units (github issue 99).
- Corrected a bug in `IgtMix.cv()` that gave incorrect values (github issue 102).
