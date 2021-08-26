"""Unit conversion module for PYroMat

To print a summary of all conversions supplied, call
>>> units.show()
          force : lb lbf kN N oz kgf 
         energy : BTU kJ J cal eV kcal BTU_ISO 
    temperature : K R eV C F 
       pressure : mmHg psi inHg MPa inH2O kPa Pa bar atm GPa torr mmH2O ksi 
          molar : Ncum NL Nm3 kmol scf n mol sci Ncc lbmol 
         volume : cumm cc mL L mm3 in3 gal UKgal cuin ft3 cuft USgal m3 cum 
         length : ft nm cm mm m km um mile in 
           mass : mg kg g oz lb lbm slug 
           time : s ms min hr ns year day us 

To obtain a list for the units recognized by a particular conversion
class, call the get() method,
>>> units.length.get()
['ft', 'nm', 'cm', 'mm', 'm', 'km', 'um', 'mile', 'in']

There are also special functions defined,
matter()                Converts between molar and mass
temperature_scale()     Converts between the temperature scales
gauge_to_abs()          Converts pressures between absolute and gauge
abs_to_gauge()

To alter the standards by which units are defined, call the setup()
function.
>>> setup(Tstd = 300.)

The setup() function re-defines the various constants used to calculate
these unit conversions.  Many of them are physical constants, but some
have arbitrary definitions.  The constants that were used to formulate
the conversions are exposed as module members beginning with 'const_'. 

Constants that are calculated in terms of other more fundamental 
constants are shown with their formulae in parentheses.  Constants that
can be adjusted by a call to the setup() function are shown with their
default values in parentheses.

const_k     1.38064852e-23      J/K             Boltzmann's constant
const_Na    6.022140857e23      count/mol       Avagardo's number
const_q     1.602176634e-19     coulomb/e       Elementary charge
const_Nc    6.24150934e18       count/coulomb   Definition of a Coulomb (1/e)
const_Ru    (k*Na)              J/mol/K         Ideal gas constant
const_Tstd  (273.15)            K               Standard temperature
const_pstd  (1.01325)           bar             Standard pressure
const_dstd  (pstd*1e5/Ru/Tstd)  mol/m3          Standard concentration
const_g     (9.80665)           m/s2            Accel. d.t. gravity
const_dh2o  (999.9720)          kg/m3           Water column density
const_dhg   (13595.1)           kg/m3           Mercury column density
"""

import numpy as np
import pyromat as pyro
import sys


class Conversion:
    """CONVERSION CLASS
The unit conversion class simulates a function that converts a value of 
one set of units to another.  

  new_value = conversion_object( old_value, from_units, to_units )

for example

  value_ft = length(6., 'in', 'ft)

When unit specifiers are omitted, the conversion object will throw an 
error unless conversion_object.config_default is a string.  If so, it
specifies a parameter that users can set in their PYroMat configuration
files.

There is an optional 'exponent' parameter.  This is to allow support for
cases like

  volume_ft = length(volume_in, 'in', 'ft', exponent=3)

or

  stress_psi = length(stress_psf, 'ft', 'in', exponent=-2)
"""



    def __init__(self, table, config_default=None):
        """Conversion(table, config_default)
Conversion objects are initialized with a conversion table in the format
of a dictionary.  The dictionary keys are the strings that identify the
various units, and the values are used to perform the conversion.

>>> inft = Conversion({'in':12., 'ft':1.})
>>> inft(6., from_units='in', to_units='ft')
.5

Note that the values in the table are chosen so that
>>> new_value = old_value * table[to_units] / table[from_units]

The __call__ method has two additional optional parameters; exponent
and inplace.  In the "inft" example above, the "exponent" parameter 
might be set to 3 to make the conversion volumetric instead of length.
The inplace parameter is intended for operating on numpy arrays 
"in-place" so that instead of returning a new array, the original value
array will have its values converted "in-place".

There is an optional config_default parameter that can be used to specify
the behavior of the 'to_units' and 'from_units' strings if they are not
specified when the conversion is evoked.  The config_default does NOT
define the default unit string; rather it specifies the parameter name
to query using pyro.utility.get_config() for the default units.  This
allows users to completely reconfigure the PYroMat unit system on the
fly.
"""
        self.table = table
        self.config_default = config_default
        

    def __contains__(self, unit):
        """Test whether a particular unit string is supported"""
        return (unit in self.table)

    def __call__(self, value=1., from_units=None, to_units=None, exponent=None, inplace=False):
        """Execultes a conversion from [from_units]**[exponent] to 
[to_units]**[exponent].  By default, [value] is 1., so that the value 
returned is the appropriate conversion factor to convert any set of 
units.  Otherwise, [value] should be the quantity to be converted to the
new system of units.
"""
        # Check to be certain from_units is specified
        if from_units is None:
            # If not, retrieve the configured default units
            if self.config_default is None:
                raise pyro.utility.PMParamError('Missing from_units, and no default specified')
            else:
                from_units = pyro.config[self.config_default]
        # Repeat for to_units.  Are they specified?
        if to_units is None:
            # If not, retrieve the configured default units
            if self.config_default is None:
                raise pyro.utility.PMParamError('Missing from_units, and no default specified')
            else:
                to_units = pyro.config[self.config_default]
        # Do not do the conversion if it is not necessary
        if from_units == to_units:
            return value

        conv = self.table[from_units] / self.table[to_units]
        if exponent:
            conv **= exponent

        if inplace:
            # Point to the original array if possible
            return np.multiply(value, conv,out=value)
        
        return np.multiply(value, conv)
        


    def __getitem__(self,item):
        return self.table.__getitem__(item)

    def __setitem__(self,item,value):
        return self.table.__setitem__(item,value)


    def get(self):
        """Return an unorded list of the units supported"""
        return self.table.keys()




def setup(Tstd=273.15, pstd=1.01325, g=9.80665, dh2o=999.9720, dhg=13595.1):
    """Set up the conversion functions
    setup(Tstd=273.15, pstd=1.01325, g=9.80665, dh2o=999.9720, dhg=13595.1)

The arguments to setup() allow users to re-define certain conventional 
constants from which unit conversions are derived.  For example, the 
standard temperature and pressure have no universal convention.  
Similarly, the density of mercury and water used for column pressure can
obey contradictory conventions, and the acceleration due to gravity might
be adjusted for operation at different points on the globe.

The parameters and their defaults are defined as follows:
Tstd = 273.15 Kelvin (0 degC)
pstd = 1.01325 Bar (The US Standard atmosphere)
g = 9.80665 m/s^2 (Acceleration due to gravity at 45 deg latitude)
dh2o = 999.9720 kg/m^3 (Density of water at 4 degC)
dhg = 13595.1 kg/m^3 (Density of mercury at 0 degC)

Changes to these values must remain in the units shown above.

This function writes constants,
# Plank constant
const_h = 6.62607015e-34

#Boltzman
# Used to convert between eV and temperature
const_k = 1.38064852e-23 # J/K

#Avagadro's number
# Used by the standard volume and molar units
const_Na = 6.022140857e23 # count/mol

#Coulomb
# Used to convert between eV and Joules
const_q = 1.602176634e-19 # coulomb
const_Nc = 1/const_q # count/coulomb

#Universal gas constant
# Used by the standard volume molar units
const_Ru = const_k * const_Na # J/mol/K

#Standard Temperature
# Used by the standard volume molar units
const_Tstd = Tstd #K

#Standard Pressure (US Std atmosphere)
# Used by the standard volume molar units
# and the atm pressure unit
const_Pstd = pstd #bar

#Molecules in a cubic meter at standard conditions
# The authoritative constant used for standard volume
const_dstd = const_pstd*1e5/const_Ru/const_Tstd  # 1/m^3

#Acceleration due to gravity on Earth
# Used by column height pressure units
const_g = g   # m/s^2

#Density of water (for water column pressure)
# Used by column water pressures
const_dh2o = dh2o

#Density of mercury (for mercury column pressure)
# Used by column mercury pressures
const_dhg = dhg

To change the standard conditions, call this function, and all 
unit conversion routines will be updated.
"""
    # Set up the exports as globals
    global const_h, const_k, const_Na, const_Nc, const_Ru, const_Tstd, \
            const_pstd, const_dstd, const_g, const_dh2o, const_dhg, const_q,\
            energy, length, time, mass, molar, temperature, \
            volume, pressure, force\
    # Plank: fundamental definition of kilogram from BIPM
    const_h = 6.62607015e-34 # Js
    #Boltzmann: fundamental definition of Kelvin from BIPM
    const_k = 1.380649e-23 # J/K
    #Avagadro: Fundamental definition of the mole from BIPM
    const_Na = 6.02214076e23 # count/mol
    #Coulumb: Fundamental definition of Ampere from BIPM
    const_q = 1.602176634e-19 # Fundamental charge in coulombs
    const_Nc = 1/const_q # count/coulomb
    #Universal gas
    const_Ru = const_k * const_Na # J/mol/K
    #Standard Temperature
    const_Tstd = Tstd #K
    # Standard Pressure (US Std atmosphere)
    const_pstd = pstd #bar
    # Standard molar concentration
    const_dstd = const_pstd*1e5/const_Ru/const_Tstd  # mol/m^3
    # Acceleration due to gravity on Earth
    const_g = g   # m/s^2
    # Density of water (for water column pressure)
    const_dh2o = dh2o
    # Density of mercury (for mercury column pressure)
    const_dhg = dhg


    # Validated 11/18/2017
    # Added angstroms 7/5/2021
    length = Conversion({
        'km':1e3,          # Kilometer
        'm':1.,             # Meter
        'cm':.01,          # Centimeter
        'mm':1e-3,         # Millimeter
        'um':1e-6,           # Micrometer
        'nm':1e-9,           # Nanometer
        'A':1e-10,           # Angstrom
        'in':.0254,         # Inch
        'nmi':1852,         # Nautical Mile
    }, 'unit_length')
    length['ft'] = length['in'] * 12.   # Foot
    length['yd'] = length['ft'] * 3.    # Yard
    length['mile'] = length['ft'] * 5280. # statute mile
    length['mi'] = length['mile']

    # Validated 11/18/2017    
    time = Conversion({
        'ns':1e9,           # Nanosecond
        'us':1e6,           # Microsecond
        'ms':1000.,         # Millisecond
        's':1.,             # Second
        'min':60.,       # Minute
        'hr':3600.,      # Hour
        'day':86400.,    # Day
        'year':3153600. # Julian year
    }, 'unit_time')

    # Validated 11/18/2017
    # Added u and amu, and changed to fundamental definition of lbm
    #    7/5/2021
    mass = Conversion({
        'kg':1.,            # Kilogram
        'g':1e-3,           # Gram
        'mg':1e-6,          # Milligram
        'lbm': 0.45359237,  # avoirdupois pound mass (precise by definition)
    }, 'unit_mass')
    mass['lb'] = mass['lbm']        # pound-mass
    mass['oz'] = mass['lbm']/16.    # ounce-mass
    # The slug is one pound adjusted by gc
    mass['slug'] = mass['lb'] * length['ft'] / (const_g * length['m'])
    # Define the atomic unit based on avagadro's number
    mass['u'] = 1./(const_Na * 1000)
    mass['amu'] = mass['u']

    # Validated 11/18/2017
    force = Conversion({
        'N':1.,             # Newton
        'kN':1000,          # Kilonewton
        'lb': const_g * mass['lb'], # Pounds
        'kgf':const_g       # Kilogram-force
    }, 'unit_force')
    force['lbf'] = force['lb']
    force['oz'] = force['lb']/16.

    # Define the molar to be consistently scaled with mass
    # This is essential for the matter function to operate correctly
    # Validated 11/18/2017
    molar = Conversion({
        'kmol':mass['kg'],              # Kilomole
        'mol':mass['g'],                # Mole
        'lbmol':mass['lbm'],            # Pound-mole
    }, 'unit_molar')
    molar['n'] = molar['mol']/const_Na  # Molecular count
    molar['Nm3'] = molar['mol']*const_dstd # Normal cubic meters
    molar['Ncum'] = molar['Nm3']        # Normal cubic meters
    molar['NL'] = molar['Nm3']*1e-3      # Normal liters
    molar['Ncc'] = molar['Nm3']*1e-6     # Normal cubic centimeters
    molar['scf'] = molar['Nm3'] * (length['ft'] / length['m'])**3
    molar['sci'] = molar['scf'] * (length['in'] / length['ft'])**3

    # Validated 11/18/2017
    # Adjusted so that Kelvin is 1.8 instead of Rankine 5/9
    # 7/6/2021
    temperature = Conversion({
        'K':1.8,                # Kelvin 
        'C':1.8,                # Celsius
        'F':1.,                 # Farenheit
        'R':1.,                 # Rankine
        'eV':1.8*const_q/const_k    # Electron-Volts
    }, 'unit_temperature')

    # Validated 11/18/2017
    energy = Conversion({
        'J':1.,             # Joule
        'kJ':1000.,         # Kilojoule
        'cal':4.184,        # Calorie (thermochemical)
        'kcal':4184.,       # Kilo Calorie (thermochemical)
        'eV':const_q        # Electron-Volt
    }, 'unit_energy')
    # Define the BTU as the energy 
    # British thermal unit (thermochemical)
    energy['BTU'] = energy['cal'] * (mass['lbm']/mass['g']) * (temperature['F']/temperature['C'])

    # Validated 11/18/2017
    volume = Conversion({
        'm3':length['m']**3,        # Cubic meters
        'mm3':length['mm']**3,      # Cubic millimeters
        'cm3':length['cm']**3,       # Cubic centimeter
        'in3':length['in']**3,      # Cubic inches
        'ft3':length['ft']**3,      # Cubic feet
    }, 'unit_volume')
    volume['L'] = volume['m3'] / 1000
    volume['mL'] = volume['L'] / 1000
    volume['uL'] = volume['L'] / 1e6
    volume['cum'] = volume['m3']
    volume['cc'] = volume['cm3']
    volume['cumm'] = volume['mm3']
    volume['cuin'] = volume['in3']
    volume['cuft'] = volume['ft3']
    volume['gal'] = volume['in3'] * 231.
    volume['USgal'] = volume['gal']
    volume['UKgal'] = volume['L'] * 4.54609 # Imperial gallons
    volume['qt'] = volume['gal'] / 4        # US liquid quart
    volume['pt'] = volume['gal'] / 8        # US liquid pint
    
    # Validated 11/18/2017
    pressure = Conversion({
        'Pa':1.,                # Pascals
        'kPa':1e3,              # Kilopascals
        'MPa':1e6,              # Megapascals
        'GPa':1e9,              # Gigapascals
        'bar':1e5,              # Bar
        'atm':1e5*const_pstd,           # Standard atmosphere
        'Torr':const_pstd * 1e5 / 760,  # Torr
        'mmHg':const_dhg*const_g*length['mm'],   # mm Mercury column height
        'mmH2O':const_dh2o*const_g*length['mm'], # mm column height water
        'psi':(force['lb']/force['N']) * (length['m'] / length['in'])**2,
    }, 'unit_pressure')
    pressure['psf'] = pressure['psi'] * (length['in'] / length['ft'])**2
    pressure['ksi'] = pressure['psi'] * 1000.   # kips per square inch
    pressure['inHg'] = pressure['mmHg'] * length['in'] / length['mm']  # in Mercury column height
    pressure['inH2O'] = pressure['mmH2O'] * length['in'] / length['mm'] # in water column height


# Execute the setup function
# This is executed BEFORE defining the other functions so the conversions
# will already be available
setup()



# Validated 11/18/2017
def temperature_scale(value, from_units=None, to_units=None, inplace=False):
    """Convert between tempertures scales
new_value = temperature_scale(old_value, from_units=None, to_units=None)

If the from_units or to_units are neglected, they will be pulled from the
'unit_temperature' pyro configuration parameter.  The temperature_scale()
function is used to convert values between scales.  The temperature() 
conversion routine converts differential or relative temperatures.  For 
example, temperature_scale() would be used to assert that 300K is 26.85C,
but temperature() could be used to show that specific heat is expressed
identically in J/kg/K and J/kg/C.

Instead of to_units or from_units, the string 'abs' may be passed.  This
implies that the temperature scale should be converted to or from its 
respective absolute scale.  This is useful to ensure an absolute 
temperature when the actual system of units is uncertain.  For example,
    T = temperature_scale(T, to_units = 'abs')
    [... some code ...]
    T = temperature_scale(T, from_units = 'abs')
This sequence of operations will convert a temperature from the system 
default to an absolute scale, perform some operation, and return to the
original default scale.
"""

    # Check to be certain from_units is specified
    if from_units is None:
        from_units = pyro.config['unit_temperature']
    # Repeat for to_units.  Are they specified?
    if to_units is None:
        to_units = pyro.config['unit_temperature']

    if from_units == to_units:
        return value
    elif from_units == 'abs':
        if to_units == 'C':
            from_units = 'K'
        elif to_units == 'F':
            from_units = 'R'
        else:
            return value
    elif to_units == 'abs':
        if from_units == 'C':
            to_units = 'K'
        elif from_units == 'F':
            to_units = 'R'
        else:
            return value

    if inplace:
        # point to the original if possible
        out = value
    else:
        # copy the original
        out = np.array(value, dtype=float)
        
    if from_units == 'C':
        np.add(out, 273.15, out=out)
    elif from_units == 'R':
        np.divide(out, 1.8, out=out)
    elif from_units == 'F':
        np.divide(out, 1.8, out=out)
        np.add(out, 255.3722222222222, out=out)
    elif from_units == 'eV':
        np.multiply(out, const_q / const_k, out=out)

    if to_units == 'C':
        np.subtract(out, 273.15, out=out)
    elif to_units == 'R':
        np.multiply(out, 1.8, out=out)
    elif to_units == 'F':
        np.subtract(out, 255.3722222222222, out=out)
        np.multiply(out, 1.8, out=out)
    elif to_units == 'eV':
        np.divide(out,const_q / const_k, out=out)
    return out

# Validated 11/18/2017
# Modified 7/4/2018 without validation -- added inplace operation
def gauge_to_abs(value, units=None, patm=const_pstd, inplace=False):
    """Adjust a gauge pressure to be in absolute units
new_value = gauge_to_abs(old_value, units=None, patm=const_pstd)

The 'units' keyword specifies the pressure units.  If it is omitted,
the 'units_pressure' PYroMat configuration parameter is used instead.

The 'patm' keyword specifies the atmospheric pressure IN BAR.  By 
default, the standard pressure is used.  If the atmospheric pressure
is already known in the same units as the gauge pressure, do not use
gauge_to_absolute(); simply add it to the old value.
"""
    if inplace:
        return np.add(
                value, 
                pressure(patm,from_units='bar',to_units=units),
                out=value)
    return np.add(
            value, 
            pressure(patm,from_units='bar',to_units=units))

# Validated 11/18/2017
# Modified 7/4/2018 without validation -- added inplace operation
def abs_to_gauge(value, units=None, patm=const_pstd):
    """Adjust an absolute pressure to be in gauge units
new_value = gauge_to_abs(old_value, units=None, patm=const_pstd)

The 'units' keyword specifies the pressure units.  If it is omitted,
the 'units_pressure' PYroMat configuration parameter is used instead.

The 'patm' keyword specifies the atmospheric pressure IN BAR.  By 
default, the standard pressure is used.  If the atmospheric pressure
is already known in the same units as the gauge pressure, do not use
gauge_to_absolute(); simply add it to the old value.
"""
    if inplace:
        return np.subtract(value, 
                pressure(patm,from_units='bar',to_units=units),
                out=value)
    return np.subtract(
            value, 
            pressure(patm,from_units='bar',to_units=units))

# Validated 11/18/2017
def matter(value, mw, from_units=None, to_units=None, exponent=None, inplace=False):
    """Convert between molar and mass units
new_value = matter(old_value, mw, 
                from_units=None, to_units=None, exponent=None)

This is a wrapper function for the mass() and molar() conversion routines
that will permit conversion between mass and molar units.  As a result, 
the molecular weight is required.

If the from_units or the to_units values are not specified, the pyromat
'unit_matter' parameter will be used.
"""

    # Check to be certain from_units is specified
    if from_units is None:
        from_units = pyro.utility.get_config('unit_matter')
    # Repeat for to_units.  Are they specified?
    if to_units is None:
        to_units = pyro.utility.get_config('unit_matter')

    if from_units == to_units:
        return value

    conv = 1.
    if from_units in mass:
        if to_units in mass:
            return mass(value, from_units, to_units, exponent, inplace=inplace)
        else:
            conv = mass[from_units] / mw / molar[to_units]
    else:
        if to_units in molar:
            return molar(value, from_units, to_units, exponent, inplace=inplace)
        else:
            conv = molar[from_units] * mw / mass[to_units]
    if exponent:
        conv**=exponent

    if inplace:
        return np.multiply(value, conv, out=value)
    return np.multiply(value, conv)


def show():
    """Print a summary of all unit conversions available"""
    for name,conv in pyro.units.__dict__.iteritems():
        if isinstance(conv,Conversion):
            sys.stdout.write( '%15s : '%name)
            for item in conv.get():
                sys.stdout.write(item + ' ')
            sys.stdout.write('\n')
    sys.stdout.write('See also...\n'
        '  abs_to_gauge()        Absolute to gauge pressure\n'+
        '  gauge_to_abs()        Gauge to absolute pressure\n'+
        '  matter()              Moles and mass conversions\n'+
        '  temperature_scale()   Correct handling of non-absolute temperatures\n')

