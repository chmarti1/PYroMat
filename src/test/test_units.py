import pyromat as pm
import pytest
import numpy as np
from pytest import approx

"""
Format of Reference Dicts is based upon return values from google.com of the 
form "1s in ns". Canonical values used where clealy available (e.g. seconds in
minute/hour). All values in each dict should be equivalent quantities after the
conversions are performed. This reflects how we would naturally form the 
equalities that we would use to represent units.
That is, we would represent the equivalency between m & cm as follows:
lengths['m'] == lengths['cm']
         1   ==         100         
"""

# ns us ms s min hr day year
times = {'s': 1,
         'ns': 1e9,
         'us': 1e6,
         'ms': 1e3,
         'min': 1/60,
         'hr': 1/(60*60),
         'day': 1/(60*60*24),
         'year': 1/(60*60*24*365)
         }
@pytest.fixture(params=list(times.keys()))
def timepair(request):
    return {'s': times['s'], request.param: times[request.param]}


# km m cm mm um nm A in nmi ft yd mile mi
lengths = {'m': 1,
           'km': 1e-3,
           'cm': 1e2,
           'mm': 1e3,
           'um': 1e6,
           'nm': 1e9,
           'A': 1e10,
           'in': 39.3701, #1/0.0254,
           'nmi': 1/1852,
           'ft': 39.3701/12, #1/0.0254/12,
           'yd': 39.3701/36, #1/0.0254/12/3,
           'mile': 0.000621371, #1/0.0254/12/5280,
           'mi': 0.000621371, #1/0.0254/12/5280,
           }
@pytest.fixture(params=list(lengths.keys()))
def lengthpair(request):
    return {'m': lengths['m'], request.param: lengths[request.param]}


# kg g mg lbm lb oz slug u amu
masses = {'kg': 1,
          'g': 1e3,
          'mg': 1e6,
          'lbm': 2.204623,
          'lb': 2.204623,
          'oz': 2.204623*16,
          'slug': 0.06852177,
          'u': pm.units.const_Na * 1000,
          'amu': pm.units.const_Na * 1000
          }
@pytest.fixture(params=list(masses.keys()))
def masspair(request):
    return {'kg': masses['kg'], request.param: masses[request.param]}


# N kN lb kgf lbf oz
forces = {'N': 1,
          'kN': 1e-3,
          'lb': 0.224809,
          'kgf': 1/pm.units.const_g,
          'lbf': 0.224809,
          'oz': 0.224809*16,
          }
@pytest.fixture(params=list(forces.keys()))
def forcepair(request):
    return {'N': forces['N'], request.param: forces[request.param]}


# kmol mol lbmol n Nm3 Ncum NL Ncc scf sci
molars = {'mol': 1,
          'kmol': 1e-3,
          'lbmol': 0.0022046226,
          'n': 1 * pm.units.const_Na,
          'Nm3': 1 / pm.units.const_dstd,
          'Ncum': 1 / pm.units.const_dstd,  # meters cubed occupied by 1 mol
          'NL': 1 / pm.units.const_dstd * 1000,
          'Ncc': 1 / pm.units.const_dstd * 1e6,
          'scf': 1 / pm.units.const_dstd * 35.3147,
          'sci': 1 / pm.units.const_dstd * 61023.7,
          }
@pytest.fixture(params=list(molars.keys()))
def molarpair(request):
    return {'mol': molars['mol'], request.param: molars[request.param]}


# K C F R eV
temperatures = {'K': 300,
                'C': 26.85,
                'F': 80.33,
                'R': 540,
                'eV': 0.02585198444922,
                }
@pytest.fixture(params=list(temperatures.keys()))
def temppair(request):
    return {'K': temperatures['K'], request.param: temperatures[request.param]}


# J kJ cal kcal eV BTU
energies = {'J': 1,
            'kJ': 1e-3,
            'cal': 1/4.184,
            'kcal': 1/4.184/1000,
            'eV': 1/pm.units.const_q,
            'BTU': 1/1054.35,  # defined via calories, not steam table val.
            }
@pytest.fixture(params=list(energies.keys()))
def energypair(request):
    return {'J': energies['J'], request.param: energies[request.param]}


# m3 mm3 cm3 in3 ft3 L mL uL cum cc cumm cuin cuft gal USgal UKgal qt pt
volumes = {'m3': 1,
           'mm3': 1e9,
           'cm3': 1e6,
           'in3': 61023.7,
           'ft3': 35.3147,
           'L': 1e3,
           'mL': 1e6,
           'uL': 1e9,
           'cum': 1,
           'cc': 1e6,
           'cumm': 1e9,
           'cuin': 61023.7,
           'cuft': 35.3147,
           'gal': 264.172,
           'USgal': 264.172,
           'UKgal': 219.969204701183,
           'qt': 1056.6879999969914934,
           'pt': 2113.3759999939829868,
           }
@pytest.fixture(params=list(volumes.keys()))
def volumepair(request):
    return {'m3': volumes['m3'], request.param: volumes[request.param]}


# Pa kPa MPa GPa bar atm Torr mmHg mmH2O psi psf ksi inHg inH2O
# https://www.justintools.com/unit-conversion/pressure.php
pressures = {'Pa': 1,
             'kPa': 1e-3,
             'MPa': 1e-6,
             'GPa': 1e-9,
             'bar': 1e-5,
             'atm': 1/101325,
             'Torr': 0.0075006168270417,
             'mmHg': 0.0075006168270417,
             'mmH2O': 0.10197442889221,
             'psi': 0.00014503773772954,
             'psf': 0.00014503773772954 * 144,
             'ksi': 0.00014503773772954 / 1000,
             'inHg': 0.0002953,
             'inH2O': 0.00401474,
             }
@pytest.fixture(params=list(pressures.keys()))
def pressurepair(request):
    return {'Pa': pressures['Pa'], request.param: pressures[request.param]}


# Manually concoct the matter values
mw = 12.2  # kg/kmol
matters = masses.copy()
for k in molars:
    matters[k] = molars[k] / mw / molars['kmol']  # manually convert
@pytest.fixture(params=list(matters.keys()))
def matterpair_kg(request):
    return {'kg': matters['kg'], request.param: matters[request.param]}
@pytest.fixture(params=list(matters.keys()))
def matterpair_kmol(request):
    return {'kmol': matters['kmol'], request.param: matters[request.param]}


# a function to wrap the manupulation that's needed for every unit type
def compare_both_ways(conv, pair, inplace=False):
    if len(pair.keys()) == 1:
        [unit_a] = list(pair.keys())
        unit_b = unit_a
        [val_a] = list(pair.values())
        val_b = val_a
    else:
        [unit_a, unit_b] = list(pair.keys())
        [val_a, val_b] = list(pair.values())

    if not inplace:
        assert conv(val_a, from_units=unit_a, to_units=unit_b) == approx(val_b)
        assert conv(val_b, from_units=unit_b, to_units=unit_a) == approx(val_a)
    else:
        assert conv(np.array(val_a, dtype=float), from_units=unit_a, to_units=unit_b, inplace=inplace) == approx(val_b)
        assert conv(np.array(val_b, dtype=float), from_units=unit_b, to_units=unit_a, inplace=inplace) == approx(val_a)
        # Could error due to lack of numpy type
        assert conv(val_a, from_units=unit_a, to_units=unit_b, inplace=inplace) == approx(val_b)


def test_time(timepair):
    compare_both_ways(pm.units.time, timepair)


def test_length(lengthpair):
    compare_both_ways(pm.units.length, lengthpair)


def test_mass(masspair):
    compare_both_ways(pm.units.mass, masspair)


def test_mass_inplace(masspair):
    compare_both_ways(pm.units.mass, masspair, inplace=True)


def test_force(forcepair):
    compare_both_ways(pm.units.force, forcepair)


def test_molar(molarpair):
    compare_both_ways(pm.units.molar, molarpair)


@pytest.mark.skip
def test_temperature():
    pass


def test_temperature_scale(temppair):
    compare_both_ways(pm.units.temperature_scale, temppair)


def test_temperature_scale_inplace(temppair):
    compare_both_ways(pm.units.temperature_scale, temppair, inplace=True)


def test_energy(energypair):
    compare_both_ways(pm.units.energy, energypair)


def test_volume(volumepair):
    compare_both_ways(pm.units.volume, volumepair)


def test_pressure(pressurepair):
    compare_both_ways(pm.units.pressure, pressurepair)


# Important to note that pressure conversions require patm in bar!
def test_abs_to_gauge():
    assert pm.units.abs_to_gauge(100, units='kPa', patm=1) == approx(0)
    assert pm.units.abs_to_gauge(0, units='kPa', patm=1) == approx(-100)
    assert pm.units.abs_to_gauge(500, units='kPa', patm=1) == approx(400)
    assert pm.units.abs_to_gauge(np.array(100, dtype=float), units='kPa', patm=1, inplace=True) == approx(0)
    # Could error due to lack of numpy type
    assert pm.units.abs_to_gauge(100, units='kPa', patm=1, inplace=True) == approx(0)


def test_gauge_to_abs():
    assert pm.units.gauge_to_abs(100, units='kPa', patm=1) == approx(200)
    assert pm.units.gauge_to_abs(-100, units='kPa', patm=1) == approx(0)
    assert pm.units.gauge_to_abs(500, units='kPa', patm=1) == approx(600)
    assert pm.units.gauge_to_abs(np.array(100, dtype=float), units='kPa', patm=1, inplace=True) == approx(200)
    # Could error due to lack of numpy type
    assert pm.units.gauge_to_abs(100, units='kPa', patm=1, inplace=True) == approx(200)


# We have to go manually here, because matter is a more complicated function call
def test_matter_kg(matterpair_kg):
    if len(matterpair_kg.keys()) == 1:
        [unit_a] = list(matterpair_kg.keys())
        unit_b = unit_a
        [val_a] = list(matterpair_kg.values())
        val_b = val_a
    else:
        [unit_a, unit_b] = list(matterpair_kg.keys())
        [val_a, val_b] = list(matterpair_kg.values())
    assert pm.units.matter(val_a, mw, from_units=unit_a, to_units=unit_b) == approx(val_b)
    assert pm.units.matter(val_b, mw, from_units=unit_b, to_units=unit_a) == approx(val_a)


def test_matter_kmol(matterpair_kmol):
    if len(matterpair_kmol.keys()) == 1:
        [unit_a] = list(matterpair_kmol.keys())
        unit_b = unit_a
        [val_a] = list(matterpair_kmol.values())
        val_b = val_a
    else:
        [unit_a, unit_b] = list(matterpair_kmol.keys())
        [val_a, val_b] = list(matterpair_kmol.values())
    assert pm.units.matter(val_a, mw, from_units=unit_a, to_units=unit_b) == approx(val_b)
    assert pm.units.matter(val_b, mw, from_units=unit_b, to_units=unit_a) == approx(val_a)


def test_matter_inplace(matterpair_kg):
    if len(matterpair_kg.keys()) == 1:
        [unit_a] = list(matterpair_kg.keys())
        unit_b = unit_a
        [val_a] = list(matterpair_kg.values())
        val_b = val_a
    else:
        [unit_a, unit_b] = list(matterpair_kg.keys())
        [val_a, val_b] = list(matterpair_kg.values())
    assert pm.units.matter(np.array(val_a, dtype=float), mw, from_units=unit_a, to_units=unit_b, inplace=True) == approx(val_b)
    assert pm.units.matter(np.array(val_b, dtype=float), mw, from_units=unit_b, to_units=unit_a, inplace=True) == approx(val_a)
    # Could error due to lack of numpy type
    assert pm.units.matter(val_a, mw, from_units=unit_a, to_units=unit_b, inplace=True) == approx(val_b)


def test_matter_exponent():
    mymw = 12.2  # kg/kmol
    unit_a = 'kmol'
    unit_b = 'kg'
    exp = -1
    val_a = 8314  # Ru in J/kmol K
    val_b = 8314 / mymw  # Rbar in J/kg K

    assert pm.units.matter(val_a, mymw, from_units=unit_a, to_units=unit_b, exponent=exp) == approx(val_b)
    assert pm.units.matter(val_b, mymw, from_units=unit_b, to_units=unit_a, exponent=exp) == approx(val_a)



