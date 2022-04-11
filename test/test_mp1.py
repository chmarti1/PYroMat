import pyromat as pm
import numpy as np
from pytest import approx, fixture, raises


@fixture
def water():
    w = pm.get('mp.H2O')
    return w


@fixture
def r134a():
    r = pm.get('mp.C2H2F4')
    return r


# Four reference points for H2O (NIST Webbook), and two saturation states
# Note that spec volumes are computed as 1/d for accuracy. NIST vals listed.
@fixture
def ref_a():
    # liquid
    vals = {'T': 300,        # K
            'p': 5,          # bar
            'd': 996.74,     # kg/m3
            'v': 1/996.74,   # 0.0010033,  # m3/kg
            'e': 112.52,     # kJ/kg
            'h': 113.02,     # kJ/kg
            's': 0.39295,    # kJ/kg K
            'x': -1
            }
    return vals


@fixture
def ref_b():
    # vapor
    vals = {'T': 600,        # K
            'p': 5,          # bar
            'd': 1.8242,     # kg/m3
            'v': 1/1.8242,   # 0.54820,    # m3/kg
            'e': 2846.0,     # kJ/kg
            'h': 3120.1,     # kJ/kg
            's': 7.5561,     # kJ/kg K
            'x': -1
            }
    return vals


@fixture
def ref_c():
    # mixture
    vals = {'T': 424.98,     # K
            'p': 5,          # bar
            'd': 5.321,      # kg/m3
            'v': 1/5.321,    # 0.18795,    # m3/kg
            'e': 1600.12,    # kJ/kg
            'h': 1694.095,   # kJ/kg
            's': 4.3406,     # kJ/kg K
            'x': 0.5,
            }
    return vals


@fixture
def ref_d():
    # supercritical
    vals = {'T': 800,        # K
            'p': 250,        # bar
            'd': 83.132,     # kg/m3
            'v': 1/83.132,   # 0.012029,   # m3/kg
            'e': 2961.5,     # kJ/kg
            'h': 3262.2,     # kJ/kg
            's': 6.0867,     # kJ/kg K
            'x': -1
            }
    return vals


@fixture
def ref_array(ref_a, ref_b, ref_c, ref_d):
    arr = {
        'T': [ref_a['T'], ref_b['T'], ref_c['T'], ref_d['T']],
        'p': [ref_a['p'], ref_b['p'], ref_c['p'], ref_d['p']],
        'd': [ref_a['d'], ref_b['d'], ref_c['d'], ref_d['d']],
        'v': [ref_a['v'], ref_b['v'], ref_c['v'], ref_d['v']],
        's': [ref_a['s'], ref_b['s'], ref_c['s'], ref_d['s']],
        'h': [ref_a['h'], ref_b['h'], ref_c['h'], ref_d['h']],
        'e': [ref_a['e'], ref_b['e'], ref_c['e'], ref_d['e']],
        'x': [ref_a['x'], ref_b['x'], ref_c['x'], ref_d['x']],
    }
    return arr


@fixture
def ref_sat_T():
    vals = {'T': 500,                    # K
            'p': 26.392,                 # bar
            'd': (831.31, 13.199),      # kg/m3
            'v': (0.012029, 0.075764),  # m3/kg
            'e': (972.26, 2602.5),      # kJ/kg
            'h': (975.43, 2802.5),      # kJ/kg
            's': (2.5810, 6.2351),      # kJ/kg K
            }
    return vals


@fixture
def ref_sat_p():
    vals = {'T': 424.98,                # K
            'p': 5.,                    # bar
            'd': (915.29, 2.6680),      # kg/m3
            'v': (0.0010925, 0.37481),  # m3/kg
            'e': (639.54, 2560.7),      # kJ/kg
            'h': (640.09, 2748.1),      # kJ/kg
            's': (1.8604, 6.8207),      # kJ/kg K
            }
    return vals


# ##### Test general functionality #####

def test_get():
    # Verify no errors
    mp1obj = pm.get('mp.H2O')
    mp1obj = pm.get('mp.C2H2F4')


# ##### Test accuracy of general calculations #####

def test_Tlim(water):
    lim = (273.16, 1273)
    assert water.Tlim() == approx(lim)
    assert water.Tlim(p=1) == approx(lim)


def test_plim(water):
    lim = (0, 10000)
    assert water.plim() == approx(lim)
    assert water.plim(T=1000) == approx(lim)


def test_critical(water):
    cr = (647.096, 220.64, 322.0)
    assert water.critical() == approx(cr[0:2])
    assert water.critical(density=True) == approx(cr)


def test_triple(water):
    tp = (273.16, 0.00611655)
    assert water.triple() == approx(tp)


# ##### Test accuracy to reference saturation data ##### #

def test_ref_sat_p_Ts(water, ref_sat_p):
    assert water.Ts(p=ref_sat_p['p']) == approx(ref_sat_p['T'], abs=1e-1)
    assert water.Ts(p=np.tile(ref_sat_p['p'],3)) == approx(np.tile(ref_sat_p['T'],3), abs=1e-1)


def test_ref_sat_p_ds(water, ref_sat_p):
    assert water.ds(p=ref_sat_p['p']) == approx(ref_sat_p['d'], abs=1e-1)
    assert water.ds(p=np.tile(ref_sat_p['p'], 3))[0] == approx(np.tile(ref_sat_p['d'][0], 3), abs=1e-1)
    assert water.ds(p=np.tile(ref_sat_p['p'], 3))[1] == approx(np.tile(ref_sat_p['d'][1], 3), abs=1e-1)


def test_ref_sat_p_es(water, ref_sat_p):
    assert water.es(p=ref_sat_p['p']) == approx(ref_sat_p['e'], abs=1e-1)
    assert water.es(p=np.tile(ref_sat_p['p'], 3))[0] == approx(np.tile(ref_sat_p['e'][0], 3), abs=1e-1)
    assert water.es(p=np.tile(ref_sat_p['p'], 3))[1] == approx(np.tile(ref_sat_p['e'][1], 3), abs=1e-1)


def test_ref_sat_p_hs(water, ref_sat_p):
    assert water.hs(p=ref_sat_p['p']) == approx(ref_sat_p['h'], abs=1e-1)
    assert water.hs(p=np.tile(ref_sat_p['p'], 3))[0] == approx(np.tile(ref_sat_p['h'][0], 3), abs=1e-1)
    assert water.hs(p=np.tile(ref_sat_p['p'], 3))[1] == approx(np.tile(ref_sat_p['h'][1], 3), abs=1e-1)


def test_ref_sat_p_ss(water, ref_sat_p):
    assert water.ss(p=ref_sat_p['p']) == approx(ref_sat_p['s'], abs=1e-1)
    assert water.ss(p=np.tile(ref_sat_p['p'], 3))[0] == approx(np.tile(ref_sat_p['s'][0], 3), abs=1e-1)
    assert water.ss(p=np.tile(ref_sat_p['p'], 3))[1] == approx(np.tile(ref_sat_p['s'][1], 3), abs=1e-1)


def test_ref_sat_T_ps(water, ref_sat_T):
    assert water.ps(T=ref_sat_T['T']) == approx(ref_sat_T['p'], abs=1e-1)
    assert water.ps(T=np.tile(ref_sat_T['T'],3)) == approx(np.tile(ref_sat_T['p'],3), abs=1e-1)


def test_ref_sat_T_ds(water, ref_sat_T):
    assert water.ds(T=ref_sat_T['T']) == approx(ref_sat_T['d'], abs=1e-1)
    assert water.ds(T=np.tile(ref_sat_T['T'], 3))[0] == approx(np.tile(ref_sat_T['d'][0], 3), abs=1e-1)
    assert water.ds(T=np.tile(ref_sat_T['T'], 3))[1] == approx(np.tile(ref_sat_T['d'][1], 3), abs=1e-1)


def test_ref_sat_T_es(water, ref_sat_T):
    assert water.es(T=ref_sat_T['T']) == approx(ref_sat_T['e'], abs=1e-1)
    assert water.es(T=np.tile(ref_sat_T['T'], 3))[0] == approx(np.tile(ref_sat_T['e'][0], 3), abs=1e-1)
    assert water.es(T=np.tile(ref_sat_T['T'], 3))[1] == approx(np.tile(ref_sat_T['e'][1], 3), abs=1e-1)


def test_ref_sat_T_hs(water, ref_sat_T):
    assert water.hs(T=ref_sat_T['T']) == approx(ref_sat_T['h'], abs=1e-1)
    assert water.hs(T=np.tile(ref_sat_T['T'], 3))[0] == approx(np.tile(ref_sat_T['h'][0], 3), abs=1e-1)
    assert water.hs(T=np.tile(ref_sat_T['T'], 3))[1] == approx(np.tile(ref_sat_T['h'][1], 3), abs=1e-1)


def test_ref_sat_T_ss(water, ref_sat_T):
    assert water.ss(T=ref_sat_T['T']) == approx(ref_sat_T['s'], abs=1e-1)
    assert water.ss(T=np.tile(ref_sat_T['T'], 3))[0] == approx(np.tile(ref_sat_T['s'][0], 3), abs=1e-1)
    assert water.ss(T=np.tile(ref_sat_T['T'], 3))[1] == approx(np.tile(ref_sat_T['s'][1], 3), abs=1e-1)


# ##### Test accuracy to reference state data ##### #


# ##### PRESSURE - REF A ####################################################################

def test_ref_a_p_Td(water, ref_a):
    assert water.p(T=ref_a['T'], d=ref_a['d']) == approx(ref_a['p'], abs=1e-1)
def test_ref_a_p_Tv(water, ref_a):
    assert water.p(T=ref_a['T'], v=ref_a['v']) == approx(ref_a['p'], abs=1e-1)
def test_ref_a_p_Ts(water, ref_a):
    assert water.p(T=ref_a['T'], s=ref_a['s']) == approx(ref_a['p'], abs=1e-1)
def test_ref_a_p_ds(water, ref_a):
    assert water.p(d=ref_a['d'], s=ref_a['s']) == approx(ref_a['p'], abs=1e-1)
def test_ref_a_p_vs(water, ref_a):
    assert water.p(v=ref_a['v'], s=ref_a['s']) == approx(ref_a['p'], abs=1e-1)
def test_ref_a_p_Th(water, ref_a):
    with raises(pm.utility.PMAnalysisError):
        assert water.p(T=ref_a['T'], h=ref_a['h']) == approx(ref_a['p'], abs=1e-1)
def test_ref_a_p_dh(water, ref_a):
    assert water.p(d=ref_a['d'], h=ref_a['h']) == approx(ref_a['p'], abs=1e-1)
def test_ref_a_p_vh(water, ref_a):
    assert water.p(v=ref_a['v'], h=ref_a['h']) == approx(ref_a['p'], abs=1e-1)
def test_ref_a_p_Te(water, ref_a):
    with raises(pm.utility.PMAnalysisError):
        assert water.p(T=ref_a['T'], e=ref_a['e']) == approx(ref_a['p'], abs=1e-1)
def test_ref_a_p_de(water, ref_a):
    assert water.p(d=ref_a['d'], e=ref_a['e']) == approx(ref_a['p'], abs=1e-1)
def test_ref_a_p_ve(water, ref_a):
    assert water.p(v=ref_a['v'], e=ref_a['e']) == approx(ref_a['p'], abs=1e-1)
def test_ref_a_p_Tx(water, ref_a):
    with raises(pm.utility.PMParamError):
        water.p(T=ref_a['T'], x=ref_a['x'])
def test_ref_a_p_dx(water, ref_a):
    with raises(pm.utility.PMParamError):
        water.p(d=ref_a['d'], x=ref_a['x'])
def test_ref_a_p_vx(water, ref_a):
    with raises(pm.utility.PMParamError):
        water.p(v=ref_a['v'], x=ref_a['x'])

# ##### PRESSURE - REF B ####################################################################

def test_ref_b_p_Td(water, ref_b):
    assert water.p(T=ref_b['T'], d=ref_b['d']) == approx(ref_b['p'], abs=1e-1)
def test_ref_b_p_Tv(water, ref_b):
    assert water.p(T=ref_b['T'], v=ref_b['v']) == approx(ref_b['p'], abs=1e-1)
def test_ref_b_p_Ts(water, ref_b):
    assert water.p(T=ref_b['T'], s=ref_b['s']) == approx(ref_b['p'], abs=1e-1)
def test_ref_b_p_ds(water, ref_b):
    assert water.p(d=ref_b['d'], s=ref_b['s']) == approx(ref_b['p'], abs=1e-1)
def test_ref_b_p_vs(water, ref_b):
    assert water.p(v=ref_b['v'], s=ref_b['s']) == approx(ref_b['p'], abs=1e-1)
def test_ref_b_p_Th(water, ref_b):
    with raises(pm.utility.PMAnalysisError):
        assert water.p(T=ref_b['T'], h=ref_b['h']) == approx(ref_b['p'], abs=1e-1)
def test_ref_b_p_dh(water, ref_b):
    assert water.p(d=ref_b['d'], h=ref_b['h']) == approx(ref_b['p'], abs=1e-1)
def test_ref_b_p_vh(water, ref_b):
    assert water.p(v=ref_b['v'], h=ref_b['h']) == approx(ref_b['p'], abs=1e-1)
def test_ref_b_p_Te(water, ref_b):
    with raises(pm.utility.PMAnalysisError):
        assert water.p(T=ref_b['T'], e=ref_b['e']) == approx(ref_b['p'], abs=1e-1)
def test_ref_b_p_de(water, ref_b):
    assert water.p(d=ref_b['d'], e=ref_b['e']) == approx(ref_b['p'], abs=1e-1)
def test_ref_b_p_ve(water, ref_b):
    assert water.p(v=ref_b['v'], e=ref_b['e']) == approx(ref_b['p'], abs=1e-1)
def test_ref_b_p_Tx(water, ref_b):
    with raises(pm.utility.PMParamError):
        water.p(T=ref_b['T'], x=ref_b['x'])
def test_ref_b_p_dx(water, ref_b):
    with raises(pm.utility.PMParamError):
        water.p(d=ref_b['d'], x=ref_b['x'])
def test_ref_b_p_vx(water, ref_b):
    with raises(pm.utility.PMParamError):
        water.p(v=ref_b['v'], x=ref_b['x'])

# ##### PRESSURE - REF C ####################################################################

def test_ref_c_p_Td(water, ref_c):
    assert water.p(T=ref_c['T'], d=ref_c['d']) == approx(ref_c['p'], abs=1e-1)
def test_ref_c_p_Tv(water, ref_c):
    assert water.p(T=ref_c['T'], v=ref_c['v']) == approx(ref_c['p'], abs=1e-1)
def test_ref_c_p_Ts(water, ref_c):
    assert water.p(T=ref_c['T'], s=ref_c['s']) == approx(ref_c['p'], abs=1e-1)
def test_ref_c_p_ds(water, ref_c):
    assert water.p(d=ref_c['d'], s=ref_c['s']) == approx(ref_c['p'], abs=1e-1)
def test_ref_c_p_vs(water, ref_c):
    assert water.p(v=ref_c['v'], s=ref_c['s']) == approx(ref_c['p'], abs=1e-1)
def test_ref_c_p_Th(water, ref_c):
    with raises(pm.utility.PMAnalysisError):
        assert water.p(T=ref_c['T'], h=ref_c['h']) == approx(ref_c['p'], abs=1e-1)
def test_ref_c_p_dh(water, ref_c):
    assert water.p(d=ref_c['d'], h=ref_c['h']) == approx(ref_c['p'], abs=1e-1)
def test_ref_c_p_vh(water, ref_c):
    assert water.p(v=ref_c['v'], h=ref_c['h']) == approx(ref_c['p'], abs=1e-1)
def test_ref_c_p_Te(water, ref_c):
    with raises(pm.utility.PMAnalysisError):
        assert water.p(T=ref_c['T'], e=ref_c['e']) == approx(ref_c['p'], abs=1e-1)
def test_ref_c_p_de(water, ref_c):
    assert water.p(d=ref_c['d'], e=ref_c['e']) == approx(ref_c['p'], abs=1e-1)
def test_ref_c_p_ve(water, ref_c):
    assert water.p(v=ref_c['v'], e=ref_c['e']) == approx(ref_c['p'], abs=1e-1)
def test_ref_c_p_Tx(water, ref_c):
    assert water.p(T=ref_c['T'], x=ref_c['x']) == approx(ref_c['p'], abs=1e-1)
def test_ref_c_p_dx(water, ref_c):
    with raises(pm.utility.PMParamError):
        assert water.p(d=ref_c['d'], x=ref_c['x'])
def test_ref_c_p_vx(water, ref_c):
    with raises(pm.utility.PMParamError):
        assert water.p(v=ref_c['v'], x=ref_c['x'])

# ##### PRESSURE - REF D ####################################################################

def test_ref_d_p_Td(water, ref_d):
    assert water.p(T=ref_d['T'], d=ref_d['d']) == approx(ref_d['p'], abs=1e-1)
def test_ref_d_p_Tv(water, ref_d):
    assert water.p(T=ref_d['T'], v=ref_d['v']) == approx(ref_d['p'], abs=1e-1)
def test_ref_d_p_Ts(water, ref_d):
    assert water.p(T=ref_d['T'], s=ref_d['s']) == approx(ref_d['p'], abs=1e-1)
def test_ref_d_p_ds(water, ref_d):
    assert water.p(d=ref_d['d'], s=ref_d['s']) == approx(ref_d['p'], abs=1e-1)
def test_ref_d_p_vs(water, ref_d):
    assert water.p(v=ref_d['v'], s=ref_d['s']) == approx(ref_d['p'], abs=1e-1)
def test_ref_d_p_Th(water, ref_d):
    with raises(pm.utility.PMAnalysisError):
        assert water.p(T=ref_d['T'], h=ref_d['h']) == approx(ref_d['p'], abs=1e-1)
def test_ref_d_p_dh(water, ref_d):
    assert water.p(d=ref_d['d'], h=ref_d['h']) == approx(ref_d['p'], abs=1e-1)
def test_ref_d_p_vh(water, ref_d):
    assert water.p(v=ref_d['v'], h=ref_d['h']) == approx(ref_d['p'], abs=1e-1)
def test_ref_d_p_Te(water, ref_d):
    with raises(pm.utility.PMAnalysisError):
        assert water.p(T=ref_d['T'], e=ref_d['e']) == approx(ref_d['p'], abs=1e-1)
def test_ref_d_p_de(water, ref_d):
    assert water.p(d=ref_d['d'], e=ref_d['e']) == approx(ref_d['p'], abs=1e-1)
def test_ref_d_p_ve(water, ref_d):
    assert water.p(v=ref_d['v'], e=ref_d['e']) == approx(ref_d['p'], abs=1e-1)
def test_ref_d_p_Tx(water, ref_d):
    with raises(pm.utility.PMParamError):
        water.p(T=ref_d['T'], x=ref_d['x'])
def test_ref_d_p_dx(water, ref_d):
    with raises(pm.utility.PMParamError):
        water.p(d=ref_d['d'], x=ref_d['x'])
def test_ref_d_p_vx(water, ref_d):
    with raises(pm.utility.PMParamError):
        water.p(v=ref_d['v'], x=ref_d['x'])

# ##### TEMPERATURE - REF A ####################################################################

def test_ref_a_T_pd(water, ref_a):
    assert water.T(p=ref_a['p'], d=ref_a['d']) == approx(ref_a['T'], abs=1e-1)
def test_ref_a_T_pv(water, ref_a):
    assert water.T(p=ref_a['p'], v=ref_a['v']) == approx(ref_a['T'], abs=1e-1)
def test_ref_a_T_ps(water, ref_a):
    assert water.T(p=ref_a['p'], s=ref_a['s']) == approx(ref_a['T'], abs=1e-1)
def test_ref_a_T_ds(water, ref_a):
    assert water.T(s=ref_a['s'], d=ref_a['d']) == approx(ref_a['T'], abs=1e-1)
def test_ref_a_T_vs(water, ref_a):
    assert water.T(s=ref_a['s'], v=ref_a['v']) == approx(ref_a['T'], abs=1e-1)
def test_ref_a_T_ph(water, ref_a):
    assert water.T(h=ref_a['h'], p=ref_a['p']) == approx(ref_a['T'], abs=1e-1)
def test_ref_a_T_dh(water, ref_a):
    assert water.T(h=ref_a['h'], d=ref_a['d']) == approx(ref_a['T'], abs=1e-1)
def test_ref_a_T_vh(water, ref_a):
    assert water.T(h=ref_a['h'], v=ref_a['v']) == approx(ref_a['T'], abs=1e-1)
def test_ref_a_T_pe(water, ref_a):
    assert water.T(e=ref_a['e'], p=ref_a['p']) == approx(ref_a['T'], abs=1e-1)
def test_ref_a_T_de(water, ref_a):
    assert water.T(e=ref_a['e'], d=ref_a['d']) == approx(ref_a['T'], abs=1e-1)
def test_ref_a_T_ve(water, ref_a):
    assert water.T(e=ref_a['e'], v=ref_a['v']) == approx(ref_a['T'], abs=1e-1)
def test_ref_a_T_px(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.T(x=ref_a['x'], p=ref_a['p'])
def test_ref_a_T_dx(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.T(x=ref_a['x'], d=ref_a['d'])
def test_ref_a_T_vx(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.T(x=ref_a['x'], v=ref_a['v'])


# ##### TEMPERATURE - REF B ####################################################################

def test_ref_b_T_pd(water, ref_b):
    assert water.T(p=ref_b['p'], d=ref_b['d']) == approx(ref_b['T'], abs=1e-1)
def test_ref_b_T_pv(water, ref_b):
    assert water.T(p=ref_b['p'], v=ref_b['v']) == approx(ref_b['T'], abs=1e-1)
def test_ref_b_T_ps(water, ref_b):
    assert water.T(p=ref_b['p'], s=ref_b['s']) == approx(ref_b['T'], abs=1e-1)
def test_ref_b_T_ds(water, ref_b):
    assert water.T(s=ref_b['s'], d=ref_b['d']) == approx(ref_b['T'], abs=1e-1)
def test_ref_b_T_vs(water, ref_b):
    assert water.T(s=ref_b['s'], v=ref_b['v']) == approx(ref_b['T'], abs=1e-1)
def test_ref_b_T_ph(water, ref_b):
    assert water.T(h=ref_b['h'], p=ref_b['p']) == approx(ref_b['T'], abs=1e-1)
def test_ref_b_T_dh(water, ref_b):
    assert water.T(h=ref_b['h'], d=ref_b['d']) == approx(ref_b['T'], abs=1e-1)
def test_ref_b_T_vh(water, ref_b):
    assert water.T(h=ref_b['h'], v=ref_b['v']) == approx(ref_b['T'], abs=1e-1)
def test_ref_b_T_pe(water, ref_b):
    assert water.T(e=ref_b['e'], p=ref_b['p']) == approx(ref_b['T'], abs=1e-1)
def test_ref_b_T_de(water, ref_b):
    assert water.T(e=ref_b['e'], d=ref_b['d']) == approx(ref_b['T'], abs=1e-1)
def test_ref_b_T_ve(water, ref_b):
    assert water.T(e=ref_b['e'], v=ref_b['v']) == approx(ref_b['T'], abs=1e-1)
def test_ref_b_T_px(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.T(x=ref_b['x'], p=ref_b['p'])
def test_ref_b_T_dx(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.T(x=ref_b['x'], d=ref_b['d'])
def test_ref_b_T_vx(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.T(x=ref_b['x'], v=ref_b['v'])

# ##### TEMPERATURE - REF C ####################################################################

def test_ref_c_T_pd(water, ref_c):
    assert water.T(p=ref_c['p'], d=ref_c['d']) == approx(ref_c['T'], abs=1e-1)
def test_ref_c_T_pv(water, ref_c):
    assert water.T(p=ref_c['p'], v=ref_c['v']) == approx(ref_c['T'], abs=1e-1)
def test_ref_c_T_ps(water, ref_c):
    assert water.T(p=ref_c['p'], s=ref_c['s']) == approx(ref_c['T'], abs=1e-1)
def test_ref_c_T_ds(water, ref_c):
    assert water.T(s=ref_c['s'], d=ref_c['d']) == approx(ref_c['T'], abs=1e-1)
def test_ref_c_T_vs(water, ref_c):
    assert water.T(s=ref_c['s'], v=ref_c['v']) == approx(ref_c['T'], abs=1e-1)
def test_ref_c_T_ph(water, ref_c):
    assert water.T(h=ref_c['h'], p=ref_c['p']) == approx(ref_c['T'], abs=1e-1)
def test_ref_c_T_dh(water, ref_c):
    assert water.T(h=ref_c['h'], d=ref_c['d']) == approx(ref_c['T'], abs=1e-1)
def test_ref_c_T_vh(water, ref_c):
    assert water.T(h=ref_c['h'], v=ref_c['v']) == approx(ref_c['T'], abs=1e-1)
def test_ref_c_T_pe(water, ref_c):
    assert water.T(e=ref_c['e'], p=ref_c['p']) == approx(ref_c['T'], abs=1e-1)
def test_ref_c_T_de(water, ref_c):
    assert water.T(e=ref_c['e'], d=ref_c['d']) == approx(ref_c['T'], abs=1e-1)
def test_ref_c_T_ve(water, ref_c):
    assert water.T(e=ref_c['e'], v=ref_c['v']) == approx(ref_c['T'], abs=1e-1)
def test_ref_c_T_px(water, ref_c):
    assert water.T(x=ref_c['x'], p=ref_c['p']) == approx(ref_c['T'], abs=1e-1)
def test_ref_c_T_dx(water, ref_c):
    with raises(pm.utility.PMParamError):
        assert water.T(x=ref_c['x'], d=ref_c['d'])
def test_ref_c_T_vx(water, ref_c):
    with raises(pm.utility.PMParamError):
        assert water.T(x=ref_c['x'], v=ref_c['v'])

# ##### TEMPERATURE - REF D ####################################################################

def test_ref_d_T_pd(water, ref_d):
    assert water.T(p=ref_d['p'], d=ref_d['d']) == approx(ref_d['T'], abs=1e-1)
def test_ref_d_T_pv(water, ref_d):
    assert water.T(p=ref_d['p'], v=ref_d['v']) == approx(ref_d['T'], abs=1e-1)
def test_ref_d_T_ps(water, ref_d):
    assert water.T(p=ref_d['p'], s=ref_d['s']) == approx(ref_d['T'], abs=1e-1)
def test_ref_d_T_ds(water, ref_d):
    assert water.T(s=ref_d['s'], d=ref_d['d']) == approx(ref_d['T'], abs=1e-1)
def test_ref_d_T_vs(water, ref_d):
    assert water.T(s=ref_d['s'], v=ref_d['v']) == approx(ref_d['T'], abs=1e-1)
def test_ref_d_T_ph(water, ref_d):
    assert water.T(h=ref_d['h'], p=ref_d['p']) == approx(ref_d['T'], abs=1e-1)
def test_ref_d_T_dh(water, ref_d):
    assert water.T(h=ref_d['h'], d=ref_d['d']) == approx(ref_d['T'], abs=1e-1)
def test_ref_d_T_vh(water, ref_d):
    assert water.T(h=ref_d['h'], v=ref_d['v']) == approx(ref_d['T'], abs=1e-1)
def test_ref_d_T_pe(water, ref_d):
    assert water.T(e=ref_d['e'], p=ref_d['p']) == approx(ref_d['T'], abs=1e-1)
def test_ref_d_T_de(water, ref_d):
    assert water.T(e=ref_d['e'], d=ref_d['d']) == approx(ref_d['T'], abs=1e-1)
def test_ref_d_T_ve(water, ref_d):
    assert water.T(e=ref_d['e'], v=ref_d['v']) == approx(ref_d['T'], abs=1e-1)
def test_ref_d_T_px(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.T(x=ref_d['x'], p=ref_d['p'])
def test_ref_d_T_dx(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.T(x=ref_d['x'], d=ref_d['d'])
def test_ref_d_T_vx(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.T(x=ref_d['x'], v=ref_d['v'])

# ##### DENSITY - REF A ####################################################################

def test_ref_a_d_pT(water, ref_a):
    assert water.d(p=ref_a['p'], T=ref_a['T']) == approx(ref_a['d'], abs=1e-2)
def test_ref_a_d_ps(water, ref_a):
    assert water.d(p=ref_a['p'], s=ref_a['s']) == approx(ref_a['d'], abs=1e-2)
def test_ref_a_d_Ts(water, ref_a):
    assert water.d(s=ref_a['s'], T=ref_a['T']) == approx(ref_a['d'], abs=1e-2)
def test_ref_a_d_ph(water, ref_a):
    assert water.d(p=ref_a['p'], h=ref_a['h']) == approx(ref_a['d'], abs=1e-2)
def test_ref_a_d_Th(water, ref_a):
    with raises(pm.utility.PMAnalysisError):
        assert water.d(T=ref_a['T'], h=ref_a['h']) == approx(ref_a['d'], abs=1e-2)
def test_ref_a_d_pe(water, ref_a):
    assert water.d(p=ref_a['p'], e=ref_a['e']) == approx(ref_a['d'], abs=1e-2)
def test_ref_a_d_Te(water, ref_a):
    with raises(pm.utility.PMAnalysisError):
        assert water.d(T=ref_a['T'], e=ref_a['e']) == approx(ref_a['d'], abs=1e-2)
def test_ref_a_d_px(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.d(p=ref_a['p'], x=ref_a['x'])
def test_ref_a_d_Tx(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.d(T=ref_a['T'], x=ref_a['x'])

# ##### DENSITY - REF B ####################################################################

def test_ref_b_d_pT(water, ref_b):
    assert water.d(p=ref_b['p'], T=ref_b['T']) == approx(ref_b['d'], abs=1e-2)
def test_ref_b_d_ps(water, ref_b):
    assert water.d(p=ref_b['p'], s=ref_b['s']) == approx(ref_b['d'], abs=1e-2)
def test_ref_b_d_Ts(water, ref_b):
    assert water.d(s=ref_b['s'], T=ref_b['T']) == approx(ref_b['d'], abs=1e-2)
def test_ref_b_d_ph(water, ref_b):
    assert water.d(p=ref_b['p'], h=ref_b['h']) == approx(ref_b['d'], abs=1e-2)
def test_ref_b_d_Th(water, ref_b):
    with raises(pm.utility.PMAnalysisError):
        assert water.d(T=ref_b['T'], h=ref_b['h']) == approx(ref_b['d'], abs=1e-2)
def test_ref_b_d_pe(water, ref_b):
    assert water.d(p=ref_b['p'], e=ref_b['e']) == approx(ref_b['d'], abs=1e-2)
def test_ref_b_d_Te(water, ref_b):
    with raises(pm.utility.PMAnalysisError):
        assert water.d(T=ref_b['T'], e=ref_b['e']) == approx(ref_b['d'], abs=1e-2)
def test_ref_b_d_px(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.d(p=ref_b['p'], x=ref_b['x'])
def test_ref_b_d_Tx(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.d(T=ref_b['T'], x=ref_b['x'])

# ##### DENSITY - REF C ####################################################################

def test_ref_c_d_pTx(water, ref_c):
    assert water.d(p=ref_c['p'], T=ref_c['T'], x=ref_c['x']) == approx(ref_c['d'], abs=1e-2)
def test_ref_c_d_ps(water, ref_c):
    assert water.d(p=ref_c['p'], s=ref_c['s']) == approx(ref_c['d'], abs=1e-2)
def test_ref_c_d_Ts(water, ref_c):
    assert water.d(s=ref_c['s'], T=ref_c['T']) == approx(ref_c['d'], abs=1e-2)
def test_ref_c_d_ph(water, ref_c):
    assert water.d(p=ref_c['p'], h=ref_c['h']) == approx(ref_c['d'], abs=1e-2)
def test_ref_c_d_Th(water, ref_c):
    with raises(pm.utility.PMAnalysisError):
        assert water.d(T=ref_c['T'], h=ref_c['h']) == approx(ref_c['d'], abs=1e-2)
def test_ref_c_d_pe(water, ref_c):
    assert water.d(p=ref_c['p'], e=ref_c['e']) == approx(ref_c['d'], abs=1e-2)
def test_ref_c_d_Te(water, ref_c):
    with raises(pm.utility.PMAnalysisError):
        assert water.d(T=ref_c['T'], e=ref_c['e']) == approx(ref_c['d'], abs=1e-2)
def test_ref_c_d_px(water, ref_c):
    assert water.d(p=ref_c['p'], x=ref_c['x']) == approx(ref_c['d'], abs=1e-2)
def test_ref_c_d_Tx(water, ref_c):
    assert water.d(T=ref_c['T'], x=ref_c['x']) == approx(ref_c['d'], abs=1e-2)

# ##### DENSITY - REF D ####################################################################

def test_ref_d_d_pT(water, ref_d):
    assert water.d(p=ref_d['p'], T=ref_d['T']) == approx(ref_d['d'], abs=1e-2)
def test_ref_d_d_ps(water, ref_d):
    assert water.d(p=ref_d['p'], s=ref_d['s']) == approx(ref_d['d'], abs=1e-2)
def test_ref_d_d_Ts(water, ref_d):
    assert water.d(s=ref_d['s'], T=ref_d['T']) == approx(ref_d['d'], abs=1e-2)
def test_ref_d_d_ph(water, ref_d):
    assert water.d(p=ref_d['p'], h=ref_d['h']) == approx(ref_d['d'], abs=1e-2)
def test_ref_d_d_Th(water, ref_d):
    with raises(pm.utility.PMAnalysisError):
        assert water.d(T=ref_d['T'], h=ref_d['h']) == approx(ref_d['d'], abs=1e-2)
def test_ref_d_d_pe(water, ref_d):
    assert water.d(p=ref_d['p'], e=ref_d['e']) == approx(ref_d['d'], abs=1e-2)
def test_ref_d_d_Te(water, ref_d):
    with raises(pm.utility.PMAnalysisError):
        assert water.d(T=ref_d['T'], e=ref_d['e']) == approx(ref_d['d'], abs=1e-2)
def test_ref_d_d_px(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.d(p=ref_d['p'], x=ref_d['x'])
def test_ref_d_d_Tx(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.d(T=ref_d['T'], x=ref_d['x'])

# ##### ENTROPY - REF A ####################################################################
def test_ref_a_s_pT(water, ref_a):
    assert water.s(p=ref_a['p'], T=ref_a['T']) == approx(ref_a['s'], abs=1e-1)
def test_ref_a_s_pd(water, ref_a):
    assert water.s(p=ref_a['p'], d=ref_a['d']) == approx(ref_a['s'], abs=1e-1)
def test_ref_a_s_pv(water, ref_a):
    assert water.s(p=ref_a['p'], v=ref_a['v']) == approx(ref_a['s'], abs=1e-1)
def test_ref_a_s_Td(water, ref_a):
    assert water.s(T=ref_a['T'], d=ref_a['d']) == approx(ref_a['s'], abs=1e-1)
def test_ref_a_s_Tv(water, ref_a):
    assert water.s(T=ref_a['T'], v=ref_a['v']) == approx(ref_a['s'], abs=1e-1)
def test_ref_a_s_ph(water, ref_a):
    assert water.s(h=ref_a['h'], p=ref_a['p']) == approx(ref_a['s'], abs=1e-1)
def test_ref_a_s_Th(water, ref_a):
    with raises(pm.utility.PMAnalysisError):
        assert water.s(h=ref_a['h'], T=ref_a['T']) == approx(ref_a['s'], abs=1e-1)
def test_ref_a_s_dh(water, ref_a):
    assert water.s(h=ref_a['h'], d=ref_a['d']) == approx(ref_a['s'], abs=1e-1)
def test_ref_a_s_vh(water, ref_a):
    assert water.s(h=ref_a['h'], v=ref_a['v']) == approx(ref_a['s'], abs=1e-1)
def test_ref_a_s_pe(water, ref_a):
    assert water.s(e=ref_a['e'], p=ref_a['p']) == approx(ref_a['s'], abs=1e-1)
def test_ref_a_s_Te(water, ref_a):
    with raises(pm.utility.PMAnalysisError):
        assert water.s(e=ref_a['e'], T=ref_a['T']) == approx(ref_a['s'], abs=1e-1)
def test_ref_a_s_de(water, ref_a):
    assert water.s(e=ref_a['e'], d=ref_a['d']) == approx(ref_a['s'], abs=1e-1)
def test_ref_a_s_ve(water, ref_a):
    assert water.s(e=ref_a['e'], v=ref_a['v']) == approx(ref_a['s'], abs=1e-1)
def test_ref_a_s_px(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.s(x=ref_a['x'], p=ref_a['p'])
def test_ref_a_s_Tx(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.s(x=ref_a['x'], T=ref_a['T'])
def test_ref_a_s_dx(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.s(x=ref_a['x'], d=ref_a['d'])
def test_ref_a_s_vx(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.s(x=ref_a['x'], v=ref_a['v'])
# ##### ENTROPY - REF B ####################################################################
def test_ref_b_s_pT(water, ref_b):
    assert water.s(p=ref_b['p'], T=ref_b['T']) == approx(ref_b['s'], abs=1e-1)
def test_ref_b_s_pd(water, ref_b):
    assert water.s(p=ref_b['p'], d=ref_b['d']) == approx(ref_b['s'], abs=1e-1)
def test_ref_b_s_pv(water, ref_b):
    assert water.s(p=ref_b['p'], v=ref_b['v']) == approx(ref_b['s'], abs=1e-1)
def test_ref_b_s_Td(water, ref_b):
    assert water.s(T=ref_b['T'], d=ref_b['d']) == approx(ref_b['s'], abs=1e-1)
def test_ref_b_s_Tv(water, ref_b):
    assert water.s(T=ref_b['T'], v=ref_b['v']) == approx(ref_b['s'], abs=1e-1)
def test_ref_b_s_ph(water, ref_b):
    assert water.s(h=ref_b['h'], p=ref_b['p']) == approx(ref_b['s'], abs=1e-1)
def test_ref_b_s_Th(water, ref_b):
    with raises(pm.utility.PMAnalysisError):
        assert water.s(h=ref_b['h'], T=ref_b['T']) == approx(ref_b['s'], abs=1e-1)
def test_ref_b_s_dh(water, ref_b):
    assert water.s(h=ref_b['h'], d=ref_b['d']) == approx(ref_b['s'], abs=1e-1)
def test_ref_b_s_vh(water, ref_b):
    assert water.s(h=ref_b['h'], v=ref_b['v']) == approx(ref_b['s'], abs=1e-1)
def test_ref_b_s_pe(water, ref_b):
    assert water.s(e=ref_b['e'], p=ref_b['p']) == approx(ref_b['s'], abs=1e-1)
def test_ref_b_s_Te(water, ref_b):
    with raises(pm.utility.PMAnalysisError):
        assert water.s(e=ref_b['e'], T=ref_b['T']) == approx(ref_b['s'], abs=1e-1)
def test_ref_b_s_de(water, ref_b):
    assert water.s(e=ref_b['e'], d=ref_b['d']) == approx(ref_b['s'], abs=1e-1)
def test_ref_b_s_ve(water, ref_b):
    assert water.s(e=ref_b['e'], v=ref_b['v']) == approx(ref_b['s'], abs=1e-1)
def test_ref_b_s_px(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.s(x=ref_b['x'], p=ref_b['p'])
def test_ref_b_s_Tx(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.s(x=ref_b['x'], T=ref_b['T'])
def test_ref_b_s_dx(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.s(x=ref_b['x'], d=ref_b['d'])
def test_ref_b_s_vx(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.s(x=ref_b['x'], v=ref_b['v'])
# ##### ENTROPY - REF C ####################################################################
def test_ref_c_s_pTx(water, ref_c):
    assert water.s(p=ref_c['p'], T=ref_c['T'], x=ref_c['x']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_pd(water, ref_c):
    assert water.s(p=ref_c['p'], d=ref_c['d']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_pv(water, ref_c):
    assert water.s(p=ref_c['p'], v=ref_c['v']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_Td(water, ref_c):
    assert water.s(T=ref_c['T'], d=ref_c['d']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_Tv(water, ref_c):
    assert water.s(T=ref_c['T'], v=ref_c['v']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_ph(water, ref_c):
    assert water.s(h=ref_c['h'], p=ref_c['p']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_Th(water, ref_c):
    with raises(pm.utility.PMAnalysisError):
        assert water.s(h=ref_c['h'], T=ref_c['T']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_dh(water, ref_c):
    assert water.s(h=ref_c['h'], d=ref_c['d']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_vh(water, ref_c):
    assert water.s(h=ref_c['h'], v=ref_c['v']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_pe(water, ref_c):
    assert water.s(e=ref_c['e'], p=ref_c['p']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_Te(water, ref_c):
    with raises(pm.utility.PMAnalysisError):
        assert water.s(e=ref_c['e'], T=ref_c['T']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_de(water, ref_c):
    assert water.s(e=ref_c['e'], d=ref_c['d']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_ve(water, ref_c):
    assert water.s(e=ref_c['e'], v=ref_c['v']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_px(water, ref_c):
    assert water.s(p=ref_c['p'], x=ref_c['x']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_Tx(water, ref_c):
    assert water.s(T=ref_c['T'], x=ref_c['x']) == approx(ref_c['s'], abs=1e-1)
def test_ref_c_s_dx(water, ref_c):
    with raises(pm.utility.PMParamError):
        assert water.s(d=ref_c['d'], x=ref_c['x'])
def test_ref_c_s_vx(water, ref_c):
    with raises(pm.utility.PMParamError):
        assert water.s(v=ref_c['v'], x=ref_c['x'])
# ##### ENTROPY - REF D ####################################################################
def test_ref_d_s_pT(water, ref_d):
    assert water.s(p=ref_d['p'], T=ref_d['T']) == approx(ref_d['s'], abs=1e-1)
def test_ref_d_s_pd(water, ref_d):
    assert water.s(p=ref_d['p'], d=ref_d['d']) == approx(ref_d['s'], abs=1e-1)
def test_ref_d_s_pv(water, ref_d):
    assert water.s(p=ref_d['p'], v=ref_d['v']) == approx(ref_d['s'], abs=1e-1)
def test_ref_d_s_Td(water, ref_d):
    assert water.s(T=ref_d['T'], d=ref_d['d']) == approx(ref_d['s'], abs=1e-1)
def test_ref_d_s_Tv(water, ref_d):
    assert water.s(T=ref_d['T'], v=ref_d['v']) == approx(ref_d['s'], abs=1e-1)
def test_ref_d_s_ph(water, ref_d):
    assert water.s(h=ref_d['h'], p=ref_d['p']) == approx(ref_d['s'], abs=1e-1)
def test_ref_d_s_Th(water, ref_d):
    with raises(pm.utility.PMAnalysisError):
        assert water.s(h=ref_d['h'], T=ref_d['T']) == approx(ref_d['s'], abs=1e-1)
def test_ref_d_s_dh(water, ref_d):
    assert water.s(h=ref_d['h'], d=ref_d['d']) == approx(ref_d['s'], abs=1e-1)
def test_ref_d_s_vh(water, ref_d):
    assert water.s(h=ref_d['h'], v=ref_d['v']) == approx(ref_d['s'], abs=1e-1)
def test_ref_d_s_pe(water, ref_d):
    assert water.s(e=ref_d['e'], p=ref_d['p']) == approx(ref_d['s'], abs=1e-1)
def test_ref_d_s_Te(water, ref_d):
    with raises(pm.utility.PMAnalysisError):
        assert water.s(e=ref_d['e'], T=ref_d['T']) == approx(ref_d['s'], abs=1e-1)
def test_ref_d_s_de(water, ref_d):
    assert water.s(e=ref_d['e'], d=ref_d['d']) == approx(ref_d['s'], abs=1e-1)
def test_ref_d_s_ve(water, ref_d):
    assert water.s(e=ref_d['e'], v=ref_d['v']) == approx(ref_d['s'], abs=1e-1)
def test_ref_d_s_px(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.s(x=ref_d['x'], p=ref_d['p'])
def test_ref_d_s_Tx(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.s(x=ref_d['x'], T=ref_d['T'])
def test_ref_d_s_dx(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.s(x=ref_d['x'], d=ref_d['d'])
def test_ref_d_s_vx(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.s(x=ref_d['x'], v=ref_d['v'])



# ##### ENTHALPY - REF A ####################################################################
def test_ref_a_h_pT(water, ref_a):
    assert water.h(p=ref_a['p'], T=ref_a['T']) == approx(ref_a['h'], abs=1e-1)
def test_ref_a_h_pd(water, ref_a):
    assert water.h(p=ref_a['p'], d=ref_a['d']) == approx(ref_a['h'], abs=1e-1)
def test_ref_a_h_pv(water, ref_a):
    assert water.h(p=ref_a['p'], v=ref_a['v']) == approx(ref_a['h'], abs=1e-1)
def test_ref_a_h_Td(water, ref_a):
    assert water.h(T=ref_a['T'], d=ref_a['d']) == approx(ref_a['h'], abs=1e-1)
def test_ref_a_h_Tv(water, ref_a):
    assert water.h(T=ref_a['T'], v=ref_a['v']) == approx(ref_a['h'], abs=1e-1)
def test_ref_a_h_ps(water, ref_a):
    assert water.h(s=ref_a['s'], p=ref_a['p']) == approx(ref_a['h'], abs=1e-1)
def test_ref_a_h_Ts(water, ref_a):
    assert water.h(s=ref_a['s'], T=ref_a['T']) == approx(ref_a['h'], abs=1e-1)
def test_ref_a_h_ds(water, ref_a):
    assert water.h(s=ref_a['s'], d=ref_a['d']) == approx(ref_a['h'], abs=1e-1)
def test_ref_a_h_vs(water, ref_a):
    assert water.h(s=ref_a['s'], v=ref_a['v']) == approx(ref_a['h'], abs=1e-1)
def test_ref_a_h_pe(water, ref_a):
    assert water.h(e=ref_a['e'], p=ref_a['p']) == approx(ref_a['h'], abs=1e-1)
def test_ref_a_h_Te(water, ref_a):
    with raises(pm.utility.PMAnalysisError):
        assert water.h(e=ref_a['e'], T=ref_a['T']) == approx(ref_a['h'], abs=1e-1)
def test_ref_a_h_de(water, ref_a):
    assert water.h(e=ref_a['e'], d=ref_a['d']) == approx(ref_a['h'], abs=1e-1)
def test_ref_a_h_ve(water, ref_a):
    assert water.h(e=ref_a['e'], v=ref_a['v']) == approx(ref_a['h'], abs=1e-1)
def test_ref_a_h_px(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.h(x=ref_a['x'], p=ref_a['p'])
def test_ref_a_h_Tx(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.h(x=ref_a['x'], T=ref_a['T'])
def test_ref_a_h_dx(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.h(x=ref_a['x'], d=ref_a['d'])
def test_ref_a_h_vx(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.h(x=ref_a['x'], v=ref_a['v'])

# ##### ENTHALPY - REF B ####################################################################
def test_ref_b_h_pT(water, ref_b):
    assert water.h(p=ref_b['p'], T=ref_b['T']) == approx(ref_b['h'], abs=1e-1)
def test_ref_b_h_pd(water, ref_b):
    assert water.h(p=ref_b['p'], d=ref_b['d']) == approx(ref_b['h'], abs=1e-1)
def test_ref_b_h_pv(water, ref_b):
    assert water.h(p=ref_b['p'], v=ref_b['v']) == approx(ref_b['h'], abs=1e-1)
def test_ref_b_h_Td(water, ref_b):
    assert water.h(T=ref_b['T'], d=ref_b['d']) == approx(ref_b['h'], abs=1e-1)
def test_ref_b_h_Tv(water, ref_b):
    assert water.h(T=ref_b['T'], v=ref_b['v']) == approx(ref_b['h'], abs=1e-1)
def test_ref_b_h_ps(water, ref_b):
    assert water.h(s=ref_b['s'], p=ref_b['p']) == approx(ref_b['h'], abs=1e-1)
def test_ref_b_h_Ts(water, ref_b):
    assert water.h(s=ref_b['s'], T=ref_b['T']) == approx(ref_b['h'], abs=1e-1)
def test_ref_b_h_ds(water, ref_b):
    assert water.h(s=ref_b['s'], d=ref_b['d']) == approx(ref_b['h'], abs=1e-1)
def test_ref_b_h_vs(water, ref_b):
    assert water.h(s=ref_b['s'], v=ref_b['v']) == approx(ref_b['h'], abs=1e-1)
def test_ref_b_h_pe(water, ref_b):
    assert water.h(e=ref_b['e'], p=ref_b['p']) == approx(ref_b['h'], abs=1e-1)
def test_ref_b_h_Te(water, ref_b):
    with raises(pm.utility.PMAnalysisError):
        assert water.h(e=ref_b['e'], T=ref_b['T']) == approx(ref_b['h'], abs=1e-1)
def test_ref_b_h_de(water, ref_b):
    assert water.h(e=ref_b['e'], d=ref_b['d']) == approx(ref_b['h'], abs=1e-1)
def test_ref_b_h_ve(water, ref_b):
    assert water.h(e=ref_b['e'], v=ref_b['v']) == approx(ref_b['h'], abs=1e-1)
def test_ref_b_h_px(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.h(x=ref_b['x'], p=ref_b['p'])
def test_ref_b_h_Tx(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.h(x=ref_b['x'], T=ref_b['T'])
def test_ref_b_h_dx(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.h(x=ref_b['x'], d=ref_b['d'])
def test_ref_b_h_vx(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.h(x=ref_b['x'], v=ref_b['v'])
# ##### ENTHALPY - REF C ####################################################################
def test_ref_c_h_pTx(water, ref_c):
    assert water.h(p=ref_c['p'], T=ref_c['T'], x=ref_c['x']) == approx(ref_c['h'], abs=1e-1)
def test_ref_c_h_pd(water, ref_c):
    assert water.h(p=ref_c['p'], d=ref_c['d']) == approx(ref_c['h'], abs=5e-1)
def test_ref_c_h_pv(water, ref_c):
    assert water.h(p=ref_c['p'], v=ref_c['v']) == approx(ref_c['h'], abs=5e-1)
def test_ref_c_h_Td(water, ref_c):
    assert water.h(T=ref_c['T'], d=ref_c['d']) == approx(ref_c['h'], abs=5e-1)
def test_ref_c_h_Tv(water, ref_c):
    assert water.h(T=ref_c['T'], v=ref_c['v']) == approx(ref_c['h'], abs=5e-1)
def test_ref_c_h_ps(water, ref_c):
    assert water.h(s=ref_c['s'], p=ref_c['p']) == approx(ref_c['h'], abs=1e-1)
def test_ref_c_h_Ts(water, ref_c):
    assert water.h(s=ref_c['s'], T=ref_c['T']) == approx(ref_c['h'], abs=1e-1)
def test_ref_c_h_ds(water, ref_c):
    assert water.h(s=ref_c['s'], d=ref_c['d']) == approx(ref_c['h'], abs=1e-1)
def test_ref_c_h_vs(water, ref_c):
    assert water.h(s=ref_c['s'], v=ref_c['v']) == approx(ref_c['h'], abs=1e-1)
def test_ref_c_h_pe(water, ref_c):
    assert water.h(e=ref_c['e'], p=ref_c['p']) == approx(ref_c['h'], abs=1e-1)
def test_ref_c_h_Te(water, ref_c):
    with raises(pm.utility.PMAnalysisError):
        assert water.h(e=ref_c['e'], T=ref_c['T']) == approx(ref_c['h'], abs=1e-1)
def test_ref_c_h_de(water, ref_c):
    assert water.h(e=ref_c['e'], d=ref_c['d']) == approx(ref_c['h'], abs=1e-1)
def test_ref_c_h_ve(water, ref_c):
    assert water.h(e=ref_c['e'], v=ref_c['v']) == approx(ref_c['h'], abs=1e-1)
def test_ref_c_h_px(water, ref_c):
    assert water.h(p=ref_c['p'], x=ref_c['x']) == approx(ref_c['h'], abs=1e-1)
def test_ref_c_h_Tx(water, ref_c):
    assert water.h(T=ref_c['T'], x=ref_c['x']) == approx(ref_c['h'], abs=1e-1)
def test_ref_c_h_dx(water, ref_c):
    with raises(pm.utility.PMParamError):
        assert water.h(d=ref_c['d'], x=ref_c['x'])
def test_ref_c_h_vx(water, ref_c):
    with raises(pm.utility.PMParamError):
        assert water.h(v=ref_c['v'], x=ref_c['x'])
# ##### ENTHALPY - REF D ####################################################################
def test_ref_d_h_pT(water, ref_d):
    assert water.h(p=ref_d['p'], T=ref_d['T']) == approx(ref_d['h'], abs=1e-1)
def test_ref_d_h_pd(water, ref_d):
    assert water.h(p=ref_d['p'], d=ref_d['d']) == approx(ref_d['h'], abs=1e-1)
def test_ref_d_h_pv(water, ref_d):
    assert water.h(p=ref_d['p'], v=ref_d['v']) == approx(ref_d['h'], abs=1e-1)
def test_ref_d_h_Td(water, ref_d):
    assert water.h(T=ref_d['T'], d=ref_d['d']) == approx(ref_d['h'], abs=1e-1)
def test_ref_d_h_Tv(water, ref_d):
    assert water.h(T=ref_d['T'], v=ref_d['v']) == approx(ref_d['h'], abs=1e-1)
def test_ref_d_h_ps(water, ref_d):
    assert water.h(s=ref_d['s'], p=ref_d['p']) == approx(ref_d['h'], abs=1e-1)
def test_ref_d_h_Ts(water, ref_d):
    assert water.h(s=ref_d['s'], T=ref_d['T']) == approx(ref_d['h'], abs=1e-1)
def test_ref_d_h_ds(water, ref_d):
    assert water.h(s=ref_d['s'], d=ref_d['d']) == approx(ref_d['h'], abs=1e-1)
def test_ref_d_h_vs(water, ref_d):
    assert water.h(s=ref_d['s'], v=ref_d['v']) == approx(ref_d['h'], abs=1e-1)
def test_ref_d_h_pe(water, ref_d):
    assert water.h(e=ref_d['e'], p=ref_d['p']) == approx(ref_d['h'], abs=1e-1)
def test_ref_d_h_Te(water, ref_d):
    with raises(pm.utility.PMAnalysisError):
        assert water.h(e=ref_d['e'], T=ref_d['T']) == approx(ref_d['h'], abs=1e-1)
def test_ref_d_h_de(water, ref_d):
    assert water.h(e=ref_d['e'], d=ref_d['d']) == approx(ref_d['h'], abs=1e-1)
def test_ref_d_h_ve(water, ref_d):
    assert water.h(e=ref_d['e'], v=ref_d['v']) == approx(ref_d['h'], abs=1e-1)
def test_ref_d_h_px(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.h(x=ref_d['x'], p=ref_d['p'])
def test_ref_d_h_Tx(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.h(x=ref_d['x'], T=ref_d['T'])
def test_ref_d_h_dx(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.h(x=ref_d['x'], d=ref_d['d'])
def test_ref_d_h_vx(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.h(x=ref_d['x'], v=ref_d['v'])


# ##### ENERGY - REF A ####################################################################
def test_ref_a_e_pT(water, ref_a):
    assert water.e(p=ref_a['p'], T=ref_a['T']) == approx(ref_a['e'], abs=1e-1)
def test_ref_a_e_pd(water, ref_a):
    assert water.e(p=ref_a['p'], d=ref_a['d']) == approx(ref_a['e'], abs=1e-1)
def test_ref_a_e_pv(water, ref_a):
    assert water.e(p=ref_a['p'], v=ref_a['v']) == approx(ref_a['e'], abs=1e-1)
def test_ref_a_e_Td(water, ref_a):
    assert water.e(T=ref_a['T'], d=ref_a['d']) == approx(ref_a['e'], abs=1e-1)
def test_ref_a_e_Tv(water, ref_a):
    assert water.e(T=ref_a['T'], v=ref_a['v']) == approx(ref_a['e'], abs=1e-1)
def test_ref_a_e_ps(water, ref_a):
    assert water.e(s=ref_a['s'], p=ref_a['p']) == approx(ref_a['e'], abs=1e-1)
def test_ref_a_e_Ts(water, ref_a):
    assert water.e(s=ref_a['s'], T=ref_a['T']) == approx(ref_a['e'], abs=1e-1)
def test_ref_a_e_ds(water, ref_a):
    assert water.e(s=ref_a['s'], d=ref_a['d']) == approx(ref_a['e'], abs=1e-1)
def test_ref_a_e_vs(water, ref_a):
    assert water.e(s=ref_a['s'], v=ref_a['v']) == approx(ref_a['e'], abs=1e-1)
def test_ref_a_e_ph(water, ref_a):
    assert water.e(h=ref_a['h'], p=ref_a['p']) == approx(ref_a['e'], abs=1e-1)
def test_ref_a_e_Th(water, ref_a):
    with raises(pm.utility.PMAnalysisError):
        assert water.e(h=ref_a['h'], T=ref_a['T']) == approx(ref_a['e'], abs=1e-1)
def test_ref_a_e_dh(water, ref_a):
    assert water.e(h=ref_a['h'], d=ref_a['d']) == approx(ref_a['e'], abs=1e-1)
def test_ref_a_e_vh(water, ref_a):
    assert water.e(h=ref_a['h'], v=ref_a['v']) == approx(ref_a['e'], abs=1e-1)
def test_ref_a_e_px(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.e(x=ref_a['x'], p=ref_a['p'])
def test_ref_a_e_Tx(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.e(x=ref_a['x'], T=ref_a['T'])
def test_ref_a_e_dx(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.e(x=ref_a['x'], d=ref_a['d'])
def test_ref_a_e_vx(water, ref_a):
    with raises(pm.utility.PMParamError):
        assert water.e(x=ref_a['x'], v=ref_a['v'])


# ##### ENERGY - REF B ####################################################################
def test_ref_b_e_pT(water, ref_b):
    assert water.e(p=ref_b['p'], T=ref_b['T']) == approx(ref_b['e'], abs=1e-1)
def test_ref_b_e_pd(water, ref_b):
    assert water.e(p=ref_b['p'], d=ref_b['d']) == approx(ref_b['e'], abs=1e-1)
def test_ref_b_e_pv(water, ref_b):
    assert water.e(p=ref_b['p'], v=ref_b['v']) == approx(ref_b['e'], abs=1e-1)
def test_ref_b_e_Td(water, ref_b):
    assert water.e(T=ref_b['T'], d=ref_b['d']) == approx(ref_b['e'], abs=1e-1)
def test_ref_b_e_Tv(water, ref_b):
    assert water.e(T=ref_b['T'], v=ref_b['v']) == approx(ref_b['e'], abs=1e-1)
def test_ref_b_e_ps(water, ref_b):
    assert water.e(s=ref_b['s'], p=ref_b['p']) == approx(ref_b['e'], abs=1e-1)
def test_ref_b_e_Ts(water, ref_b):
    assert water.e(s=ref_b['s'], T=ref_b['T']) == approx(ref_b['e'], abs=1e-1)
def test_ref_b_e_ds(water, ref_b):
    assert water.e(s=ref_b['s'], d=ref_b['d']) == approx(ref_b['e'], abs=1e-1)
def test_ref_b_e_vs(water, ref_b):
    assert water.e(s=ref_b['s'], v=ref_b['v']) == approx(ref_b['e'], abs=1e-1)
def test_ref_b_e_ph(water, ref_b):
    assert water.e(h=ref_b['h'], p=ref_b['p']) == approx(ref_b['e'], abs=1e-1)
def test_ref_b_e_Th(water, ref_b):
    with raises(pm.utility.PMAnalysisError):
        assert water.e(h=ref_b['h'], T=ref_b['T']) == approx(ref_b['e'], abs=1e-1)
def test_ref_b_e_dh(water, ref_b):
    assert water.e(h=ref_b['h'], d=ref_b['d']) == approx(ref_b['e'], abs=1e-1)
def test_ref_b_e_vh(water, ref_b):
    assert water.e(h=ref_b['h'], v=ref_b['v']) == approx(ref_b['e'], abs=1e-1)
def test_ref_b_e_px(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.e(x=ref_b['x'], p=ref_b['p'])
def test_ref_b_e_Tx(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.e(x=ref_b['x'], T=ref_b['T'])
def test_ref_b_e_dx(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.e(x=ref_b['x'], d=ref_b['d'])
def test_ref_b_e_vx(water, ref_b):
    with raises(pm.utility.PMParamError):
        assert water.e(x=ref_b['x'], v=ref_b['v'])

# ##### ENERGY - REF C ####################################################################
def test_ref_c_e_pTx(water, ref_c):
    assert water.e(p=ref_c['p'], T=ref_c['T'], x=ref_c['x']) == approx(ref_c['e'], abs=1e-1)
def test_ref_c_e_pd(water, ref_c):
    assert water.e(p=ref_c['p'], d=ref_c['d']) == approx(ref_c['e'], abs=5e-1)
def test_ref_c_e_pv(water, ref_c):
    assert water.e(p=ref_c['p'], v=ref_c['v']) == approx(ref_c['e'], abs=5e-1)
def test_ref_c_e_Td(water, ref_c):
    assert water.e(T=ref_c['T'], d=ref_c['d']) == approx(ref_c['e'], abs=5e-1)
def test_ref_c_e_Tv(water, ref_c):
    assert water.e(T=ref_c['T'], v=ref_c['v']) == approx(ref_c['e'], abs=5e-1)
def test_ref_c_e_ps(water, ref_c):
    assert water.e(s=ref_c['s'], p=ref_c['p']) == approx(ref_c['e'], abs=1e-1)
def test_ref_c_e_Ts(water, ref_c):
    assert water.e(s=ref_c['s'], T=ref_c['T']) == approx(ref_c['e'], abs=1e-1)
def test_ref_c_e_ds(water, ref_c):
    assert water.e(s=ref_c['s'], d=ref_c['d']) == approx(ref_c['e'], abs=1e-1)
def test_ref_c_e_vs(water, ref_c):
    assert water.e(s=ref_c['s'], v=ref_c['v']) == approx(ref_c['e'], abs=1e-1)
def test_ref_c_e_ph(water, ref_c):
    assert water.e(h=ref_c['h'], p=ref_c['p']) == approx(ref_c['e'], abs=1e-1)
def test_ref_c_e_Th(water, ref_c):
    with raises(pm.utility.PMAnalysisError):
        assert water.e(h=ref_c['h'], T=ref_c['T']) == approx(ref_c['e'], abs=1e-1)
def test_ref_c_e_dh(water, ref_c):
    assert water.e(h=ref_c['h'], d=ref_c['d']) == approx(ref_c['e'], abs=1e-1)
def test_ref_c_e_vh(water, ref_c):
    assert water.e(h=ref_c['h'], v=ref_c['v']) == approx(ref_c['e'], abs=1e-1)
def test_ref_c_e_px(water, ref_c):
    assert water.e(p=ref_c['p'], x=ref_c['x']) == approx(ref_c['e'], abs=1e-1)
def test_ref_c_e_Tx(water, ref_c):
    assert water.e(T=ref_c['T'], x=ref_c['x']) == approx(ref_c['e'], abs=1e-1)
def test_ref_c_e_dx(water, ref_c):
    with raises(pm.utility.PMParamError):
        assert water.e(d=ref_c['d'], x=ref_c['x'])
def test_ref_c_e_vx(water, ref_c):
    with raises(pm.utility.PMParamError):
        assert water.e(v=ref_c['v'], x=ref_c['x'])
# ##### ENERGY - REF D ####################################################################
def test_ref_d_e_pT(water, ref_d):
    assert water.e(p=ref_d['p'], T=ref_d['T']) == approx(ref_d['e'], abs=1e-1)
def test_ref_d_e_pd(water, ref_d):
    assert water.e(p=ref_d['p'], d=ref_d['d']) == approx(ref_d['e'], abs=1e-1)
def test_ref_d_e_pv(water, ref_d):
    assert water.e(p=ref_d['p'], v=ref_d['v']) == approx(ref_d['e'], abs=1e-1)
def test_ref_d_e_Td(water, ref_d):
    assert water.e(T=ref_d['T'], d=ref_d['d']) == approx(ref_d['e'], abs=1e-1)
def test_ref_d_e_Tv(water, ref_d):
    assert water.e(T=ref_d['T'], v=ref_d['v']) == approx(ref_d['e'], abs=1e-1)
def test_ref_d_e_ps(water, ref_d):
    assert water.e(s=ref_d['s'], p=ref_d['p']) == approx(ref_d['e'], abs=1e-1)
def test_ref_d_e_Ts(water, ref_d):
    assert water.e(s=ref_d['s'], T=ref_d['T']) == approx(ref_d['e'], abs=1e-1)
def test_ref_d_e_ds(water, ref_d):
    assert water.e(s=ref_d['s'], d=ref_d['d']) == approx(ref_d['e'], abs=1e-1)
def test_ref_d_e_vs(water, ref_d):
    assert water.e(s=ref_d['s'], v=ref_d['v']) == approx(ref_d['e'], abs=1e-1)
def test_ref_d_e_ph(water, ref_d):
    assert water.e(h=ref_d['h'], p=ref_d['p']) == approx(ref_d['e'], abs=1e-1)
def test_ref_d_e_Th(water, ref_d):
    with raises(pm.utility.PMAnalysisError):
        assert water.e(h=ref_d['h'], T=ref_d['T']) == approx(ref_d['e'], abs=1e-1)
def test_ref_d_e_dh(water, ref_d):
    assert water.e(h=ref_d['h'], d=ref_d['d']) == approx(ref_d['e'], abs=1e-1)
def test_ref_d_e_vh(water, ref_d):
    assert water.e(h=ref_d['h'], v=ref_d['v']) == approx(ref_d['e'], abs=1e-1)
def test_ref_d_e_px(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.e(x=ref_d['x'], p=ref_d['p'])
def test_ref_d_e_Tx(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.e(x=ref_d['x'], T=ref_d['T'])
def test_ref_d_e_dx(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.e(x=ref_d['x'], d=ref_d['d'])
def test_ref_d_e_vx(water, ref_d):
    with raises(pm.utility.PMParamError):
        assert water.e(x=ref_d['x'], v=ref_d['v'])


# ##### Test combinations and different data formats ##### #

# ##### PRESSURE ####################################################################
def test_multi_P_Td(water, ref_array):
    assert water.p(T=ref_array['T'], d=ref_array['d']) == approx(ref_array['p'], abs=1e-1)
def test_multi_P_Tv(water, ref_array):
    assert water.p(T=ref_array['T'], v=ref_array['v']) == approx(ref_array['p'], abs=1e-1)
def test_multi_P_Ts(water, ref_array):
    assert water.p(T=ref_array['T'], s=ref_array['s']) == approx(ref_array['p'], abs=2.5e-1)
def test_multi_P_ds(water, ref_array):
    assert water.p(d=ref_array['d'], s=ref_array['s']) == approx(ref_array['p'], abs=1e-1)
def test_multi_P_vs(water, ref_array):
    assert water.p(v=ref_array['v'], s=ref_array['s']) == approx(ref_array['p'], abs=1e-1)
def test_multi_P_Th(water, ref_array):
    with raises(pm.utility.PMAnalysisError):
        assert water.p(T=ref_array['T'], h=ref_array['h']) == approx(ref_array['p'], abs=1e-1)
def test_multi_P_dh(water, ref_array):
    assert water.p(d=ref_array['d'], h=ref_array['h']) == approx(ref_array['p'], abs=1e-1)
def test_multi_P_vh(water, ref_array):
    assert water.p(v=ref_array['v'], h=ref_array['h']) == approx(ref_array['p'], abs=1e-1)
def test_multi_P_Te(water, ref_array):
    with raises(pm.utility.PMAnalysisError):
        assert water.p(T=ref_array['T'], e=ref_array['e']) == approx(ref_array['p'], abs=1e-1)
def test_multi_P_de(water, ref_array):
    assert water.p(d=ref_array['d'], e=ref_array['e']) == approx(ref_array['p'], abs=1e-1)
def test_multi_P_ve(water, ref_array):
    assert water.p(v=ref_array['v'], e=ref_array['e']) == approx(ref_array['p'], abs=1e-1)
def test_multi_P_Tx(water, ref_array):
    with raises(pm.utility.PMParamError):  # Some x values are illegal
        water.p(T=ref_array['T'], x=ref_array['x'])

# ##### TEMPERATURE ####################################################################
def test_multi_T_pd(water, ref_array):
    assert water.T(p=ref_array['p'], d=ref_array['d']) == approx(ref_array['T'], abs=1e-1)
def test_multi_T_pv(water, ref_array):
    assert water.T(p=ref_array['p'], v=ref_array['v']) == approx(ref_array['T'], abs=1e-1)
def test_multi_T_ps(water, ref_array):
    assert water.T(p=ref_array['p'], s=ref_array['s']) == approx(ref_array['T'], abs=1e-1)
def test_multi_T_ds(water, ref_array):
    assert water.T(d=ref_array['d'], s=ref_array['s']) == approx(ref_array['T'], abs=1e-1)
def test_multi_T_vs(water, ref_array):
    assert water.T(v=ref_array['v'], s=ref_array['s']) == approx(ref_array['T'], abs=1e-1)
def test_multi_T_ph(water, ref_array):
    assert water.T(p=ref_array['p'], h=ref_array['h']) == approx(ref_array['T'], abs=1e-1)
def test_multi_T_dh(water, ref_array):
    assert water.T(d=ref_array['d'], h=ref_array['h']) == approx(ref_array['T'], abs=1e-1)
def test_multi_T_vh(water, ref_array):
    assert water.T(v=ref_array['v'], h=ref_array['h']) == approx(ref_array['T'], abs=1e-1)
def test_multi_T_pe(water, ref_array):
    assert water.T(p=ref_array['p'], e=ref_array['e']) == approx(ref_array['T'], abs=1e-1)
def test_multi_T_de(water, ref_array):
    assert water.T(d=ref_array['d'], e=ref_array['e']) == approx(ref_array['T'], abs=1e-1)
def test_multi_T_ve(water, ref_array):
    assert water.T(v=ref_array['v'], e=ref_array['e']) == approx(ref_array['T'], abs=1e-1)
def test_multi_T_px(water, ref_array):
    with raises(pm.utility.PMParamError):  # Some x values are illegal
        water.T(p=ref_array['p'], x=ref_array['x'])


# ##### DENSITY ####################################################################
def test_multi_d_pTx(water, ref_array):
    assert water.d(p=ref_array['p'], T=ref_array['T'], x=ref_array['x']) == approx(ref_array['d'], abs=1e-1)
def test_multi_d_ps(water, ref_array):
    assert water.d(p=ref_array['p'], s=ref_array['s']) == approx(ref_array['d'], abs=1e-1)
def test_multi_d_Ts(water, ref_array):
    assert water.d(T=ref_array['T'], s=ref_array['s']) == approx(ref_array['d'], abs=1e-1)
def test_multi_d_ph(water, ref_array):
    assert water.d(p=ref_array['p'], h=ref_array['h']) == approx(ref_array['d'], abs=1e-1)
def test_multi_d_Th(water, ref_array):
    with raises(pm.utility.PMAnalysisError):
        assert water.d(T=ref_array['T'], h=ref_array['h']) == approx(ref_array['d'], abs=1e-1)
def test_multi_d_pe(water, ref_array):
    assert water.d(p=ref_array['p'], e=ref_array['e']) == approx(ref_array['d'], abs=1e-1)
def test_multi_d_Te(water, ref_array):
    with raises(pm.utility.PMAnalysisError):
        assert water.d(T=ref_array['T'], e=ref_array['e']) == approx(ref_array['d'], abs=1e-1)

# ##### ENTROPY ####################################################################
def test_multi_s_pTx(water, ref_array):
    assert water.s(p=ref_array['p'], T=ref_array['T'], x=ref_array['x']) == approx(ref_array['s'], abs=1e-1)
def test_multi_s_pd(water, ref_array):
    assert water.s(p=ref_array['p'], d=ref_array['d']) == approx(ref_array['s'], abs=1e-1)
def test_multi_s_pv(water, ref_array):
    assert water.s(p=ref_array['p'], v=ref_array['v']) == approx(ref_array['s'], abs=1e-1)
def test_multi_s_Td(water, ref_array):
    assert water.s(T=ref_array['T'], d=ref_array['d']) == approx(ref_array['s'], abs=1e-1)
def test_multi_s_Tv(water, ref_array):
    assert water.s(T=ref_array['T'], v=ref_array['v']) == approx(ref_array['s'], abs=1e-1)
def test_multi_s_ph(water, ref_array):
    assert water.s(p=ref_array['p'], h=ref_array['h']) == approx(ref_array['s'], abs=1e-1)
def test_multi_s_Th(water, ref_array):
    with raises(pm.utility.PMAnalysisError):
        assert water.s(T=ref_array['T'], h=ref_array['h']) == approx(ref_array['s'], abs=1e-1)
def test_multi_s_dh(water, ref_array):
    assert water.s(d=ref_array['d'], h=ref_array['h']) == approx(ref_array['s'], abs=1e-1)
def test_multi_s_vh(water, ref_array):
    assert water.s(v=ref_array['v'], h=ref_array['h']) == approx(ref_array['s'], abs=1e-1)
def test_multi_s_pe(water, ref_array):
    assert water.s(p=ref_array['p'], e=ref_array['e']) == approx(ref_array['s'], abs=1e-1)
def test_multi_s_Te(water, ref_array):
    with raises(pm.utility.PMAnalysisError):
        assert water.s(T=ref_array['T'], e=ref_array['e']) == approx(ref_array['s'], abs=1e-1)
def test_multi_s_de(water, ref_array):
    assert water.s(d=ref_array['d'], e=ref_array['e']) == approx(ref_array['s'], abs=1e-1)
def test_multi_s_ve(water, ref_array):
    assert water.s(v=ref_array['v'], e=ref_array['e']) == approx(ref_array['s'], abs=1e-1)


# ##### ENTHALPY ####################################################################
def test_multi_h_pTx(water, ref_array):
    assert water.h(p=ref_array['p'], T=ref_array['T'], x=ref_array['x']) == approx(ref_array['h'], abs=1e-1)
def test_multi_h_pd(water, ref_array):
    assert water.h(p=ref_array['p'], d=ref_array['d']) == approx(ref_array['h'], abs=1e-1)
def test_multi_h_pv(water, ref_array):
    assert water.h(p=ref_array['p'], v=ref_array['v']) == approx(ref_array['h'], abs=1e-1)
def test_multi_h_Td(water, ref_array):
    assert water.h(T=ref_array['T'], d=ref_array['d']) == approx(ref_array['h'], abs=1e-1)
def test_multi_h_Tv(water, ref_array):
    assert water.h(T=ref_array['T'], v=ref_array['v']) == approx(ref_array['h'], abs=1e-1)
def test_multi_h_ps(water, ref_array):
    assert water.h(p=ref_array['p'], s=ref_array['s']) == approx(ref_array['h'], abs=1e-1)
def test_multi_h_Ts(water, ref_array):
    assert water.h(T=ref_array['T'], s=ref_array['s']) == approx(ref_array['h'], abs=1e-1)
def test_multi_h_ds(water, ref_array):
    assert water.h(d=ref_array['d'], s=ref_array['s']) == approx(ref_array['h'], abs=1e-1)
def test_multi_h_vs(water, ref_array):
    assert water.h(v=ref_array['v'], s=ref_array['s']) == approx(ref_array['h'], abs=1e-1)
def test_multi_h_pe(water, ref_array):
    assert water.h(p=ref_array['p'], e=ref_array['e']) == approx(ref_array['h'], abs=1e-1)
def test_multi_h_Te(water, ref_array):
    with raises(pm.utility.PMAnalysisError):
        assert water.h(T=ref_array['T'], e=ref_array['e']) == approx(ref_array['h'], abs=1e-1)
def test_multi_h_de(water, ref_array):
    assert water.h(d=ref_array['d'], e=ref_array['e']) == approx(ref_array['h'], abs=1e-1)
def test_multi_h_ve(water, ref_array):
    assert water.h(v=ref_array['v'], e=ref_array['e']) == approx(ref_array['h'], abs=1e-1)


# ##### INTERNAL ENERGY ####################################################################
def test_multi_e_pTx(water, ref_array):
    assert water.e(p=ref_array['p'], T=ref_array['T'], x=ref_array['x']) == approx(ref_array['e'], abs=1e-1)
def test_multi_e_pd(water, ref_array):
    assert water.e(p=ref_array['p'], d=ref_array['d']) == approx(ref_array['e'], abs=1e-1)
def test_multi_e_pv(water, ref_array):
    assert water.e(p=ref_array['p'], v=ref_array['v']) == approx(ref_array['e'], abs=1e-1)
def test_multi_e_Td(water, ref_array):
    assert water.e(T=ref_array['T'], d=ref_array['d']) == approx(ref_array['e'], abs=1e-1)
def test_multi_e_Tv(water, ref_array):
    assert water.e(T=ref_array['T'], v=ref_array['v']) == approx(ref_array['e'], abs=1e-1)
def test_multi_e_ps(water, ref_array):
    assert water.e(p=ref_array['p'], s=ref_array['s']) == approx(ref_array['e'], abs=1e-1)
def test_multi_e_Ts(water, ref_array):
    assert water.e(T=ref_array['T'], s=ref_array['s']) == approx(ref_array['e'], abs=1e-1)
def test_multi_e_ds(water, ref_array):
    assert water.e(d=ref_array['d'], s=ref_array['s']) == approx(ref_array['e'], abs=1e-1)
def test_multi_e_vs(water, ref_array):
    assert water.e(v=ref_array['v'], s=ref_array['s']) == approx(ref_array['e'], abs=1e-1)
def test_multi_e_ph(water, ref_array):
    assert water.e(p=ref_array['p'], h=ref_array['h']) == approx(ref_array['e'], abs=1e-1)
def test_multi_e_Th(water, ref_array):
    with raises(pm.utility.PMAnalysisError):
        assert water.e(T=ref_array['T'], h=ref_array['h']) == approx(ref_array['e'], abs=1e-1)
def test_multi_e_dh(water, ref_array):
    assert water.e(d=ref_array['d'], h=ref_array['h']) == approx(ref_array['e'], abs=1e-1)
def test_multi_e_vh(water, ref_array):
    assert water.e(v=ref_array['v'], h=ref_array['h']) == approx(ref_array['e'], abs=1e-1)



# ##### KNOWN FORMER BUG CONDITIONS #######################################################


def test_d_s_crash(water):
    # A documented crash condition that reaches s_ and I with different sizes
    # within d_s
    s = 8
    T = [275, 335]
    d = water.d_s(T=T, s=s)
    assert not np.isnan(d).any()


def test_multipoint_nan_failure(water):
    h = [4372, 4373]  # kJ/kg
    p = 1000  # bar
    ans = water.T_h(h=h, p=p)  # expect [1272.80748982, np.nan] as T exceeds Tmax
    assert not np.isnan(ans[0])
    assert np.isnan(ans[1])


def test_error_condition(water):
    # A documented crash condition that reaches parts of hsd with different
    # sizes of I.
    water.hsd(p=[1.611654e-2, 2.8095e-2], d=[999.056])
