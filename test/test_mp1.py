import pyromat as pm
import numpy as np
from pytest import approx, fixture


@fixture
def water():
    w = pm.get('mp.H2O')
    return w


@fixture
def r134a():
    r = pm.get('mp.C2H2F4')
    return r


# Four reference points for H2O (NIST Webbook), and two saturation states
@fixture
def ref_a():
    # liquid
    vals = {'T': 300,        # K
            'p': 5,          # bar
            'd': 996.74,     # kg/m3
            'v': 0.0010033,  # m3/kg
            'e': 112.52,     # kJ/kg
            'h': 113.02,     # kJ/kg
            's': 0.39295,    # kJ/kg K
            }
    return vals


@fixture
def ref_b():
    # vapor
    vals = {'T': 600,        # K
            'p': 5,          # bar
            'd': 1.8242,     # kg/m3
            'v': 0.54820,    # m3/kg
            'e': 2846.0,     # kJ/kg
            'h': 3120.1,     # kJ/kg
            's': 7.5561,     # kJ/kg K
            }
    return vals


@fixture
def ref_c():
    # mixture
    vals = {'T': 424.98,     # K
            'p': 5,          # bar
            'd': 5.321,      # kg/m3
            'v': 0.18795,    # m3/kg
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
            'v': 0.012029,   # m3/kg
            'e': 2961.5,     # kJ/kg
            'h': 3262.2,     # kJ/kg
            's': 6.0867,     # kJ/kg K
            }
    return vals


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
    # Verify no errors  # TODO run built-in test
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

def test_ref_a_p(water, ref_a):
    assert water.p(T=ref_a['T'], d=ref_a['d']) == approx(ref_a['p'], abs=1e-1)


def test_ref_b_p(water, ref_b):
    assert water.p(T=ref_b['T'], d=ref_b['d']) == approx(ref_b['p'], abs=1e-1)


def test_ref_c_p(water, ref_c):
    assert water.p(T=ref_c['T'], d=ref_c['d']) == approx(ref_c['p'], abs=1e-1)
    assert water.p(T=ref_c['T'], d=ref_c['d'], quality=True) == approx((ref_c['p'], ref_c['x']), abs=1e-1)
    assert water.p(T=ref_c['T'], x=ref_c['x']) == approx(ref_c['p'], abs=1e-1)


def test_ref_d_p(water, ref_d):
    assert water.p(T=ref_d['T'], d=ref_d['d']) == approx(ref_d['p'], abs=1e-1)


def test_ref_a_T(water, ref_a):
    assert water.T(p=ref_a['p'], d=ref_a['d']) == approx(ref_a['T'], abs=1e-1)


def test_ref_b_T(water, ref_b):
    assert water.T(p=ref_b['p'], d=ref_b['d']) == approx(ref_b['T'], abs=1e-1)


def test_ref_c_T(water, ref_c):
    assert water.T(p=ref_c['p'], d=ref_c['d']) == approx(ref_c['T'], abs=1e-1)
    assert water.T(p=ref_c['p'], x=ref_c['x']) == approx(ref_c['T'], abs=1e-1)


def test_ref_d_T(water, ref_d):
    assert water.T(p=ref_d['p'], d=ref_d['d']) == approx(ref_d['T'], abs=1e-1)


def test_ref_a_T_s(water, ref_a):
    assert water.T_s(p=ref_a['p'], s=ref_a['s']) == approx(ref_a['T'], abs=1e-1)
    assert water.T_s(s=ref_a['s'], d=ref_a['d']) == approx(ref_a['T'], abs=1e-1)


def test_ref_b_T_s(water, ref_b):
    assert water.T_s(p=ref_b['p'], s=ref_b['s']) == approx(ref_b['T'], abs=1e-1)
    assert water.T_s(s=ref_b['s'], d=ref_b['d']) == approx(ref_b['T'], abs=1e-1)


def test_ref_c_T_s(water, ref_c):
    assert water.T_s(p=ref_c['p'], s=ref_c['s']) == approx(ref_c['T'], abs=1e-1)
    assert water.T_s(p=ref_c['p'], s=ref_c['s'], quality=True) == approx((ref_c['T'], ref_c['x']), abs=1e-1)
    assert water.T_s(s=ref_c['s'], d=ref_c['d']) == approx(ref_c['T'], abs=1e-1)


def test_ref_d_T_s(water, ref_d):
    assert water.T_s(p=ref_d['p'], s=ref_d['s']) == approx(ref_d['T'], abs=1e-1)
    assert water.T_s(s=ref_d['s'], d=ref_d['d']) == approx(ref_d['T'], abs=1e-1)


def test_ref_a_T_h(water, ref_a):
    assert water.T_h(p=ref_a['p'], h=ref_a['h']) == approx(ref_a['T'], abs=1e-1)
    assert water.T_h(h=ref_a['h'], d=ref_a['d']) == approx(ref_a['T'], abs=1e-1)


def test_ref_b_T_h(water, ref_b):
    assert water.T_h(p=ref_b['p'], h=ref_b['h']) == approx(ref_b['T'], abs=1e-1)
    assert water.T_h(h=ref_b['h'], d=ref_b['d']) == approx(ref_b['T'], abs=1e-1)


def test_ref_c_T_h(water, ref_c):
    assert water.T_h(p=ref_c['p'], h=ref_c['h']) == approx(ref_c['T'], abs=1e-1)
    assert water.T_h(p=ref_c['p'], h=ref_c['h'], quality=True) == approx((ref_c['T'], ref_c['x']), abs=1e-1)
    assert water.T_h(h=ref_c['h'], d=ref_c['d']) == approx(ref_c['T'], abs=1e-1)


def test_ref_d_T_h(water, ref_d):
    assert water.T_h(p=ref_d['p'], h=ref_d['h']) == approx(ref_d['T'], abs=1e-1)
    assert water.T_h(h=ref_d['h'], d=ref_d['d']) == approx(ref_d['T'], abs=1e-1)


def test_ref_a_d(water, ref_a):
    assert water.d(p=ref_a['p'], T=ref_a['T']) == approx(ref_a['d'], abs=1e-2)


def test_ref_b_d(water, ref_b):
    assert water.d(p=ref_b['p'], T=ref_b['T']) == approx(ref_b['d'], abs=1e-2)


def test_ref_c_d(water, ref_c):
    assert water.d(p=ref_c['p'], x=ref_c['x']) == approx(ref_c['d'], abs=1e-2)
    assert water.d(T=ref_c['T'], x=ref_c['x']) == approx(ref_c['d'], abs=1e-2)


def test_ref_d_d(water, ref_d):
    assert water.d(p=ref_d['p'], T=ref_d['T']) == approx(ref_d['d'], abs=1e-2)


def test_ref_a_d_s(water, ref_a):
    assert water.d_s(s=ref_a['s'], T=ref_a['T']) == approx(ref_a['d'], abs=1e-2)


def test_ref_b_d_s(water, ref_b):
    assert water.d_s(s=ref_b['s'], T=ref_b['T']) == approx(ref_b['d'], abs=1e-2)


def test_ref_c_d_s(water, ref_c):
    assert water.d_s(s=ref_c['s'], T=ref_c['T']) == approx(ref_c['d'], abs=1e-2)


def test_ref_d_d_s(water, ref_d):
    assert water.d_s(s=ref_d['s'], T=ref_d['T']) == approx(ref_d['d'], abs=1e-2)


def test_ref_a_s(water, ref_a):
    assert water.s(p=ref_a['p'], T=ref_a['T']) == approx(ref_a['s'], abs=1e-4)
    assert water.s(p=ref_a['p'], d=ref_a['d']) == approx(ref_a['s'], abs=1e-3)
    assert water.s(T=ref_a['T'], d=ref_a['d']) == approx(ref_a['s'], abs=1e-3)


def test_ref_b_s(water, ref_b):
    assert water.s(p=ref_b['p'], T=ref_b['T']) == approx(ref_b['s'], abs=1e-4)
    assert water.s(p=ref_b['p'], d=ref_b['d']) == approx(ref_b['s'], abs=1e-3)
    assert water.s(T=ref_b['T'], d=ref_b['d']) == approx(ref_b['s'], abs=1e-3)


def test_ref_c_s(water, ref_c):
    assert water.s(p=ref_c['p'], x=ref_c['x']) == approx(ref_c['s'], abs=1e-4)
    assert water.s(T=ref_c['T'], x=ref_c['x']) == approx(ref_c['s'], abs=1e-4)
    assert water.s(p=ref_c['p'], d=ref_c['d']) == approx(ref_c['s'], abs=1e-3)
    assert water.s(T=ref_c['T'], d=ref_c['d']) == approx(ref_c['s'], abs=1e-3)
    assert water.s(T=ref_c['T'], d=ref_c['d'], quality=True) == approx((ref_c['s'], ref_c['x']), abs=1e-3)


def test_ref_d_s(water, ref_d):
    assert water.s(p=ref_d['p'], T=ref_d['T']) == approx(ref_d['s'], abs=1e-4)
    assert water.s(p=ref_d['p'], d=ref_d['d']) == approx(ref_d['s'], abs=1e-3)
    assert water.s(T=ref_d['T'], d=ref_d['d']) == approx(ref_d['s'], abs=1e-3)


def test_ref_a_h(water, ref_a):
    assert water.h(p=ref_a['p'], T=ref_a['T']) == approx(ref_a['h'], abs=1e-1)
    assert water.h(p=ref_a['p'], d=ref_a['d']) == approx(ref_a['h'], abs=1e-1)
    assert water.h(T=ref_a['T'], d=ref_a['d']) == approx(ref_a['h'], abs=1e-1)


def test_ref_b_h(water, ref_b):
    assert water.h(p=ref_b['p'], T=ref_b['T']) == approx(ref_b['h'], abs=1e-1)
    assert water.h(p=ref_b['p'], d=ref_b['d']) == approx(ref_b['h'], abs=1e-1)
    assert water.h(T=ref_b['T'], d=ref_b['d']) == approx(ref_b['h'], abs=1e-1)


def test_ref_c_h(water, ref_c):
    assert water.h(p=ref_c['p'], x=ref_c['x']) == approx(ref_c['h'], abs=1e-0)
    assert water.h(T=ref_c['T'], x=ref_c['x']) == approx(ref_c['h'], abs=1e-0)
    assert water.h(p=ref_c['p'], d=ref_c['d']) == approx(ref_c['h'], abs=1e-0)
    assert water.h(T=ref_c['T'], d=ref_c['d']) == approx(ref_c['h'], abs=1e-0)
    assert water.h(T=ref_c['T'], d=ref_c['d'], quality=True) == approx((ref_c['h'], ref_c['x']), abs=1e-0)


def test_ref_d_h(water, ref_d):
    assert water.h(p=ref_d['p'], T=ref_d['T']) == approx(ref_d['h'], abs=1e-1)
    assert water.h(p=ref_d['p'], d=ref_d['d']) == approx(ref_d['h'], abs=1e-1)
    assert water.h(T=ref_d['T'], d=ref_d['d']) == approx(ref_d['h'], abs=1e-1)


def test_ref_a_e(water, ref_a):
    assert water.e(p=ref_a['p'], T=ref_a['T']) == approx(ref_a['e'], abs=1e-1)
    assert water.e(p=ref_a['p'], d=ref_a['d']) == approx(ref_a['e'], abs=1e-1)
    assert water.e(T=ref_a['T'], d=ref_a['d']) == approx(ref_a['e'], abs=1e-1)


def test_ref_b_e(water, ref_b):
    assert water.e(p=ref_b['p'], T=ref_b['T']) == approx(ref_b['e'], abs=1e-1)
    assert water.e(p=ref_b['p'], d=ref_b['d']) == approx(ref_b['e'], abs=1e-1)
    assert water.e(T=ref_b['T'], d=ref_b['d']) == approx(ref_b['e'], abs=1e-1)


def test_ref_c_e(water, ref_c):
    assert water.e(p=ref_c['p'], x=ref_c['x']) == approx(ref_c['e'], abs=1e-0)
    assert water.e(T=ref_c['T'], x=ref_c['x']) == approx(ref_c['e'], abs=1e-0)
    assert water.e(p=ref_c['p'], d=ref_c['d']) == approx(ref_c['e'], abs=1e-0)
    assert water.e(T=ref_c['T'], d=ref_c['d']) == approx(ref_c['e'], abs=1e-0)
    assert water.e(T=ref_c['T'], d=ref_c['d'], quality=True) == approx((ref_c['e'], ref_c['x']), abs=1e-0)


def test_ref_d_e(water, ref_d):
    assert water.e(p=ref_d['p'], T=ref_d['T']) == approx(ref_d['e'], abs=1e-1)
    assert water.e(p=ref_d['p'], d=ref_d['d']) == approx(ref_d['e'], abs=1e-1)
    assert water.e(T=ref_d['T'], d=ref_d['d']) == approx(ref_d['e'], abs=1e-1)


# ##### Test combinations and different data formats ##### #

def test_multi_P(water, ref_a, ref_b, ref_c, ref_d):
    T = [ref_a['T'], ref_b['T'], ref_c['T'], ref_d['T']]
    p = [ref_a['p'], ref_b['p'], ref_c['p'], ref_d['p']]
    d = [ref_a['d'], ref_b['d'], ref_c['d'], ref_d['d']]
    s = [ref_a['s'], ref_b['s'], ref_c['s'], ref_d['s']]
    h = [ref_a['h'], ref_b['h'], ref_c['h'], ref_d['h']]
    e = [ref_a['e'], ref_b['e'], ref_c['e'], ref_d['e']]
    assert water.p(T=T, d=d) == approx(p, abs=1e-1)
    assert water.p(T=np.array(T), d=np.array(d)) == approx(p, abs=1e-1)


def test_multi_T(water, ref_a, ref_b, ref_c, ref_d):
    T = [ref_a['T'], ref_b['T'], ref_c['T'], ref_d['T']]
    p = [ref_a['p'], ref_b['p'], ref_c['p'], ref_d['p']]
    d = [ref_a['d'], ref_b['d'], ref_c['d'], ref_d['d']]
    s = [ref_a['s'], ref_b['s'], ref_c['s'], ref_d['s']]
    h = [ref_a['h'], ref_b['h'], ref_c['h'], ref_d['h']]
    e = [ref_a['e'], ref_b['e'], ref_c['e'], ref_d['e']]
    assert water.T(p=p, d=d) == approx(T, abs=1e-1)
    assert water.T(p=np.array(p), d=np.array(d)) == approx(T, abs=1e-1)
    assert water.T_s(p=p, s=s) == approx(T, abs=1e-1)
    assert water.T_s(s=s, d=d) == approx(T, abs=1e-1)
    assert water.T_h(h=h, d=d) == approx(T, abs=1e-1)
    assert water.T_h(h=h, p=p) == approx(T, abs=1e-1)


def test_multi_d(water, ref_a, ref_b, ref_c, ref_d):
    T = [ref_a['T'], ref_b['T'], ref_c['T'], ref_d['T']]
    p = [ref_a['p'], ref_b['p'], ref_c['p'], ref_d['p']]
    d = [ref_a['d'], ref_b['d'], ref_c['d'], ref_d['d']]
    s = [ref_a['s'], ref_b['s'], ref_c['s'], ref_d['s']]
    h = [ref_a['h'], ref_b['h'], ref_c['h'], ref_d['h']]
    e = [ref_a['e'], ref_b['e'], ref_c['e'], ref_d['e']]
    x = [-1, -1, ref_c['x'], -1]
    assert water.d(p=p, T=T, x=x) == approx(d, abs=1e-1)
    assert water.d(p=np.array(p), T=np.array(T), x=np.array(x)) == approx(d, abs=1e-1)
    assert water.d_s(s=s, T=T) == approx(d, abs=1e-1)


def test_multi_s(water, ref_a, ref_b, ref_c, ref_d):
    T = [ref_a['T'], ref_b['T'], ref_c['T'], ref_d['T']]
    p = [ref_a['p'], ref_b['p'], ref_c['p'], ref_d['p']]
    d = [ref_a['d'], ref_b['d'], ref_c['d'], ref_d['d']]
    s = [ref_a['s'], ref_b['s'], ref_c['s'], ref_d['s']]
    h = [ref_a['h'], ref_b['h'], ref_c['h'], ref_d['h']]
    e = [ref_a['e'], ref_b['e'], ref_c['e'], ref_d['e']]
    x = [-1, -1, ref_c['x'], -1]
    assert water.s(p=p, T=T, x=x) == approx(s, abs=1e-1)
    assert water.s(p=p, d=d) == approx(s, abs=1e-1)
    assert water.s(T=T, d=d) == approx(s, abs=1e-1)
    assert water.s(p=np.array(p), T=np.array(T), x=np.array(x)) == approx(s, abs=1e-1)


def test_multi_h(water, ref_a, ref_b, ref_c, ref_d):
    T = [ref_a['T'], ref_b['T'], ref_c['T'], ref_d['T']]
    p = [ref_a['p'], ref_b['p'], ref_c['p'], ref_d['p']]
    d = [ref_a['d'], ref_b['d'], ref_c['d'], ref_d['d']]
    s = [ref_a['s'], ref_b['s'], ref_c['s'], ref_d['s']]
    h = [ref_a['h'], ref_b['h'], ref_c['h'], ref_d['h']]
    e = [ref_a['e'], ref_b['e'], ref_c['e'], ref_d['e']]
    x = [-1, -1, ref_c['x'], -1]
    assert water.h(p=p, T=T, x=x) == approx(h, abs=1e-0)
    assert water.h(p=p, d=d) == approx(h, abs=1e-0)
    assert water.h(T=T, d=d) == approx(h, abs=1e-0)
    assert water.h(p=np.array(p), T=np.array(T), x=np.array(x)) == approx(h, abs=1e-0)


def test_multi_e(water, ref_a, ref_b, ref_c, ref_d):
    T = [ref_a['T'], ref_b['T'], ref_c['T'], ref_d['T']]
    p = [ref_a['p'], ref_b['p'], ref_c['p'], ref_d['p']]
    d = [ref_a['d'], ref_b['d'], ref_c['d'], ref_d['d']]
    s = [ref_a['s'], ref_b['s'], ref_c['s'], ref_d['s']]
    h = [ref_a['h'], ref_b['h'], ref_c['h'], ref_d['h']]
    e = [ref_a['e'], ref_b['e'], ref_c['e'], ref_d['e']]
    x = [-1, -1, ref_c['x'], -1]
    assert water.e(p=p, T=T, x=x) == approx(e, abs=1e-0)
    assert water.e(p=p, d=d) == approx(e, abs=1e-0)
    assert water.e(T=T, d=d) == approx(e, abs=1e-0)
    assert water.e(p=np.array(p), T=np.array(T), x=np.array(x)) == approx(e, abs=1e-0)


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
