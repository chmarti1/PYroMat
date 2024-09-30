import pyromat as pm
import numpy as np
from pytest import approx, raises
import pytest


class TestInputErrors:
    @pytest.fixture(params=('mp.H2O', 'mp.C2H2F4'),
                    ids=('water', 'R134a'))
    def subst(self, request):
        o = pm.get(request.param)
        return o

    def test_T_only_varg(self, subst):
        # Based on how .T is written, it interprets this as a temperature given
        # So just make sure that works
        assert subst.T(300) == subst.T(T=300, p=pm.config['def_p'])

    def test_T_double_specified(self, subst):
        with raises(pm.utility.PMParamError):
            subst.T(300, T=300)

    def test_p_double_specified(self, subst):
        with raises(pm.utility.PMParamError):
            subst.T(300, 1, p=1)

    def test_varg_specified(self, subst):
        assert subst.s(300, 1) == subst.s(p=1, T=300)

    def test_varg_toomany(self, subst):
        with raises(pm.utility.PMParamError):
            subst.T(300, 1, 5)

    def test_any_arg_toomany(self, subst):
        with raises(pm.utility.PMParamError):
            subst.T(300, 1, s=5)
        with raises(pm.utility.PMParamError):
            subst.T(T=300, p=1, s=5)

    def test_default_usage(self, subst):
        assert subst.s() == subst.s(T=pm.config['def_T'], p=pm.config['def_p'])
        assert subst.s(T=pm.config['def_T']) == subst.s(T=pm.config['def_T'], p=pm.config['def_p'])
        assert subst.s(p=pm.config['def_p']) == subst.s(T=pm.config['def_T'], p=pm.config['def_p'])
        # And for saturated
        assert subst.ss() == subst.ss(T=pm.config['def_T'])
        assert subst.Ts() == subst.Ts(p=pm.config['def_p'])
        assert subst.ps() == subst.ps(T=pm.config['def_T'])

    def test_invarg_toomany(self, subst):
        with raises(pm.utility.PMParamError):
            subst.T(h=300, s=5)
        with raises(pm.utility.PMParamError):
            subst.T(e=300, s=5)

    def test_d_v_collision(self, subst):
        with raises(pm.utility.PMParamError):
            subst.T(d=1, v=1)

    def test_illegal_prop(self, subst):
        with raises(pm.utility.PMParamError):
            subst.T(T=300, fakeprop=10000)

    def test_T_oob_all(self, subst):
        with raises(pm.utility.PMParamError):
            subst.T(T=30000, p=1)

    def test_T_oob_some(self, subst):
        vals = subst.T(T=[300, 30000], p=1)
        assert vals[0] == 300
        assert np.isnan(vals[1])

    def test_p_oob_all(self, subst):
        with raises(pm.utility.PMParamError):
            subst.T(T=300, p=1e100)

    def test_p_oob_some(self, subst):
        vals = subst.p(T=300, p=[1, 1e100])
        assert vals[0] == approx(1)
        assert np.isnan(vals[1])

    def test_x_oob(self, subst):
        with raises(pm.utility.PMParamError):
            subst.T(p=1, x=2)
        with raises(pm.utility.PMParamError):
            subst.T(p=1, x=-2)

    def test_sat_oob_T(self, subst):
        Tt, pt = subst.triple()
        Tc, pc = subst.critical()
        with raises(pm.utility.PMParamError):
            subst.ps(T=Tt-1)
        with raises(pm.utility.PMParamError):
            subst.ps(T=Tc+1)
        with raises(pm.utility.PMParamError):
            subst.ds(T=Tt - 1)
        with raises(pm.utility.PMParamError):
            subst.ds(T=Tc + 1)

    def test_sat_oob_p(self, subst):
        Tt, pt = subst.triple()
        Tc, pc = subst.critical()
        with raises(pm.utility.PMParamError):
            subst.Ts(p=pt-1)
        with raises(pm.utility.PMParamError):
            subst.Ts(p=pc+1)
        with raises(pm.utility.PMParamError):
            subst.ds(p=pt-1)
        with raises(pm.utility.PMParamError):
            subst.ds(p=pc+1)

    def test_sat_toomany(self, subst):
        with raises(pm.utility.PMParamError):
            subst.ds(T=300, p=1)


class TestGeneral:

    @pytest.fixture
    def water(self):
        return pm.get('mp.H2O')

    @pytest.mark.parametrize('sub', ('mp.H2O', 'mp.C2H2F4'), ids=('water', 'r134a'))
    def test_get(self, sub):
        # Just confirm no error
        mp1obj = pm.get(sub)

    def test_Tlim(self, water):
        lim = (273.16, 1273)
        assert water.Tlim() == approx(lim)
        assert water.Tlim(p=1) == approx(lim)

    def test_plim(self, water):
        lim = (0, 10000)
        assert water.plim() == approx(lim)
        assert water.plim(T=1000) == approx(lim)

    def test_critical(self, water):
        cr = (647.096, 220.64, 322.0)
        assert water.critical() == approx(cr[0:2])
        assert water.critical(density=True) == approx(cr)

    def test_triple(self, water):
        tp = (273.16, 0.00611655)
        assert water.triple() == approx(tp)

    @pytest.mark.parametrize('sub, expect',
                             [('mp.H2O', 18.015268),
                              ('mp.C2H2F4', 102.032)])
    def test_mw(self, sub, expect):
        assert pm.get(sub).mw() == approx(expect)

    @pytest.mark.parametrize('sub, expect',
                             [('mp.H2O', pm.units.const_Ru / 18.015268),
                              ('mp.C2H2F4', pm.units.const_Ru / 102.032)])
    def test_R(self, sub, expect):
        assert pm.get(sub).R() == approx(expect)

    @pytest.mark.parametrize('sub, expect', [
        ('mp.H2O', {'H': 2, 'O': 1}),
        ('mp.C2H2F4', {'C': 2, 'H': 2, 'F': 4})
    ])
    def test_atoms(self, sub, expect):
        assert pm.get(sub).atoms() == approx(expect)

    @pytest.mark.parametrize('sub, expect', [
        ('mp.H2O', 4.1813595),
        ('mp.C2H2F4', 0.85124365)
    ])
    def test_cp(self, sub, expect):
        assert pm.get(sub).cp() == approx(expect)
        assert pm.get(sub).cp(x=0.5) == approx(np.inf)
        cp, x = pm.get(sub).cp(quality=True)
        assert cp == approx(expect)

    @pytest.mark.parametrize('sub, expect, expect05', [
        ('mp.H2O', 4.13760893, 67.6735),
        ('mp.C2H2F4', 0.7603562, 3.39377)
    ])
    def test_cv(self, sub, expect, expect05):
        assert pm.get(sub).cv() == approx(expect)
        assert pm.get(sub).cv(x=0.5) == approx(expect05)
        cv, x = pm.get(sub).cv(quality=True)
        assert cv == approx(expect)

    @pytest.mark.parametrize('sub, expect', [
        ('mp.H2O', 1.01057388),
        ('mp.C2H2F4', 1.11953272)
    ])
    def test_gam(self, sub, expect):
        assert pm.get(sub).gam() == approx(expect)
        assert pm.get(sub).gam(x=0.5) == approx(np.inf)
        gam, x = pm.get(sub).gam(quality=True)
        assert gam == approx(expect)


class TestSat:

    @pytest.fixture
    def water(self):
        return pm.get('mp.H2O')

    @pytest.fixture
    def ref_sat_T(self):  # Some slight deviation from NIST. Keeping both vals listed
        vals = {'T': 500,  # K
                'p': 26.392226747183663,  # bar # NIST 26.392
                'd': (831.35713738, 13.19861192),  # kg/m3  # NIST 831.31, 13.199
                'v': (0.00120285, 0.07576554),  # m3/kg  # NIST 0.0012029, 0.07564
                'e': (972.22040467, 2602.55303444),  # kJ/kg # NIST 972.26, 2602.5
                'h': (975.45086017, 2802.5115447),  # kJ/kg  # NIST 975.43, 2.5810
                's': (2.58098058, 6.23522113),  # kJ/kg K  # NIST 2.5810, 6.2351
                }
        return vals

    @pytest.fixture
    def ref_sat_p(self):
        vals = {'T': 424.98,  # K
                'p': 5.,  # bar
                'd': (915.29133902, 2.66790101),  # kg/m3  # NIST 915.29, 2.6680
                'v': (0.00109255, 0.3748265),  # m3/kg  # NIST 0.0010925, 0.37481
                'e': (639.54607876, 2560.73598375),  # kJ/kg  # NIST 639.54, 2560.7
                'h': (640.09563612, 2748.14164556),  # kJ/kg  # NIST 640.09, 2748.1
                's': (1.86040396, 6.82076331),  # kJ/kg K  # NIST 1.8604, 6.8207
                }
        return vals

    # ##### Test accuracy to reference saturation data ##### #

    def test_ref_sat_p_Ts(self, water, ref_sat_p):
        assert water.Ts(p=ref_sat_p['p']) == approx(ref_sat_p['T'], rel=1e-5, abs=1e-2)
        assert water.Ts(p=np.tile(ref_sat_p['p'], 3)) == approx(
            np.tile(ref_sat_p['T'], 3), rel=1e-5, abs=1e-2)

    def test_ref_sat_p_ds(self, water, ref_sat_p):
        assert water.ds(p=ref_sat_p['p']) == approx(ref_sat_p['d'], rel=1e-5, abs=1e-2)
        assert water.ds(p=np.tile(ref_sat_p['p'], 3))[0] == approx(
            np.tile(ref_sat_p['d'][0], 3), rel=1e-5, abs=1e-2)
        assert water.ds(p=np.tile(ref_sat_p['p'], 3))[1] == approx(
            np.tile(ref_sat_p['d'][1], 3), rel=1e-5, abs=1e-2)

    def test_ref_sat_p_vs(self, water, ref_sat_p):
        assert water.vs(p=ref_sat_p['p']) == approx(ref_sat_p['v'], rel=1e-5, abs=1e-2)
        assert water.vs(p=np.tile(ref_sat_p['p'], 3))[0] == approx(
            np.tile(ref_sat_p['v'][0], 3), rel=1e-5, abs=1e-2)
        assert water.vs(p=np.tile(ref_sat_p['p'], 3))[1] == approx(
            np.tile(ref_sat_p['v'][1], 3), rel=1e-5, abs=1e-2)

    def test_ref_sat_p_es(self, water, ref_sat_p):
        assert water.es(p=ref_sat_p['p']) == approx(ref_sat_p['e'], rel=1e-5, abs=1e-2)
        assert water.es(p=np.tile(ref_sat_p['p'], 3))[0] == approx(
            np.tile(ref_sat_p['e'][0], 3), rel=1e-5, abs=1e-2)
        assert water.es(p=np.tile(ref_sat_p['p'], 3))[1] == approx(
            np.tile(ref_sat_p['e'][1], 3), rel=1e-5, abs=1e-2)

    def test_ref_sat_p_hs(self, water, ref_sat_p):
        assert water.hs(p=ref_sat_p['p']) == approx(ref_sat_p['h'], rel=1e-5, abs=1e-2)
        assert water.hs(p=np.tile(ref_sat_p['p'], 3))[0] == approx(
            np.tile(ref_sat_p['h'][0], 3), rel=1e-5, abs=1e-2)
        assert water.hs(p=np.tile(ref_sat_p['p'], 3))[1] == approx(
            np.tile(ref_sat_p['h'][1], 3), rel=1e-5, abs=1e-2)

    def test_ref_sat_p_ss(self, water, ref_sat_p):
        assert water.ss(p=ref_sat_p['p']) == approx(ref_sat_p['s'], rel=1e-5, abs=1e-2)
        assert water.ss(p=np.tile(ref_sat_p['p'], 3))[0] == approx(
            np.tile(ref_sat_p['s'][0], 3), rel=1e-5, abs=1e-2)
        assert water.ss(p=np.tile(ref_sat_p['p'], 3))[1] == approx(
            np.tile(ref_sat_p['s'][1], 3), rel=1e-5, abs=1e-2)

    def test_ref_sat_T_ps(self, water, ref_sat_T):
        assert water.ps(T=ref_sat_T['T']) == approx(ref_sat_T['p'], rel=1e-5, abs=1e-2)
        assert water.ps(T=np.tile(ref_sat_T['T'], 3)) == approx(
            np.tile(ref_sat_T['p'], 3), rel=1e-5, abs=1e-2)

    def test_ref_sat_T_ds(self, water, ref_sat_T):
        assert water.ds(T=ref_sat_T['T']) == approx(ref_sat_T['d'], rel=1e-5, abs=1e-2)
        assert water.ds(T=np.tile(ref_sat_T['T'], 3))[0] == approx(
            np.tile(ref_sat_T['d'][0], 3), rel=1e-5, abs=1e-2)
        assert water.ds(T=np.tile(ref_sat_T['T'], 3))[1] == approx(
            np.tile(ref_sat_T['d'][1], 3), rel=1e-5, abs=1e-2)

    def test_ref_sat_T_es(self, water, ref_sat_T):
        assert water.es(T=ref_sat_T['T']) == approx(ref_sat_T['e'], rel=1e-5, abs=1e-2)
        assert water.es(T=np.tile(ref_sat_T['T'], 3))[0] == approx(
            np.tile(ref_sat_T['e'][0], 3), rel=1e-5, abs=1e-2)
        assert water.es(T=np.tile(ref_sat_T['T'], 3))[1] == approx(
            np.tile(ref_sat_T['e'][1], 3), rel=1e-5, abs=1e-2)

    def test_ref_sat_T_hs(self, water, ref_sat_T):
        assert water.hs(T=ref_sat_T['T']) == approx(ref_sat_T['h'], rel=1e-5, abs=1e-2)
        assert water.hs(T=np.tile(ref_sat_T['T'], 3))[0] == approx(
            np.tile(ref_sat_T['h'][0], 3), rel=1e-5, abs=1e-2)
        assert water.hs(T=np.tile(ref_sat_T['T'], 3))[1] == approx(
            np.tile(ref_sat_T['h'][1], 3), rel=1e-5, abs=1e-2)

    def test_ref_sat_T_ss(self, water, ref_sat_T):
        assert water.ss(T=ref_sat_T['T']) == approx(ref_sat_T['s'], rel=1e-5, abs=1e-2)
        assert water.ss(T=np.tile(ref_sat_T['T'], 3))[0] == approx(
            np.tile(ref_sat_T['s'][0], 3), rel=1e-5, abs=1e-2)
        assert water.ss(T=np.tile(ref_sat_T['T'], 3))[1] == approx(
            np.tile(ref_sat_T['s'][1], 3), rel=1e-5, abs=1e-2)


def make_props_array(refs, subid):
    """
    Take a property array formatted like refs{}
    """
    propkeys = None
    myrefs = []

    # Find which ref states match this substance
    for refkey in refs:
        if refs[refkey]['sub'] == subid:
            myrefs.append(refs[refkey])
            if propkeys is None:
                propkeys = refs[refkey]['props']

    # Build a new dict from the refs
    newdict= {
        'props': {
            key: [ref['props'][key] for ref in myrefs] for key in propkeys
        },
        'sub': subid,
        'comment': 'multi'
    }
    return newdict


# Four reference points for H2O (NIST Webbook), and two saturation states
# Note that spec volumes are computed as 1/d for accuracy. NIST vals listed.
refs = {}
refs["water_a"] = {
    'props': {
        # values generated by pyromat using T&p. Checks with NIST
        'T': 300,  # K
        'p': 5,  # bar
        'd': 996.73584742,  # kg/m3
        'v': 1/996.73584742,  # m3/kg
        'e': 112.5214752,  # kJ/kg
        'h': 113.02311262,  # kJ/kg
        's': 0.39295625,  # kJ/kg K
        'x': -1,
        'a': 1502.20206096,  # m/s
    },
    'NIST': {
        'T': 300,  # K
        'p': 5,  # bar
        'd': 996.74,  # kg/m3
        'v': 0.0010033,  # m3/kg
        'e': 112.52,  # kJ/kg
        'h': 113.02,  # kJ/kg
        's': 0.39295,  # kJ/kg K
        'x': -1,
        'a': 1502.2,  # m/s
    },
    'sub': 'mp.H2O',
    'comment': 'compressed liquid'
}
refs["water_b"] = {
    'props': {
        # values generated by pyromat using T&p. Checks with NIST
        'T': 600,  # K
        'p': 5,  # bar
        'd': 1.82413427,  # kg/m3
        'v': 1/1.82413427,  # m3/kg
        'e': 2846.01833204,  # kJ/kg
        'h': 3120.12094234,  # kJ/kg
        's': 7.55619253,  # kJ/kg K
        'x': -1,
        'a': 595.97412567,  # m/s
    },
    'NIST': {
        'T': 600,  # K
        'p': 5,  # bar
        'd': 1.8242,  # kg/m3
        'v': 0.54820,  # m3/kg
        'e': 2846.0,  # kJ/kg
        'h': 3120.1,  # kJ/kg
        's': 7.5561,  # kJ/kg K
        'x': -1,
        'a': 595.97,  # m/s
    },
    'sub': 'mp.H2O',
    'comment': 'superheated vapor'
}
refs["water_c"] = {
    'props': {
        # values generated by pyromat using p&x. Checks with NIST
        'T': 424.98151125,  # K
        'p': 5,  # bar
        'd': 5.32029437,  # kg/m3
        'v': 1/5.32029437,  # m3/kg
        'e': 1600.14103126,  # kJ/kg
        'h': 1694.11864084,  # kJ/kg
        's': 4.34058363,  # kJ/kg K
        'x': 0.5,
        'a': np.nan,  # m/s
    },
    'NIST': {
        'T': 424.98,  # K
        'p': 5,  # bar
        'd': 5.321,  # kg/m3
        'v': 0.18795,  # m3/kg
        'e': 1600.12,  # kJ/kg
        'h': 1694.095,  # kJ/kg
        's': 4.3406,  # kJ/kg K
        'x': 0.5,
        'a': np.nan,  # m/s
    },
    'sub': 'mp.H2O',
    'comment': 'saturated mixture'
}
refs["water_d"] = {
    'props': {
        # values generated by pyromat using p&x. Checks with NIST
        'T': 800,  # K
        'p': 250,  # bar
        'd': 83.13115268,  # kg/m3
        'v': 1/83.13115268,  # m3/kg
        'e': 2961.48937738,  # kJ/kg
        'h': 3262.21899789,  # kJ/kg
        's': 6.08676043,  # kJ/kg K
        'x': -1,
        'a': 626.78786445,  # m/s
    },
    'NIST': {
        'T': 800,  # K
        'p': 250,  # bar
        'd': 83.132,  # kg/m3
        'v': 0.012029,  # m3/kg
        'e': 2961.5,  # kJ/kg
        'h': 3262.2,  # kJ/kg
        's': 6.0867,  # kJ/kg K
        'x': -1,
        'a': 626.78,  # m/s
    },
    'sub': 'mp.H2O',
    'comment': 'supercritical fluid'
}

for ref in refs:
    if 'water' in ref:
        refs[ref]['props']['f'] = refs[ref]['props']['e'] - refs[ref]['props']['T'] * refs[ref]['props']['s']
        refs[ref]['props']['g'] = refs[ref]['props']['h'] - refs[ref]['props']['T'] * refs[ref]['props']['s']
        refs[ref]['NIST']['f'] = refs[ref]['NIST']['e'] - refs[ref]['NIST']['T'] * refs[ref]['NIST']['s']
        refs[ref]['NIST']['g'] = refs[ref]['NIST']['h'] - refs[ref]['NIST']['T'] * refs[ref]['NIST']['s']

refs['water_array'] = make_props_array(refs, 'mp.H2O')

refs["r134a_a"] = {
    'props': {
        # values generated by pyromat using T&p. Checks with NIST
        'T': 200,  # K
        'p': 5,  # bar
        'd': 1511.25811875,  # kg/m3
        'v': 1/1511.25811875,  # m3/kg
        'e': 107.27245971,  # kJ/kg
        'h': 107.60330988,  # kJ/kg
        's': 0.606735,  # kJ/kg K
        'x': -1,
        'a': 969.8363053,  # m/s
    },
    'NIST': {
        'T': 300,  # K
        'p': 5,  # bar
        'd': 1511.3,  # kg/m3
        'v': 0.00066170	,  # m3/kg
        'e': 107.27,  # kJ/kg
        'h': 107.60,  # kJ/kg
        's': 0.60674,  # kJ/kg K
        'x': -1,
        'a': 969.84,  # m/s
    },
    'sub': 'mp.C2H2F4',
    'comment': 'compressed liquid'
}

refs["r134a_b"] = {
    'props': {
        # values generated by pyromat using T&p. Checks with NIST
        'T': 300,  # K
        'p': 5,  # bar
        'd': 22.90875049,  # kg/m3
        'v': 1/22.90875049,  # m3/kg
        'e': 396.33598821,  # kJ/kg
        'h': 418.16170929,  # kJ/kg
        's': 1.75600437,  # kJ/kg K
        'x': -1,
        'a': 150.76818444,  # m/s
    },
    'NIST': {
        'T': 300,  # K
        'p': 5,  # bar
        'd': 22.909,  # kg/m3
        'v': 0.043651,  # m3/kg
        'e': 396.34,  # kJ/kg
        'h': 418.16,  # kJ/kg
        's': 1.7560,  # kJ/kg K
        'x': -1,
        'a': 150.77,  # m/s
    },
    'sub': 'mp.C2H2F4',
    'comment': 'superheated vapor'
}

refs["r134a_c"] = {
    'props': {
        # values generated by pyromat using p&x. Checks with NIST
        'T': 288.8850165,  # K
        'p': 5,  # bar
        'd': 47.67918444,  # kg/m3
        'v': 1/47.67918444,  # m3/kg
        'e': 303.99268139,  # kJ/kg
        'h': 314.49195319,  # kJ/kg
        's': 1.3977908,  # kJ/kg K
        'x': 0.5,
        'a': np.nan,  # m/s
    },
    'NIST': {
        'T': 288.88,  # K
        'p': 5,  # bar
        'd': 47.700,  # kg/m3
        'v': 0.020964,  # m3/kg
        'e': 304.01,  # kJ/kg
        'h': 314.49,  # kJ/kg
        's': 1.3978,  # kJ/kg K
        'x': 0.5,
        'a': np.nan,  # m/s
    },
    'sub': 'mp.C2H2F4',
    'comment': 'saturated mixture'
}

refs["r134a_d"] = {
    'props': {
        # values generated by pyromat using T&p. Checks with NIST
        'T': 400,  # K
        'p': 50,  # bar
        'd': 285.05284058,  # kg/m3
        'v': 1/285.05284058,  # m3/kg
        'e': 439.61684889,  # kJ/kg
        'h': 457.15745307,  # kJ/kg
        's': 1.7310427,  # kJ/kg K
        'x': -1,
        'a': 124.51396721,  # m/s
    },
    'NIST': {
        'T': 400,  # K
        'p': 50,  # bar
        'd': 285.05,  # kg/m3
        'v': 0.0035081,  # m3/kg
        'e': 439.62,  # kJ/kg
        'h': 457.16,  # kJ/kg
        's': 1.7310,  # kJ/kg K
        'x': -1,
        'a': 124.51,  # m/s
    },
    'sub': 'mp.C2H2F4',
    'comment': 'supercritical'
}

for ref in refs:
    if 'r134a' in ref:
        refs[ref]['props']['f'] = refs[ref]['props']['e'] - refs[ref]['props']['T'] * refs[ref]['props']['s']
        refs[ref]['props']['g'] = refs[ref]['props']['h'] - refs[ref]['props']['T'] * refs[ref]['props']['s']
        refs[ref]['NIST']['f'] = refs[ref]['NIST']['e'] - refs[ref]['NIST']['T'] * refs[ref]['NIST']['s']
        refs[ref]['NIST']['g'] = refs[ref]['NIST']['h'] - refs[ref]['NIST']['T'] * refs[ref]['NIST']['s']

refs['r134a_array'] = make_props_array(refs, 'mp.C2H2F4')


class TestRefs:

    @pytest.fixture(params=(refs[key] for key in refs), ids=(refs.keys()))
    def refdat(self, request):
        return {'sub': pm.get(request.param['sub']),
                'data': request.param['props']}

    @pytest.fixture(params=('p', 'T', 'd', 'v', 'e', 'h', 's', 'x', 'f', 'g', 'a'))
    def param(self, request):
        return request.param

    def test_pT(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        if (np.array(ref['x']) > 0).any():
            pass
        else:
            assert fn(p=ref['p'], T=ref['T']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_pTx(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(p=ref['p'], T=ref['T'], x=ref['x']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_pd(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(p=ref['p'], d=ref['d']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_pd_qualityflag(self, param, refdat):
        if param == 'x':
            pass  # quality flag not supported for substance.x()
        else:
            sub, ref = refdat['sub'], refdat['data']
            fn = getattr(sub, param)
            parval, x = fn(p=ref['p'], d=ref['d'], quality=True)
            assert parval == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)
            assert x == approx(ref['x'], rel=1e-5, abs=1e-2)

    def test_pv(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(p=ref['p'], v=ref['v']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_Td(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(T=ref['T'], d=ref['d']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_Tv(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(T=ref['T'], v=ref['v']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_ps(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(p=ref['p'], s=ref['s']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_Ts(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(T=ref['T'], s=ref['s']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_ds(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(d=ref['d'], s=ref['s']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_vs(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(v=ref['v'], s=ref['s']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_ph(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(p=ref['p'], h=ref['h']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    @pytest.mark.skip(reason="T&h and T&e aren't viable combos for most cases")
    def test_Th(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMAnalysisError):
            assert fn(T=ref['T'], h=ref['h']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_dh(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(d=ref['d'], h=ref['h']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_vh(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(v=ref['v'], h=ref['h']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_pe(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(p=ref['p'], e=ref['e']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    @pytest.mark.skip(reason="T&h and T&e aren't viable combos for most cases")
    def test_Te(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMAnalysisError):
            assert fn(T=ref['T'], e=ref['e']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_de(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(d=ref['d'], e=ref['e']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_ve(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(v=ref['v'], e=ref['e']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_px(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        if (np.array(ref['x']) < 0).any():
            # Illegal
            with raises(pm.utility.PMParamError):
                fn(p=ref['p'], x=ref['x'])
        else:
            assert fn(p=ref['p'], x=ref['x']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    def test_Tx(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        if (np.array(ref['x']) < 0).any():
            # Illegal
            with raises(pm.utility.PMParamError):
                fn(T=ref['T'], x=ref['x'])
        else:
            assert fn(T=ref['T'], x=ref['x']) == approx(ref[param], rel=1e-5, abs=1e-2, nan_ok=True)

    # Always Unsupported Cases
    def test_dv(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMParamError):
            assert fn(d=ref['d'], v=ref['v'])

    def test_hs(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMParamError):
            assert fn(h=ref['h'], s=ref['s'])

    def test_es(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMParamError):
            assert fn(e=ref['e'], s=ref['s'])

    def test_he(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMParamError):
            assert fn(h=ref['h'], e=ref['e'])

    def test_hx(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMParamError):
            assert fn(h=ref['h'], x=ref['x'])

    def test_ex(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMParamError):
            assert fn(e=ref['e'], x=ref['x'])

    def test_sx(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMParamError):
            assert fn(s=ref['s'], x=ref['x'])

    def test_dx(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        if (np.array(ref['x']) < 0).any():
            # Illegal
            with raises(pm.utility.PMParamError):
                fn(d=ref['d'], x=ref['x'])
        else:
            # Legal but not supported
            with raises(pm.utility.PMParamError):
                fn(d=ref['d'], x=ref['x'])

    def test_vx(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        if (np.array(ref['x']) < 0).any():
            # Illegal
            with raises(pm.utility.PMParamError):
                fn(v=ref['v'], x=ref['x'])
        else:
            # Legal but not supported
            with raises(pm.utility.PMParamError):
                fn(d=ref['d'], x=ref['x'])


class TestState:

    @pytest.fixture(params=(refs[key] for key in refs), ids=(refs.keys()))
    def refdat(self, request):
        return {'sub': pm.get(request.param['sub']),
                'data': request.param['props']}

    @pytest.fixture(params=('p', 'T', 'd', 'v', 'e', 'h', 's', 'x', 'f', 'g'))
    def prop(self, request):
        return request.param

    def test_state_pTx(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(T=ref['T'], p=ref['p'], x=ref['x'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_pd(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(p=ref['p'], d=ref['d'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_pv(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(p=ref['p'], v=ref['v'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_Td(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(T=ref['T'], d=ref['d'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_Tv(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(T=ref['T'], v=ref['v'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_ps(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(p=ref['p'], s=ref['s'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_Ts(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(T=ref['T'], s=ref['s'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_ds(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(d=ref['d'], s=ref['s'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_vs(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(v=ref['v'], s=ref['s'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_ph(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(p=ref['p'], h=ref['h'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    @pytest.mark.skip(reason="T&h and T&e aren't viable combos for most cases")
    def test_state_Th(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        with raises(pm.utility.PMAnalysisError):
            state = sub.state(T=ref['T'], h=ref['h'])

    def test_state_dh(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(d=ref['d'], h=ref['h'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_vh(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(v=ref['v'], h=ref['h'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_pe(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(p=ref['p'], e=ref['e'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    @pytest.mark.skip(reason="T&h and T&e aren't viable combos for most cases")
    def test_state_Te(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        with raises(pm.utility.PMAnalysisError):
            state = sub.state(T=ref['T'], e=ref['e'])

    def test_state_de(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(d=ref['d'], e=ref['e'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_ve(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(v=ref['v'], e=ref['e'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)
