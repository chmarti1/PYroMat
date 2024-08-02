import pyromat as pm
import numpy as np
from pytest import approx, raises
import pytest


class TestInputErrors:
    @pytest.fixture(params=('ig.O2', 'ig.BH3O3', 'ig.air'),
                    ids=('oxygen-ig2', 'BH3O3-ig', 'air-ig'))
    def gas(self, request):
        o = pm.get(request.param)
        return o

    def test_T_only_varg(self, gas):
        # Based on how .T is written, it interprets this as a temperature given
        # So just make sure that works
        assert gas.T(300) == gas.T(T=300, p=pm.config['def_p'])

    def test_T_double_specified(self, gas):
        with raises(pm.utility.PMParamError):
            gas.T(300, T=300)

    def test_p_double_specified(self, gas):
        with raises(pm.utility.PMParamError):
            gas.T(300, 1, p=1)

    def test_varg_specified(self, gas):
        assert gas.s(300, 1) == gas.s(p=1, T=300)

    def test_varg_toomany(self, gas):
        with raises(pm.utility.PMParamError):
            gas.T(300, 1, 5)

    def test_any_arg_toomany(self, gas):
        with raises(pm.utility.PMParamError):
            gas.T(300, 1, s=5)
        with raises(pm.utility.PMParamError):
            gas.T(T=300, p=1, s=5)

    def test_default_usage(self, gas):
        assert gas.s() == gas.s(T=pm.config['def_T'], p=pm.config['def_p'])
        assert gas.s(T=pm.config['def_T']) == gas.s(T=pm.config['def_T'], p=pm.config['def_p'])
        assert gas.s(p=pm.config['def_p']) == gas.s(T=pm.config['def_T'], p=pm.config['def_p'])

    def test_invarg_toomany(self, gas):
        with raises(pm.utility.PMParamError):
            gas.T(h=300, s=5)
        with raises(pm.utility.PMParamError):
            gas.T(e=300, s=5)

    def test_d_v_collision(self, gas):
        with raises(pm.utility.PMParamError):
            gas.T(d=1, v=1)

    def test_illegal_prop(self, gas):
        with raises(pm.utility.PMParamError):
            gas.T(T=300, fakeprop=10000)

    def test_T_oob_all(self, gas):
        with raises(pm.utility.PMParamError):
            gas.T(T=30000, p=1)

    def test_T_oob_some(self, gas):
        vals = gas.T(T=[300, 30000], p=1)
        assert vals[0] == 300
        assert np.isnan(vals[1])


class TestGeneral:

    @pytest.mark.parametrize('sub', ('ig.O2', 'ig.BH3O3', 'ig.air'), ids=('oxygen-ig2', 'boric acid-ig', 'air-ig'))
    def test_get(self, sub):
        # Just confirm no error
        igobj = pm.get(sub)

    @ pytest.mark.parametrize('gas, expect', [('ig.O2',(200, 6000)), ('ig.BH3O3',(298, 6000)), ('ig.air',(200, 6000))])
    def test_Tlim(self, gas, expect):
        assert pm.get(gas).Tlim() == approx(expect)

    @ pytest.mark.parametrize('gas, expect', [('ig.O2', 31.9988), ('ig.BH3O3', 61.833), ('ig.air', 28.96493)])
    def test_mw(self, gas, expect):
        assert pm.get(gas).mw() == approx(expect)

    @pytest.mark.parametrize('gas, expect', [('ig.O2', pm.units.const_Ru/31.9988), ('ig.BH3O3', pm.units.const_Ru/61.833), ('ig.air', pm.units.const_Ru/28.96493)])
    def test_R(self, gas, expect):
        assert pm.get(gas).R() == approx(expect)

    @pytest.mark.parametrize('gas, expect', [
        ('ig.O2', {'O': 2}),
        ('ig.BH3O3', {'B': 1, 'H': 3, 'O': 3}),
        ('ig.air', {'Ar': 0.009350187003740077, 'C': 0.00031400628012560254, 'O': 0.4195883917678354, 'N': 1.5617112342246846})
    ])
    def test_atoms(self, gas, expect):
        assert pm.get(gas).atoms() == approx(expect)

    @pytest.mark.parametrize('gas, expect', [
        ('ig.air', {'ig.Ar': 0.009350187003740077, 'ig.CO2': 0.00031400628012560254, 'ig.N2': 0.7808556171123423, 'ig.O2': 0.2094801896037921})
    ])
    def test_X(self, gas, expect):
        assert pm.get(gas).X() == approx(expect)

    @pytest.mark.parametrize('gas, expect', [
        ('ig.air', {'ig.Ar': 0.012895634840195168, 'ig.CO2': 0.0004771062632750561, 'ig.N2': 0.7552055804206431, 'ig.O2': 0.23142167847588652})
    ])
    def test_Y(self, gas, expect):
        assert pm.get(gas).Y() == approx(expect)

def complete_props_theory(propdict, subid):
    """
    Fill in some properties based on theoretical Ideal Gas values
    """
    sub = pm.get(subid)

    newdict = propdict.copy()

    newdict['d'] = newdict['p'] * 1e5 / (sub.R() * 1000 * newdict['T'])
    newdict['v'] = 1/newdict['d']
    newdict['e'] = newdict['h'] - sub.R()*newdict['T']
    newdict['cv'] = newdict['cp'] - sub.R()
    newdict['gam'] = newdict['cp']/newdict['cv']

    return newdict


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


# Four reference points for oxygen computed using pyromat, due to variability
# of reference states.
# Note that specific volumes are computed as 1/d for accuracy.
refs = {}
refs["oxygen_a"] = {
    'props': {
        # values generated by pyromat using T&p
        'T': 300,  # K
        'p': 1,  # bar
        'h': 1.69877563,  # kJ/kg
        's': 6.41680485,  # kJ/kg K
        'cp': 0.91841166  # kJ/kg K
    },
    'sub': 'ig.O2',
    'comment': 'Room Temp and Pressure'
}
refs['oxygen_a']['props'] = complete_props_theory(refs['oxygen_a']['props'], refs['oxygen_a']['sub'])

refs["oxygen_b"] = {
    'props': {
        # values generated by pyromat using T&p
        'T': 800,  # K
        'p': 1,  # bar
        'h': 494.96011806,  # kJ/kg
        's': 7.37301205,  # kJ/kg K
        'cp': 1.05471624  # kJ/kg K
    },
    'sub': 'ig.O2',
    'comment': 'Increased Temp, Room Pressure'
}
refs['oxygen_b']['props'] = complete_props_theory(refs['oxygen_b']['props'], refs['oxygen_b']['sub'])

refs["oxygen_c"] = {
    'props': {
        # values generated by pyromat using T&p
        'T': 900,  # K
        'p': 1*(900/800)**(1.32282754/(1.32282754-1)),  # bar, gamma of 850K = 1.32282754
        'h': 601.41481211,  # kJ/kg
        's': 7.37301205,  # kJ/kg K
        'cp': 1.07371185  # kJ/kg K
    },
    'sub': 'ig.O2',
    'comment': 'Isentropic from state b'
}
refs['oxygen_c']['props'] = complete_props_theory(refs['oxygen_c']['props'], refs['oxygen_c']['sub'])

refs["oxygen_d"] = {
    'props': {
        # values generated by pyromat using T&p
        'T': 1800,  # K
        'p': 1000,  # bar
        'h': 1614.04575249,  # kJ/kg
        's': 6.47977787,  # kJ/kg K
        'cp': 1.16705079  # kJ/kg K
    },
    'sub': 'ig.O2',
    'comment': 'High T & P'
}
refs['oxygen_d']['props'] = complete_props_theory(refs['oxygen_d']['props'], refs['oxygen_c']['sub'])

# Array
refs['oxygen_array'] = make_props_array(refs, refs['oxygen_a']['sub'])


refs["boricacid_a"] = {
    'props': {
        # values generated by pyromat using T&p
        'T': 300,  # K
        'p': 1,  # bar
        'h': -16045.78193466,  # kJ/kg
        's': 4.78111299,  # kJ/kg K
        'cp': 1.05979064  # kJ/kg K
    },
    'sub': 'ig.BH3O3',
    'comment': 'Room Temp and Pressure'
}
refs['boricacid_a']['props'] = complete_props_theory(refs['boricacid_a']['props'], refs['boricacid_a']['sub'])

refs["boricacid_b"] = {
    'props': {
        # values generated by pyromat using T&p
        'T': 800,  # K
        'p': 1,  # bar
        'h': -15328.26168563,  # kJ/kg
        's': 6.13679451,  # kJ/kg K
        'cp': 1.70773773  # kJ/kg K
    },
    'sub': 'ig.BH3O3',
    'comment': 'Increased Temp, Room Pressure'
}
refs['boricacid_b']['props'] = complete_props_theory(refs['boricacid_b']['props'], refs['boricacid_b']['sub'])

refs["boricacid_c"] = {
    'props': {
        # values generated by pyromat using T&p
        'T': 900,  # K
        'p': 1*(900/800)**(1.08349905/(1.08349905-1)),  # bar, gamma of 850K = 1.08349905
        'h': -15153.82697422,  # kJ/kg
        's': 6.13679451,  # kJ/kg K
        'cp': 1.7789004  # kJ/kg K
    },
    'sub': 'ig.BH3O3',
    'comment': 'Isentropic from state b'
}
refs['boricacid_c']['props'] = complete_props_theory(refs['boricacid_c']['props'], refs['boricacid_c']['sub'])

refs["boricacid_d"] = {
    'props': {
        # values generated by pyromat using T&p
        'T': 1800,  # K
        'p': 1000,  # bar
        'h': -13367.19402821,  # kJ/kg
        's': 6.77542321,  # kJ/kg K
        'cp': 2.12425583  # kJ/kg K
    },
    'sub': 'ig.BH3O3',
    'comment': 'High T & P'
}
refs['boricacid_d']['props'] = complete_props_theory(refs['boricacid_d']['props'], refs['boricacid_d']['sub'])

# Array
refs['boricacid_array'] = make_props_array(refs, refs['boricacid_a']['sub'])


refs["air_a"] = {
    'props': {
        # values generated by pyromat using T&p
        'T': 300,  # K
        'p': 1,  # bar
        'h': -2.4071345,  # kJ/kg
        # 's': 6.7077026,  # kJ/kg K  pre smix value
        's': 6.87040994,  #kJ/kg K
        'cp': 1.00483493  # kJ/kg K
    },
    'sub': 'ig.air',
    'comment': 'Room Temp and Pressure'
}
refs['air_a']['props'] = complete_props_theory(refs['air_a']['props'], refs['air_a']['sub'])

refs["air_b"] = {
    'props': {
        # values generated by pyromat using T&p
        'T': 800,  # K
        'p': 1,  # bar
        'h': 519.47733875,  # kJ/kg
        # 's': 7.72402011,  # kJ/kg K  pre smix value
        's': 7.88672745,  # kJ/kg K
        'cp': 1.09862262  # kJ/kg K
    },
    'sub': 'ig.air',
    'comment': 'Increased Temp, Room Pressure'
}
refs['air_b']['props'] = complete_props_theory(refs['air_b']['props'], refs['air_b']['sub'])

refs["air_c"] = {
    'props': {
        # values generated by pyromat using T&p
        'T': 900,  # K
        'p': 1*(900/800)**(1.34859507/(1.34859507-1)),  # bar, gamma of 850K = 1.34859507
        'h': 630.51678945,  # kJ/kg
        # 's': 7.72397993,  # kJ/kg K  pre smix value
        's': 7.88672745,  # kJ/kg K
        'cp': 1.12170972  # kJ/kg K
    },
    'sub': 'ig.air',
    'comment': 'Isentropic from state b'
}
refs['air_c']['props'] = complete_props_theory(refs['air_c']['props'], refs['air_c']['sub'])

refs["air_d"] = {
    'props': {
        # values generated by pyromat using T&p
        'T': 1800,  # K
        'p': 1000,  # bar
        'h': 1699.26825827,  # kJ/kg
        # 's': 6.69042738,  # kJ/kg K  pre smix value
        's': 6.85313473,  # kJ/kg K
        'cp': 1.23698868  # kJ/kg K
    },
    'sub': 'ig.air',
    'comment': 'High T & P'
}
refs['air_d']['props'] = complete_props_theory(refs['air_d']['props'], refs['air_d']['sub'])

# Array
refs['air_array'] = make_props_array(refs, refs['air_a']['sub'])


class TestRefs:

    @pytest.fixture(params=(refs[key] for key in refs), ids=(refs.keys()))
    def refdat(self, request):
        return {'sub': pm.get(request.param['sub']),
                'data': request.param['props']}

    @pytest.fixture(params=('p', 'T', 'd', 'v', 'e', 'h', 's', 'cp', 'cv', 'gam'))
    def param(self, request):
        return request.param

    def test_pT(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(p=ref['p'], T=ref['T']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_pd(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(p=ref['p'], d=ref['d']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_pv(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(p=ref['p'], v=ref['v']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_Td(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(T=ref['T'], d=ref['d']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_Tv(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(T=ref['T'], v=ref['v']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_ps(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(p=ref['p'], s=ref['s']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_Ts(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(T=ref['T'], s=ref['s']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_ds(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(d=ref['d'], s=ref['s']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_vs(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(v=ref['v'], s=ref['s']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_ph(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(p=ref['p'], h=ref['h']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_dh(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(d=ref['d'], h=ref['h']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_vh(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(v=ref['v'], h=ref['h']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_pe(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(p=ref['p'], e=ref['e']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_de(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(d=ref['d'], e=ref['e']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_ve(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(v=ref['v'], e=ref['e']) == approx(ref[param], rel=1e-5, abs=1e-1)

    # Always Unsupported Cases
    def test_Th(self, param, refdat):  # Can't work for ideal gases
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMParamError):
            fn(T=ref['T'], h=ref['h'])

    def test_Te(self, param, refdat):  # Can't work for ideal gases
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMParamError):
            fn(T=ref['T'], e=ref['e'])

    def test_dv(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMParamError):
            fn(d=ref['d'], v=ref['v'])

    def test_hs(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(h=ref['h'], s=ref['s']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_es(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(v=ref['e'], e=ref['s']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_he(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMParamError):
            fn(h=ref['h'], e=ref['e'])


class TestState:

    @pytest.fixture(params=(refs[key] for key in refs), ids=(refs.keys()))
    def refdat(self, request):
        return {'sub': pm.get(request.param['sub']),
                'data': request.param['props']}

    @pytest.fixture(params=('p', 'T', 'd', 'v', 'e', 'h', 's', 'cp', 'cv', 'gam'))
    def prop(self, request):
        return request.param

    def test_state_pT(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(T=ref['T'], p=ref['p'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_pd(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(p=ref['p'], d=ref['d'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_pv(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(p=ref['p'], v=ref['v'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_Td(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(T=ref['T'], d=ref['d'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_Tv(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(T=ref['T'], v=ref['v'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_ps(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(p=ref['p'], s=ref['s'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_Ts(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(T=ref['T'], s=ref['s'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_ds(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(d=ref['d'], s=ref['s'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_vs(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(v=ref['v'], s=ref['s'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_ph(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(p=ref['p'], h=ref['h'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_Th(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        with raises(pm.utility.PMParamError):
            state = sub.state(T=ref['T'], h=ref['h'])

    def test_state_dh(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(d=ref['d'], h=ref['h'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_vh(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(v=ref['v'], h=ref['h'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_pe(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(p=ref['p'], e=ref['e'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_Te(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        with raises(pm.utility.PMParamError):
            state = sub.state(T=ref['T'], e=ref['e'])

    def test_state_de(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(d=ref['d'], e=ref['e'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)

    def test_state_ve(self, refdat, prop):
        sub, ref = refdat['sub'], refdat['data']
        state = sub.state(v=ref['v'], e=ref['e'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-1)
