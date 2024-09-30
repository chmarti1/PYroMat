import pyromat as pm
from pyromat.igtools import IgtMix
from pyromat import igtools as igt
import numpy as np
from pytest import approx, raises
import pytest

class TestDefinitions:
    @pytest.fixture(params=('ig.O2', 'ig.BH3O3', 'ig.air'),
                    ids=('oxygen-ig2', 'BH3O3-ig', 'air-ig'))
    def gas(self, request):
        o = pm.get(request.param)
        return o

    def test_string_def(self):
        air = IgtMix('O2 + 3.76 N2')
        assert air.nsubst() == 2
        assert air['N2'] == 3.76
        assert air['O2'] == 1

    def test_string_def_molar(self):
        air = IgtMix('O2 + 3.76 N2', units='kmol')
        assert air.molar() == approx(4.76)

    def test_kwd_def(self):
        air = IgtMix(N2=3.76, O2=1)
        assert air.nsubst() == 2
        assert air['N2'] == 3.76
        assert air['O2'] == 1

    def test_empty_def(self):
        air = IgtMix(['N2', 'O2'])
        assert air.nsubst() == 2
        assert air.mass() == 0

    def test_empty_def_addvals(self):
        air = IgtMix(['N2', 'O2'])
        air['N2'] = 3.76
        air['O2'] = 1
        assert air.mass() == 4.76

    def test_dict_def(self):
        air = IgtMix({'N2': 3.76, 'O2': 1})
        assert air.nsubst() == 2
        assert air['N2'] == 3.76
        assert air['O2'] == 1

    @pytest.mark.skip("Currently failing")
    def test_summation_def(self):
        o2 = IgtMix('1 O2')
        n2 = IgtMix('3.76 N2')
        air = o2 + n2
        assert air.nsubst() == 2
        assert air['N2'] == 3.76
        assert air['O2'] == 1

    def test_from_igmix_def(self):
        air = igt.fromigmix(pm.get('ig.air'))
        print(air)
        assert air.nsubst() == 4

#
# def complete_props_theory(propdict, subid):
#     """
#     Fill in some properties based on theoretical Ideal Gas values
#     """
#     sub = pm.get(subid)
#
#     newdict = propdict.copy()
#
#     newdict['d'] = newdict['p'] * 1e5 / (sub.R() * 1000 * newdict['T'])
#     newdict['v'] = 1/newdict['d']
#     newdict['e'] = newdict['h'] - sub.R()*newdict['T']
#     newdict['cv'] = newdict['cp'] - sub.R()
#     newdict['gam'] = newdict['cp']/newdict['cv']
#     newdict['f'] = newdict['e'] - newdict['T'] * newdict['s']
#     newdict['g'] = newdict['h'] - newdict['T'] * newdict['s']
#
#     return newdict
#
# # Define some air properties
# refs = {}
# refs["air_a"] = {
#     'props': {
#         # values generated by pyromat using T&p
#         'T': 300,  # K
#         'p': 1,  # bar
#         'h': -2.4071345,  # kJ/kg
#         # 's': 6.7077026,  # kJ/kg K  pre smix value
#         's': 6.87040994,  #kJ/kg K
#         'cp': 1.00483493  # kJ/kg K
#     },
#     'sub': 'ig.air',
#     'comment': 'Room Temp and Pressure'
# }
# refs['air_a']['props'] = complete_props_theory(refs['air_a']['props'], refs['air_a']['sub'])
#
# refs["air_b"] = {
#     'props': {
#         # values generated by pyromat using T&p
#         'T': 800,  # K
#         'p': 1,  # bar
#         'h': 519.47733875,  # kJ/kg
#         # 's': 7.72402011,  # kJ/kg K  pre smix value
#         's': 7.88672745,  # kJ/kg K
#         'cp': 1.09862262  # kJ/kg K
#     },
#     'sub': 'ig.air',
#     'comment': 'Increased Temp, Room Pressure'
# }
# refs['air_b']['props'] = complete_props_theory(refs['air_b']['props'], refs['air_b']['sub'])
#
# refs["air_c"] = {
#     'props': {
#         # values generated by pyromat using T&p
#         'T': 900,  # K
#         'p': 1*(900/800)**(1.34859507/(1.34859507-1)),  # bar, gamma of 850K = 1.34859507
#         'h': 630.51678945,  # kJ/kg
#         # 's': 7.72397993,  # kJ/kg K  pre smix value
#         's': 7.88672745,  # kJ/kg K
#         'cp': 1.12170972  # kJ/kg K
#     },
#     'sub': 'ig.air',
#     'comment': 'Isentropic from state b'
# }
# refs['air_c']['props'] = complete_props_theory(refs['air_c']['props'], refs['air_c']['sub'])
#
# refs["air_d"] = {
#     'props': {
#         # values generated by pyromat using T&p
#         'T': 1800,  # K
#         'p': 1000,  # bar
#         'h': 1699.26825827,  # kJ/kg
#         # 's': 6.69042738,  # kJ/kg K  pre smix value
#         's': 6.85313473,  # kJ/kg K
#         'cp': 1.23698868  # kJ/kg K
#     },
#     'sub': 'ig.air',
#     'comment': 'High T & P'
# }
# refs['air_d']['props'] = complete_props_theory(refs['air_d']['props'], refs['air_d']['sub'])
#

from test_ig import refs
refs = {key: refs[key] for key in ['air_a', 'air_b', 'air_c', 'air_d']}

class TestRefs:

    @pytest.fixture(params=(refs[key] for key in refs), ids=(refs.keys()))
    def refdat(self, request):
        return {'sub': IgtMix('0.01289563 Ar + 0.00047711 CO2 + 0.75520558 N2 + 0.23142168 O2'),
                'data': request.param['props']}

    @pytest.fixture(params=('p', 'T', 'd', 'e', 'h', 's', 'cp', 'cv', 'gam', 'f', 'g'))
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

    @pytest.mark.skip("d-s not implemented")
    def test_ds(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        assert fn(d=ref['d'], s=ref['s']) == approx(ref[param], rel=1e-5, abs=1e-1)

    @pytest.mark.skip("v not implemented")
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

    @pytest.mark.skip("v not implemented")
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

    @pytest.mark.skip("v not implemented")
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
        assert fn(e=ref['e'], s=ref['s']) == approx(ref[param], rel=1e-5, abs=1e-1)

    def test_he(self, param, refdat):
        sub, ref = refdat['sub'], refdat['data']
        fn = getattr(sub, param)
        with raises(pm.utility.PMParamError):
            fn(h=ref['h'], e=ref['e'])