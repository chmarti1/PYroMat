import pyromat as pm
import numpy as np
from pytest import approx, raises
import pytest


class TestDevBugs:

    def test_bug_git42(self):
        sub = pm.get('mp.H2O')
        assert not np.isnan(sub.Ts(p=110))

    def test_bug_git43(self):
        # Cases that failed due to the density being under the dome and yielding
        # insane pressure results in the limits of the temperature iterator.
        sub = pm.get('mp.H2O')
        print(sub.T(p=2497, d=525))  # errors with the old method, works with the new, expect around 1026
        print(sub.T(p=2497, d=1090))  # density doesn't intersect the dome, should work for old or new method, expect around 291
        with raises(pm.utility.PMParamError):
            print(sub.T(p=2497, d=1100))  # really out of range, sub.d(p=2497, T=Tmin is around 1097), expect errors

    def test_bug_git44(self):
        # Cases that failed due to the pressure being out of gamut on some of the density
        # limits. Corrected by giving better density limits.
        sub = pm.get('mp.CH4')
        print(sub.d(p=0.11707696, T=624.375))
        print(sub.d(p=9999.0, T=624.375))

        # Cases where low P&T yield an extremely high density that is legitimately out
        # of bounds on the density limits. These cases really should fail.
        # sub = pm.get('mp.CO2')
        # # print(sub.d(p=7959, T=216.8))
        # # In this case our density exceeds dmax, so p&T in combo put us out of gamut
        # print(sub.p(d=sub.data['dlim'][1], T=216.8))
        # sub = pm.get('mp.C3H2F4_1')
        # # print(sub.d(p=999, T=169.2))
        # print(sub.p(d=sub.data['dlim'][1], T=169.2))
        # sub = pm.get('mp.C2H2F4')
        # # print(sub.d(p=699.3, T=170))
        # print(sub.p(d=sub.data['dlim'][1], T=170))
        # sub = pm.get('mp.H2O')
        # print(sub.p(d=sub.data['dlim'][1], T=sub.data['Tlim'][0]))

    def test_bug_git45(self):
        # p vs T curve has two values. Iteration for _p() from T fails
        sub = pm.get('mp.N2')
        # Input state ref = sub.state(p=1.6e4, T=68)
        sub.T(p=1.6e4, d=1355)

        # s vs d curve has two values, iteration for _d() from s fails.
        sub = pm.get('mp.CH4')
        # ref = sub.state(p=2475, T=223) basis reference state
        ref = sub.d(T=223, s=-5.5)  # 460 or so

    def test_bug_git46(self):
        co2 = pm.get('mp.CO2')
        assert co2.h(T=190., p=7900.) == approx(-64.14)  # raises error

    def test_bug_git47(self):
        w = pm.get('mp.H2O')

        # Here's a case where we have some points that should be np.nan
        h = [1694.095, 4372, 4373]  # kJ/kg
        p = [5, 100, 1000]  # bar
        ans = w.state(h=h, p=p)  # expect ans['x'] = [0.5, -1, np.nan]
        for prop in ans:
            print(prop)
            assert not np.isnan(ans[prop][0])
            assert not np.isnan(ans[prop][1])
            assert np.isnan(ans[prop][2])

        # It can also happen to the temperature
        T = w.state(T=300, p=[1, 1e10])['T']
        assert not np.isnan(T[0])
        assert np.isnan(T[1])

    def test_bug_git52(self):
        # ig, ig2, igmix, mp1
        subs = [pm.get("ig.BH3O3"), pm.get("ig.O2"), pm.get('ig.air'), pm.get('mp.H2O')]

        for sub in subs:
            p = 1
            T = 500
            h = sub.h(T=T, p=p)
            s = sub.s(T=T, p=p)
            assert sub.T_s(s, p) == approx(sub.T_s(s, p=p))
            assert sub.T_s(s, p) == approx(sub.T_s(s=s, p=p))
            assert sub.T_h(h, p) == approx(sub.T_h(h, p=p))
            assert sub.T_h(h, p) == approx(sub.T_h(h=h, p=p))
            if hasattr(sub, "p_s"):
                assert sub.p_s(s, T) == approx(sub.p_s(s, T=T))
                assert sub.p_s(s, T) == approx(sub.p_s(s=s, T=T))
            if hasattr(sub, "d_s"):
                assert sub.d_s(s, T) == approx(sub.d_s(s, T=T))
                assert sub.d_s(s, T) == approx(sub.d_s(s=s, T=T))
