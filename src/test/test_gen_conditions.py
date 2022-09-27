import numpy as np
import pytest
from pytest import approx, raises

import pyromat as pm


# Note that this generates an extremely large number of tests, rougly 800k
# for ig. So we're currently running things for a single PM class at a time.
# We should look for ways to manage this or separate things into different
# runs.

# The purpose of these tests is really to hit every single substance to make
# sure that our calculations are reliable across all, rather than the few
# targeted substances listed previously.

n = 5  # Number of points per temperature/pressure range
pmclass = 'mp1'  # The class to test: one of igmix, ig, ig2, mp1


# Which parameters to test the outputs of based on class
params = {'mp1': ('p', 'T', 'd', 'v', 'e', 'h', 's', 'x'),
          'ig': ('p', 'T', 'd', 'v', 'e', 'h', 's'),
          'ig2': ('p', 'T', 'd', 'v', 'e', 'h', 's'),
          'igmix': ('p', 'T', 'd', 'v', 'e', 'h', 's')
          }


def getT(sub, i, n):
    """
    Get a temperature from the standard list by index using a uniform
    methodology by substance. This is necessary because each substance has its
    own legal temperatures, but the pytest.fixtures require that everything
    be prespecified. Thus we feed the fixtures a list of indices, and allow
    this function to generate the actual data.
    """
    if sub.data['id'].startswith('mp'):
        Tmin, Tmax = sub.Tlim()
        # Tt, pt = sub.triple()
        # Tmin = max(Tmin, Tt)
        Tmax *= 0.999
        Tmin *= 1.001
        Tlist = np.linspace(Tmin, Tmax, n)
    else:
        Tmin, Tmax = sub.Tlim()
        Tmin *= 1.001
        Tmax *= 0.999
        Tlist = np.linspace(Tmin, Tmax, n)
    return Tlist[i]


def getp(sub, i, n):
    """
    Get a pressure from the standard list by index using a uniform
    methodology by substance. This is necessary because each substance has its
    own legal temperatures, but the pytest.fixtures require that everything
    be prespecified. Thus we feed the fixtures a list of indices, and allow
    this function to generate the actual data.
    """
    if sub.data['id'].startswith('mp'):
        pmin, pmax = sub.plim()
        # at low T and high p, we can actually cause density to go out of range
        pmax2 = sub.p(d=sub.data['dlim'][1], T=sub.data['Tlim'][0])
        pmax = min(pmax, pmax2)
        pmax *= 0.999
        Tt, pt = sub.triple()
        pt *= 1.001
        plist = np.linspace(pt, pmax, n)
    else:
        Tmin, Tmax = sub.Tlim()
        pmin, pmax = np.array([sub.p(T=Tmin, d=0.01), sub.p(T=Tmax, d=1000)]).flatten()
        plist = np.linspace(pmin, pmax, n)
    return plist[i]


class TestArbitary:

    @pytest.fixture(params=[item.data['id'] for item in pm.search(pmclass=pmclass)], ids=[item.data['id'] for item in pm.search(pmclass=pmclass)])
    def sub(self, request):
        sub = pm.get(request.param)
        return sub

    @pytest.fixture(params=range(n), ids=["p"+str(i) for i in range(n)])
    def p(self, request, sub):
        return getp(sub, request.param, n)

    @pytest.fixture(params=range(n), ids=["T"+str(i) for i in range(n)])
    def T(self, request, sub):
        return getT(sub, request.param, n)

    @pytest.fixture(params=params[pmclass])
    def prop(self, request):
        return request.param

    def test_state_pTx(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        if sub.data['id'].startswith('mp'):
            state = sub.state(T=ref['T'], p=ref['p'], x=ref['x'])
        else:
            state = sub.state(T=ref['T'], p=ref['p'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_pd(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(p=ref['p'], d=ref['d'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_pv(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(p=ref['p'], v=ref['v'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_Td(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(T=ref['T'], d=ref['d'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_Tv(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(T=ref['T'], v=ref['v'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_ps(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(p=ref['p'], s=ref['s'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_Ts(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(T=ref['T'], s=ref['s'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_ds(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(d=ref['d'], s=ref['s'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_vs(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(v=ref['v'], s=ref['s'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_ph(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(p=ref['p'], h=ref['h'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    # def test_state_Th(self, sub, p, T, prop):
    #     ref = sub.state(p=p, T=T)
    #     if sub.data['id'].startswith('mp'):
    #         err = pm.utility.PMAnalysisError
    #     else:
    #         err = pm.utility.PMParamError
    #     with raises(err):
    #         state = sub.state(T=ref['T'], h=ref['h'])

    def test_state_dh(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(d=ref['d'], h=ref['h'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_vh(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(v=ref['v'], h=ref['h'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_pe(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(p=ref['p'], e=ref['e'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    # def test_state_Te(self, sub, p, T, prop):
    #     ref = sub.state(p=p, T=T)
    #     if sub.data['id'].startswith('mp'):
    #         err = pm.utility.PMAnalysisError
    #     else:
    #         err = pm.utility.PMParamError
    #     with raises(err):
    #         state = sub.state(T=ref['T'], e=ref['e'])

    def test_state_de(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(d=ref['d'], e=ref['e'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)

    def test_state_ve(self, sub, p, T, prop):
        ref = sub.state(p=p, T=T)
        state = sub.state(v=ref['v'], e=ref['e'])
        assert state[prop] == approx(ref[prop], rel=1e-5, abs=1e-2)