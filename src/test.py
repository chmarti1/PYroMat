#!bin/python

import pyromat as pm
import numpy as np
import matplotlib.pyplot as plt
import time

S = pm.get('mp.N2')
Tmin,Tmax = S.Tlim()
T = np.random.random(10000) * (Tmax - Tmin) + Tmin
d = np.random.random(10000) * S.data['dlim'][1]
_,d1,d2,x,I = S._argparse(T=T,d=d)
s = S.s(T=T,d=d)
h = S.h(T=T,d=d)
p = S.p(T=T,d=d)
e = S.e(T=T,d=d)

def generate(subst, N=10000):
    Tmin, Tmax = subst.Tlim()
    _,pmax = S.plim()
    dmax = S.d(T=Tmin, p=pmax)
    T = np.random.random(N) * (Tmax - Tmin) + Tmin
    d = np.random.random(N) * dmax
    p = S.p(T=T, d=d)
    I = p > pmax
    if I.any():
        dd = np.random.random(np.sum(I)) * dmax
        d[I] = dd
        p = S.p(T=T[I], d=dd)
        I[I] = p > pmax
        
    return S.state(T=T, d=d)


def benchmark(subst, N=10, **kwarg):
    Targ = 0
    Tstate = 0
    for count in range(N):
        print(count)
        # Time argparse()
        start = time.time()
        subst._argparse(**kwarg)
        this = time.time()-start
        Targ += this
        # Time the state method
        start = time.time()
        subst.state(**kwarg)
        this = time.time()-start
        Tstate += this
    Targ /= N
    Tstate /= N
    Tstate -= Targ
    return Targ, Tstate
    
def auto_benchmark(subst, prop1, prop2, N=[10,100,1000,10000,100000]):
    Targ = []
    Tstate = []
    for n in N:
        print('==>', n, '<==')
        state = generate(subst, n)
        args = {prop1:state[prop1], prop2:state[prop2]}
        targ, tstate = benchmark(subst, **args)
        Targ.append(targ)
        Tstate.append(tstate)
    return N, Targ, Tstate

#   7.1: x is specified
#       7.1.1: x,T,p        GOOD
#       7.1.2: x,T          
#       7.1.3: x,p
#   7.2: T,?
#       7.2.1: T,d          GOOD
#       7.2.2: T,p          GOOD -- precise equality with saturation states give inconsistent (but stable results)
#       7.2.3: T + inverse  GOOD -- occasionally misidentifies points numerically close to the saturation line
#   7.3: d,?
#       7.3.1: d + inverse  GOOD -- LOOKS LIKE d,g IS INVALID
#   7.4: p,?
#       p + inverse         GOOD -- not allowed with g or f
#   7.5: g,?
#       7.5: g + inverse    CRASHES
#   7.6: ?,?
#       Any two remaining inverse       CRASHES

# In N2 (mp2)
# (s,h) fails with density errors at
# T = array([664.27513038, 582.74124236, 663.40093498, 591.07042537,
#       742.3821707 , 563.61851251, 616.06854254, 553.1869754 ,
#       597.64432506, 855.97881845, 645.293752  , 586.47720309,
#       444.9969719 , 612.90227   , 472.89750766, 678.3283894 ,
#       611.38801314])
# d = array([ 0.94504185, 11.19105334,  2.19517483,  9.71148211,  2.15269466,
#        1.97054596,  4.18058009, 11.80185871,  7.9976818 ,  2.82713949,
#        3.11305451,  9.33314563,  6.87374511,  8.05231862, 65.24656384,
#        5.52651253,  2.07673124])
#
# (s,h) fails with NaN at
# T = array([91.21508767, 86.71596801, 67.78798953, 71.02038424, 89.64740762,
#       89.13277341, 80.20867577, 89.7364947 , 76.1229656 ])
# d = array([740.16405429, 778.96422903, 857.63224626, 850.49609392,
#       761.65375645, 756.48328942, 804.08683827, 761.64942785,
#       818.30321754])
#
# (s,e) fails with density errors at
# T = array([663.40093498, 742.3821707 ])
# d = array([2.19517483, 2.15269466])
#
# (s,e) fails with NaN at
# T = array([63.18018277])
# d = array([945.12284745])
# 
# (d,p) fails with NaN at
# T = array([70.04344752, 72.03489332, 69.09605407])
# d = array([1414.2209862 , 1414.61209767, 1401.26925118])
#
# There are persistent out-of-bounds errors at high-density and low 
# temperature.  The same errors also occur at low temperature liquid just
# outside of the dome.  Are these errors actually out-of-bounds!?



def test(subst, state, p1, p2):
    d = state['d']
    T, d1, d2, x, I = subst._argparse(T=state['T'], d=d)
    args = {p1:state[p1], p2:state[p2]}
    T_, d1_, d2_, x_, I_ = subst._argparse(**args)
    eT = np.abs(T_ - T)
    e1 = np.abs(d1_ - d1)
    e2 = np.abs(d2_ - d2)
    
    maxT = np.max(eT)
    max1 = np.max(e1)
    max2 = np.max(e2)
    
    print(f'Max errors: T: {maxT}, d1: {max1}, d2: {max2}')
    
    IT = np.abs(eT)>1e-6*T
    I1 = np.abs(e1)>1e-6*d1
    I2 = np.abs(e2)>1e-6*d2
    IN = np.isnan(T_) + np.isnan(d1_) + np.isnan(d2_) + np.isnan(x_)

    NT = np.sum(IT)
    N1 = np.sum(I1)
    N2 = np.sum(I2)
    NN = np.sum(IN)

    print(f'Potential Problems: T: {NT}, d1: {N1}, d2: {N2}, NaN: {NN}')
    
    fig,ax = plt.subplots(3,1)
    ax[0].semilogy(T, eT, '.', ms=1)
    ax[1].semilogy(d1, e1, '.', ms=1)
    ax[2].semilogy(d2, e2, '.', ms=1)
    
    fig,ax = plt.subplots(1,1)
    ax.plot(d[IT], T[IT], 'g^', ms=3)
    ax.plot(d[I1], T[I1], 'bx', ms=3)
    ax.plot(d[I2], T[I2], 'ro', ms=2)
    ax.plot(d[IN], T[IN], 'ms', ms=3)
    ax.plot(subst._sattable['dL'], subst._sattable['T'], 'k')
    ax.plot(subst._sattable['dV'], subst._sattable['T'], 'k')
    
    plt.show()
    
    out = {'IT':IT, 'I1':I1, 'I2':I2, 'IN':IN, 'T':T_, 'd1': d1_, 'd2':d2_, 'x':x_, 'I':I_}
    
    return out



def failsearch(**kwarg):
    f0str, f1str = kwarg.keys()
    f0value = kwarg[f0str]
    f1value = kwarg[f1str]
    a = 0
    b = f0value.size
    while b-a > 1:
        c = int(0.5 * (a+b))
        print(a,c)
        try:
            S._argparse(**{f0str:f0value[a:c], f1str:f1value[a:c]})
            a = c
        except:
            b = c
    return {'T':T[a], 'd':d[a], 'p':p[a], 's':s[a], 'e':e[a], 'h':h[a], 'd1':d1[a], 'd2':d2[a]}

def renew():
    global S
    pm.reg.regload()
    pm.dat.load()
    S = pm.get('mp.N2')

