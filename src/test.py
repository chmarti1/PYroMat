#!bin/python

import pyromat as pm
import numpy as np
import matplotlib.pyplot as plt

S = pm.get('mp.N2')
Tmax,Tmin = S.Tlim()
T = np.random.random(10000) * (Tmax - Tmin) + Tmin
d = np.random.random(10000) * S.data['dlim'][1]
_,d1,d2,x,I = S._argparse(T=T,d=d)
s = S.s(T=T,d=d)
h = S.h(T=T,d=d)
p = S.p(T=T,d=d)
e = S.e(T=T,d=d)

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

def test(**kwarg):
    T_, d1_, d2_, x_, I_ = S._argparse(**kwarg)
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

