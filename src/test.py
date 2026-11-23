#!bin/python

import pyromat as pm
import numpy as np
import matplotlib.pyplot as plt
import time

#
# p, h produces NaN at points very near the saturated liquid line.
# In Water
#   h = np.array([ 166.2512912 , 1101.49739185,  530.50773555,  449.94259603])
#   p = np.array([ 55.52310365, 112.25935846,  86.23674416,  77.39399102])
# In Water
#   Convergence errors near the saturated liquid line with p,h
#   state(p=2., h=505)
#   returns h = 503.906...

def generate(subst, N=10000):
    Tmin, Tmax = subst.Tlim()
    _,pmax = subst.plim()
    dmax = subst.d(T=Tmin, p=pmax)
    T = np.random.random(N) * (Tmax - Tmin) + Tmin
    d = np.random.random(N) * dmax
    p = subst.p(T=T, d=d)
    I = p > pmax
    while I.any():
        N = np.sum(I)
        print('Adjusting densities:', N)
        dd = np.random.random(N) * dmax
        d[I] = dd
        p = S.p(T=T[I], d=dd)
        I[I] = p > pmax
        
    return subst.state(T=T, d=d)


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





def test(subst, state, p1, p2):
    d = state['d']
    T, d1, d2, x, I = subst._argparse(T=state['T'], d=d)
    args = {p1:state[p1], p2:state[p2]}
    T_, d1_, d2_, x_, I_ = subst._argparse(**args)
    eT = np.abs(T_ - T)
    e1 = np.abs(d1_ - d1)
    e2 = np.abs(d2_ - d2)
    
    IT = np.zeros_like(T_, dtype=bool)
    I1 = np.zeros_like(T_, dtype=bool)
    I2 = np.zeros_like(T_, dtype=bool)
    
    IN = np.logical_not(np.isnan(T_) + np.isnan(d1_) + np.isnan(d2_) + np.isnan(x_))
    IT[IN] = eT[IN]>1e-6*T[IN]
    I1[IN] = e1[IN]>1e-6*d1[IN]
    I2[IN] = e2[IN]>1e-6*d2[IN]
    Itest = IT + I1 + I2
    
    maxT = np.max(eT[IN])
    max1 = np.max(e1[IN])
    max2 = np.max(e2[IN])
    
    print(f'Max errors: T: {maxT}, d1: {max1}, d2: {max2}')
    
    IN = np.logical_not(IN)

    NT = np.sum(IT)
    N1 = np.sum(I1)
    N2 = np.sum(I2)
    Ntest = np.sum(Itest)
    NN = np.sum(IN)

    print(f'Potential Problems: T: {NT}, d1: {N1}, d2: {N2}, Unique: {Ntest}, NaN: {NN}')
    
    # Re-calculate p1 and p2 to see if the "bad" solution agrees -- is this
    # a case where there is more than one solution that matches the inputs?
    args = {p1:state[p1][Itest], p2:state[p2][Itest]}
    state_ = subst.state(**args)
    
    p1e = np.abs(state[p1][Itest] - state_[p1])
    p2e = np.abs(state[p2][Itest] - state_[p2])
    
    Ip1 = np.zeros_like(T_, dtype=bool)
    Ip2 = np.zeros_like(T_, dtype=bool)
    
    Ip1[Itest] = p1e > state[p1][Itest] * 1e-6
    Ip2[Itest] = p2e > state[p2][Itest] * 1e-6
    Ibad = Ip1 + Ip2
    
    Np1 = np.sum(Ip1)
    Np2 = np.sum(Ip2)
    Nbad = np.sum(Ibad)
    
    print(f'Incorrect solutions: {p1}: {Np1}, {p2}: {Np2}, Unique: {Nbad}')
    
    fig,ax = plt.subplots(3,1)
    ax[0].semilogy(T, eT, '.', ms=1)
    ax[1].semilogy(d1, e1, '.', ms=1)
    ax[2].semilogy(d2, e2, '.', ms=1)
    
    fig,ax = plt.subplots(1,1)
    ax.plot(d[IT], T[IT], 'g^', ms=3, label='T')
    ax.plot(d[I1], T[I1], 'bx', ms=4, label='d1')
    ax.plot(d[I2], T[I2], 'o', mec='r', mfc='none', ms=4, label='d2')
    ax.plot(d[IN], T[IN], 'ms', ms=3, label='nan')
    ax.plot(d[Ibad], T[Ibad], 's', mec='c', mfc='none', ms=4, label='bad')
    
    ax.legend(loc=0)
    ax.plot(subst._sattable['dL'], subst._sattable['T'], 'k')
    ax.plot(subst._sattable['dV'], subst._sattable['T'], 'k')
    
    plt.show()
    
    out = {'IT':IT, 'I1':I1, 'I2':I2, 'IN':IN, 'T':T_, 'd1': d1_, 'd2':d2_, 'x':x_, 'I':I_, 'Ip1':Ip1, 'Ip2':Ip2}
    
    return out



def failsearch(subst, state, p1, p2):
    p1value = state[p1]
    p2value = state[p2]
    a = 0
    b = p1value.size
    failed = False
    while b-a > 1:
        c = int(0.5 * (a+b))
        print(a,b)
        try:
            subst._argparse(**{p1:p1value[a:c], p2:p2value[a:c]})
            a = c
        except:
            failed = True
            b = c
    if failed:
        return a
    return -1

def renew():
    global S
    pm.reg.regload()
    pm.dat.load()
    S = pm.get('mp.N2')

