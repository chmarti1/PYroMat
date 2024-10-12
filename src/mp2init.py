#!/usr/bin/python3

import pyromat as pm
import numpy as np

def satiterT(subst, T, dL, dV):
    for count in range(50):
        gL,gLT,gLd = subst._g(T=T, d=dL, diff=1)
        gV,gVT,gVd = subst._g(T=T, d=dV, diff=1)
        pL,pLT,pLd = subst._p(T=T, d=dL, diff=1)
        pV,pVT,pVd = subst._p(T=T, d=dV, diff=1)

        e = np.array([gL - gV, pL - pV])
        J = np.array([[gLd, -gVd], [pLd, -pVd]])
        
        dd = np.linalg.solve(J,e)
        
        if np.abs(dd[0]) < dL*1e-6 and np.abs(dd[1]) < dV*1e-6:
            return T, pV, dL, dV
            
        dL -= dd[0]
        dV -= dd[1]
    raise Exception('Failed to converge.')
            
def satiterL(subst, T, dL, dV):
    for count in range(50):
        print("  ", T, dV)
        gL,gLT,gLd = subst._g(T=T, d=dL, diff=1)
        gV,gVT,gVd = subst._g(T=T, d=dV, diff=1)
        pL,pLT,pLd = subst._p(T=T, d=dL, diff=1)
        pV,pVT,pVd = subst._p(T=T, d=dV, diff=1)

        e = np.array([gL - gV, pL - pV])
        J = np.array([[gLT-gVT, -gVd], [pLT-pVT, -pVd]])
        
        dd = np.linalg.solve(J,e)
        
        if np.abs(dd[0]) < T*1e-6 and np.abs(dd[1]) < dV*1e-6:
            return T, pV, dL, dV

        T -= dd[0]
        dV -= dd[1]
    raise Exception('Failed to converge.')

def satline(subst, xtrans=0.25, dstep = 0.01, ntstep=100):
    Tt = subst.data['Tt']
    Tc = subst.data['Tc']
    dc = subst.data['dc']

    # Initialize the saturation points
    dsL = dc
    dsV = dc
    Ts = Tc
    ps = subst._p(T=Tc, d=dc)[0]

    dsL_array = [dc]
    dsV_array = [dc]
    Ts_array = [Tc]
    ps_array = [ps]
    
    # Incrementally increase the saturated liquid density until the 
    # temperature is 25% of the way to the triple point
    Ttrans = xtrans * Tt + (1-xtrans) * Tc
    while Ts > Ttrans:
        # Increase dsL
        dsL *= (1+dstep)
        dsV *= (1-dstep)
        print(Ts)
        Ts,ps,dsL,dsV = satiterL(subst,Ts,dsL,dsV)
        dsL_array.append(dsL)
        dsV_array.append(dsV)
        Ts_array.append(Ts)
        ps_array.append(ps)
    Ttrans = Ts
    for Ts in np.linspace(Ttrans, Tt, ntstep):
        print("*", Ts)
        Ts,ps,dsL,dsV = satiterT(subst,Ts,dsL,dsV)
        dsL_array.append(dsL)
        dsV_array.append(dsV)
        Ts_array.append(Ts)
        ps_array.append(ps)
    
    return np.flip(Ts_array), np.flip(ps_array), np.flip(dsL_array), np.flip(dsV_array)

if __name__ == '__main__':
    
