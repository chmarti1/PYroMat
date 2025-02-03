#!/usr/bin/python3

import pyromat as pm
import numpy as np
import os,sys
import matplotlib.pyplot as plt
import json

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
    """Establish points along the liquid-vapor saturation line
    T,p,dL,dV = satline(subst, xtrans=.25, dstep=.01, nstep=100)
    
The saturation line is calculated iteratively using Maxwell's criteria.
In general, the iteration is unstable unless the initial guess is close
to the correct value, so the iteration progresses along the saturation
line in small steps from the previous valid solution; starting at the
critical point.

Very close to the critical point, changes in the saturation temperature
are very nearly zero, so small changes in density are made instead.

*Arguments*
subst   (mandatory)
The multi-phase substance instance.

xtrans  (optional)
Near the critical point, 
"""
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
        Ts,ps,dsL,dsV = satiterL(subst,Ts,dsL,dsV)
        dsL_array.append(dsL)
        dsV_array.append(dsV)
        Ts_array.append(Ts)
        ps_array.append(ps)
    Ttrans = Ts
    for Ts in np.linspace(Ttrans, Tt, ntstep):
        Ts,ps,dsL,dsV = satiterT(subst,Ts,dsL,dsV)
        dsL_array.append(dsL)
        dsV_array.append(dsV)
        Ts_array.append(Ts)
        ps_array.append(ps)
    
    return np.flip(Ts_array), np.flip(ps_array), np.flip(dsL_array), np.flip(dsV_array)

if __name__ == '__main__':
    targets = pm.search(pmclass='mp2')
    for this in targets:
        Ts,ps,dsL,dsV = satline(this, xtrans=0.1, dstep=.01, ntstep=100)
        
        pt = ps[0]
        print('Tt:', Ts[0])
        print('pt:', pt)
        
        
        fig,ax = plt.subplots(1,1)
        ax.set_title(this.hill())
        ax.plot(dsL,Ts,'b.')
        ax.plot(dsV,Ts,'r.')
        plt.show()
        
        ui = input('Accept? (y/n):')
        if ui == 'y':
            data = this.data
            data['pt'] = pt
            data['SATGroup'] = {
                'Ts': list(Ts),
                'ps': list(ps),
                'dsL': list(dsL),
                'dsV': list(dsV)
            }
            fromfile = data['fromfile']
            with open(fromfile,'w') as fd:
                json.dump(data, fd, indent=2)
            
