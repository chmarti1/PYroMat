#!../../bin/python

import pyromat as pm
import numpy as np


substances = pm.search(pmclass='mp2')

for subst in substances:
    print(subst)
    Tmin,Tmax = subst.Tlim()
    _,pmax = subst.plim()
    dmax = subst.data['dlim'][1]
    print('  Registered maximum density:', dmax)
    dmax = subst.d(T=Tmin, p=pmax)[0]
    print('  Calculated maximum density:', dmax)
    dmax = np.ceil(dmax)
    response = None
    while response not in {'Y', 'n'}:
        response = input(f'  Update to {dmax}? (Y/n):')
    if response == 'Y':
        subst.data['dlim'][1] = dmax
    
pm.dat.updatefiles()
