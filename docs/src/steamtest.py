#
#   Steam phase change demo
#   By Chris Martin (c) 2016
#   GPL v3.0
#   Enjoy!
#
# In this code, we generate a surface plot that spans the saturation
# curve from triple point to critical point.  We'll add two red lines
# to show where the liquid- and gas-phase saturation properties are in 
# each plot.  

import pyromat as pyro
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


# Start with a blank slate
plt.close('all')

# Get the steam object
S = pyro.get('mp.H2O')

# Get the critical and triple point properties
Tt,pt = S.triple()
Tc,pc = S.critical()

# Explore the temperatures between Tt and Tc in 5K increments
Ts = np.arange(Tt,Tc,5.)

# Now, obtain the saturation properties
# Note that when a saturation property is called with the "tp" flag True
# it also returns the saturation temperature and pressure.
# This is more efficient than separate calls to Ts and ps.
Ts,ps,dL,dV = S.ds(T=Ts,tp=True)
# It is faster to explicitly pass both Ts and ps to the saturation methods.
# This saves them a redundant call to Ts or ps since both are required anyway.
hL,hV = S.hs(T=Ts,p=ps)
sL,sV = S.ss(T=Ts,p=ps)
eL,eV = S.es(T=Ts,p=ps)
# Note that this code results in warnings that the accuracy of saturation
# properties is reduced above 623.15K.  This is normal; PYro is simply 
# making sure the user is aware that saturation properties in this region 
# can be difficult to provide accurately.

# Now, crank out some surface plots
T,p = np.meshgrid(Ts,ps)

# Calculate density, enthalpy, entropy, and internal energy
d = S.d(T,p)
h = S.h(T,p)
s = S.s(T,p)
e = S.e(T,p)


f = plt.figure(1)
ax = f.add_subplot(111,projection='3d')
ax.plot_surface(T,p,d)
ax.plot(Ts,ps,dL,'r')
ax.plot(Ts,ps,dV,'r')
ax.set_xlabel('T')
ax.set_ylabel('p')
ax.set_zlabel('density')

f = plt.figure(2)
ax = f.add_subplot(111,projection='3d')
ax.plot_surface(T,p,h)
ax.plot(Ts,ps,hL,'r')
ax.plot(Ts,ps,hV,'r')
ax.set_xlabel('T')
ax.set_ylabel('p')
ax.set_zlabel('enthalpy')

f = plt.figure(3)
ax = f.add_subplot(111,projection='3d')
ax.plot_surface(T,p,s)
ax.plot(Ts,ps,sL,'r')
ax.plot(Ts,ps,sV,'r')
ax.set_xlabel('T')
ax.set_ylabel('p')
ax.set_zlabel('entropy')

f = plt.figure(4)
ax = f.add_subplot(111,projection='3d')
ax.plot_surface(T,p,e)
ax.plot(Ts,ps,eL,'r')
ax.plot(Ts,ps,eV,'r')
ax.set_xlabel('T')
ax.set_ylabel('p')
ax.set_zlabel('energy')

plt.show(block=False)
