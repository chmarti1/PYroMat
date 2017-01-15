#
#   Steam phase change demo
#   By Chris Martin (c) 2016
#   GPL v3.0
#   Enjoy!
#

import pyromat as pyro
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Start with a blank slate
plt.close('all')

# Call the steam object
S = pyro.get('steam')
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
NN = T.shape
# Reshape T and p to be 1D vectors
# As of v1.2, the steam class chokes on higher dimensional vectors
p.shape = (T.size,)
T.shape = (T.size,)

# Calculate density, enthalpy, entropy, and internal energy
d = S.d(T,p)
h = S.h(T,p)
s = S.s(T,p)
e = S.e(T,p)

# Reshape the property vectors for ploting
p.shape = NN
T.shape = NN
d.shape = NN
h.shape = NN
s.shape = NN
e.shape = NN

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
