import pyromat as pm
import numpy as np
import matplotlib.pylab as pylab

#get steam
mp1obj = pm.get('mp.H2O')

#Public methods
#Tlim
mp1obj.Tlim()
mp1obj.Tlim(p=1)

#Plim
mp1obj.plim()
mp1obj.plim(T=1000)

#critical
mp1obj.critical()
mp1obj.critical(density=True)

#triple
mp1obj.triple()

#ps
mp1obj.ps(T=373) #should be around 1 bar
mp1obj.ps(T=[373,373])
mp1obj.ps(T=np.asarray(373))
mp1obj.ps(T=np.asarray([373,373]))
#print(mp1obj.ps(T=373))

#Ts
mp1obj.Ts(p=1) #should be around 1 bar
mp1obj.Ts(p=[1,1])
mp1obj.Ts(p=np.asarray(1))
mp1obj.Ts(p=np.asarray([1,1]))

#ds
mp1obj.ds(p=1)
mp1obj.ds(p=[1,1])
mp1obj.ds(p=np.asarray(1))
mp1obj.ds(p=np.asarray([1,1]))
mp1obj.ds(T=373)
mp1obj.ds(T=[373,373])
mp1obj.ds(T=np.asarray(373))
mp1obj.ds(T=np.asarray([373,373]))

#es
mp1obj.es(p=1)
mp1obj.es(p=[1,1])
mp1obj.es(p=np.asarray(1))
mp1obj.es(p=np.asarray([1,1]))
mp1obj.es(T=373)
mp1obj.es(T=[373,373])
mp1obj.es(T=np.asarray(373))
mp1obj.es(T=np.asarray([373,373]))

#hs
mp1obj.hs(p=1)
mp1obj.hs(p=[1,1])
mp1obj.hs(p=np.asarray(1))
mp1obj.hs(p=np.asarray([1,1]))
mp1obj.hs(T=373)
mp1obj.hs(T=[373,373])
mp1obj.hs(T=np.asarray(373))
mp1obj.hs(T=np.asarray([373,373]))

#ss
mp1obj.ss(p=1)
mp1obj.ss(p=[1,1])
mp1obj.ss(p=np.asarray(1))
mp1obj.ss(p=np.asarray([1,1]))
mp1obj.ss(T=373)
mp1obj.ss(T=[373,373])
mp1obj.ss(T=np.asarray(373))
mp1obj.ss(T=np.asarray([373,373]))

#Four reference points for H2O: A - liquid, B - vapor, C - mixture, D - supercritical (NIST Webbook)
#A: T = 300 K, P = 5 bar, d = 996.74 kg/m3, v = 0.0010033 m3/kg, h = 113.02 kJ/kg, s = 0.39295 kJ/kg K
#B: T = 600 K, P = 5 bar, d = 1.8242 kg/m3, v = 0.54820 m3/kg, h = 3120.1 kJ/kg, s = 7.5561 kJ/kg K
#C: T = 424.98 K, P = 5 bar, d = 5.321 kg/m3, v = 0.18795 m3/kg, h = 1694.095 kJ/kg, s = 4.3406 kJ/kgK, x = 0.5
#D: T = 800 K, P = 250 bar, d = 83.132 kg/m3, v = 0.012029 m3/kg, h = 3262.2 kJ/kg, s = 6.0867 kJ/kgK


#p
mp1obj.p(T=300,d=996.74) 
mp1obj.p(T=[373,373],d=996.74)
mp1obj.p(T=np.asarray(373),d=996.74)
mp1obj.p(T=np.asarray([373,373]),d=996.74)
print(mp1obj.p(T=300,d=996.74))
print(mp1obj.p(T=600,d=1.8242))
print(mp1obj.p(T=424.98,x=0.5))
print(mp1obj.p(T=800,d=83.132))

# Performing computations at ref A
#print(mp1obj.T(p=5,d=996.74)) #error
