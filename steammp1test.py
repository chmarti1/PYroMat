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
mp1obj.p(T=[300,300],d=996.74)
mp1obj.p(T=np.asarray(300),d=996.74)
mp1obj.p(T=np.asarray([300,300]),d=996.74)
mp1obj.p(T=np.asarray([300,300]),d=996.74,quality=True)
print('Pressure')
print(mp1obj.p(T=300,d=996.74))
print(mp1obj.p(T=600,d=1.8242))
print(mp1obj.p(T=424.98,x=0.5,quality=True))
print(mp1obj.p(T=800,d=83.132))

#d
mp1obj.d(T=300,p=5) 
mp1obj.d(T=[300,300],p=5)
mp1obj.d(T=np.asarray(300),p=5)
mp1obj.d(T=np.asarray([300,300]),p=5)
print('Density')
print(mp1obj.d(T=300,p=5))
print(mp1obj.d(T=600,p=5))
print(mp1obj.d(T=424.98,x=0.5))
print(mp1obj.d(T=800,p=250))

#e
mp1obj.e(T=300,d=996.74) 
mp1obj.e(T=[300,300],p=5)
mp1obj.e(T=np.asarray(300),p=5)
mp1obj.e(T=np.asarray([300,300]),p=5)
print('Energy')
print(mp1obj.e(T=300,p=5))
print(mp1obj.e(T=600,d=1.8242))
print(mp1obj.e(T=424.98,x=0.5))
print(mp1obj.e(T=800,p=250))

#h
mp1obj.h(T=300,d=996.74) 
mp1obj.h(T=[300,300],p=5)
mp1obj.h(T=np.asarray(300),p=5)
mp1obj.h(T=np.asarray([300,300]),p=5)
print('Enthalpy')
print(mp1obj.h(T=300,p=5))
print(mp1obj.h(T=600,d=1.8242))
print(mp1obj.h(T=424.98,x=0.5))
print(mp1obj.h(T=800,p=250))

#s
mp1obj.s(T=300,d=996.74) 
mp1obj.s(T=[300,300],p=5)
mp1obj.s(T=np.asarray(300),p=5)
mp1obj.s(T=np.asarray([300,300]),p=5)
print('Entropy')
print(mp1obj.s(T=300,p=5))
print(mp1obj.s(T=600,d=1.8242))
print(mp1obj.s(T=424.98,x=0.5))
print(mp1obj.s(T=800,p=250))

#hsd
mp1obj.hsd(T=300,d=996.74) 
mp1obj.hsd(T=[300,300],p=5)
mp1obj.hsd(T=np.asarray(300),p=5)
mp1obj.hsd(T=np.asarray([300,300]),p=5)
print('HSD')
print(mp1obj.hsd(T=300,p=5))
print(mp1obj.hsd(T=600,d=1.8242))
print(mp1obj.hsd(T=424.98,x=0.5))
print(mp1obj.hsd(T=800,p=250))

#cp
mp1obj.cp(T=300,d=996.74) 
mp1obj.cp(T=[300,300],p=5)
mp1obj.cp(T=np.asarray(300),p=5)
mp1obj.cp(T=np.asarray([300,300]),p=5)
print('cp')
print(mp1obj.cp(T=300,p=5))
print(mp1obj.cp(T=600,d=1.8242))
print(mp1obj.cp(T=424.98,x=0.5))
print(mp1obj.cp(T=800,p=250))

#cv
mp1obj.cv(T=300,d=996.74) 
mp1obj.cv(T=[300,300],p=5)
mp1obj.cv(T=np.asarray(300),p=5)
mp1obj.cv(T=np.asarray([300,300]),p=5)
print('cv')
print(mp1obj.cv(T=300,p=5))
print(mp1obj.cv(T=600,d=1.8242))
print(mp1obj.cv(T=424.98,x=0.5))
print(mp1obj.cv(T=800,p=250))

#T_s (has bugs)
mp1obj.T_s(p=5,s=0.39295) 
mp1obj.T_s(p=[5,5],s=0.39295)
mp1obj.T_s(p=np.asarray(5),s=0.39295)
mp1obj.T_s(p=np.asarray([5,5]),s=0.39295)
mp1obj.T_s(p=5,s=[0.39295,0.39295]) 
mp1obj.T_s(p=5,s=np.asarray(0.39295))
mp1obj.T_s(p=np.asarray(5),s=np.asarray([0.39295,0.39295]))
mp1obj.T_s(p=np.asarray(5),s=np.asarray(0.39295))
mp1obj.T_s(p=np.asarray([5,5]),s=np.asarray([0.39295,0.39295]))
print('T_s')
print(mp1obj.T_s(p=[5],s=0.39295))
print(mp1obj.T_s(p=[5],s=7.5561))
print(mp1obj.T_s(p=[5],s=4.3406,quality=True))
print(mp1obj.T_s(p=[250],s=6.0867))

#T_h (has bugs)
mp1obj.T_h(p=5,h=113.02) 
mp1obj.T_h(p=[5,5],h=113.02)
mp1obj.T_h(p=np.asarray(5),h=113.02)
mp1obj.T_h(p=np.asarray([5,5]),h=113.02)
mp1obj.T_h(p=5,h=[113.02,113.02]) 
mp1obj.T_h(p=5,h=np.asarray(113.02))
mp1obj.T_h(p=np.asarray(5),h=np.asarray([113.02,113.02]))
mp1obj.T_h(p=np.asarray(5),h=np.asarray(113.02))
mp1obj.T_h(p=np.asarray([5,5]),h=np.asarray([113.02,113.02]))
print('T_h')
print(mp1obj.T_h(p=[5],h=113.02))
print(mp1obj.T_h(p=[5],h=3120.1))
print(mp1obj.T_h(p=[5],h=1694.095,quality=True))
print(mp1obj.T_h(p=[250],h=3262.2))

# Performing computations at ref A
#print(mp1obj.T(p=5,d=996.74)) #error
