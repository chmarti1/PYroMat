import pyromat as pyro
N2 = pyro.get('N2')
T_h = pyro.solve.solve1n('T',f=N2.h, df=N2.cp, param_lim=(N2.data['Tmin'],N2.data['Tmax']))
T_h(100.)
