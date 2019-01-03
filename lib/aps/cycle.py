"""PYroMat thermodynamic cycle solvers

"""

import pyromat as pm
import numpy as np
import matplotlib.pyplot as plt


class PMCycleError(Exception):
    pass



class Cycle(object):
    """The base Cycle class defines members:
    param       a dictionary of cycle parameters
    meta        a dictionary of cycle analysis results

The following members describe the fluid properties at the discrete
cycle states.  They are lists with an element for each state.
    T   Temperature
    p   Pressure
    d   Density
    x   Quality (if applicable)
    s   Entropy
    h   Enthalpy

The following members describe the processes that separate the states. 
They are lists with a member for each process.  They are ordered so that
the 0-element of each list corresponds to the process between states 0 
and 1.
    q   Heat per unit matter added during the process
    w   Work extracted per unit mater during the process
    
The __repr__() method constructs a printed summary of the parameters, 
states, and processes.  It is up to the individual cycle classes to 
enforce that the states and processes are lists of dictionaries of the 
same length.  Each state dictionary contains keyword and value pairs
identifying properties and their values.  Each process dictionary should
contain 'w' and 'q' keywords with corresponding values for work and heat
respectively.

The _test() method is a utility for conducting basic data integrity
checks.
"""
    def __init__(self, N=None):
        if N is None:
            N = 0
        self.param = {}
        self.meta = {}
        self.T = [0] * N
        self.p = [0] * N
        self.d = [0] * N
        self.x = [0] * N
        self.s = [0] * N
        self.h = [0] * N
        self.w = [0] * N
        self.q = [0] * N
        
    def __repr__(self):
        if self._test() or len(self.T)<1:
            return '<Incomplete ' + str(self.__class__.__name__) + '>'
        
        # What parameters?
        out = str(self.__class__.__name__) + '\nParameters:\n'
        kk = self.param.keys()
        kk.sort()
        for this in kk:
            out += '  ' + this + ' = ' + repr(self.param[this]) + '\n'
        
        if isinstance(self.T[0],np.ndarray) and self.T[0].size>1:
            out += '[Performance data are arrays]\n'
            return out
        
        # 
        out += 'State Table\n'
        out += '    %12s %12s %12s %12s %12s %12s\n'%('T', 'p', 'd', 'x', 'h', 's')
        for index in range(len(self.T)):
            out += '%3d %12f %12f %12f %12f %12f %12f\n'%(\
                    index+1, self.T[index], self.p[index], self.d[index],
                    self.x[index], self.h[index], self.s[index])
        
        out += '\nProcess Table\n'
        out += '      %12s %12s\n'%('q', 'w')
        for index in range(len(self.w)):
            out += '%3d-%1d %12f %12f\n'%(index+1, index+2, self.q[index], self.w[index])
        return out
        
        
    def _test(self, N=None, require=None):
        """Returns a non-zero if there is a problem with the cycle data.
    Cycle.test()
    
This is a basic data integrity test.  When none of the optional 
parameters are passed, _test() only ensures that the number of states 
matches the number of processes, and that all process and state dicts
have compatible keyword elements.

Optional keyword arguments:

N           [positive integer]  
Number of states and processes required by the cycle.  If this test 
fails, 

require     [tuple or list of strings]
A tuple or list of string parameter names to require.  Return a code if 
these required string keywords are not found in the params dictionary.

Return codes:
0:  There were no errors
1:  The state properties or process parameter lengths do not match
3:  The require keywords are not found in params
-1: There was an unhandled exception either due to the parameters passed
        to _test() or due to an unhandled case in the object's data
"""
        # If N is not given, learn it from the number of states
        if N is None:
            N = len(self.T)
            
        if  len(self.T) != N or\
                len(self.p) != N or\
                len(self.d) != N or\
                len(self.x) != N or\
                len(self.s) != N or\
                len(self.h) != N or\
                len(self.w) != N:
            return 1
        
        if require is not None:
            for this in require:
                if this not in self.param:
                    return 3
            
        return 0
        
    def _initplot(self, xlabel, ylabel, ax=None, fig=None):
        """Initialize a figure axes
    ax = _initplot(xlabel, ylabel, ax=None, fig=None)
    
This helper funciton processes the axes and figure optional parameters
for the cycle plotting functions.  It returns a properly formatted axes
for the plot.  By default, it creates a new figure with a blank axes and
returns the axes object.

If the AX parameter is specified, then it is formatted and returned.
If the FIG parameter is specified, it is interpreted as a figure object.
If that figure has current axes, the first one is selected and treated
as the plot axes.  If the figure has no axes, then axes are created and
those will be returned.
"""
        if ax is None:
            if fig is None:
                fig = plt.figure()
            # Recover a list of all axes in the active figure
            ax = fig.get_axes()
            # If it is not empty, use the first in the list
            if ax:
                ax = ax[0]
            else:
                ax = fig.add_subplot(111)

        # Now, we have an axis; format it.
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True)
        
        return ax
        
        
        
class BraytonCycle(Cycle):
    """The Brayton Cycle
"""
    def __init__(self):
        # Call the super init
        Cycle.__init__(self, N=4)
        # Initialize the parameters to generally safe values
        self.param = {
                'p1':1.01325, 
                'p2':20, 
                'fluid':'ig.air', 
                'T1':300., 
                'T3':2000.}
        

    def update(self):
        """Update the cycle states and processes
"""
        if self._test(require=('fluid', 'p1', 'p2', 'T1', 'T3')):
            raise PMCycleError('There was a problem with the Brayton Cycle configuration')

        fluid = self.param['fluid']
        p1 = self.param['p1']
        p2 = self.param['p2']
        T1 = self.param['T1']
        T3 = self.param['T3']

        # Parse the subsance input
        if isinstance(fluid, str):
            fluid = pm.get(fluid)

        # Format the pressure inputs
        p1 = np.asarray(p1,dtype=float)
        if p1.ndim==0:
            p1 = np.reshape(p1, (1,))
        
        p2 = np.asarray(p2,dtype=float)
        if p2.ndim==0:
            p2 = np.reshape(p2, (1,))
            
        T1 = np.asarray(T1,dtype=float)
        if T1.ndim==0:
            T1 = np.reshape(T1, (1,))
        
        T3 = np.asarray(T3,dtype=float)
        if T3.ndim==0:
            T3 = np.reshape(T3, (1,))
        
        p1,p2,T1,T3 = np.broadcast_arrays(p1,p2,T1,T3)
            
        # p1 < p2 < pc
        if (p2 <= p1).any():
            raise PMCycleError('BraytonCycle requires p2 to be greater than p1.')
            
        if (T3 <= T1).any():
            raise PMCycleError('BraytonCycle requires T3 to be less than T1.')
        
        
        # Calculate state 1
        self.p[0] = p1
        self.T[0] = T1
        self.s[0] = fluid.s(T=T1,p=p1)
        self.h[0] = fluid.h(T=T1,p=p1)
        self.d[0] = fluid.d(T=T1,p=p1)

        # Calculate state 2
        self.p[1] = p2
        self.s[1] = self.s[0]
        self.T[1] = fluid.T_s(s=self.s[0],p=p2)
        self.d[1] = fluid.d(T=self.T[1],p=p2)
        self.h[1] = fluid.h(T=self.T[1],p=p2)
        
        # Calculate state 3
        self.p[2] = p2
        self.T[2] = T3
        self.s[2] = fluid.s(T=T3, p=p2)
        self.h[2] = fluid.h(T=T3, p=p2)
        self.d[2] = fluid.d(T=T3, p=p2)
        
        # Calculate state 4
        self.p[3] = p1
        self.s[3] = self.s[2]
        self.T[3] = fluid.T_s(s=self.s[3], p=p1)
        self.h[3] = fluid.h(T=self.T[3], p=p1)
        self.d[3] = fluid.d(T=self.d[3], p=p1)
        
        self.w = [
                self.h[0] - self.h[1],
                0.,
                self.h[2] - self.h[3],
                0.]
        
        self.q = [
                0.,
                self.h[2] - self.h[1],
                0.,
                self.h[0] - self.h[3]]
        
        self.meta['eta'] = (self.w[2]+self.w[0]) / self.q[1]
    
    
class RankineCycle(Cycle):
    """The Rankine cycle, also known as the steam cycle, once powered the world's
trains and factories, and produces most of the world's electricity today.
Steam is brought to a boil under pressure and expanded through a turbine
or a piston to produce shaft work.  In the basic Rankine cycle, sat-
urated steam is drawn directly from a boiler into the turbine or piston.
All modern systems use at least one additional super-heat process to 
push the steam into the super-heated region.  To model a super-heated
Rankine cycle, see the RankineSHCycle class.

Parameters:
    p1          The reservoir pressure
    p2          The boiler pressure
    fluid       The working fluid

The basic Rankine cycle consists of four processes and four states.

(State 1)       [Liquid Reservoir at pressure p1]

(Process 1-2)   A feed-water-pump performs approximately isentropic 
                compression of liquid water.  Treating this as an 
                isothermal pressure increase is a practical but 
                technically imprecise approximation.
                
(State 2)       [Boiler inlet at p2]
                
(Process 2-3)   A boiler produces saturated steam under pressure.

(State 3)       [Boiler outlet at p2]

(Process 3-4)   Saturated steam is expanded isentropically to extract
                work.
                
(State 4)       [Condenser inlet at p1]
                
(Process 4-1)   The liquid-vapor mixture is condensed back into a sat-
                urated liquid by extracting waste heat from a condenser.

"""

    def __init__(self):
        # Call the super init
        Cycle.__init__(self, N=4)
        # Initialize the parameters to generally safe values
        self.param = {
                'p1':1.01325, 
                'p2':10, 
                'fluid':'mp.H2O',
                'eta12':1.,
                'eta23':1.,
                'eta34':1.}
        


    def update(self):
        """Update the cycle states and processes
"""
        if self._test(require=('fluid', 'p1', 'p2', 'eta12', 'eta23', 'eta34')):
            raise PMCycleError('There was a problem with the Rankine Cycle configuration')

        fluid = self.param['fluid']
        p1 = self.param['p1']
        p2 = self.param['p2']
        eta12 = self.param['eta12']
        eta23 = self.param['eta23']
        eta34 = self.param['eta34']
        

        # Parse the subsance input
        if isinstance(fluid, str):
            fluid = pm.get(fluid)
        if not isinstance(fluid, pm.reg.registry['mp1']):
            raise PMCycleError('RankineCycle requires substances of class mp1')

        # Format the pressure inputs
        p1 = np.asarray(p1,dtype=float)
        if p1.ndim==0:
            p1 = np.reshape(p1, (1,))
        
        p2 = np.asarray(p2,dtype=float)
        if p2.ndim==0:
            p2 = np.reshape(p2, (1,))
        
        eta12 = np.asarray(eta12, dtype=float)
        if eta12.ndim==0:
            eta12 = np.reshape(eta12, (1,))
            
        eta23 = np.asarray(eta23, dtype=float)
        if eta23.ndim==0:
            eta23 = np.reshape(eta23, (1,))
            
        eta34 = np.asarray(eta34, dtype=float)
        if eta12.ndim==0:
            eta34 = np.reshape(eta34, (1,))
        
        p1,p2,eta12,eta23,eta34 = np.broadcast_arrays(p1,p2,eta12,eta23,eta34)
        
            
        # Get the critical and triple points
        Tc,pc = fluid.critical()
        Tt,pt = fluid.triple()
        
        # Test the inputs
        # pt < p1 < pc
        if (p1 <= pt).any():
            raise PMCycleError('RankineCycle requires p1 to be greater than the triple-point pressure')
        elif (p1 >= pc).any():
            raise PMCycleError('Rankine Cycle requires p1 to be less than the critical-point pressure')
            
        # p1 < p2 < pc
        if (p2 <= p1).any():
            raise PMCycleError('RankineCycle requires p2 to be greater than p1.')
        elif (p2 >= pc).any():
            raise PMCycleerror('RankineCycle requires p1 to be less than the critical-point pressure')
            
        
        # Calculate state 1
        self.p[0] = p1
        self.T[0] = fluid.Ts(p=p1)
        self.s[0],_ = fluid.ss(T=self.T[0])
        self.h[0],_ = fluid.hs(T=self.T[0])
        self.d[0],_ = fluid.ds(T=self.T[0])
        self.x[0] = 0.

        # Calculate state 2
        self.p[1] = p2
        # Start with the isentropic performance
        T2s = fluid.T_s(s=self.s[0], p=p2)
        h2s = fluid.h(T=T2s, p=p2)
        # Adjust the enthalpy rise by the inverse of pump efficiency
        self.h[1] = self.h[0] + (h2s - self.h[0])/eta12
        self.T[1],self.x[1] = fluid.T_h(h=self.h[1],p=p2,quality=True)
        self.d[1] = fluid.d(T=self.T[1], p=p2, x=self.x[1])
        # Now that we have density, we can use it directly to solve for entropy
        self.s[1] = fluid.s(T=self.T[1], d=self.d[1])
        
        # Calculate state 3
        self.p[2] = p2
        self.T[2] = fluid.Ts(p=p2)
        _,self.s[2] = fluid.ss(T=self.T[2])
        _,self.h[2] = fluid.hs(T=self.T[2])
        _,self.d[2] = fluid.ds(T=self.T[2])
        self.x[2] = 1.
        
        # Calculate state 4
        self.p[3] = p1
        T4s,x4s = fluid.T_s(self.s[2], p=p1, quality=True)
        h4s = fluid.h(T=T4s, x=x4s)
        # Adjust enthlpy fall by the turbine efficiency
        self.h[3] = self.h[2] + (h4s - self.h[2])*eta34
        self.T[3],self.x[3] = fluid.T_h(h=self.h[3], p=p1, quality=True)
        self.d[3] = fluid.d(T=self.T[3], p=p1, x=self.x[3])
        # Now that we have density, we can use it directly to solve for entropy
        self.s[3] = fluid.s(T=self.T[3], d=self.d[3])
        
        self.w = [
                self.h[0] - self.h[1],
                0.,
                self.h[2] - self.h[3],
                0.]
        
        # Adjust the required boiler heat by the efficiency
        self.q = [
                0.,
                (self.h[2] - self.h[1]) / eta23,
                0.,
                self.h[0] - self.h[3]]
        
        self.meta['eta'] = (self.w[2]+self.w[0]) / self.q[1]
        self.meta['uE'] = pm.config['unit_energy']
        self.meta['uT'] = pm.config['unit_temperature']
        self.meta['uP'] = pm.config['unit_pressure']
        self.meta['uM'] = pm.config['unit_matter']
        self.meta['uV'] = pm.config['unit_volume']
        self.meta['fluid'] = fluid
    
    
    
    def tsplot(self, ax=None, fig=None, slabels=True, 
            satstyle=None, procstyle=None, statestyle=None, boxstyle=None):
        """
"""
        # If the properties are empty, abort
        if not self.T:
            raise PMCycleError('The plot cannot be generated for a cycle that has not yet been updated.  Run update().')
        if isinstance(self.T, np.ndarray) and self.T.size>1:
            raise PMCycleError('Plotting arrays of cycle data is not supported.')
        
        ax = self._initplot(
                's (%s/%s%s)'%(self.meta['uE'], self.meta['uM'], self.meta['uT']), 
                'T (%s)'%(self.meta['uT']), ax=ax, fig=fig)
        
        if satstyle is None:
            satstyle = {'lw':2, 'c':'k', 'ls':'solid', 'marker':'None'}
        if procstyle is None:
            procstyle = {'lw':2, 'c':'r', 'ls':'solid', 'marker':'None'}
        if statestyle is None:
            statestyle = {'ls':'None', 'marker':'o', 'mec':'k', 'mew':1, 'mfc':'r', 'ms':6}
        if boxstyle is None:
            boxstyle = {'fc':'w', 'ec':'k', 'boxstyle':'square,pad=.25'}
            
        # Get the working fluid object
        fluid = self.meta['fluid']
            
        # Plot the dome
        Tc,pc = fluid.critical()
        Tt,pt = fluid.triple()
        T = np.linspace(Tt, 0.99999*Tc, 100)
        sL, sV = fluid.ss(T=T)
        ax.plot(sL, T, **satstyle)
        ax.plot(sV, T, **satstyle)
            
        # Polytropic pump
        T = np.linspace(self.T[0], self.T[1], 21)
        d = np.linspace(self.d[0], self.d[1], 21)
        s = fluid.s(T=T,d=d)
        ax.plot( s, T, **procstyle)
        
        # Isobaric boiler
        s = np.linspace(self.s[1], self.s[2], 51)
        T = fluid.T_s(s, p=self.p[1])
        ax.plot(s, T, **procstyle)
        
        # Polytropic turbine
        # Use a linear approximation
        ax.plot(
                [self.s[2], self.s[3]],
                [self.T[2], self.T[3]], 
                **procstyle)
        
        # Isobaric condensation
        s = np.linspace(self.s[3], self.s[0], 51)
        T = fluid.T_s(s, p=self.p[0])
        ax.plot(s, T, **procstyle)
        
        # deterine offsets in T and s space
        ds = ax.get_xlim()
        ds = .02 * (ds[1] - ds[0])
        dT = ax.get_ylim()
        dT = .02 * (dT[1] - dT[0])
        
        # Deal with the state labels
        ax.plot(self.s[0], self.T[0], **statestyle)
        tt = ax.text(self.s[0]-ds, self.T[0]-dT, '1', ha='right', va='top')
        tt.set_bbox(boxstyle)
        ax.plot(self.s[1], self.T[1], **statestyle)
        tt = ax.text(self.s[1]-ds, self.T[1]+dT, '2', ha='right', va='bottom')
        tt.set_bbox(boxstyle)
        ax.plot(self.s[2], self.T[2], **statestyle)
        tt = ax.text(self.s[2]+ds, self.T[2]+dT, '3', ha='left', va='bottom')
        tt.set_bbox(boxstyle)
        ax.plot(self.s[3], self.T[3], **statestyle)
        tt = ax.text(self.s[3]+ds, self.T[3]-dT, '4', ha='left', va='top')
        tt.set_bbox(boxstyle)
        
        return ax
    
    
class RankineSHCycle(Cycle):
    """The Rankine cycle, also known as the steam cycle, once powered the world's
trains and factories, and produces most of the world's electricity today.
Steam is brought to a boil under pressure and expanded through a turbine
or a piston to produce shaft work.  In the Super-Heat Rankine cycle, 
saturated steam is drawn from a boiler and further heated to form super-
heated steam.  Unlike the basic Rankine cycle, the turbine or piston
receives super-heated steam with the goal of increasing power, improving
efficiency, and mitigating the problems that come from condensing liquid 
in a turbine or piston.  Of course, this comes at the cost of higher 
operating temperatures, which stress the superheater materials.

Parameters:
    p1          The reservoir pressure
    p2          The boiler pressure
    T4          Super-heat temperature (optional)
    d34         Super-heater heat (optional)
    fluid       The working fluid

The Rankine Super-Heated Cycle consists of five processes and states.

(State 1)       [Liquid Reservoir at pressure p1]

(Process 1-2)   A feed-water-pump performs approximately isentropic 
                compression of liquid water.  Treating this as an 
                isothermal pressure increase is a practical but 
                technically imprecise approximation.
                
(State 2)       [Boiler inlet at p2]
                
(Process 2-3)   A boiler produces saturated steam under pressure.

(State 3)       [Saturated Boiler outlet at p2]

(Process 3-4)   Heat is added in a superheater.

(State 4)       [Superheated turbine inlet at p2]

(Process 4-5)   Saturated steam is expanded isentropically to extract
                work.
                
(State 5)       [Condenser inlet at p1]
                
(Process 5-1)   The exhaust of the piston/turbine is condensed back into
                a saturated liquid by extracting waste heat from a 
                condenser.

In addition to the reservoir and boiler pressures p1 and p2, RankineSH
requires some means of specifying the heat added by the superheater.
"""
    def __init__(self):
        # Call the super init
        Cycle.__init__(self, N=5)
        # Initialize the parameters to generally safe values
        self.param = {
                'p1':1.01325, 
                'p2':10., 
                'fluid':'mp.H2O'}
        

    def update(self):
        """Update the cycle states and processes
"""
        if self._test(require=('fluid', 'p1', 'p2')):
            raise PMCycleError('There was a problem with the Rankine Cycle configuration')

        fluid = self.param['fluid']
        p1 = self.param['p1']
        p2 = self.param['p2']
        q34 = 0
        T4 = 0
        mode = 0
        if 'q34' in self.param:
            q34 = self.param['q34']
            mode = 2
        elif 'T4' in self.param:
            T4 = self.param['T4']
            mode = 1

        # Parse the subsance input
        if isinstance(fluid, str):
            fluid = pm.get(fluid)
        if not isinstance(fluid, pm.reg.registry['mp1']):
            raise PMCycleError('RankineCycle requires substances of class mp1')

        # Format the pressure inputs
        p1 = np.asarray(p1,dtype=float)
        if p1.ndim==0:
            p1 = np.reshape(p1, (1,))
        
        p2 = np.asarray(p2,dtype=float)
        if p2.ndim==0:
            p2 = np.reshape(p2, (1,))
        
        T4 = np.asarray(T4,dtype=float)
        if T4.ndim==0:
            T4 = np.reshape(T4, (1,))
        
        q34 = np.asarray(q34,dtype=float)
        if q34.ndim==0:
            q34 = np.reshape(q34, (1,))
        
        p1,p2,T4,q34 = np.broadcast_arrays(p1,p2,T4,q34)
            
            
        # Get the critical and triple points
        Tc,pc = fluid.critical()
        Tt,pt = fluid.triple()
        
        # Test the inputs
        # pt < p1 < pc
        if (p1 <= pt).any():
            raise PMCycleError('RankineCycle requires p1 to be greater than the triple-point pressure')
        elif (p1 >= pc).any():
            raise PMCycleError('Rankine Cycle requires p1 to be less than the critical-point pressure')
            
        # p1 < p2 < pc
        if (p2 <= p1).any():
            raise PMCycleError('RankineCycle requires p2 to be greater than p1.')
        elif (p2 >= pc).any():
            raise PMCycleerror('RankineCycle requires p1 to be less than the critical-point pressure')
            
        
        # Calculate state 1
        self.p[0] = p1
        self.T[0] = fluid.Ts(p=p1)
        self.s[0],_ = fluid.ss(T=self.T[0])
        self.h[0],_ = fluid.hs(T=self.T[0])
        self.d[0],_ = fluid.ds(T=self.T[0])
        self.x[0] = 0.

        # Calculate state 2
        self.p[1] = p2
        self.s[1] = self.s[0]
        self.T[1] = fluid.T_s(s=self.s[0],p=p2)
        self.d[1] = fluid.d(T=self.T[1],p=p2)
        self.h[1] = fluid.h(T=self.T[1], d=self.d[1])
        self.x[1] = -1.
        
        # Calculate state 3
        self.p[2] = p2
        self.T[2] = fluid.Ts(p=p2)
        _,self.s[2] = fluid.ss(T=self.T[2])
        _,self.h[2] = fluid.hs(T=self.T[2])
        _,self.d[2] = fluid.ds(T=self.T[2])
        self.x[2] = 1.
        
        # Calculate states 4 and 5
        # If state 5 is saturated vapor
        if mode == 0:
            self.p[4] = p1
            self.T[4] = fluid.Ts(p=p1)
            _,self.d[4] = fluid.ds(T=self.T[4])
            _,self.h[4] = fluid.hs(T=self.T[4])
            _,self.s[4] = fluid.ss(T=self.T[4])
            self.x[4] = 1.
            
            self.p[3] = p2
            self.s[3] = self.s[4]
            self.T[3] = fluid.T_s(s=self.s[4], p=p2)
            self.d[3] = fluid.d(T=self.T[3], p=p2)
            self.h[3] = fluid.h(T=self.T[3],d=self.d[3])
            self.x[3] = -1.
                
        # If state 4 is temperature limited
        elif mode == 1:
            # Test the temperature
            if (T4 < self.T[2]).any():
                raise PMCycleError('RankineSHCycle requires T4 to be greater than T3.')
            
            self.p[3] = p2
            self.T[3] = T4
            self.h[3], self.s[3], self.d[3] = fluid.hsd(T=T4, p=p2)
            self.x[3] = -1.
            
            # Isentropic expansion
            self.p[4] = p1
            self.T[4], self.x[4] = fluid.T_s(s=self.s[3], p=p1, quality=True)
            # Initialize h,s, and d results
            self.h[4] = np.zeros_like(self.T[4])
            self.s[4] = np.zeros_like(self.T[4])
            self.d[4] = np.zeros_like(self.T[4])
            # Select points that are super-heated
            I = self.x[4] < 0
            self.h[4][I], self.s[4][I], self.d[4][I] = fluid.hsd(T=self.T[4][I], p=p1)
            I = np.logical_not(I)
            self.h[4][I], self.s[4][I], self.d[4][I] = fluid.hsd(T=self.T[4][I], x=self.x[4][I])
            
        # If state 4 is determined by heat addition
        else:
            self.p[3] = p2
            self.h[3] = self.h[2] + q34
            self.T[3] = fluid.T_h(h=self.h[3], p=p2)
            self.d[3] = fluid.d(T=self.T[3], p=p2)
            self.s[3] = fluid.s(T=self.T[3], d=self.d[3])
            self.x[3] = -1
            
            # Isentropic expansion
            self.p[4] = p1
            self.T[4], self.x[4] = fluid.T_s(s=self.s[3], p=p1, quality=True)
            # Initialize h,s, and d results
            self.h[4] = np.zeros_like(self.T[4])
            self.s[4] = np.zeros_like(self.T[4])
            self.d[4] = np.zeros_like(self.T[4])
            # Select points that are super-heated
            I = self.x[4] < 0
            self.h[4][I], self.s[4][I], self.d[4][I] = fluid.hsd(T=self.T[4][I], p=p1)
            I = np.logical_not(I)
            self.h[4][I], self.s[4][I], self.d[4][I] = fluid.hsd(T=self.T[4][I], x=self.x[4][I])
            
        
        self.w = [
                self.h[0] - self.h[1],
                0.,
                0.,
                self.h[3] - self.h[4],
                0.]
        
        self.q = [
                0.,
                self.h[2] - self.h[1],
                self.h[3] - self.h[2],
                0.,
                self.h[0] - self.h[4]]
        
        self.meta['eta'] = (self.w[3]+self.w[0]) / (self.q[1] + self.q[2])
        self.meta['uE'] = pm.config['unit_energy']
        self.meta['uT'] = pm.config['unit_temperature']
        self.meta['uP'] = pm.config['unit_pressure']
        self.meta['uM'] = pm.config['unit_matter']
        self.meta['uV'] = pm.config['unit_volume']
        self.meta['fluid'] = fluid
        
        
    def tsplot(self, ax=None, fig=None, slabels=True, 
            satstyle=None, procstyle=None, statestyle=None):
        """
"""
        # If the properties are empty, abort
        if not self.T:
            raise PMCycleError('The plot cannot be generated for a cycle that has not yet been updated.  Run update().')
        if isinstance(self.T, np.ndarray) and self.T.size>1:
            raise PMCycleError('Plotting arrays of cycle data is not supported.')
        
        ax = self._initplot(
                's (%s/%s%s)'%(self.meta['uE'], self.meta['uM'], self.meta['uT']), 
                'T (%s)'%(self.meta['uT']), ax=ax, fig=fig)
        
        if satstyle is None:
            satstyle = {'lw':2, 'c':'k', 'ls':'solid', 'marker':'None'}
        if procstyle is None:
            procstyle = {'lw':2, 'c':'r', 'ls':'solid', 'marker':'None'}
        if statestyle is None:
            statestyle = {'ls':'None', 'marker':'o', 'mec':'k', 'mew':1, 'mfc':'r', 'ms':6}
        boxstyle = {'fc':'w', 'ec':'k', 'boxstyle':'square,pad=.25'}
            
        # Get the working fluid object
        fluid = self.meta['fluid']
            
        # Plot the dome
        Tc,pc = fluid.critical()
        Tt,pt = fluid.triple()
        T = np.linspace(Tt, 0.99999*Tc, 100)
        sL, sV = fluid.ss(T=T)
        ax.plot(sL, T, **satstyle)
        ax.plot(sV, T, **satstyle)
            
        # Isentropic pump is a straight vertical line
        ax.plot( 
                [self.s[0], self.s[1]], 
                [self.T[0], self.T[1]],
                **procstyle)
                
        # Isobaric boiler
        s = np.linspace(self.s[1], self.s[2], 51)
        T = fluid.T_s(s, p=self.p[1])
        ax.plot(s, T, **procstyle)
        
        # Isobaric super-heater
        s = np.linspace(self.s[2], self.s[3], 21)
        T = fluid.T_s(s, p=self.p[1])
        ax.plot(s, T, **procstyle)
        
        # Isentropic turbine is a straight vertical line
        ax.plot(
                [self.s[3], self.s[4]],
                [self.T[3], self.T[4]],
                **procstyle)
                
        # Isobaric condensation
        s = np.linspace(self.s[4], self.s[0], 51)
        T = fluid.T_s(s, p=self.p[0])
        ax.plot(s, T, **procstyle)
        
        # deterine offsets in T and s space
        ds = ax.get_xlim()
        ds = .02 * (ds[1] - ds[0])
        dT = ax.get_ylim()
        dT = .02 * (dT[1] - dT[0])
        
        # Deal with the state labels
        ax.plot(self.s[0], self.T[0], **statestyle)
        tt = ax.text(self.s[0]-ds, self.T[0]-dT, '1', ha='right', va='top')
        tt.set_bbox(boxstyle)
        ax.plot(self.s[1], self.T[1], **statestyle)
        tt = ax.text(self.s[1]-ds, self.T[1]+dT, '2', ha='right', va='bottom')
        tt.set_bbox(boxstyle)
        ax.plot(self.s[2], self.T[2], **statestyle)
        tt = ax.text(self.s[2]-ds, self.T[2]+dT, '3', ha='right', va='bottom')
        tt.set_bbox(boxstyle)
        ax.plot(self.s[3], self.T[3], **statestyle)
        tt = ax.text(self.s[3]+ds, self.T[3]+dT, '4', ha='left', va='bottom')
        tt.set_bbox(boxstyle)
        ax.plot(self.s[4], self.T[4], **statestyle)
        tt = ax.text(self.s[4]+ds, self.T[4]-dT, '5', ha='left', va='top')
        tt.set_bbox(boxstyle)
        
        return ax
    
    
class RefrigerationCycle(Cycle):
    pass
    
class OttoCycle(Cycle):
    pass
    
class DieselCycle(Cycle):
    pass
