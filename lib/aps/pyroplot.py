"""PYroPlot: The PYroMat plotting module

ts()        Temperature-entropy diagrams
ph()        Pressure-enthalpy diagrams

"""

import matplotlib.pyplot as plt
import matplotlib.figure
import numpy as np
import pyromat as pm



config = {
    'sat_color':'k',
    'sat_width':2,
    'sat_style':'-',
    'p_color':[1.0,0.5,0.5],
    'p_width':1,
    'p_style':'-',
    'd_color':[0.5,0.5,1.0],
    'd_width':1,
    'd_style':'-',
}


def _interval(start, stop, count, inc=[1,2,5]):
    """Auto-generate conveniently spaced values in a range
    array = _interval(start, stop, count)
    
Generates an array of approximately uniformly spaced values beginning 
after start, and stopping before stop.  Unlike arange() or linspace(),
this function chooses values that are rounded to a nearest convenient 
interval.  This is commonly done for plotting.

The optional 'inc' parameter indicates the intervals that may be used
in any given decade.  They should be values larger than or equal to 1
and smaller than 10.  By default, inc=[1,2,5], which indicates that 
valid step sizes might be .01, .02, .05, .1, .2, .5, 1, 2, 5, 10, 20, 50
etc... 

To produce elegantly roundable numbers, the start value will be rounded
up to the nearest integer multiple of the selected incrementer.
"""
    # Calculate an increment that is the next largerst value in inc
    # scaled to the relevant decade
    dx = abs(stop-start) / count
    power = np.floor(np.log10(dx))
    dx /= 10**power
    inc.sort(reverse=True)
    for ii in inc:
        if ii<=dx:
            break
    dx = ii * 10**power
    
    # Round to the nearest incrementer
    start = dx * np.ceil(start / dx)
    return np.arange(start, stop, dx)
    
def _log_interval(start, stop, count, inc=[1,2,5]):
    """Auto-generate conveniently spaced values in a range
    array = _log_interval(start, stop, count)
    
Generates an array of approximately count values beginning after start, 
and stopping before stop.  Just like _interval(), this function does not
strictly respect the start, stop, or count parameters, but rounds to
the nearest conveniently represented decimal value. 

Unlike _interval(), _log_interval() adjusts the incrementer as the 
series progresses, so that it will always match the decade of the 
current value.

The optional 'inc' parameter indicates the intervals that may be used
in any given decade.  They should be values larger than or equal to 1
and smaller than 10.  By default, inc=[1,2,5], which indicates that 
valid step sizes might be .01, .02, .05, .1, .2, .5, 1, 2, 5, 10, 20, 50
etc... 

To produce elegantly roundable numbers, the start value will be rounded
up to the nearest integer multiple of the selected incrementer.
"""
    inc = np.asarray(inc)
    # transform the problem to a log-space
    log_start = np.log(start)
    log_stop = np.log(stop)
    log_dx = abs(log_stop - log_start) / count
    # produce an evenly spaced logarithmic series
    out = np.exp(np.arange(log_start, log_stop, log_dx))
    dx = out * log_dx
    for ii in range(out.size):
        # Calculate a nominal spacing
        dx = np.abs(out[ii] * log_dx)
        power = np.floor(np.log10(dx))
        dx /= 10**power
        # Adjust it to the nearest inc value
        dx = inc[np.argmin(np.abs(inc-dx))] * 10 ** power
        # Round out to the nearest integer multiple of the dx value
        out[ii] = dx * np.round(out[ii]/dx)
        
    # eliminate any points outside of start and stop
    out = out[np.logical_and(
            out <= max(start,stop),
            out >= min(start,stop)
            )]
    return out


def Tp(mpobj, fig=None, ax=None, Tlim=None, plim=None, dlines=None):
    """Temperature-pressure diagram
"""

    # Select a figure
    if fig is None:
        if ax is not None:
            fig = ax.get_figure()
        else:
            fig = plt.figure()
    elif isinstance(fig, matplotlib.figure.Figure):
        pass
    else:
        fig = plt.figure(fig)
    
    # Select an axes
    if ax is None:
        fig.clf()
        ax = fig.add_subplot(111)
    
    if Tlim is None:
        Tlim = mpobj.Tlim()
        
    if plim is None:
        plim = mpobj.plim()
    
    Tc,pc,dc = mpobj.critical(density=True)
    Tt,pt = mpobj.triple()
    
    plim[0] = max(plim[0], pt)
    
    if dlines is None:
        dlim = [dc/1000., dc]
        DLINES = np.flip(_log_interval(dlim[0], dlim[1], 10), 0)
    else:
        DLINES = np.asarray(dlines)
    
    Tn = (Tc - Tt) / 1000.

    # Generate lines
    T = np.linspace(Tlim[0]+Tn, Tlim[1]-Tn, 151)
    
    # Lines of constant density
    for d in DLINES:
        p = mpobj.p(T=T,d=d)
        ax.loglog(p,T,
                config['d_style'],
                color=config['d_color'],
                lw=config['d_width'])
                
    # Generate the dome
    T = np.linspace(Tt+Tn,Tc-Tn,101)
    p = mpobj.ps(T)

    ax.loglog(p,T,
            ls=config['sat_style'],
            color=config['sat_color'],
            lw=config['sat_width'])
        
    # Label the s-axis
    ax.set_xlabel('p [%s]'%(pm.config['unit_pressure']))
            
    # Label the T-axis
    ax.set_ylabel('T [%s]'%(
            pm.config['unit_temperature']))
        
    # Label the figure
    ax.set_title('%s T-p Diagram'%(mpobj.data['id']))
        
    plt.show(block=False)
    return ax


def Ts(mpobj, fig=None, ax=None, satlines=True, Tlim=None, dlines=None, plines=None):
    """Temperature-enthalpy diagram
    ax = TS(mpobj)
    
"""
    # Select a figure
    if fig is None:
        if ax is not None:
            fig = ax.get_figure()
        else:
            fig = plt.figure()
    elif isinstance(fig, matplotlib.figure.Figure):
        pass
    else:
        fig = plt.figure(fig)
    
    # Select an axes
    if ax is None:
        fig.clf()
        ax = fig.add_subplot(111)
    
    if Tlim is None:
        Tlim = mpobj.Tlim()
    
    Tc,pc,dc = mpobj.critical(density=True)
    Tt,pt = mpobj.triple()
    
    if dlines is None:
        dlim = [dc/1000., 2*dc]
        DLINES = np.flip(_log_interval(dlim[0], dlim[1], 10), 0)
    else:
        DLINES = np.asarray(dlines)
    
    if plines is None:
        plim = mpobj.plim()
        # Force plim to be greater than pt
        plim[0] = max(pt, plim[0])
        PLINES = _log_interval(plim[0], plim[1], 10)
    else:
        PLINES = np.asarray(plines)

    Tn = (Tc - Tt) / 1000.

    # Generate lines
    T = np.linspace(Tlim[0]+Tn, Tlim[1]-Tn, 151)
    
    # Lines of constant pressure
    for p in PLINES:
        s = mpobj.s(T=T,p=p)
        ax.plot(s,T,
                config['p_style'],
                color=config['p_color'], 
                lw=config['p_width'])
    
    # Lines of constant density
    for d in DLINES:
        s = mpobj.s(T=T,d=d)
        ax.plot(s,T,
                config['d_style'],
                color=config['d_color'],
                lw=config['d_width'])
                
    # Generate the dome
    T = np.linspace(Tt+Tn,Tc-Tn,101)
    ssL,ssV = mpobj.ss(T)

    ax.plot(ssL,T,
            ls=config['sat_style'],
            color=config['sat_color'],
            lw=config['sat_width'])
    ax.plot(ssV,T,
            ls=config['sat_style'],
            color=config['sat_color'],
            lw=config['sat_width'])
                
    # LABELS of constant pressure
    dT = 0.8*(Tc-Tt)/len(PLINES)
    for p in PLINES:
        if p<pc:
            T = mpobj.Ts(p)
            ss = mpobj.ss(T)
            s = 0.7*ss[0]+0.3*ss[1]
        else:
            T = Tc + (Tlim[1]-Tc) * (p-pc)/(plim[1]-pc)
            s = mpobj.s(T=T,p=p)
        label = '%s%s'%(str(p),pm.config['unit_pressure'])
        ax.text(s, T, label,
                color=config['p_color'],
                ha='center',
                va='center',
                backgroundcolor='w')
        
    # LABELS of constant density
    T = Tlim[1] - 0.05*(Tlim[1] - Tlim[0])
    dT = 0.8*(Tc - Tt)/len(DLINES)
    for d in DLINES:
        s = mpobj.s(T, d=d)
        label = '%s%s/%s'%(str(d),pm.config['unit_matter'],pm.config['unit_volume'])
        ax.text(s, T, label,
                color=config['d_color'],
                ha='center',
                va='center',
                backgroundcolor='w')
        T-=dT
        
    # Label the s-axis
    ax.set_xlabel('s [%s/(%s%s)]'%(
            pm.config['unit_energy'],
            pm.config['unit_matter'],
            pm.config['unit_temperature']))
            
    # Label the T-axis
    ax.set_ylabel('T [%s]'%(
            pm.config['unit_temperature']))
        
    # Label the figure
    ax.set_title('%s T-S Diagram'%(mpobj.data['id']))
        
    plt.show(block=False)
    return ax
