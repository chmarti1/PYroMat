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
    'p_color':[0.8,0.8,0.8],
    'p_width':1,
    'p_style':'-',
    'd_color':[0.8,0.8,0.8],
    'd_width':1,
    'd_style':'--',
}


def _interval(start, stop, count, inc=[1,2,5]):
    """Auto-generate conveniently spaced values in a range
    array = _interval(start, stop, count)
    
Generates an array of approximately count values beginning after start, 
and stopping before stop.  Unlike arange() or linspace(), this function
chooses values that are rounded to a nearest convenient interval.  This
is commonly done for plotting.

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



def ts(mpobj, fig=None, ax=None, Tlim=None, slim=None):
    """Temperature-enthalpy diagram
    ax = TS(mpobj)
    
"""
    if Tlim is None:
        Tlim = mpobj.Tlim()
        
    if slim is None:
        slim = mpobj.ss(T=Tlim[0])
    
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
    
    Tc,pc,dc = mpobj.critical(density=True)
    Tt,pt = mpobj.triple()
    
    # Generate the dome
    T = np.linspace(Tt,Tc,101)
    ssL,ssV = mpobj.ss(T)

    ax.plot(ssL,T,
            ls=config['sat_style'],
            color=config['sat_color'],
            lw=config['sat_width'])
    ax.plot(ssV,T,
            ls=config['sat_style'],
            color=config['sat_color'],
            lw=config['sat_width'])
    
    PLINES = _log_interval(pt, pc, 10)
    DLINES = _log_interval(dc/1000., dc, 10)
    
    # Lines of constant pressure
    for p in PLINES:
        s = mpobj.s(T=T,p=p)
        ax.plot(s,T,
                config['p_style'],
                color=config['p_color'], 
                lw=config['p_width'])
                
    for d in DLINES:
        s = mpobj.s(T=T,d=d)
        ax.plot(s,T,
                config['d_style'],
                color=config['d_color'],
                lw=config['d_width'])
                
    # LABELS of constant pressure
    for p in PLINES:
        Ts = mpobj.Ts(p)
        ss = mpobj.ss(Ts)
        label = 'p=%s%s'%(str(p),pm.config['unit_pressure'])
        ax.text(0.5*(ss[0]+ss[1]), Ts, label,
                color=config['p_color'],
                ha='center',
                va='center',
                backgroundcolor='w')
        
    ax.set_xlim(mpobj.ss(Tt))
        
    plt.show(block=False)
    return ax
