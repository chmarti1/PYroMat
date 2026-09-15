import pyromat as pm
import numpy as np

def isovalues(subst, xprop, xlim, yprop, ylim, isoprop, xlog=False, ylog=False, isolog=None):
    """Recommend isoline values on an x,y property plot
    value_array = isolvalues(subst, xprop, xlim, yprop, ylim, isoprop)
    
Recommended values will appear on a plot as evenly spaced as is 
practical with round numbers, minimizing the number of significant 
figures that need to be displayed on the plot.
    
Arguments:
subst   
    The PYroMat substance instance
    
xprop, yprop, isoprop
    Strings (one character eacg) indicating which property is being
    plotted on the x- and y-axes and which property's isoline is being
    generated.  Accepted characters are
        'e'     Internal energy
        'd'     Density
        'f'     Helmholtz free energy
        'g'     Gibbs energy
        'h'     Enthalpy
        'p'     Pressure
        's'     Entropy
        'T'     Temperature
        'v'     Specific volume
        
xlim, ylim
    Two-element list or tuple representing the [lower, upper] limits of
    the plot.
    
Optional arguments are
xlog, ylog  (default False)
    Indicate whether the plot is using logarithmic x- or y-axes.
    
isolog      (default None)
    Indicate whether the isovalues should progress linearly or
    logarythmically.  If set to None, isovalues() decides automatically.
"""
    pass


