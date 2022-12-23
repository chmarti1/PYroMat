import pyromat as pm
import numpy as np



class IGTMix(object):
    """Ideal Gas Tools dynamic Mixture
    
    The dynamic mixture object is a wrapper for the various ideal gas
classes.  It deals in extensive quantities, and quantities can be 
expressed in mass or molar units.  Once defined, gas mixtures objects 
can then be mixed with each other at the command line.

** DEFINING A MIXTURE **

There are five ways to define a mixture:

(1) A single-constituent mixture...
This is most useful as a means to build mixtures algebraicaly (see #2) 
from the command line.

    m = IGTMix('ig.N2')

(2) A list of constituents with no quantities...
All of the constituents will be initialized containing one unit_matter
in the mixture.
 
    m = IGTMix(['ig.N2', 'ig.O2', 'ig.Ar'])
    
(3) A dictionary of constituents with their quantities as values
By default, quantities will be understood in the config['unit_matter']
units.  See UNITS below for more information.

    air = IGTMix({'ig.N2':0.76, 'ig.O2':0.23, 'ig.Ar':.01})

(4) From keyword arguments...
Keyword arguments accept abbreviated substance ID strings as keywords 
(see below) under "SUBSTANCE ID".
    
    air = IGTMix(N2=0.76, O2=0.23, Ar=0.01)

(5) Algebraicaly from other mixtures...
    n2 = IGTMix('ig.N2')
    o2 = IGTMix('ig.O2')
    ar = IGTMix('ig.Ar')
    air = 0.76*n2 + 0.23*o2 + 0.01*ar

** SUBSTANCE ID **

In the examples above, constituent gases are specified by their 
substance ID string, but they may be specified three ways:
(1) By their full substance ID string
    m = IGTMix('ig.N2')
(2) By their abbreviated substance ID string (with no leading ig.)
    m = IGTMix('N2')
(3) Or by their full data instance
    n2 = pm.get('ig.N2')
    m = IGTMix(n2)
    
** UNTIS **

By default, the IGTM class will interpret quantities in the units in
config['unit_matter'], but units may alternately be specified by the
optional "units" keyword.

    m = IGTMix({'ig.N2':0.76, 'ig.O2':0.23, 'ig.Ar':.01}, units = 'kg')
    
** ALGEBRA WITH MIXTURES **

Because mixtures work in absolute quantities (as opposed to mass or 
mole fractions), they can be incrased, decreased, subtracted from, or
added to using basic math operations at the command line.

** REPRESENTATION **

The mixture is represented in the back-end by attributes that are not 
intended to be accessed directly.  Instead, these attributes should only
be read or modified by the appropriate methods.

_q      is a numpy array of quantities; the first index in the array
        identifies the substance.  All other dimensions may be freely
        broadcast to eachother.
_c      is a numpy array of IGMix or ideal gas instances that make up
        the mixture.
_units  is a string identifying the matter units used in the _q array.
"""

    def __init__(self, contents=None, units=None, **kwarg):
        
        # Before we do anything, let's deal with the cases that call for
        # recursion.  All cases should route to a call to __init__ with
        # a dictionary assigned to contents.
        # If no contents have been specified...
        if contents is None:
            # Check the kwarg dict
            if len(kwarg)==0:
                raise pm.utility.PMParamError('Cannot initialize an IGTMix instance without contents.')
            return self.__init__(contents=kwarg, units=units)
        # If the contents is a string or data instance, build a dummy dict
        elif isinstance(contents, (str, IGTMix, pm.reg.__basedata__)):
            return self.__init__(contents={contents:1}, units=units)
        # if the contents is a list or tuple of constituents
        elif isinstance(contents, (list, tuple)):
            return self.__init__(dict.fromkeys(contents,1), units=units)
        # Finally, if the contents is not a dict, raise an error
        elif not isinstance(contents, dict):
            raise Exception(pm.utility.PMParamError('Unexpected data type for IGTMix contents: ' + repr(type(contents))))
        
        
        # At this point, contents is always a dict with substance 
        # identifiers and quantities. We'll loop over the contents twice.
        # In the first iteration, we'll broadcast the quatity array 
        # shapes, and we'll gather the substance instances.
        shape = ()
        self._c = []
        for sid,qty in contents.items():
            # First, deal with the substance
            # If this is an ID string
            if isinstance(sid,str):
                if '.' not in sid:
                    sid = 'ig.' + sid
                # Go find the substance
                # This will raise a meaningful error if it is not found
                self._c.append(pm.get(sid))
            # If this is a PM data class
            elif isinstance(sid,pm.reg.__basedata__):
                # Verify that it's an ideal gas first
                if sid.collection() != 'ig':
                    raise pm.utility.PMParamError('IGTMix constituent was not an ideal gas: ' + sid.sid())
                self._c.append(sid)
            # IGTMix components are included automatically
            elif isinstance(sid,IGTMix):
                self._c.append(sid)
            else:
                raise pm.utility.PMParamError('IGTMix: Unrecognized constituent: ' + repr(sid))

            # Next, process the quantity as an array
            qty = np.atleast_1d(np.asarray(qty))
            try:
                shape = np.broadcast_shapes(qty.shape, shape)
            except:
                raise pm.utility.PMParamError('IGTMix: Failed to broadcast quantity array for: ' + repr(sid))
            
        # Initialize the quantity array
        N = len(self._c)
        self._q = np.empty((N,) + shape, dtype=float)

        # In the second loop, we'll broadcast all of the quantity arrays
        # to the correct shape.
        for index,qty in enumerate(contents.values()):
            self._q[index,:] = np.broadcast_to(qty, shape)

            
    def __iadd__(self, b):
        """Add mixture b to mixture a
    a += b
"""
        for sid,(subst, qty) in self.items():
            if sid in self._contents:
                # If the substances don't match raise an error
                if self._contents[sid][0] != subst:
                    pm.print_error('While adding two mixtures, substances with the same ID strings were found pointing to contradictory data instances. This means there are two contradictory definitions for the same substance ID.')
                    raise pm.utility.PMDataError('Contradictory data instances for SID: ' + str(sid))
                self._contents[sid][1] += qty
            else:
                self._contents[sid] = [subst,qty]
                
    def __add__(self, b):
        """Combine two mixtures to form a third unique mixture
    c = a + b
"""
        c = IGTMix(self)
        c += b
        return c
        
    def _update(self):
        """_UPDATE - a back-end helper function that updates all parameters

_update() should be called any time there is a change to the mixture 
composition.  It re-calculates the molecular weight, mole, and mass 
fractions.  The IGMIX class has to use the _bootstrap() method for this
"""
        
    def sety(self, y, sid=None):
        pass
        
    def setx(self, x, sid=None):
        pass
        



def get(identifier):
    """Return an 
    
"""
    pass
