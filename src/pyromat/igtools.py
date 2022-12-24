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
            qty = np.atleast_1d(qty)
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
            
        # Finally, stash the units string
        if units is None:
            self._units = pm.config['unit_matter']
        elif units in pm.units.mass or units in pm.units.molar:
            self._units = units
        else:
            raise pm.utility.PMParamError('IGTMix: Unrecognized unit matter: ' + repr(units))
            
    def __iadd__(self, b):
        """Add mixture b to mixture a
    a += b
"""
        # Detect the number of unique substances and add them 
        # to the result list
        N = len(self._c)
        NN = N
        for bc in b._c:
            if bc not in self._c:
                self._c.append(bc)
                NN += 1
                
        # Grow the q array appropariately
        # Check for broadcastable quantity shapes
        ashape = self.shape()
        cshape = np.broadcast_shapes(ashape,b.shape())
        # If broadcasting is necessary, do it now
        if cshape != ashape:
            self._q = np.broadcast_to(self._q, (N,)+cshape)
        # If there are new substances, grow the quantity array
        if NN>N:
            self._q = np.concatenate((self._q,np.zeros((NN-N,)+cshape)),axis=0)
        
        for bi,subst in enumerate(b._c):
            ai = self._c.index(subst)
            self._q[ai,:] += pm.units.matter(subst.mw(),b._q[bi,:],b._units,self._units)
        return self
                
                
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
        
    def _sindex(self, sid):
        """SINDEX - retrieve the index corresponding to a substance
    index = _sindex(sid)
    
SID can be a substance ID string or a data instance to be matched.

Returns None if not found.
"""
        # If we're looking for an SID string
        if isinstance(sid, str):
            if '.' not in sid:
                sid = 'ig.' + sid
            for index,this in enumerate(self._c):
                if isinstance(this,pm.reg.__basedata__) and this.sid()==sid:
                    return index
            return None
        # If we're looking for anything else...
        else:
            for index,this in enumerate(self._c):
                if sid is this:
                    return index
            return None
        
    #
    # Length and size operations
    #
    def __len__(self):
        """LEN - detect the number of mixtures represented here
"""
        return np.product(self.shape())
            
    def shape(self):
        """SHAPE - return the dimensions of the quantity arrays
"""
        return self._q.shape[1:]


    def nsubst(self):
        """NSUBST - return the number of constituent substances
"""
        return len(self._c)
        
    def sety(self, y, sid=None):
        pass
        
    def setx(self, x, sid=None):
        pass
        



def get(identifier):
    """Return an 
    
"""
    pass
