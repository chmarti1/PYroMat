"""IGTOOLS - A module for ideal gas tools

::Tools for ideal gas mixtures::



"""


import pyromat as pm
import numpy as np


def parse_mixstr(mixstr):
    """PARSE_MIXSTR - construct a mixture dictionary from a mixture string
    mdict = parse_mixstr( ... )
"""
    result = {}
    parts = mixstr.split('+')
    for this in parts:
        coef = 1.0
        subst = None
        # Scan for the first alphabetical character
        for ii,c in enumerate(this):
            if c.isalpha():
                break
        coefstr = this[:ii].strip()
        if coefstr:
            try:
                coef = float(coefstr)
            except:
                raise pm.utility.PMParamError('PARSE_MIXSTR: Failed to parse mixture fragment: ' + this)
        subst = this[ii:].strip()
        
        result[subst] = coef
    return result
            
class IGTMix(object):
    """Ideal Gas Tools dynamic Mixture
    
    The dynamic mixture object is a wrapper for the various ideal gas
classes.  It deals in extensive quantities, and quantities can be 
expressed in mass or molar units.  Once defined, gas mixtures objects 
can then be mixed with each other at the command line.

** DEFINING A MIXTURE **

There are six ways to define a mixture:

(1) From a string...
Strings are expected in the format: 'qty0 subst0 + qty1 subst1 + ...'
All whitespace is ignored.  Omitted quantities are presumed to be unity,
and substances can be specified with or without their 'ig.' collection
prefix.

    mymix = IGTMix('10 ig.N2')
    mymix = IGTMix('2.4 N2 + Ar')

*note* This method cannot be used to specify an ion like 'Ne+'. Instead,
use the dictionary or list methods below.

(2) From keyword arguments...
Keyword arguments accept abbreviated substance ID strings as keywords 
(see below) under "SUBSTANCE ID".
    
    air = IGTMix(N2=0.76, O2=0.23, Ar=0.01)

*note* This method cannot be used to specify an ion like 'Ne+'. Instead,
use the dictionary or list methods below.

(3) A list of constituents with no quantities...
All of the constituents will be initialized containing one unit_matter
in the mixture.
 
    mymix = IGTMix(['ig.N2', 'ig.O2', 'ig.Ar'])
    
(4) A dictionary of constituents with their quantities as values
By default, quantities will be understood in the config['unit_matter']
units.  See UNITS below for more information.

    air = IGTMix({'ig.N2':0.76, 'ig.O2':0.23, 'ig.Ar':.01})

(5) Algebraicaly...
    air = IGTMix('.76 N2 + .23 O2 + .01 Ar')
    fuel = IGTMix('.23 CH4 + .44 C3H8')
    reactants = air + 0.4*fuel

(6) From an existing igmix instance, using the fromigmix() function.
    air = fromigmix(pm.get('ig.air'))
    
When the fromigmix() function is used, the igmix instance is split into
its constituent gases to form the IGTMix instance, so the example above
results in a mixture with argon, carbon dioxide, nitrogen, and oxygen.  
If an igmix instance is passed directly to the IGTMix class, it is 
treated the same as any other constituent gas, so the example below 
produces a mixture that only has one gas component.
    air = IGTMix('ig.air')

** SPECIFYING SUBSTANCES **

In the examples above, constituent gases are specified by their 
substance ID string, but they may be specified three ways:
(1) By their full substance ID string
    mymix = IGTMix('ig.N2')
    
(2) By their abbreviated substance ID string (with no leading ig.)
    mymix = IGTMix('N2')
    
(3) Or by their full data instance
    n2 = pm.get('ig.N2')
    mymix = IGTMix(n2)
    
Any of these -- 'N2', 'ig.N2', or pm.get('ig.N2') -- may be used as the
substance identifier in the mixture list or dictionary methods above.
    
** QUANTITIES AND UNTIS **

The IGTMix class interprets quantities in the matter units configured in
config['unit_matter'].  This can be overridden by the 'units' keyword, 
but the mixture will be immediately converted into configured units.
For example, this code segment defines a mixture with 5.2kmol of water 
and 4.7kmol of carbon dioxide, but when queried, the answer is returned
in kilograms.
    pm.config['unit_matter'] = 'kg'
    mymix = IGTMix({'H2O':5.2, 'CO2':4.7}, units='kmol')
    mymix['H2O']
        array([93.679456])

In cases where no quantity is specifed, IGTMix always assumes one unit 
of matter.  This lets the algebraic methods above (in method 5) work.

If the units are changed after an IGTMix instance has already been 
defined, they will automatically be converted into the new units the 
next time one of the instance methods is called.  To ensure this process
is correctly handled, users should never interact with an IGTMix instance
through members or methods with a leading underscore -- these are 
intended for internal use only.

Quantities are always handled as arrays, so a single IGTMix instance can
actually manage arrays of mixtures with an arbitrary shape.  In this 
example the mixture instance is an array of mixtures of neon and its 
first ion.
    x = np.linspace(0,0.1,21)
    mymix = IGTMix{{'Ne':1-x, 'Ne+':x})
    mymix.shape()
        (21,)

** INDEXING **

Indexing a mixture allows users to access and change the mixture 
quantities.  The first index can be an integer or a string with the 
substance ID.  If there are additional indices, they indicate the 
element of the array to access.
    mymix = IGTMix({'C3H8':[2,1], 'CH4':[5,10]})
    mymix['CH4']
        array([ 5., 10.])
    mymix['C3H8',1]
        1.0
    mymix['CH4', 0] = 4.

** ITERATING **

The items() method works like the Python built-in dict.items().  It 
returns an iterator that produces an ordered pair of each substance 
instance and its quanitty array.
    for subst,qty in mymix.items():
        ... do some stuff ...

** ALGEBRA WITH MIXTURES **

Because mixtures work in absolute quantities (as opposed to mass or 
mole fractions), they can be incrased, decreased, subtracted from, or
added to using basic math operations at the command line.  In this 
example, mix3 has 1kmol of water, 2kmol carbon dioxide, and 2kmol each
of nitrogen and argon.
    pm.config['unit_matter'] = 'kmol'
    mix1 = IGTMix({'H2O':2, 'CO2':4})
    mix2 = IGTMix({'N2':1, 'Ar':1})
    mix3 = 0.5*mix1 + 2*mix2

** REPRESENTATION **

The mixture is represented in the back-end by attributes that are not 
intended to be accessed directly.  Instead, these attributes should only
be read or modified by the appropriate methods.

_q      is a numpy array of quantities; the first index in the array
        identifies the substance.  All other dimensions may be freely
        broadcast to eachother.
_c      is a numpy array of IGMix or ideal gas instances that make up
        the mixture.
_units  is a string identifying the matter units that were in use when
        the mixture was defined.
"""

    def __init__(self, contents=None, units=None, **kwarg):
        
        # Before we do anything, let's deal with the cases that call for
        # recursion.  All cases should route to a call to __init__ with
        # a dictionary assigned to contents.
        # If no contents have been specified...
        if contents is None:
            return self.__init__(contents=kwarg, units=units)
        # If this is a copy operation
        elif isinstance(contents, IGTMix):
            self._c = contents._c.copy()
            self._q = contents._q.copy()
            self._units = contents._units
            self._update()
            return
        # If the contents is a string or data instance, build a dummy dict
        elif isinstance(contents, str):
            return self.__init__(contents=parse_mixstr(contents), units=units)
        # If contents is a base PYroMat class
        elif isinstance(contents, pm.reg.__basedata__):
            return self.__init__(contents={contents:1}, units=units)
        # if the contents is a list or tuple of constituents
        elif isinstance(contents, (list, tuple)):
            return self.__init__(dict.fromkeys(contents,1), units=units)
        # Finally, if the contents is not a dict, raise an error
        elif not isinstance(contents, dict):
            raise Exception(pm.utility.PMParamError('Unexpected data type for IGTMix contents: ' + repr(type(contents))))
        
        
        # At this point, contents is always a dict with substance 
        # identifiers and quantities. We'll loop over the contents twice.
        # In the first iteration, we'll determine the quatity array 
        # shape and we'll gather the substance instances.
        shape = None
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
                raise pm.utility.PMParamError('IGTMix instances cannot be nested. Use the + operator instead to combine mixtures.')
            else:
                raise pm.utility.PMParamError('IGTMix: Unrecognized constituent: ' + repr(sid))

            # Next, process the quantity as an array
            qty = np.atleast_1d(qty)
            try:
                shape = np.broadcast(qty, shape)
            except:
                raise pm.utility.PMParamError('IGTMix: Failed to broadcast quantity array for: ' + repr(sid))
            
        # Initialize the quantity array
        N = len(self._c)
        self._q = np.empty((N,) + shape.shape, dtype=float)

        # In the second loop, we'll broadcast all of the quantity arrays
        # to the correct shape.
        for index,qty in enumerate(contents.values()):
            self._q[index,:] = np.broadcast_to(qty, shape.shape)
            
        # Finally, stash the units string
        if units is None:
            self._units = pm.config['unit_matter']
        elif units in pm.units.mass or units in pm.units.molar:
            self._units = units
        else:
            raise pm.utility.PMParamError('IGTMix: Unrecognized unit matter: ' + repr(units))
        
        self._update()    
        
        
    def __iadd__(self, b):
        """Add mixture b to self
    a += b
    
This combines the quantities of gas in mixture b into mixture a.
"""
        # If this is not an IGTMix, see if we can convert it to a mixture
        if not isinstance(b, IGTMix):
            return self.__iadd__(IGTMix(b))
        # Detect the number of unique substances and add them 
        # to the result list
        N = len(self._c)
        NN = N
        for bc in b._c:
            if bc not in self._c:
                self._c.append(bc)
                NN += 1
        # Grow the q array appropariately
        # First, detect the starting and finishing array shapes
        sshape = self.shape()
        bshape = b.shape()
        bb = np.broadcast(self._q[0,:], b._q[0,:])
        
        # If the existing in-place array needs to grow, do that first
        if bb.shape != sshape:
            self._q = np.broadcast_to(self._q, (N,)+bb.shape)
        # If there are new substances, grow the quantity array
        if NN>N:
            self._q = np.concatenate((self._q,np.zeros((NN-N,)+bb.shape)),axis=0)
        
        for bi,subst in enumerate(b._c):
            ai = self._c.index(subst)
            self._q[ai,:] += pm.units.matter(b._q[bi,:],subst.mw(),b._units,self._units)
        
        self._update()
        
        return self
        
    
    def __add__(self, b):
        """Combine two mixtures to form a third unique mixture
    c = a + b
    
Makes a copy of a and calls __iadd__() to execute the operation.
"""
        # Use in-place addition
        # Make a copy of self and return the in-place addition result
        c = IGTMix(self)
        c.__iadd__(b)
        return c
        
    def __radd__(self, b):
        c = IGTMix(self)
        c.__iadd__(b)
        return c
      
        
    def __imul__(self, b):
        """Scale the quantities of a mixture
    a *= qty
"""
        # Force b into a shape that can be broadcast over the entire
        # _q array.
        shape = self.shape()
        b = np.atleast_1d(b)
        # If the shapes are already matched, no broadcasting will happen
        # In-place multiplication can happen
        if b.shape == shape:
            self._q *= b
        # Otherwise, there will be a new _q array
        else:
            self._q = b * self._q
        self._update()
        return self
        
    def __mul__(self, b):
        """Scale the quantity of a mixture
    c = a * qty
"""
        c = IGTMix(self)
        c *= b
        return c
        
    def __rmul__(self,b):
        return self.__mul__(b)
        

    def _update(self):
        # Verify that _units agrees with the current matter setting
        # If not, convert!
        if self._units != pm.config['unit_matter']:
            for subst,qty in self.items():
                pm.units.matter(qty, subst.mw(), from_units=self._units, inplace=True)        
            self._units = pm.config['unit_matter']
        
    def _peval(self, prop, *varg, **kwarg):
        """PEVAL - evaluate a property for all constituents
    parray = m._peval(prop, ...)

The property to be evaluated is evaluated by a string, prop. For example 
a call to evaluate enthalpy might appear:
    parray = m_peval('h', T=300, p=1.01325)
    
All arguments, with or without keywords, are passed verbatim to the 
property method of each constituent substance.
"""
        out = None
        for ii,subst in enumerate(self._c):
            v = getattr(subst,prop)(*varg, **kwarg)
            # If this is the first result, we'll rely on the substance's
            # broadcasting to determine the dimension of the result
            if out is None:
                shape = np.atleast_1d(v).shape
                out = np.zeros((len(self._c),) + shape, dtype=float)
            out[ii,:] = v
        return out
        
    def _sindex(self, index):
        """SINDEX - convert an item index into an array index
    sindex = m._sindex(index)
    
When the index is an integer, string, slice, or a data instance it is 
interpreted as indexing the constituent substances.  Strings are 
converted into substance ID strings, and data instances are matched 
against their 
"""
        # If we're looking for an SID string
        if isinstance(index, tuple):
            # Force the first index to be an integer
            return (self._sindex(index[0]),) + index[1:]
        # Integers and slices pass thru
        elif isinstance(index, (int,slice)):
            return index
        # Strings are searched for
        elif isinstance(index, str):
            if '.' not in index:
                index = 'ig.' + index
            for ii,this in enumerate(self._c):
                if isinstance(this,pm.reg.__basedata__) and this.sid()==index:
                    return ii
            raise pm.utility.PMParamError('Not in this mixture: ' + index)
        # If we're looking for anything else...
        else:
            try:
                return self._c.index(index)
            except ValueError:
                raise pm.utility.PMParamError('Not in this mixture: ' + repr(index))
            except:
                raise

    #
    # The property interface
    #
    def _argparse(self, T=None, p=None, d=None, v=None):
        """ARGPARSE - provides the standard property argument interface
    _argparse(self, T=None, p=None, d=None, v=None)
The IGTMix class is currently implemented with a much simpler argument
interface than the core PYroMat classes.  It expects
"""

    #
    # Indexing
    #
    def __getitem__(self, index):
        return self._q[self._sindex(index)]
        
    def __setitem__(self, index, value):
        self._q[self._sindex(index)] = value
        
    #
    # Length and size operations
    #
    def __len__(self):
        """LEN - detect the number of mixtures represented here
"""
        return np.product(self.shape())
    
    def __str__(self):
        """STR - pretty print of the mixture
"""
        out = ''
        if self.__len__() > 1:
            for sid in self._c[:-1]:
                out += f'[...] {sid} + '
            out += f'[...] {self._c[-1]}'
        else:
            for sid,qty in zip(self._c[:-1], self._q[:-1]):
                out += f'{qty} {sid} + '
            out += f'{self._q[-1]} {self._c[-1]}'
        return out
    
    def shape(self):
        """SHAPE - return the dimensions of the quantity arrays
"""
        return self._q.shape[1:]


    def nsubst(self):
        """NSUBST - return the number of constituent substances
"""
        return len(self._c)
        
    #
    # Type manipulation
    #
    def toigmix(self,sid):
        """TOIGMIX - convert the mixture to a static IGMIX instance
    [... igmix list ...] = igtm.toigmix(sid)
    
Returns a list of static igmix (low-level PYroMat data instances) built
from the IGTMix instance quantity arrays. These are less flexible, but
substantially faster than the IGTMix instances for property calculation.

sid     substance id string to use as the igmix name

If the IGTMix instance only defines one mxiture (shape == (1,)), the
sid string is used verbatim to name the new igmix instance.  If there is
more than one mixture, each is appended with "_{index}" in the same 
manner as the PYroMat collections' deduplication scheme.

It is good practice to choose an SID string that is unique (not in the
PYroMat collection), but it is not strictly required unless this is to
be saved in the permanent collection for later access.

The simpler lower-level igmix class does not allow dynamic manipulation
of mixtures, but it is substantially faster for a number of reasons:
    1) igmix instances store important intermediate calculations
    2) Unlike IGTMix instances, igmix instances call the ig and ig2 
       inner routines directly, so the top-level unit conversions and
       argument parsing is not performed redundantly.
    3) igmix instances have a faster and more flexible property 
       interface.
    4) igmix instances can be saved and added to the permanent 
       collection, while IGTMix instances cannot.

"""
        out = []
        multiple_f = len(self)>1
        for count, index in enumerate(np.ndindex(self.shape())):
            data = {'class':'igmix', 
                'doc':'Auto-generated using IGTMix.toigmix()',
                'fromfile':__name__}
            if multiple_f:
                data['id'] = sid + '_' + repr(count)
            else:
                data['id'] = sid
                
            slice_index = (slice(0,len(self._c)),) + index
            contents = {subst.sid():qty for subst,qty in zip(self._c, self._q[slice_index])}
            data['contents'] = contents
            data['bymass' ] = pm.units.ismass(self._units)
            out.append(pm.reg.registry['igmix'](data))
        return out
        
    #
    # Measures of quantity
    #
    
    def atoms(self):
        """ATOMS - count the number of each atom in the mixture
    atoms = m.atoms()
    
returns a dictionary with the atomic symbol as the key and the quantity
array as the value.  The values are always in configured 'unit_molar'.

Note that many users may be in the habit of setting the 'unit_matter' 
setting without updating 'unit_molar' or 'unit_mass' accordingly.  
This will yield confusing results!
"""
        out = {}
        for subst,qty in self.items():
            qty = pm.units.matter(qty, subst.mw(), from_units=self._units, to_units=pm.config['unit_molar'])
            for atom,aq in subst.atoms().items():
                working = out.get(atom)
                if working is None:
                    working = np.zeros(qty.shape, dtype=float)
                    out[atom] = working
                working += qty * aq
        return out
            
    def X(self):
        """X - calculate the mole fractions of the substances
    x = m.X()
    
returns a dictionary with the substance id strings as the keys and the
fraction array as values.  Values are always between zero and one.
"""
        x = {}
        # If this is in mass-based units, we'll need to convert to molar
        if pm.units.ismass(self._units):
            # First, total all the moles
            total = np.zeros(self.shape())
            for subst,qty in self.items():
                temp = qty / subst.mw()
                total += temp
                x[subst.sid()] = temp
            # Finally, normalize by the total
            for key,qty in x.items():
                x[key] /= total
        else:
            total = np.sum(self._q, axis=0)
            for subst,qty in self.items():
                x[subst.sid()] = qty/total
        return x
        
    def Y(self):
        """Y - calculate the mass fractions of the substances
    y = m.Y()
    
returns a dictionary with the substance id strings as the keys and the
fraction array as values.  Values are always between zero and one.
"""
        y = {}
        # If this is in molar-based units, we'll need to convert to mass
        if pm.units.ismass(self._units):
            total = np.sum(self._q, axis=0)
            for subst,qty in self.items():
                y[subst.sid()] = qty/total
        else:
            # First, total all the mass
            total = np.zeros(self.shape())
            for subst,qty in self.items():
                temp = qty * subst.mw()
                total += temp
                y[subst.sid()] = temp
            # Finally, normalize by the total
            for key,qty in y.items():
                y[key] /= total
        return y
        
        
    def mass(self):
        """MASS - calculate the total mass of the mixture
    mass = m.mass()
    
Returns an array with m.shape() shape.  Always reports units in the 
configured 'unit_mass' units.  Note that many users may be in the habit
of setting the 'unit_matter' setting without updating 'unit_molar' or 
'unit_mass' accordingly.  This will yield confusing results!
"""
        if pm.units.ismass(self._units):
            temp = pm.units.mass(self._q, from_units=self._units, to_units=pm.config['unit_mass'])
            return np.sum(temp, axis=0)
        else:
            total = np.zeros(self.shape())
            for subst,qty in self.items():
                total += pm.units.matter(qty, subst.mw(), from_units=self._units, to_units=pm.config['unit_mass'])
            return total
        raise pm.utility.PMAnalysisError('PANIC! Unhandled exception.')
        
        
    def molar(self):
        """MOLAR - calculate the total mole count of the mixture
    molar = m.molar()
    
Returns an array with m.shape() shape.  Always reports units in the 
configured 'unit_molar' units.  Note that many users may be in the habit
of setting the 'unit_matter' setting without updating 'unit_molar' or 
'unit_mass' accordingly.  This will yield confusing results!
"""
        if pm.units.ismass(self._units):
            total = np.zeros(self.shape())
            for subst,qty in self.items():
                total += pm.units.matter(qty, subst.mw(), from_units=self._units, to_units=pm.config['unit_molar'])
            return total
        else:
            temp = pm.units.molar(self._q, from_units=self._units, to_units=pm.config['unit_molar'])
            return np.sum(temp, axis=0)
        raise pm.utility.PMAnalysisError('PANIC! Unhandled exception.')
        
    #
    # Iteration methods
    #
    
    def items(self):
        """ITEMS - return an iterator on the contents of the mixture
        
This works like the items() method in Python dictionaries
    for subst,qty in mix.items():
        ...
        
In this loop format, subst will be a PYroMat data instance from the ideal
gas collection, or (if there is a recursively defined mixture) it will
be another IGTMix instance.  qty is the quantity array corresponding to
each substance.
"""
        return zip(self._c, self._q)

    #
    # Properties
    #
    
    def mw(self):
        """MW - molecular weight

    mw = m.mw()
    
Returns the molecular weight in [unit_mass] / [unit_molar].
"""
        self._update()
        mw = self._peval('mw')
        total = np.sum(self._q, axis=0)
        out = 0.
        if pm.units.ismass(self._units):
            for ii,subst in enumerate(self._c):
                out += self._q[ii] / mw[ii]
            out = total / out
        else:
            for ii,subst in enumerate(self._c):
                out += self._q[ii] * mw[ii]
            out /= total
        return out
        
    def R(self):
        """R - ideal gas constant
        
Returns the gas constant in [unit_energy]/[unit_matter][unit_temperature]
"""
        # No need to call _update() - it will be called by mw()
        R = np.broadcast_to(pm.units.const_Ru, self.shape())
        R = pm.units.energy(R, from_units='J')
        R = pm.units.temperature(R, from_units='K', exponent=-1)
        R = pm.units.matter(R, self.mw(), from_units='mol', exponent=-1)
        return R

    
        



def fromigmix(source):
    """FROMIGMIX - create an IGTMix instance from and igmix instance
    igtmix = fromigmix(igm)
    
Creates an IGTMix instance from a low-level igmix instance.  The IGTMix
instance will automatically contain a single unit of matter (with the 
currently configured unit matter).  
"""
    if pm.units.ismass():
        return IGTMix(source.Y())
    else:
        return IGTMix(source.X())
