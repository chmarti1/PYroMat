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

        OR JUST AS A STRING!
        
    reactants = '.76 N2 + .23 O2 + .01 Ar' + 0.4*fuel

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

IGTMix instances act much like an ordered dictionary with PYroMat ideal
gas class isntances as keys and quantitiy arrays as values.  The items() 
method allows simultaneous iteration over the substances and quantities.

    for subst,qty in mymix.items():
        # subst is now a PYroMat substance instance
        # qty is now a Numpy array of the substance's quantity
        
    for subst in mymix:
        # subst is now a PYroMat substance instance

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

The addition algorithm attempts to convert non-IGTMix instances to 
IGTMix instances.  That means that strings, lists, dictionaries, and
ordinary PYroMat instances may be folded into mixtures using plain 
command-line algebra.
    mix4 = mix2 + '0.5 CO'
    mix5 = mix2 + pm.get('ig.H2O')
    mix6 = mix2 + {'H2O':[0.5,0.12], 'C2H2':0.8}

** REPRESENTATION **

The mixture is represented in the back-end by attributes that are not 
intended to be accessed directly.  Instead, these attributes should only
be read or modified by the appropriate methods.

_q      is a numpy array of quantities; the first index in the array
        identifies the substance.  All other dimensions may be freely
        broadcast to eachother.
_c      is a numpy array of IGMix or ideal gas instances that make up
        the mixture.
_units  is a string identifying the matter units in which the quanitity
        array is currently expressed.
"""

    def __init__(self, contents=None, units=None, **kwarg):
        
        # Before we do anything, let's deal with the cases that call for
        # recursion.  All cases should route to a call to __init__ with
        # a dictionary assigned to contents.
        # If no contents have been specified...
        if contents is None:
            if kwarg:
                return self.__init__(contents=kwarg, units=units)
            else:
                contents = {}
        # If this is a copy operation
        elif isinstance(contents, IGTMix):
            self._c = list(contents._c)
            self._q = np.array(contents._q)
            self._units = str(contents._units)
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
        shape = (1,)
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
            qshape = np.shape(qty)
            shape = self._broadcast_shape(qshape, shape)
            
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
        
        self.update()
        
    def __iter__(self):
        return self._c.__iter__()
        
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
        bshape = b.shape()
        shape = self._broadcast_shape(b.shape())
        
        # If the q-array does not change shape, then this can be done
        # in-place.  Otherwise, we're going to need to declare a new
        # q-array.
        if NN!=N or self.shape() != shape:
            q = np.zeros((NN,) + shape)
            q[:N] = self._q
            self._q = q
                
        # If units are not equal.
        if self._units != b._units:
            for bi,subst in enumerate(b._c):
                ai = self._c.index(subst)
                self._q[ai,:] += + pm.units.matter(b._q[bi,:],subst.mw(),b._units,self._units)
        else:
            for bi,subst in enumerate(b._c):
                ai = self._c.index(subst)
                self._q[ai,:] += b._q[bi,:]
        
        return self
        
    
    def __add__(self, b):
        """Combine two mixtures to form a third unique mixture
    c = a + b
    
Makes a copy of a and calls __iadd__() to execute the operation.
"""
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
        

    def update(self, units=None):
        """UPDATE - bring the mixture into agreement with the unit matter config entry
    Returns nothing.
If the mixture's contents are currently listed in a unit matter that 
disagrees with PYroMat's pm.config['unit_matter'] setting, the contents
are converted to the current unit matter.
"""
        if units is None:
            units = pm.config['unit_matter']
            
        if self._units != units:
            for subst,qty in self.items():
                pm.units.matter(qty, subst.mw(), from_units=self._units, to_units=units, inplace=True)        
            self._units = units
        

    def _sindex(self, index):
        """SINDEX - convert an item index into an array index
    sindex = m._sindex(index)
    
When the index is an integer, string, slice, or a data instance it is 
interpreted as indexing the constituent substances.  Strings are 
converted into substance ID strings, and data instances are matched 
against their id hash.

When the index is a tuple, the first element is interpreted as a 
substance index and the trailing indices are interpreted as quantity
indices.  The substance index is passed recursively to _sindex() to 
convert it to an integer.

If the substance does not appear to be in this mixture, _sindex() 
returns None.
"""
        # If the index also contains a quantity index, only operate on
        # the first (substance) index.
        if isinstance(index, tuple):
            # Force the first index to be an integer
            si = self._sindex(index[0])
            if si is None:
                return None
            return (si,) + index[1:]
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
            return None
        # If we're looking for anything else...
        else:
            try:
                return self._c.index(index)
            except ValueError:
                return None
            except:
                raise

    def _broadcast_shape(self, qshape, sshape=None):
        """BROADCAST_SHAPE - determine a new quantity shape
    shape = _broadcast_shape(qshape)
        OR
    shape = _broadcast_shape(qshape, sshape)
   
If the mixture is interacting with a new quantity, either due to insertion
or addition, _broadcast_shape() compares the new quantity's shape tuple
to the existing mixture quantity array to determine how they should be
broadcast together.  

Returns the new mixture shape tuple if the operation is to be completed.

Raises an error if the shapes cannot be broadcast.

By default, the current mixture shape is used, but this can be overridden
by manually passing a second shape tuple to sshape.

Note that this method replicates the behavior of numpy's broadcast_shapes
function.  Writing our own prevents needing to bump the numpy version 
requirement.
"""
        if sshape is None:
            sshape = self._q.shape[1:]
        Ns = len(sshape)
        Nq = len(qshape)
        
        if Ns > Nq:
            qshape += (Ns - Nq) * (1,)
        elif Nq > Ns:
            sshape += (Nq - Ns) * (1,)
            
        out = tuple()
        for s,q in zip(sshape,qshape):
            if s == q or s == 1:
                out += (q,)
            elif q == 1:
                out += (s,)
            else:
                raise pm.utility.PMParamError(f'IGTMix: cannot broadcast shapes: {qshape}, {sshape}')
        return out
        
    #
    # The property interface
    #
    def _argparse(self, *varg, **kwarg):
        """Parse the arguments supplied to an IGTMix property method
    T,p,d,X,mw = _argparse(*varg, **kwarg)

_ARGPARSE automatically applies the default temperature and pressure,
def_T or def_p, from the pyromat.config system to deal with unspecified
parameters.  All inputs are re-cast as numpy arrays of at least one 
dimension and inputs are automatically converted from the configured 
user units into kJ, kmol, m^3.

The returned variables are arrays of temperature, T, pressure, p, and 
the density, d.  Temperature will always be returned, and at least one
of p and d will be populated as well, but one of them may be None.  
_argparse decides which to populate based on what is most efficient.
"""
        # 1) Handle varg and kward and their defaults
        # 2) Apply the argument rules...
        #   2.1: All arguments must be legal
        #   2.2: Only 2 arguments
        #   2.3: T, e, and h may not be specified together
        #   2.4: d and v may not be specified together 
        # 3) Convert the arguments to arrays with dim 1 or greater
        # 4) Convert to standard units
        # 5) Process equivalent arguments
        #   5.1: replace v with d
        #   5.2: replace e or h with T
        # 6) Check for out-of-bounds on basic arguments (This was shifted to last)
        # 7) Case out the possible combinations
        # 8) Broadcast the arrays appropriately
        # 9) Calculate T,p,d
        

        # Fancy tool for tracking iteration issues
        debug = False
        
        # Calculate some preliminaries
        # We're going to need the temperature limits in Kelvin
        Tlow, Thigh = self.Tlim()
        Tlow = pm.units.temperature_scale(Tlow, to_units='K')
        Thigh = pm.units.temperature_scale(Thigh, to_units='K')
        # Let's go ahead and get the mole fraction array
        X = self.X(asarray=True)
        # Get the molecular weight array
        mw = self.mw(X=X)

        # 1) Handle varg and kwarg and apply defaults

        # If varg is specified, assign its values to T,p
        if len(varg) > 0:
            if 'T' in kwarg:
                raise pm.utility.PMParamError('T was specified both positionally and with a keyword.')
            kwarg['T'] = varg[0]
        if len(varg) > 1:
            if 'p' in kwarg:
                raise pm.utility.PMParamError('p was specified both positionally and with a keyword.')
            kwarg['p'] = varg[1]
        if len(varg) > 2:
            raise pm.utility.PMParamError('There are only two positional arguments: T, p.')

        # Count the number of arguments
        nargs = len(kwarg)
        if nargs == 1:
            if 'T' not in kwarg:
                kwarg['T'] = pm.config.def_T()
            else:
                kwarg['p'] = pm.config.def_p()
        elif nargs == 0:
            kwarg['T'] = pm.config.def_T()
            kwarg['p'] = pm.config.def_p()
        
        # 2) Apply the argument rules
        # Re-measure the number of arguments and use sets to enforce
        # the remaining rules
        nargs = len(kwarg)
        args = set(kwarg.keys())
        # A set of the legal arguments
        legal_args = set(['T','p','d','v','e','h','s'])
        
        # 2.1: There may only be 2 arguments
        if nargs>2:
            raise pm.utility.PMParamError(
                    'Specifying more than two simultaneous parameters is illegal.')
        
        # 2.2: All arguments must be "legal" recognized arguments
        these_args = args - legal_args
        if these_args:
            message = 'Unrecognized propert(y/ies):'
            prefix = '  '
            for name in these_args:
                message += prefix + name
                prefix = ', '
            raise pm.utility.PMParamError(message)
        
        # 2.3: T, e, and h may not be specified together
        these_args = set(['T', 'e', 'h']).intersection(args)
        if len(these_args) > 1:
            message = 'Properties may not be specified together:'
            prefix = ' '
            for name in inverse_args:
                message += prefix + name
                prefix = ', '
            raise pm.utility.PMParamError(message)
        
        # 2.4: Density and specific volume cannot be specified together
        if 'v' in args and 'd' in args:
            raise pm.utility.PMParamError('Density (d) and specific volume (v) cannot be specified together.')
        
        # 3) Convert all arguments to numpy arrays
        #    The asarray function does NOT copy the array if it is already
        #    a numpy array.
        for name,value in kwarg.items():
            value = np.asarray(value, dtype=float)
            if value.ndim == 0:
                value = np.reshape(value, (1,))
            kwarg[name] = value
        
        # 4) Convert the units appropriately
        #   This step will only make a copy of the array if the units need
        #   to be converted.  Otherwise, the array is passed through verbatim
        #   As a result, the input array will ONLY be copied if it needs to
        #   be reshaped, converted, or retyped.
        # 5) Replace v with d and e/h with T
        if 'T' in kwarg:
            kwarg['T'] = pm.units.temperature_scale(kwarg['T'], to_units='K')
        if 'p' in kwarg:
            kwarg['p'] = pm.units.pressure(kwarg['p'], to_units='Pa')
        if 'd' in kwarg:
            value = pm.units.volume(kwarg['d'], to_units='m3', exponent=-1)
            kwarg['d'] = pm.units.matter(value, mw, to_units='kmol')
        if 'v' in kwarg:
            # Convert and replace with d at the same time
            value = pm.units.volume(kwarg['v'], to_units='m3')
            kwarg['d'] = 1./pm.units.matter(value, mw, to_units='kmol', exponent=-1)
        if 'h' in kwarg:
            value = kwarg['h']
            value = pm.units.energy(value, to_units='kJ')
            value = pm.units.matter(value, mw, to_units='kmol', exponent=-1)
            shape = self._broadcast_shape(value.shape)
            h = np.broadcast_to(value,shape)
            T = np.full(shape, 0.5*(Tlow+Thigh))
            I = np.ones(shape, dtype=bool)
            self._iter1('_h', h, T, I, Tlow, Thigh, X=X)
            kwarg['T'] = T
        if 'e'  in kwarg:
            value = kwarg['e']
            value = pm.units.energy(value, to_units='kJ')
            value = pm.units.matter(value, mw, to_units='kmol', exponent=-1)
            shape = self._broadcast_shape(value.shape)
            e = np.broadcast_to(value,shape)
            T = np.full(shape, 0.5*(Tlow+Thigh))
            I = np.ones(shape, dtype=bool)
            self._iter1('_e', e, T, I, Tlow, Thigh, X=X)
            kwarg['T'] = T
        if 's' in kwarg:
            value = kwarg['s']
            value = pm.units.energy(value, to_units='kJ')
            value = pm.units.matter(value, mw, to_units='kmol', exponent=-1)
            value = pm.units.temperature(value, to_units='K', exponent=-1)
            kwarg['s'] = value

        # Convert R into J/kmol/K - use this for p = dRT
        # Do NOT use this for s, e, and h relationships
        R = 1000 * pm.units.const_Ru

        T = p = d = None
        # Entropy requires special iteration
        if 's' in kwarg:
            # If density is specified
            if 'd' in args:
                raise pm.utility.PMParamError('(s,d) is not currently supported only for IGTMix classes')
            # If pressure is specified
            elif 'p' in kwarg:
                shape = self._broadcast_shape(kwarg['s'].shape)
                shape = self._broadcast_shape(kwarg['p'].shape, shape)
                #s = np.broadcast_to(kwarg['s'], shape)
                p = np.broadcast_to(kwarg['p'], shape)
                # adjust entropy to the reference pressure
                s = kwarg['s'] + pm.units.const_Ru * np.log(p / self._pref()) - self._smix()
                T = np.full_like(s, 0.5*(Tlow+Thigh))
                I = np.ones_like(s,dtype=bool)
                self._iter1('_s', s, T, I, Tlow, Thigh, X=X)
            # If temperature is specified
            elif 'T' in kwarg:
                shape = self._broadcast_shape(kwarg['s'].shape)
                shape = self._broadcast_shape(kwarg['T'].shape, shape)
                s = np.broadcast_to(kwarg['s'], shape)
                T = np.broadcast_to(kwarg['T'], shape)
                # Calculate the reference entropy at the specified temperature
                s0 = self._propeval('_s', T)[0] + self._smix()
                p = self._pref() * np.exp((s0-s)/pm.units.const_Ru)
                # Otherwise, this is an illegal combination!
            else:
                raise pm.utility.PMParamError(f'Cannot simultaneously specify parameters: {args}')
        # If temperature is specified
        elif 'T' in kwarg:
            # There isn't much work to do
            if 'p' in kwarg:
                shape = self._broadcast_shape(kwarg['T'].shape)
                shape = self._broadcast_shape(kwarg['p'].shape, shape)
                T = np.broadcast_to(kwarg['T'], shape)
                p = np.broadcast_to(kwarg['p'], shape)
            elif 'd' in kwarg:
                shape = self._broadcast_shape(kwarg['T'].shape)
                shape = self._broadcast_shape(kwarg['d'].shape, shape)
                T = np.broadcast_to(kwarg['T'], shape)
                d = np.broadcast_to(kwarg['d'], shape)
            else:
                message = 'Please report a bug: Unhandled event [T] in IGTMix._argparse with args:'
                prefix = ' '
                for name in kwarg:
                    message += prefix + name
                    prefix = ', '
                raise pm.utility.PMParamError(message)
        # If pressure is specified
        elif 'p' in kwarg:
            if 'd' in kwarg:
                shape = self._broadcast_shape(kwarg['p'].shape)
                shape = self._broadcast_shape(kwarg['d'].shape, shape)
                p = np.broadcast_to(kwarg['p'], shape)
                d = np.broadcast_to(kwarg['d'], shape)
                T = p / (R * d)
            else:
                message = 'Please report a bug: Unhandled event [p] in IGTMix._argparse with args:'
                prefix = ' '
                for name in kwarg:
                    message += prefix + name
                    prefix = ', '
                raise pm.utility.PMParamError(message)
        else:
            message = 'Please report a bug: Unhandled event [MASTER] in IGTMix._argparse with args:'
            prefix = ' '
            for name in kwarg:
                message += prefix + name
                prefix = ', '
            raise pm.utility.PMParamError(message)
            
        # Test the temperatures for out-of-bounds
        I = np.logical_or(T < Tlow, T > Thigh)
        if I.all():
            raise pm.utility.PMParamError('All of the specified states were out-of-bounds.  '  
                    'Legal temperatures are between {} and {} Kelvin.'.format(Tlow, Thigh))
        elif I.any():
            T[I] = pm.config['def_oob']
            pm.utility.print_warning('Some of the states were out of bounds - setting to config[\'def_oob\'].  '
                    'Legal temperatures are between {} and {} Kelvin.'.format(Tlow, Thigh))
        return T,p,d,X,mw
        
        
    #
    # Indexing
    #
    def __getitem__(self, sid):
        index = self._sindex(sid)
        if index is None:
            raise pm.utility.PMParamError('Substance is not in the mixture: ' + repr(index))
        return self._q[index]
        
    def __setitem__(self, sid, value):
        index = self._sindex(sid)
        # Are we inserting a new substance?
        if index is None:
            self.insert(sid,value)
            return
        # Otherwise, use the numpy broadcasting rules
        self._q[index] = value
        
    def __contains__(self, index):
        return self._sindex(index) is not None
        
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
                out += f'[...]{sid.hill()} + '
            out += f'[...]' + self._c[-1].hill()
        else:
            for sid,qty in zip(self._c[:-1], self._q[:-1]):
                out += f'{qty}{sid.hill()} + '
            out += f'{self._q[-1]}' + self._c[-1].hill()
        out += f' ({self._units})'
        return out
    
    def shape(self):
        """SHAPE - return the dimensions of the quantity arrays
    shape = m.shape()
    
Returns a tuple indicating the shape of each of the substance quantity 
arrays.  This also indicates the total number of individual mixtures
represented in this mixture instance.
"""
        return self._q.shape[1:]

    def reshape(self, newshape):
        """RESHAPE - change the shape of the mixture quantity arrays
    m.reshape( shape_tuple )
    
Like the Numpy reshape() method, the tuple indicates the dimensions of
the new mixture array.
"""
        newshape = (self._q.shape[0],) + tuple(newshape)
        self._q = self._q.reshape(newshape)

    def nsubst(self):
        """NSUBST - return the number of constituent substances
    N = m.nsubst()
    
The integer number of substances in the mixture (may be zero).
"""
        return len(self._c)
        
    def insert(self, sid, qty):
        """INSERT - add a new substance to the mixture
    m.insert(sid, qty)
    
sid     The substance identifier may be a string ('H2' or 'ig.H2') or a
        PYroMat substance instance.
        
qty     The quantity array may be a scalar or array-like that is 
        broadcastable to the existing mixture's quantity array.
"""
        if isinstance(sid,str):
            if '.' not in sid:
                sid = 'ig.' + sid
            # Go find the substance
            # This will raise a meaningful error if it is not found
            subst = pm.get(sid)
            
        # If this is a PM data class
        elif isinstance(sid,pm.reg.__basedata__):
            # Verify that it's an ideal gas first
            if sid.collection() != 'ig':
                raise pm.utility.PMParamError('IGTMix constituent was not an ideal gas: ' + sid.sid())
            subst = sid
        # IGTMix components cannot be nested.
        elif isinstance(sid,IGTMix):
            raise pm.utility.PMParamError('IGTMix instances cannot be nested. Use the + operator instead to combine mixtures.')
        else:
            raise pm.utility.PMParamError('IGTMix.insert(): Unrecognized constituent: ' + repr(sid))
        
        # Next, test for prior existence
        for this in self._c:
            if this is subst:
                raise pm.utility.PMParamError('IGTMix.insert(): substance is already in the mixture. Use the + operator to modify.')
        
        # Broadcast the quantity
        qty = np.atleast_1d(qty)
        newshape = self._broadcast_shape(qty.shape)
        newq = np.empty((self.nsubst() + 1,) + newshape, dtype=float)
        newq[-1] = qty
        newq[:-1] = self._q
        
        # All is well, so append
        self._c.append(subst)
        self._q = newq
    
    def remove(self, sid):
        """REMOVE - remove a substance from the mixture
    m.remove(sid)
    
sid     The substance identifier may be a string ('H2' or 'ig.H2') or a
        PYroMat substance instance.  If present, this substance will be
        removed from the mixture.  Otherwise, an error is raised.
"""
        index = self._sindex(sid)
        if index is None:
            raise pm.utility.PMParamError('IGTMix.remove(): substance is not in the mixture: ' + repr(sid))
        del self._c[index]
        self._q = np.delete(self._q, index, axis=0)
    
    #
    # Type manipulation
    #
    def toigmix(self,sid):
        """TOIGMIX - convert the mixture to a static IGMIX instance
    [... igmix list ...] = m.toigmix(sid)
    
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
        
    def tolist(self):
        """TOLIST - generate a list of the mixture contents
    mixlist = m.tolist()

The list entries are the substance id strings for each of the mixture's
contents.
"""
        return [this.sid() for this in self._c]
        
    def todict(self):
        """TODICT - generate a dictionary representing the mixture
    mixdict = m.todict()
    
The dictionary keys are the substance id strings for each of the mixture's
contents.  The values are arrays representing the quantities of each 
substance.  Note that these are exported in the units currently used by
the mixture - not necessarily the currently configured PYroMat units.
"""
        return {subst.sid():qty for subst,qty in self.items()}
        
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
            
    def X(self, asarray=False):
        """X - calculate the mole fractions of the substances
    x = m.X()
    
returns a dictionary with the substance id strings as the keys and the
fraction array as values.  Values are always between zero and one.
"""
        # Make a copy of _q
        x = np.array(self._q)
        # If this is in mass-based units, we'll need to convert to molar
        if pm.units.ismass(self._units):
            # Broadcasting rules should work just fine 
            mw = np.array([subst.mw() for subst in self])
            x /= mw.reshape((mw.size,) + (x.ndim-1)*(1,))
        x /= np.sum(x,axis=0)
        # If we want an array, we're done
        if asarray:
            return x
        # Convert to dict
        return {subst.sid():xx for subst,xx in zip(self, x)}
        
    def Y(self, asarray=False):
        """Y - calculate the mass fractions of the substances
    y = m.Y()
    
returns a dictionary with the substance id strings as the keys and the
fraction array as values.  Values are always between zero and one.
"""
        # Make a copy of _q
        y = np.array(self._q)
        # If this is in mass-based units, we'll need to convert to molar
        if not pm.units.ismass(self._units):
            # Broadcasting rules should work just fine 
            mw = np.array([subst.mw() for subst in self])
            y *= mw.reshape((mw.size,) + (x.ndim-1)*(1,))
        y /= np.sum(y,axis=0)
        # If we want an array, we're done
        if asarray:
            return y
        # Convert to dict
        return {subst.sid():yy for subst,yy in zip(self, y)}
        
        
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
    # Internal Method for Properties
    #
    def _propeval(self, prop, T, X=None, diff=False):
        """_PROPEVAL - generic evaluator for internal properties
    pval, pvalT = m._propeval(T, diff=False)
        OR
    pval, pvalT = m._propeval(T, X=X)
    
T       A numpy array of temperatures in units Kelvin.  
prop    The string name of the property method to call on members
diff    If True, the property derivative with respect to temperature is
        also evaluated.  If set to None, then the property has no 
        option for a differential output.
X       Optional argument allows the mole fraction array (not dict) to
        be passed to prevent redundant calculations.

Returns pval, the substance property values averaged by mole fraction,
and pvalT, the substance property's derivative with respect to 
temperature.
"""
        pval = 0
        pvalT = None;
        if diff:
            pvalT=0
        if X is None:
            X = self.X(asarray=True)
        if diff is None:
            for subst,x in zip(self, X):
                pmethod = getattr(subst, prop)
                this = pmethod(T)
                pval += x*this
        else:
            for subst,x in zip(self, X):
                pmethod = getattr(subst, prop)
                this, thisT = pmethod(T,diff=diff)
                pval += x*this
                if diff:
                    pvalT += x*thisT
        return pval,pvalT

    def _pref(self, X=None):
        """_PREF - Calculate the reference pressure in Pa
    pref = m._pref()
        OR
    pref = m._pref(X=X)
    
Optional keyword, X, allows the mole fraction array (not dict) to be 
passed to prevent redundant calculations.
    
If all substances in the mixture share the same reference pressure, this
simply returns a scalar.  Otherwise, _pref() calculates a reference 
pressure array with the same shape as the mixture array.
"""
        # Start a reference pressure array
        pref_array = []
        pref = None     # Track a scalar in case they are all the same
        same = True     # Keep a flag to indicate if the reference pressures
                        # disagree.
        for subst in self._c:
            pmclass = subst.pmclass()
            if pmclass == 'ig':
                this = subst._pref_pa
            elif pmclass == 'ig2':
                this = subst.data['pref']
            elif pmclass == 'igmix':
                subst._bootstrap()
                this = subst._pref_pa
            else:
                raise pm.utility.PMDataError('IGMIX: Cannot determine the reference pressure of PYroMat data class: ' + pmclass)
            
            pref_array.append(this)
            if pref is None:
                pref = this
            # If the values disagree, give up on the scalar
            elif pref != this:
                same = False
        
        # If all the reference pressures were the same, just return a 
        # scalar.
        if same:
            return pref
        
        # If at least one reference pressure was different, we'll need
        # the logarithmic weighted average.  Construct a broadcastable 
        # pressure reference array
        pref_array = np.array(pref_array).reshape((self.nsubst(),) + (self._q.ndim-1)*(1,))
        if X is None:
            X = self.X(asarray=True)
        pref = np.exp(np.sum(X * np.log(pref_array), axis=0))
        return pref
        
    def _smix(self, X=None):
        """_SMIX - calculate the entropy of mixing in kJ/kmol/K
    smix = m._smix()
        OR
    smix = m._smix(X=X)
    
Optional keyword, X, allows the mole fraction array (not dict) to be 
passed to prevent redundant calculations.
"""
        if X is None:
            X = self.X(asarray=True)
        temp = np.log(X)
        temp[X==0] = 0.     # Force log(0) to zero to prevent nan problems
        return -pm.units.const_Ru * np.sum(X * temp, axis=0)
        

    def _iter1(self, prop, y, T, Ids, Tmin, Tmax,
                ep=1e-6, Nmax=20, verbose=False,
                X=None):
        """Invert the _propeval method (NOT STANDARD ITER1!!!)
        
    _iter1(prop, y, T, Ids, xmin, xmax)
    
Iteration is performed in-place on the x array.

*** Required Parameters ***
prop        The string keyword index of the property to be calculated.
y           An array of N target values for fn().
T           The result array.  It should be an N-element floating point
            array that has already been initialized.
Ids         A down-select boolean index array; only x[Ids],y[Ids] will 
            be evaluated.  This allows iteration in-place on data sets 
            where only a portion of the data require iteration.  If y is
            a floating point array with N elements, Ids must be a bool
            array with N elements.  It will specify a down-selected 
            data set with M elements, where M<=N.
Tmin, Tmax  Upper and lower limit arrays for the T values.  These must
            broadcastable to match T and y.  Even values outside of the
            down-select region should have legal values.
*** Optional Parameters ***
ep          Epsilon; fractional error permitted in y (default 1e-6)
Nmax        Maximum number of iterations (default 20)
X           The mol fraction array, preventing its repeated calculation

"""
        # As the iteration progresses, the number of True elements in 
        # Ids will decrease until they are all false
        # There are some important intermediate values that will also
        # require indexable arrays
        dT = np.zeros_like(y, dtype=float)
        error = np.zeros_like(dT, dtype=float)
        IooB = np.zeros_like(Ids, dtype=bool)

        count = 0

        if verbose:
            print('x, yy, yyx, dx, Ids')
        while Ids.any():
            # Evaluate the funciton and isolate its derivative
            yy,yyT = self._propeval(prop, T, diff=True, X=X)
            # note that x[Ids], yy, yyx, and all the other floating 
            # intermediates are now in m-space; the sub-set of values
            # still under iteration.
            # Calculate the error, the linear change in x, and the new x
            error[Ids] = y[Ids] - yy
            dT[Ids] = error[Ids] / yyT
            if verbose:
                print(T, yy, yyT, dT, Ids)
            T[Ids] += dT[Ids]
            # An out-of-bounds index
            IooB = np.logical_or( T < Tmin, T > Tmax)
            count_oob = 0
            while IooB.any():
                dT[IooB] /= 2.
                T[IooB] -= dT[IooB]
                IooB = np.logical_or( T < Tmin, T > Tmax)
                # Prevent a while-loop-trap
                count_oob += 1
                if count_oob>Nmax:
                    raise pm.utility.PMAnalysisError(
                        'iter1_() failed to produce a guess that was in-bounds')
            
            # Check the iteration convergence
            Ids[Ids] = abs(error[Ids]) > abs(ep*y[Ids])
            # Prevent a while-loop-trap
            count += 1
            if count>Nmax:                
                pm.utility.print_warning(\
                    'iter1_() failed to converge for %d elements after %d attempts'%(\
                    Ids.sum(), Nmax))
                return

    #
    # Properties
    #
    
    def mw(self, *varg, X=None, **kwarg):
        """MW - molecular weight

    mw = m.mw()
        OR
    mw = m.mw(X=X)
    
Returns the molecular weight in [unit_mass] / [unit_molar].

The optional keyword, X, allows passing the mole fraction array (not dict)
to prevent redundant calculations.
"""
        if X is None:
            X = self.X(asarray=True)
        mw = np.array([subst.data['mw'] for subst in self._c]).reshape((X.shape[0],) + (X.ndim-1)*(1,))
        out = np.sum(mw * X, axis=0)
        pm.units.mass(out, from_units='kg', inplace=True)
        pm.units.molar(out, from_units='kmol', inplace=True, exponent=-1)
        return out
        
    def R(self, *varg, X=None, **kwarg):
        """R - ideal gas constant
    R = m.R()
    
Returns the gas constant in [unit_energy]/[unit_matter][unit_temperature]
"""
        # Convert the universal constant to the configured units
        R = pm.units.energy(pm.units.const_Ru, from_units='kJ')
        R = pm.units.temperature(R, from_units='K', exponent=-1)
        R = pm.units.molar(R, from_units='kmol', exponent=-1)
        # If we need to switch into mass units, call mw()
        if pm.units.ismass(pm.config['unit_matter']):
            R /= self.mw(X=X)
        return R

    def Tlim(self):
        """TLIM - Calculate the temperature limits on the model
    Tlow, Thigh = m.Tlim()
    
The high and low temperatures are scalar values that apply for all 
individual mixtures in the quantity arrays.
"""
        Tlow = float('inf')
        Thigh = float('-inf')
        for subst in self:
            this_low, this_high = subst.Tlim()
            if this_low < Tlow:
                Tlow = this_low
            if this_high > Thigh:
                Thigh = this_high
        return Tlow, Thigh
    
    def cp(self, *varg, **kwarg):
        """CP - Constant pressure specific heat
        
"""
        T,p,d,X,mw = self._argparse(*varg, **kwarg)
        cp = self._propeval('_cp', T, X=X, diff=None)[0]
        pm.units.energy(cp, from_units='kJ', inplace=True)
        pm.units.temperature(cp, from_units='K', inplace=True, exponent=-1)
        pm.units.matter(cp, mw, from_units='kmol', inplace=True, exponent=-1)
        return cp
        
    def cv(self, *varg, **kwarg):
        """CV - Constant volume specific heat
        
"""
        T,p,d,X,mw = self._argparse(*varg, **kwarg)
        cv = self._propeval('_cp', T, X=X, diff=None)[0]
        cv -= pm.units.const_Ru
        pm.units.energy(cv, from_units='kJ', inplace=True)
        pm.units.temperature(cv, from_units='K', inplace=True, exponent=-1)
        pm.units.matter(cv, mw, from_units='kmol', inplace=True, exponent=-1)
        return cv
        
    def gam(self, *varg, **kwarg):
        """GAM - Specific heat ratio
        
"""
        T,p,d,X,mw = self._argparse(*varg, **kwarg)
        cp = self._propeval('_cp', T, X=X, diff=None)[0]
        return cp / (cp - pm.units.const_Ru)
        
    def h(self, *varg, **kwarg):
        """H - Enthalpy
        
"""
        T,p,d,X,mw = self._argparse(*varg, **kwarg)
        h = self._propeval('_h', T, X=X)[0]
        pm.units.energy(h, from_units='kJ', inplace=True)
        pm.units.matter(h, mw, from_units='kmol', inplace=True, exponent=-1)
        return h
        
    def e(self, *varg, **kwarg):
        """E - Internal energy
        
"""
        T,p,d,X,mw = self._argparse(*varg, **kwarg)
        e = self._propeval('_e', T, X=X)[0]
        pm.units.energy(e, from_units='kJ', inplace=True)
        pm.units.matter(e, mw, from_units='kmol', inplace=True, exponent=-1)
        return e
        
    def s(self, *varg, **kwarg):
        """S - Entropy
        
"""
        T,p,d,X,mw = self._argparse(*varg, **kwarg)
        # We need pressure
        if p is None:
            p = d * (1000*pm.units.const_Ru * T)
        s = self._propeval('_s', T, X=X)[0]
        s += self._smix(X=X)
        s -= pm.units.const_Ru * np.log(p / self._pref(X=X))
        pm.units.energy(s, from_units='kJ', inplace=True)
        pm.units.matter(s, mw, from_units='kmol', inplace=True, exponent=-1)
        pm.units.temperature(s, from_units='K', inplace=True, exponent=-1)
        return s
        
    def g(self, *varg, **kwarg):
        """G - Gibbs energy
        
"""
        T,p,d,X,mw = self._argparse(*varg, **kwarg)
        # We need pressure
        if p is None:
            p = d * (1000*pm.units.const_Ru * T)
        g = self._propeval('_g', T, X=X)[0]
        g -= T * (self._smix(X=X) - pm.units.const_Ru * np.log(p / self._pref(X=X)))
        pm.units.energy(g, from_units='kJ', inplace=True)
        pm.units.matter(g, mw, from_units='kmol', inplace=True, exponent=-1)
        return g

    def f(self, *varg, **kwarg):
        """F - Free (Helmholtz) energy
        
"""
        T,p,d,X,mw = self._argparse(*varg, **kwarg)
        # We need pressure
        if p is None:
            p = d * (1000*pm.units.const_Ru * T)
        f = self._propeval('_f', T, X=X)[0]
        f -= T * (self._smix(X=X) - pm.units.const_Ru * np.log(p / self._pref(X=X)))
        pm.units.energy(f, from_units='kJ', inplace=True)
        pm.units.matter(f, mw, from_units='kmol', inplace=True, exponent=-1)
        return f

    def T(self, *varg, **kwarg):
        """T - Temperature
        
"""
        T,p,d,X,mw = self._argparse(*varg, **kwarg)            
        pm.units.temperature_scale(T, from_units='K', inplace=True)
        return T
    
    def p(self, *varg, **kwarg):
        """P - Pressure
        
"""
        T,p,d,X,mw = self._argparse(*varg, **kwarg)
        if p is None:
            p = d * (1000*pm.units.const_Ru * T)
            
        pm.units.pressure(p, from_units='bar', inplace=True)
        return p

    def d(self, *varg, **kwarg):
        """D - Density
        
"""
        T,p,d,X,mw = self._argparse(*varg, **kwarg)
        if d is None:
            d = p / (1000*pm.units.const_Ru * T)
            
        pm.units.matter(d, mw, from_units='kmol', inplace=True)
        pm.units.volume(d, from_units='m3', inplace=True, exponent=-1)
        return d
        
        
        

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
