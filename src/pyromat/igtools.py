"""IGTOOLS - A module for ideal gas tools

::Tools for ideal gas mixtures::

BETA RELEASE v2.4.5
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

class IgtMix(object):
    """Ideal Gas Tools dynamic Mixture
    
    m = IgtMix( ... )
    
The IgtMix class instances allow users to dynamically define mixtures of
component ideal gases with a flexible interface convenient for the 
command line and scripts alike.  IgtMix instances have property methods 
just like the rest of the PYroMat classes, but the total quanitity 
(the extent) of the mixture is specified.  Therefore, properties are
calculate in total and not per unit metter (not intensive).

In this example, we quickly compute the enthalpy of a mixture 35% 
hydrogen and balance argon by volume.
    
    mymix = IgtMix('0.35 H2 + .65 Ar', units='kmol')
    mymix.h(T=300)
        array([43.67153007])

** METHODS **

IgtMix instances provide the same properties as the core PYroMat classes

  cp() spec. heat       (unit_energy / unit_temperature)
  cv() spec. heat       (unit_energy / unit_temperature)
  d()  density          (unit_matter / unit_volume)
  e()  internal energy  (unit_energy)
  f()  free energy      (unit_energy)
  g()  Gibbs energy     (unit_energy)
  gam()  spec. heat ratio (dless)
  h()  enthalpy         (unit_energy)
  mw() molecular weight (unit_mass / unit_molar)
  T()  temperature      (unit_temperature)
  p()  pressure         (unit_pressure)
  R()  gas constant     (unit_energy / unit_matter / unit_temperature)
  s()  entropy          (unit_energy / unit_temperature)
  v()  specific volume  (unit_volume / unit_matter)
There is not currently a state() method, but there are plans to add one.

For information about the mixture, see
  atoms()   qty. of each element    (unit_molar)
  mass()    total mass              (unit_mass)
  molar()   total number of moles   (unit_molar)
  Tlim()    lower & upper limits    (unit_temperature)
  X()       mole fractions          (dless)
  Y()       mass fractions          (dless)

To operate on the mixture like a Numpy array, see
  insert()      Add a new substance
  nsubst()      Number of substances in the mixture
  ndim()        Mixture array dimension
  remove()      Remove a substance
  reshape()     Change the mixture array's shape
  shape()       Mixture array shape
See also the indexing features described below.

To convert to other data types, see
  todict()      A dictionary with keys=subst, values=quantities
  toigmix()     List of igmix instances
  tolist()      Return a list of the substances
There is also afromigmix() method.


** DEFINING A MIXTURE **

There are six ways to define a mixture:

(1) From a string...
Strings are expected in the format: 'qty0 subst0 + qty1 subst1 + ...'
Quantities specify an amount of each substance in the units configured
in pm.config['unit_matter'] unless the optional 'units' keyword is set.
Omitted quantities are presumed to be unity, and substances can be 
specified with or without their 'ig.' collection prefix.  All whitespace
is ignored.  

    mymix = IgtMix('10 ig.N2')
    mymix = igt.IgtMix('2.4 N2 + Ar')
    print(mymix)
        [2.4]N2 + [1.]Ar

*note* This method cannot be used to specify an ion like 'Ne+'. Instead,
use the dictionary or list methods below.

(2) From keyword arguments...
Keyword arguments accept abbreviated substance ID strings as keywords 
(see below) under "SUBSTANCE ID".
    
    air = IgtMix(N2=0.76, O2=0.23, Ar=0.01)

*note* This method cannot be used to specify an ion like 'Ne+'. Instead,
use the dictionary or list methods below.

(3) A list of constituents with no quantities...
This assigns zero matter to each mixture.  Then, users are free to 
modify the mixture contents by indexing.

    mymix = IgtMix(['ig.N2', 'ig.O2', 'ig.Ar'])
    mymix['N2'] = 0.76
    mymix['O2'] = 0.23
    mymix['Ar'] = 0.01
    
(4) A dictionary of constituents with their quantities as values...
Dictionary keys are interpreted to be substance identifiers, and the 
corresponding values are interpreted as quantities.

    air = IgtMix({'ig.N2':0.76, 'ig.O2':0.23, 'ig.Ar':.01})

(5) Algebraicaly...
IgtMix instances can be combined with each other using addition and 
multiplication at the command line.

    air = IgtMix('.76 N2 + .23 O2 + .01 Ar')
    fuel = IgtMix('.23 CH4 + .44 C3H8')
    reactants = air + 0.4*fuel
    print(reactants)
[0.76]N2 + [0.23]O2 + [0.01]Ar + [0.092]CH4 + [0.176]C3H8

Any other data type in addition with an IgtMix is interpreted as a 
mixture definition as well.  For example,

    reactants = '.76 N2 + .23 O2 + .01 Ar' + 0.4*fuel       
    reactants = {'N2':.75, 'O2':.23, 'Ar':.01} + 0.4*fuel   
    mymix = pm.get('ig.H2O') + 0.4 * air                    

(6) From an existing igmix instance, using the fromigmix() function.

    air = fromigmix(pm.get('ig.air'))
    print(air)
[0.00935019]Ar + [0.00031401]CO2 + [0.78085562]N2 + [0.20948019]O2
    
This is a useful trick if users want to make quick changes to an 
existing mixture.  When the fromigmix()function is used, the igmix 
instance is split into its constituent gases to form the IgtMix 
instance, so the example above results in a mixture with argon, carbon 
dioxide, nitrogen, and oxygen.  If an igmix instance is passed directly 
to the IgtMix class, it is treated the same as any other constituent 
gas, so the example below produces a mixture that only has one gas 
component.

    air = igt.IgtMix(pm.get('ig.air'))
    print(air)
[1.]air

** SPECIFYING SUBSTANCES **

In the examples above, constituent gases are specified by their 
substance ID string, but they may be specified three ways:
(1) By their full substance ID string

    mymix = IgtMix('ig.N2')
    
(2) By their abbreviated substance ID string (with no leading ig.)

    mymix = IgtMix('N2')
    
(3) Or by their full data instance

    n2 = pm.get('ig.N2')
    mymix = IgtMix(n2)
    
Any of these -- 'N2', 'ig.N2', or pm.get('ig.N2') -- may be used as the
substance identifier in the mixture list or dictionary methods above.

** ARITHMETIC WITH MIXTURES **

Because mixtures work in absolute quantities (as opposed to mass or 
mole fractions), they can be incrased, decreased, subtracted from, or
added to using basic math operations at the command line.  In this 
example, mix3 has 1kmol of water, 2 kmol carbon dioxide, and 2 kmol each
of nitrogen and argon.

    pm.config['unit_matter'] = 'kmol'
    mix1 = IgtMix({'H2O':2, 'CO2':4})
    mix2 = IgtMix({'N2':1, 'Ar':1})
    mix3 = 0.5*mix1 + 2*mix2

Addition and subtraction attempt to convert non-IgtMix instances to 
IgtMix instances.  That means that strings, lists, dictionaries, and
ordinary PYroMat instances may be folded into mixtures using plain 
command-line arithmetic.

    mix4 = mix2 + '0.5 CO'
    mix5 = mix2 + pm.get('ig.H2O')
    mix6 = mix2 + {'H2O':[0.5,0.12], 'C2H2':0.8}

The last examples show how mixture arrays can be defined.  Array broad-
casting is inherently supported, so even though H2O was the only 
substance with multiple values, all other values are broadcast,

    print(mix6)
        N2   : array([1., 1.])
        Ar   : array([1., 1.])
        H2O  : array([0.5 , 0.12])
        C2H2 : array([0.8, 0.8])

** QUANTITIES AND UNTIS **

When an IgtMix instance is initialized, the quantities are interpreted
in PYroMat's currently configured 'unit_matter' unless the behavior is
overridden by the optional units keyword.

    pm.config['unit_matter'] = 'kg'
    mymix = IgtMix({'H2O':5.2, 'CO2':4.7}, units='kmol')
        [93.679456]H2O + [206.84606]CO2

Changes to the 'unit_matter' configuration parameter do not affect the
mixture's definition.  Instead, it only changes the way the mixture's 
quantities are represented in a print() operation.

** ARRAYS AND INDEXING **

Quantities are always handled as arrays, so a single IgtMix instance can
actually manage arrays of mixtures with an arbitrary shape.  The 
indexing rules are constructed around the idea that the IgtMix instance
is actually an array of many individual mixtures.  In this example, 
the mixture instance is an array of mixtures of neon, its first ion, and
the free electron.

    x = np.linspace(0,1,11)
    mymix = IgtMix{{'Ne':1-x, 'Ne+':x, 'e':x}, units='kmol')
    print(mymix)
Ne  : array([1. , 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0. ])
Ne+ : array([0. , 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. ])
e-  : array([0. , 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. ])

The shape() method returns a tuple like Numpy's shape attribute, but
it describes only the array of mixtures.  It excludes the axis with the
constituent substances.  To obtain the number of substances, use the 
nsubst() method.

    mymix.shape()
        (11,)
    mymix.nsubst()
        3

See also reshape().

IgtMix instances may be indexed like a normal Numpy array, but with the
special rule that the first dimension of the array is always indexed by
a substance ID.  From the example above,

    mymix['Ne']
        array([1. , 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0. ])
    mymix['Ne',4]
        0.6

When the first index is a slice, the value returned is a sub-mixture
instead of an array

    print(mymix[:,4])
        0.6Ne + 0.4Ne+ + 0.4e-
    print(mymix[:,3:6])
        Ne  : array([0.7 , 0.65, 0.5 ])
        Ne+ : array([0.3 , 0.35, 0.5 ])
        e-  : array([0.3 , 0.35, 0.5 ])

Assignment works as well.  The example below overwrites the fifth mix-
ture.

    mymix[:,4] = {'Ne':0.65, 'Ne+':0.35, 'e':0.35}
    print(mymix)
        Ne  : array([1.  , 0.9 , 0.8 , 0.7 , 0.65, 0.5 , 0.4 , 0.3 , ...
        Ne+ : array([0.  , 0.1 , 0.2 , 0.3 , 0.35, 0.5 , 0.6 , 0.7 , ...
        e-  : array([0.  , 0.1 , 0.2 , 0.3 , 0.35, 0.5 , 0.6 , 0.7 , ...

** ITERATING **

IgtMix instances act much like an ordered dictionary with PYroMat ideal
gas class isntances as keys and quantitiy arrays as values.  The items() 
method allows simultaneous iteration over the substances and quantities.

    for subst,qty in mymix.items():
        # subst is now a PYroMat substance instance
        # qty is now a Numpy array of the substance's quantity
        
    for subst in mymix:
        # subst is now a PYroMat substance instance

** REPRESENTATION **

The mixture is represented in the back-end by attributes that are not 
intended to be accessed directly.  Instead, these attributes should only
be read or modified by the appropriate methods.

_q      is a numpy array of quantities in kmol. The first index in the 
        array corresponds to the substance.  All other dimensions allow 
        the IgtMix instance to specify an arbitrary number of individual 
        mixtures.
_c      is a numpy array of ideal gas instances that make up the mixture.
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
        elif isinstance(contents, IgtMix):
            self._c = list(contents._c)
            self._q = np.array(contents._q)
            return
        # If the contents is a string or data instance, build a dummy dict
        elif isinstance(contents, str):
            return self.__init__(contents=parse_mixstr(contents), units=units)
        # If contents is a base PYroMat class
        elif isinstance(contents, pm.reg.__basedata__):
            return self.__init__(contents={contents:1}, units=units)
        # if the contents is a list or tuple of constituents
        elif isinstance(contents, (list, tuple)):
            return self.__init__(dict.fromkeys(contents,0), units=units)
        # Finally, if the contents is not a dict, raise an error
        elif not isinstance(contents, dict):
            raise Exception(pm.utility.PMParamError('Unexpected data type for IgtMix contents: ' + repr(type(contents))))
        
        
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
                    raise pm.utility.PMParamError('IgtMix constituent was not an ideal gas: ' + sid.sid())
                self._c.append(sid)
            # IgtMix components are included automatically
            elif isinstance(sid,IgtMix):
                raise pm.utility.PMParamError('IgtMix instances cannot be nested. Use the + operator instead to combine mixtures.')
            else:
                raise pm.utility.PMParamError('IgtMix: Unrecognized constituent: ' + repr(sid))

            # Next, process the quantity as an array
            qshape = np.shape(qty)
            shape = self._broadcast_shape(qshape, shape)
            
        # Initialize the quantity array
        N = len(self._c)
        self._q = np.empty((N,) + shape, dtype=float)

        # In the second loop, we'll broadcast all of the quantity arrays
        # to the correct shape and with the correct units
        for index,qty in enumerate(contents.values()):
            self._q[index,:] = pm.units.matter(qty, self._c[index].mw(), from_units=units, to_units='kmol')
            
        
    def __iter__(self):
        return self._c.__iter__()
        
    def __iadd__(self, b):
        """Add mixture b to self
    a += b
    
This combines the quantities of gas in mixture b into mixture a.
"""
        # If this is not an IgtMix, see if we can convert it to a mixture
        if not isinstance(b, IgtMix):
            return self.__iadd__(IgtMix(b))
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
             
        for bi,subst in enumerate(b._c):
            ai = self._c.index(subst)
            self._q[ai,:] += b._q[bi,:]
    
        return self
        
    
    def __add__(self, b):
        """Combine two mixtures to form a third unique mixture
    c = a + b
    
Makes a copy of a and calls __iadd__() to execute the operation.
"""
        c = IgtMix(self)
        c.__iadd__(b)
        return c
        
        
    def __radd__(self, b):
        c = IgtMix(self)
        c.__iadd__(b)
        return c
      
    def __neg__(self):
        """Create a mixture with negative quantities
    b = -a
"""
        b = IgtMix()
        b._c = self._c.copy()
        b._q = -self._q
        return b

      
    def __sub__(self, b):
        """Deduct the contents of one mixture from another
    c = a - b
    
Negates b and uses __iadd__() to perform the operation.
"""
        c = b.__neg__()
        c.__iadd__(self)
        return c
        
    def __isub__(self, b):
        self.__iadd__(b.__neg__())
        return self
        
    def __rsub__(self, b):
        c = self.__neg__()
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
        return IgtMix(self).__imul__(b)
        
    def __rmul__(self,b):
        return self.__mul__(b)
        
    def __itruediv__(self, b):
        return self.__imul__(1./b)
            
    def __truediv__(self, b):
        return IgtMix(self).__itruediv__(b)

    def _sindex(self, index):
        """SINDEX - convert an item index into an array index (inner routine)
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
            # Convert the first index, pass the others through
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
        """BROADCAST_SHAPE - determine a new quantity shape (inner routine)
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
        
        # In the special case that both are empty, force an empty 1d array
        if Ns == 0 and Nq == 0:
            return (0,)
        
        # If the shapes have dissimilar numbers of dimension, add the
        # necessary unity dimensions
        if Ns > Nq:
            qshape += (Ns - Nq) * (1,)
        elif Nq > Ns:
            sshape += (Nq - Ns) * (1,)
        
        # Check each dimension for compatibility
        out = tuple()
        for s,q in zip(sshape,qshape):
            if s == q or s == 1:
                out += (q,)
            elif q == 1:
                out += (s,)
            else:
                raise pm.utility.PMParamError(f'IgtMix: cannot broadcast shapes: {qshape}, {sshape}')
                
        return out
        
    #
    # The property interface
    #
    def _argparse(self, *varg, **kwarg):
        """Parse the arguments supplied to an IgtMix property method (inner routine)
    T,p,d = _argparse(*varg, **kwarg)

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
        legal_args = {'T','p','d','v','e','h','s'}
        
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
            for name in these_args:
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
            kwarg['d'] = pm.units.matter(value, self._mw(), to_units='kmol')
        if 'v' in kwarg:
            # Convert and replace with d at the same time
            value = pm.units.volume(kwarg['v'], to_units='m3')
            kwarg['d'] = 1./pm.units.matter(value, self._mw(), to_units='kmol', exponent=-1)
        if 'h' in kwarg:
            value = kwarg['h']
            value = pm.units.energy(value, to_units='kJ')
            shape = self._broadcast_shape(value.shape)
            h = np.broadcast_to(value,shape)
            T = np.full(shape, 0.5*(Tlow+Thigh))
            I = np.ones(shape, dtype=bool)
            self._iter1('_h', h, T, I, Tlow, Thigh)
            kwarg['T'] = T
        if 'e'  in kwarg:
            value = kwarg['e']
            value = pm.units.energy(value, to_units='kJ')
            shape = self._broadcast_shape(value.shape)
            e = np.broadcast_to(value,shape)
            T = np.full(shape, 0.5*(Tlow+Thigh))
            I = np.ones(shape, dtype=bool)
            self._iter1('_e', e, T, I, Tlow, Thigh)
            kwarg['T'] = T
        if 's' in kwarg:
            value = kwarg['s']
            value = pm.units.energy(value, to_units='kJ')
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
                raise pm.utility.PMParamError('(s,d) is not currently supported only for IgtMix classes')
            # If pressure is specified
            elif 'p' in kwarg:
                shape = self._broadcast_shape(kwarg['s'].shape)
                shape = self._broadcast_shape(kwarg['p'].shape, shape)
                #s = np.broadcast_to(kwarg['s'], shape)
                p = np.broadcast_to(kwarg['p'], shape)
                # adjust entropy to the reference pressure
                s = kwarg['s'] - self._sp(p)[0]
                T = np.full_like(s, 0.5*(Tlow+Thigh))
                I = np.ones_like(s,dtype=bool)
                self._iter1('_s', s, T, I, Tlow, Thigh)
            # If temperature is specified
            elif 'T' in kwarg:
                shape = self._broadcast_shape(kwarg['s'].shape)
                shape = self._broadcast_shape(kwarg['T'].shape, shape)
                s = np.broadcast_to(kwarg['s'], shape)
                T = np.broadcast_to(kwarg['T'], shape)
                # Calculate the pressure portion of entropy by deducting
                # the standard entropy and the "entropy of mixing"
                sp = kwarg['s'] - self._propeval('_s', T)[0]
                # Deduct the entropy of mixing
                qt = np.sum(self._q, axis=0)
                temp = self._q / qt
                # Force values with nonpositive q to have a log of zero
                temp[self._q <= 0] = 1.
                sp += pm.units.const_Ru * np.sum(self._q * np.log(temp), axis=0)
                # Calculate pressure
                p = self._pref() * np.exp( -sp / pm.units.const_Ru / qt)
                
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
                message = 'Please report a bug: Unhandled event [T] in IgtMix._argparse with args:'
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
                message = 'Please report a bug: Unhandled event [p] in IgtMix._argparse with args:'
                prefix = ' '
                for name in kwarg:
                    message += prefix + name
                    prefix = ', '
                raise pm.utility.PMParamError(message)
        else:
            message = 'Please report a bug: Unhandled event [MASTER] in IgtMix._argparse with args:'
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
        return T,p,d
        
    def _mw(self):
        """MW - molecular weight (inner routine)

    mw = m._mw()
    
Returns the molecular weight in kg / kmol.
"""
        mw = np.array([subst.data['mw'] for subst in self._c]).reshape((self._q.shape[0],) + (self._q.ndim-1)*(1,))
        return np.sum(mw * self._q, axis=0) / np.sum(self._q, axis=0)
        
        
    def _sp(self, p, diff=False):
        """SP - The pressure portion of entropy (inner routine)
    ss, ssp = m._sp(p, diff=False)
    
Accepts pressure in Pa.
"""
        # Calculate the partial pressure for each constituent
        qt = np.sum(self._q, axis=0)
        sp = np.zeros_like(p, dtype=float)
        spp = None
        if diff:
            spp = 0.
        # The different classes deal with reference pressure differently
        for subst,q in zip(self._c, self._q):
            if subst.data['class'] == 'ig':
                pref = subst._pref_pa
            elif subst.data['class'] == 'ig2':
                pref = subst.data['pref']
            elif subst.data['class'] == 'igmix':
                subst._bootstrap()
                pref = subst._pref_pa
            else:
                raise pm.utility.PMDataError('Substance is not a recongized IG class: ' + subst.data['id'] + ', ' + subst.data['class'])

            # Detect elements with negative or zero quantities
            I = q<=0
            # Partial pressure is p * q / qt
            temp = p * q / qt / pref
            temp[I] = 1.    # Force to a value that will not cause an error in log
            # Force q<=0 values to 
            sp -= q * np.log(temp)

        sp *= pm.units.const_Ru
        if diff:
            spp = pm.units.const_Ru * qt / p
        return sp, spp
        
    #
    # Indexing
    #
    def __getitem__(self, index):
        ii = self._sindex(index)
        if ii is None:
            raise pm.utility.PMParamError('Could not index mixture with argument: ' + repr(index))
        # Case out the four ways to index the mixture
        # If the index is a slice of multiple substances
        if isinstance(ii, slice):
            out = IgtMix()
            out._c = self._c[ii]
            out._q = self._q[ii]
            return out
        # If the index is a tuple of indices
        elif isinstance(ii, tuple):
            # If the first (substance) index is a slice
            if isinstance(ii[0],slice):
                out = IgtMix()
                out._c = self._c[ii[0]]
                out._q = self._q[ii]
                return out
            # If a specific species is indexed
            else:
                subst = self._c[ii[0]]
        # If the index is merely an integer
        else:
            subst = self._c[ii]

        return pm.units.matter(self._q[ii], subst.mw(), from_units='kmol')
        
        
    def __setitem__(self, sid, value):
        index = self._sindex(sid)
        # Attempt to insert a new substance?
        if index is None:
            raise pm.utility.PMParamError('Could not index the substance with: ' + repr(sid))
        # If we're slicing substances
        elif isinstance(index, slice) or\
                isinstance(index, tuple) and isinstance(index[0],slice):
            # Force the value to be a mixture
            if not isinstance(value, IgtMix):
                try:
                    value = IgtMix(value)
                except:
                    raise pm.utility.PMParamError('Could not interpret the assigned values to a mixture: ' + repr(value))
                # Recursively set the items one row at a time
            for ii,subst in enumerate(value._c):
                # Build a subindex
                if isinstance(index, tuple):
                    subindex = (subst,) + index[1:]
                else:
                    subindex = subst
                self.__setitem__(subindex, value._q[ii])

        else:
            # Otherwise, use the numpy broadcasting rules
            value = np.atleast_1d(value)
            shape = self._broadcast_shape(value.shape)
            # If the existing array needs to be reshaped
            if shape != self.shape():
                # Create a new array and broadcast the old data into it
                _q = np.empty((self._q.shape[0],) + shape)
                _q[:] = self._q
                self._q = _q
            # If the assigned array needs to be reshaped, broadcasting will
            # do it automatically
            self._q[index] = pm.units.matter(value, self._c[index].mw(), to_units='kmol')
        
    def __contains__(self, index):
        return self._sindex(index) is not None
        
    #
    # Length and size operations
    #
    def __len__(self):
        """LEN - detect the number of mixtures represented here
"""
        return np.prod(self.shape())
    
    def __str__(self):
        """STR - pretty print of the mixture
"""
        out = ''
        if self.__len__() > 1:
            # Measure the Hill string length
            hilllen = 0
            hills = []
            for sid in self._c:
                this = sid.hill()
                hilllen = max(hilllen, len(this))
                hills.append(this)
            
            for ii,this in enumerate(hills):
                out += this + (hilllen + 1 - len(this))*' '
                out += ': ' + repr(self._q[ii]) + '\n'
        else:
            for subst,qty in zip(self._c[:-1], self._q[:-1]):
                qty = pm.units.matter(qty, subst.mw(), from_units='kmol')
                out += f'{qty}{subst.hill()} + '
            subst = self._c[-1]
            qty = pm.units.matter(self._q[-1], subst.mw(), from_units='kmol')
            out += f'{qty}{subst.hill()}'
        return out
        
        
    def __repr__(self):
        """REPR - short command-line representation
"""
        return f'<IgtMix, {len(self._c)} subst, {self.shape()}>'
        
    
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
the new mixture array.  Reshape may change the number of dimensions of
the mixture array, but it may not change the total number of mixtures.
"""
        newshape = (self._q.shape[0],) + tuple(newshape)
        self._q = self._q.reshape(newshape)

    def ndim(self):
        """NDIM - return the number of dimensions of the quantity arrays
    ndim = m.ndim()

Returns an integer indicating the dimension of the mixture array.
"""
        return self._q.ndim-1

    def nsubst(self):
        """NSUBST - return the number of constituent substances
    N = m.nsubst()
    
The integer number of substances in the mixture (may be zero).
"""
        return len(self._c)
        
    def insert(self, sid, qty=0., units=None):
        """INSERT - add a new substance to the mixture
    m.insert(sid)
        OR
    m.insert(sid, qty)
        OR
    m.insert(sid, qty, units=unit_str)
    
sid     The substance identifier may be a string ('H2' or 'ig.H2') or a
        PYroMat substance instance.
        
qty     If no quantity is specified, then the substance is added with 
        zero quantity.  Otherwise, the quantity must be a scalar or an 
        array-like entity that is broadcastable to the mixture shape.
        
units   Optional units string specifies the units of the quantity being
        inserted.  When unspecified, uses the existing mixture units.
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
                raise pm.utility.PMParamError('IgtMix constituent was not an ideal gas: ' + sid.sid())
            subst = sid
        # IgtMix components cannot be nested.
        elif isinstance(sid,IgtMix):
            raise pm.utility.PMParamError('IgtMix instances cannot be nested. Use the + operator instead to combine mixtures.')
        else:
            raise pm.utility.PMParamError('IgtMix.insert(): Unrecognized constituent: ' + repr(sid))
        
        # Next, test for prior existence
        for this in self._c:
            if this is subst:
                raise pm.utility.PMParamError('IgtMix.insert(): substance is already in the mixture. Use the + operator to modify.')
        
        # Broadcast the quantity
        # First, force qty to be an array
        qty = np.atleast_1d(qty)
        # Handle the units
        if units is not None:
            qty = pm.units.matter(qty, subst.mw(), from_units=units, to_units='kmol')
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
            raise pm.utility.PMParamError('IgtMix.remove(): substance is not in the mixture: ' + repr(sid))
        del self._c[index]
        self._q = np.delete(self._q, index, axis=0)
    
    #
    # Type manipulation
    #
    def toigmix(self,sid):
        """TOIGMIX - convert the mixture to a static IGMIX instance
    [... igmix list ...] = m.toigmix(sid)
    
Returns a list of static igmix (low-level PYroMat data instances) built
from the IgtMix instance quantity arrays. These are less flexible, but
substantially faster than the IgtMix instances for property calculation.

sid     substance id string to use as the igmix name

If the IgtMix instance only defines one mxiture (shape == (1,)), the
sid string is used verbatim to name the new igmix instance.  If there is
more than one mixture, each is appended with "_{index}" in the same 
manner as the PYroMat collections' deduplication scheme.

It is good practice to choose an SID string that is unique (not in the
PYroMat collection), but it is not strictly required unless this is to
be saved in the permanent collection for later access.

The simpler lower-level igmix class does not allow dynamic manipulation
of mixtures, but it is substantially faster for a number of reasons:
    1) igmix instances store important intermediate calculations
    2) Unlike IgtMix instances, igmix instances call the ig and ig2 
       inner routines directly, so the top-level unit conversions and
       argument parsing is not performed redundantly.
    3) igmix instances have a faster and more flexible property 
       interface.
    4) igmix instances can be saved and added to the permanent 
       collection, while IgtMix instances cannot.

"""
        out = []
        multiple_f = len(self)>1
        for count, index in enumerate(np.ndindex(self.shape())):
            data = {'class':'igmix', 
                'doc':'Auto-generated using IgtMix.toigmix()',
                'fromfile':__name__}
            if multiple_f:
                data['id'] = sid + '_' + repr(count)
            else:
                data['id'] = sid
                
            slice_index = (slice(0,len(self._c)),) + index
            
            contents = {subst.sid():qty for subst,qty in zip(self._c, self._q[slice_index])}
            data['contents'] = contents
            data['bymass' ] = False
            out.append(pm.reg.registry['igmix'](data))
        return out
        
    def tolist(self, mode='sid'):
        """TOLIST - generate a list of the mixture contents
    mixlist = m.tolist()
        OR
    mixlist = m.tolist(mode_str)

The list entries are the substance id strings for each of the mixture's
contents.  The behavior of tolist() can be changed by setting the mode 
string to one of the following:

'sid'   (default)
Builds a list of the full substance ID strings for each substance.  For
example, ['ig.H2', 'ig.Ar'].

'hill'
Builds a list of the hill notation for each substance.  For example,
['H2', 'Ar']

'instance'
Builds a list of the ig class instances for each substance.  For example,
[<ig2, ig.H2>, <ig2, ig.O2>]

"""
        if mode == 'sid':
            return [this.sid() for this in self._c]
        elif mode == 'hill':
            return [this.hill() for this in self._c]
        elif mode == 'instance':
            return self._c.copy()
        raise pm.utility.PMParamError('IgtMix.tolist(): Unrecognized mode: ' + repr(mode))
        
        
    def todict(self, units=None):
        """TODICT - generate a dictionary representing the mixture
    mixdict = m.todict()
        OR
    mixdict = m.todict(units=unit_matter)
    
The dictionary keys are the substance id strings for each of the mixture's
contents.  The values are arrays representing the quantities of each 
substance.  By default, the quantities are reported in config['unit_matter']
but this can be overridden by setting the optional units keyword.
"""
        return {\
            subst.sid():pm.units.matter(qty, subst.mw(), from_units='kmol', to_units=units)\
            for subst,qty in self.items()}
        
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
            qty = pm.units.molar(qty, from_units='kmol', to_units=pm.config['unit_molar'])
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
        # Use broadcasting rules to convert to a mass-based quantity
        # Any will do
        mw = np.array([subst.mw() for subst in self])
        y *= mw.reshape((mw.size,) + (y.ndim-1)*(1,))
        y /= np.sum(y,axis=0)
        # If we want an array, we're done
        if asarray:
            return y
        # Convert to dict
        return {subst.sid():yy for subst,yy in zip(self, y)}
        
        
    def mass(self):
        """MASS - calculate the total mixture quantity in [unit_mass]
    mass = m.mass()
    
Returns an array with m.shape() shape.  Always reports units in the 
configured 'unit_mass' units.  Note that many users may be in the habit
of setting the 'unit_matter' setting without updating 'unit_molar' or 
'unit_mass' accordingly.  This will yield confusing results!
"""
        total = np.zeros(self.shape())
        for subst,qty in self.items():
            total += pm.units.matter(qty, subst.mw(), from_units='kmol', to_units=pm.config['unit_mass'])
        return total
        
        
    def molar(self):
        """MOLAR - calculate the total mixture quantity in [unit_molar]
    molar = m.molar()
    
Returns an array with m.shape() shape.  Always reports units in the 
configured 'unit_molar' units.  Note that many users may be in the habit
of setting the 'unit_matter' setting without updating 'unit_molar' or 
'unit_mass' accordingly.  This will yield confusing results!
"""
        n = np.sum(self._q, axis=0)
        return pm.units.molar(n, from_units='kmol', to_units=pm.config['unit_molar'])
        
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
be another IgtMix instance.  qty is the quantity array corresponding to
each substance.
"""
        return zip(self._c, self._q)

    #
    # Internal Method for Properties
    #
    def _propeval(self, prop, T, diff=False):
        """_PROPEVAL - generic evaluator for internal properties
    pval, pvalT = m._propeval(T, diff=False)
        OR
    pval, pvalT = m._propeval(T, X=X)
    
T       A numpy array of temperatures in units Kelvin.  
prop    The string name of the property method to call on members
diff    If True, the property derivative with respect to temperature is
        also evaluated.  If set to None, then the property has no 
        option for a differential output.

Returns pval, the substance property values averaged by mole fraction,
and pvalT, the substance property's derivative with respect to 
temperature.
"""
        pval = 0
        pvalT = None;
        if diff:
            pvalT=0
    
        for subst,q in zip(self._c, self._q):
            pmethod = getattr(subst, prop)
            if diff is None:
                this = pmethod(T)
            else:
                this, thisT = pmethod(T,diff=diff)
            pval += q*this
            if diff:
                pvalT += q*thisT
        return pval,pvalT
        
    def _pref(self):
        """_PREF - Calculate the reference pressure for the mixture
    pref = m._pref()
    
Returns the reference pressure in Pa.  If all the substances have the 
same reference pressure (almost always the case) then a scalar is 
returned.  Otherwise, an array with the same shape as the mixture array
is returned.
"""
        # Form an array of pressures with singleton dimensions in all
        # but the first.
        shape = (self._q.shape[0],) + (self._q.ndim-1)*(1,)
        pref = np.empty(shape, dtype=float)
        for ii,subst in enumerate(self._c):
            if subst.data['class'] == 'ig':
                pref[ii] = subst._pref_pa
            elif subst.data['class'] == 'ig2':
                pref[ii] = subst.data['pref']
            elif subst.data['class'] == 'igmix':
                subst._bootstrap()
                pref[ii] = subst._pref_pa
            else:
                raise pm.utility.PMDataError('Substance is not a recongized IG class: ' + subst.data['id'] + ', ' + subst.data['class'])
        # If all values are equal
        if np.all(pref == pref[0]):
            return pref[0].squeeze()
        # Otherwise, perform a logarithmic mean
        out = np.exp(np.sum(self._q * np.log(pref), axis=0) / np.sum(self._q, axis=0))
        return out
        

    def _iter1(self, prop, y, T, Ids, Tmin, Tmax,
                ep=1e-6, Nmax=20, verbose=False):
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
            yy,yyT = self._propeval(prop, T, diff=True)
            # note that x[Ids], yy, yyx, and all the other floating 
            # intermediates are now in m-space; the sub-set of values
            # still under iteration.
            # Calculate the error, the linear change in x, and the new x
            error[Ids] = y[Ids] - yy[Ids]
            dT[Ids] = error[Ids] / yyT[Ids]
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
    
    def mw(self, *varg, **kwarg):
        """MW - molecular weight

    mw = m.mw()
    
Returns the molecular weight in [unit_mass] / [unit_molar].
"""
        out = self._mw()
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
        R = pm.units.matter(R, self._mw(), from_units='kmol', exponent=-1)
        return R

    def Tlim(self):
        """TLIM - Calculate the temperature limits on the model
    Tlow, Thigh = m.Tlim()
    
The upper and lower temperature limits are established from the 
temperature limits of the constituent substances.  The range reported
is the widest over which all constituents report valid data.
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
    cp = m.cp(...)
    
Specific heat is calculated in units [unit_energy / unit_temperature]

All IgtMix instance properties accept flexible keyword arguments in 
pairs to specify these properties:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume
    e   internal energy     unit_energy
    h   enthalpy            unit_energy
    s   entropy             unit_energy / unit_temperature
"""
        T,p,d = self._argparse(*varg, **kwarg)
        cp = self._propeval('_cp', T, diff=None)[0]
        pm.units.energy(cp, from_units='kJ', inplace=True)
        pm.units.temperature(cp, from_units='K', inplace=True, exponent=-1)
        return cp
        
    def cv(self, *varg, **kwarg):
        """CV - Constant volume specific heat
    cv = m.cv(...)
    
Specific heat is calculated in units [unit_energy / unit_temperature] 

All IgtMix instance properties accept flexible keyword arguments in 
pairs to specify these properties:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume
    e   internal energy     unit_energy
    h   enthalpy            unit_energy
    s   entropy             unit_energy / unit_temperature
"""
        T,p,d = self._argparse(*varg, **kwarg)
        cv = self._propeval('_cp', T, diff=None)[0]
        # const_Ru is in J/mol/K == kJ/kmol/K
        cv -= pm.units.const_Ru * np.sum(self._q, axis=0)
        pm.units.energy(cv, from_units='kJ', inplace=True)
        pm.units.temperature(cv, from_units='K', inplace=True, exponent=-1)
        return cv
        
    def gam(self, *varg, **kwarg):
        """GAM - Specific heat ratio
    gam = m.gam(...)

Specific heat ratio, gamma, is dimensionless.

All IgtMix instance properties accept flexible keyword arguments in 
pairs to specify these properties:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume
    e   internal energy     unit_energy
    h   enthalpy            unit_energy
    s   entropy             unit_energy / unit_temperature        
"""
        T,p,d = self._argparse(*varg, **kwarg)
        cp = self._propeval('_cp', T, diff=None)[0] / np.sum(self._q, axis=0)
        return cp / (cp - pm.units.const_Ru)
        
    def h(self, *varg, **kwarg):
        """H - Enthalpy
    h = m.h(...)
    
Enthalpy is calculated with units [unit_energy]
    
All IgtMix instance properties accept flexible keyword arguments in 
pairs to specify these properties:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume
    e   internal energy     unit_energy
    h   enthalpy            unit_energy
    s   entropy             unit_energy / unit_temperature
"""
        T,p,d = self._argparse(*varg, **kwarg)
        h = self._propeval('_h', T)[0]
        pm.units.energy(h, from_units='kJ', inplace=True)
        return h
        
    def e(self, *varg, **kwarg):
        """E - Internal energy
    e = m.e(...)
    
Internal energy is calculated with units [unit_energy]

All IgtMix instance properties accept flexible keyword arguments in 
pairs to specify these properties:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume
    e   internal energy     unit_energy
    h   enthalpy            unit_energy
    s   entropy             unit_energy / unit_temperature
"""
        T,p,d = self._argparse(*varg, **kwarg)
        e = self._propeval('_e', T)[0]
        pm.units.energy(e, from_units='kJ', inplace=True)
        return e
        
    def s(self, *varg, **kwarg):
        """S - Entropy
    s = m.s(...)
    
Entropy is calculated with units [unit_energy / unit_temperature]

All IgtMix instance properties accept flexible keyword arguments in 
pairs to specify these properties:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume
    e   internal energy     unit_energy
    h   enthalpy            unit_energy
    s   entropy             unit_energy / unit_temperature
"""
        T,p,d = self._argparse(*varg, **kwarg)
        # We need pressure
        if p is None:
            p = d * (1000*pm.units.const_Ru * T)
        s = self._propeval('_s', T)[0] + self._sp(p=p)[0]
        pm.units.energy(s, from_units='kJ', inplace=True)
        pm.units.temperature(s, from_units='K', inplace=True, exponent=-1)
        return s
        
    def g(self, *varg, **kwarg):
        """G - Gibbs energy
    g = m.g(...)
    
Gibbs energy is calculated with units [unit_energy]

All IgtMix instance properties accept flexible keyword arguments in 
pairs to specify these properties:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume
    e   internal energy     unit_energy
    h   enthalpy            unit_energy
    s   entropy             unit_energy / unit_temperature
"""
        T,p,d = self._argparse(*varg, **kwarg)
        # We need pressure
        if p is None:
            p = d * (1000*pm.units.const_Ru * T)
        g = self._propeval('_g', T)[0] - T * self._sp(p=p)[0]
        pm.units.energy(g, from_units='kJ', inplace=True)
        return g

    def f(self, *varg, **kwarg):
        """F - Free (Helmholtz) energy
    f = m.f(...)
    
Free (Helmholtz) energy is calculated with units [unit_energy]

All IgtMix instance properties accept flexible keyword arguments in 
pairs to specify these properties:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume
    e   internal energy     unit_energy
    h   enthalpy            unit_energy
    s   entropy             unit_energy / unit_temperature
"""
        T,p,d = self._argparse(*varg, **kwarg)
        # We need pressure
        if p is None:
            p = d * (1000*pm.units.const_Ru * T)
        f = self._propeval('_f', T)[0] - T*self._sp(p=p)[0]
        pm.units.energy(f, from_units='kJ', inplace=True)
        return f

    def T(self, *varg, **kwarg):
        """T - Temperature
        
    T = m.T(...)
    
Temperature is calculated with units [unit_temperature]

All IgtMix instance properties accept flexible keyword arguments in 
pairs to specify these properties:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume
    e   internal energy     unit_energy
    h   enthalpy            unit_energy
    s   entropy             unit_energy / unit_temperature
"""
        T,p,d = self._argparse(*varg, **kwarg)
        if T.base is not None:
            T = np.copy(T)
        pm.units.temperature_scale(T, from_units='K', inplace=True)
        return T
    
    def p(self, *varg, **kwarg):
        """P - Pressure
    p = m.p(...)
    
Pressure is calculated with units [unit_pressure]

All IgtMix instance properties accept flexible keyword arguments in 
pairs to specify these properties:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume
    e   internal energy     unit_energy
    h   enthalpy            unit_energy
    s   entropy             unit_energy / unit_temperature
"""
        T,p,d = self._argparse(*varg, **kwarg)
        if p is None:
            p = d * (1000*pm.units.const_Ru * T)
        # Test to see if this is a view
        if p.base is not None:
            p = np.copy(p)
        pm.units.pressure(p, from_units='Pa', inplace=True)
        return p

    def d(self, *varg, **kwarg):
        """D - Density
    d = m.d(...)
    
Density is calculated with units [unit_matter / unit_volume]

All IgtMix instance properties accept flexible keyword arguments in 
pairs to specify these properties:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume
    e   internal energy     unit_energy
    h   enthalpy            unit_energy
    s   entropy             unit_energy / unit_temperature
"""
        T,p,d = self._argparse(*varg, **kwarg)
        if d is None:
            d = p / (1000*pm.units.const_Ru * T)
        # Test to see if this is a view
        if d.base is not None:
            d = np.copy(d)
        scale = pm.units.matter(1., self._mw(), from_units='kmol')
        scale = pm.units.volume(scale, from_units='m3', exponent=-1)
        np.multiply(d, scale, out=d)
        return d
        
        
    def v(self, *varg, **kwarg):
        """V - Specific Volume
    v = m.v(...)
    
Specific volume is calculated with units [unit_volume / unit_matter]

All IgtMix instance properties accept flexible keyword arguments in 
pairs to specify these properties:
    T   temperature         unit_temperature
    p   pressure            unit_pressure
    d   density             unit_matter / unit_volume
    v   specific volume     unit_volume
    e   internal energy     unit_energy
    h   enthalpy            unit_energy
    s   entropy             unit_energy / unit_temperature
"""
        return 1. / self.d(*varg, **kwarg)
        
        

def fromigmix(source):
    """FROMIGMIX - create an IgtMix instance from and igmix instance
    igtmix = fromigmix(igm)
    
Creates an IgtMix instance from a low-level igmix instance.  The IgtMix
instance will automatically contain a single unit of matter (with the 
currently configured unit matter).  
"""
    if pm.units.ismass():
        return IgtMix(source.Y())
    else:
        return IgtMix(source.X())

