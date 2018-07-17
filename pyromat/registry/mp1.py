# MP1
#   PYroMat Multi-phase generalist class
#   Calculates physical properties from a fit for the helmholtz free 
#   energy in terms of density and temperature.

import numpy as np
import pyromat as pm


def _f1_x(x,param,diff=2):
    """f1 pass-thru
    x
"""
    if diff>0:
        return x, 1., 0.
    return x, 0., 0.

def _f1_lin(x,param,diff=2):
    """f1 linear
    param[0] + param[1]*x
"""
    f = param[0] + param[1]*x
    if diff>0:
        return f, param[1], 0.
    return f, 0., 0.

def _f1_pow(x,param,diff=2):
    """f1 power
    x**param
"""
    f = x**param
    if diff>0:
        fx = param*f/x
    else:
        fx = 0.
    if diff>1:
        fxx = (param-1.)*fx/x
    else:
        fxx = 0.
    return f,fx,fxx
    
def _f1_inv(x,param,diff=2):
    """f1 inverse
    1/x
"""
    f = 1./x
    if diff>0:
        fx = -f/x
    else:
        fx = 0.
    if diff>1:
        fxx = -2*fx/x
    else:
        fxx = 0.
    return f,fx,fxx

def _f1_exp(x,param,diff=2):
    """f1 exponent
    exp(param*x)
"""
    f = np.exp(x*param)
    if diff>0:
        fx = param*f
    else:
        fx = 0.
    if diff>1:
        fxx = param*fx
    else:
        fxx = 0.
    return f,fx,fxx

def _f1_log(x,param,diff=2):
    """f1 natural log
    log(x)
"""
    f = np.log(x)
    if diff>0:
        fx = 1./x
    else:
        fx = 0.
    if diff>1:
        fxx = -fx/x
    else:
        fxx = 0.
    return f,fx,fxx

def _f1_poly(x,coef,diff=2):
    """Polynomial evaluation
(p, px, pxx) = _f1_poly(x,coef,diff=2)

Evaluates a polynomial on x and its derivatives
x       x value
coef    the coefficient list
diff    the highest order derivative to evaluate (0,1, or 2)

Returns
p       polynomial value at p(x)
px      dp/dx
pxx     d2p/dx2

Each element of coef is a two-element list defining a term in the 
polynomial; the exponent on x and the corresponding coefficient.  It 
must be sorted in descending order by the first column.

The list,
[[4, 0.1], [2, 0.2], [1, 1.2], [0, 0.5]]
corrsponds to the polynomial
p(x) = .5 + 1.2*x + 0.2*x**2 + 0.1*x**4

This approach assumes that most polynomials used will be "sparse";
that most of the coefficients will be zero, so they need not be 
stored.
"""
    # initialize the final polynomial and its derivatives
    p = 0.  # total polynomial
    px = 0.
    pxx = 0.
    # From here, we loop over terms of the form a*(x**ii)
    # If a particular ii is not found in the data, then
    # its coefficient is treated as zero.
    # What is the largest ii?
    II = coef[0][0]
    
    # On which coefficient are we currently operating?
    index = 0
    # This is a flag that indicates the active index was used in 
    # the last loop, so it needs to be incremented.
    
    for ii in range(II,-1,-1):
        if index<len(coef) and coef[index][0] == ii:
            # Fold the coefficient values into the p expansion
            if diff > 1:
                pxx = 2*px + x * pxx
            if diff > 0:
                px = p + x * px
            p = coef[index][1] + x * p
            # Move on to the next coefficient
            index += 1
        # If the current x exponent is not represented, execute a
        # p-expansion with zero q.
        else:
            if diff > 1:
                pxx = 2*px + x * pxx
            if diff > 0:
                px = p + x * px
            p = x * p
            
    return p,px,pxx
    
    
def _f2_x(x,y,param,diff=2):
    f,fx,fxx = _f1_x(x,param,diff)
    return f,fx,0.,fxx,0.,0.
    
def _f2_y(x,y,param,diff=2):
    f,fy,fyy = _f1_x(y,param,diff)
    return f,0.,fy,0.,0.,fyy

def _f2_lin(x,y,param,diff=2):
    """f2 linear
    param[0] + param[1]*x + param[2]*y
"""
    f = param[0] + param[1]*x + param[2]*y
    if diff>0:
        return f, param[1], param[2], 0., 0., 0.
    return f, 0., 0., 0., 0., 0.
    
def _f2_linx(x,y,param,diff=2):
    f,fx,fxx = _f1_lin(x,param,diff)
    return f, fx, 0., fxx, 0., 0.

def _f2_liny(x,y,param,diff=2):
    f,fy,fyy = _f1_lin(y,param,diff)
    return f, 0., fy, 0., 0., fyy

def _f2_powx(x,y,param,diff=2):
    f,fx,fxx = _f1_pow(x,param,diff)
    return f, fx, 0., fxx, 0., 0.

def _f2_powy(x,y,param,diff=2):
    f,fy,fyy = _f1_pow(y,param,diff)
    return f, 0., fy, 0., 0., fyy
    
def _f2_invx(x,y,param,diff=2):
    f,fx,fxx = _f1_inv(x,param,diff)
    return f, fx, 0., fxx, 0., 0.

def _f2_invy(x,y,param,diff=2):
    f,fy,fyy = _f1_inv(y,param,diff)
    return f, 0., fy, 0., 0., fyy
    
def _f2_expx(x,y,param,diff=2):
    f,fx,fxx = _f1_exp(x,param,diff)
    return f, fx, 0., fxx, 0., 0.

def _f2_expy(x,y,param,diff=2):
    f,fy,fyy = _f1_exp(y,param,diff)
    return f, 0., fy, 0., 0., fyy

def _f2_logx(x,y,param,diff=2):
    f,fx,fxx = _f1_log(x,param,diff)
    return f, fx, 0., fxx, 0., 0.

def _f2_logy(x,y,param,diff=2):
    f,fy,fyy = _f1_log(y,param,diff)
    return f, 0., fy, 0., 0., fyy

def _f2_poly(x,y,coef,diff=2):
    """Polynomial evaluation
(p, px, py, pxx, pxy, pyy) = _f2_poly(x,y,coef,diff=2)

Evaluates a polynomial on x and y and its derivatives.
x       x value
y       y value
coef    coefficient list
diff    the highest order derivative to evaluate (0,1, or 2)

Returns
p       polynomial value at p(x,y)
px      dp/dx
py      dp/dy
pxx     d2p/dx2
pxy     d2p/dxdy
pyy     d2p/dy2

Each element of coef is a three-element list defining a term in the 
polynomial; the x-exponent, the y-exponent, and the corresponding
coefficient.  It must be sorted in descending order by the first column
and then the second column.

The list,
[[1, 1, 0.1], [0, 2, 0.2], [0, 1, 1.2], [0, 0, 0.5]]
corrsponds to the polynomial
p(x,y) = .5 + 1.2y + .2y**2 + 0.1xy

This approach assumes that most polynomials used will be "sparse";
that most of the coefficients will be zero, so they need not be 
stored.
"""
    
    # initialize the final polynomial and its derivatives
    p = 0.  # total polynomial
    px = 0.
    py = 0.
    pxx = 0.
    pxy = 0.
    pyy = 0.
    # From here, we loop over terms of the form a*(x**ii)*(y**jj)
    # If a particular ii,jj combination is not found in the data, then
    # its coefficient is treated as zero.
    # What is the largest ii?
    II = coef[0][0]
    
    # On which coefficient are we currently operating?
    index = 0
    # This is a flag that indicates the active index was used in 
    # the last loop, so it needs to be incremented.
    
    for ii in range(II,-1,-1):
        # If the current x-exponent is the same one represented in
        # the active coefficient row, then calculate q.
        if index<len(coef) and coef[index][0] == ii:
            # For this value of ii, what is the largest jj?
            JJ = coef[index][1]
            # q is a sub-polynomial on y that represents the 
            # variation on y of all terms that share the same
            # power in x.  This inner loop is much like the loop
            # on ii, except that it looks at both the x and y 
            # exponents.
            q = 0
            qy = 0
            qyy = 0
            for jj in range(JJ,-1,-1):
                if diff > 1:
                    qyy = 2*qy + y*qyy
                if diff > 0:
                    qy = q + y*qy
                # If the current y-exponent is represented in the 
                # active coefficient row, then fold it into the q
                # expansion.
                if index<len(coef) and coef[index][0] == ii and coef[index][1] == jj:
                    q = coef[index][2] + y*q
                    # increment the active index
                    index += 1
                else:
                    q *= y
                
            # Fold the current q values into the p expansion
            # Update the highest derivatives first since they depend
            # on the historical values of the lower derivatives
            if diff > 1:
                pyy = qyy + x * pyy
                pxx = 2*px + x * pxx
                pxy = py + x * pxy
            if diff > 0:
                px = p + x * px
                py = qy + x * py
            p = q + x * p
        # If the current x exponent is not represented, execute a
        # p-expansion with zero q.
        else:
            if diff > 1:
                pyy = x * pyy
                pxx = 2*px + x * pxx
                pxy = py + x * pxy
            if diff > 0:
                px = p + x * px
                py = x * py
            p = x * p
            
    return p,px,py,pxx,pxy,pyy
    

def _f2_spoly(x,y,coef,diff):
    "Sparse polynomial"
    p = 0.
    px = 0.
    py = 0.
    pxx = 0.
    pxy = 0.
    pyy = 0.
    for row in coef:
        t = x**row[0] * y**row[1] * row[2]
        if diff>0:
            tx = t * row[0] / x
            ty = t * row[1] / y
        else:
            tx = 0.
            ty = 0.
        if diff>1:
            txx = tx * (row[0]-1) / x
            tyy = ty * (row[1]-1) / y
            txy = tx * row[1] / y
        else:
            txx = 0.
            tyy = 0.
            txy = 0.
        p += t
        if diff>0:
            px += tx
            py += ty
        if diff>1:
            pxx += txx
            pxy += txy
            pyy += tyy
    return p,px,py,pxx,pxy,pyy
    
    
_f1_dict = {
    'x':_f1_x,
    'lin':_f1_lin,
    'pow':_f1_pow,
    'inv':_f1_inv,
    'exp':_f1_exp,
    'log':_f1_log,
    'poly':_f1_poly
}


def _f1eval(x, fdict, diff=2):
    """Evaluates a group funciton element recursively
f, fx, fxx = _f1eval(x, fdict, diff=2)

fdict is the dictionary that defines the function type and behavior.  
The dictionary must contain a 'type' key, which identifies the function
to be evaluated.

There are three keys that define the function's behavior:
    KEY     DESCRIPTION
    type    The name of the function to be evaluated
    param   Parameters to define the behavior of the function
    arg     specifies a group or function for recursion

Each function uses the 'param' key differently to define its behavior.  
The table below lists the function types, and shows how the param key 
is used.
    TYPE        FORM
    x           x
    lin         param[0] + param[1]*x
    poly        _p1eval(x, param)
    pow         x**param
    inv         1/x
    exp         exp(param*x)
    log         log(x)
    
If the 'arg' key is defined, then its value will be treated as a group
or function used to evaluate the argument to the present function.

For example,
    p( x**0.5 )
    
Might be evaluated by
{
    'type':'poly', 
    'param':[...coef...], 
    'arg':{
        'type':'pow', 
        'param':0.5,
    }
}
"""
    # If there is nesting
    if 'arg' in fdict:
        if isinstance(fdict['arg'],dict):
            x_0,x_1,x_2 = _f1eval(x, fdict['arg'], diff=diff)
        else:
            x_0,x_1,x_2 = _g1eval(x, fdict['arg'], diff=diff)
    else:
        x_0 = x
        x_1 = 1.
        x_2 = 0.
    
    f,fx,fxx = _f1_dict[fdict['type']](x_0, fdict['param'], diff)
    if diff>1:
        fxx = fxx*x_1*x_1 + fx*x_2
    if diff>0:
        fx *= x_1
        
    return f,fx,fxx



_f2_dict = {
    'x':_f2_x,
    'y':_f2_y,
    'lin':_f2_lin,
    'linx':_f2_linx,
    'liny':_f2_liny,
    'powx':_f2_powx,
    'powy':_f2_powy,
    'invx':_f2_invx,
    'invy':_f2_invy,
    'expx':_f2_expx,
    'expy':_f2_expy,
    'logx':_f2_logx,
    'logy':_f2_logy,
    'poly':_f2_poly,
    'spoly':_f2_spoly
}

def _f2eval(x, y, fdict, diff=2):
    """Evaluates a group funciton element recursively
f, fx, fy, fxx, fxy, fyy = _f2eval(x, y, fdict, diff=2)

fdict is the dictionary that defines the function type and behavior.  
The dictionary must contain a 'type' key, which identifies the function
to be evaluated.

There are three keys that define the function's behavior:
    KEY     DESCRIPTION
    type    The name of the function to be evaluated
    param   Parameters to define the behavior of the function
    argx    specifies a group or function for recursion on x only
    argy    specifies a group or function for recursino on y only

Each function uses the 'param' key differently to define its behavior.  
The table below lists the function types, and shows how the param key 
is used.
    TYPE        FORM
    x           x
    y           y
    poly        _p1eval(x, y, param)
    invx        1/x
    invy        1/y
    linx        param[0] + param[1]*x
    liny        param[0] + param[1]*y
    lin         param[0] + param[1]*x + param[2]*y
    powx        x**param
    powy        y**param
    expx        exp(param*x)
    expy        exp(param*y)
    logx        log(x)
    logy        log(y)
    
If the 'arg' key is defined, then its value will be treated as a group
or function used to evaluate the argument to the present function.

For example,
    p( x**0.5, y )
    
Might be evaluated by
{
    'type':'poly', 
    'param':[...coef...], 
    'argx':{
        'type':'pow', 
        'param':0.5,
    }
}
"""
    f = 0.
    fx = 0.
    fy = 0.
    fxx = 0.
    fxy = 0.
    fyy = 0.

    # If there is nesting
    if 'argx' in fdict:
        if isinstance(fdict['argx'],dict):
            x_0,x_1,x_2 = _f1eval(x, fdict['argx'], diff)
        else:
            x_0,x_1,x_2 = _g1eval(x, fdict['argx'], diff)
    else:
        x_0 = x
        x_1 = 1.
        x_2 = 0.
        
    if 'argy' in fdict:
        if isinstance(fdict['argy'],dict):
            y_0,y_1,y_2 = _f1eval(y, fdict['argy'], diff)
        else:
            y_0,y_1,y_2 = _g1eval(y, fdict['argy'], diff)
    else:
        y_0 = y
        y_1 = 1.
        y_2 = 0.

    f,fx,fy,fxx,fxy,fyy = _f2_dict[fdict['type']](x_0, y_0, fdict['param'], diff)
    if diff>1:
        fxx = fxx*x_1*x_1 + fx*x_2
        fyy = fyy*y_1*y_1 + fy*y_2
        fxy = fxy*x_1*y_1
    if diff>0:
        fx = fx*x_1
        fy = fy*y_1
    return f,fx,fy,fxx,fxy,fyy


def _g1eval(x, group, diff=2):
    """Evaluates a group of terms and their derivatives
g, gx, gxx = _g1eval(x, group, diff=2)

Returns
g       The value of the group at x
gx      The derivative of the group with respect to x
gxx     The second derivative of the group with respect to x
Arguments
x       The floating point value for x
group   The group to be evaluated at x (see below)
diff    The highest derivative to evaluate (0,1, or 2)

Groups are expansions of the form 
g = sum_k  f1(x) * f2(x) * ...
Each group is a sum of individual terms, each of which is a 
multiplication of standard functions on x.  Each function, f1, f2....
must be selected from a bank of familiar functions, and x may be 
modified before it is passed to the individual functions.

Groups are lists of terms.  Each term is a list of individual functions.
Therefore, groups are lists of lists.  Each element of the innermost 
list must either be a dictionary that is evaluated by _f1eval() or 
another nested list that can be evaluated recursively by _g1eval().

For example, the group
    g(x) = exp(-x) * (0.5*x**2 - 0.1*x) + sqrt(x)
    
might be represented by the group
    group = [
        [
            {'type':'exp', 'param':-1},
            {'type':'poly': 'coef':[[2,0.5], [1,-0.1]]}
        ],
        [
            {'type':'pow', 'param':0.5}
        ]
    ]
"""
    g = 0.
    gx = 0.
    gxx = 0.

    # Evaluate the terms one-by-one
    for term in group:
        # Let t be the term value and its derivatives
        t = 1.
        tx = 0.
        txx = 0.
        
        for function in term:
            if isinstance(function, dict):
                f,fx,fxx = _f1eval(x, function, diff)
                if diff>1:
                    txx = txx*f + 2*tx*fx + t*fxx
                if diff>0:
                    tx = f*tx + fx*t
                t *= f
            elif isinstance(function, list):
                f,fx,fxx = _g1eval(x, function, diff)
                if diff>1:
                    txx = txx*f + 2*tx*fx + t*fxx
                if diff>0:
                    tx = f*tx + fx*t
                t *= f
            else:
                t *= function
                if diff>0:
                    tx*=function
                if diff>1:
                    txx*=function
    
        # Fold in the last term to the group
        g += t
        gx += tx
        gxx += txx
    
    # All done!
    return g,gx,gxx




def _g2eval(x, y, group, diff=2):
    """Evaluates a group of terms and their derivatives
g, gx, gy, gxx, gxy, gyy = _g2eval(x, y, group, diff=2)

Returns
g       The value of the group at x
gx      The derivative of the group with respect to x
gy
gxx     The second derivative of the group with respect to x
gxy
gyy
Arguments
x       The floating point value for x
y
group   The group to be evaluated at x (see below)
diff    The highest derivative to evaluate (0,1, or 2)

Groups are expansions of the form 
g = sum  f1(x,y) * f2(x,y) * ...
Each group is a sum of individual terms, each of which is a 
multiplication of standard functions on x.  Each function, f1, f2....
must be selected from a bank of familiar functions, and x may be 
modified before it is passed to the individual functions.

Groups are lists of terms.  Each term is a list of individual functions.
Therefore, groups are lists of lists.  Each element of the innermost 
list must either be a dictionary that is evaluated by _f2eval() or 
another nested list that can be evaluated recursively by _g2eval().

For example, the group
    g(x) = exp(-x) * (0.5*x**2 - 0.1*x) + 2.5*sqrt(x)
    
might be represented by the group
    group = [
        [
            {'type':'exp', 'param':-1},
            {'type':'poly': 'coef':[[2,0.5], [1,-0.1]]}
        ],
        [
            2.5,
            {'type':'pow', 'param':0.5}
        ]
    ]
"""
    g = 0.
    gx = 0.
    gy = 0.
    gxx = 0.
    gxy = 0.
    gyy = 0.

    # Evaluate the terms one-by-one
    for term in group:
        # Let t be the term value and its derivatives
        t = 1.
        tx = 0.
        ty = 0.
        txx = 0.
        txy = 0.
        tyy = 0.
        
        for function in term:
            if isinstance(function, dict):
                f,fx,fy,fxx,fxy,fyy = _f2eval(x, y, function, diff)
                if diff>1:
                    txx = txx*f + 2*tx*fx + t*fxx
                    tyy = tyy*f + 2*ty*fy + t*fyy
                    txy = txy*f + tx*fy + ty*fx + t*fxy
                if diff>0:
                    tx = f*tx + fx*t
                    ty = f*ty + fy*t
                t *= f
            elif isinstance(function, list):
                f,fx,fy,fxx,fxy,fyy = _g2eval(x, y, function, diff)
                if diff>1:
                    txx = txx*f + 2*tx*fx + t*fxx
                    tyy = tyy*f + 2*ty*fy + t*fyy
                    txy = txy*f + tx*fy + ty*fx + t*fxy
                if diff>0:
                    tx = f*tx + fx*t
                    ty = f*ty + fy*t
                t *= f
            else:
                t *= function
                if diff>0:
                    tx*=function
                    ty*=function
                if diff>1:
                    txx*=function
                    tyy*=function
                    txy*=function
            #print "  -->",f,fx,fy,fxx,fxy,fyy
        #print "==>",t,tx,ty,txx,txy,tyy
    
        # Fold in the last term to the group
        g += t
        gx += tx
        gy += ty
        gxx += txx
        gxy += txy
        gyy += tyy
    
    # All done!
    return g,gx,gy,gxx,gxy,gyy


class mp1(pm.reg.__basedata__):
    """The PYroMat multi-phase generalist class 1
    
Models thermo-physical properties of a liquid-gas system using a general
fit for helmholtz free energy.  These "Spand & Wagner" fits are 
evaluated in a polynomial form with exponential post factors.

The primary property methods calculate their respective properties from
temperature and pressure, and appear with no subscripts.  However, these
may not always be the most convenient or the fastest algorithms.  There
are a number of "secondary" methods that provide the same properties as
functions of other parameters.

The MP1 class is divided into three layers of methods (routines).  

--- USER ROUTINES ---
Accept data in any format (array or scalar) and in whatever units are
configured in the PYroMat configuration object.  These routines are
responsible for handling the necessary unit conversions and arranging
all the arguments in appropriately typed numpy arrays.

Values from these methods are returned in appropriately broadcast arrays
in the correctly configured units.

--- INNER ROUTINES ---
These methods presume that all arguments are numpy arrays and that they
are in a common unit system.  This prevents redundant calls to the unit
conversion functions as MP1 methods call one another.
    Energy:     J
    Mass:       kg
    Moles:      kmol
    Pressure:   Pa
    Temperature:K

--- PRIMATIVE ROUTINES ---
Methods that have been labeled as primative routines should UNDER NO
CIRCUMSTANCES be called by the user.  They accept non-dimensionalized
arguments and return non-dimensional parameters.  These are encapsulated
as independent methods either because they are complicated and need to 
be called by a number of other methods, or because separating them made
sense for numerical efficiency.  In summary: these aren't the methods
you're looking for.

--- DATA DICTIONARY ---
The MP1 data dictionary must have certain data "groups" to define the 
various empirical fits.  Each group is a dictionary (within the 
dictionary) that defines the various parameters necessary for at least
one of the inner methods.

PSgroup         Saturated pressure data group
    Tscale      Temperature scale for normalizing T in the fit
    pscale      Pressure scale for re-scaling the result
    coef        a coefficient group to be passed to _poly1()
If tt = T/Tscale
    tt*log(ps/pscale) = p(1-tt)
where p() is the polynomial defined by the coef list.

DSLgroup        Saturated liquid density data group
    Tscale      Temperature scale for normalizing T in the fit
    dscale      Density scale for re-scaling the result
    coef        a coefficient group to be passed to _poly1()
If tt = T/Tscale
    dsl / dscale = p(1-tt)
where p() is the polynomial defined by the coef list.

DSVgroup        Saturated vapor density data group
    Tscale      Temperature scale for normalizing T in the fit
    dscale      Density scale for re-scaling the result
    coef        a coefficient group to be passed to _poly1()
If tt = T/Tscale
    log(dsv/dscale) = p(1-tt)
where p() is the polynomial defined by the coef list.

AOgroup         Helmholtz free energy ideal gas group
    Tscale      Temperature scale for normalizing T
    dscale      density scale for normalizing d
    logt        a scalar coefficient of a log(tt) term
    coef        a coefficient list to be passed to _poly1()
If tt = Tscale/T    <=== INVERSE!
and dd = d/dscale
    ao = log(d) + LOGT*log(tt) + p(tt)
    Ao = ao * R * T
where LOGT is the coefficient defined by the 'logt' parameter, and p is
the polynomial defined by the coef list

ARgroup         Helmholtz free energy residual group
    Tscale      Temperature scale for normalizing T
    dscale      density scale for normalizing d
    coef        a nested list of coefficient lists
Each element of coef is, itself a coefficient list intended to be passed
to _poly2().  After the first element, each individual polynomial is 
multiplied by exp(-dd**k) where k is the index in the coef list.

If tt = Tscale/T    <=== INVERSE!
and dd = d/dscale
    ar = p0(tt,dd) + exp(-dd)*p1(tt,dd) + exp(-dd**2)*p2(tt,dd) + ...
    Ar = ar * R * T

Additionally, there are a number of parameters that define static 
properties

Tlim            A two-element list of the upper and lower temperatures
                for which the data set is valid.
plim            A two-element list of the upper and lower pressures for
                which the data set is valid.
dlim            A two-element list the represent practical maximum and
                minimum densities over the entire data set.  These are 
                NOT guaranteed limits of validity.
Tc, pc, dc      Critical temperature, pressure, and density
Tt, pt          Triple-point temperature and pressure
R               Ideal gas constant 8.314 / mw
mw              Molecular weight
contents        A dictionary with a key for each atom and a value for 
                its count in the molecule.  For example, CO2 would 
                have content = {'C':1, 'O':2}
                
There are also the typical mandatory PYroMat meta data elements:
id              What substance is this?
doc             Where did it come from?
class           What class should be used to evaluate the data?
"""

    def _poly2(self,x,y,group,diff=2):    
        """Polynomial evaluation
(p, px, py, pxx, pxy, pyy) = _poly(x,y,coef,diff=2)

Evaluates a polynomial on x and y and its derivatives.
x       x value
y       y value
param   coefficient list
diff    the highest order derivative to evaluate (0,1, or 2)

Returns
p       polynomial value at p(x,y)
px      dp/dx
py      dp/dy
pxx     d2p/dx2
pxy     d2p/dxdy
pyy     d2p/dy2

The behavior of poly2 is very much the same as poly1, but for funcitons
of two variables.  The pre- and post- exponents for poly2 are expected
to be lists or tuples: [prex, prey], [postx, posty]
The innermost lists defining the terms must contain three elements:
[powxN, powyN, coefN].  The coefficients must be sorted by x-power and
then by y-power in descending order.

A coefficient list might appear

coef = [
    [
        [prex, prey], 
        [postx, posty],
        [powxN, powyN, coefN],
        ...
        [powx0, powy0, coef0]
    ],
    [
        [prex, prey], 
        [postx, posty],
        [powxN, powyN, coefN],
        ...
        [powx0, powy0, coef0]
    ]
]
    
The pre-exponents are applied to the arguments to the polynomial, and
the post-exponents are applied after the polynomial is evaluated, so 
that the value returned is
    x**postx * y**posty * p(x**prex, y**prey)

Starting with the third element (element 2) of the coefficient list,
each element of coef is a three-element list defining a term in the 
polynomial; the x-exponent, the y-exponent, and the corresponding
coefficient.  It must be sorted in descending order by the first column
and then the second column.

For example, the list,
[[[1,1], [0,0], [1, 1, 0.1], [0, 2, 0.2], [0, 1, 1.2], [0, 0, 0.5]]]
corresponds to the polynomial
p(x,y) = .5 + 1.2y + .2y**2 + 0.1xy

Efficient polynomial evaluation algorithms are normally restricted to
positive integer exponents, but many thermodynamic property models use 
much more interesting polynomials.  The pre- and post- exponents can be
used to acheive a much wider range of functions.

For example,
    p(x,y) = x**(-1.5) + x**(3.5)
might be expressed as a coefficient list
    [[[0.5,1], [-1.5, 0], [10, 0, 1], [0, 0, 1]]]
The pre- and post- exponents make this equivalent to
    xx = x**0.5
    p(x,y) = x**(-1.5) (xx**10 + 1)
which is equivalent to the original polynomial, except that the core of
the evaluation algorithm only operates on positive integers.
"""
    
        g = 0.  # total group
        gx = 0.
        gy = 0.
        gxx = 0.
        gxy = 0.
        gyy = 0.

        for coef in group:
            # initialize the final polynomial and its derivatives
            p = 0.  # total polynomial
            px = 0.
            py = 0.
            pxx = 0.
            pxy = 0.
            pyy = 0.
            
            # collect the pre- and post-exponentials
            prex,prey = coef[0]
            postx,posty = coef[1]
            
            # Apply the pre-exponentials
            if prex!=1.:
                x_0 = x**prex
                if diff>0:
                    x_1 = x_0*prex/x
                else:
                    x_1 = 0.
                if diff>1:
                    x_2 = x_1*(prex-1.)/x
                else:
                    x_2 = 0.
            else:
                x_0 = x
                x_1 = 1.
                x_2 = 0.
                
            if prey!=1.:
                y_0 = y**prey
                if diff>0:
                    y_1 = y_0*prey/y
                else:
                    y_1 = 0.
                if diff>1:
                    y_2 = y_1*(prey-1.)/y
                else:
                    y_2 = 0.
            else:
                y_0 = y
                y_1 = 1.
                y_2 = 0.

            # From here, we loop over terms of the form a*(x**ii)*(y**jj)
            # If a particular ii,jj combination is not found in the data, then
            # its coefficient is treated as zero.
            # What is the largest ii?
            II = coef[2][0]
            
            # On which coefficient are we currently operating?
            index = 2
            # This is a flag that indicates the active index was used in 
            # the last loop, so it needs to be incremented.
            
            for ii in range(II,-1,-1):
                # If the current x-exponent is the same one represented in
                # the active coefficient row, then calculate q.
                if index<len(coef) and coef[index][0] == ii:
                    # For this value of ii, what is the largest jj?
                    JJ = coef[index][1]
                    # q is a sub-polynomial on y that represents the 
                    # variation on y of all terms that share the same
                    # power in x.  This inner loop is much like the loop
                    # on ii, except that it looks at both the x and y 
                    # exponents.
                    q = 0
                    qy = 0
                    qyy = 0
                    for jj in range(JJ,-1,-1):
                        if diff > 1:
                            qyy = 2*qy + y_0*qyy
                        if diff > 0:
                            qy = q + y_0*qy
                        # If the current y-exponent is represented in the 
                        # active coefficient row, then fold it into the q
                        # expansion.
                        if index<len(coef) and coef[index][0] == ii and coef[index][1] == jj:
                            q = coef[index][2] + y_0*q
                            # increment the active index
                            index += 1
                        else:
                            q *= y_0
                        
                    # Fold the current q values into the p expansion
                    # Update the highest derivatives first since they depend
                    # on the historical values of the lower derivatives
                    if diff > 1:
                        pyy = qyy + x_0 * pyy
                        pxx = 2*px + x_0 * pxx
                        pxy = py + x_0 * pxy
                    if diff > 0:
                        px = p + x_0 * px
                        py = qy + x_0 * py
                    p = q + x_0 * p
                # If the current x exponent is not represented, execute a
                # p-expansion with zero q.
                else:
                    if diff > 0:
                        if diff > 1:
                            pyy = x_0 * pyy
                            pxx = 2*px + x_0 * pxx
                            pxy = py + x_0 * pxy
                        px = p + x_0 * px
                        py = x_0 * py
                    p = x_0 * p
                    
            # Modify the derivatives for the pre-exponnetials
            if prex!=1.:
                if diff>0:
                    if diff>1:
                        pxx = pxx*x_1*x_1 + px*x_2
                        pyy = pyy*y_1*y_1 + py*y_2
                        pxy = pxy*x_1*y_1
                    px *= x_1
                    py *= y_1
                
            # Apply the post-exponentials
            if postx!=0:
                f = x**postx
                if diff>0:
                    fx = postx*f/x
                    if diff>1:
                        fxx = fx*(postx-1)/x
                        pxx = pxx*f + 2.*px*fx + p*fxx
                        pyy = pyy*f
                        pxy = pxy*f + py*fx
                    px = px*f + p*fx
                    py = py*f
                p *= f
            if posty!=0:
                f = y**posty
                if diff>0:
                    fy = posty*f/y
                    if diff>1:
                        fyy = fy*(posty-1)/y
                        pyy = pyy*f + 2.*py*fy + p*fyy
                        pxx = pxx*f
                        pxy = pxy*f + px*fy
                    py = py*f + p*fy
                    px = px*f
                p *= f
            
            # If the group has only one coefficient set, just return
            if len(group) == 1:
                return p,px,py,pxx,pxy,pyy
                
            g += p
            if diff>0:
                gx += px
                gy += py
                if diff>1:
                    gxx += pxx
                    gxy += pxy
                    gyy += pyy
                    
        return g,gx,gy,gxx,gxy,gyy


    def _poly1(self,x,group,diff=2):    
        """Polynomial evaluation
(p, px, pxx) = _poly1(x,coef,diff=2)

Evaluates a polynomial on x and y and its derivatives.
x       x value
coef    coefficient list
diff    the highest order derivative to evaluate (0,1, or 2)

Returns
p       polynomial value at p(x)
px      dp/dx
pxx     d2p/dx2

When diff is less than 2, the corresponding values of px and pxx are
returned as 0.  The default for diff is 2 to protect against careless
treatment as if these values ARE zero, but reducing diff will make poly1
execute more efficiently.

The coefficient list represents a nested list structure defining a poly-
nomial series.  The inner-most lists contain two elements, specifying
an integer exponent and a coefficient: [power, coefficient] for each
term.  They are contained in a list that represents groups of terms.
The groups must be sorted by power from highest to lowest.

Each group is lead by two elements that define a pre- and post- 
exponents.
    [pre, post, [powN, coefN], ... , [pow0, coef0]]

This defines a polynomial of the form
    x**post * p(x**pre)
    
The powers must be integers, but no such restriciton exists on the pre-
and post- exponents.  This permits efficient evaluaiton of polynomials
with rational exponents.

The highest level list contains a list of these groups, so that separate
pre- and post- exponentials may be applied to certain terms of the 
polynomial.

coef = [
    [
        pre, 
        post,
        [
            [powN, coefN],
            ...,
            [pow0, coef0]
        ]
    ],
    [
        pre, 
        post,
        [
            [powN, coefN],
            ...,
            [pow0, coef0]
        ]
    ]
]

In a simple example, the polynomial,
    p(x) = 2*x**-1.5 - x**0.5
    
might be specified
[[  0.5, -1.5, [0, 2.], [4, -1.]]]
"""
        g = 0.
        gx = 0.
        gxx = 0.
    
        for coef in group:
            # initialize the final polynomial and its derivatives
            p = 0.  # total polynomial
            px = 0.
            pxx = 0.

            # collect the pre- and post-exponentials
            pre = coef[0]
            post = coef[1]
            
            # Apply the pre-exponentials
            if pre!=1.:
                x_0 = x**pre
                if diff>0:
                    x_1 = x_0*pre/x
                else:
                    x_1 = 0.
                if diff>1:
                    x_2 = x_1*(pre-1.)/x
                else:
                    x_2 = 0.
            else:
                x_0 = x
                x_1 = 1.
                x_2 = 0.

            # From here, we loop over terms of the form a*(x**ii)
            # If a particular ii,jj combination is not found in the data, then
            # its coefficient is treated as zero.
            # What is the largest ii?
            II = coef[2][0]
            
            # On which coefficient are we currently operating?
            index = 2
            # This is a flag that indicates the active index was used in 
            # the last loop, so it needs to be incremented.
            
            for ii in range(II,-1,-1):
                # If the current x-exponent is the same one represented in
                # the active coefficient row, then calculate q.
                if index<len(coef) and coef[index][0] == ii:
                    # Fold the current coefficient into the p expansion
                    # Update the highest derivatives first since they depend
                    # on the historical values of the lower derivatives
                    if diff>0:
                        if diff > 1:
                            pxx = 2*px + x_0 * pxx
                        px = p + x_0 * px
                    p = coef[index][1] + x_0 * p
                    index += 1
                # If the current x exponent is not represented, execute a
                # p-expansion with zero q.
                else:
                    if diff > 0:
                        if diff > 1:
                            pxx = 2*px + x_0 * pxx
                        px = p + x_0 * px
                    p = x_0 * p
                    
            # Modify the derivatives for the pre-exponnetials
            if pre!=1.:
                if diff>0:
                    if diff>1:
                        pxx = pxx*x_1*x_1 + px*x_2
                    px *= x_1
                
            # Apply the post-exponentials
            if post!=0:
                f = x**post
                if diff>0:
                    fx = post*f/x
                    if diff>1:
                        fxx = fx*(post-1)/x
                        pxx = pxx*f + 2.*px*fx + p*fxx
                    px = px*f + p*fx
                p *= f

            if len(group)==1:
                return p,px,pxx

            g += p
            if diff>0:
                gx += px
                if diff>1:
                    gxx += pxx
                    
        return g,gx,gxx
    
    
    def _iter1(self, fn, prop, y, x, Ids, xmin, xmax,
                ep=1e-6, Nmax=20, fx_index=1, verbose=False,
                 **kwarg):
        """Invert an inner routine.
        
    _iter1(fn, y, x, xmin, xmax, Ids)
    
Iteration is performed in-place on the x array.

*** Required Parameters ***
fn          The inner routine (method) to be inverted.  It must have a 
            call signature 
                f, fx0, ... = fn(x0, x1, ..., diff)
            where f is the value of fn, and fx0 is the derivative of fn
            with respect to prop0. The fx_index keyword can be used to
            change where fx is found in the returned tuple.  By default
            it is 1.
prop        The string keyword index of the property to be calculated.
y           An array of N target values for fn().
x           The result array.  It should be an N-element floating point
            array that has already been initialized.
Ids         A down-select boolean index array; only x[Ids],y[Ids] will 
            be evaluated.  This allows iteration in-place on data sets 
            where only a portion of the data require iteration.  If y is
            a floating point array with N elements, Ids must be a bool
            array with N elements.  It will specify a down-selected 
            data set with M elements, where M<=N.
xmin, xmax  Upper and lower limit arrays for the x values.  These must
            broadcastable to match x and y.  Even values outside of the
            down-select region should have legal values.
*** Optional Parameters ***
ep          Epsilon; fractional error permitted in y (default 1e-6)
Nmax        Maximum number of iterations (default 20)
fx_index    The location of the property derivative in the call 
            signature (default 1)
**kwarg     All other keyword arguments are passed directly to the 
            inner routine being inverted.

"""
        # As the iteration progresses, the number of True elements in 
        # Ids will decrease until they are all false
        # There are some important intermediate values that will also
        # require indexable arrays
        dx = np.zeros_like(y, dtype=float)
        error = np.zeros_like(dx, dtype=float)
        IooB = np.zeros_like(Ids, dtype=bool)

        arg = kwarg.copy()
        count = 0
        while Ids.any():
            # Build the new argument list
            for kw in kwarg:
                arg[kw] = kwarg[kw][Ids]
            arg[prop] = x[Ids]
            # Evaluate the funciton and isolate its derivative
            FF = fn( diff=1, **arg)
            yy = FF[0]
            yyx = FF[fx_index]
            # note that x[Ids], yy, yyx, and all the other floating 
            # intermediates are now in m-space; the sub-set of values
            # still under iteration.
            # Calculate the error, the linear change in x, and the new x
            error[Ids] = y[Ids] - yy
            dx[Ids] = error[Ids] / yyx
            if verbose:
                print x, yy, yyx, dx, Ids
            x[Ids] += dx[Ids]
            # An out-of-bounds index
            IooB = np.logical_or( x < xmin, x > xmax)
            count_oob = 0
            while IooB.any():
                dx[IooB] /= 2.
                x[IooB] -= dx[IooB]
                IooB = np.logical_or( x < xmin, x > xmax)
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


    def _ao(self, tt, dd, diff=2):
        """Dimensionless ideal gas helmholtz free energy (primative routine)
Evaluates an ideal gas equation of the form
    a = log(dd) + AOlogt*log(tt) + p(t)
    
where
    dd = d / dscale
    tt = Tscale / T
    AOlt = logt
    
This is a PRIMATIVE ROUTINE.  The arguments must already be 
nondimensionalized, and the returned values are non-dimensionalzied.
"""
        
        
        # Evaluate the temperature contributions
        A,At,Att = self._poly1(tt,self.data['AOgroup']['coef'],diff)
        Ad = 0.
        Atd = 0.
        Add = 0.
        
        # Add in the logarithmic term
        A += self.data['AOgroup']['logt'] * np.log(tt)
        # Add the logarithmic density term
        A += np.log(dd)
        if diff>0:
            dlog = self.data['AOgroup']['logt']/tt
            At += dlog
            if diff>1:
                Att += -dlog/tt
                
            dlog = 1./dd
            Ad += dlog
            if diff>1:
                Add += -dlog/tt

        return A, At, Ad, Att, Atd, Add


    def _ar(self, tt, dd, diff=2):
        """Dimensionless residual helmhotz free energy (primative routine)
Each fit in the group is of the form
    a = exp(-dd**k) * pk(tt, dd)
    
when dd = d / dscale, tt = Tscale / T
    
    A,Ad,At,Add,Adt,Att = _ar(dd, tt, order=2)

Returns the Helmholtz free energy and its derivatives up to diff.

This is a PRIMATIVE ROUTINE.  The arguments must already be 
nondimensionalized, and the returned values are non-dimensionalzied.
"""
        # Sparse polynomial evaluation is roughly 2x as fast as dense
        # polynomial evaluation for the R134a polynomials.  The 
        # explicitly defined algorithm is also roughly 2x as fast.  the
        # p_d() algorithm for R134a evaluated in about 80us on a 4 core
        # AMD A10-9700B R7

        # First evaluate the polynomial without an exponential coefficient
        A,At,Ad,Att,Atd,Add = self._poly2(tt,dd,self.data['ARgroup']['coef'][0],diff)
        
        k=0
        ddk = 1.
        for term in self.data['ARgroup']['coef'][1:]:
            p,pt,pd,ptt,ptd,pdd = self._poly2(tt, dd, term, diff)
            
            k += 1
            ddk *= dd
            e = np.exp(-ddk)
            if diff>0:
                # Calculate the derivative of exp(-dd**k) without
                # the exponential; it will be multiplied through next.
                ed = -k*ddk/dd

                if diff>1:
                    # Calculate the second derivative of exp(-dd**k) without
                    # the exponential; it will be multiplied through next.
                    edd = ed*((k-1.)/dd + ed)
                    # Multiply the exponential into all the terms
                    pdd = e*(p*edd + 2.*pd*ed + pdd)
                    ptd = e*(ptd + pt*ed)
                    ptt *= e
                    
                    Att += ptt
                    Atd += ptd
                    Add += pdd
                    
                pd = e*(p*ed + pd)
                pt *= e
                
                At += pt
                Ad += pd
                
            p *= e
            
            A += p
        
        return A,At,Ad,Att,Atd,Add
        
        
    def _dsv(self,T,diff=0):
        """Saturated vapor density
"""
        DSVt = self.data['DSVgroup']['Tscale']
        DSVd = self.data['DSVgroup']['dscale']
        tt = 1. - (T/DSVt)
        
        d,dt,dtt = self._poly1(tt, self.data['DSVgroup']['coef'], diff)
        d = DSVd * np.exp(d)
        if diff>0:
            dt *= -d / DSVt
            if diff>1:
                dtt *= d / (DSVt * DSVt)
        return d, dt, dtt
        
        
    def _dsl(self,T,diff=0):
        """Saturated liquid density
"""
        DSLt = self.data['DSLgroup']['Tscale']
        DSLd = self.data['DSLgroup']['dscale']
        tt = 1. - T/DSLt
        
        d,dt,dtt = self._poly1(tt, self.data['DSLgroup']['coef'], diff)
        d *= DSLd
        if diff>0:
            dt *= -DSLd / DSLt
            if diff>1:
                dtt  *= DSLd / (DSLt * DSLt)
        return d, dt, dtt
        
        
    def _ps(self,T,diff=0):
        """Saturation pressure (inner routine)
    ps, ps_T, ps_TT = _ps(T, diff=0)
    
Presumes temperature is in Kelvin, reports pressure in Pa
"""
        ps = 0.
        pst = 0.
        pstt = 0.
        
        PSt = self.data['PSgroup']['Tscale']
        PSp = self.data['PSgroup']['pscale']
        tt = T/PSt
        p,pt,ptt = self._poly1(1. - tt, self.data['PSgroup']['coef'], diff)
        
        invt = 1. / tt
        ps = np.exp(p*invt) * PSp
        if diff>0:
            pst = -invt*(pt + invt*p)
            if diff>1:
                pstt = pst*pst
                pstt += invt*(ptt + invt*(2.*pt + invt*(2.*p))) 
                pstt *= ps / (PSt * PSt)
            pst *= ps / PSt
        return ps, pst, pstt
        

    def _Ts(self,p):
        
        # Initialize the result array
        T = np.ones_like(p, dtype=float) * \
                0.5*(self.data['Tt'] + self.data['Tc'])
        # Create a down-select array
        Ids = np.logical_and(
                p >= self.data['pt'],
                p <= self.data['pc'])
        # Execute the iteration
        self._iter1(
                self._ps,           # Invert the saturation pressure
                'T',                # Solve for temperature
                p,                  # such that _ps(T) = p
                T,                  # The initial T values
                Ids,                # The down-select array
                self.data['Tt'],    # Minimum at the triple temp.
                self.data['Tc'])    # Maximum at the critical temp.
        return T

        
    def _p(self, T, d, diff=0):
        """Calculate pressure from (T,d) (inner routine)
    p, pt, pd = _p
    
_p() does NOT handle cases where d is "under the dome."  _p() expects
sub-critical densities to be either purely liquid or vapor.
"""
        p = 0.
        pt = 0.
        pd = 0.

        At = self.data['ARgroup']['Tscale']
        Ad = self.data['ARgroup']['dscale']
        R = self.data['R']
        # Calculate dimensionless arrays
        tt = At/T
        dd = d/Ad
        # Calculate the Helmholtz free energy
        _,_,ard,_,artd,ardd = self._ar(tt,dd,diff+1)
        temp = R*(1. + dd*ard)
        p = temp*T*d
        if diff>0:
            pt = d*(temp - R*dd*artd*tt)
            pd = T*(temp + dd*R*(ard + dd*ardd))

        return p,pt,pd
        
        
    def _d(self,T,p):
        """Density iterator - calculate density from T,p (inner routine)
T and p MUST be ndarrays
"""
        # Benchmarking shows that calls to _p() with fewer than 100
        # data points are all equivalently expensive; even when 
        # utilizing only a single thread.  As a result, iterations must
        # under no circumstances be conducted in series.  This bisection
        # algorithm acts on all valid data in parallel.
        
        # Create a down-select array
        I = np.ones_like(T, dtype=bool)
        # And initialize a solution array
        d = np.zeros_like(T, dtype=float)
        # Initialize upper and lower iteration densities
        da = np.zeros_like(T, dtype=float)
        db = np.zeros_like(T, dtype=float)
        
        # Separate out sub-critical and super-critical values for 
        # initial conditions.  For temperatures that are super-critical, 
        # use the extreme density limits of the data set.
        Itest = T>=self.data['Tc']
        da[Itest] = self.data['dlim'][0]
        db[Itest] = self.data['dlim'][1]
        d[Itest] = 0.5*(self.data['dlim'][0] + self.data['dlim'][1])
        # For temperatures that are sub-critical, detect whether the 
        # state is liquid or gaseous.  Set Itest to sub-critical.  
        Itest = np.logical_not(Itest)
        if Itest.any():
            # Now, isolate the vapor points; set the upper density to the
            # saturated vapor density FORCE Istate to be an ndarray
            Istate = np.zeros_like(T, dtype=bool)
            Istate[Itest] = p[Itest] < self._ps(T[Itest], 0)[0]
            da[Istate] = self.data['dlim'][0]
            db[Istate] = self._dsv(T[Istate], 0)[0]
            d[Istate] = 0.5*(da[Istate] + db[Istate])
            # Now, isolate the liquid points; set the lower density to the
            # saturated liquid density
            Istate[Itest] = np.logical_not(Istate[Itest])
            da[Istate] = self._dsl(T[Istate], 0)[0]
            db[Istate] = self.data['dlim'][1]
            d[Istate] = 0.5*(da[Istate] + db[Istate])
        
            # Release the memory from Istate
            del Istate

        # perform the iteration
        self._iter1(
                self._p,
                'd',
                p,
                d,
                I,
                da,
                db,
                fx_index = 2,
                T = T)
                
        return d
        
        
    def _T(self,d,p,sat=False):
        """Temperature iterator - calculate temperature from d,p (inner routine)
d and p MUST be ndarrays

    T = _T(d,p,sat=False)

Unlike _p(), _T() DOES handle cases where d is "under the dome."  These
calculations are relatively expensive, but they are necessary to the _T
inversion process.  When sat is set to True, these intermediate 
calculations are returned to prevent redundent saturation property calls

    T,dsL,dsV,Isat = _T(d,p,sat=True)
    
dsL and dsV are the saturation densities at p
Isat is a boolean index array that is True at points where d is between
    dsL and dsV.

Calling _T() should be avoided when possible, since it is one of the
more expensive iterators.  It requires iterative steps to calculate
the saturation properties in terms of pressure AND the EOS has to be
inverted to calculate T
"""
        # Benchmarking shows that calls to _p() with fewer than 100
        # data points are all equivalently expensive; even when 
        # utilizing only a single thread.  As a result, iterations must
        # under no circumstances be conducted in series.  This bisection
        # algorithm acts on all valid data in parallel.
        
        # Initialize a down-select array
        I = np.ones_like(d, dtype=bool)
        # Initialize a saturation index array
        Isat = np.zeros_like(I, dtype=bool)
        # Initialize a result array
        T = np.zeros_like(d, dtype=float)
        # Initialize upper and lower iteration densities
        Ta = np.zeros_like(d, dtype=float)
        Tb = np.zeros_like(d, dtype=float)
        # Saturaiton density arrays
        dsL = np.zeros_like(d, dtype=float)
        dsV = np.zeros_like(d, dtype=float)
        
        # Separate out sub-critical and super-critical values for 
        # initial conditions.  For pressures that are super-critical, 
        # use the extreme temperature limits of the data set.
        Itest = np.logical_and(p>=self.data['pc'], I)
        Ta[Itest] = self.data['Tlim'][0]
        Tb[Itest] = self.data['Tlim'][1]
        T[Itest] = 0.5*(self.data['Tlim'][0] + self.data['Tlim'][1])
        
        # For pressures that are sub-critical, detect whether the 
        # state is liquid or gaseous.  Set Itest to sub-critical.  
        Itest = np.logical_not(Itest)
        if Itest.any():
            # Now, identify the points in liquid, vapor, and mixed states
            # First, we'll need the saturation temperatures... this is 
            # a numerically expensive process since Ts() is iterative.
            # Let Ta temporarily be the saturation temperature
            Ta[Itest] = self._Ts(p[Itest])
            dsL[Itest] = self._dsl(Ta[Itest], 0)[0]
            dsV[Itest] = self._dsv(Ta[Itest], 0)[0]
        
            # Now, identify the liquid points
            Isat[Itest] = d[Itest] > dsL[Itest]
            # Shift the saturation temperature to Tb
            Tb[Isat] = Ta[Isat]
            Ta[Isat] = self.data['Tlim'][0]
            T[Isat] = 0.5*(Ta[Isat] + Tb[Isat])
            
            # Now, identify the vapor points
            Isat[Itest] = d[Itest] < dsV[Itest]
            # Leave Ta as the saturation temperature
            Tb[Isat] = self.data['Tlim'][1]
            T[Isat] = 0.5*(Ta[Isat] + Tb[Isat])
            
            # Now, get the saturated states
            Isat[Itest] = np.logical_and(
                    d[Itest] >= dsV[Itest],
                    d[Itest] <= dsL[Itest])
            # We now have the solution at these points.
            # Assign the value to T
            T[Isat] = Ta[Isat]
            # Put safe values in Ta and Tb
            Tb[Isat] = self.data['Tlim'][1]
            Ta[Isat] = self.data['Tlim'][0]
            # Eliminate these from the down-select array - no iteraiton required.
            I[Isat] = False
        
        self._iter1(
                self._p,
                'T',
                p,
                T,
                I,
                Ta,
                Tb,
                d=d)
        
        if sat:
            return T, dsL, dsV, Isat
        return T
        
 
    def _Tdsat_test(self, T, d):
        """Determine whether a T,d pair are under the dome (inner routine)
    sat, dsL, dsV = _Tdsat_test(T,d)
    
sat         is True for points where T,d is saturated
dsL         Saturated liquid density
dsV         Saturated vapor density

dsL and dsV are required as intermediate calculations, so they are 
returned to prevent the need for redundant saturation density funciton
calls.
"""
        dsL = np.zeros_like(T)
        dsV = np.zeros_like(T)
        I = np.asarray(T<self.data['Tc'])
        dsL[I] = self._dsl(T[I],0)[0]
        dsV[I] = self._dsv(T[I],0)[0]
        np.logical_and(I, dsL>d, out=I)
        np.logical_and(I, dsV<d, out=I)
        return I, dsL, dsV
        
        
    def _argparse(self, T=None, p=None, d=None, x=None):
        """Present a standard argument scheme for all top-level property methods
    T,d1,d2,x,I = _argparse(self, T=None, p=None, d=None, x=None)

T the temperature array in K.
d1 and d2 are densities in kg/m3.  If the conditions are under the dome,
    then d1 is the liquid, and d2 is the vapor density.
x is the quality.  If conditions are not under the dome, it can be 
    ignored.
I is a boolean array whose elements will be True for conditions that are
    under the dome.  Its values will be False otherwise.

returns temperature, density and quality at which the property is to be
evaluated as independent arrays.  When the conditions are under the dome
d1 and d2 represent the liquid and vapor densities respectively.  At all 
other conditions, x<0 and d1 == d2.
"""
        # Case out the possible combinations
        # There are 11 possible pairs: Unspecified, T, p, d, x
        # 
        # 1) Convert the inputs into the correct units
        # 2) Broadcast the arrays appropriately
        # 3) Calculate T,d1,d2,x, and I
        if T is not None:
            T = pm.units.temperature_scale(np.asarray(T), to_units='K')
            # Force T to have at least one dimension
            if T.ndim == 0:
                T.resize((1,))
            # Ensure that T is a legal value
            if ((T<self.data['Tlim'][0]).any() or 
                    (T>self.data['Tlim'][1]).any()):
                raise pm.utility.PMParamError('MP1: Pressure is out-of-bounds.')
            
            # T,p
            # If p is defined, then none of the conditions are under the
            # dome.  Use _d() to invert into density units.
            if p is not None:
                # Convert pressure
                p = pm.units.pressure(np.asarray(p), to_units='Pa')
                # Force compatible arrays
                T,p = np.broadcast_arrays(T,p)
                d1 = self._d(T,p)
                d2 = d1
                x = np.broadcast_to(-1, T.shape)
                I = np.broadcast_to(False, T.shape)
                return T, d1, d2, x, I
            
            # T,d
            # If T,d then points CAN be under the dome
            elif d is not None:
                # Convert density
                d1 = pm.units.matter(d, self.data['mw'],
                        to_units='kg')
                pm.units.volume(d1, to_units='m3', exponent=-1, inplace=True)
                # broadcast the arrays
                T,d1 = np.broadcast_arrays(T,d1)
                # Isolate the sub-critical temperatures
                I = np.asarray(T<self.data['Tc'])
                # Calculate the saturation densities
                dsL = self._dsl(T[I])[0]
                dsV = self._dsv(T[I])[0]
                # Identify the densities that are under the dome
                Isat = np.logical_and(
                            d1[I] < dsL, d1[I] > dsV)
                # Modify I to include only the points that are saturated
                I[I] = Isat
                # If there are any densities under the dome
                if Isat.any():
                    # Calculate the quality
                    x = -np.ones_like(T, dtype=float)
                    x[I] = (d1[I] - dsL[Isat]) / (dsV[Isat] - dsL[Isat])
                    # Update the densities
                    d2 = d1.copy()
                    d1[I] = dsL[Isat]
                    d2[I] = dsV[Isat]
                else:
                    d2 = d1
                    x = np.broadcast_to(-1,T.shape)
                return T, d1, d2, x, I
                
            # T,x
            # If quality is defined, then the points MUST be saturated
            elif x is not None:
                T,x,I = np.broadcast_arrays(T,x,True)
                if (T>self.data['Tc']).any():
                    raise pm.utility.PMParamError(
                        'Quality cannot be specified above the critical temperature.')
                d1 = self._dsl(T,0)[0]
                d2 = self._dsv(T,0)[0]
                return T,d1,d2,x,I
                
            # If only T is specified
            else:
                return self._argparse(T=T, p=pm.config['def_p'])
                
        # p
        # If p is the primary parameter
        elif p is not None:
            # Convert p to the correct units
            p = pm.units.pressure(p, to_units='Pa')
            # Force p to have dimension 1 or greater
            if p.ndim==0:
                p.resize((1,))
            # Ensure that p is a legal value
            if ((p<self.data['plim'][0]).any() or 
                    (p>self.data['plim'][1]).any()):
                raise pm.utility.PMParamError('MP1: Pressure is out-of-bounds.')
                
            # p,d
            # Pressure, density is an expensive combination since it
            # involves iteration to determine the saturation properties
            # AND to recover temperature.
            if d is not None:
                # Convert to kg/m3
                d1 = pm.units.matter(d, self.data['mw'], to_units='kg')
                pm.units.volume(d1, to_units='m3', exponent=-1, inplace=True)
                # Broadcast the arrays
                d1,p = np.broadcast_arrays(d1,p)
                # This one's an expensive funciton call
                # Get temperature and the saturation densities
                T,dsL,dsV,Isat = self._T(d1,p,sat=True)
                # If there are any saturated points
                if Isat.any():
                    # Calculate the quality
                    x = -np.ones_like(p, dtype=float)
                    x[Isat] = (d1[Isat] - dsL[Isat]) / (dsV[Isat] - dsL[Isat])
                    # Update the densities
                    d2 = d1.copy()
                    d1[Isat] = dsL[Isat]
                    d2[Isat] = dsV[Isat]
                else:
                    d2 = d1
                    x = np.broadcast_to(-1,T.shape)
                return T,d1,d2,x,Isat
            
            # p,x
            # If quality is defined, we are saturated
            elif x is not None:
                # Ensure that p is sub-critical
                if (p>self.data['pc']).any():
                    raise pm.utility.PMParamError('Quality cannot be specified at pressures above the critical point.')
                    
                p,x,I = np.broadcast_arrays(p,x,True)
                # T is just the saturation temperature
                T = self._Ts(p)
                d1 = self._dsl(T,0)[0]
                d2 = self._dsv(T,0)[0]
                return T, d1, d2, x, I
                
            # If only p is specified
            else:
                return self._argparse(T=pm.config['def_T'], p=p)
                
        # d
        elif d is not None:
            # d,x
            # This combination is not supported!
            # This represents an expensive inversion problem that is
            # not certain to have a solution, and it is highly unusual 
            # to specify quality AND density.
            if x is not None:
                raise pm.utility.PMParamError(
                    'Specifying properties by density and quality is not currently supported.')
                
            # If only d is specified
            else:
                return self._argparse(T=pm.config['def_T'], d=d)
                
        # If only x is specified
        elif x is not None:
            return self._argparse(T=pm.config['def_T'], x=x)
        return self._argparse(T=pm.config['def_T'], p=pm.config['def_p'])
            


    def _e(self,T,d,diff=0):
        """Internal energy (inner routine)
    e,eT,ed = _e(T,d,diff=0)
"""
        eT = 0.
        ed = 0.

        # The IG part        
        R = self.data['R']
        Tscale = self.data['AOgroup']['Tscale']
        dscale = self.data['AOgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,at,_,att,atd,_ = self._ao(tt,dd,diff+1)
        
        e = at
        if diff>0:
            eT = tt*tt*att
            ed = atd
        
        # The residual part
        Tscale = self.data['ARgroup']['Tscale']
        dscale = self.data['ARgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,at,ad,att,atd,add = self._ar(tt,dd,diff+1)
        e += at
        e *= R*Tscale
        if diff>0:
            eT += tt*tt*att
            eT *= -R
            ed += atd
            ed *= R*Tscale/dscale

        return e,eT,ed


    def _h(self,T,d,diff=0):
        """enthalpy (inner routine)
    h,hT,hd = _h(T,d,diff=0)
"""
        hT = 0.
        hd = 0.

        # The IG part        
        R = self.data['R']
        Tscale = self.data['AOgroup']['Tscale']
        dscale = self.data['AOgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,at,_,att,atd,_ = self._ao(tt,dd,diff+1)
        
        h = 1. + tt*at
        if diff>0:
            hT = 1. - tt*tt*att
            hd = tt*atd
        
        # The residual part
        Tscale = self.data['ARgroup']['Tscale']
        dscale = self.data['ARgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,at,ad,att,atd,add = self._ar(tt,dd,diff+1)
        h += dd*ad + tt*at
        h *= R*T
        if diff>0:
            hT += dd*ad - tt*(tt*att + dd*atd)
            hT *= R
            hd += ad + dd*add + tt*atd
            hd *= R*T/dscale

        return h,hT,hd


    def _s(self,T,d,diff=0):
        """entropy (inner routine)
    s,sT,sd = _s(T,d,diff=0)
"""
        sT = 0.
        sd = 0.

        # The IG part        
        R = self.data['R']
        Tscale = self.data['AOgroup']['Tscale']
        dscale = self.data['AOgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        a,at,ad,att,atd,_ = self._ao(tt,dd,diff+1)
        
        s = tt*at - a
        if diff>0:
            sT = tt*tt*att
            sd = tt*atd - ad
        
        # The residual part
        Tscale = self.data['ARgroup']['Tscale']
        dscale = self.data['ARgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        a,at,ad,att,atd,_ = self._ar(tt,dd,diff+1)
        s += tt*at - a
        s *= R
        if diff>0:
            sT += tt*tt*att
            sT *= -R/T
            sd += tt*atd - ad
            sd *= R/dscale

        return s,sT,sd
        
        
    def _cp(self,T,d):
        """Isobaric specific heat (inner routine)
    cp = _cp(T,d)
"""

        # The IG part        
        R = self.data['R']
        Tscale = self.data['AOgroup']['Tscale']
        dscale = self.data['AOgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,_,att,_,_ = self._ao(tt,dd,2)
        
        cp = -tt*tt*att
        
        # The residual part
        Tscale = self.data['ARgroup']['Tscale']
        dscale = self.data['ARgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,ad,att,atd,add = self._ar(tt,dd,2)

        temp = 1.+dd*(ad-tt*atd)
        cp += -tt*tt*att + temp*temp/(1.+dd*(2.*ad+dd*add))
        cp *= R
        return cp
        
        
    def _cv(self,T,d):
        """Isochoric specific heat (inner routine)
    cv = _cv(T,d)
"""

        # The IG part        
        R = self.data['R']
        Tscale = self.data['AOgroup']['Tscale']
        dscale = self.data['AOgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,_,att,_,_ = self._ao(tt,dd,2)
        
        cv = tt*tt*att
        
        # The residual part
        Tscale = self.data['ARgroup']['Tscale']
        dscale = self.data['ARgroup']['dscale']
        tt = Tscale / T
        dd = d / dscale
        _,_,_,att,_,_ = self._ar(tt,dd,2)

        cv += tt*tt*att
        cv *= -R
        return cv


    #               #
    # USER ROUTINES #
    #               #
    
    #               #
    # Data limits   #
    #               #
    def Tlim(self, p=None):
        """Return the temperature limits for the data set
    Tmin, Tmax = Tlim(p=None)
    
Tlim accepts pressure as an argument for extensibility, but the MP1 
class has homogeneous temperature limits.

Returns the temperature limits in [unit_temperature].
"""
        return pm.units.temperature_scale(
            np.asarray(self.data['Tlim']),
            from_units='K')
        
        
    def plim(self, T=None):
        """Returns the pressure limits for the data set
    pmin, pmax = plim(T=None)
    
plim accepts temperature as an argument for extensibility, but the MP1 
class has homogeneous pressure limits.

Returns the pressure limits in [unit_pressure]
"""
        return pm.units.pressure(
            np.asarray(self.data['plim']),
            from_units='Pa')
        
    #                               #
    # Critical and triple points    #
    #                               #
        
    def critical(self, density=False):
        """Critical point
    Tc, pc = critical()
    
To also return the density, set the 'density' keyword to True

    Tc, pc, dc = critical(density=True)
    
Returns the critical temperature, pressure, and density in 
[unit_temperature], [unit_pressure], [unit_matter/unit_volume]
"""
        if density:
            return  pm.units.temperature_scale( \
                        self.data['Tc'], from_units='K'),\
                    pm.units.pressure( \
                        self.data['pc'], from_units='Pa'), \
                    pm.units.volume(\
                        pm.units.matter( \
                            self.data['dc'], \
                            self.data['mw'], \
                            from_units='kg'),\
                        from_units='m3', exponent=-1)
                    
        return  pm.units.temperature_scale( \
                    self.data['Tc'], from_units='K'),\
                pm.units.pressure( \
                    self.data['pc'], from_units='Pa')
        
        
    def triple(self):
        """Triple point
    Tt, pt = triple()
    
Returns the triple temperature and pressure in a tuple pair in
[unit_temperature], [unit_pressure]
"""
        return  pm.units.temperature_scale( \
                    self.data['Tt'], from_units='K'),\
                pm.units.pressure( \
                    self.data['pt'], from_units='Pa')
        
    #                       #
    # Saturaiton properties #
    #                       #
    
    def ps(self, T=None):
        """Saturation pressure
    psat = ps(T)
    
Returns the saturaiton pressure in [unit_pressure]

Calls to ps() are MUCH faster than calls to Ts(), so when given a choice
ps() should always be preferred.  The MP1 class exposes ps() as a 
fundamental empirical relationship, while Ts() has to perform iterative 
numerical inversion.
"""
        if T is None:
            T = pm.config['def_T']

        # Replace T with an array of the correct units
        T = pm.units.temperature_scale(
                np.asarray(T, dtype=float), 
                to_units='K')
        # Exclude points outside the triple-critical range
        if np.logical_or( T<self.data['Tt'], T>self.data['Tc'] ).any():
            raise pm.utility.PMParamError(
                'Saturation properties are not ' +
                'available above the critical point Tc=%f K or below the '%self.data['Tc'] +
                'triple point Tt=%f K.'%self.data['Tt'] )

        # Initialize 
        p = np.ndarray(T.shape, dtype=float)
        # Set illegal values to -1
        p[I] = -1.
        # Get to work on the legal values
        I = np.logical_not(I)
        p[I] = self._ps(T[I])[0]
        return pm.units.pressure(p, from_units='Pa')
        
        
    def Ts(self,p=None):
        """Saturation temperature
    Tsat = Ts(p)
    
Uses Newton iteration to calculate Ts from the _ps() inner method
"""
        if p is None:
            p = pm.config['def_p']

        # Replace p with an array of the correct units
        p = pm.units.pressure(
                np.asarray(p, dtype=float), 
                to_units='Pa')
        # Force p to have at least 1 dimension
        if p.ndim==0:
            p.resize((1,))
                
        # Exclude points outside the triple-critical range
        if np.logical_or( p<self.data['pt'], p>self.data['pc'] ).any():
            raise pm.utility.PMParamError(
                'Saturation properties are not ' +
                'available above the critical point pc=%f bar or below the '%(self.data['pc']/1e5) +
                'triple point pt=%f bar.'%(self.data['pt']/1e5) )
        
        return pm.units.temperature_scale( \
            self._Ts(p), from_units='Pa')
        
        
    def ds(self, T=None):
        """Saturation density
    dsL, dsV = ds(T)
    
Returns the liquid (dsL) and vapor (dsV) saturation density as a 
function of temperature.  To determine saturation densities from 
pressure, first use the Ts funciton to determine the saturation 
temperature.  

    dsL, dsV = ds(Ts(p))
"""
        if T is None:
            T = pm.config['def_T']
            
        # Replace T with an array of the correct units
        T = pm.units.temperature_scale(
                np.asarray(T, dtype=float), 
                to_units='K')
        # Exclude points outside the triple-critical range
        I = np.logical_or( T<self.data['Tt'], T>self.data['Tc'] )
        if I.any():
            pm.utility.print_warning('Saturation properties are not ' +
                'available above the critical point Tc=%f K or below the '%self.data['Tc'] +
                'triple point Tt=%f K.'%self.data['Tt'] )
        I = np.logical_not(I)
                
        dsL = np.zeros_like(T)
        dsV = np.zeros_like(T)
        dsL[I] = self._dsl(T[I],diff=0)[0]
        dsV[I] = self._dsv(T[I],diff=0)[0]
        # Get a conversion factor
        conv = pm.units.matter(1., self.data['mw'],
                from_units='kg')
        conv = pm.units.volume(conv, from_units='m3', exponent=-1)
        dsL[I] *= conv
        dsV[I] *= conv
        return dsL, dsV

    def es(self, T=None, p=None):
        pass

    def hs(self, T=None, p=None):
        pass
        
    def ss(self, T=None, p=None):
        pass


    #                       #
    # EOS properties T,p,d  #
    #                       #
    
    def p(self, T=None, d=None, quality=False):
        """Pressure   p_d(T, d, quality=False)
    p = p(T, d)
    
If the optional "quality" keyword is set to True, then the two-phase
mixture quality will also be returned.

    p,x = p_d(T,d,quality=True)
    
For points that are not "under the dome" quality will be computed to be
-1.
"""
        # If d is None, then just return the default pressure
        if d is None:
            if T is None:
                return np.array(pm.config['def_p'])
            T = np.asarray(T)
            return np.broadcast_to(pm.config['def_p'], T.shape)
        
        if T is None:
            T = pm.config['def_T']
            
        # Replace T and d with arrays of the correct units
        # If T and d are already the correct units, then they will not
        # be modified or copied.
        
        T = pm.units.temperature_scale(
                T, to_units='K')
        d = pm.units.matter( 
                d, self.data['mw'], to_units='kg')
        pm.units.volume(d, to_units='m3', exponent=-1, inplace=True)
        
        # Throw an error if T is out of bounds
        if (T < self.data['Tlim'][0]).any() or \
                (T > self.data['Tlim'][1]).any():
            raise pm.utility.PMParamError('MP1: Temperature was out-of-bounds')
        
        # Broadcast
        T,d = np.broadcast_arrays(T,d)
        
        # Initialize the results
        p = np.ndarray(T.shape, dtype=float)
        if quality:
            x = np.ndarray(T.shape, dtype=float)

        # Calculate any values under the dome
        I,dsL,dsV = self._Tdsat_test(T,d)
        if I.any():
            p[I] = self._ps(T[I],diff=0)[0]
            if quality:
                x[I] = (d[I]-dsL[I])/(dsV[I]-dsL[I])

        # Now do the rest
        I = np.logical_not(I)
        p[I] = self._p(T[I],d[I],diff=0)[0]
        if quality:
            x[I] = -1
            return p,x
        return p
        
        
    def e(self, *varg, **kwarg):
        """Internal energy  e(T=None, p=None, d=None, x=None)
From any two of the provided primary properties
    
e   Int. Energy [unit_energy / unit_matter]
T   Temperature [unit_temperature]
p   Pressure    [unit_pressure]
d   Density     [unit_matter / unit_volume]
x   Quality     [dimensionless]
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        e = self._e(T,d1,0)[0]
        if I.any():
            e[I] *= (1.-x[I])
            e[I] += self._e(T,d2,0)[0] * x[I]
        # Convert the units back to user space
        pm.units.energy(e, from_units='J', inplace=True)
        pm.units.matter(e, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        return e
        
        
    def h(self, *varg, **kwarg):
        """Enthalpy  h(T=None, p=None, d=None, x=None)
From any two of the provided primary properties
    
h   Enthalpy    [unit_energy / unit_matter]
T   Temperature [unit_temperature]
p   Pressure    [unit_pressure]
d   Density     [unit_matter / unit_volume]
x   Quality     [dimensionless]
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        h = self._h(T,d1,0)[0]
        if I.any():
            h[I] *= (1.-x[I])
            h[I] += self._h(T,d2,0)[0] * x[I]
        # Convert the units back to user space
        pm.units.energy(h, from_units='J', inplace=True)
        pm.units.matter(h, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        return h


    def s(self, *varg, **kwarg):
        """Entropy  s(T=None, p=None, d=None, x=None)
From any two of the provided primary properties
    
s   Entropy     [unit_energy / unit_matter / unit_temperature]
T   Temperature [unit_temperature]
p   Pressure    [unit_pressure]
d   Density     [unit_matter / unit_volume]
x   Quality     [dimensionless]
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        s = self._s(T,d1,0)[0]
        if I.any():
            s[I] *= (1.-x[I])
            s[I] += self._s(T,d2,0)[0] * x[I]
        # Convert the units back to user space
        pm.units.energy(s, from_units='J', inplace=True)
        pm.units.matter(s, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        pm.units.temperature(s, from_units='K', 
                exponent=-1, inplace=True)
        return s


    def cp(self, *varg, **kwarg):
        """Isobaric Specific Heat  cp(T=None, p=None, d=None, x=None)
From any two of the provided primary properties
    
cp  Sp. heat    [unit_energy / unit_matter / unit_temperature]
T   Temperature [unit_temperature]
p   Pressure    [unit_pressure]
d   Density     [unit_matter / unit_volume]
x   Quality     [dimensionless]
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        cp = self._cp(T,d1)
        if I.any():
            cp[I] *= (1.-x[I])
            cp[I] += self._cp(T,d2) * x[I]
        # Convert the units back to user space
        pm.units.energy(cp, from_units='J', inplace=True)
        pm.units.matter(cp, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        pm.units.temperature(cp, from_units='K', 
                exponent=-1, inplace=True)
        return cp


    def cv(self, *varg, **kwarg):
        """Isobaric Specific Heat  cp(T=None, p=None, d=None, x=None)
From any two of the provided primary properties
    
cv  Sp. heat    [unit_energy / unit_matter / unit_temperature]
T   Temperature [unit_temperature]
p   Pressure    [unit_pressure]
d   Density     [unit_matter / unit_volume]
x   Quality     [dimensionless]
"""
        T,d1,d2,x,I = self._argparse(*varg, **kwarg)
        cv = self._cv(T,d1)
        if I.any():
            cv[I] *= (1.-x[I])
            cv[I] += self._cv(T,d2) * x[I]
        # Convert the units back to user space
        pm.units.energy(cv, from_units='J', inplace=True)
        pm.units.matter(cv, self.data['mw'], 
                from_units='kg', exponent=-1, inplace=True)
        pm.units.temperature(cv, from_units='K', 
                exponent=-1, inplace=True)
        return cv
