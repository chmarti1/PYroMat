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

The first two elements of the coef list indicate pre- and post-
exponents for x and y.
    coef = [[prex, prey], [postx, posty], ... ]
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
[[1,1], [0,0], [1, 1, 0.1], [0, 2, 0.2], [0, 1, 1.2], [0, 0, 0.5]]
corrsponds to the polynomial
p(x,y) = .5 + 1.2y + .2y**2 + 0.1xy

Efficient polynomial evaluation algorithms are normally restricted to
positive integer exponents, but many thermodynamic property models use 
much more interesting polynomials.  The pre- and post- exponents can be
used to acheive a much wider range of functions.

For example,
    p(x,y) = x**(-1.5) + x**(3.5)
might be expressed as a coefficient list
    [[0.5,1], [-1.5, 0], [10, 0, 1], [0, 0, 1]]
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


    def __ao(self, tt, dd, diff=2):
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
        A,At,Att = self._poly1(tt,dd,self['AOgroup']['coef'],diff)
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
        if diff>0:
            dlog = 1./dd
            Ad += dlog
            if diff>1:
                Att += -dlog/tt

        return A, At, Ad, Att, Atd, Add


    def __ar(self, tt, dd, diff=2):
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
        
        
    def _dsv(self,T,diff=2):
        """Saturated vapor density
"""
        DSVt = self.data['DSVgroup']['Tscale']
        DSVd = self.data['DSVgroup']['dscale']
        tt = 1. - (T/DSVt)
        
        d,dt,dtt = self._poly1(tt, self.data['DSVgroup']['coef'], diff)
        d = DSVd * np.exp(d)
        if diff>0:
            dt = -dt * d / DSVt
            if diff>1:
                dtt *= d / (DSVt * DSVt)
        return d, dt, dtt
        
        
    def _dsl(self,T,diff=2):
        """Saturated liquid density
"""
        DSLt = self.data['DSLgroup']['Tscale']
        DSLd = self.data['DSLgroup']['dscale']
        tt = 1. - T/DSLt
        
        d,dt,dtt = self._poly1(tt, self.data['DSLgroup']['coef'], diff)
        d *= DSLd
        if diff>0:
            dt = -dt * DSLd / DSLt
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
        
        
    def _p(self, T, d, diff=0):
        """Calculate pressure from (T,d) (inner routine)
    p, pt, pd = _p
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
        _,_,ard,_,artd,ardd = self.__ar(tt,dd,diff+1)

        temp = R * (1. + dd*ard)
        p = temp * T * d
        if diff>0:
            pt = d*(temp - R*dd*artd*At/T)
            pd = T*(temp + d*R*(ard + ardd/Ad)/Ad)

        return p,pt,pd
        
    def _d(self,T,p):
        """Density iterator - calculate density from T,p (inner routine)
"""
        # Benchmarking shows that calls to p_d() with fewer than 100
        # data points are all equivalently expensive; even when 
        # utilizing only a single thread.  As a result, iterations must
        # under no circumstances be conducted in series.
        
        # To decide on safe initial conditions, determine whether we are
        # super-critical, liquid, or vapor for each point
        dd = np.zeros_like(T)
        
        tt = self.data['ARgroup']['Tscale'] / T
        
        
        
        
    def Tlim(self, p=None):
        """Return the temperature limits for the data set
    Tmin, Tmax = Tlim(p=None)
    
Tlim accepts pressure as an argument for extensibility, but the MP1 
class has homogeneous temperature limits.
"""
        return pm.units.temperature_scale(
            np.asarray(self.data['Tlim']),
            from_units='K')
        
    def plim(self, T=None):
        """Returns the pressure limits for the data set
    pmin, pmax = plim(T=None)
    
plim accepts temperature as an argument for extensibility, but the MP1 
class has homogeneous pressure limits.
"""
        return pm.units.pressure(
            np.asarray(self.data['plim']),
            from_units='Pa')
        
    def p_d(self, T,d):
        # Replace T and d with arrays of the correct units
        # If T and d are already the correct units, then they will not
        # be modified or copied.
        T = pm.units.temperature_scale(T, to_units='K')
        d = pm.units.matter(d, self.data['mw'], to_units='kg')
        pm.units.length(d, to_units='m', exponent=-3, inplace=True)




