
def _peval(x,y,coef,diff=2):
    """Polynomial evaluation
(p, px, py, pxx, pxy, pyy) = _peval(x,y,A,diff=2)

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
