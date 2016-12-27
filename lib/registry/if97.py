

class if97(__basedata__):
    
    def _peval(self,x,y,A,order=2):
        """Polynomial evaluation
    (p, dpdx, dpdy) = _peval(x,y,A)
    
Evaluates a polynomial on x and y and its derivatives.
x   x value
y   y value
A   coefficient list

Each element of A is a three-element list defining a term in the 
polynomial; the x-exponent, the y-exponent, and the corresponding
coefficient.

The list,
    [[0, 0, .5], [0, 1, 1.2], [0, 2, 0.2], [1, 1, 0.1]]
corrsponds to the polynomial
    p(x,y) = .5 + 1.2y + .2y**2 + 0.1xy
"""
        # _peval() performs sanity checks to make sure it doesn't
        # lock in an infinite loop if the data is corrupted
        # LARGE is the largest reasonable exponent
        LARGE = 100        
        
        # start at the highest order term
        index = len(A)-1
        # initialize the final polynomial and its derivatives
        p = 0.  # total polynomial
        px = 0.
        py = 0.
        pxx = 0.
        pxy = 0.
        pyy = 0.
        # From here, we loop over terms of the form a*(x**m)*(y**n)
        # If a particular m,n combination is not found in the data, then
        # its coefficient is treated as zero.
        # Start off with the largest x-exponent found in the data
        m = A[index][0]
        if m>LARGE:
            raise Exception()
        # loop through the x terms so long as there is data left to process
        while index>=0:
            # For each value of m, a "sub-polynomial" q is constructed
            # on y.  As soon as the x index changes in the data, q is re-
            # initialized and the process begins again.
            q = 0.
            dq = 0.
            ddq = 0.
            # if the current x exponent corresponds to the next entry
            if A[index][0]==m:
                # grab the largest y exponent for this sub-polynomial
                n = A[index][1]
                # sanity check
                if n>LARGE:
                    raise Exception()
                # loop through the y terms so long as the x exponent is
                # unchanged.
                while index>=0 and A[index][0]==m:
                    # if the current y exponent corresponds to the next entry
                    # pull it into the polynomial
                    if n==A[index][1]:
                        if order>1:
                            ddq = ddq*y + 2.*dq
                        if order>0:
                            dq = dq*y + q
                        q = q*y + A[index][2]
                        index-=1
                    # if the next entry doesn't match the y exponent, then
                    # there is no coefficient for this term
                    else:
                        if order>1:
                            ddq = ddq*y + 2.*dq
                        if order>0:
                            dq = dq*y + q
                        q *= y
                    # decrement the exponent
                    n-=1
                    # sanity check (this one is important!)
                    # if the data are unordered, an infinite loop
                    # can occur.  This addresses the problem.
                    if n<-LARGE:
                        raise Exception()
                # restore the y exponent to the last value used
                n+=1
                # if the last y exponent wasn't zero
                if n:
                    # there is a special algorithm for n==1
                    if n==1:
                        # modify q and its derivatives
                        ddyn = 0.
                        dyn = 1.
                        yn = y
                    else:
                        ddyn = y**(n-2)
                        dyn = ddyn*y
                        yn = dyn*y
                        dyn*=n
                        ddyn*=(n-1)*n
                    if order>1:
                        ddq = ddq*yn + 2.*dq*dyn + q*ddyn
                    if order>0:
                        dq = q*dyn + dq*yn
                    q *= yn
                # update p and derivatives
                if order>1:
                    pyy = pyy*x + ddq
                    pxy = pxy*x + py
                    pxx = pxx*x + 2.*px
                if order>0:
                    py = py*x + dq
                    px = px*x + p
                p = p*x + q
            # if the next entry doesn't match the x exponent, then
            # there are no coefficients for this term
            else:
                if order>1:
                    pyy *= x
                    pxy = pxy*x + py
                    pxx = pxx*x + 2.*px
                if order>0:
                    py *= x
                    px = px*x + p
                p *= x
            
            m-=1
            # sanity check (this one is important!)
            if m<-LARGE:
                raise Exception()
        # restore the x exponent to the last value used
        m+=1
        # if the last x exponent wasn't zero
        if m:
            # modify p and its derivatives
            # there is a special algorithm for m==1
            if m==1:
                # modify q and its derivatives
                ddxm = 0.
                dxm = 1.
                xm = x
            else:
                ddxm = x**(m-2)
                dxm = ddxm*x
                xm = dxm*x
                dxm*=m
                ddxm*=(m-1)*m
            if order>1:
                pxx = pxx*xm + 2.*px*dxm + p*ddxm
                pxy = pxy*xm + py*dxm
                pyy *= xm
            if order>0:
                px = px*xm + p*dxm
                py *= xm
            p *= xm

        return p,px,py,pxx,pxy,pyy



    def _g1(self,T,p,order=2):
        """Gibbs energy in region 1
    g,gp,gt,gpp,gpt,gtt = _g1(T,p,order=2)

Calculates the dimensionless gibbs free energy and its derivatives in
region 1.  The 'order' keyword indicates the number of derivatives to
calculate.
"""
        # apply the region 1 scaling
        t = 1386./T
        pi = p/165.3
        g,gp,gt,gpp,gpt,gtt = self._peval(7.1-pi,t-1.222,self.data['r1'],order=order)
        # gp is negative because pi appears as a negative in the argument to peval
        # The effect cancels in gpp, but also appears in gpt.
        return (pi,t,g,-gp,gt,gpp,-gpt,gtt)



    def _g2(self,T,p,order=2):
        """Gibbs energy in region 2
    g,gp,gt,gpp,gpt,gtt = _g2(T,p,order=2)

Calculates the dimensionless gibbs free energy and its derivatives in
region 1.  The 'order' keyword indicates the number of derivatives to
calculate.
"""
        np = pyro.utility.np
        # apply the region 2 scaling
        t = 540./T
        pi = p/10.
        g,gp,gt,gpp,gpt,gtt = self._peval(pi,t,self.data['r2o'],order=order)
        gr,grp,grt,grpp,grpt,grtt = self._peval(pi,t-0.5,self.data['r2r'],order=order)
        LP = np.log(pi)
        g+=(gr + LP)
        if order>0:
            LP = 1./pi
            gp+=(grp + LP)
            gt+=grt
        if order>1:
            LP = -LP/pi
            gpp+=(grpp+LP)
            gpt+=grpt
            gtt+=grtt
        return pi,t,g,gp,gt,gpp,gpt,gtt



    def _f3(self,T,p,order=2,epsilon=1e-6):
        """Helmholtz free energy for region 3
    f,fx,fy,fxx,fxy,fyy,delta = _d3(T,p,order=2,epsilon=1e-6)

"""
        np = pyro.utility.np
        if T.ndim>0:
            n = np.zeros(T.shape)
            t = np.zeros(T.shape)
            f = np.zeros(T.shape)
            fx = np.zeros(T.shape)
            fy = np.zeros(T.shape)
            fxx = np.zeros(T.shape)
            fxy = np.zeros(T.shape)
            fyy = np.zeros(T.shape)
            for index in range(T.size):
                (n[index], t[index], f[index],fx[index],fy[index],fxx[index],
                 fxy[index], fyy[index]) = self._f3(T[index],p[index],epsilon=epsilon)
            return n,t,f,fx,fy,fxx,fxy,fyy
            
        # static configuration parameters
        Niter = 20
        R = self.data['R']      # ideal gas constant
        dc = self.data['dc']    # critical density
        Tc = self.data['Tc']    # critical temperature
        r3 = self.data['r3'][1:]    # poly coefficients
        A = self.data['r3'][0][2]   # natural log coefficient
        
        # initialization        
        N = 20
        # nondimensionalize parameters
        pp = p * 1e2 / (dc * R * T)   # dimensionless target pressure
        t = Tc / T              # dimensionless temperature inverse
        # create a helper funciton to calculate the pressure
        # _pfromd() quietly updates N, pt, dpt, and enew
        # it is also responsible for evaluating the curve fit
        # and its derivatives in the third region
        def _pfromd(nnew):
            # Evaluate the curve fit polynomial terms
            f,fx,fy,fxx,fxy,fyy = self._peval(nnew,t,r3)
            # Modify the function and its derivatives to include the
            # logarithmic terms.  
            f += A*np.log(nnew)
            DLN = A/nnew
            fx += DLN
            fxx -= DLN/nnew
            # Calculate the dimensionless pressure
            pc = nnew*nnew*fx
            # aggregate the values
            values = (f,fx,fy,fxx,fxy,fyy)
            return pc, values

        # use bisection to find dimensionless density
        na = 0.10       # minimum dimensionless density
        nb = 2.33       # maximum dimensionless density
        # continue until we exceed the iteration limit
        for index in range(N):
            nc = 0.5*(na+nb)
            pc,values = _pfromd(nc)
            if pc>pp:
                nb = nc
            else:
                na = nc
            
        return (nc,t) + tuple(values)



    def _g5(self,T,p,order=2):
        """Gibbs energy in region 5
    g,gp,gt,gpp,gpt,gtt = _g2(T,p,order=2)

Calculates the dimensionless gibbs free energy and its derivatives in
region 5.  The 'order' keyword indicates the number of derivatives to
calculate.
"""
        np = pyro.utility.np
        # apply the region 5 scaling
        t = 1000./T
        pi = p/10.
        g,gp,gt,gpp,gpt,gtt = self._peval(pi,t,self.data['r5o'],order=order)
        gr,grp,grt,grpp,grpt,grtt = self._peval(pi,t,self.data['r5r'],order=order)
        LP = np.log(pi)
        g+=(gr + LP)
        if order>0:
            LP = 1./pi
            gp+=(grp + LP)
            gt+=grt
        if order>1:
            LP = -LP/pi
            gpp+=(grpp+LP)
            gpt+=grpt
            gtt+=grtt
        return pi,t,g,gp,gt,gpp,gpt,gtt



    def _region(self,T,p,root=True):
        """Identify the region in the IF97 model
    r,d = mps._region(T,p)
    
Returns an array, r, matching the shape of T and identifying the region
of each element of the T,p pair.  Implicitly, T and p must be numpy 
arrays with compatible sizes.  _region() doesn't test this condition, so
failure to comply may give unpredictable results.
"""
        np = pyro.utility.np
        nan = -1
        # test for multi-element arrays
        if root:
            N = max(T.size,p.size)
            r = np.zeros((N,),dtype=int)
            for index in range(T.size):
                r[index] = self._region(T[index],p[index],root=False)
            return r
        
        T13 = 623.15
        T32 = 863.15
        T25 = 1073.15
        Tmin = 273.15
        Tmax = 2273.15
        pmax = 1000.
        p5max = 500.

        pc = self.data['pc']
        Tc = self.data['Tc']
        
        if p<0.:
            return nan
        elif T>Tmax:
            return nan
        elif T>T25:
            if p>p5max:
                return nan
            else:
                return 5
        # at all temperatures below T25, the pressure must be less than pmax
        if p>pmax:
            return nan
        elif T>T32:
            return 2
        elif T>T13:
            # these coefficients were modified from the original IF97 report's
            # B23 equation (section 4 page 5 equation 5) to produce a pressure 
            # in bar rather than MPa
            n1 = 0.34805185628969e4
            n2 = -0.11671859879975e2
            n3 = 0.10192970039326e-1
            ptest = (n3*T + n2)*T + n1
    
            if p<ptest:
                return 2
            else:
                # we're in region 3.  Calculate the density using a 
                # bisection algorithm
                pa = pmax
                pb = ptest
                return 3
        else:
            if p<self.ps(T):
                return 2
            else:
                return 1
            



    def critical(self):
        """Returns the critical point
    (Tc,pc) = critical()
"""
        return self.data['Tc'],self.data['pc']



    def triple(self):
        """Returns the tripple point
    (Tt,pt) = triple()
"""
        return self.data['Tt'],self.data['pt']



    def ps(self,T=None):
        """Saturation pressure
    ps(T)
Return the saturation pressure as a function of temperature.
"""
        np = pyro.utility.np
        if T is None:
            T = pyro.utility.get_config('def_T')
        if not isinstance(T,np.ndarray):
            T = np.array(T)
        if (T < self.data['Tt']).any():
            raise pyro.utility.PyroInputError(
            'Saturation properties are not available below the triple point.')
        if (T > self.data['Tc']).any():
            raise pyro.utility.PyroInputError(
            'Saturation properties are not available above the critical point.')
        # revert to a float if possible
        if T.size==1:
            T = np.float(T)
        
        r4 = self.data['r4']
        n10 = r4[9]
        n9 = r4[8]
        # calculate the scaled temperature
        t = T + n9/(T-n10)
        # calculate the quadratic coefficients
        a = (t + r4[0])*t + r4[1]
        b = (r4[2]*t + r4[3])*t + r4[4]
        c = (r4[5]*t + r4[6])*t + r4[7]
        pmpa = (2.*c / (-b + np.sqrt(b*b-4.*a*c)))**4
        # IF97 reports pressure in MPa
        # convert to bar
        return pmpa*10.



    def Ts(self,p=None):
        """Saturation temperature
    Ts(p)
Returns the saturation temperature as a function of pressure.
"""
        np = pyro.utility.np
        if p is None:
            p = pyro.utility.get_config('def_p')
        if not isinstance(p,np.ndarray):
            p = np.array(p)
        if (p < self.data['pt']).any():
            raise pyro.utility.PyroInputError(
            'Saturation properties are not available below the triple point.')
        if (p > self.data['pc']).any():
            raise pyro.utility.PyroInputError(
            'Saturation properties are not available above the critical point.')
        # revert to a float if possible
        if p.size==1:
            p = np.float(p)

        r4 = self.data['r4']
        n10 = r4[9]
        n9 = r4[8]
        pp = (p / 10.)**.25
        # calculate the quadratic coefficients
        a = (pp + r4[2])*pp + r4[5]
        b = (r4[0]*pp + r4[3])*pp + r4[6]
        c = (r4[1]*pp + r4[4])*pp + r4[7]
        # compute the scaled temperature
        t = 2.*c / (-b - np.sqrt(b*b-4.*a*c))
        # de-scale the temperature
        tt = n10 + t
        T = 0.5*(tt - np.sqrt(tt*tt - 4*(n9+n10*t)))
        return T


    def hs(self, T=None, p=None, tp=False):
        """Saturation enthalpy (kJ/kg)
    (hL, hV) = hs(...)
    
Saturation properties are calculated in the liquid (L) and vapor (V) 
states at the saturation condition.  Saturation properties can be 
calculated from either the temperature or the pressure or both.  If one
is provided, the other will be calculated.  If both are provided, hs 
assumes that the user has already called Ts() or ps(), and skips the 
calculation.  Most users will want to pass only one property to avoid 
the potential for error.

If neither T or p are provided, the saturation properties will use the
the 'def_p' config parameter.  To prompt any of the saturation property
functions to return their values for temperature and pressure, set the
'tp' keyword to True.

    (T,p,hL,hV) = hs(..., tp=True)
"""
        np = pyro.utility.np
        R = self.data['R']
        # ensure that one property is defined
        # use the default pressure when in doubt
        if T is None and p is None:
            p = pyro.utility.get_config('def_p')

        if T is None:
            if not isinstance(p,np.ndarray):
                p = np.array(p)
            T = self.Ts(p)
        if p is None:
            if not isinstance(T,np.ndarray):
                T = np.array(T)
            p = self.ps(T)
        if (T>623.15).any():
            print (
    "Warning::Accuracy of saturation properties above 623.15K is reduced.")
            
        pi,t,_,_,gt,_,_,_ = self._g1(T,p,order=1)
        hL = R * T * t * gt
        pi,t,_,_,gt,_,_,_ = self._g2(T,p,order=1)
        hV = R * T * t * gt
        if tp:
            return (T,p,hL,hV)
        return hL,hV



    def es(self, T=None, p=None, tp=False):
        """Saturation internal energy (kJ/kg)
    (eL, eV) = es(...)
    
Saturation properties are calculated in the liquid (L) and vapor (V) 
states at the saturation condition.  Saturation properties can be 
calculated from either the temperature or the pressure or both.  If one
is provided, the other will be calculated.  If both are provided, hs 
assumes that the user has already called Ts() or ps(), and skips the 
calculation.  Most users will want to pass only one property to avoid 
the potential for error.

If neither T or p are provided, the saturation properties will use the
the 'def_p' config parameter.  To prompt any of the saturation property
functions to return their values for temperature and pressure, set the
'tp' keyword to True.

    (T,p,eL,eV) = es(..., tp=True)
"""
        np = pyro.utility.np
        R = self.data['R']
        # ensure that one property is defined
        # use the default pressure when in doubt
        if T is None and p is None:
            p = np.array(pyro.utility.get_config('def_p'))

        if T is None:
            if not isinstance(p,np.ndarray):
                p = np.array(p)
            T = self.Ts(p)
        if p is None:
            if not isinstance(T,np.ndarray):
                T = np.array(T)
            p = self.ps(T)
        if (T>623.15).any():
            print (
    "Warning::Accuracy of saturation properties above 623.15K is reduced.")
        pi,t,_,gp,gt,_,_,_ = self._g1(T,p,order=1)
        eL = T * R * (t*gt - pi*gp)
        pi,t,_,gp,gt,_,_,_ = self._g2(T,p,order=1)
        eV = T * R * (t*gt - pi*gp)
        if tp:
            return (T,p,eL,eV)
        return eL,eV
        
        
        

    def ds(self, T=None, p=None, tp=False):
        """Saturation density (kg/m**3)
    (dL, dV) = ds(...)

Saturation properties are calculated in the liquid (L) and vapor (V) 
states at the saturation condition.  Saturation properties can be 
calculated from either the temperature or the pressure or both.  If one
is provided, the other will be calculated.  If both are provided, hs 
assumes that the user has already called Ts() or ps(), and skips the 
calculation.  Most users will want to pass only one property to avoid 
the potential for error.

If neither T or p are provided, the saturation properties will use the
the 'def_p' config parameter.  To prompt any of the saturation property
functions to return their values for temperature and pressure, set the
'tp' keyword to True.

    (T,p,dL,dV) = ds(..., tp=True)
"""
        np = pyro.utility.np
        R = self.data['R']
        # ensure that one property is defined
        # use the default pressure when in doubt
        if T is None and p is None:
            p = np.array(pyro.utility.get_config('def_p'))

        if T is None:
            if not isinstance(p,np.ndarray):
                p = np.array(p)
            T = self.Ts(p)
        if p is None:
            if not isinstance(T,np.ndarray):
                T = np.array(T)
            p = self.ps(T)
        if (T>623.15).any():
            print (
    "Warning::Accuracy of saturation properties above 623.15K is reduced.")
        pi,t,_,gp,_,_,_,_ = self._g1(T,p,order=1)
        dL = p * 100 / (R * T * pi * gp)
        pi,t,_,gp,_,_,_,_ = self._g2(T,p,order=1)
        dV = p * 100 / (R * T * pi * gp)
        if tp:
            return (T,p,dL,dV)
        return dL,dV
        
        
    def ss(self, T=None, p=None, tp=False):
        """Saturation entropy (kJ/kg/K)
    (sL, sV) = ss(...)

Saturation properties are calculated in the liquid (L) and vapor (V) 
states at the saturation condition.  Saturation properties can be 
calculated from either the temperature or the pressure or both.  If one
is provided, the other will be calculated.  If both are provided, hs 
assumes that the user has already called Ts() or ps(), and skips the 
calculation.  Most users will want to pass only one property to avoid 
the potential for error.

If neither T or p are provided, the saturation properties will use the
the 'def_p' config parameter.  To prompt any of the saturation property
functions to return their values for temperature and pressure, set the
'tp' keyword to True.

    (T,p,dL,dV) = ds(..., tp=True)
"""
        np = pyro.utility.np
        R = self.data['R']
        # ensure that one property is defined
        # use the default pressure when in doubt
        if T is None and p is None:
            p = np.array(pyro.utility.get_config('def_p'))

        if T is None:
            if not isinstance(p,np.ndarray):
                p = np.array(p)
            T = self.Ts(p)
        if p is None:
            if not isinstance(T,np.ndarray):
                T = np.array(T)
            p = self.ps(T)
        if (T>623.15).any():
            print (
    "Warning::Accuracy of saturation properties above 623.15K is reduced.")
        pi,t,g,_,gt,_,_,_ = self._g1(T,p,order=1)
        sL = R * (t*gt - g)
        pi,t,g,_,gt,_,_,_ = self._g2(T,p,order=1)
        sV = R * (t*gt - g)
        if tp:
            return (T,p,sL,sV)
        return sL,sV        
        
# To do...
# Saturation properties:
# hs()
# es()
# ds()
# ss()
# cps()
# cvs()
# General use functions:
# psolve() <--- this one will take some care



    def hsd(self, T=None, p=None, x=None):
        """Calculate enthalpy entropy and density
    (h,s,d) = hsd(T,p)

This funciton calculates these three common properties together to save
the substantial computational overhead in cases where multiple 
properties are needed.  
"""
        np = pyro.utility.np
        R = self.data['R']
        if x is None:
            T,p,h = self._vectorize(T,p,out_init=True,allow_scalar=False)
            s = h.copy()
            d = h.copy()
            r = self._region(T,p)
            # start with region 1
            I = (r==1)
            if any(I):
                pi,t,g,gp,gt,_,_,_ = self._g1(T[I],p[I],order=1)
                h[I] = R * T[I] * t * gt
                s[I] = R * (t*gt - g)
                d[I] = p[I] * 100 / (R * T[I] * pi * gp)
            # now region 2
            I = (r==2)
            if any(I):
                pi,t,g,gp,gt,_,_,_ = self._g2(T[I],p[I],order=1)
                h[I] = R*T[I] * t * gt
                s[I] = R * (t*gt - g)
                d[I] = p[I] * 100 / (R * T[I] * pi * gp)
            I = (r==3)
            if any(I):
                n,t,f,fn,ft,_,_,_ = self._f3(T[I],p[I])
                h[I] = R*T[I] * (n*fn + t*ft)
                s[I] = R * (t*ft - f)
                d[I] = self.data['dc'] * n
            I = (r==5)
            if any(I):
                pi,t,g,gp,gt,_,_,_ = self._g5(T[I],p[I],order=1)
                h[I] = R*T[I] * t * gt
                s[I] = R * (t*gt - g)
                d[I] = p[I] * 100 / (R * T[I] * pi * gp)
            if h.size==1:
                return np.float(h),np.float(s),np.float(d)
        else:
            if T is None:
                x,p = self._vectorize(x,p)
                T = self.Ts(p=p)
            if p is None:
                T,x = self._vectorize(T,x)
                p = self.ps(T=T)
            pi,t,g,gp,gt,_,_,_ = self._g1(T,p,order=1)
            hL = R * T * t * gt
            sL = R * (t*gt - g)
            dL = p * 100 / (R * T * pi * gp)
            pi,t,g,gp,gt,_,_,_ = self._g2(T,p,order=1)
            hV = R*T * t * gt
            sV = R * (t*gt - g)
            dV = p * 100 / (R * T * pi * gp)
            h = hL + (hV-hL)*x
            s = sL + (sV-sL)*x
            d = dL + (dV-dL)*x
        return h,s,d
                

    def h(self,T=None,p=None,x=None):
        """Enthalpy (kJ/kg)
    h(T,p)
"""
        np = pyro.utility.np
        # if x is unspecified, proceed as normal
        if x is None:
            T,p,h = self._vectorize(T,p,out_init=True,allow_scalar=False)
            r = self._region(T,p)
            R = self.data['R']
            # start with region 1
            I = (r==1)
            if any(I):
                pi,t,_,_,gt,_,_,_ = self._g1(T[I],p[I],order=1)
                h[I] = R * T[I] * t * gt
            # now region 2
            I = (r==2)
            if any(I):
                pi,t,_,_,gt,_,_,_ = self._g2(T[I],p[I],order=1)
                h[I] = R*T[I] * t * gt
            I = (r==3)
            if any(I):
                n,t,_,fn,ft,_,_,_ = self._f3(T[I],p[I])
                h[I] = R*T[I] * (n*fn + t*ft)
            I = (r==5)
            if any(I):
                pi,t,_,_,gt,_,_,_ = self._g5(T[I],p[I],order=1)
                h[I] = R*T[I] * t * gt
            if h.size==1:
                return np.float(h)
            return h
        # if T is unspecified
        else:
            if T is None:
                x,p = self._vectorize(x,p)
                (L,V) = self.hs(p=p)
            else:
                T,x = self._vectorize(T,x)
                (L,V) = self.hs(T=T)
            out = L + (V-L)*x
            if out.size==1:
                return np.float(out)
            return out


    def d(self,T=None,p=None,x=None):
        """Density (kg/m**3)
    d(T,p)
"""
        np = pyro.utility.np
        if x is None:
            T,p,d = self._vectorize(T,p,out_init=True,allow_scalar=False)
            r = self._region(T,p)
            R = self.data['R']
            # calculate density as specific volume
            # start with region 1
            I = (r==1)
            if any(I):
                pi,t,_,gp,_,_,_,_ = self._g1(T[I],p[I],order=1)
                d[I] = p[I] * 100 / (R * T[I] * pi * gp)
            # now region 2
            I = (r==2)
            if any(I):
                pi,t,_,gp,_,_,_,_ = self._g2(T[I],p[I],order=1)
                d[I] = p[I] * 100 / (R * T[I] * pi * gp)
            I = (r==3)
            if any(I):
                n,t,_,_,_,_,_,_ = self._f3(T[I],p[I])
                d[I] = self.data['dc'] * n
            I = (r==5)
            if any(I):
                pi,t,_,gp,_,_,_,_ = self._g5(T[I],p[I], order=1)
                d[I] = p[I] * 100 / (R * T[I] * pi * gp)
            if d.size==1:
                return np.float(d)
            return d
        else:
            if T is None:
                x,p = self._vectorize(x,p)
                (L,V) = self.ds(p=p)
            else:
                T,x = self._vectorize(T,x)
                (L,V) = self.ds(T=T)
            out = L + (V-L)*x
            if out.size==1:
                return np.float(out)
            return out


        
    def s(self,T=None,p=None,x=None):
        """Entropy (kJ/kg/K)
        s(T,p)
"""
        np = pyro.utility.np
        if x is None:
            T,p,s = self._vectorize(T,p,out_init=True,allow_scalar=False)
            r = self._region(T,p)
            R = self.data['R']
            # start with region 1
            I = (r==1)
            if any(I):
                pi,t,g,_,gt,_,_,_ = self._g1(T[I],p[I],order=1)
                s[I] = R * (t*gt - g)
            # now region 2
            I = (r==2)
            if any(I):
                pi,t,g,_,gt,_,_,_ = self._g2(T[I],p[I],order=1)
                s[I] = R * (t*gt - g)
            I = (r==3)
            if any(I):
                n,t,f,_,ft,_,_,_ = self._f3(T[I],p[I])
                s[I] = R * (t*ft - f)
            I = (r==5)
            if any(I):
                pi,t,g,_,gt,_,_,_ = self._g5(T[I],p[I], order=1)
                s[I] = R * (t*gt - g)
            if s.size==1:
                return np.float(s)
            return s
        else:
            if T is None:
                x,p = self._vectorize(x,p)
                (L,V) = self.ss(p=p)
            else:
                T,x = self._vectorize(T,x)
                (L,V) = self.ss(T=T)
            out = L + (V-L)*x
            if out.size==1:
                return np.float(out)
            return out

            
                
    
    def e(self,T=None,p=None,x=None):
        """Internal Energy (kJ/kg)
        e(T,p)
"""
        np = pyro.utility.np
        if x is None:
            T,p,e = self._vectorize(T,p,out_init=True,allow_scalar=False)
            r = self._region(T,p)
            R = self.data['R']
            # start with region 1
            I = (r==1)
            if any(I):
                pi,t,_,gp,gt,_,_,_ = self._g1(T[I],p[I],order=1)
                e[I] = T[I] * R * (t*gt - pi*gp)
            # now region 2
            I = (r==2)
            if any(I):
                pi,t,_,gp,gt,_,_,_ = self._g2(T[I],p[I],order=1)
                e[I] = T[I] * R * (t*gt - pi*gp)
            I = (r==3)
            if any(I):
                n,t,_,_,ft,_,_,_ = self._f3(T[I],p[I])
                e[I] = T[I] * R * t*ft
            I = (r==5)
            if any(I):
                pi,t,_,gp,gt,_,_,_ = self._g5(T[I],p[I], order=1)
                e[I] = T[I] * R * (t*gt - pi*gp)
            if e.size==1:
                return np.float(e)
            return e
        else:
            if T is None:
                x,p = self._vectorize(x,p)
                (L,V) = self.es(p=p)
            else:
                T,x = self._vectorize(T,x)
                (L,V) = self.es(T=T)
            out = L + (V-L)*x
            if out.size==1:
                return np.float(out)
            return out


    def cp(self,T=None,p=None):
        """Constant pressure specific heat (kJ/kg)
        cp(T,p)
"""
        np = pyro.utility.np
        T,p,cp = self._vectorize(T,p,out_init=True,allow_scalar=False)
        r = self._region(T,p)
        R = self.data['R']
        # start with region 1
        I = (r==1)
        if any(I):
            pi,t,_,_,_,_,_,gtt = self._g1(T[I],p[I],order=2)
            cp[I] = -R * t*t*gtt
        # now region 2
        I = (r==2)
        if any(I):
            pi,t,_,_,_,_,_,gtt = self._g2(T[I],p[I],order=2)
            cp[I] = -R * t*t*gtt
        I = (r==3)
        if any(I):
            n,t,_,fp,ft,fpp,fpt,ftt = self._f3(T[I],p[I])
            temp = n*fp - n*t*fpt
            temp = temp*temp/(2*n*fp + n*n*fpp)
            cp[I] = R * (-t*t*ftt + temp)
        I = (r==5)
        if any(I):
            pi,t,_,_,_,_,_,gtt = self._g5(T[I],p[I], order=2)
            cp[I] = -R * t*t*gtt
        if cp.size==1:
            return np.float(cp)
        return cp
        
        
    def cv(self,T=None,p=None):
        """Constant volume specific heat (kJ/kg)
        cv(T,p)
"""
        np = pyro.utility.np
        T,p,cv = self._vectorize(T,p,out_init=True,allow_scalar=False)
        r = self._region(T,p)
        R = self.data['R']
        # start with region 1
        I = (r==1)
        if any(I):
            pi,t,_,gp,gt,gpp,gpt,gtt = self._g1(T[I],p[I],order=2)
            temp = gp - t*gpt
            temp = temp*temp/gpp
            cv[I] = R * (temp - t*t*gtt)
        # now region 2
        I = (r==2)
        if any(I):
            pi,t,_,gp,gt,gpp,gpt,gtt = self._g2(T[I],p[I],order=2)
            temp = 1. + pi*gp - t*pi*gpt
            temp = temp*temp/(1. - pi*pi*gpp)
            cv[I] = -R * (temp + t*t*gtt)
        I = (r==3)
        if any(I):
            n,t,_,_,_,_,_,ftt = self._f3(T[I],p[I])
            cv[I] = -R * t*t*ftt
        I = (r==5)
        if any(I):
            pi,t,_,gp,gt,gpp,gpt,gtt = self._g5(T[I],p[I], order=2)
            temp = gp - t*gpt
            temp = temp*temp/gpp
            cv[I] = R * (temp - t*t*gtt)
        if cv.size==1:
            return np.float(cv)
        return cv



    def mw(self,T=None,p=None):
        """Molecular weight (kg/kmol)
    mw()
"""
        return self.data['mw']
        


