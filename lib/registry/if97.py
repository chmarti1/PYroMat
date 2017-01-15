import pyromat as pyro
import numpy as np

class if97(pyro.reg.__basedata__):


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


    def _th1(self,h,p):
        """Temperature from enthalpy and pressure
    T = _th1(h,p)

Applies the inverse relations for region 1 to calculate temperature from
enthalpy and pressure.
"""
        eta = h/2500.
        pi = p/10.
        T,_,_,_,_,_ = self._peval(pi, eta+1.,self.data['th1'],order=0)
        return T


    def _ts1(self,s,p):
        """Temperature from entropy and pressure
    T = _th1(s,p)

Applies the inverse relations for region 1 to calculate temperature from
entropy and pressure.
"""
        pi = p/10.
        T,_,_,_,_,_ = self._peval(pi, s+2.,self.data['ts1'],order=0)
        return T


    def _g2(self,T,p,order=2):
        """Gibbs energy in region 2
    g,gp,gt,gpp,gpt,gtt = _g2(T,p,order=2)

Calculates the dimensionless gibbs free energy and its derivatives in
region 1.  The 'order' keyword indicates the number of derivatives to
calculate.
"""
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


    def _th2a(self,h,p):
        """Temperature from enthalpy and pressure
    T = _th2a(h,p)

Applies the inverse relations for region 2a to calculate temperature from
enthalpy and pressure.
"""
        eta = h/2000.
        pi = p/10.
        T,_,_,_,_,_ = self._peval(pi, eta-2.1,self.data['th2a'],order=0)
        return T

    def _th2b(self,h,p):
        """Temperature from enthalpy and pressure
    T = _th2b(h,p)

Applies the inverse relations for region 2b to calculate temperature from
enthalpy and pressure.
"""
        eta = h/2000.
        pi = p/10.
        T,_,_,_,_,_ = self._peval(pi-2., eta-2.6,self.data['th2b'],order=0)
        return T

    def _th2c(self,h,p):
        """Temperature from enthalpy and pressure
    T = _th2c(h,p)

Applies the inverse relations for region 2c to calculate temperature from
enthalpy and pressure.
"""
        eta = h/2000.
        pi = p/10.
        T,_,_,_,_,_ = self._peval(pi+25., eta-1.8,self.data['th2c'],order=0)
        return T

    def _ts2a(self,s,p):
        """Temperature from entropy and pressure
    T = _ts2a(s,p)

Applies the inverse relations for region 2a to calculate temperature from
enthalpy and pressure.
"""
        sigma = s/2.
        pi = p/10.
        T,_,_,_,_,_ = self._peval(pi**.25, sigma-2.,self.data['ts2a'],order=0)
        return T

    def _ts2b(self,s,p):
        """Temperature from entropy and pressure
    T = _ts2b(s,p)

Applies the inverse relations for region 2b to calculate temperature from
enthalpy and pressure.
"""
        sigma = s/.7853
        pi = p/10.
        T,_,_,_,_,_ = self._peval(pi, 10.-sigma,self.data['ts2b'],order=0)
        return T

    def _ts2c(self,s,p):
        """Temperature from entropy and pressure
    T = _ts2c(s,p)

Applies the inverse relations for region 2c to calculate temperature from
enthalpy and pressure.
"""
        sigma = s/2.9251
        pi = p/10.
        T,_,_,_,_,_ = self._peval(pi, 2.-sigma,self.data['ts2c'],order=0)
        return T


    def _f3(self,T,p,order=2):
        """Helmholtz free energy for region 3
    f,fx,fy,fxx,fxy,fyy,delta = _d3(T,p,order=2)

"""
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
                 fxy[index], fyy[index]) = self._f3(T[index],p[index])
            return n,t,f,fx,fy,fxx,fxy,fyy
            
        # static configuration parameters
        R = self.data['R']      # ideal gas constant
        dc = self.data['dc']    # critical density
        Tc = self.data['Tc']    # critical temperature
        r3 = self.data['r3']
        A = self.data['r3ln']   # natural log coefficient
        
        # initialization        
        N = 14
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



    def _th3(self,h,p,Tinit,dinit=500.):
        """Temperature from enthalpy and pressure in region 3
    T = _th3(h,p,Tinit,dinit)

Unlike the other region evaluation functions, _th3 does NOT accept 
arrays.  It requires initial values for temperature (Tinit) and 
density (dinit).

The IF-97 document does not supply inverse relationships in regime 3.  
Instead, _th3 uses Newton iteration to match enthalpy and pressure.
In order to know that the point lies in region 3, the controlling 
algorithm will already need to have evaluated the enthalpy at the 
region 3 boundary with region 1 and region 2 with pressure p.
"""
        # Define some important constants
        R = self.data['R']      # ideal gas constant
        dc = self.data['dc']    # critical density
        Tc = self.data['Tc']    # critical temperature
        A = self.data['r3ln']   # natural log coefficient
        maxiter = 30  # maximum iterations
        epsilon = 1e-6

        # nondimensionalize parameters
        pp = p * 1e2 / (dc * R * Tc)   # dimensionless target pressure
        hh = h / (R * Tc)   # dimensionless target enthalpy (sortof)
        n = dinit/dc        # dimensionless density (delta)
        t = Tc/Tinit        # dimensionless temperature (tau)

        for count in range(maxiter):
            f,fx,fy,fxx,fxy,fyy = self._peval(n,t,self.data['r3'])
            # Modify the function and its derivatives to include the
            # logarithmic terms.  
            f += A*np.log(n)
            DLN = A/n
            fx += DLN
            fxx -= DLN/n

            ptest = n*n*fx/t - pp
            htest = n*fx/t + fy - hh
            if abs(ptest)<epsilon*pp and abs(htest)<epsilon*hh:
                return Tc/t
            dpdn = n/t * (2.*fx + n*fxx)
            dpdt = n*n/t * (fxy - fx/t)
            dhdn = fxy + (fx + n*fxx)/t
            dhdt = fyy + n/t*(fxy - fx/t)
            dx = np.linalg.solve(
                [[dpdn, dpdt],[dhdn, dhdt]], [-ptest, -htest])
            n += dx[0]
            t += dx[1]
        raise pyro.utility.PMAnalysisError('Steam _TH3 failed to converge.')



    def _ts3(self,s,p,Tinit,dinit=500.):
        """Temperature from entropy and pressure in region 3
    T = _ts3(h,p,Tinit,dinit)

Unlike the other region evaluation functions, _ts3 does NOT accept 
arrays.  It requires initial values for temperature (Tinit) and 
density (dinit).

The IF-97 document does not supply inverse relationships in regime 3.  
Instead, _th3 uses Newton iteration to match enthalpy and pressure.
In order to know that the point lies in region 3, the controlling 
algorithm will already need to have evaluated the enthalpy at the 
region 3 boundary with region 1 and region 2 with pressure p.
"""
        # Define some important constants
        R = self.data['R']      # ideal gas constant
        dc = self.data['dc']    # critical density
        Tc = self.data['Tc']    # critical temperature
        A = self.data['r3ln']   # natural log coefficient
        maxiter = 30  # maximum iterations
        epsilon = 1e-6

        # nondimensionalize parameters
        pp = p * 1e2 / (dc * R * Tc)   # dimensionless target pressure
        ss = s / R          # dimensionless target enthalpy
        n = dinit/dc        # dimensionless density (delta)
        t = Tc/Tinit        # dimensionless temperature (tau)

        for count in range(maxiter):
            f,fx,fy,fxx,fxy,fyy = self._peval(n,t,self.data['r3'])
            # Modify the function and its derivatives to include the
            # logarithmic terms.  
            f += A*np.log(n)
            DLN = A/n
            fx += DLN
            fxx -= DLN/n

            ptest = n*n*fx/t - pp
            stest = t*fy - f - ss
            if abs(ptest)<epsilon*pp and abs(stest)<epsilon*ss:
                return Tc/t
            dpdn = n/t * (2.*fx + n*fxx)
            dpdt = n*n/t * (fxy - fx/t)
            dsdn = t*fxy - fx
            dsdt = t*fyy
            dx = np.linalg.solve(
                [[dpdn, dpdt],[dsdn, dsdt]], [-ptest, -stest])
            n += dx[0]
            t += dx[1]
        raise pyro.utility.PMAnalysisError('Steam _TH3 failed to converge.')


    def _g5(self,T,p,order=2):
        """Gibbs energy in region 5
    g,gp,gt,gpp,gpt,gtt = _g2(T,p,order=2)

Calculates the dimensionless gibbs free energy and its derivatives in
region 5.  The 'order' keyword indicates the number of derivatives to
calculate.
"""
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


    def _th5(self,h,p,Tinit):
        """Temperature from enthalpy and pressure in region 5
    T = _th5(h,p,Tinit)

Unlike the other region evaluation functions, _th5 does NOT accept 
arrays.  It requires an initial value for temperature (Tinit).

The IF-97 document does not supply inverse relationships in regime 5.  
Instead, _th5 uses Newton iteration to match enthalpy and pressure.
In order to know that the point lies in region 5, the controlling 
algorithm will already need to have evaluated the enthalpy at the 
region 5 boundary with region 2 with pressure p.
"""
        # Define some important constants
        R = self.data['R']      # ideal gas constant
        ps = 10.        # pressure scale
        Ts = 1000.      # temperature scale
        maxiter = 30    # maximum iterations
        epsilon = 1e-6
        # nondimensional terms
        t = Ts/Tinit
        hh = h/R/Ts

        for count in range(maxiter):
            _,_,g,gp,gt,gpp,gpt,gtt = self._g5(T=Ts/t, p=p, order=2)
            htest = gt - hh
            if abs(htest)<epsilon*hh:
                return Ts/t
            dhdt = gtt
            t -= htest/dhdt
        raise pyro.utility.PMAnalysisError('Steam _th5() failed to converge.')


    def _ts5(self,s,p,Tinit):
        """Temperature from entropy and pressure in region 5
    T = _ts5(s,p,Tinit)

Unlike the other region evaluation functions, _ts5 does NOT accept 
arrays.  It requires an initial value for temperature (Tinit).

The IF-97 document does not supply inverse relationships in regime 5.  
Instead, _th5 uses Newton iteration to match entropy and pressure.
In order to know that the point lies in region 5, the controlling 
algorithm will already need to have evaluated the enthalpy at the 
region 5 boundary with region 2 with pressure p.
"""
        # Define some important constants
        R = self.data['R']      # ideal gas constant
        ps = 10.        # pressure scale
        Ts = 1000.      # temperature scale
        maxiter = 30    # maximum iterations
        epsilon = 1e-6
        # nondimensional terms
        t = Ts/Tinit
        ss = s/R

        for count in range(maxiter):
            _,_,g,gp,gt,gpp,gpt,gtt = self._g5(T=Ts/t, p=p, order=2)
            stest = t*gt - g - ss
            if abs(stest)<epsilon*ss:
                return Ts/t
            dsdt = t*gtt
            t -= stest/dsdt
        raise pyro.utility.PMAnalysisError('Steam _ts5() failed to converge.')


    def _b23(self,T=None,p=None):
        """Calculate the 2-3 T,p boundary
    p = _b23(T=T)
        or
    T = _b23(p=p)
Which ever is supplied (T or p), _b23 supplies the other.  Uses the B23
equations 5 and 6 modified for pressure in bar.
"""
        if T is not None:
            n = self.data['b23']
            return (n[2]*T + n[1])*T + n[0]
        elif p is not None:
            n = self.data['b23']
            return n[3] + np.sqrt((p-n[4])/n[2])
        else:
            raise Exception('_b23 requires either T or p')


    def _b2bc(self,h=None,p=None):
        """Calculate the 2b-2c h,p boundary
    h = _b2bc(p=p)
        or
    p = _b2bc(h=h)
Which ever is supplied (h or p), _b2bc supplies the other.  Uses the B2bc
equations 20 and 21 modified for pressure in bar.
"""
        if h is not None:
            n = self.data['b2bc']
            return (n[2]*h + n[1])*h + n[0]
        elif p is not None:
            n = self.data['b2bc']
            return n[3] + np.sqrt((p-n[4])/n[2])
        else:
            raise Exception('_b2bc requires either h or p')


    def _region(self,T,p,root=True):
        """Identify the region in the IF97 model
    r = mps._region(T,p)
    
Returns an array, r, matching the shape of T and identifying the region
of each element of the T,p pair.  Implicitly, T and p must be numpy 
arrays with compatible sizes.  _region() doesn't test this condition, so
failure to comply may give unpredictable results.
"""
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
            # Test pressure against the 2-3 boundary
            if p<self._b23(T=T):
                return 2
            else:
                # we're in region 3.
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
        if T is None:
            T = pyro.utility.get_config('def_T')
        if not isinstance(T,np.ndarray):
            T = np.array(T)
        if (T < self.data['Tt']).any():
            raise pyro.utility.PMParamError(
            'Saturation properties are not available below the triple point.')
        if (T > self.data['Tc']).any():
            raise pyro.utility.PMParamError(
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
        if p is None:
            p = pyro.utility.get_config('def_p')
        if not isinstance(p,np.ndarray):
            p = np.array(p)
        if (p < self.data['pt']).any():
            raise pyro.utility.PMParamError(
            'Saturation properties are not available below the triple point.')
        if (p > self.data['pc']).any():
            raise pyro.utility.PMParamError(
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
        


    def hsd(self, T=None, p=None, x=None):
        """Calculate enthalpy entropy and density
    (h,s,d) = hsd(T,p)

This funciton calculates these three common properties together to save
the substantial computational overhead in cases where multiple 
properties are needed.  
"""
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
        



    def T_h(self,h,p=None,quality=False):
        """Temperature calculated from enthalpy
    T = T_h(h)
        or
    T = T_h(h,p)
        or
    T,x = T_h(h,p=None,quality=False)

Calculates temperature from the enthalpy and pressure.

When quality is set to True, T_h also returns the quality of saturated 
conditions.  Non-saturated conditions are assigned quality -1.
"""
        h,p,T = self._vectorize(h,p,out_init=True,allow_scalar=False,def_T=-1.)
        N = h.size
        x = np.zeros(T.shape)
        # First, figure out what region we're in.  We need to use a custom 
        # region-finding routine since we're in h,p coordinates and not T,p
        # Comments describe the process relative to the inverse region diagram 
        # in the IF-97 report on page 21 (Figure 2)
        T13 = 623.15    # vertical boundary between regions 2 and 3
        p2ab = 40.      # horizontal boundary between regions 2a and 2b
        p3 = 165.292    # The bottom of region 3
        pmax = 1000.    # The top of the IF-97 region
        p5max = 500.       # The top of region 5
        T25 = 1073.15    # Left edge of region 5
        T5max = 2273.15    # Right edge of region 5

        for index in range(N):
            thish = h[index]
            thisp = p[index]
            thisT = -1.
            thisx = -1.
            if thisp < p3:
                Ts,_,hL,hV = self.hs(p=thisp,tp=True)
                if thish<hL:
                    # Region 1
                    thisT = self._th1(h=thish,p=thisp)
                elif thish<hV:
                    # Saturation
                    thisT = Ts
                    thisx = (thish-hL)/(hV-hL)
                elif thisp<p2ab:
                    h25 = self.h(T=T25,p=thisp)
                    if thish<h25:
                        # Region 2a
                        thisT = self._th2a(h=thish,p=thisp)
                    else:
                        # Region 5
                        h5 = self.h(T=T5max,p=thisp)
                        if thish>h5:
                            raise pyro.utility.PMParamError(
                            'Steam T_h(): the state is not in the IF-97 domain.')
                        thisT = T25 + (T5max-T25)*(thish-h25)/(h5-h25)
                        thisT = self._th5(h=thish, p=thisp, Tinit=thisT)
                elif thish<self._b2bc(p=thisp):
                    # Region 2c
                    thisT = self._th2c(h=thish, p=thisp)
                else:
                    h25 = self.h(T=T25,p=thisp)
                    if thish<h25:
                        # Region 2b
                        thisT = self._th2b(h=thish, p=thisp)
                    else:
                        # Region 5
                        h5 = self.h(T=T5max,p=thisp)
                        if thish>h5:
                            raise pyro.utility.PMParamError(
                            'Steam T_h(): the state is not in the IF-97 domain.')
                        thisT = T25 + (T5max-T25)*(thish-h25)/(h5-h25)
                        thisT = self._th5(h=thish, p=thisp, Tinit=thisT)
            elif thisp < pmax:
                h13 = self.h(T=T13, p=thisp)
                if thish<h13:
                    # Region 1
                    thisT = self._th1(h=thish,p=thisp)
                else:
                    T23 = self._b23(p=thisp)
                    h23 = self.h(T=T23,p=thisp)
                    if thish<h23:
                        # Region 3
                        thisT = T13 + (T23-T13)*(thish-h13)/(h23-h13)
                        thisT = self._th3(h=thish, p=thisp, Tinit=thisT)
                    elif thish<self._b2bc(p=thisp):
                        # Region 2c
                        thisT = self._th2c(h=thish,p=thisp)
                    else:
                        h25 = self.h(T=T25,p=thisp)
                        if thish<h25:
                            # Region 2b
                            thisT = self._th2b(h=thish,p=thisp)
                        elif thisp<p5max:
                            # Region 5
                            h5 = self.h(T=T5max,p=thisp)
                            if thish>h5:
                                raise pyro.utility.PMParamError(
                                'Steam T_h(): the state is not in the IF-97 ' +
                                'domain.')
                            thisT = T25 + (T5max-T25)*(thish-h25)/(h5-h25)
                            thisT = self._th5(h=thish, p=thisp, Tinit=thisT)
                        else:
                            raise pyro.utility.PMParamError(
                            'Steam T_h(): the state is not in the IF-97 domain.')
            else:
                raise pyro.utility.PMParamError(
                'Steam T_h(): pressure is above the IF-97 maximum (1000bar)')
  
            T[index] = thisT
            x[index] = thisx
        if T.size==1:
            T = np.float(T)
            x = np.float(x)
        if quality:
            return T,x
        return T



    def T_s(self,s,p=None,quality=False):
        """Temperature calculated from entropy
    T = T_s(h)
        or
    T = T_s(h,p)
        or
    T,x = T_s(s,p=None,quality=False)

Calculates temperature from the entropy and pressure.

When quality is set to True, T_s also returns the quality of saturated 
conditions.  Non-saturated conditions are assigned quality -1.
"""
        s,p,T = self._vectorize(s,p,out_init=True,allow_scalar=False,def_T=-1.)
        N = s.size
        x = np.zeros(T.shape)
        # First, figure out what region we're in.  We need to use a custom 
        # region-finding routine since we're in s,p coordinates and not T,p
        # Comments describe the process relative to the inverse region diagram 
        # in the IF-97 report on page 21 (Figure 2)
        T13 = 623.15    # vertical boundary between regions 2 and 3
        p2ab = 40.      # horizontal boundary between regions 2a and 2b
        s2bc = 5.85     # entropy boundary between regions 2b and 2c
        p3 = 165.292    # The bottom of region 3
        pmax = 1000.    # The top of the IF-97 region
        p5max = 500.       # The top of region 5
        T25 = 1073.15    # Left edge of region 5
        T5max = 2273.15    # Right edge of region 5

        for index in range(N):
            thiss = s[index]
            thisp = p[index]
            thisT = -1.
            thisx = -1.
            if thisp < p3:
                Ts,_,sL,sV = self.ss(p=thisp,tp=True)
                if thiss<sL:
                    # Region 1
                    thisT = self._ts1(s=thiss,p=thisp)
                elif thiss<sV:
                    # Saturation
                    thisT = Ts
                    thisx = (thiss-sL)/(sV-sL)
                elif thisp<p2ab:
                    s25 = self.s(T=T25,p=thisp)
                    if thiss<s25:
                        # Region 2a
                        thisT = self._ts2a(s=thiss,p=thisp)
                    else:
                        # Region 5
                        s5 = self.s(T=T5max,p=thisp)
                        if thiss>s5:
                            raise pyro.utility.PMParamError(
                            'Steam T_h(): the state is not in the IF-97 domain.')
                        thisT = T25 + (T5max-T25)*(thiss-s25)/(s5-s25)
                        thisT = self._ts5(s=thiss, p=thisp, Tinit=thisT)
                elif thiss<s2bc:
                    # Region 2c
                    thisT = self._ts2c(s=thiss, p=thisp)
                else:
                    s25 = self.s(T=T25,p=thisp)
                    if thiss<s25:
                        # Region 2b
                        thisT = self._ts2b(s=thiss, p=thisp)
                    else:
                        # Region 5
                        s5 = self.s(T=T5max,p=thisp)
                        if thiss>s5:
                            raise pyro.utility.PMParamError(
                            'Steam T_h(): the state is not in the IF-97 domain.')
                        thisT = T25 + (T5max-T25)*(thiss-s25)/(s5-s25)
                        thisT = self._ts5(s=thiss, p=thisp, Tinit=thisT)
            elif thisp < pmax:
                s13 = self.s(T=T13, p=thisp)
                if thiss<s13:
                    # Region 1
                    thisT = self._ts1(s=thiss,p=thisp)
                else:
                    T23 = self._b23(p=thisp)
                    s23 = self.s(T=T23,p=thisp)
                    if thiss<s23:
                        # Region 3
                        thisT = T13 + (T23-T13)*(thiss-s13)/(s23-s13)
                        thisT = self._ts3(s=thiss, p=thisp, Tinit=thisT)
                    elif thiss<s2bc:
                        # Region 2c
                        thisT = self._ts2c(s=thiss,p=thisp)
                    else:
                        s25 = self.s(T=T25,p=thisp)
                        if thiss<s25:
                            # Region 2b
                            thisT = self._ts2b(s=thiss,p=thisp)
                        elif thisp<p5max:
                            # Region 5
                            s5 = self.s(T=T5max,p=thisp)
                            if thiss>s5:
                                raise pyro.utility.PMParamError(
                                'Steam T_h(): the state is not in the IF-97 ' +
                                'domain.')
                            thisT = T25 + (T5max-T25)*(thiss-s25)/(s5-s25)
                            thisT = self._ts5(s=thiss, p=thisp, Tinit=thisT)
                        else:
                            raise pyro.utility.PMParamError(
                            'Steam T_h(): the state is not in the IF-97 domain.')
            else:
                raise pyro.utility.PMParamError(
                'Steam T_h(): pressure is above the IF-97 maximum (1000bar)')
  
            T[index] = thisT
            x[index] = thisx
        if T.size==1:
            T = np.float(T)
            x = np.float(x)
        if quality:
            return T,x
        return T
