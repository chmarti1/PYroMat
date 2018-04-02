import pyromat as pyro
import numpy as np

class if97(pyro.reg.__basedata__):
    """The IF-97 class

Based on the International Association for the Properties of Water and 
Steam's 1997 Industrial Formulation (IF-97).

http://www.iapws.org/relguide/IF97-Rev.html

Offers property methods:
  cp() spec. heat       (unit_energy / unit_temperature / unit_matter)
  cv() spec. heat       (unit_energy / unit_temperature / unit_matter)
  d()  density          (unit_matter / unit_volume)
  e()  internal energy  (unit_energy / unit_matter)
  h()  enthalpy         (unit_energy / unit_matter)
  gam() sp. heat ratio  (dless)
  mw() molecular weight (unit_mass / unit_molar)
  s()  entropy          (unit_energy / unit_temperature / unit_matter)

There is also a special function for simultaneously evaluating enthalpy,
entropy and density.
  hsd()  enthalpy, entropy, and density

And special methods for calculating properties along the saturation curve
  critical()    The critical point (T,p) as a tuple
  triple()      The triple point (T,p) as a tuple
  ds()          Saturation density as a tuple (dL, dV)
  es()          Saturation energy as a tuple (eL, eV)
  hs()          Saturation enthalpy as a tuple (hL, hV)
  ss()          Saturation entropy as a tuple (sL, sV)
  ps()  Saturation pressure as a function of temperature
  Ts()  Saturation temperature as a function of pressure
  
There are also routines to invert properties; e.g. calculating 
temperature from enthalpy or from entropy and pressure.
  T_h()  temperature from enthalpy
  T_s()  temperature from entropy and pressure

The limits on the IF-97 validity are reported by the Tlim() and plim()
methods.
"""



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
        nb = 2.4        # maximum dimensionless density
        # continue until we exceed the iteration limit
        for index in range(N):
            nc = 0.5*(na+nb)
            pc,values = _pfromd(nc)
            if pc>pp:
                nb = nc
            else:
                na = nc
            
        return (nc,t) + tuple(values)



    def _th3(self,h,p,Tinit=None,dinit=500.):
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
        maxiter = 200  # maximum iterations
        epsilon = 1e-5
        # This is the convergence threshold for error^2
        # Computing the square of error saves the square root
        threshold = epsilon*epsilon

        # nondimensionalize parameters
        pp = p * 1e2 / (dc * R * Tc)   # dimensionless target pressure
        hh = h / (R * Tc)   # dimensionless target enthalpy (sortof)
        n = dinit/dc        # dimensionless density (delta)
        t = Tc/Tinit        # dimensionless temperature (tau)

        n_old = n
        t_old = t
        error_old = float('inf')
        dx = np.array([0.,0.])
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

            # calculate the new square of error
            # This is the 2D distance from the target point.
            error = ptest*ptest + htest*htest

            # Test for convergence
            if error < threshold:
                return Tc/t
            # If the error in both variables is not reduced, cut the step size
            # in half and repeat: this is back-tracking along the line of 
            # descent to look for a valley.
            elif error >= error_old:
                #print t,n,htest,ptest,'*'
                dx /= 2.
                n = n_old + dx[0]
                t = t_old + dx[1]
            # Otherwise, use the typical Newton algorithm
            else:
                #print t,n,htest,ptest
                # Retire the current error and n,t values
                error_old = error
                n_old = n
                t_old = t

                dpdn = n/t * (2.*fx + n*fxx)
                dpdt = n*n/t * (fxy - fx/t)
                dhdn = fxy + (fx + n*fxx)/t
                dhdt = fyy + n/t*(fxy - fx/t)
                dx = np.linalg.solve(
                    [[dpdn, dpdt],[dhdn, dhdt]], [-ptest, -htest])
                n += dx[0]
                t += dx[1]
        raise pyro.utility.PMAnalysisError('Steam _th3 failed to converge. h=%f, p=%f'%(h,p))



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
        maxiter = 200  # maximum iterations
        epsilon = 1e-5
        # This is the convergence threshold for error^2
        # Computing the square of error saves the square root
        threshold = epsilon*epsilon

        # nondimensionalize parameters
        pp = p * 1e2 / (dc * R * Tc)   # dimensionless target pressure
        ss = s / R          # dimensionless target enthalpy
        n = dinit/dc        # dimensionless density (delta)
        t = Tc/Tinit        # dimensionless temperature (tau)

        n_old = n
        t_old = t
        error_old = float('inf')
        dx = np.array([0.,0.])
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

            # calculate the new square of error
            # This is the 2D distance from the target point.
            error = ptest*ptest + stest*stest

            # Test for convergence
            if error < threshold:
                return Tc/t
            # If the error in both variables is not reduced, cut the step size
            # in half and repeat: this is back-tracking along the line of 
            # descent to look for a valley.
            elif error >= error_old:
                #print t,n,htest,ptest,'*'
                dx /= 2.
                n = n_old + dx[0]
                t = t_old + dx[1]
            # Otherwise, use the typical Newton algorithm
            else:
                #print t,n,htest,ptest
                # Retire the current error and n,t values
                error_old = error
                n_old = n
                t_old = t

                dpdn = n/t * (2.*fx + n*fxx)
                dpdt = n*n/t * (fxy - fx/t)
                dsdn = t*fxy - fx
                dsdt = t*fyy
                dx = np.linalg.solve(
                    [[dpdn, dpdt],[dsdn, dsdt]], [-ptest, -stest])
                n += dx[0]
                t += dx[1]
        raise pyro.utility.PMAnalysisError('Steam _ts3 failed to converge. s=%f, p=%f'%(s,p))


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


    def _region(self,T,p):
        """Identify the region in the IF97 model
    r = mps._region(T,p)
    
For a scalar T and p, returns the IF-97 region index.  If a valid region
is not found, _region() returns -1.

Accepts K, bar
Returns dimensionless
"""
        nan = -1
        
        T13 = 623.15
        T32 = 863.15
        T25 = 1073.15
        Tmin = 273.15
        Tmax = 2273.15
        pmax = 1000.
        p5max = 500.

        pc = self.data['pc']
        Tc = self.data['Tc']
        
        if p<0. or p>pmax:
            return nan
        elif T>Tmax:
            return nan
        elif T>T25:
            if p>p5max:
                return nan
            else:
                return 5
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
            if p<self._ps(T):
                return 2
            else:
                return 1
        raise Exception('Unhandled exception in IF97 _region()')


    def _ps(self, T):
        """Saturation pressure
Accepts K
Returns bar
"""
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


    def _Ts(self, p):
        """Saturation temperature
Accepts bar
Returns K
"""
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



    def Tlim(self,p=None):
        """Returns the upper and lower temperature limit
    (Tmin, Tmax) = Tlim(p)

The IF-97 report has a piece-wise limit based on pressure.  If pressure
is omitted, it obeys the same defaults as with any other property.
If all p are in the same range, then Tmax will be compressed into a 
scalar.  Otherwise, Tmax will be an array with the same dimension as p.
Tmin is always the same value, regardless of p, so it will always be
a scalar.
"""
        if p is None:
            p = pyro.config['def_p']
        if not isinstance(p,np.ndarray):
            p = np.array(p)
        p = pyro.units.pressure(p,to_units='bar')

        # First, parse out the Tmax limit        
        # Test the range for each element of p
        RR = (p > 500.)

        bigT = pyro.units.temperature_scale(2273.15, from_units='K')
        smallT = pyro.units.temperature_scale(1073.15, from_units='K')
        # If all of the pressures are in one of the ranges, let Tmax be a scalar
        if RR.all():
            Tmax = smallT
        elif not RR.any():
            Tmax = bigT
        else:
            # Otherwise, iterate element by element
            # This is more complicated than boolean indexing, but it handles
            # the scalar case gracefully
            it = np.nditer((None,RR), 
				    op_flags=[['readwrite','allocate'],['readonly','copy']], 
					    op_dtypes='float')

            for Tmax_,RR_ in it:
                if RR_:
                    Tmax_[...] = smallT
                else:
                    Tmax_[...] = bigT

            Tmax = it.operands[0]
        # Tmin is always the same regardless, so it is always a scalar
        Tmin = pyro.units.temperature_scale(273.15, from_units='K')
        return Tmin, Tmax


    def plim(self,T=None):
        """Returns the upper and lower pressure limit
    (pmin, pmax) = plim(T)

The IF-97 report has a piece-wise limit based on Temperature.  If T
is omitted, it obeys the same defaults as with any other property.
If all T are in the same range, then pmax will be compressed into a 
scalar.  Otherwise, pmax will be an array with the same dimension as T.
pmin is always the same value, regardless of p, so it will always be
a scalar.
"""
        if T is None:
            T = pyro.config['def_T']
        if not isinstance(T,np.ndarray):
            T = np.array(T)
        T = pyro.units.temperature_scale(T,to_units='K')

        # First, parse out the Tmax limit        
        # Test the range for each element of p
        RR = (T > 1073.15)

        bigp = pyro.units.pressure(1000., from_units='bar')
        smallp = pyro.units.pressure(500., from_units='bar')
        # If all of the pressures are in one of the ranges, let Tmax be a scalar
        if RR.all():
            pmax = smallp
        elif not RR.any():
            pmax = bigp
        else:
            # Otherwise, iterate element by element
            # This is more complicated than boolean indexing, but it handles
            # the scalar case gracefully
            it = np.nditer((None,RR), 
				    op_flags=[['readwrite','allocate'],['readonly','copy']], 
					    op_dtypes='float')

            for pmax_,RR_ in it:
                if RR_:
                    pmax_[...] = smallp
                else:
                    pmax_[...] = bigp

            pmax = it.operands[0]
        # pmin is always the same regardless, so it is always a scalar
        pmin = pyro.units.pressure(self.data['pt'], from_units='bar')
        return pmin, pmax


    def critical(self):
        """Returns the critical point
    (Tc,pc) = critical()

Accepts None
Returns unit_temperature, unit_pressure
"""
        T = pyro.units.temperature_scale(self.data['Tc'], from_units='K')
        p = pyro.units.pressure(self.data['pc'], from_units='bar')
        return T,p



    def triple(self):
        """Returns the tripple point
    (Tt,pt) = triple()

Accepts None
Returns unit_temperature, unit_pressure
"""
        T = pyro.units.temperature_scale(self.data['Tt'], from_units='K')
        p = pyro.units.pressure(self.data['pt'], from_units='bar')
        return T,p



    def ps(self,T=None):
        """Saturation pressure
    ps(T)
Return the saturation pressure as a function of temperature.

Accepts unit_temperature
Returns unit_pressure
"""
        if T is None:
            T = pyro.config['def_T']
        if not isinstance(T,np.ndarray):
            T = np.array(T)

        T = pyro.units.temperature_scale(T, to_units='K')

        if (T < self.data['Tt']).any():
            raise pyro.utility.PMParamError(
            'Saturation properties are not available below the triple point.')
        if (T > self.data['Tc']).any():
            raise pyro.utility.PMParamError(
            'Saturation properties are not available above the critical point.')
        
        return pyro.units.pressure(self._ps(T), from_units='bar')



    def Ts(self,p=None):
        """Saturation temperature
    Ts(p)
Returns the saturation temperature as a function of pressure.

Accepts unit_pressure
Returns unit_temperature
"""
        if p is None:
            p = pyro.config['def_p']
        if not isinstance(p,np.ndarray):
            p = np.array(p)

        p = pyro.units.pressure(p, to_units='bar')

        if (p < self.data['pt']).any():
            raise pyro.utility.PMParamError(
            'Saturation properties are not available below the triple point.')
        if (p > self.data['pc']).any():
            raise pyro.utility.PMParamError(
            'Saturation properties are not available above the critical point.')

        T = self._Ts(p)
        return pyro.units.temperature_scale(T, from_units='K')


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

Accepts unit_temperature
        unit_pressure
Returns (unit_temperature)
        (unit_pressure)
        unit_energy / unit_matter 
"""
        R = self.data['R']
        # ensure that one property is defined
        # use the default pressure when in doubt
        if T is None and p is None:
            p = pyro.config['def_p']

        if T is not None:
            TT = pyro.units.temperature_scale(T, to_units='K')
            if (TT < self.data['Tt']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available below the triple point.')
            if (TT > self.data['Tc']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available above the critical point.')

        if p is not None:
            pp = pyro.units.pressure(p, to_units='bar')
            if (pp < self.data['pt']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available below the triple point.')
            if (pp > self.data['pc']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available above the critical point.')
        else:
            pp = self._ps(TT)
            p = pyro.units.pressure(pp, from_units='bar')

        if T is None:
            TT = self._Ts(pp)
            T = pyro.units.temperature_scale(TT, from_units='K')

        if (TT>623.15).any():
            pyro.utility.print_warning(
    "Accuracy of steam saturation properties above 623.15K is reduced.")

        scale = pyro.units.energy(from_units='kJ')
        scale = pyro.units.matter(scale, self.data['mw'], from_units='kg', exponent=-1)

        pi,t,_,_,gt,_,_,_ = self._g1(TT,pp,order=1)
        hL = scale * R * TT * t * gt
        pi,t,_,_,gt,_,_,_ = self._g2(TT,pp,order=1)
        hV = scale * R * TT * t * gt
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

Accepts unit_temperature
        unit_pressure
Returns (unit_temperature)
        (unit_pressure)
        unit_energy / unit_matter 
"""
        R = self.data['R']
        # ensure that one property is defined
        # use the default pressure when in doubt
        if T is None and p is None:
            p = pyro.config['def_p']

        if T is not None:
            TT = pyro.units.temperature_scale(T, to_units='K')
            if (TT < self.data['Tt']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available below the triple point.')
            if (TT > self.data['Tc']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available above the critical point.')

        if p is not None:
            pp = pyro.units.pressure(p, to_units='bar')
            if (pp < self.data['pt']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available below the triple point.')
            if (pp > self.data['pc']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available above the critical point.')
        else:
            pp = self._ps(TT)
            p = pyro.units.pressure(pp, from_units='bar')

        if T is None:
            TT = self._Ts(pp)
            T = pyro.units.temperature_scale(TT, from_units='K')

        if (TT>623.15).any():
            pyro.utility.print_warning(
    "Accuracy of steam saturation properties above 623.15K is reduced.")

        scale = pyro.units.energy(from_units='kJ')
        scale = pyro.units.matter(scale, self.data['mw'], from_units='kg', exponent=-1)

        pi,t,_,gp,gt,_,_,_ = self._g1(TT,pp,order=1)
        eL = scale * TT * R * (t*gt - pi*gp)
        pi,t,_,gp,gt,_,_,_ = self._g2(TT,pp,order=1)
        eV = scale * TT * R * (t*gt - pi*gp)
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

Accepts unit_temperature
        unit_pressure
Returns (unit_temperature)
        (unit_pressure)
        unit_matter / unit_volume
"""
        R = self.data['R']
        # ensure that one property is defined
        # use the default pressure when in doubt
        if T is None and p is None:
            p = pyro.config['def_p']

        if T is not None:
            TT = pyro.units.temperature_scale(T, to_units='K')
            if (TT < self.data['Tt']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available below the triple point.')
            if (TT > self.data['Tc']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available above the critical point.')

        if p is not None:
            pp = pyro.units.pressure(p, to_units='bar')
            if (pp < self.data['pt']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available below the triple point.')
            if (pp > self.data['pc']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available above the critical point.')
        else:
            pp = self._ps(TT)
            p = pyro.units.pressure(pp, from_units='bar')

        if T is None:
            TT = self._Ts(pp)
            T = pyro.units.temperature_scale(TT, from_units='K')

        if (TT>623.15).any():
            pyro.utility.print_warning(
    "Accuracy of steam saturation properties above 623.15K is reduced.")

        scale = pyro.units.volume(from_units='m3',exponent=-1)
        scale = pyro.units.matter(scale, self.data['mw'], from_units='kg')

        pi,t,_,gp,_,_,_,_ = self._g1(TT,pp,order=1)
        dL = scale * pp * 100 / (R * TT * pi * gp)
        pi,t,_,gp,_,_,_,_ = self._g2(TT,pp,order=1)
        dV = scale * pp * 100 / (R * TT * pi * gp)
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

Accepts unit_temperature
        unit_pressure
Returns (unit_temperature)
        (unit_pressure)
        unit_energy / unit_temperature / unit_matter
"""
        R = self.data['R']
        # ensure that one property is defined
        # use the default pressure when in doubt
        if T is None and p is None:
            p = pyro.config['def_p']

        if T is not None:
            TT = pyro.units.temperature_scale(T, to_units='K')
            if (TT < self.data['Tt']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available below the triple point.')
            if (TT > self.data['Tc']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available above the critical point.')

        if p is not None:
            pp = pyro.units.pressure(p, to_units='bar')
            if (pp < self.data['pt']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available below the triple point.')
            if (pp > self.data['pc']).any():
                raise pyro.utility.PMParamError(
                'Saturation properties are not available above the critical point.')
        else:
            pp = self._ps(TT)
            p = pyro.units.pressure(pp, from_units='bar')

        if T is None:
            TT = self._Ts(pp)
            T = pyro.units.temperature_scale(TT, from_units='K')

        if (TT>623.15).any():
            pyro.utility.print_warning(
    "Accuracy of steam saturation properties above 623.15K is reduced.")

        scale = pyro.units.energy(from_units='kJ')
        scale = pyro.units.matter(scale, self.data['mw'], from_units='kg', exponent=-1)

        pi,t,g,_,gt,_,_,_ = self._g1(TT,pp,order=1)
        sL = scale * R * (t*gt - g)
        pi,t,g,_,gt,_,_,_ = self._g2(TT,pp,order=1)
        sV = scale * R * (t*gt - g)
        if tp:
            return (T,p,sL,sV)
        return sL,sV        
        


    def hsd(self, T=None, p=None, x=None):
        """Calculate enthalpy entropy and density
    (h,s,d) = hsd(T,p)
    (h,s,d) = hsd(T,x)
    (h,s,d) = hsd(p,x)

This funciton calculates these three common properties together to save
the substantial computational overhead in cases where multiple 
properties are needed.

T   temperature
p   pressure
x   quality (mass fraction in vapor phase)

Accepts unit_temperature
        unit_pressure
        dimensionless
Returns unit_energy / unit_matter
        unit_energy / unit_matter / unit_temperature
        unit_matter / unit_volume
"""
        def_T = False
        if T is None:
            T = pyro.config['def_T']
            def_T = True

        def_p = False
        if p is None:
            p = pyro.config['def_p']
            def_p = True
        pscale = pyro.units.pressure(to_units='bar')

        if x is None:
            x = -1.

        R = self.data['R']

        it = np.nditer((None,None,None,T,p,x), 
				op_flags=[['readwrite','allocate'],['readwrite','allocate'],['readwrite','allocate'],
					['readonly','copy'],['readonly','copy'],['readonly','copy']], 
					op_dtypes='float')

        for h_,s_,d_,T_,p_,x_ in it:
            TT = pyro.units.temperature_scale(T_, to_units='K')
            pp = p_ * pscale

            # If x is unspecified
            if x_<0.:

                # Detect the region
                r = self._region(TT,pp)
                # Case out the region indices
                if r==1:
                    pi,t,g,gp,gt,_,_,_ = self._g1(TT,pp,order=1)
                    h_[...] = R * TT * t * gt
                    s_[...] = R * (t*gt - g)
                    d_[...] = pp * 100 / (R * TT * pi * gp)
                elif r==2:
                    pi,t,g,gp,gt,_,_,_ = self._g2(TT,pp,order=1)
                    h_[...] = R*TT * t * gt
                    s_[...] = R * (t*gt - g)
                    d_[...] = pp * 100 / (R * TT * pi * gp)
                elif r==3:
                    n,t,f,fn,ft,_,_,_ = self._f3(TT,pp)
                    h_[...] = R*TT * (n*fn + t*ft)
                    s_[...] = R * (t*ft - f)
                    d_[...] = self.data['dc'] * n
                elif r==5:
                    pi,t,g,gp,gt,_,_,_ = self._g5(TT,pp,order=1)
                    h_[...] = R*TT * t * gt
                    s_[...] = R * (t*gt - g)
                    d_[...] = pp * 100 / (R * TT * pi * gp)
                else:
                    raise pyro.utility.PMParamError('Invalid property combination T=%f K, p=%f bar'%(TT,pp))

            else:
                # If T was unspecified
                if def_T:
                    # Override the default T with the saturation T
                    TT = self._Ts(pp)
                # If pressure was unspecified, but temperature WAS
                elif pp<0.:
                    # Override the default p with the saturation p
                    pp = self._ps(TT)

                pi,t,g,gp,gt,_,_,_ = self._g1(TT,pp,order=1)
                hL = R * TT * t * gt
                sL = R * (t*gt - g)
                vL = (R * TT * pi * gp) / (pp * 100)
                pi,t,g,gp,gt,_,_,_ = self._g2(TT,pp,order=1)
                hV = R * TT * t * gt
                sV = R * (t*gt - g)
                vV = (R * TT * pi * gp) / (pp * 100)

                h_[...] = hL + (hV-hL)*x_
                s_[...] = sL + (sV-sL)*x_
                d_[...] = 1./(vL + (vV-vL)*x_)

        h,s,d = it.operands[0:3]

        # Convert the results
        hscale = pyro.units.energy(from_units='kJ')
        hscale = pyro.units.matter(hscale,self.data['mw'],from_units='kg',exponent=-1)
        sscale = hscale
        sscale = pyro.units.temperature(sscale,from_units='K')
        dscale = pyro.units.volume(from_units='m3',exponent=-1)
        dscale = pyro.units.matter(dscale, self.data['mw'], from_units='kg')
        
        return hscale*h, sscale*s, dscale*d


    def h(self,T=None,p=None,x=None):
        """Enthalpy
    h(T=None,p=None,x=None)

T   Temperature
p   pressure
x   quality

Accepts the temperature, pressure, or quality in any of the three 
possible combinations.  If none of the arguments are supplied, then 
def_T,def_p are used, and x is ignored.  If only x is supplied, p is 
presumed to be def_p, and T is calculated to be the saturation 
temperature.

Accepts unit_temperature
        unit_pressure
        dimensionless
Returns unit_energy / unit_matter
"""
        def_T = False
        if T is None:
            T = pyro.config['def_T']
            def_T = True

        def_p = False
        if p is None:
            p = pyro.config['def_p']
            def_p = True
        pscale = pyro.units.pressure(to_units='bar')

        if x is None:
            x = -1.

        R = self.data['R']

        it = np.nditer((None,T,p,x), 
                op_flags=[['readwrite','allocate'],['readonly','copy'],
                    ['readonly','copy'],['readonly','copy']], op_dtypes='float')

        for h_,T_,p_,x_ in it:
            TT = pyro.units.temperature_scale(T_, to_units='K')
            pp = p_ * pscale

            # If x is unspecified
            if x_<0.:

                # Detect the region
                r = self._region(TT,pp)
                # Case out the region indices
                if r==1:
                    pi,t,g,gp,gt,_,_,_ = self._g1(TT,pp,order=1)
                    h_[...] = R * TT * t * gt
                elif r==2:
                    pi,t,g,gp,gt,_,_,_ = self._g2(TT,pp,order=1)
                    h_[...] = R*TT * t * gt
                elif r==3:
                    n,t,f,fn,ft,_,_,_ = self._f3(TT,pp)
                    h_[...] = R*TT * (n*fn + t*ft)
                elif r==5:
                    pi,t,g,gp,gt,_,_,_ = self._g5(TT,pp,order=1)
                    h_[...] = R*TT * t * gt
                else:
                    raise pyro.utility.PMParamError('Invalid property combination T=%f K, p=%f bar'%(TT,pp))


            else:
                # If T was unspecified
                if def_T:
                    # Override the default T with the saturation T
                    TT = self._Ts(pp)
                # If pressure was unspecified, but temperature WAS
                elif def_p:
                    # Override the default p with the saturation p
                    pp = self._ps(TT)

                pi,t,g,gp,gt,_,_,_ = self._g1(TT,pp,order=1)
                hL = R * TT * t * gt
                pi,t,g,gp,gt,_,_,_ = self._g2(TT,pp,order=1)
                hV = R * TT * t * gt

                h_[...] = hL + (hV-hL)*x_

        h = it.operands[0]

        # Convert the results
        hscale = pyro.units.energy(from_units='kJ')
        hscale = pyro.units.matter(hscale,self.data['mw'],from_units='kg',exponent=-1)

        return hscale*h



    def d(self,T=None,p=None,x=None):
        """Density
    d(T=None,p=None,x=None)

T   Temperature
p   pressure
x   quality

Accepts the temperature, pressure, or quality in any of the three 
possible combinations.  If none of the arguments are supplied, then 
def_T,def_p are used, and x is ignored.  If only x is supplied, p is 
presumed to be def_p, and T is calculated to be the saturation 
temperature.

Accepts unit_temperature
        unit_pressure
        dimensionless
Returns unit_matter / unit_volume
"""
        def_T = False
        if T is None:
            T = pyro.config['def_T']
            def_T = True

        def_p = False
        if p is None:
            p = pyro.config['def_p']
            def_p = True
        pscale = pyro.units.pressure(to_units='bar')

        if x is None:
            x = -1.

        R = self.data['R']

        it = np.nditer((None,T,p,x), 
                op_flags=[['readwrite','allocate'],['readonly','copy'],
                    ['readonly','copy'],['readonly','copy']], op_dtypes='float')

        for d_,T_,p_,x_ in it:
            TT = pyro.units.temperature_scale(T_, to_units='K')
            pp = p_ * pscale

            # If x is unspecified
            if x_<0.:

                # Detect the region
                r = self._region(TT,pp)
                # Case out the region indices
                if r==1:
                    pi,t,g,gp,gt,_,_,_ = self._g1(TT,pp,order=1)
                    d_[...] = pp * 100 / (R * TT * pi * gp)
                elif r==2:
                    pi,t,g,gp,gt,_,_,_ = self._g2(TT,pp,order=1)
                    d_[...] = pp * 100 / (R * TT * pi * gp)
                elif r==3:
                    n,t,f,fn,ft,_,_,_ = self._f3(TT,pp)
                    d_[...] = self.data['dc'] * n
                elif r==5:
                    pi,t,g,gp,gt,_,_,_ = self._g5(TT,pp,order=1)
                    d_[...] = pp * 100 / (R * TT * pi * gp)
                else:
                    raise pyro.utility.PMParamError('Invalid property combination T=%f K, p=%f bar'%(TT,pp))


            else:
                # If T was unspecified
                if def_T:
                    # Override the default T with the saturation T
                    TT = self._Ts(pp)
                # If pressure was unspecified, but temperature WAS
                elif def_p:
                    # Override the default p with the saturation p
                    pp = self._ps(TT)

                pi,t,g,gp,gt,_,_,_ = self._g1(TT,pp,order=1)
                vL = (R * TT * pi * gp) / (pp * 100)
                pi,t,g,gp,gt,_,_,_ = self._g2(TT,pp,order=1)
                vV = (R * TT * pi * gp) / (pp * 100)

                d_[...] = 1./(vL + (vV-vL)*x_)

        d = it.operands[0]

        # Convert the results
        dscale = pyro.units.volume(from_units='m3',exponent=-1)
        dscale = pyro.units.matter(dscale, self.data['mw'], from_units='kg')
        
        return dscale*d



    def s(self,T=None,p=None,x=None):
        """Entropy
    s(T=None,p=None,x=None)

T   Temperature
p   pressure
x   quality

Accepts the temperature, pressure, or quality in any of the three 
possible combinations.  If none of the arguments are supplied, then 
def_T,def_p are used, and x is ignored.  If only x is supplied, p is 
presumed to be def_p, and T is calculated to be the saturation 
temperature.

Accepts unit_temperature
        unit_pressure
        dimensionless
Returns unit_energy / unit_matter / unit_temperature
"""
        def_T = False
        if T is None:
            T = pyro.config['def_T']
            def_T = True

        def_p = False
        if p is None:
            p = pyro.config['def_p']
            def_p = True
        pscale = pyro.units.pressure(to_units='bar')

        if x is None:
            x = -1.

        R = self.data['R']

        it = np.nditer((None,T,p,x), 
                op_flags=[['readwrite','allocate'],['readonly','copy'],
                    ['readonly','copy'],['readonly','copy']], op_dtypes='float')

        for s_,T_,p_,x_ in it:
            TT = pyro.units.temperature_scale(T_, to_units='K')
            pp = p_ * pscale

            # If x is unspecified
            if x_<0.:

                # Detect the region
                r = self._region(TT,pp)
                # Case out the region indices
                if r==1:
                    pi,t,g,gp,gt,_,_,_ = self._g1(TT,pp,order=1)
                    s_[...] = R * (t*gt - g)
                elif r==2:
                    pi,t,g,gp,gt,_,_,_ = self._g2(TT,pp,order=1)
                    s_[...] = R * (t*gt - g)
                elif r==3:
                    n,t,f,fn,ft,_,_,_ = self._f3(TT,pp)
                    s_[...] = R * (t*ft - f)
                elif r==5:
                    pi,t,g,gp,gt,_,_,_ = self._g5(TT,pp,order=1)
                    s_[...] = R * (t*gt - g)
                else:
                    raise pyro.utility.PMParamError('Invalid property combination T=%f K, p=%f bar'%(TT,pp))


            else:
                # If T was unspecified
                if def_T:
                    # Override the default T with the saturation T
                    TT = self._Ts(pp)
                # If pressure was unspecified, but temperature WAS
                elif def_p:
                    # Override the default p with the saturation p
                    pp = self._ps(TT)

                pi,t,g,gp,gt,_,_,_ = self._g1(TT,pp,order=1)
                sL = R * (t*gt - g)
                pi,t,g,gp,gt,_,_,_ = self._g2(TT,pp,order=1)
                sV = R * (t*gt - g)

                s_[...] = sL + (sV-sL)*x_

        s = it.operands[0]

        # Convert the results
        sscale = pyro.units.energy(from_units='kJ')
        sscale = pyro.units.matter(sscale,self.data['mw'],from_units='kg',exponent=-1)
        sscale = pyro.units.temperature(sscale,from_units='K')
        
        return sscale*s



    def e(self,T=None,p=None,x=None):
        """Internal Energy
        e(T=None,p=None,x=None)

T   Temperature
p   pressure
x   quality

Accepts the temperature, pressure, or quality in any of the three 
possible combinations.  If none of the arguments are supplied, then 
def_T,def_p are used, and x is ignored.  If only x is supplied, p is 
presumed to be def_p, and T is calculated to be the saturation 
temperature.

Accepts unit_temperature
        unit_pressure
        dimensionless
Returns unit_energy / unit_matter
"""
        def_T = False
        if T is None:
            T = pyro.config['def_T']
            def_T = True

        def_p = False
        if p is None:
            p = pyro.config['def_p']
            def_p = True
        pscale = pyro.units.pressure(to_units='bar')

        if x is None:
            x = -1.

        R = self.data['R']

        it = np.nditer((None,T,p,x), 
                op_flags=[['readwrite','allocate'],['readonly','copy'],
                    ['readonly','copy'],['readonly','copy']], op_dtypes='float')

        for e_,T_,p_,x_ in it:
            TT = pyro.units.temperature_scale(T_,to_units='K')
            pp = p_ * pscale

            # If x is unspecified
            if x_<0.:

                # Detect the region
                r = self._region(TT,pp)
                # Case out the region indices
                if r==1:
                    pi,t,_,gp,gt,_,_,_ = self._g1(TT,pp,order=1)
                    e_[...] = TT * R * (t*gt - pi*gp)
                elif r==2:
                    pi,t,_,gp,gt,_,_,_ = self._g2(TT,pp,order=1)
                    e_[...] = TT * R * (t*gt - pi*gp)
                elif r==3:
                    n,t,_,_,ft,_,_,_ = self._f3(TT,pp)
                    e_[...] = TT * R * t*ft
                elif r==5:
                    pi,t,_,gp,gt,_,_,_ = self._g5(TT,pp,order=1)
                    e_[...] = TT * R * (t*gt - pi*gp)
                else:
                    raise pyro.utility.PMParamError('Invalid property combination T=%f K, p=%f bar'%(TT,pp))


            else:
                # If T was unspecified
                if def_T:
                    # Override the default T with the saturation T
                    TT = self._Ts(pp)
                # If pressure was unspecified, but temperature WAS
                elif def_p:
                    # Override the default p with the saturation p
                    pp = self._ps(TT)

                pi,t,_,gp,gt,_,_,_ = self._g1(TT,pp,order=1)
                eL = TT * R * (t*gt - pi*gp)
                pi,t,g,gp,gt,_,_,_ = self._g2(TT,pp,order=1)
                eV = TT * R * (t*gt - pi*gp)

                e_[...] = eL + (eV-eL)*x_

        e = it.operands[0]

        # Convert the results
        scale = pyro.units.energy(from_units='kJ')
        scale = pyro.units.matter(scale,self.data['mw'],from_units='kg',exponent=-1)

        return e * scale



    def cp(self,T=None,p=None,x=None):
        """Constant pressure specific heat
        cp(T=None,p=None,x=None)

T   Temperature
p   pressure
x   quality

Accepts the temperature, pressure, or quality in any of the three 
possible combinations.  If none of the arguments are supplied, then 
def_T,def_p are used, and x is ignored.  If only x is supplied, p is 
presumed to be def_p, and T is calculated to be the saturation 
temperature.

Accepts unit_temperature
        unit_pressure
        dimensionless
Returns unit_energy / unit_matter / unit_temperature
"""
        def_T = False
        if T is None:
            T = pyro.config['def_T']
            def_T = True

        def_p = False
        if p is None:
            p = pyro.config['def_p']
            def_p = True
        pscale = pyro.units.pressure(to_units='bar')

        if x is None:
            x = -1.

        R = self.data['R']

        it = np.nditer((None,T,p,x), 
                op_flags=[['readwrite','allocate'],['readonly','copy'],
                    ['readonly','copy'],['readonly','copy']], op_dtypes='float')

        for cp_,T_,p_,x_ in it:
            TT = pyro.units.temperature_scale(T_, to_units='K')
            pp = p_ * pscale
            # If x is unspecified
            if x_<0.:

                # Detect the region
                r = self._region(TT,pp)
                # Case out the region indices
                if r==1:
                    pi,t,_,_,_,_,_,gtt = self._g1(TT,pp,order=2)
                    cp_[...] = -R * t*t*gtt
                elif r==2:
                    pi,t,_,_,_,_,_,gtt = self._g2(TT,pp,order=2)
                    cp_[...] = -R * t*t*gtt
                elif r==3:
                    n,t,_,fp,ft,fpp,fpt,ftt = self._f3(TT,pp)
                    temp = n*fp - n*t*fpt
                    temp = temp*temp/(2*n*fp + n*n*fpp)
                    cp_[...] = R * (-t*t*ftt + temp)
                elif r==5:
                    pi,t,_,_,_,_,_,gtt = self._g5(TT,pp,order=2)
                    cp_[...] = -R * t*t*gtt
                else:
                    raise pyro.utility.PMParamError('Invalid property combination T=%f K, p=%f bar'%(T_,pp))


            else:
                # If T was unspecified
                if def_T:
                    # Override the default T with the saturation T
                    TT = self._Ts(pp)
                # If pressure was unspecified, but temperature WAS
                elif def_p:
                    # Override the default p with the saturation p
                    pp = self._ps(TT)

                pi,t,_,_,_,_,_,gtt = self._g1(TT,pp,order=2)
                cpL = -R * t*t*gtt
                pi,t,_,_,_,_,_,gtt = self._g2(TT,pp,order=2)
                cpV = -R * t*t*gtt

                cp_[...] = cpL + (cpV-cpL)*x_

        cp = it.operands[0]

        # Convert the results
        scale = pyro.units.energy(from_units='kJ')
        scale = pyro.units.matter(scale,self.data['mw'],from_units='kg',exponent=-1)
        scale = pyro.units.temperature(scale,from_units='K',exponent=-1)

        return cp * scale

        
        
    def cv(self,T=None,p=None,x=None):
        """Constant volume specific heat (kJ/kg)
        cv(T=None,p=None,x=None)

T   Temperature
p   pressure
x   quality

Accepts the temperature, pressure, or quality in any of the three 
possible combinations.  If none of the arguments are supplied, then 
def_T,def_p are used, and x is ignored.  If only x is supplied, p is 
presumed to be def_p, and T is calculated to be the saturation 
temperature.

Accepts unit_temperature
        unit_pressure
        dimensionless
Returns unit_energy / unit_matter / unit_temperature
"""
        def_T = False
        if T is None:
            T = pyro.config['def_T']
            def_T = True

        def_p = False
        if p is None:
            p = pyro.config['def_p']
            def_p = True
        pscale = pyro.units.pressure(to_units='bar')

        if x is None:
            x = -1.

        R = self.data['R']

        it = np.nditer((None,T,p,x), 
                op_flags=[['readwrite','allocate'],['readonly','copy'],
                    ['readonly','copy'],['readonly','copy']], op_dtypes='float')

        for cv_,T_,p_,x_ in it:
            TT = pyro.units.temperature_scale(T_,to_units='K')
            pp = p_ * pscale

            # If x is unspecified
            if x_<0.:

                # Detect the region
                r = self._region(TT,pp)
                # Case out the region indices
                if r==1:
                    pi,t,_,gp,gt,gpp,gpt,gtt = self._g1(TT,pp,order=2)
                    temp = gp - t*gpt
                    temp = temp*temp/gpp
                    cv_[...] = R * (temp - t*t*gtt)
                elif r==2:
                    pi,t,_,gp,gt,gpp,gpt,gtt = self._g2(TT,pp,order=2)
                    temp = gp - t*gpt
                    temp = temp*temp/gpp
                    cv_[...] = R * (temp - t*t*gtt)
                elif r==3:
                    n,t,_,_,_,_,_,ftt = self._f3(TT,pp)
                    cv_[...] = -R * t*t*ftt
                elif r==5:
                    pi,t,_,gp,gt,gpp,gpt,gtt = self._g5(TT,pp,order=2)
                    temp = gp - t*gpt
                    temp = temp*temp/gpp
                    cv_[...] = R * (temp - t*t*gtt)
                else:
                    raise pyro.utility.PMParamError('Invalid property combination T=%f K, p=%f bar'%(T_,pp))


            else:
                # If T was unspecified
                if def_T:
                    # Override the default T with the saturation T
                    TT = self._Ts(pp)
                # If pressure was unspecified, but temperature WAS
                elif def_p:
                    # Override the default p with the saturation p
                    pp = self._ps(TT)

                pi,t,_,gp,gt,gpp,gpt,gtt = self._g1(TT,pp,order=2)
                temp = gp - t*gpt
                temp = temp*temp/gpp
                cvL = R * (temp - t*t*gtt)
                pi,t,_,gp,gt,gpp,gpt,gtt = self._g2(TT,pp,order=2)
                temp = gp - t*gpt
                temp = temp*temp/gpp
                cvV = R * (temp - t*t*gtt)

                cv_[...] = cvL + (cvV-cvL)*x_

        cv = it.operands[0]

        # Convert the results
        scale = pyro.units.energy(from_units='kJ')
        scale = pyro.units.matter(scale,self.data['mw'],from_units='kg',exponent=-1)
        scale = pyro.units.temperature(scale,from_units='K',exponent=-1)

        return cv * scale


    def gam(self,T=None,p=None,x=None):
        """Specific heat ratio
    gam(T=None,p=None,x=None)

T   Temperature
p   pressure
x   quality

Accepts the temperature, pressure, or quality in any of the three 
possible combinations.  If none of the arguments are supplied, then 
def_T,def_p are used, and x is ignored.  If only x is supplied, p is 
presumed to be def_p, and T is calculated to be the saturation 
temperature.

Accepts unit_temperature
        unit_pressure
        dimensionless
Returns unit_energy / unit_matter / unit_temperature
"""
        return cp(T,p,x) / cv(T,p,x)


    def mw(self,T=None,p=None,x=None):
        """Molecular weight (kg/kmol)
    mw(T=None,p=None,x=None)

T   Temperature
p   pressure
x   quality

Accepts the temperature, pressure, or quality, but ignores all values.

Accepts unit_temperature
        unit_pressure
        dimensionless
Returns unit_mass / unit_mol
"""
        mw = pyro.units.mass(self.data['mw'],from_units='g')
        mw = pyro.units.molar(mw,from_units='mol',exponent=-1)
        return mw
        



    def T_h(self,h,p=None,quality=False):
        """Temperature calculated from enthalpy
    T = T_h(h)
        or
    T = T_h(h,p)
        or
    T,x = T_h(h,p=None,quality=True)

Calculates temperature from the enthalpy and pressure.

When quality is set to True, T_h also returns the quality of saturated 
conditions.  Non-saturated conditions are assigned quality -1.

Accepts unit_energy / unit_matter
        unit_pressure
        (boolean)
Returns unit_temperature
        (dimensionless)
"""

        # First, figure out what region we're in.  We need to use a custom 
        # region-finding routine since we're in h,p coordinates and not T,p
        # Comments describe the process relative to the inverse region diagram 
        # in the IF-97 report on page 21 (Figure 2)
        T13 = 623.15    # vertical boundary between regions 1 and 3
        p2ab = 40.      # horizontal boundary between regions 2a and 2b
        p2c = 45.257578905948   # the bottom of the 2c region
        p3 = 165.292    # The bottom of region 3
        pmax = 1000.    # The top of the IF-97 region
        p5max = 500.       # The top of region 5
        T25 = 1073.15    # Left edge of region 5
        T5max = 2273.15    # Right edge of region 5
        R = self.data['R']

        # Prepare h and p
        if p is None:
            p = pyro.config['def_p']
        if not isinstance(p,np.ndarray):
            p = np.array(p)
        pscale = pyro.units.pressure(to_units='bar')

        scale = pyro.units.energy(to_units='kJ')
        scale = pyro.units.matter(scale,self.data['mw'],to_units='kg')
		
        it = np.nditer((None,None,h,p),
                op_flags=[['readwrite','allocate'],['readwrite','allocate'],
					['readonly','copy'], ['readonly','copy']],
					op_dtypes='float')

        for T_,x_,h_,p_ in it:
            T_[...] = -1.
            x_[...] = -1.
            pp = p_ * pscale
            # Scale the enthalpy to kJ/kg
            hh = h_ * scale
            if p_ < p3:
                # Calculate the saturation temperature and enthalpies
                Ts = self._Ts(p_)
                pi,t,_,_,gt,_,_,_ = self._g1(Ts,p_,order=1)
                hL = R * Ts * t * gt
                pi,t,_,_,gt,_,_,_ = self._g2(Ts,p_,order=1)
                hV = R * Ts * t * gt
                # If h is below the liquid enthalpy, use region 1
                if hh<hL:
                    T_[...] = self._th1(h=hh,p=p_)
                # If h is below the vapor enthalpy, this is a saturated mixture
                elif hh<hV:
                    T_[...] = Ts
                    x_[...] = (hh-hL)/(hV-hL)
                # If p is below the a-b boundary
                elif p_<p2ab:
                    # Calculate the enthalpy at the 2-5 boarder
                    pi,t,g,gp,gt,_,_,_ = self._g5(T25,p_,order=1)
                    h25 = R*T25 * t * gt
                    # If h is below the 2-5 boundary, this is region 2a
                    if hh<h25:
                        T_[...] = self._th2a(h=hh,p=p_)
                    else:
                        # Calculate the enthalpy at the upper bound of r5
                        pi,t,g,gp,gt,_,_,_ = self._g5(T5max,p_,order=1)
                        h5 = R*T5max * t * gt
                        if hh>h5:
                            raise pyro.utility.PMParamError(
                            'Steam T_h(): the state is not in the IF-97 domain.')
                        Tinit = T25 + (T5max-T25)*(hh-h25)/(h5-h25)
                        T_[...] = self._th5(h=hh, p=p_, Tinit=Tinit)
                # If h is below (left of) the b-c boarder
                elif p_>p2c and hh<self._b2bc(p=p_):
                    # Region 2c
                    T_[...] = self._th2c(h=hh, p=p_)
                else:
                    # Calculate the enthalpy at the 2-5 boarder
                    pi,t,g,gp,gt,_,_,_ = self._g5(T25,p_,order=1)
                    h25 = R*T25 * t * gt
                    # All that's left is either 2b or 5
                    if hh<h25:
                        # Region 2b
                        T_[...] = self._th2b(h=hh, p=p_)
                    else:
                        # Calculate the enthalpy at the upper bound of r5
                        pi,t,g,gp,gt,_,_,_ = self._g5(T5max,p_,order=1)
                        h5 = R*T5max * t * gt
                        if hh>h5:
                            raise pyro.utility.PMParamError(
                            'Steam T_h(): the state is not in the IF-97 domain.')
                        Tinit = T25 + (T5max-T25)*(hh-h25)/(h5-h25)
                        T_[...] = self._th5(h=hh, p=p_, Tinit=Tinit)
            elif p_ <= pmax:
                # Calculate the enthalpy at the 1-3 boarder
                pi,t,g,gp,gt,_,_,_ = self._g1(T13,p_,order=1)
                h13 = R*T13 * t * gt
                if hh<h13:
                    # Region 1
                    T_[...] = self._th1(h=hh,p=p_)
                else:
                    # Calculate T and h at the 2-3 boarder
                    T23 = self._b23(p=p_)
                    h23 = self.h(T=T23,p=p_)
                    pi,t,g,gp,gt,_,_,_ = self._g2(T23,p_,order=1)
                    h23 = R*T23 * t * gt
                    if hh<h23:
                        # Region 3
                        Tinit = T13 + (T23-T13)*(hh-h13)/(h23-h13)
                        T_[...] = self._th3(h=hh, p=p_, Tinit=Tinit)
                    elif hh<self._b2bc(p=p_):
                        # Region 2c
                        T_[...] = self._th2c(h=hh,p=p_)
                    else:
                        # Calculate the enthalpy at the 2-5 boarder
                        pi,t,g,gp,gt,_,_,_ = self._g5(T25,p_,order=1)
                        h25 = R*T25 * t * gt
                        if hh<h25:
                            # Region 2b
                            T_[...] = self._th2b(h=hh,p=p_)
                        elif p_<p5max:
                            # Calculate the enthalpy at the upper bound of r5
                            pi,t,g,gp,gt,_,_,_ = self._g5(T5max,p_,order=1)
                            h5 = R*T5max * t * gt
                            if hh>h5:
                                raise pyro.utility.PMParamError(
                                'Steam T_h(): the state is not in the IF-97 ' +
                                'domain.')
                            Tinit = T25 + (T5max-T25)*(hh-h25)/(h5-h25)
                            T_[...] = self._th5(h=hh, p=p_, Tinit=Tinit)
                        else:
                            raise pyro.utility.PMParamError(
                            'Steam T_h(): the state is not in the IF-97 domain.')
            else:
                raise pyro.utility.PMParamError(
                'Steam T_h(): pressure is above the IF-97 maximum (1000bar)')

            T_[...] = pyro.units.temperature_scale(T_,from_units='K')
  
        if quality:
            return it.operands[0:2]
        return it.operands[0]



    def T_s(self,s,p=None,quality=False):
        """Temperature calculated from entropy
    T = T_s(s)
        or
    T = T_s(s,p)
        or
    T,x = T_s(s,p=None,quality=False)

Calculates temperature from the entropy and pressure.

When quality is set to True, T_s also returns the quality of saturated 
conditions.  Non-saturated conditions are assigned quality -1.

Accepts unit_energy / unit_matter / unit_temperature
        unit_pressure
        (boolean)
Returns unit_temperature
        (dimensionless)
"""
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

        R = self.data['R']

        # Prepare s and p
        if p is None:
            p = pyro.config['def_p']
        if not isinstance(p,np.ndarray):
            p = np.array(p)
        p = pyro.units.pressure(p, to_units='bar')

        scale = pyro.units.energy(to_units='kJ')
        scale = pyro.units.matter(scale,self.data['mw'],to_units='kg')
        scale = pyro.units.temperature(scale,to_units='K')

        it = np.nditer((None,None,s,p),
		            op_flags=[['readwrite','allocate'],['readwrite','allocate'],
					['readonly','copy'], ['readonly','copy']],
					op_dtypes='float')

        for T_,x_,s_,p_ in it:
            T_[...] = -1.
            x_[...] = -1.
            # Scale the entropy to kJ/kg/K
            ss = s_ * scale
            if p_ < p3:
                Ts = self.Ts(p=p_)
                pi,t,g,_,gt,_,_,_ = self._g1(Ts,p_,order=1)
                sL = R * (t*gt - g)
                pi,t,g,_,gt,_,_,_ = self._g2(Ts,p_,order=1)
                sV = R * (t*gt - g)
                if ss<sL:
                    # Region 1
                    T_[...] = self._ts1(s=ss,p=p_)
                elif ss<sV:
                    # Saturation
                    T_[...] = Ts
                    x_[...] = (ss-sL)/(sV-sL)
                elif p_<p2ab:
                    s25 = self.s(T=T25,p=p_)
                    if ss<s25:
                        # Region 2a
                        T_[...] = self._ts2a(s=ss,p=p_)
                    else:
                        # Region 5
                        s5 = self.s(T=T5max,p=p_)
                        if ss>s5:
                            raise pyro.utility.PMParamError(
                            'Steam T_h(): the state is not in the IF-97 domain.')
                        Tinit = T25 + (T5max-T25)*(ss-s25)/(s5-s25)
                        T_[...] = self._ts5(s=ss, p=p_, Tinit=Tinit)
                elif ss<s2bc:
                    # Region 2c
                    T_[...] = self._ts2c(s=ss, p=p_)
                else:
                    s25 = self.s(T=T25,p=p_)
                    if ss<s25:
                        # Region 2b
                        T_[...] = self._ts2b(s=ss, p=p_)
                    else:
                        # Region 5
                        s5 = self.s(T=T5max,p=p_)
                        if ss>s5:
                            raise pyro.utility.PMParamError(
                            'Steam T_h(): the state is not in the IF-97 domain.')
                        Tinit = T25 + (T5max-T25)*(ss-s25)/(s5-s25)
                        T_[...] = self._ts5(s=ss, p=p_, Tinit=Tinit)
            elif p_ <= pmax:
                s13 = self.s(T=T13, p=p_)
                if ss<s13:
                    # Region 1
                    T_[...] = self._ts1(s=ss,p=p_)
                else:
                    T23 = self._b23(p=p_)
                    s23 = self.s(T=T23,p=p_)
                    if ss<s23:
                        # Region 3
                        Tinit = T13 + (T23-T13)*(ss-s13)/(s23-s13)
                        T_[...] = self._ts3(s=ss, p=p_, Tinit=Tinit)
                    elif ss<s2bc:
                        # Region 2c
                        T_[...] = self._ts2c(s=ss,p=p_)
                    else:
                        s25 = self.s(T=T25,p=p_)
                        if ss<s25:
                            # Region 2b
                            T_[...] = self._ts2b(s=ss,p=p_)
                        elif p_<=p5max:
                            # Region 5
                            s5 = self.s(T=T5max,p=p_)
                            if ss>s5:
                                raise pyro.utility.PMParamError(
                                'Steam T_s(): the state is not in the IF-97 ' +
                                'domain.')
                            Tinit = T25 + (T5max-T25)*(ss-s25)/(s5-s25)
                            T_[...] = self._ts5(s=ss, p=p_, Tinit=Tinit)
                        else:
                            raise pyro.utility.PMParamError(
                            '*Steam T_s(): the state is not in the IF-97 domain.')
            else:
                raise pyro.utility.PMParamError(
                'Steam T_h(): pressure is above the IF-97 maximum (1000bar)')

            T_[...] = pyro.units.temperature_scale(T_,from_units='K')
  
        if quality:
            return it.operands[0:2]
        return it.operands[0]
