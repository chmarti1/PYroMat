import pyromat as pyro
import os, sys
import numpy as np

with open('test.log','w+') as writeto:

    writeto.write("PYroMat Validation Report\n")
    writeto.write("  Version: {:s}\n".format(pyro.config['version']))
    writeto.write("  Installation: {:s}\n".format(pyro.config['install_dir']))
    writeto.write("  Found {:d} species\n".format(len(pyro.dat.data)))

    # Python version
    writeto.write("Python version: {:s}\n".format(sys.version.split()[0]))
    # Operating system
    writeto.write("Running in a {:s} environment\n".format(os.name))
    if os.name=='posix':
        writeto.write("  {0:s} {2:s} {4:s}\n".format(*os.uname()))
    writeto.write('\n')

    # Test the igfit class
    species = 'O2'
    print "Testing " + species
    test = pyro.get(species)
    T = 500.
    p = 20.
    # Reference values from the NIST tables
    mw_ref = 31.9988
    R_ref = 8.314 / mw_ref  # kJ/kg
    cp_ref = 31.091 / mw_ref
    h_ref = 6.084 * 1000. / mw_ref
    s_ref = 220.693 / mw_ref - R_ref * np.log(p/1.)
    d_ref = p*1e2 / R_ref / T   # 1e5/1e3->1e2 conversion from bar->Pa and kJ->J

    mw = test.mw()
    R = test.R()
    cp = test.cp(T=T,p=p)
    h = test.h(T=T,p=p)
    s = test.s(T=T,p=p)
    d = test.d(T=T,p=p)

    writeto.write( 
    "Using " + species + " to test the IGFIT class\n" +
    "  Reference values found at (http://kinetics.nist.gov/janaf/html/O-029.html)\n" +
    "  p = {:.5e} bar\n".format(p) +
    "  T = {:.5e} K\n".format(T) +
    "  mw[kg/kmol] pyro={:.5e} ref={:.5e} err={:.3%}\n".format(mw,mw_ref,(mw-mw_ref)/mw_ref) +
    "  R[kJ/kg/K]  pyro={:.5e} ref={:.5e} err={:.3%}\n".format(R,R_ref,(R-R_ref)/R_ref) +
    "  cp[kg/kg/K] pyro={:.5e} ref={:.5e} err={:.3%}\n".format(cp,cp_ref,(cp-cp_ref)/cp_ref) +
    "  h[kJ/kg]    pyro={:.5e} ref={:.5e} err={:.3%}\n".format(h,h_ref,(h-h_ref)/h_ref) + 
    "  s[kJ/kg/K]  pyro={:.5e} ref={:.5e} err={:.3%}\n".format(s,s_ref,(s-s_ref)/s_ref) + 
    "  d[kg/m**3]  pyro={:.5e} ref={:.5e} err={:.3%}\n".format(d,d_ref,(d-d_ref)/d_ref) +
    "\n"
    )

    # Test the igtab class
    species = 'I2'
    print "Testing " + species
    test = pyro.get(species)
    T = 600.
    p = 20.
    # Reference values from the NIST tables
    mw_ref = 253.80894
    R_ref = 8.314 / mw_ref  # kJ/kg/K
    cp_ref = 37.612 / mw_ref
    h_ref = 73.690 * 1000. / mw_ref
    s_ref = 286.764 / mw_ref - R_ref * np.log(p/1.)
    d_ref = p*1e2 / R_ref / T   # 1e5/1e3->1e2 conversion from bar->Pa and kJ->J

    mw = test.mw()
    R = test.R()
    cp = test.cp(T=T,p=p)
    h = test.h(T=T,p=p)
    s = test.s(T=T,p=p)
    d = test.d(T=T,p=p)

    writeto.write( 
    "Using " + species + " to test the IGTAB class\n" +
    "  Reference values found at (http://kinetics.nist.gov/janaf/html/Xe-001.html)\n" +
    "  p = {:.5e} bar\n".format(p) +
    "  T = {:.5e} K\n".format(T) +
    "  mw[kg/kmol] pyro={:.5e} ref={:.5e} err={:.3%}\n".format(mw,mw_ref,(mw-mw_ref)/mw_ref) +
    "  R[kJ/kg/K]  pyro={:.5e} ref={:.5e} err={:.3%}\n".format(R,R_ref,(R-R_ref)/R_ref) +
    "  cp[kg/kg/K] pyro={:.5e} ref={:.5e} err={:.3%}\n".format(cp,cp_ref,(cp-cp_ref)/cp_ref) +
    "  h[kJ/kg]    pyro={:.5e} ref={:.5e} err={:.3%}\n".format(h,h_ref,(h-h_ref)/h_ref) + 
    "  s[kJ/kg/K]  pyro={:.5e} ref={:.5e} err={:.3%}\n".format(s,s_ref,(s-s_ref)/s_ref) + 
    "  d[kg/m**3]  pyro={:.5e} ref={:.5e} err={:.3%}\n".format(d,d_ref,(d-d_ref)/d_ref) +
    "\n"
    )


    # Test the if97 class
    species = 'steam'
    print "Testing " + species
    test = pyro.get(species)

    T = np.array([300., 300., 500.])
    p = np.array([30., 800., 30.])

    # Reference values from IF-97
    mw_ref = 18.0141
    cp_ref = np.array([4.17301218, 4.01008987, 4.65580682])
    h_ref = np.array([115.331273, 184.142828, 975.542239])
    s_ref = np.array([.392294792, .368563852, 2.58041912])
    e_ref = np.array([112.324818, 106.448356, 971.934985])
    d_ref = 1./np.array([.100215168e-2, .971180894e-3, .120241800e-2])

    mw = test.mw()
    cp = test.cp(T=T,p=p)
    h = test.h(T=T,p=p)
    s = test.s(T=T,p=p)
    e = test.e(T=T,p=p)
    d = test.d(T=T,p=p)

    writeto.write( 
    "Using " + species + " to test the IF97 class in region 1\n" +
    "  Reference values found at (http://www.iapws.org/relguide/IF97-Rev.pdf)\n" +
    "       Table 5 pg 9\n" +
    "  mw[kg/kmol] pyro={:.5e} ref={:.5e} err%={:.5f}\n".format(mw,mw_ref,(mw-mw_ref)/mw_ref) +
    "            p[bar] =" + ("{:13.5e}"*3).format(*p) + "\n" +
    "              T[K] =" + ("{:13.5e}"*3).format(*T) + "\n" +
    "  cp[kg/kg/K] pyro =" + ("{:13.5e}"*3).format(*cp) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*cp_ref) + "\n" +
    "               err =" + ("{:13.3%}"*3).format(*((cp-cp_ref)/cp_ref)) + "\n" + 
    "     h[kJ/kg] pyro =" + ("{:13.5e}"*3).format(*h) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*h_ref) + "\n" +
    "               err =" + ("{:13.3%}"*3).format(*((h-h_ref)/h_ref)) + "\n" + 
    "   s[kJ/kg/K] pyro =" + ("{:13.5e}"*3).format(*s) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*s_ref) + "\n" +
    "          err frac =" + ("{:13.3%}"*3).format(*((s-s_ref)/s_ref)) + "\n" + 
    "     e[kJ/kg] pyro =" + ("{:13.5e}"*3).format(*e) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*e_ref) + "\n" +
    "          err frac =" + ("{:13.3%}"*3).format(*((e-e_ref)/e_ref)) + "\n" + 
    "   d[kg/m**3] pyro =" + ("{:13.5e}"*3).format(*d) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*d_ref) + "\n" +
    "          err frac =" + ("{:13.3%}"*3).format(*((d-d_ref)/d_ref)) + "\n" +
    "\n")


    T = np.array([300., 700., 700.])
    p = np.array([.035, .035, 300.])

    # Reference values from IF-97
    mw_ref = 18.0141
    cp_ref = np.array([1.91300162, 2.08141274, 10.3505092])
    h_ref = np.array([.254991145e4, .333568375e4, .263149474e4])
    s_ref = np.array([.852238967e1, .101749996e2, .517540298e1])
    e_ref = np.array([.241169160e4, .301262819e4, .246861076e4])
    d_ref = 1./np.array([.100215168e-2, .971180894e-3, .120241800e-2])

    mw = test.mw()
    cp = test.cp(T=T,p=p)
    h = test.h(T=T,p=p)
    s = test.s(T=T,p=p)
    e = test.e(T=T,p=p)
    d = test.d(T=T,p=p)

    writeto.write( 
    "Using " + species + " to test the IF97 class in region 2\n" +
    "  Reference values found at (http://www.iapws.org/relguide/IF97-Rev.pdf)\n" +
    "       Table 15 pg 17\n" +
    "  mw[kg/kmol] pyro={:.5e} ref={:.5e} err%={:.5f}\n".format(mw,mw_ref,(mw-mw_ref)/mw_ref) +
    "            p[bar] =" + ("{:13.5e}"*3).format(*p) + "\n" +
    "              T[K] =" + ("{:13.5e}"*3).format(*T) + "\n" +
    "  cp[kg/kg/K] pyro =" + ("{:13.5e}"*3).format(*cp) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*cp_ref) + "\n" +
    "               err =" + ("{:13.5%}"*3).format(*((cp-cp_ref)/cp_ref)) + "\n" + 
    "     h[kJ/kg] pyro =" + ("{:13.5e}"*3).format(*h) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*h_ref) + "\n" +
    "               err =" + ("{:13.5%}"*3).format(*((h-h_ref)/h_ref)) + "\n" + 
    "   s[kJ/kg/K] pyro =" + ("{:13.5e}"*3).format(*s) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*s_ref) + "\n" +
    "               err =" + ("{:13.5%}"*3).format(*((s-s_ref)/s_ref)) + "\n" + 
    "     e[kJ/kg] pyro =" + ("{:13.5e}"*3).format(*e) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*e_ref) + "\n" +
    "               err =" + ("{:13.5%}"*3).format(*((e-e_ref)/e_ref)) + "\n" + 
    "   d[kg/m**3] pyro =" + ("{:13.5e}"*3).format(*e) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*e_ref) + "\n" +
    "               err =" + ("{:13.5%}"*3).format(*((e-e_ref)/e_ref)) + "\n" +
    "\n")


    T = np.array([650., 650., 750.])
    p = np.array([.255837018e3, .222930643e3, .783095639e3])

    # Reference values from IF-97
    mw_ref = 18.0141
    cp_ref = np.array([.138935717e2, .446579342e2, .634165359e1])
    h_ref = np.array([.186343019e4, .237512401e4, .225868845e4])
    s_ref = np.array([.405427273e1, .485438792e1, .446971906e1])
    e_ref = np.array([.181226279e4, .226365868e4, .210206932e4])
    d_ref = np.array([500., 200., 500.])

    mw = test.mw()
    cp = test.cp(T=T,p=p)
    h = test.h(T=T,p=p)
    s = test.s(T=T,p=p)
    e = test.e(T=T,p=p)
    d = test.d(T=T,p=p)

    writeto.write( 
    "Using " + species + " to test the IF97 class in region 3\n" +
    "  Reference values found at (http://www.iapws.org/relguide/IF97-Rev.pdf)\n" +
    "       Table 33 pg 32\n" +
    "  mw[kg/kmol] pyro={:.5e} ref={:.5e} err%={:.5f}\n".format(mw,mw_ref,(mw-mw_ref)/mw_ref) +
    "            p[bar] =" + ("{:13.5e}"*3).format(*p) + "\n" +
    "              T[K] =" + ("{:13.5e}"*3).format(*T) + "\n" +
    "  cp[kg/kg/K] pyro =" + ("{:13.5e}"*3).format(*cp) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*cp_ref) + "\n" +
    "               err =" + ("{:13.3%}"*3).format(*((cp-cp_ref)/cp_ref)) + "\n" + 
    "     h[kJ/kg] pyro =" + ("{:13.5e}"*3).format(*h) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*h_ref) + "\n" +
    "               err =" + ("{:13.3%}"*3).format(*((h-h_ref)/h_ref)) + "\n" + 
    "   s[kJ/kg/K] pyro =" + ("{:13.5e}"*3).format(*s) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*s_ref) + "\n" +
    "               err =" + ("{:13.3%}"*3).format(*((s-s_ref)/s_ref)) + "\n" + 
    "     e[kJ/kg] pyro =" + ("{:13.5e}"*3).format(*e) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*e_ref) + "\n" +
    "               err =" + ("{:13.3%}"*3).format(*((e-e_ref)/e_ref)) + "\n" + 
    "   d[kg/m**3] pyro =" + ("{:13.5e}"*3).format(*e) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*e_ref) + "\n" +
    "               err =" + ("{:13.3%}"*3).format(*((e-e_ref)/e_ref)) + "\n" +
    "\n")



    T = np.array([1500., 1500., 2000.])
    p = np.array([5., 300., 300.])

    # Reference values from IF-97
    mw_ref = 18.0141
    cp_ref = np.array([.261609445e1, .272724317e1, .288569882e1])
    h_ref = np.array([.521976855e4, .516723514e4, .657122604e4])
    s_ref = np.array([.965408875e1, .772970133e1, .853640623e1])
    e_ref = np.array([.452749310e4, .447495124e4, .563707038e4])
    d_ref = 1./np.array([.138455090e1, .230761299e-1, .311385219e-1])

    mw = test.mw()
    cp = test.cp(T=T,p=p)
    h = test.h(T=T,p=p)
    s = test.s(T=T,p=p)
    e = test.e(T=T,p=p)
    d = test.d(T=T,p=p)

    writeto.write( 
    "Using " + species + " to test the IF97 class in region 5\n" +
    "  Reference values found at (http://www.iapws.org/relguide/IF97-Rev.pdf)\n" +
    "       Table 15 pg 17\n" +
    "  mw[kg/kmol] pyro={:.5e} ref={:.5e} err%={:.5f}\n".format(mw,mw_ref,(mw-mw_ref)/mw_ref) +
    "            p[bar] =" + ("{:13.5e}"*3).format(*p) + "\n" +
    "              T[K] =" + ("{:13.5e}"*3).format(*T) + "\n" +
    "  cp[kg/kg/K] pyro =" + ("{:13.5e}"*3).format(*cp) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*cp_ref) + "\n" +
    "               err =" + ("{:13.3%}"*3).format(*((cp-cp_ref)/cp_ref)) + "\n" + 
    "     h[kJ/kg] pyro =" + ("{:13.5e}"*3).format(*h) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*h_ref) + "\n" +
    "               err =" + ("{:13.3%}"*3).format(*((h-h_ref)/h_ref)) + "\n" + 
    "   s[kJ/kg/K] pyro =" + ("{:13.5e}"*3).format(*s) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*s_ref) + "\n" +
    "               err =" + ("{:13.3%}"*3).format(*((s-s_ref)/s_ref)) + "\n" + 
    "     e[kJ/kg] pyro =" + ("{:13.5e}"*3).format(*e) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*e_ref) + "\n" +
    "               err =" + ("{:13.3%}"*3).format(*((e-e_ref)/e_ref)) + "\n" + 
    "   d[kg/m**3] pyro =" + ("{:13.5e}"*3).format(*e) + "\n" +
    "               ref =" + ("{:13.5e}"*3).format(*e_ref) + "\n" +
    "               err =" + ("{:13.3%}"*3).format(*((e-e_ref)/e_ref)) + "\n" +
    "\n")





    T = np.array([300., 500., 600.])
    ps = test.ps(T=T)
    ps_ref = np.array([.353658941e-1, .263889776e2, .123443146e3])
    
    p = np.array([1., 10., 100.])
    Ts = test.Ts(p=p)
    Ts_ref = np.array([.372755919e3, .453035632e3, .584149488e3])


    writeto.write(
    "Using " + species + " to test the saturation conditions for the IF97 class\n" +
    "  Reference values found at (http://www.iapws.org/relguide/IF97-Rev.pdf)\n" +
    "       Table 35 pg 34, Table 36 pg 36\n" +
    "       T[K] =" + ("{:13.5e}"*3).format(*T) + "\n" + 
    "    ps[bar] =" + ("{:13.5e}"*3).format(*ps) + "\n" +
    "        ref =" + ("{:13.5e}"*3).format(*ps_ref) + "\n" + 
    "        err =" + ("{:13.3%}"*3).format(*((ps-ps_ref)/ps_ref)) + "\n\n" +
    "     p[bar] =" + ("{:13.5e}"*3).format(*p) + "\n" +
    "      Ts[K] =" + ("{:13.5e}"*3).format(*Ts) + "\n" +
    "        ref =" + ("{:13.5e}"*3).format(*Ts_ref) + "\n" +
    "        err =" + ("{:13.3%}"*3).format(*((Ts-Ts_ref)/Ts_ref)) + "\n\n"
)
