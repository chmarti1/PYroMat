# This script validates all IG class data against their embedded 

import pyromat as pyro
import os, sys
import numpy as np
import time

logdir = 'log'

def runargtest(outfile, species, args, reference, error_threshold = .001):
    """Test properties for a species object given flexible arguments
    runargtest(outfile, species, args, reference)

outfile
    An open writable file to which to stream the results summary text.

species
    The pyromat object to test or the string id of the object to get

args
    A dictionary containing keyword arguments to pass to the methods identified
by reference.

reference
    A dictionary with keys that must correspond to property methods belonging
to the species.  The corresponding value to each key are the expected results
from the Tp test.  Examples are in test.py.
"""
    if isinstance(species,str):
        test = pyro.get(species)
    elif issubclass(type(species),pyro.reg.__basedata__):
        test = species
    else:
        raise Exception("Species needs to be a PYroMat ID string or a PYroMat object\n")
    sys.stdout.write("Testing " + str(test) + "..")
    outfile.write("Validation of {:s}\n".format(str(test)))
    outfile.write("  Error failure threshold {:%}\n".format(error_threshold))
    error = False

    for argname in args:
        N = np.array(args[argname]).size
        break

    if N>1:
        linefmt = "{:>15s} = " + "{:<15.5e}"*N + "\n"
        errfmt = "{:>15s} = " + "{:<15.5%}"*N + "\n"
        for argname in args:
            outfile.write(linefmt.format(argname,*args[argname]))
        outfile.write(72*'-' + '\n')
        for thisparam in reference:
            # Evaluate the parameter under test
            val = getattr(test,thisparam)(**args)
            val_ref = reference[thisparam]
            val_error = np.abs((val - val_ref)/val_ref)
            outfile.write(linefmt.format(thisparam,*val))
            outfile.write(linefmt.format('ref',*val_ref))
            outfile.write(errfmt.format('error',*val_error))
            if np.any(val_error > error_threshold):
                error = True
                outfile.write("    [FAILED]\n")
            else:
                outfile.write("    [passed]\n")
    else:
        linefmt = "{:>15s} = {:<15.5e}\n"
        errfmt = "{:>15s} = {:<15.5%}\n"
        for argname in args:
            outfile.write(linefmt.format(argname,float(args[argname])))
        for thisparam in reference:
            # Evaluate the parameter under test
            val = getattr(test,thisparam)(**args)
            val_ref = reference[thisparam]
            val_error = np.abs((val - val_ref)/val_ref)
            outfile.write(linefmt.format(thisparam,float(val)))
            outfile.write(linefmt.format('ref',float(val_ref)))
            outfile.write(errfmt.format('error',float(val_error)))
            if np.any(val_error > error_threshold):
                error = True
                outfile.write("    [FAILED]\n")
            else:
                outfile.write("    [passed]\n")
    outfile.write("\n")

    if error:
        sys.stdout.write("[FAILED]\n")
    else:
        sys.stdout.write("[passed]\n")
    return error


# Create the log directory if it doesn't already exist
if not os.path.isdir(logdir):
    os.mkdir(logdir)

logfile = os.path.join(logdir,'steam.log')
with open(logfile,'w+') as writeto:

    writeto.write("PYroMat Validation Report\n")
    writeto.write(time.strftime("%Y-%m-%d\n"))
    writeto.write("  PYroMat Version: {:s}\n".format(pyro.config['version']))
    writeto.write("  Installation: {:s}\n".format(pyro.config['install_dir']))
    writeto.write("  Found {:d} species\n".format(len(pyro.dat.data)))

    # Python version
    writeto.write("Python version: {:s}\n".format(sys.version.split()[0]))
    writeto.write("Numpy version: {:s}\n".format(np.version.version))
    # Operating system
    writeto.write("Running in a {:s} environment\n".format(os.name))
    if os.name=='posix':
        writeto.write("  {0:s} {2:s} {4:s}\n".format(*os.uname()))
    writeto.write('\n')

    # STEAM
    T = np.array([300., 300., 500.])
    p = np.array([30., 800., 30.])
    args = {'T':T, 'p':p}

    reference = {'cp':np.array([4.17301218, 4.01008987, 4.65580682]),
        'h':np.array([115.331273, 184.142828, 975.542239]),
        's':np.array([.392294792, .368563852, 2.58041912]),
        'e':np.array([112.324818, 106.448356, 971.934985]),
        'd':1./np.array([.100215168e-2, .971180894e-3, .120241800e-2])
    }

    writeto.write( "Steam validation values for region 1 found at \n" +
    "http://www.iapws.org/relguide/IF97-Rev.pdf\n" +
    "   Table 5 pg 9\n" )
    runargtest(writeto,'steam',args,reference)

    T = np.array([300., 700., 700.])
    p = np.array([.035, .035, 300.])
    args = {'T':T, 'p':p}

    # Reference values from IF-97
    reference = {'cp':np.array([1.91300162, 2.08141274, 10.3505092]),
        'h':np.array([.254991145e4, .333568375e4, .263149474e4]),
        's':np.array([.852238967e1, .101749996e2, .517540298e1]),
        'e':np.array([.241169160e4, .301262819e4, .246861076e4]),
        'd':1./np.array([.394913866e2, .923015898e2, .542946619e-2])
    }

    writeto.write( "Steam validation values for region 2 found at \n" +
    "http://www.iapws.org/relguide/IF97-Rev.pdf\n" +
    "   Table 15 pg 17\n" )
    runargtest(writeto,'steam',args,reference)

    T = np.array([650., 650., 750.])
    p = np.array([.255837018e3, .222930643e3, .783095639e3])
    args = {'T':T, 'p':p}

    # Reference values from IF-97
    reference = {'cp':np.array([.138935717e2, .446579342e2, .634165359e1]),
        'h':np.array([.186343019e4, .237512401e4, .225868845e4]),
        's':np.array([.405427273e1, .485438792e1, .446971906e1]),
        'e':np.array([.181226279e4, .226365868e4, .210206932e4]),
        'd':np.array([500., 200., 500.])
    }

    writeto.write( "Steam validation values for region 3 found at \n" +
    "http://www.iapws.org/relguide/IF97-Rev.pdf\n" +
    "   Table 33 pg 32\n" )
    runargtest(writeto,'steam',args,reference)

    T = np.array([1500., 1500., 2000.])
    p = np.array([5., 300., 300.])
    args = {'T':T, 'p':p}

    # Reference values from IF-97
    reference = {'cp':np.array([.261609445e1, .272724317e1, .288569882e1]),
        'h':np.array([.521976855e4, .516723514e4, .657122604e4]),
        's':np.array([.965408875e1, .772970133e1, .853640623e1]),
        'e':np.array([.452749310e4, .447495124e4, .563707038e4]),
        'd':1./np.array([.138455090e1, .230761299e-1, .311385219e-1])
    }

    writeto.write( "Steam validation values for region 5 found at \n" +
    "http://www.iapws.org/relguide/IF97-Rev.pdf\n" +
    "   Table 42 pg 40\n" )
    runargtest(writeto,'steam',args,reference)


    test = pyro.get('steam')
    T = np.array([300., 500., 600.])
    args = {'T':T}
    reference = {'ps':np.array([.353658941e-1, .263889776e2, .123443146e3])}

    writeto.write( "Steam validation values for saturation (region 4) found at \n" +
    "http://www.iapws.org/relguide/IF97-Rev.pdf\n" +
    "   Table 35 pg 34 and Table 36 pg 36\n" )
    runargtest(writeto,'steam',args,reference)

    p = np.array([1., 10., 100.])
    args = {'p':p}
    reference = {'Ts': np.array([.372755919e3, .453035632e3, .584149488e3])}
    runargtest(writeto,'steam',args,reference)


    # Inverse relations
    writeto.write( "Steam validation values for the inverse relations\n" +
    "   Values are borrowed from the above tests, but run in reverse.\n" )

    T = []
    p = []
    s = []
    h = []
    # Values from Region 1
    T += [300., 300., 500.]
    p += [30., 800., 30.]
    s += [.392294792, .368563852, 2.58041912]
    h += [115.331273, 184.142828, 975.542239]
    # Values from Region 2
    T += [300., 700., 700.]
    p += [.035, .035, 300.]
    s += [.852238967e1, .101749996e2, .517540298e1]
    h += [.254991145e4, .333568375e4, .263149474e4]
    # Values from Region 3
    T += [650., 650., 750.]
    p += [.255837018e3, .222930643e3, .783095639e3]
    s += [.405427273e1, .485438792e1, .446971906e1]
    h += [.186343019e4, .237512401e4, .225868845e4]
    # Values from Region 5
    T += [1500., 1500., 2000.]
    p += [5., 300., 300.]
    s += [.965408875e1, .772970133e1, .853640623e1]
    h += [.521976855e4, .516723514e4, .657122604e4]

    reference = {'T_h':T}
    args = {'h':h, 'p':p}
    runargtest(writeto,'steam',args,reference)

    reference = {'T_s':T}
    args = {'s':s, 'p':p}
    runargtest(writeto,'steam',args,reference)

