# This script validates all IG class data against their embedded 

import pyromat as pyro
import os, sys
import numpy as np
import time

logdir = os.path.abspath('log')
masterlog = os.path.join(logdir,'ig.log')

# Create the log directory if it doesn't already exist
if not os.path.isdir(logdir):
    os.mkdir(logdir)

with open(masterlog, 'w+') as ff:
    for this in pyro.dat.data:
        if isinstance(this, pyro.reg.registry['ig']):
            with open(os.path.join(logdir,
