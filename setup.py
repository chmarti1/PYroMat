#!/usr/bin/env python
"""
Installation file for PYro
Chris Martin (c) 2016
"""

from distutils.core import setup
import os

# Where is the installation?
install_from = 'lib'
# Where is the master __init__ file?
install_init = os.path.join(install_from, '__init__.py')


# scans __init__.py for the __version__ string
def get_version():
    lookfor = '__version__'
    localtemp = {lookfor:'0.0'}
    with open(install_init,'r') as init_file:
        for line in init_file:
            if line[:len(lookfor)] == lookfor:
                exec(line, {}, localtemp)
                break
    return localtemp[lookfor]


setup(
    name='PYroMat',
    version=get_version(),
    description="Thermodynamic property calculations for Python.",
    long_description=\
"""Provides a Python interface for thermo-physical properties of a wide range of 
species. PYroMat is distinguished from other excellent thermodynamic resources 
in a number of ways; first among them that it is implemented entirely in Python. 
PYroMat combines a wide variety of data in incompatible formats in one seamsless
interface that strives to expose only what the user needs while keeping the nuts
and bolts to the background.""",
    author="Chris Martin",
    author_email="crm28@psu.edu",
    url="https://chmarti1.github.io/PYroMat/",
    packages=['pyromat'],
    package_dir={'pyromat':install_from},
    package_data={'pyromat':['registry/*.py','data/*.hpd']},
    requires=['numpy','json','distutils']
    )
