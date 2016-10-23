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
    __version__='0.0'
    with open(install_init,'r') as init_file:
        for line in init_file:
            if line[:len(lookfor)] == lookfor:
                exec(line)
                break
    return __version__


setup(
    name='PYroMat',
    version=get_version(),
    description="Thermodynamic property calculations for Python.",
    long_description="Provides a Python interface for thermo-physical properties of a wide range of species. PYro is distinguished from other excellent thermodynamic resources (like Cantera) by its emphasis on being application agnostic. To do that, it provides base classes and data within a highly flexible framework that allows users to define their own species, new properties, and even entirely new data classes.",
    author="Chris Martin",
    author_email="crm28@psu.edu",
    url="https://chmarti1.github.io/PYroMat/",
    packages=['pyro'],
    package_dir={'pyro':install_from},
    package_data={'pyro':['registry/*.py','data/*.hpd']},
    requires=['numpy','json','distutils']
    )
