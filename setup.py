#!/usr/bin/env python
"""
Installation file for PYro
Chris Martin (c) 2016,2017
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
"""PYroMat provides a Python interface for thermo-physical properties of a 
wide range of species. PYroMat is distinguished from other excellent 
thermodynamic resources in a number of ways;
    (1) It takes full advantage of Python's unique features to 
        streamline the user interface (no long hard-to-remember
        function calls)
    (2) There is no attempt to standardize the data at the back-end,
        but there is one standard interface for a wide range of data
        types and sources.
    (3) Data sources are cited with academic rigor.
    (4) The standard package includes tools for users to generate their
        own data and data classes.
The result is a system that is easy to use from the command line, but
with the performance needed for advanced modeling.""",
    author="Chris Martin",
    author_email="crm28@psu.edu",
    url="https://chmarti1.github.io/PYroMat/",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Education',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Manufacturing',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Natural Language :: English',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3'],
    license = 'GNU General Public License v3 (GPLv3)',
    keywords = 'thermodynamic properties',
    packages=['pyromat'],
    package_dir={'pyromat':install_from},
    package_data={'pyromat':['registry/*.py','data/*.hpd']},
    data_files=[('.',['README.md','CHANGELOG.md','test.py'])],
    requires=['numpy','json','distutils'],
    provides=['pyromat']
    )
