#!/usr/bin/env python
"""
Installation file for PYro
Chris Martin (c) 2016,2018
"""

from setuptools import setup
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
    description="Thermodynamic properties in Python.",
    long_description=\
"""PYroMat provides a Python interface for thermo-physical properties of a 
wide range of species. Visit the PYroMat homepage for more documentation, 
installation options, and more (https://chmarti1.github.io/PYroMat/).""",
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
    package_data={'':['data/mp/*.hpd', 'data/ig/*.hpd','data/ig2/*.hpd','registry/*.py']},
    data_files=[('.',[  'README.md', 'CHANGELOG.md', 'LICENSE.txt',
                        'test.py'])],
    install_requires=['numpy>=1.7'],
    provides=['pyromat']
    )
