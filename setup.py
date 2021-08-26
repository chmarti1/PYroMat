#!/usr/bin/env python
"""
Installation file for PYro
Chris Martin (c) 2016,2018
"""

import setuptools
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

# Load the README for the long description
with open('./README.md', 'r') as ff:
    readme = ff.read()

setuptools.setup(
    name='PYroMat',
    version=get_version(),
    description="Thermodynamic properties in Python.",
    long_description=readme,
    long_description_content_type = "text/markdown",
    author="Chris Martin",
    author_email="crm28@psu.edu",
    url="http://pyromat.org",
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
    packages=setuptools.find_packages(where=install_from),
    package_dir={'':install_from},
    py_modules=['dat', 'reg', 'utility', 'units'],
    package_data={'':['config.py', 'data/mp/*.hpd', 'data/ig/*.hpd','data/ig2/*.hpd','registry/*.py','aps/*.py']},
    license_files=['LICENSE.txt'],
    install_requires=['numpy>=1.7'],
    provides=['pyromat']
    )
