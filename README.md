# PYroMat

Thermodynamic tools for Python

By Chris Martin [crm28@psu.edu](mailto:crm28@psu.edu)

PYroMat is a flexible platform for conveniently working with thermodynamic data.  The expanding collection of substances includes data for the properties people need most, exposed in an intuitively designed object interface [Come read more.](http://www.pyromat.org)

## Installation from the Python Package Index
If you have pip installed, you can install PYroMat with a single command.
```
$ pip install pyromat 
```
If you are upgrading, you can always use
```
$ pip install pyromat --upgrade
```

## Getting started
```python
>>> import pyromat as pm
>>> O2 = pm.get('ig.O2')
>>> h = O2.h(492,1.01)  # enthalpy at 492K, 1.01bar
>>> pm.info('O2')     # where did these data come from?
>>> pm.config['unit_pressure'] = 'psi'  # Don't like working in bar?
```

Happy calculating!

## License
Copyright (c) 2015-2021 Christopher R. Martin

PYroMat is released under the GNU [General Public License v3.0](http://www.gnu.org/licenses/gpl-3.0.en.html).

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

