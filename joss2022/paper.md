---
title: 'PYroMat: A Python package for thermodynamic properties'
tags:
  - Python
  - thermodynamics
  - properties
  - ideal gas
  - multi-phase
authors:
  - name: Christopher Martin
    orcid: 0000-0002-7129-9367
    corresponding: true
    affiliation: 1
  - name: Joseph Ranalli
    orcid: 0000-0002-8184-9895
    affiliation: 2
  - name: Jacob Moore
    orcid: 0000-0001-7513-5979
    affiliation: 3
affiliations:
 - name: Penn State Altoona, Altoona, PA, USA
   index: 1
 - name: Penn State Hazleton, Hazleton, PA, USA
   index: 2
 - name: Penn State Mont Alto, Mont Alto, PA, USA 
   index: 3
date: 13 August 2022
bibliography: paper.bib

---

# Summary

Lookup of thermodynamic properties of substances is a common task in numerous engineering disciplines and is also an important skill taught to students in several common engineering disciplines. While properties were traditionally found using tables of thermodynamic data, availability of modern software has essentially eliminated the need for laborious use of tables, and has created opportunities to support quick iteration and revision of calculations. As thermodyanmic data are published in a variety of formats, it is necessary for software to provide a standardized interface to operate on a number of different data source formats. Additionally, thermodynamic software needs to be numerically efficient and robust, to deal with the diverse computational methodologies required to produce property calculations. PYroMat seeks to create a tool that meets these needs and can grow with users at different levels of expertise. It offers a low barrier to entry and is easy enough for students to use, but also powerful and flexible enough to serve the needs of engineering professionals and scientific researchers.

# Statement of need

As far back as the birth of the industrial revolution, engineers and scientists have needed precise calculations for the thermodynamic properties of fluids to predict the behaviors of systems of global importance.  Today, people working in aerospace propulsion, electrical power generation, plasma physics, refrigeration, building heating and cooling, combusion, and countless other fields still rely on decades of excellent data for these calculations.  It is exceedingly rare to find either the original data or software that performs these important calculations in the public domain.

The current industry standard, NIST's `REFPROP` [@Lemmon:2018], is neither free nor open.  Another excellent alternative is `coolprop`[@Bell:2014], which is an NSF-funded flexible interface for an impressive variety of languages, which focuses heavily on multi-phase substances.  There are also other less widely embraced alternatives, but `PYroMat` distinguishes itself because:
a) every aspect of it (and its dependencies) is open source, and 
b) it provides a standard Pythonic interface simultaneously simple enough for novices and powerful enough for professionals. 
Several conference publications have already described early versions of the software
[@Martin:2016], and its application in the undergraduate thermodyanmics 
classroom [@Martin:2017; @Ranalli:2019]. 

# Features

As of version 2.2.1, `PYroMat` provided property models for nearly 1,000 substances, and there are plans to continue expanding.  Each substance (whether pure or mixture) is represented by an instance of a `PYroMat` class.  Its methods are responsible for calculating its properties using a standard interface, so users need not be aware of the nuances of the back-end models.  For example, this segement of code retrieves the diatomic nitrogen model and calculates its enthalpy at 372.15 Kelvin and 1.4 bar.  Then, it returns the critical point (Kelvin, bar) of nitrogen.
```python
>>> import pyromat as pm
>>> n2 = pm.get('mp.N2')
>>> n2.h(T=372.15, p=1.4)
array([386.32246521])
>>> n2.critical()
(126.192, 33.958000000000006)
```

`PYroMat` can be configured to work in virtually any system of units in common use, and users can specify states with almost any pair of available properties.  The interface is well documented using Python's in-line help system, on the PYroMat website, [pyromat.org](http://pyromat.org), and in the regularly updated PYroMat User and Developer's Handbook [@pmhandbook].

Substances are organized into two collections: ideal gas (`ig`) and multi-phase (`mp`).  They are further organized by the specific underlying data model (class), so users always know what assumptions are implicit in the property calculations.  All of the classes implement a common flexible call signature, so the user can concentrate on the application instead of the nuances of the model.  As of version 2.2.1, `PYroMat` implements ideal gas data using the Shomate model (used by NIST) with the `ig` class, the time tested `NASA` polynomials with the `ig2` class, and ideal gas mixtures using the `igmix` class.  The so-called Span and Wagner multi-phase models are implemented with the `mp1` class.  There are plans to add more classes in later releases.

Users can conveniently browse for substances that interest them by name, chemical formula, atomic composition, molecular weight, InChIE or CAS identifier string, `PYroMat` class, or by collection.  For example, this search returns a table of all substances that contain two nitrogen atoms.

```python
>>> pm.info(contains={'N':2})
  PYroMat
Thermodynamic computational tools for Python
version: 2.2.1
---------------------------------------------------------------------------------------
 ID         : class : name                         : properties
---------------------------------------------------------------------------------------
 ig.C2H6N2  :  ig2  :                              : T p d v cp cv gam e h s mw R    
 ig.C2K2N2  :  ig2  : Potassium cyanide dimer      : T p d v cp cv gam e h s mw R    
 ig.C2N2    :  ig2  : Cyanogen                     : T p d v cp cv gam e h s mw R    
 ig.C2N2Na2 :  ig2  : Sodium cyanide               : T p d v cp cv gam e h s mw R    
 ig.C4N2    :  ig2  : 2-Butynedinitrile            : T p d v cp cv gam e h s mw R    
 ig.CN2     :  ig2  : CNN Radical                  : T p d v cp cv gam e h s mw R    
 ig.CN2_1   :  ig   : Amidogen, methanetetraylbis- : T p d v cp cv gam e h s mw R    
 ig.D2N2    :  ig2  : Diazene-d2, cis              : T p d v cp cv gam e h s mw R    
 ig.F2N2    :  ig2  : Nitrogen fluoride, (E)-      : T p d v cp cv gam e h s mw R    
 ig.F2N2_1  :  ig   : (Z)-Difluorodiazene          : T p d v cp cv gam e h s mw R    
 ig.F4N2    :  ig2  : Tetrafluorohydrazine         : T p d v cp cv gam e h s mw R    
 ig.H2N2    :  ig2  : (Z)-Diazene                  : T p d v cp cv gam e h s mw R    
 ig.H2N2O2  :  ig2  :                              : T p d v cp cv gam e h s mw R    
 ig.H4N2    :  ig2  : Hydrazine                    : T p d v cp cv gam e h s mw R    
 ig.N2      :  ig2  : Nitrogen                     : T p d v cp cv gam e h s mw R    
---------------------------------------------------------------------------------------
 ID         : class : name                         : properties
---------------------------------------------------------------------------------------
 ig.N2+     :  ig2  :                              : T p d v cp cv gam e h s mw R    
 ig.N2-     :  ig2  :                              : T p d v cp cv gam e h s mw R    
 ig.N2O     :  ig2  : Nitrous oxide                : T p d v cp cv gam e h s mw R    
 ig.N2O+    :  ig2  :                              : T p d v cp cv gam e h s mw R    
 ig.N2O3    :  ig2  : Dinitrogen trioxide          : T p d v cp cv gam e h s mw R    
 ig.N2O4    :  ig2  : Dinitrogen tetroxide         : T p d v cp cv gam e h s mw R    
 ig.N2O5    :  ig2  : Dinitrogen pentoxide         : T p d v cp cv gam e h s mw R    
 mp.N2      :  mp1  : Nitrogen                     : T p d v cp cv gam e h s mw R    
```

Using the Numpy library for compiled back-end numerical efficiency, `PYroMat` has been designed to efficiently handle large arrays of input data.  The iterative back-end operations on hundred-thousand-element arrays consistently execute in a few seconds on a single core of an old laptop.  

# Limitations and Future Plans

Developing a graphical web interface that exposes most (if not all) of `PYroMat`'s functionality is the next developmental objective for the project.  The goal is to provide a free interface that can be used without needing to know Python.  The majority of all corresponding users are either graduate students or thermodynamics instructors, so it seems this is the greatest barrier to expanding the project's user base.

The collection of multi-phase substances still lags significantly behind the likes of `coolprop` and `REFPROP`.  At present, it includes, methane, carbon dioxide, water/steam, nitrogen, oxygen, R-134a, and R-1234ze.  These models are time consuming to port to the package, so authors are choosing to strategically target substances of common industrial and scientific interest, like argon, hydrogen, ammonia, propane, R-22, R-12, and other refrigerants.

There is not currently support for transport physics (like viscosity or diffusivity), and users who want to calculate tertiary properties like Gibbs energy or Helmholtz energy are required to do so from more fudnamental properties.  These additions are planned for the distant future, but users have not yet asked for these features.

# Acknowledgements

PYroMat development was partially funded through an Affordable Course Transformation grant by Penn State Teaching and Learning with Technology.

# References
