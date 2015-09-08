******************
Background of PCSE
******************

The `Python Crop Simulation Environment` was developed because of a need to re-implement crop simulation
models that were developed in Wageningen. Many of the Wageningen crop simulation models were originally developed in
FORTRAN77 or using the `FORTRAN Simulation Translator (FST)`. Although this approach has yielded high quality models
with high numerical performance, the inherent limitations of models written in FORTRAN is also becoming increasingly
evident:

* The structure of the models is often rather monolithic and the different parts are very tightly coupled.
  Replacing parts of the model with another simulation approach is not easy.

* The models rely on file-based I/O which is difficult to change. For example, interfacing with databases
  is complicated in FORTRAN.

* In general, with low-level languages like FORTRAN, simple things already take many lines of code and mistakes
  are easily made, particularly by agronomists and crop scientist that have limited experience in developing or
  adapting software.

To overcome many of the limitations above, the Python Crop Simulation Environment (PCSE) was developed which provides
an environment for developing simulation models as well as a number of implementations of crop simulation models.
PCSE is written in pure python code which makes it more flexible, easier to modify and extensible allowing easy
interfacing with databases, graphical user interfaces, visualization tools and numerical/statistical packages. PCSE has
several interesting features:

* Implementation in pure python and the core system has a small number of dependencies outside the python standard
  library. Most of these can be automatically installed from the Python Package Index (PyPI) (`SQLAlchemy`, `PyYAML`,
  `tabulate`, `xlwt`, `xlrd`) although some additional modules rely on `NumPy`.

* Modular design allowing you to add or change components relatively quickly with a simple but powerful approach
  to communicate variables between modules.

* Similar to FST, it enforces good model design by explicitly separating parameters, rate variables and state
  variables. Moreover PCSE takes care of the module initialization, calculation of rates of changes, updating
  of state variables and actions needed to finalize the simulation.

* Input/Output is completely separated from the simulation model itself. Therefore PCSE models can easily
  read from and write to text files, databases and scientific formats such as HDF or NetCDF.

* Built-in testing of program modules ensuring integrity of the system

Why python
==========
PCSE was first and foremost developed from a scientific need, to be able to quickly adapt models and test ideas.
In science, python is quickly becoming a tool for implementing algorithms, visualization and explorative analysis
due to its clear syntax and ease of use. An additional advantage is that the C implementation of python
can be easily interfaces with routines written in FORTRAN and therefore many FORTRAN routines can be reused by
simulation models written with PCSE.

Many packages exist for numeric analysis (e.g. NumPy, SciPy),
visualisation (e.g. MatPlotLib, Chaco), distributed computing (e.g. iPython, pyMPI) and interfacing with databases
(e.g. SQLAlchemy). Moreover, for statistical analyses an interface with R-project can be established through
Rpy or Rserve. Even large companies like Google and YouTube are heavy python users, MicroSoft has implemented python
in its .NET framework (IronPython) and ESRI uses python as the main scripting language in its
Geographical Information System (ArcGIS). Finally, python is an Open Source interpreted programming language that
runs on almost any hardware and operating system.

Given the above considerations, it was quickly recognized that python was a good choice. Although, PCSE was
developed for scientific purposes, it has already been implemented for some tasks in a production environment.

History of PCSE
===============

Up till version 4.1, PCSE was called "PyWOFOST" as its primary goal was to provide an implementation of the
WOFOST crop simulation model in python.
However, as the system has grown it has become evident that the system can be used to implement, extend or
hybridize (crop) simulation models. Therefore, the name "PyWOFOST" became too narrow and the name Python Crop
Simulation Environment was selected in analog with the FORTRAN Simulation Environment (FSE).


Limitations of PCSE
===================

PCSE also has its limitations, in fact there are several:

* Speed: flexibility comes a at a price; PCSE is considerably slower than equivalent models written in FORTRAN or
  another compiled language.

* The simulation approach in PCSE is currently limited to rectangular (Euler) integration with a fixed daily
  time-step. Although the internal time-step of modules can be made more fine-grained if needed.

* No graphical user interface. However the lack of a user interface is partly compensated by using PCSE with the
  `pandas <http://pandas.pydata.org/>`_ package and the `IPython notebook <http://ipython.org/notebook.html>`_.
  PCSE output can be easily converted to a pandas `DataFrame` which can be used to display charts in an IPython
  notebook.

License
=======

The source code of PCSE is made available under the European Union
Public License (EUPL), Version 1.1 or as soon they will be approved by the
European Commission - subsequent versions of the EUPL (the "Licence").
You may not use this work except in compliance with the Licence. You may obtain
a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl

The PCSE package contains some modules that have been taken and/or modified
from other open source projects:

* the `pydispatch` module obtained from http://pydispatcher.sourceforge.net/
  which is distributed under a BSD style license.

* The `traitlets` and `import string` modules which were taken and adapted from the
  `ipython` project (http://ipython.org/) which are distributed under a
  BSD style license.

See the project pages of both projects for exact license terms.