Python Crop Simulation Environment - PCSE
=========================================

PCSE is a framework developed for implementing crop simulation models developed in
Wageningen. Many of the Wageningen crop simulation models were originally developed
in FORTRAN77 or using the FORTRAN Simulation Translator (`FST`_). Although this
approach has yielded high quality models with high numerical performance, the
inherent limitations of models written in FORTRAN is also becoming increasingly evident:

* The structure of the models is often rather monolithic and the different parts are
  very tightly coupled. Replacing parts of the model with another simulation approach
  is not easy.
* The models rely on file-based I/O which is difficult to change. For example,
  interfacing with databases is complicated in FORTRAN.
* In general, with low-level languages like FORTRAN, simple things already take many
  lines of code and mistakes are easily made, particularly by agronomists and crop
  scientist that have limited experience in developing or adapting software.

To overcome many of the limitations above, the Python Crop Simulation Environment
(PCSE) was developed which provides an environment for developing simulation models
as well as a number of implementations of crop simulation models. PCSE is written
in pure python code which makes it more flexible, easier to modify and extensible
allowing easy interfacing with databases, graphical user interfaces, visualization
tools and numerical/statistical packages. PCSE has several interesting features:

* Implementation in pure python and the core system has a small number of dependencies outside the python standard
  library. Most of these can be automatically installed from the Python Package Index (PyPI) (`SQLAlchemy`, `PyYAML`,
  `tabulate`, `xlwt`, `xlrd`) although some additional modules rely on `NumPy`.

* Modular design allowing you to add or change components relatively quickly with
  a simple but powerful approach to communicate variables between modules.

* Similar to `FST`_, it enforces good model design by explicitly separating parameters,
  rate variables and state variables. Moreover PCSE takes care of the module
  initialization, calculation of rates of changes, updating of state variables
  and actions needed to finalize the simulation.

* Input/Output is completely separated from the simulation model itself. Therefore
  PCSE models can easily read from and write to text files, databases and scientific
  formats such as HDF or NetCDF.

* Tools are available for reading parameter and weather files from existing models to
  have as much backward compatibility as possible.

* An `AgroManager` module which allows to define the agromanagement actions that
  happen on a farmers field. Such actions can be specified as events based on
  time or model state.

* Built-in testing of program modules ensuring integrity of the system.

To contribute to PCSE, you can fork your own copy at https://github.com/ajwdewit/pcse

Full documentation is available on http://pcse.readthedocs.io

.. _FST: http://models.pps.wur.nl/sites/models.pps.wur.nl/files/FST%203.pdf