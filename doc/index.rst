.. PCSE documentation master file, created by
   sphinx-quickstart on Sun Jul  1 23:03:43 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: abbreviations.txt

PCSE: The Python Crop Simulation Environment
============================================

PCSE (Python Crop Simulation Environment) is a Python package for building crop simulation models,
in particular the crop models developed in Wageningen (Netherlands). PCSE provides the
environment to implement crop simulation models, the tools for reading ancillary
data (weather, soil, agromanagement) and the components for simulating biophysical
processes such as phenology, respiration and evapotranspiration. PCSE also
includes implementations of the
`WOFOST <http://www.wageningenur.nl/wofost>`_ `LINGRA <https://edepot.wur.nl/336784>`_ and
`LINTUL3 <https://models.pps.wur.nl/system/files/LINTUL-N-Shibu-article_1.pdf>`_
crop and grassland simulation models
which have been widely used around the world. For example, WOFOST has been implemented in
the MARS crop yield forecasting system which is used operationally for crop monitoring and
yield prediction in Europe and beyond.

Originally, models developed in Wageningen were often written using FORTRAN or the
FORTRAN Simulation Environment (FSE). Both are very good tools, but they have become
somewhat outdated and are difficult to integrate with many of the great tools that are available
nowadays (XML, databases, web, etc).
Like so many other software packages, PCSE was developed to facilitate my own research work. I wanted something
that was more easy to work with, more interactive and more flexible while still implementing
the sound computational approach of FSE. For this reason PCSE was developed in Python,
which has become an important programming language for scientific purposes.

Traditionally, crop simulation models in Wageningen have been provided including the
full source code. PCSE is no exception and its source code is open and licensed under the
European Union Public License. PCSE runs on Python 2.7+ and 3.2+ and has a decent
test coverage of the implementation of the biophysical processes.

What's new
----------

.. toctree::
   :maxdepth: 2

   whatsnew.rst

Crop models Available in PCSE
-----------------------------

.. toctree::
   :maxdepth: 2

   available_models.rst

User guide
----------
.. toctree::
   :maxdepth: 2
   
   user_guide.rst

Reference guide
---------------
.. toctree::
   :maxdepth: 2

   reference_guide.rst

Code documentation
------------------

.. toctree::
   :maxdepth: 2

   code.rst

This document was generated on |date|/|time|.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

