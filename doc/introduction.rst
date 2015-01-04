**********************
Introduction on WOFOST
**********************

WOFOST (WOrld FOod STudies) is a simulation model for quantitative analysis of growth and production of annual field crops. It is a mechanistic model that explains crop growth and development based on the underlying processes, such as phenology, photosynthesis, respiration and how these processes are influenced by environmental conditions.

WOFOST calculates the attainable crop production, biomass, water use, etc. for a location given knowledge about soil type, crop type, weather data and crop management factors (e.g. sowing date). It has been used  by many researchers world wide and has been applied for various crops over a large range of climatic and management conditions. WOFOST is also implemented in the Crop Growth Monitoring System which is used operationally to monitor arable crops in Europe and to forecast crop yields for the current growing season.

WOFOST originated in the framework of interdisciplinary studies on world food security and on the potential world
food production by the Center for World Food Studies (CWFS) in cooperation with the Wageningen Agricultural
University, Department of Theoretical Production Ecology (WAU-TPE) and the DLO-Center for Agrobiological Research and
 Soil Fertility (AB-DLO), Wageningen, the Netherlands. After cessation of CWFS in 1988,
 the DLO Winand Staring Centre (SC-DLO) has continued the development co-operation with AB-DLO and WAU-TPE.

Currently, the WOFOST model is maintained and further developed by Alterra in co-operation with the Plant Production Systems Group of Wageningen University (http://www.pps.wur.nl/UK) and the Agri4Cast unit of the Joint Research Centre in Italy (http://mars.jrc.it/mars/About-us/AGRI4CAST).

About PyWOFOST and other WOFOST implementations
===============================================

Several WOFOST implementations exist. The original implementation dates back to the early 1990-ies and is written in FORTRAN77. This version has been well tested and used in many studies, but the implementation in FORTRAN77 limits the interfacing with other components such as databases and visualization applications. Next, WOFOST is implemented within the Crop Growth Monitoring System (CGMS) of the European Commission's Joint Research Centre (JRC). Although CGMS can be downloaded from the JRC webpage, the source code of the system is proprietary. Moreover, CGMS/WOFOST reads/writes to a fairly complex database and is implemented in the C++ language which is hard to use for agronomists and students who do not have an IT background. Finally, there is the implementation of WOFOST in the BioMa platform of JRC (http://bioma.jrc.ec.europa.eu/) but this version is not open source although the individual components are available in binary form.

PyWOFOST is a specific WOFOST implementation which is completely written in the python programming language. It is more flexible, easier to modify and also applies a modern language which allows for easy interfacing with databases, graphical interfaces, visualization tools and numerical/statistical packages.  PyWOFOST has several interesting features:

* Implementation in pure python with limited dependencies (only NumPy and SQLAlchemy).

* Modular design allowing you to add or change components relatively quickly

* Easily read/write from databases as well as files

* Built-in testing of program modules ensuring integrity of the system

Why python
==========
PyWOFOST was first and foremost developed from a scientific need, to be able to quickly adapt models and test ideas. In science, python is quickly becoming a tool for implementing algorithms, visualization and explorative analysis due to its clear syntax and ease of use. Many packages exist for numeric analysis (e.g. NumPy, SciPy), visualisation (e.g. MatPlotLib, Chaco), distributed computing (e.g. iPython, pyMPI) and interfacing with databases (e.g. SQLAlchemy). Moreover, for statistical analysis an interface with R-project can be established through Rpy or Rserve. Even large companies like Google and YouTube are heavy python users, MicroSoft has implemented python in its .NET framework (IronPython) and ESRI uses python as the main scripting language in its Geographical Information System (ArcGIS). Finally, python is an Open Source interpreted programming language that runs on almost any hardware and operating system.

Given the above considerations, it was quickly recognized that python was a good choice. Although, PyWOFOST was developed for scientific purposes, PyWOFOST can be used in a production environment and can for example easily replace the WOFOST implementation in CGMS.

History of PyWOFOST
===================

The version 5.0.0 of the PyWOFOST model is the first pure python implementation of PyWOFOST. Earlier implementations were mere wrappers around the FORTRAN77 processing core. These earlier versions worked well and have been used in a couple of studies. Nevertheless, they had some inherent limitations, particularly the biophysical core was much more difficult to modify. Also from a design point of view, the internal model structure had some problems due to the inherent difficulties in interfacing FORTRAN and python.

Requirements and dependencies
=============================

PyWOFOST was developed on Ubuntu Linux 10.04 using python 2.6.5, but is known to work with the 2.7 series (3.x series not tested yet). On Windows7, PyWOFOST is tested using version 7 of the Enthought Python Distribution 7.0 (EPD70), available from http://www.enthought.com/products/epd.php and free for academic use. The following screen dump shows the minimal requirements::

    Python 2.6.5 (r265:79063, Apr 16 2010, 13:57:41) 
    [GCC 4.4.3] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import numpy as np
    >>> np.__version__
    '1.3.0'
    >>> import sqlalchemy as sa
    >>> sa.__version__
    '0.5.8'
    >>> 

License
=======

The source code of PyWOFOST is made available under the European Union
Public License (EUPL), Version 1.1 or as soon they will be approved by the
European Commission - subsequent versions of the EUPL (the "Licence").
You may not use this work except in compliance with the Licence. You may obtain
a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl

The PyWOFOST package contains some modules that have been taken and/or modified
from other open source projects:

* the `pydispatch` module obtained from http://pydispatcher.sourceforge.net/
  which is distributed under a BSD style license.

* The `traitlets` and `import string` modules which were taken from the
  `ipython` project (http://ipython.org/) which are distributed under a
  BSD style license.

See the project pages of both projects for exact license terms.