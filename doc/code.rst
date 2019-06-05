.. include:: abbreviations.txt

How to read
===========

The API documentation provides a description of the interface and internals of 
all SimulationObjects, AncillaryObjects and utility routines available in the
PCSE source distribution. All SimulationObjects and AncillaryObjects are
described using the same structure:
    
    1. A short description of the object
    2. The positional parameters and keywords specified in the interface.
    3. A table specifying the simulation parameters needed for the simulation
    4. A table specifying the state variables of the SimulationObject
    5. A table specifying the rate variables of the SimulationObject
    6. Signals sent or received by the SimulationObject
    7. External dependencies on state/rate variables of other SimulationObjects.
    8. The exceptions that are raised under which conditions.
    
One or more of these sections may be excluded when they are not appropriate
for the SimulationObject that is described.

The table specifying the simulation parameters has the following columns:

    1. The name of the parameter.
    2. A description of the parameter.
    3. The type of the parameter. This is provided as a three-character code
       with the following interpretation. The first character indicates of the
       parameter is a scalar **(S)** or table **(T)** parameter. The second and
       third
       character indicate whether this parameter should be present in the
       timerdata '**Ti**', cropdata '**Cr**', soildata '**So**' or
       sitedata '**Si**' dictionary.
    4. The physical unit of the parameter.

The tables specifying state/rate variables have the following columns:
    1. The name of the variable.
    2. A description of the variable.
    3. Whether the variable is published in the kiosk or not: Y|N
    4. The physical unit of the variable.
    
Finally, all public methods of all objects are described as well.

Engine and models
=================

.. automodule:: pcse.engine
    :members:

.. automodule:: pcse.models
    :members:


.. _AgromanagementCode:

Agromanagement modules
======================

The routines below implement the agromanagement system in PCSE including crop calendars, rotations,
state and timed events. For reading agromanagement data from a file or a database structure see the sections
on the :ref:`reading file input <FileInput>` and the :ref:`database tools <DBtools>`.

.. autoclass:: pcse.agromanager.AgroManager
    :members:

.. autoclass:: pcse.agromanager.CropCalendar
    :members:

.. autoclass:: pcse.agromanager.TimedEventsDispatcher
    :members:

.. autoclass:: pcse.agromanager.StateEventsDispatcher
    :members:


The Timer
=========

.. autoclass:: pcse.timer.Timer
    :members:

The waterbalance
================

The PCSE distribution provides several waterbalance modules:
    1. WaterbalancePP which is used for simulation under non-water-limited
       production
    2. WaterbalanceFD which is used for simulation of water-limited production
       under conditions of freely draining soils
    3. The `SnowMAUS` for simulation the build-up and melting of the snow cover.
    4. A multi-layer waterbalance implementing simulations for potential
       conditions, water-limited free drainage conditions and
       water-limited groundwater conditions (in case of shallow ground
       water tables). This waterbalance is in a prototype stage and not yet
       usable, although the source code is available in PCSE.

.. autoclass:: pcse.soil.WaterbalancePP
    :members:

.. autoclass:: pcse.soil.WaterbalanceFD
    :members:

.. autoclass:: pcse.soil.SnowMAUS
    :members:

Crop simulation processes
=========================

Phenology
---------

.. autoclass:: pcse.crop.phenology.DVS_Phenology
    :members:

.. autoclass:: pcse.crop.phenology.Vernalisation
    :members:

Partitioning
------------

.. autoclass:: pcse.crop.partitioning.DVS_Partitioning
    :members:


|CO2| Assimilation
------------------

.. autoclass:: pcse.crop.assimilation.WOFOST_Assimilation
    :members:
    
Maintenance respiration
-----------------------
.. autoclass:: pcse.crop.respiration.WOFOST_Maintenance_Respiration
    :members:
    
Evapotranspiration
------------------
.. autoclass:: pcse.crop.evapotranspiration.Evapotranspiration
    :members:

.. autofunction:: pcse.crop.evapotranspiration.SWEAF

    
Leaf dynamics
-------------
.. autoclass:: pcse.crop.leaf_dynamics.WOFOST_Leaf_Dynamics
    :members:

Root dynamics
-------------
.. autoclass:: pcse.crop.root_dynamics.WOFOST_Root_Dynamics
    :members:

Stem dynamics
-------------
.. autoclass:: pcse.crop.stem_dynamics.WOFOST_Stem_Dynamics
    :members:

Storage organ dynamics
----------------------
.. autoclass:: pcse.crop.storage_organ_dynamics.WOFOST_Storage_Organ_Dynamics
    :members:

Abiotic damage
--------------
.. autoclass:: pcse.crop.abioticdamage.FROSTOL
    :members:

.. autoclass:: pcse.crop.abioticdamage.CrownTemperature
    :members:

.. .. autoclass:: pcse.crop.abioticdamage.CERES_WinterKill
..    :members:

.. _BaseClasses:

Base classes
============

The base classes define much of the functionality which is used "under the
hood" in PCSE. Except for the `VariableKiosk` and the `WeatherDataContainer`
all classes are not to be called directly but should be subclassed instead.


VariableKiosk
-------------
.. autoclass:: pcse.base.VariableKiosk
    :members:


Base classes for parameters, rates and states
---------------------------------------------

.. autoclass:: pcse.base.StatesTemplate
    :members:

.. autoclass:: pcse.base.RatesTemplate
    :members:

.. autoclass:: pcse.base.ParamTemplate
    :members:

Base and utility classes for weather data
-----------------------------------------

.. autoclass:: pcse.base.WeatherDataProvider
    :members:

.. autoclass:: pcse.base.WeatherDataContainer
    :members:

.. _Signals:

Signals defined
===============
.. automodule:: pcse.signals
    :members:
    

Utilities
=========

The utilities section deals with tools for reading weather data and parameter
values from files or databases.

.. _FileInput:

Tools for reading input files
-----------------------------

The file_input tools contain three classes: the `CABOFileReader` for reading in
parameter files in the CABO format,  the `CABOWeatherDataProvider` for reading
files from the CABO weather system and the `PCSEFileReader` for reading
files in PCSE format.

.. _CABOFileReader:
.. autoclass:: pcse.fileinput.CABOFileReader
    :members:

.. _CABOWeatherDataProvider:
.. autoclass:: pcse.fileinput.CABOWeatherDataProvider
    :members:

.. _PCSEFileReader:
.. autoclass:: pcse.fileinput.PCSEFileReader
    :members:

.. _ExcelWeatherDataProvider:
.. autoclass:: pcse.fileinput.ExcelWeatherDataProvider

.. _CSVWeatherDataProvider:
.. autoclass:: pcse.fileinput.CSVWeatherDataProvider

.. _YAMLAgroManagementReader:
.. autoclass:: pcse.fileinput.YAMLAgroManagementReader

.. _YAMLCropDataProvider:
.. autoclass:: pcse.fileinput.YAMLCropDataProvider


Simple or dummy data providers
------------------------------

This class of data providers can be used to provide parameter values in cases
where separate files or a database is not needed or not practical. An example
is the set of soil parameters for simulation of potential production conditions
where the value of the parameters does not matter but nevertheless some values
must be provided to the model.

.. _DummySoilDataProvider:
.. autoclass:: pcse.util.DummySoilDataProvider

.. _WOFOST71SiteDataProvider:
.. autoclass:: pcse.util.WOFOST71SiteDataProvider


.. _DBtools:

The database tools
------------------

The database tools contain functions and classes for retrieving agromanagement,
parameter values and weather variables from database structures implemented for
different versions of the European `Crop Growth Monitoring System <CGMS>`_.

Note that the data providers only provide functionality for *reading* data,
there are no tools here *writing* simulation results to a CGMS database. This was
done on purpose as writing data can be a complex matter and it is our
experience that this can be done more easily with dedicated database loader
tools such as `SQLLoader`_ for ORACLE or the ``load data infile`` syntax of MySQL

.. _SQLLoader: http://www.oracle.com/technetwork/database/enterprise-edition/sql-loader-overview-095816.html

.. _CGMS: https://ec.europa.eu/jrc/en/mars

.. _CGMS8tools:

The CGMS8 database
..................

The CGMS8 tools are for reading data from a database structure that is used
by CGMS executable version 9 and 10.

.. autoclass:: pcse.db.cgms8.GridWeatherDataProvider
    :members:

.. autoclass:: pcse.db.cgms8.SoilDataIterator
    :members:

.. autoclass:: pcse.db.cgms8.CropDataProvider
    :members:

.. autoclass:: pcse.db.cgms8.STU_Suitability
    :members:

.. autoclass:: pcse.db.cgms8.SiteDataProvider
    :members:

.. _CGMS12tools:

The CGMS12 database
...................

The CGMS12 tools are for reading data from a CGMS12 database structure that
is used by CGMS executable version 11 and BioMA 2014.


.. automodule:: pcse.db.cgms12
    :members:

.. autoclass:: WeatherObsGridDataProvider
    :members:

.. autoclass:: AgroManagementDataProvider
    :members:

.. autoclass:: SoilDataIterator
    :members:

.. autoclass:: CropDataProvider
    :members:

.. autoclass:: STU_Suitability
    :members:

.. autoclass:: SiteDataProvider
    :members:

.. _CGMS14tools:

The CGMS14 database
...................

The CGMS14 database is the database structure that is compatible with the 2015 BioMA implementation
of WOFOST. Note that the CGMS14 database structure is considerably different
from CGMS8 and CGMS12.

.. _CGMS14_data_providers:


The NASA POWER database
.......................

.. _NASAPowerWeatherDataProvider:

.. autoclass:: pcse.db.NASAPowerWeatherDataProvider
    :members:


Convenience routines
--------------------

These routines are there for conveniently starting a WOFOST simulation
for the demonstration and tutorials. They can serve as an example to
build your own script but have no further relevance.

.. autofunction:: pcse.start_wofost.start_wofost


Miscelaneous utilities
----------------------

Many miscelaneous for a variety of purposes such as the Arbitrary Function
Generator (AfGen) for linear interpolation and functions for calculating
Penman Penman/Monteith reference evapotranspiration,
the Angstrom equation and astronomical calculations such as day length.

.. autofunction:: pcse.util.reference_ET
.. autofunction:: pcse.util.penman_monteith
.. autofunction:: pcse.util.penman
.. autofunction:: pcse.util.check_angstromAB
.. autofunction:: pcse.util.wind10to2
.. autofunction:: pcse.util.angstrom
.. autofunction:: pcse.util.doy
.. autofunction:: pcse.util.limit
.. autofunction:: pcse.util.daylength
.. autofunction:: pcse.util.astro
.. autofunction:: pcse.util.merge_dict
.. autoclass:: pcse.util.Afgen
    :members:
.. autoclass:: pcse.util.ConfigurationLoader
    :members:
.. autofunction:: pcse.util.is_a_month
.. autofunction:: pcse.util.is_a_dekad
.. autofunction:: pcse.util.is_a_week
.. autofunction:: pcse.util.load_SQLite_dump_file
