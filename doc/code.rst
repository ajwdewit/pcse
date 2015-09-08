.. include:: abbreviations.txt

******************************
Model code documentation (API)
******************************

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
    
Finally, all public methods of all object are described as well.

Engine and models
=================

.. automodule:: pcse.engine
    :members:

.. automodule:: pcse.models
    :members:


.. _Agromanagement:

Agromanagement
==============

.. note::
   Currently two modules are available, the new `AgroManager` and the old `AgroManagementSingleCrop`.
   The old `AgroManagementSingleCrop` is still here for backward compatibility but will be removed
   in future versions of PCSE.

.. automodule:: pcse.agromanager
    :members:

.. automodule:: pcse.agromanagement
    :members:

The Timer
=========

.. autoclass:: pcse.timer.Timer
    :members:

The waterbalance
================

The PCSE distribution provides several waterbalance modules:
    1. WaterbalancePP which is used for simulation under non-waterlimited
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

Base classes
============

The base classes define much of the functionality which is used "under the
hood" in PCSE. Except for the `VariableKiosk` and the `WeatherDataContainer`
all classes are not to be called directly but should be subclassed instead.


VariableKiosk
-------------
.. autoclass:: pcse.base_classes.VariableKiosk
    :members:


Base classes for parameters, rates and states
---------------------------------------------

.. autoclass:: pcse.base_classes.StatesTemplate
    :members:

.. autoclass:: pcse.base_classes.RatesTemplate
    :members:

.. autoclass:: pcse.base_classes.ParamTemplate
    :members:

Base and utility classes for weather data
-----------------------------------------

.. autoclass:: pcse.base_classes.WeatherDataProvider
    :members:

.. autoclass:: pcse.base_classes.WeatherDataContainer
    :members:

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

The database tools
------------------

The database tools contain functions and classes for retrieving parameters
and weather database from several database structures. The database
structure is mostly derived from the database used for the Crop Growth
Monitoring System (CGMS_).

.. _CGMS: http://mars.jrc.ec.europa.eu/mars/About-us/AGRI4CAST/Models-Software-Tools/Crop-Growth-Monitoring-System-CGMS

The PCSE database
.................

The PCSE database structure is very similar to a CGMS9 structure but has some
modifications for dealing with dates in the CROP_CALENDAR table and uses
different table names and structure for model output.

.. autofunction:: pcse.db.pcse.fetch_cropdata
.. autofunction:: pcse.db.pcse.fetch_soildata
.. autofunction:: pcse.db.pcse.fetch_timerdata
.. autofunction:: pcse.db.pcse.fetch_sitedata

.. _GridWeatherDataProvider:

.. autoclass:: pcse.db.pcse.GridWeatherDataProvider
    :members:

The CGMS9 database
..................

The CGMS9 tools are for reading data from a database structure that is used
by CGMS executable version 9 and 10.

.. autoclass:: pcse.db.cgms9.GridWeatherDataProvider
    :members:

The CGMS11 database
...................

Tools for reading a CGMS11 database are under construction.

The NASA POWER database
.......................

.. _NASAPowerWeatherDataProvider:

.. autoclass:: pcse.db.NASAPowerWeatherDataProvider
    :members:


Convenience routines
--------------------

These routines are there for conveniently starting a WOFOST simulation. They
are mainly used in the tutorial and examples but can be used to further
elaborating on when writing your own scripts.

.. autofunction:: pcse.start_wofost.start_wofost

.. autofunction:: pcse.run_wofost.run_wofost


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
