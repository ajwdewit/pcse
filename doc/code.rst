.. include:: abbreviations.txt

##################
Code documentation
##################

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
    
One or more of these sections may be excluded when they are not relevant
for the SimulationObject that is described.

The table specifying the simulation parameters has the following columns:

    1. The name of the parameter.
    2. A description of the parameter.
    3. The type of the parameter. This is provided as a three-character code
       with the following interpretation. The first character indicates of the
       parameter is a scalar **(S)** or table **(T)** parameter. The second and
       third
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

Soil process modules
====================

Water balance modules
---------------------

The PCSE distribution provides several waterbalance modules:
    1. WaterbalancePP which is used for simulation under non-water-limited
       production
    2. WaterbalanceFD which is used for simulation of water-limited production
       under conditions of freely draining soils
    3. The `SnowMAUS` for simulation the build-up and melting of the snow cover.
    4. A multi-layer waterbalance implementing simulations for potential
       conditions, water-limited free drainage conditions. Currently the model
       does not support the impact of shallow ground water tables but this will
       implemented in the future.

.. autoclass:: pcse.soil.WaterbalancePP

.. autoclass:: pcse.soil.WaterbalanceFD

.. autoclass:: pcse.soil.WaterBalanceLayered

.. autoclass:: pcse.soil.soil_profile.SoilProfile

.. autoclass:: pcse.soil.soil_profile.SoilLayer

.. autoclass:: pcse.soil.SnowMAUS

Nitrogen and Carbon modules
---------------------------

PCSE contains two modules for nitrogen and carbon in the soil:
    1. The simple N_Soil_Dynamics module which only simulates N availability as a pool of available N
       without any dynamic processes like leach, volatilization, etc.
    2. The SNOMIN module (Soil Nitrogen module for Mineral and Inorganic Nitrogen) which is a layered soil
       carbon/nitrogen balance that also requires the layered soil water balance. It includes the full
       N dynamics in the soil as well as the impact of organic matter and organic amendments (manure) on the
       availability of nitrogen in the soil.

.. autoclass:: pcse.soil.N_Soil_Dynamics

.. autoclass:: pcse.soil.SNOMIN


Crop simulation processes for WOFOST
====================================

Phenology
---------

.. autoclass:: pcse.crop.phenology.DVS_Phenology


.. autoclass:: pcse.crop.phenology.Vernalisation


Partitioning
------------

.. autoclass:: pcse.crop.partitioning.DVS_Partitioning


|CO2| Assimilation
------------------

.. autoclass:: pcse.crop.assimilation.WOFOST72_Assimilation

.. autoclass:: pcse.crop.assimilation.WOFOST73_Assimilation

.. autoclass:: pcse.crop.assimilation.WOFOST81_Assimilation


Maintenance respiration
-----------------------
.. autoclass:: pcse.crop.respiration.WOFOST_Maintenance_Respiration

Evapotranspiration
------------------
.. autoclass:: pcse.crop.evapotranspiration.Evapotranspiration

.. autoclass:: pcse.crop.evapotranspiration.EvapotranspirationCO2

.. autoclass:: pcse.crop.evapotranspiration.EvapotranspirationCO2Layered

.. autofunction:: pcse.crop.evapotranspiration.SWEAF

    
Leaf dynamics
-------------
.. autoclass:: pcse.crop.leaf_dynamics.WOFOST_Leaf_Dynamics

.. autoclass:: pcse.crop.leaf_dynamics.WOFOST_Leaf_Dynamics_N

.. autoclass:: pcse.crop.leaf_dynamics.CSDM_Leaf_Dynamics

Root dynamics
-------------
.. autoclass:: pcse.crop.root_dynamics.WOFOST_Root_Dynamics

Stem dynamics
-------------
.. autoclass:: pcse.crop.stem_dynamics.WOFOST_Stem_Dynamics

Storage organ dynamics
----------------------
.. autoclass:: pcse.crop.storage_organ_dynamics.WOFOST_Storage_Organ_Dynamics

Crop N dynamics
---------------

.. autoclass:: pcse.crop.n_dynamics.N_Crop_Dynamics
.. autoclass:: pcse.crop.nutrients.N_Demand_Uptake
.. autoclass:: pcse.crop.nutrients.N_Stress


Abiotic damage
--------------
.. autoclass:: pcse.crop.abioticdamage.FROSTOL

.. autoclass:: pcse.crop.abioticdamage.CrownTemperature


Crop simulation processes for LINGRA
====================================

.. automodule:: pcse.crop.lingra

Overall grassland model
-----------------------

.. autoclass:: pcse.crop.lingra.LINGRA

Source/Sink limited growth
--------------------------

.. autoclass:: pcse.crop.lingra.SourceLimitedGrowth

.. autoclass:: pcse.crop.lingra.SinkLimitedGrowth

Nitrogen dynamics
-----------------

.. autoclass:: pcse.crop.lingra_ndynamics.N_Demand_Uptake

.. autoclass:: pcse.crop.lingra_ndynamics.N_Stress

.. autoclass:: pcse.crop.lingra_ndynamics.N_Crop_Dynamics

Crop simulation processes for LINTUL
====================================

.. autoclass:: pcse.crop.lintul3.Lintul3
    :members:


.. Crop simulation processes for the ALCEPAS model
.. ===============================================


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

Configuration loading
---------------------
.. autoclass:: pcse.base.ConfigurationLoader
    :members:

.. _Signals:

Signals defined
===============
.. automodule:: pcse.signals
    :members:
    

Ancillary code
==============

The ancillary code section deals with tools for reading weather data and parameter
values from files or databases.


.. _Input:
Data providers
--------------

The module `pcse.input` contains all classes for reading weather files,
parameter files and agromanagement files.

.. _NASAPowerWeatherDataProvider:
.. autoclass:: pcse.input.NASAPowerWeatherDataProvider

.. _OpenMeteoWeatherDataProvider:
.. autoclass:: pcse.input.OpenMeteoWeatherDataProvider

.. _CABOWeatherDataProvider:
.. autoclass:: pcse.input.CABOWeatherDataProvider

.. _ExcelWeatherDataProvider:
.. autoclass:: pcse.input.ExcelWeatherDataProvider

.. _CSVWeatherDataProvider:
.. autoclass:: pcse.input.CSVWeatherDataProvider

.. _CABOFileReader:
.. autoclass:: pcse.input.CABOFileReader

.. _PCSEFileReader:
.. autoclass:: pcse.input.PCSEFileReader

.. _YAMLAgroManagementReader:
.. autoclass:: pcse.input.YAMLAgroManagementReader

.. _YAMLCropDataProvider:
.. autoclass:: pcse.input.YAMLCropDataProvider

.. _WOFOST72SiteDataProvider:
.. autoclass:: pcse.input.WOFOST72SiteDataProvider

.. _WOFOST73SiteDataProvider:
.. autoclass:: pcse.input.WOFOST73SiteDataProvider

.. _WOFOST81SiteDataProvider_classic:
.. autoclass:: pcse.input.WOFOST81SiteDataProvider_Classic

.. _WOFOST81SiteDataProvider_SNOMIN:
.. autoclass:: pcse.input.WOFOST81SiteDataProvider_SNOMIN


Simple or dummy data providers
------------------------------

This class of data providers can be used to provide parameter values in cases
where separate files or a database is not needed or not practical. An example
is the set of soil parameters for simulation of potential production conditions
where the value of the parameters does not matter but nevertheless some values
must be provided to the model.

.. _DummySoilDataProvider:
.. autoclass:: pcse.util.DummySoilDataProvider


.. _DBtools:

The database tools
------------------

.. note::
    The dataproviders for CGMS database were removed from PCSE starting with version 6.0.10 because they were forcing
    a dependency on PCSE (SQLAlchemy) which was generating problems with other packages. Moreover SQLAlchemy is not
    required for running PCSE and the DB tools were not broadly used anyway.

The database tools contain functions and classes for retrieving agromanagement,
parameter values and weather variables from database structures implemented for
different versions of the European `Crop Growth Monitoring System <CGMS>`_.

Note that the data providers only provide functionality for *reading* data,
there are no tools here *writing* simulation results to a CGMS database. This was
done on purpose as writing data can be a complex matter and it is our
experience that this can be done more easily with dedicated database loader
tools such as `SQLLoader`_ for ORACLE or the ``load data infile`` syntax of MySQL.


.. _SQLLoader: http://www.oracle.com/technetwork/database/enterprise-edition/sql-loader-overview-095816.html

.. _CGMS: https://ec.europa.eu/jrc/en/mars

.. _CGMS8tools:



Convenience routines
--------------------

These routines are there for conveniently starting a WOFOST simulation
for the demonstration and tutorials. They can serve as an example to
build your own script but have no further relevance.

.. autofunction:: pcse.start_wofost.start_wofost


Miscelaneous utilities
----------------------

Many miscelaneous function for a variety of purposes such as the Arbitrary Function
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
.. autofunction:: pcse.util.is_a_month
.. autofunction:: pcse.util.is_a_dekad
.. autofunction:: pcse.util.is_a_week
.. autofunction:: pcse.util.load_SQLite_dump_file
