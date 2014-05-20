*************************
Getting started with PCSE
*************************

This quickstart guide will help you to install PCSE and provides
some examples to get you started with modelling. All examples are currently focused on applying
the WOFOST crop simulation model, although other crop simulations may become available within
PCSE in the future.

Installing PCSE
===============

Requirements and dependencies
-----------------------------

PCSE is being developed on Ubuntu Linux 10.04 using python 2.7.6 and is known to work with the 3.x series (using the 2to3
tool). As python is a platform independent language, PCSE works equally well on Windows or Mac OSX.  The most
straightforward approach for installing python is through one of the prepackaged python distributions
such as `Enthought Canopy`_, `Anaconda`_ or `PythonXY`_.
The following screen dump shows the version of python, numpy and SQLAlchemy that were used to develop PCSE::

    Python 2.7.6 (default, Dec 16 2013, 12:39:22)
    [GCC 4.4.3] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import numpy as np
    >>> np.__version__
    '1.8.0'
    >>> import sqlalchemy as sa
    >>> sa.__version__
    '0.8.4'

.. _Enthought Canopy: https://store.continuum.io/cshop/anaconda/
.. _Anaconda: https://www.enthought.com/products/canopy/
.. _PythonXY: https://code.google.com/p/pythonxy/wiki/Welcome

All examples in this quickstart guide were developed under Windows 7 with the Enthought Python distribution
(version 7.0) which is by now superseeded by Enthought Canopy.

Downloading PCSE
----------------

The PCSE package can be downloaded as a zip file from GitHub.com using the link `here`_. Just unzip the package
at a suitable location. Note that the top directory in the zip file is `pcse-<branchname>`. The actual PCSE is package is
inside this folder and needs to be put on your file system.

.. _here: https://github.com/ajwdewit/pcse/archive/develop.zip


Testing the PCSE package
------------------------
To guarantee its integrity, the PCSE package includes a number of self
tests that test individual components as well as the entire simulation. These tests
verify that the output produced by the different components matches with the
expected outputs. Test data for the individual components can be found
in the `pcse.tests.test_data` package, while the test data for the entire chain
is stored in an SQLite database (pcse.db). This database can be found under
`.pcse` in your home folder and will be automatically generated when importing
PCSE for the first time. When you delete the database file manually it will be
regenerated..

We assume here that PCSE is installed under 'D:\\USERDATA\\pylib\\' and
this location needs to be added to the search path of python::

    C:\>python
    Enthought Python Distribution -- www.enthought.com
    Version: 7.0-2 (32-bit)

    Python 2.7.1 |EPD 7.0-2 (32-bit)| (r271:86832, Dec  2 2010, 10:35:02) [MSC v.1500 32 bit (Intel)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import sys
    >>> sys.path.append(r"D:\USERDATA\pylib\pcse")

Next, PCSE can be imported and the tests can be executed by calling
the `test()` function at the top of the package::

    >>> import pcse
    Building PCSE demo database at: C:\Users\wit015\.pcse\pcse.db
    >>> pcse.test()
    runTest (pcse.tests.test_abioticdamage.Test_FROSTOL) ... ok
    runTest (pcse.tests.test_assimilation.Test_WOFOST_Assimilation) ... ok
    runTest (pcse.tests.test_partitioning.Test_DVS_Partitioning) ... ok
    runTest (pcse.tests.test_evapotranspiration.Test_PotentialEvapotranspiration) ... ok
    runTest (pcse.tests.test_evapotranspiration.Test_WaterLimitedEvapotranspiration1) ... ok
    runTest (pcse.tests.test_evapotranspiration.Test_WaterLimitedEvapotranspiration2) ... ok
    runTest (pcse.tests.test_respiration.Test_WOFOSTMaintenanceRespiration) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedPotato) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialGrainMaize) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialSpringBarley) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedWinterWheat) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialSunflower) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedGrainMaize) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialWinterRapeseed) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedWinterRapeseed) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialPotato) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedSpringBarley) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedSunflower) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialWinterWheat) ... ok

    ----------------------------------------------------------------------
    Ran 19 tests in 29.748s

    OK
    >>>

If the model output matches the expected output the test will report 'OK',
otherwise an error will be produced with a detailed traceback on where the
problem occurred.

Part 1: An interactive PCSE/WOFOST session
==========================================

The easiest way to demonstrate PCSE is to import WOFOST from PCSE and run it from
an interactive Python session. We will be using the `start_wofost()` script that
connects to a the demo database which contains meteorologic data, soil data
and crop data for a grid location in South-Spain.

Initializing PCSE/WOFOST and advancing model state
--------------------------------------------------
Let's start a WOFOST object for modelling winter-wheat (crop=1) on a
location in South-Spain (grid 31031) for the year 2000 under water-limited
conditions for a freely draining soil (mode='wlp')::

    >>> wofost_object = pcse.start_wofost(grid=31031, crop=1, year=2000, mode='wlp')
    >>> type(wofost_object)
    <class 'pcse.models.Wofost71_WLP_FD'>

You have just successfully initialized a PCSE/WOFOST object in the python
interpreter, which is in its initial state and waiting to do some simulation. We
can now advance the model state for example with 1 day::

    >>> wofost_object.run()

Advancing the crop simulation with only 1 day, is often not so useful so the
number of days to simulate can be specified as well::

    >>> wofost_object.run(days=10)

Getting information about state and rate variables
--------------------------------------------------
Retrieving information about the calculated model states or rates 
can be done with the `get_variable()` method on a PCSE object.
For example, to retrieve the leaf area index value in the current
model state you can do::

    >>> wofost_object.get_variable('LAI')
    0.28708095263317146 
    >>> wofost_object.run(days=25)
    >>> wofost_object.get_variable('LAI')
    1.5281215808337203

Showing that after 11 days the LAI value is 0.287. When we increase time
with another 25 days, the LAI increases to 1.528. The `get_variable` method
can retrieve any state or rate variable that is defined somewhere in the
model. Finally, we can finish the crop season by simply specifying sufficient days
and store the results to a file 'myresults.csv'::

    >>> wofost_object.run(days=300)
    >>> wofost_object.store_to_file("myresults.txt")

Which should look like this :download:`myresults.txt`

Part 2: Running PCSE/WOFOST with custom input data
==================================================

For setting up PCSE/WOFOST with your
own data sources you should understand that WOFOST uses 5 different types of
inputs: `cropdata`, `soildata`, `timerdata`, `sitedata` and `driving variables`
(e.g. weather data). The fact that these names end with 'data' is a bit of
misnomer as they contain a mixture of parameter values, boundary conditions
and events rather than data, except for the driving variables which
can be considered as (observed) data. This terminology was inherited from the 
previous WOFOST versions and it was kept because changing it would
cause more confusion.

All the input `\*data` must be provided as python dictionaries
storing key/value pairs and several tools are available in the PCSE
distribution to read these from a file or a database. Moreover,
there are several tools available for reading weather data.

For the second example we will run a simulation for sugar beet in
Wageningen (Netherlands) and we will read the input data step by step from
several different sources instead of using the pre-configured `start_wofost()`
script. For the example we will assume that data files are in the directory
`D:\\userdata\\pcse_examples`. First we will import the necessary modules and
import set the data directory::

    >>> import os
    >>> import pcse
    >>> import matplotlib.pyplot as plt
    >>> data_dir = r'D:\userdata\pcse_examples'

Cropdata
--------

Cropdata consist of parameter names (dictionary keys) and the
corresponding parameter values that are needed to parameterize the
components of the crop simulation model. These are
crop-specific values regarding phenology, assimilation, respiration,
biomass partitioning, etc. The parameter file for sugar beet can be
downloaded here: :download:`sug0601.crop` and is taken from the
crop files in the `WOFOST Control Centre`_.

.. _WOFOST Control Centre: http://www.wageningenur.nl/wofost

The crop parameter values for many models in
Wageningen are often provided in the CABO format that could be read
with the `TTUTIL <http://edepot.wur.nl/17847>`_ FORTRAN library. PCSE
tries to be backward compatible as much as possible and provides a
tool for reading parameter files in CABO format::

    >>> from pcse.fileinput import CABOFileReader
    >>> cropfile = os.path.join(data_dir, 'sug0601.crop')
    >>> cropdata = CABOFileReader(cropfile)
    >>> print cropdata

printing the cropdata dictionary gives you an listing of the header and
all parameters and their values.

Soildata
--------

The soildata dictionary must provide the parameter name/value pairs related
to the soil type and soil physical properties. The number of parameters is
variable depending on the soil water balance type that is used for the
simulation. For this example, we will use the water balance for freely
draining soils and use the soil file for medium fine sand: :download:`ec3.soil`.
This file is also taken from the soil files in the `WOFOST Control Centre`_ ::

    >>> soilfile = os.path.join(data_dir, 'ec3.soil')
    >>> soildata = CABOFileReader(soilfile)

Timerdata
---------

The timerdata dictionary provides the start date of the water balance,
the start date and type of the crop simulation, the end date and type of the crop
simulation and the maximum duration of the crop simulation. The latter is
included to avoid unrealistically long simulations for example as a results of
a too high temperature sum requirement. These values are used by the AgroManagement
unit of PCSE. Currently, there is only an AgroManagement unit for single cropping
seasons but will change in the future allowing for crop rotations. Therefore,
the approach for providing AgroManagement data (timerdata) will change.

The following list gives an overview of the parameter names, values and types that
need to be specified in the `timerdata` dictionary::

        CAMPAIGNYEAR: year of the agricultural campaign (e.g. harvest year)
          START_DATE: date of the start of the simulation
            END_DATE: date last possible day of the simulation
     CROP_START_TYPE: 'emergence' or 'sowing'
     CROP_START_DATE: date of the start of the crop simulation
       CROP_END_TYPE: 'maturity' | 'harvest' |'earliest'
       CROP_END_DATE: date of the end of the crop simulation in case of CROP_END_TYPE == 'harvest' | 'earliest'
        MAX_DURATION: maximum number of days of the crop simulation

The CABO format has no support for dates, therefore the PCSE file format was
developed that does allow to use dates. The crop calendar file for sugar beet
in Wageningen can be downloaded here: :download:`sugarbeet_calendar.pcse`::

    >>> from pcse.fileinput import PCSEFileReader
    >>> crop_calendar_file = os.path.join(data_dir, 'sugarbeet_calendar.pcse')
    >>> timerdata = PCSEFileReader(crop_calendar_file)
    >>> print timerdata
    PCSE parameter file contents loaded from:
    D:\\userdata\\pcse_examples\\sugarbeet_calendar.pcse

    CAMPAIGNYEAR: 2000 (<type 'int'>)
    CROP_START_DATE: 2000-04-05 (<type 'datetime.date'>)
    END_DATE: 2000-12-31 (<type 'datetime.date'>)
    MAX_DURATION: 300 (<type 'int'>)
    CROP_END_DATE: 2000-10-20 (<type 'datetime.date'>)
    CROP_START_TYPE: emergence (<type 'str'>)
    CROP_END_TYPE: harvest (<type 'str'>)
    START_DATE: 2000-01-01 (<type 'datetime.date'>)

Sitedata
--------

The sitedata dictionary provides ancillary parameters that are not related to
the crop, the soil or the agromanagement. Examples are the initial conditions of
the water balance such as the initial soil moisture content (WAV) and
the initial and maximum surface storage (SSI, SSMAX). For the moment, we will
define these parameters directly on the python commandline::

    >>> sitedata = {'SSMAX'  : 0.,
                    'IFUNRN' : 0,
                    'NOTINF' : 0,
                    'SSI'    : 0,
                    'WAV'    : 100,
                    'SMLIM'  : 0.03}

Driving variables (weather data)
--------------------------------

Daily weather variables are needed for running the simulation. Currently, three
options are available in PCSE for retrieving weather data:

    1. The database structure as provided by the Crop Growth Monitoring
       System. Weather data will be read from the GRID_WEATHER table which
       is implemented using `pcse.db.pcse.GridWeatherDataProvider`.
    2. The file structure as defined by the `CABO Weather System`_ which is
       implemented using `pcse.fileinput.CABOWeatherDataProvider`.
    3. The global weather data provided by the agroclimatology from the
       `NASA Power database`_ at a resolution of 1x1 degree. PCSE
       provides the `pcse.db.NASAPowerWeatherDataProvider' which retrieves
       the NASA Power data from the internet for a given latitude and
       longitude.

.. _CABO Weather System: http://edepot.wur.nl/43010
.. _NASA Power database: http://power.larc.nasa.gov

For this example we will use the weather data from the NASA Power database
for the location of Wageningen. Note that it can take around 30 seconds
to retrieve the weather data from the NASA Power server the first time::

    >>> from pcse.db import NASAPowerWeatherDataProvider
    >>> wdp = NASAPowerWeatherDataProvider(latitude=52, longitude=5)
    >>> print wdp
    Weather data provided by: NASAPowerWeatherDataProvider
    --------Description---------
    NASA/POWER Agroclimatology Daily Averaged Data
    Dates (month/day/year): 01/01/1984 through 05/10/2014
    Location: Latitude 52   Longitude 5
    Location clarification: Integer values may indicate the lower left (south and west)
    corner of the one degree lat/lon region that includes the requested locations
    Elevation (meters): Average for one degree lat/lon region = 5
    Methodology Documentation:
    *Vegetation type: "Airport": flat rough grass
    ----Site characteristics----
    Elevation:    5.0
    Latitude:  52.000
    Longitude:  5.000
    Data available for 1997-01-01 - 2014-01-31
    Number of missing days: 47

Importing, initializing and running a PCSE model
------------------------------------------------

Internally, PCSE uses a simulation `engine` to run a crop simulation. This
engine takes a configuration file that specifies the components for the crop,
the soil and the agromanagement that need to be used for the simulation.
So any PCSE model can be started by importing the `engine` and initializing
it with a given configuration file and the corresponding sitedata, cropdata,
soildata, timerdata and weather data.

However, as many users of PCSE only need a particular configuration (for
example the WOFOST model for potential production), preconfigured Engines
are provided in `pcse.models`. For the sugarbeet example we will import
the WOFOST model for water-limited simulation under freely draining soils::

    >>> from pcse.models import Wofost71_WLP_FD
    >>> wofsim = Wofost71_WLP_FD(sitedata, timerdata, soildata, cropdata, wdp)

We can then run the simulation and show some final results such as the anthesis and
harvest dates (DOA, DOH), total biomass (TAGP) and maximum LAI (LAIMAX).
Next, we retrieve the time series of daily simulation output using the `get_output()`
method on the WOFOST object::

    >>> wofsim.run(days=400)
    >>> print wofsim.get_variable("DOA")
    2000-06-09
    >>> print wofsim.get_variable("DOH")
    2000-10-20
    >>> print wofsim.get_variable("TAGP")
    22783.5023325
    >>> print wofsim.get_variable("LAIMAX")
    5.11868342855
    >>> output = wofsim.get_output()
    >>> len(output)
    294

As the output is returned as a list of dictionaries, we need to unpack these variables
from the list of output::

    >>> varnames = ["day", "DVS", "TAGP", "LAI", "SM"]
    >>> tmp = {}
    >>> for var in varnames:
    >>>     tmp[var] = [t[var] for t in output]

Finally, we can generate some figures of WOFOST variables such as the
development (DVS), total biomass (TAGP), leaf area
index (LAI) and root-zone soil moisture (SM) using the `MatPlotLib`_ plotting package::

    >>> day = tmp.pop("day")
    >>> fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,8))
    >>> for var, ax in zip(["DVS", "TAGP", "LAI", "SM"], axes.flatten()):
    >>>     ax.plot_date(day, tmp[var], 'b-')
    >>>     ax.set_title(var)
    >>> fig.autofmt_xdate()
    >>> fig.savefig('sugarbeet.png')

.. _MatPlotLib: http://matplotlib.org/

This should provide generate a figure of the simulation results as shown below.


.. image:: sugarbeet.png

