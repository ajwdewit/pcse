Getting started
===============

This guide will help you install PCSE as well as provide
some examples to get you started with modelling. The examples are currently focused on applying
the WOFOST and LINTUL3 crop simulation models, although other crop simulation models may become available within
PCSE in the future.


An interactive PCSE/WOFOST session
==================================

The easiest way to demonstrate PCSE is to import WOFOST from PCSE and run it from
an interactive Python session. We will be using the `start_wofost()` script that
connects to a the demo database which contains meteorologic data, soil data,
crop data and management data for a grid location in South-Spain.

Initializing PCSE/WOFOST and advancing model state
..................................................

Let's start a WOFOST object for modelling winter-wheat (crop=1) on a
location in South-Spain (grid 31031) for the year 2000 under water-limited
conditions for a freely draining soil (mode='wlp')::

    >>> wofost_object = pcse.start_wofost(grid=31031, crop=1, year=2000, mode='wlp')
    >>> type(wofost_object)
    <class 'pcse.models.Wofost72_WLP_FD'>

You have just successfully initialized a PCSE/WOFOST object in the Python
interpreter, which is in its initial state and waiting to do some simulation. We
can now advance the model state for example with 1 day::

    >>> wofost_object.run()

Advancing the crop simulation with only 1 day, is often not so useful so the
number of days to simulate can be specified as well::

    >>> wofost_object.run(days=10)

Getting information about state and rate variables
..................................................
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
model. Finally, we can finish the crop season by letting it run until the
model terminates because the crop reaches maturity or the harvest date::

    >>> wofost_object.run_till_terminate()

Next we retrieve the simulation results at each time-step ('output') of the
simulation::

    >>> output = wofost_object.get_output()

We can now use the pandas package to turn the simulation output into a
dataframe which is much easier to handle and can be exported to different
file types. For example an Excel file which should look like this
:download:`downloads/wofost_results.xls`::

    >>> import pandas as pd
    >>> df = pd.DataFrame(output)
    >>> df.to_excel("wofost_results.xls")

Finally, we can retrieve the results at the end of the crop cycle (summary results)
and have a look at these as well::

    >>> summary_output = wofost_object.get_summary_output()
    >>> msg = "Reached maturity at {DOM} with total biomass {TAGP} kg/ha "\
    "and a yield of {TWSO} kg/ha."
    >>> print(msg.format(**summary_output[0]))
    Reached maturity at 2000-05-31 with total biomass 15261.7521735 kg/ha and a yield of 7179.80460783 kg/ha.

    >>> summary_output
    [{'CTRAT': 22.457536342947606,
      'DOA': datetime.date(2000, 3, 28),
      'DOE': datetime.date(2000, 1, 1),
      'DOH': None,
      'DOM': datetime.date(2000, 5, 31),
      'DOS': None,
      'DOV': None,
      'DVS': 2.01745939841335,
      'LAIMAX': 6.132711275237731,
      'RD': 60.0,
      'TAGP': 15261.752173534584,
      'TWLV': 3029.3693107257263,
      'TWRT': 1546.990661062695,
      'TWSO': 7179.8046078262705,
      'TWST': 5052.578254982587}]

Running PCSE/WOFOST with custom input data
------------------------------------------

For running PCSE/WOFOST (and PCSE models in general) with your own data sources you need three different types of
inputs:

1. Model parameters that parameterize the different model components. These parameters usually
   consist of a set of crop parameters (or multiple sets in case of crop rotations), a set of soil parameters
   and a set of site parameters. The latter provide ancillary parameters that are specific for a location.
2. Driving variables represented by weather data which can be derived from various sources.
3. Agromanagement actions which specify the farm activities that will take place on the field that is simulated
   by PCSE.

For the second example we will run a simulation for sugar beet in
Wageningen (Netherlands) and we will read the input data step by step from
several different sources instead of using the pre-configured `start_wofost()`
script. For the example we will assume that data files are in the directory
`D:\\userdata\\pcse_examples` and all the parameter files needed can be
found by unpacking this zip file :download:`downloads/quickstart_part2.zip`.

First we will import the necessary modules and define the data directory::

    >>> import os
    >>> import pcse
    >>> import matplotlib.pyplot as plt
    >>> data_dir = r'D:\userdata\pcse_examples'

Crop parameters
...............

The crop parameters consist of parameter names and the
corresponding parameter values that are needed to parameterize the
components of the crop simulation model. These are
crop-specific values regarding phenology, assimilation, respiration,
biomass partitioning, etc. The parameter file for sugar beet
is taken from the crop files in the `WOFOST Control Centre`_.

.. _WOFOST Control Centre: http://www.wageningenur.nl/wofost

The crop parameters for many models in
Wageningen are often provided in the CABO format that could be read
with the `TTUTIL <http://edepot.wur.nl/17847>`_ FORTRAN library. PCSE
tries to be backward compatible as much as possible and provides the
:ref:`CABOFileReader <CABOFileReader>` for reading parameter files in CABO format.
the CABOFileReader returns a dictionary with the parameter name/value pairs::

    >>> from pcse.fileinput import CABOFileReader
    >>> cropfile = os.path.join(data_dir, 'sug0601.crop')
    >>> cropdata = CABOFileReader(cropfile)
    >>> print(cropdata)

Printing the cropdata dictionary gives you a listing of the header and
all parameters and their values.

Soil parameters
...............

The soildata dictionary provides the parameter name/value pairs related
to the soil type and soil physical properties. The number of parameters is
variable depending on the soil water balance type that is used for the
simulation. For this example, we will use the water balance for freely
draining soils and use the soil file for medium fine sand: `ec3.soil`.
This file is also taken from the soil files in the `WOFOST Control Centre`_ ::

    >>> soilfile = os.path.join(data_dir, 'ec3.soil')
    >>> soildata = CABOFileReader(soilfile)

Site parameters
...............

The site parameters provide ancillary parameters that are not related to
the crop or the soil. Examples are the initial conditions of
the water balance such as the initial soil moisture content (WAV) and
the initial and maximum surface storage (SSI, SSMAX). Also the
atmospheric CO2 concentration is a typical site parameter.
For the moment, we can define these parameters directly on the Python commandline
as a simple python dictionary. However, it is more convenient to use the
:ref:`WOFOST71SiteDataProvider <WOFOST71SiteDataProvider>` that documents the
site parameters and provides sensible defaults::

    >>> from pcse.util import WOFOST71SiteDataProvider
    >>> sitedata = WOFOST71SiteDataProvider(WAV=100, CO2=360)
    >>> print(sitedata)
    {'SMLIM': 0.4, 'NOTINF': 0, 'CO2': 360.0, 'SSI': 0.0, 'SSMAX': 0.0, 'IFUNRN': 0, 'WAV': 100.0}

Finally, we need to pack the different sets of parameters into one variable
using the `ParameterProvider`. This is needed because PCSE expects one
variable that contains all parameter values. Using this approach has the
additional advantage that parameters value can be easily overridden in case
of running multiple simulations with slightly different parameter values::

     >>> from pcse.base import ParameterProvider
     >>> parameters = ParameterProvider(cropdata=cropdata, soildata=soildata, sitedata=sitedata)

AgroManagement
..............

The agromanagement inputs provide the start date of the agricultural campaign,
the start_date/start_type of the crop simulation, the end_date/end_type of the crop
simulation and the maximum duration of the crop simulation. The latter is
included to avoid unrealistically long simulations for example as a results of
a too high temperature sum requirement.

The agromanagement inputs are defined with a special syntax called `YAML`_ which allows
to easily create more complex structures which is needed for defining the agromanagement.
The agromanagement file for sugar beet in Wageningen `sugarbeet_calendar.agro` can be read with
the :ref:`YAMLAgroManagementReader <YAMLAgroManagementReader>`::

    >>> from pcse.fileinput import YAMLAgroManagementReader
    >>> agromanagement_file = os.path.join(data_dir, 'sugarbeet_calendar.agro')
    >>> agromanagement = YAMLAgroManagementReader(agromanagement_file)
    >>> print(agromanagement)
     !!python/object/new:pcse.fileinput.yaml_agro_loader.YAMLAgroManagementReader
     listitems:
     - 2000-01-01:
         CropCalendar:
           crop_name: sugarbeet
           variety_name: sugar_beet_601
           crop_start_date: 2000-04-05
           crop_start_type: emergence
           crop_end_date: 2000-10-20
           crop_end_type: harvest
           max_duration: 300
         StateEvents: null
         TimedEvents: null

Daily weather observations
..........................

Daily weather variables are needed for running the simulation. There are several
data providers in PCSE for reading weather data, see the section on
:ref:`weather data providers <Weather data providers>` to get an overview.

For this example we will use the weather data from the NASA Power database
which provides global weather data with a spatial resolution of 0.5 degree (~50 km).
We will retrieve the data from the Power database for the location of Wageningen.
Note that it can take around 30 seconds
to retrieve the weather data from the NASA Power server the first time::

    >>> from pcse.db import NASAPowerWeatherDataProvider
    >>> wdp = NASAPowerWeatherDataProvider(latitude=52, longitude=5)
    >>> print(wdp)
    Weather data provided by: NASAPowerWeatherDataProvider
    --------Description---------
    NASA/POWER SRB/FLASHFlux/MERRA2/GEOS 5.12.4 (FP-IT) 0.5 x 0.5 Degree Daily Averaged Data
    ----Site characteristics----
    Elevation:    4.7
    Latitude:  52.000
    Longitude:  5.000
    Data available for 1983-07-01 - 2018-09-16
    Number of missing days: 8

Importing, initializing and running a PCSE model
................................................

Internally, PCSE uses a simulation `engine` to run a crop simulation. This
engine takes a configuration file that specifies the components for the crop,
the soil and the agromanagement that need to be used for the simulation.
So any PCSE model can be started by importing the `engine` and initializing
it with a given configuration file and the corresponding parameters, weather
data and agromanagement.

However, as many users of PCSE only need a particular configuration (for
example the WOFOST model for potential production), preconfigured Engines
are provided in `pcse.models`. For the sugarbeet example we will import
the WOFOST model for water-limited simulation under freely draining soil
conditions::

    >>> from pcse.models import Wofost71_WLP_FD
    >>> wofsim = Wofost71_WLP_FD(parameters, wdp, agromanagement)

We can then run the simulation and show some final results such as the anthesis and
harvest dates (DOA, DOH), total biomass (TAGP) and maximum LAI (LAIMAX).
Next, we retrieve the time series of daily simulation output using the `get_output()`
method on the WOFOST object::

    >>> wofsim.run_till_terminate()
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

This should generate a figure of the simulation results as shown below. The complete Python
script for this examples can be downloaded here :download:`downloads/quickstart_demo2.py`

.. image:: figures/sugarbeet.png


.. _RunningLINTUL3:

Running a simulation with PCSE/LINTUL3
--------------------------------------

The LINTUL model (Light INTerception and UtiLisation) is a simple generic crop model, which simulates dry
matter production as the result of light interception and utilization with a constant light use efficiency.
In PCSE the LINTUL family of models has been implemented including the LINTUL3 model which is used for
simulation of crop production under water-limited and nitrogen-limited conditions.

For the third example, we will use LINTUL3 for simulating spring-wheat in the Netherlands under water-limited
and nitrogen-limited conditions. We will again assume that data files are in the directory
`D:\\userdata\\pcse_examples` and all the parameter files needed can be
found by unpacking this zip file :download:`downloads/quickstart_part3.zip`. Note that this guide is also available
as an IPython notebook: :download:`downloads/running_LINTUL3.ipynb`.

First we will import the necessary modules and define the data directory. We also assume that you have the
`matplotlib`_, `pandas`_ and `PyYAML`_ packages installed on your system.::

    >>> import os
    >>> import pcse
    >>> import matplotlib.pyplot as plt
    >>> import pandas as pd
    >>> import yaml
    >>> data_dir = r'D:\userdata\pcse_examples'

.. _pandas: http://pandas.pydata.org
.. _PyYAML: http://pyyaml.org/wiki/PyYAML

Similar to the previous example, for running the PCSE/LINTUL3 model we need to define the tree types of inputs
(parameters, weather data and agromanagement).

Reading model parameters
........................
Model parameters can be easily read from the input files using the `PCSEFileReader` as we have seen
in the previous example::

    >>> from pcse.fileinput import PCSEFileReader
    >>> crop = PCSEFileReader(os.path.join(data_dir, "lintul3_springwheat.crop"))
    >>> soil = PCSEFileReader(os.path.join(data_dir, "lintul3_springwheat.soil"))
    >>> site = PCSEFileReader(os.path.join(data_dir, "lintul3_springwheat.site"))

However, PCSE models expect a single set of parameters and therefore they need to be combined using the
`ParameterProvider`::

    >>> from pcse.base import ParameterProvider
    >>> parameterprovider = ParameterProvider(soildata=soil, cropdata=crop, sitedata=site)

Reading weather data
....................
For reading weather data we will use the ExcelWeatherDataProvider. This WeatherDataProvider uses nearly the same
file format as is used for the CABO weather files but stores its data in an MicroSoft Excel file which makes the
weather files easier to create and update::

    >>> from pcse.fileinput import ExcelWeatherDataProvider
    >>> weatherdataprovider = ExcelWeatherDataProvider(os.path.join(data_dir, "nl1.xlsx"))
    >>> print(weatherdataprovider)
    Weather data provided by: ExcelWeatherDataProvider
    --------Description---------
    Weather data for:
    Country: Netherlands
    Station: Wageningen, Location Haarweg
    Description: Observed data from Station Haarweg in Wageningen
    Source: Meteorology and Air Quality Group, Wageningen University
    Contact: Peter Uithol
    ----Site characteristics----
    Elevation:    7.0
    Latitude:  51.970
    Longitude:  5.670
    Data available for 2004-01-02 - 2008-12-31
    Number of missing days: 32

Defining agromanagement
.......................
Defining agromanagement needs a bit more explanation because agromanagement is a relatively
complex piece of PCSE. The agromanagement definition for PCSE is written in a format called `YAML`_ and
for the current example looks like this:

.. code:: yaml

    Version: 1.0.0
    AgroManagement:
    - 2006-01-01:
        CropCalendar:
            crop_name: wheat
            variety_name: spring-wheat
            crop_start_date: 2006-03-31
            crop_start_type: emergence
            crop_end_date: 2006-08-20
            crop_end_type: earliest
            max_duration: 300
        TimedEvents:
        -   event_signal: apply_n
            name:  Nitrogen application table
            comment: All nitrogen amounts in g N m-2
            events_table:
            - 2006-04-10: {amount: 10, recovery: 0.7}
            - 2006-05-05: {amount:  5, recovery: 0.7}
        StateEvents: null

.. _YAML: http://yaml.org/

The agromanagement definition starts with `Version:` indicating the version number of the agromanagement file
while the actual definition starts after the label `AgroManagement:`. Next a date must be provided which sets the
start date of the campaign (and the start date of the simulation). Each campaign is defined by zero or one
CropCalendars and zero or more TimedEvents and/or StateEvents. The CropCalendar defines the crop name,
variety_name, date of sowing, date of harvesting, etc. while the Timed/StateEvents define actions that are
either connected to a date or to a model state.

In the current example, the campaign starts on 2006-01-01, there is a crop calendar for spring-wheat starting on
2006-03-31 with a harvest date of 2006-08-20 or earlier if the crop reaches maturity before this date.
Next there are timed events defined for applying N fertilizer at 2006-04-10 and 2006-05-05. The current example
has no state events. For a thorough description of all possibilities see the section on AgroManagement in the
Reference Guide (Chapter 3).

Loading the agromanagement definition must by done with the YAMLAgroManagementReader::

    >>> from pcse.fileinput import YAMLAgroManagementReader
    >>> agromanagement = YAMLAgroManagementReader(os.path.join(data_dir, "lintul3_springwheat.amgt"))
    >>> print(agromanagement)
    !!python/object/new:pcse.fileinput.yaml_agro_loader.YAMLAgroManagementReader
    listitems:
    - 2006-01-01:
        CropCalendar:
          crop_end_date: 2006-10-20
          crop_end_type: earliest
          crop_name: wheat
          variety_name: spring-wheat
          crop_start_date: 2006-03-31
          crop_start_type: emergence
          max_duration: 300
        StateEvents: null
        TimedEvents:
        - comment: All nitrogen amounts in g N m-2
          event_signal: apply_n
          events_table:
          - 2006-04-10:
              amount: 10
              recovery: 0.7
          - 2006-05-05:
              amount: 5
              recovery: 0.7
          name: Nitrogen application table


Starting and running the LINTUL3 model
......................................
We have now all parameters, weather data and agromanagement information available to start the LINTUL3 model::

    >>> from pcse.models import LINTUL3
    >>> lintul3 = LINTUL3(parameterprovider, weatherdataprovider, agromanagement)
    >>> lintul3.run_till_terminate()

Next, we can easily get the output from the model using the get_output() method and turn it into a pandas DataFrame::

    >>> output = lintul3.get_output()
    >>> df = pd.DataFrame(output).set_index("day")
    >>> df.tail()
                     DVS       LAI     NUPTT       TAGBM     TGROWTH  TIRRIG  \
    day
    2006-07-28  1.931748  0.384372  4.705356  560.213626  626.053663       0
    2006-07-29  1.953592  0.368403  4.705356  560.213626  626.053663       0
    2006-07-30  1.974029  0.353715  4.705356  560.213626  626.053663       0
    2006-07-31  1.995291  0.339133  4.705356  560.213626  626.053663       0
    2006-08-01  2.014272  0.326169  4.705356  560.213626  626.053663       0

                   TNSOIL  TRAIN  TRAN  TRANRF  TRUNOF      TTRAN        WC  \
    day
    2006-07-28  11.794644  375.4     0       0       0  71.142104  0.198576
    2006-07-29  11.794644  376.3     0       0       0  71.142104  0.197346
    2006-07-30  11.794644  376.3     0       0       0  71.142104  0.196293
    2006-07-31  11.794644  381.6     0       0       0  71.142104  0.198484
    2006-08-01  11.794644  381.7     0       0       0  71.142104  0.197384

                     WLVD       WLVG        WRT         WSO         WST
    day
    2006-07-28  88.548865  17.687197  16.649830  184.991591  268.985974
    2006-07-29  89.284828  16.951234  16.150335  184.991591  268.985974
    2006-07-30  89.962276  16.273785  15.665825  184.991591  268.985974
    2006-07-31  90.635216  15.600845  15.195850  184.991591  268.985974
    2006-08-01  91.233828  15.002234  14.739974  184.991591  268.985974

Finally, we can visualize the results from the pandas DataFrame with a few commands if your
environment supports plotting::

    >>> fig, axes = plt.subplots(nrows=9, ncols=2, figsize=(16,40))
    >>> for key, axis in zip(df.columns, axes.flatten()):
    >>>     df[key].plot(ax=axis, title=key)
    >>> fig.autofmt_xdate()
    >>> fig.savefig(os.path.join(data_dir, "lintul3_springwheat.png"))

.. image:: downloads/lintul3_springwheat.png
