*****************************
Getting started with PCSE
*****************************

This quickstart guide gives some examples to get you started with the Python
Crop Simulation Environment. All examples are currently focused on applying
the WOFOST crop simulation model, although other crop simulations may
become available within PCSE in the future.
All tests and examples were developed with the Enthought Python Distribution
(http://www.enthought.com/products/epd.php) under windows7 (EPD 32bit). This
product is now superseeded by Enthought Canopy.

Testing the PCSE package
========================
To guarantee its integrity, the PCSE package includes a number of self
tests that test individual components as well as the entire chain. These tests
verify that the output produced by the different components matches with the
expected outputs. Test data for the individual components can be found
in the `pcse.tests.test_data` package, while the test data for the entire chain
is stored in an SQLite database (pcse.db) that can be found under
`pcse.db.pcse`.

We assume here that PCSE is installed under 'D:\\USERDATA\\pylib\\' and
this location needs to be added to the search path of python::

    E:\temp>python
    Enthought Python Distribution -- www.enthought.com
    Version: 7.0-2 (32-bit)

    Python 2.7.1 |EPD 7.0-2 (32-bit)| (r271:86832, Dec  2 2010, 10:35:02) [MSC v.1500 32 bit (Intel)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import sys
    >>> import sys.path.append(r"D:\USERDATA\pylib\")

Next, PCSE can be imported and the tests can be executed by calling
the `test()` function at the top of the package::

    >>> import pcse
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

An interactive PCSE/WOFOST session
==================================
The easiest way to demonstrate PCSE is to import WOFOST from PCSE and run it from
an interactive Python session using the built-in demo database which contains
meteorologic data, soil data and crop data for a grid location in South-Spain.
We will use a

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
can be done with the `get_variable()` method on a PyWOFOST object.
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
model.

Finally, we can finish the crop season by simply specifying sufficient days
and store the results to a file 'myresults.csv'::

    >>> wofost_object.run(days=300)
    >>> wofost_object.store_to_file("myresults.txt")

Which should look like this :download:`myresults.txt`

Getting input data for PCSE/WOFOST
==================================

After running the examples you may be wondering where the data come
from that are used to run WOFOST. In fact, these data are retrieved from
an SQLite database 'pcse.db' that is included with the source distribution
and can be found in the pcse.db.pcse folder.

For setting up PCSE/WOFOST with your
own data sources you should understand that WOFOST uses 5 different types of
inputs: `cropdata`, `soildata`, `timerdata`, `sitedata` and driving variables
(e.g. weather data). The fact that these are called `\*data` is a bit of
misnomer as they contain a mixture of parameter values, boundary conditions
and events rather than data. Except for the driving variables which 
can be considered as (observed) data. This terminology was inherited from the 
previous WOFOST versions and it was kept because changing it would
cause more confusion.

All the input `\*data` must be provided as python dictionaries
storing key/value pairs and several tools are available in the PCSE
distribution to read these from a file or a database. Moreover,
there are several tools available for reading weather data.

For the second example we will run a simulation for sugar beet in
Wageningen (Netherlands) and we will read the input data step by step from
several different sources instead of using the pre-configured `start_wofost`
script.

Cropdata
--------

Cropdata consist of parameter names (dictionary keys) and the
corresponding parameter values that are needed to parameterize the
components of the crop simulation model. These are
crop-specific values regarding phenology, assimilation, respiration,
biomass partitioning, etc. The parameter file for sugar beet can be
downloaded here :download:`sug0601.cab` and is taken from the
crop files in the `WOFOST Control Centre`.

The crop parameter values for many models in
Wageningen are often provided in the CABO format that could be read
with the `TTUTIL <http://edepot.wur.nl/17847>`_ FORTRAN library. PCSE
tries to be backward compatible as much as possible and provides a
tool for reading parameter files in CABO format::

    >>> import pcse
    >>> from pcse.fileinput import CABOFileReader
    >>> cropdata = CABOFileReader('sug0601.crop')
    >>> print cropdata

printing the cropdata dictionary gives you an listing of the header and
all parameters and their values.

Soildata
--------

The soildata dictionary must provide the parameter name/value pairs related
to the soil type and soil physical properties. The number of parameters is
variable depending on the soil water balance type that is used for the
simulation. For this example, we will use the water balance for freely
draining soils and use the soil file for medium fine sand from the soil
files in the `WOFOST Control Centre` :download:`ec3.soil`::

    >>> soildata = CABOFileReader('ec3.soil')

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
    >>> timerdata = PCSEFileReader('sugarbeet_calendar.pcse')
    >>> print timerdata
    PCSE parameter file contents loaded from:
    /home/allard/Sources/python/pcse/doc/sugarbeet_calendar.pcse

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
the water balance such as the initial soil water in amount (WAV) and
the initial and maximum surface storage (SSI, SSMAX). For the moment, we will
define these parameters directly on the python commandline::

        >>> sitedata = {'SSMAX'  = 0.,
                        'IFUNRN' = 0,
                        'NOTINF' = 0,
                        'SSI'    = 0,
                        'WAV'    = 100,
                        'SMLIM'  = 0.03}

Driving variables (weather data)
--------------------------------

Daily weather variables are needed for running the simulation, see the section
on :ref:`DrivingVar` for reference. Currently, two options are available in
PyWOFOST for storing and retrieving weather data:

    1. The database structure as provided by the Crop Growth Monitoring
       System. Weather data will be read from the GRID_WEATHER table which
       is implemented using `db_util.GridWeatherDataProvider`.
    2. The file structure as defined by the CABO Weather System which is
       implemented using `cabo.CABOWeatherDataProvider`. For more details see
       :ref:`TheCABOtools`.

Additional WeatherDataProviders may be implemented as needed. Obvious
candidates are a WeatherDataProvider that can read the NASA POWER database for
agrometeorological modelling (http://power.larc.nasa.gov).

.. note::
    
    Recent CGMS versions (9 & 10) use the FAO Penman-Monteith approach to
    calculate reference evapotranspiration (`FAO56`_). However, the 
    `cabo.CABOWeatherDataProvider` still uses the classic modified Penman
    approach (`FAO24`_) which may result in slightly different reference
    evapotranspiration values. Future extensions to the 
    `cabo.CABOWeatherDataProvider` will most likely include the FAO56 method
    for estimating reference evapotranspiration.
    
.. _FAO56: http://www.fao.org/docrep/X0490E/x0490e00.htm#Contents
.. _FAO24: http://www4.fao.org/cgi-bin/faobib.exe?rec_id=131493&database=faobib&search_type=link&table=mona&back_path=/faobib/mona&lang=eng&format_name=EFMON


Building your own SimulationObjects
===================================

Often new modelling ideas need to be implemented and tested. Defining
new SimulationObjects is fairly simple and 
we will walk through an example SimulationObject (see :ref:`SimObjExample`)
line by line to understand the workings and the implementation details.

importing the base classes
--------------------------
First of all you need to make sure that the functionality needed is available
within your module. This is accomplished by importing the necessary
components from the pyWOFOST source distribution. This includes the template
classes for parameters, rate variables and state variables as well as the
base class `SimulationObject`. Finally, some decorators need to be imported as
well as the type definitions (in this case a Float) from the traitlets module.

The tutorial assumes that your new module is located inside the 'pywofost/'
folder otherwise you need to specify your imports like this::

    from pywofost.base_classes import ParamTemplate, StatesTemplate,
         RatesTemplate, SimulationObject
    from pywofost.decorators import prepare_rates, prepare_states
    from pywofost.traitlets import Float

Defining your SimulationObject
------------------------------
Defining a SimulationObject is demonstrated in line 6 and can be named any
valid python identifier with the base class `SimulationObject` between the
brackets. In this example, the SimulationObject is named `CropProcess`.
PyWOFOST does not enfore the use of a SimulationObject as base class
(you could to without it) but as much of the underlying logic is already
defined in your SimulationObject it would give you a hard time to do otherwise.

Next, the Parameters of the SimulationObject are defined, then the
StateVariables and the RateVariables in arbitrary order. Note that these
definitions are also class definitions, so we are defining here a class within
a class. This is nothing scary in python, it simply means that this definition
of `Parameters` is only known with the class `CropProcess`.

Defining parameters
-------------------
The Parameters section is where the simulation parameters are defined. The
parameter names are case sensitive and can be any valid python identifier,
except that it should not start with an underscore (_) or with the name 'trait'
as these will be ignored when assigning values. Various data types can
be used such as Float, Integer, Str or Instance (the
latter defines an instance of a certain class, for example a date).

A dummy value can be provided to indicate that the parameter value has not yet been
initialized (-99. in this case). The system enforces that all
parameters receive a value before the simulation can start. 

Defining rate and state variables
---------------------------------
The `StateVariables` section is where the variables are defined that describe
the state of the system at a certain time step. Similarly the `RateVariables`
section is
where the variables are defined that describe the rate of change from one time
step to the next. Regarding variable names, the same restrictions apply
as for the Parameters section.

Rate and state variables are *protected* in the sense that there value cannot
be changed outside the program sections where their values are updated. Thus,
state values cannot be changed during rate calculation and similarly, rates cannot be
changed during state updates.

.. note::

    Names of parameters, state variables and rate variables:
    
    - should *not* start with '_' or 'trait' 
      
    - are *case sensitive*, as are all variables in python.

The initialize() section
------------------------

In the `initialize` section preparations for the simulation take place.
The first compulsary positional argument is the
simulation day at which the object gets initialized. Secondly the VariableKiosk
must be specified. Next an arbitrary number of positional and keyword arguments
can be specified.

At the initialization step, first of all the kiosk is assigned to the
SimulationObject in order to have a reference to the kiosk available. Next,
the parameter values are initialized with the parvalues argument and assigned
to `self.params`. The parvalues argument is a dictionary specifying key/value
pairs where the keys are matched to parameter names. If a parameter cannot
be found in the dictionary, an error is raised. The dictionary can contain
more parameters than needed however this is ignored as these parameters may be
needed by other SimulationObjects. After initialization the `Parameters`
object is read-only and the parameter values cannot be changed anymore.

Next, the state variables must be assigned an initial value. This is done within
the `self.states = self.StateVariables(...)` section. The kiosk must be given
as first argument, furthermore an initial value must be provide for all state
variables as keywords. Neglecting this will results in
an errror. Finally, one can choose to publish one or more state variables,
which will make them available in the kiosk for other project sections.

Usually the last step is to initialize the rate variables, which is done with
`self.rates = self.RateVariables(kiosk)`. The initial value of a rate variable
does not need to be specified, it will be calculated at each time step.

The calc_rates() section
------------------------

In the `calc_rates()` section the rates of change are calculated for the
current `day` given current state variables and driving variables `drv` for the
current day. The definition `def calc_rates(day, drv)` **must** be
preceeded with the decorator `@prepare_rates` on the previous line. This
decorator prepares the rates object for updates, otherwise
you will not be able to assign new values to rate variables.

The variable `drv` is a container variable which contains the
different meteorological driving variables  which can
be accessed with a dotted notation. So `drv.TEMP` returns the daily mean
temperature (see :ref:`DrivingVar` for available driving variables).

In the rate variables section, the rates of change must be calculated
and assigned to the rate variables that are stored in `self.rates`. In
the example the rate variable RATE1 is calculated as function of PAR1 and
drv.IRRAD (radiation), while the rate variable RATE2 is calculated as a
function of PAR2 and the square of drv.TEMP.

Currently, the system does not check if all rate variables have been assigned
a value in the current time step. Although such checks could be implemented in
future versions.

The integrate() section
-----------------------

In the `integrate()` section the rates of change (that were calculated in the
prevous time step) are added to the state variables. Also here, the definition
`integrate(day)` **must** be preceeded with the decorator `@prepare_states`
in order to prepare `self.states` for states updates.

Since PyWOFOST works with daily time-steps and all rates of change are
calculated as <unit>/day, the rates of change can simply be added up to the
state variables without taking the size of the time step into account.
Similar to `calc_rates()`, there is no check that all state variables get
updated at each time step.


The finalize() section
----------------------
The `finalize()` section is optional. It can be defined when a
SimulationObject needs to 'wrap-up' calculations. Examples are the
water balance closure checks or the calculation of a harvest
index after the crop has finished. In the example, a ratio between the
state variables STATE1 and STATE2 is calculated which is assigned to the
state variable RATIO.

A definition of a finalize section must also be preceeded by the
`@prepare_states` at the preceeding line. This also implies that variables
that are calculated during the `finalize()` call *must* be defined in
StateVariables section.

When a `finalize()` section is defined the
`finalize` method on the superclass `SimulationObject` has to be called
because SimulationObjects
that are embedded in the current one need to be finalized as well. The
underlying logic is that a SimulationObject always receives a `finalize()`
method inherited from `SimulationObject`. The last line in the example
calls the `finalize()` of the base class to propagate the call to objects
lower in the hierarchy.

An alternative for finalize()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Defining a `finalize()` section can sometimes be tedious, particularly if there
is nothing to calculate and only some assignments need to be done. This occurs
when some information needs to be stored in `self.states` but you are not
inside the `integrate()` section so you are unable to assign to variables
in `self.states`.

For such simple cases a `finalize()` section is not needed.
Instead, each SimulationObject contains a dictionary variable
`self._for_finalize`. All key/value pairs that are added to this dictionary
will be assigned to the corresponding variable in `self.states` at the end of
the simulation.

For an example of the use of `self._for_finalize` see the `_on_CROP_FINISH`
method in cropsimulation.py

.. _SimObjExample:

An example SimulationObject
---------------------------

.. literalinclude:: example_simobj.py
   :linenos:
   :language: python

Embedding a new SimulationObject in PyWOFOST
--------------------------------------------


