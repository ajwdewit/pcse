*****************************
General structure of PyWOFOST
*****************************

PyWOFOST has a rather different design philosophy compared to the heritage
FORTRAN WOFOST model. Instead of using relatively large blocks of code
where many calculations are carried out, PyWOFOST uses SimulationObjects. 
Parts of the crop simulation model that form an
entity are grouped into separate SimulationObjects. In this way the crop
simulation model is grouped into logical sections that implement certain
biophysical processes such as phenology, assimilation, respiration, etc.

This approach has several advantages:

* Model code with a certain purpose is grouped together, making it easier
  to read, understand and maintain.
* A SimulationObject contains only parameters, rate and state variables
  that are needed. This is differs from FORTRAN WOFOST where a
  large number of variables are defined. It is often unclear (at
  first glance at least) what biophysical process they belong to.
* Implementation of isolated processes causing
  less dependencies which makes it easier to modify individual
  SimulationObjects.
* SimulationObjects can be tested individually by comparing output vs the
  expected output (e.g. unit testing).
* SimulationObjects can be exchanged for other objects with the same purpose
  but a different biophysical approach. For example, the WOFOST assimilation
  approach could be easily replaced by a more simple Light Use Efficiency or
  Water Use Efficiency approach, only by replacing the SimulationObject that
  handles the |CO2| assimilation.

.. |CO2| replace:: CO\ :sub:`2`\


Characteristics of SimulationObjects
====================================

The actual calculations are all organized in SimulationObjects.
Each SimulationObject is defined in the same way and has a couple of standard
sections and methods which facilitates understanding and readability (in some
cases a slightly modified approach is used, which will be discussed later).
Each SimulationObject has parameters to define the mathematical relationships,
it has state variables to define the state of the system and it has rate
variables that describe the rate of change from one time step to the next.
Moreover, a SimulationObject often contains other SimulationObjects that
together form a logical structure. Finally, driving variables (meteorology)
must be retrieved and managed to fit into the program structure. 

How a SimulationObject simulates
--------------------------------

In a SimulationObject, continuous simulation is implemented in separate code
sections that implement initialization, rate calculation and integration of the
rates of change. Retrieval and calculation of driving variables is
implemented in different program sections. Moreover, time is managed by a
specialized timer module.

The strict separation of program logic was copied by the Fortran Simulation
Environment (FSE, Rappoldt and Van Kraalingen 1996. http://edepot.wur.nl/4411
and Van Kraalingen 1995, http://edepot.wur.nl/35555) and 
is needed to ensure that the simulation results are correct.
The different calculations types (integration, driving variables and
rate calculations) should be strictly separated. In other words, first all
states should be updated, subsequently all driving variables should be calculated,
after which all rates of change should be calculated. If this rule is not
applied rigorously, some rates may pertain to states at
the current time whereas others will pertain to states from the previous time
step (Van Kraalingen, 1995). 

To implement the continuous simulation, PyWOFOST uses the same approach as
FSE: Euler integration with a
fixed time step of one day. Previous versions of WOFOST could be applied on
10-daily or even monthly time steps but this option has been discarded in
PyWOFOST. The following figures shows the principle of continuous simulation
and the execution order of various steps.

.. figure:: continuous_simulation.png
    :align: center
    :scale: 50%
    
    Order of calculations for continuous simulation using Euler integration
    (after Van Kraalingen, 1995).

From the previous figure it is clear that before the simulation can start an
initialization has to be carried out which comprises of:

1. The timer module must be initialized on the starting day;
2. All parameters of all SimulationObjects must be assigned a value;
3. All state variables of all SimulationObjects must be assigned an initial
   value;
4. The driving variables must be retrieved for the starting day;
5. The initial rates of change based on the initial states and driving
   variables must be calculated;
6. Finally, output can be collected to save the initial states and rates of
   the simulation. 

Step 1 is a system-wide initialization as there is only a single timer module.
The steps 2 and 3 are carried out by the the `initialize()` section of 
each SimulationObject. Next, the driving variables are retrieved (step 4) and
used as input into the `calc_rates()` section of each SimulationObject 
to calculate the initial rates of change (step 5). Finally, model state or
rate variables may be saved for output (step 6).

A point of attention is that the initial rates are not calculated
within the `initialize()` section of each SimulationObject.
Instead the calculation is delayed after all SimulationObject have been
initialized, because the initial rates of one object can depend on the
initial states of other objects. The execution of the initial rates
calculation is therefore handled by the SimulationObject in the top of the
hierarchy.

The next cycle in the simulation will now start with an update of the timer to
the next time step (e.g. day) and the integration of the rates of change of
the previous day onto the state variables of each SimulationObject. The latter
is carried out in the `integrate()` section of each SimulationObject. Next,
the driving variables for the current day are retrieved. 
Calculation of the rates of change based on the new driving variables and
updated model states and so forth.

The simulation loop will terminate when some finish condition has been reached.
The condition can be set by various model components. For example,
the final day of the simulation period can be reached, the `Phenology` module
can signal that the crop has reached maturity or the `Leaf_Dynamics` module
finds that the crop has died (e.g. no more living leaves).

When the simulation loop terminates, the `finalize()` section of each
SimulationObject will be called in order to allow calculations that have to be
deferred to the end of the simulation, such as a harvest index or water balance
totals to check the closure of the water balance.

Compared to the FSE system, the `initialize()`, `calc_rates()`, `integrate()`
and `finalize()` sections match with the *ITASK* numbers 1, 2, 3, 4. 

Communication between SimulationObjects
---------------------------------------

A complicating factor that arises when using modular code is how to arrange
the communication between SimulationObjects. For example, the `evapotranspiration`
SimulationObject will need information about the leaf area index from the
`leaf_dynamics` SimulationObject to calculate the crop transpiration
values.

In FORTRAN WOFOST, this was
handled by putting the variables that needed to be communicated into the 
subroutines calls. This leads to large unwieldy
routine calls with a large risk of mapping variables onto the wrong variable
names. Moreover, the parameters in the call are often a mixture of model
parameters, driving variables, rate/state variables and variables that are
purely for program logic.  See for example:

.. code-block:: fortran

      SUBROUTINE CROPSI
     &         (ITASK, IDAY  , DELT , TIME , IDEM, DOANTH, IDHALT,
     &         TERMNL, ISTATE, IWB  , IOX  ,
     &         LAT   , AVRAD , TMIN , TMAX , E0  , ES0, ET0,
     &         CRFILE, IUPL  , IUOUT, IULOG,
     &         SM    , SM0   , SMFCF, SMW  , CRAIRC,
     &         EVWMX , EVSMX , TRA  , FR   , RRI   , IAIRDU,
     &         RDI   , RDMCR)


In PyWOFOST the communication between
SimulationObjects is taken care of by the so-called `VariableKiosk`. The
metaphore kiosk is used because the SimulationObjects publish
their rate and/or state variables (or a subset)into the kiosk, other
SimulationObjects can subsequently request the variable value from the kiosk
without any knowledge about the SimulationObject that published it.
Note that the VariableKiosk is shared by all SimulationObjects.


The VariableKiosk
-----------------

The VariableKiosk is an integral part of PyWOFOST and it has several
important functions. First of all it takes 
care of communicating variables between SimulationObjects (as already
mentioned in the previous paragraph). When a variable has been published,
PyWOFOST updates the state and rate variables into the `VariableKiosk` at each
simulation cycle. Under the hood, this functionality is provided by the
`traitlets` module from which all SimulationObjects are derived.

To avoid that state or rate variables of previous time steps keep
lagging in the VariableKiosk, the kiosk contents kiosk are flushed during each
time step. After all rates have been calculated, the values of all state
variables are flushed from the VariableKiosk. Similarly, after the update of the state 
variables, the rate variables are flushed from the kiosk.

Second, to enforce that variable names are unique across the entire model
the VariableKiosk registers all variables across all SimulationObjects.
Variable names that are defined in SimulationObjects are automatically
registered in the VariableKiosk and the name is checked for uniqueness.
Moreover, it is tracked which SimulationObject registers which variable and
only the SimulationObject that registers and publishes a variable can change
its value. All other objects can retrieve that value, but not change it.
Uniqueness of variable names is important for retrieving their value; if two
SimulationObject would define the same variable than Python would only find the
first one. It is not even guaranteed that his would be the same
variable between sessions.

Finally, the VariableKiosk serves as a unique identifier for a PyWofost instance.
Since the VariableKiosk is shared across all SimulationObjects within the
hierarchy, it can be used as a unique identifier of a PyWOFOST model instance
by retrieving its identifier with python's `id()` function. This aspect is 
important when running ensembles of PyWOFOST models in combination with
sending and receiving signals (see :ref:`EventsAndSignals`).

Retrieving variables
--------------------
By calling the `get_variable(<varname>)` method of a PyWOFOST instance,
state or rate variables are retrieved. The
`get_variable` method first searches for `<varname>` within its own definitions
of state and rate variables and returns the value if it is found. If not,
it searches for other embedded SimulationObjects and calls their
`get_variable(<varname>)` methods. This way, the call to `get_variable()`
travels recursively through the hierarchy thereby returning directly when a
variable is found. If the variable is not found, `get_variable()` will return
`None`.

A side effect is that a call to get_variable() will not result in an error
when you specify a variable name that does not exist (for example due to a
typo). The reason for this behaviour is that although a variable may not exist
now, it may exist later in the simulation period. For example, as long as
there is no sowing event, there is no crop simulation object and thus
variables of the crop simulation model do not exist. However, after sowing
these variables will be defined and can be found by `get_variable()`.

A way around this is to first check in the VariableKiosk whether a variable
name is registered by calling `variable_exists(<varname>)` on the
VariableKiosk.

.. note::

    In the FORTRAN code of Wageningen crop simulation models it was customary
    to put variable and parameter names in capitals. As a result many crop
    parameter files have parameter names defined in capitals and so have many
    database tables. For this reason, parameters and variable names in PyWOFOST
    are also defined in capitals although it is free to mix upper and lower case
    characters. For convenience `get_variable('MyVar')` both searches for
    'MyVar' as well as 'MYVAR'. 

.. _EventsAndSignals:

Sending and handling signals
----------------------------
During the simulation, some events have to be notified to one or
more parts of the model. Examples of events are:

- sowing, the system has to be notified that the crop simulation
  model has to be initialized.
- finishing of the crop simulation.
- saving the states of the simulation at regular intervals for later output. 
- termination of the entire model simulation.
  
These events can be generated from different parts of the model. Particularly,
crop simulation ending can be signalled by several modules:
1) the harvest date is reached (the AgroManagement module), 2) the crop
reaches physiological maturity (the phenology module) and 3) all leaves die
(leaf dynamics module). Moreover, several
simulation objects  may want to be notified of such an event and in
some cases variables need to be passed with a given event.
It is clear that implementing this kind of communication through classical
variable passing leads to large dependencies between modules which is
undesirable.

Therefore  PyWOFOST uses the pyDispatcher module available from
http://pydispatcher.sourceforge.net/. This module allows a
SimulationObject to send signals and to register one or more
handlers which 'react' to a given signal.


Currently only a limited number of signals are
implemented in the PyWOFOST. These can be found in the `signals` module:

- "CROP_START" when the crop simulations needs to be started.
- "CROP_FINISH" when the crop simulations are finished.
- "TERMINATE" when the entire system needs to be terminated.
- "OUTPUT" when output must be generated for a given set of state variables.

The pyDispatcher module allows for great flexibility and new
events can be easily added to PyWOFOST. Examples are signals for
AgroManagement such as irrigation scheduling and nutrient application
which may be added in future versions of PyWOFOST.

Logging
-------
Besides the features described above, a SimulationObject has a standard
interface for sending log messages through `self.logger` member. This works
using the standard python logging facility, so `self.logger.info(<msg>)`
sends a message with level INFO to the log file.
See the `logging` module in the python documentation for more information.

.. _DrivingVar:

Driving variables
=================

To run the crop simulation meteorological inputs are needed as
driving variables (e.g. the meteorology drives the system from one state to the
next). Internally, PyWOFOST uses daily meteorological variables with the
following characteristics:

====== ========================================================= ===============
Name   Description                                               Unit
====== ========================================================= ===============
TMAX   Daily maximum temperature                                 |C|
TMIN   Daily minimum temperature                                 |C|
VAP    Mean daily vapour pressure                                |hPa|
WIND   Mean daily wind speed at 2 m above ground level           |msec-1|
RAIN   Precipitation (rainfall or water equivalent in case of
       snow or hail).                                            |cmday-1|
IRRAD  Daily global radiation                                    |Jm-2day-1|
LAT    Latitude of the meteorological observations.              |DD|
====== ========================================================= ===============

Strictly speaking the latitude is not a driving variable, but as it is needed
for some of the calculations involving driving variables, it is included with
the driving variables for convenience.

Further there are some meteorological variables which are derived from the
previous ones:

====== ========================================================= ===============
Name   Description                                               Unit
====== ========================================================= ===============
E0     Penman potential evaporation from a free water surface    |cmday-1|
ES0    Penman potential evaporation from a bare soil surface     |cmday-1|
ET0    Penman potential evaporation from a reference crop
       canopy                                                    |cmday-1|
TEMP   Mean daily temperature (TMIN + TMAX)/2                    |C|
DTEMP  Mean daytime temperature (TEMP + TMAX)/2                  |C|
TMINRA The 7-day running average of TMIN                         |C| 
====== ========================================================= ===============


Components with the PyWOFOST package
====================================

The figure below gives an overview of the model components that are currently
available in PyWOFOST. All green components are directly derived from the
WOFOST7.1 source distribution, while the components marked in purple are new
developments, that are based on existing models published in the literature.

On the highest level there are four main components that are composed of one
or more sub-components:
 
1. The water balance which has two sub-components, one for simulations under
   potential and one for free drainage conditions. Moreover, the water balance
   can be combined with the SnowMAUS model for accumulation of snow on the soil
   surface.
2. The crop simulation object which is composed of the many sub-components for
   different processes regarding crop growth such as phenology, assimilation,
   respiration and the dynamics for roots, stems, leaves and storage organs.
   Moreover, components are included for frost damage assessment (FROSTOL,
   CERES_Winterkill) and the estimation of the crown temperature
3. The AgroManagement module which implements management actions such as
   sowing and harvesting
4. The timer module which keeps track of time and generates model
   output in specified intervals (daily, every x days, dekadal, monthly or None)

.. figure:: components.png
    :align: center
    :scale: 45 %
    
    Graphical overview of 
    components (implemented as SimulationObjects) available in the PyWOFOST
    source distribution. The waterbalance "water-limited groundwater" is
    not yet implemented.
    
The PyWOFOST distribution contains a number of additional packages and
modules that are not displayed in the figure as they are utility packages or used
for setting up the environment. These packages are:

  * the `base_classes` module which defines functionality underlying PyWOFOST.
  * the `cabo` package for reading CABO weather and parameter files for 
    retrieving model input.
  * the `db_util` module for communicating with a CGMS database for retrieving
    model input.
  * the `database` package which contains an SQLite database with example data.
    and some utilities for setting up a PyWOFOST database.
  * the `util` module with functions such as penman, angstrom and astro.
  * the `traitlets` module for defining attributes on classes.
  * the `pydispatch` module for sending and handling signals.
  * the `signals` module which defines the used signals.
  * the `test_data` package which defines the test data for some unit tests.
  * the `tests` package which defines the unit tests for many SimulationObjects.


Package structure of PyWOFOST
=============================

The PyWOFOST package structure is:

.. literalinclude:: package_structure.txt




.. |C| replace:: :math:`^{\circ}C`
.. |hPa| replace:: :math:`hPa`
.. |msec-1| replace:: :math:`m sec^{-1}`
.. |cmday-1| replace:: :math:`cm day^{-1}`
.. |Jm-2day-1| replace:: :math:`J m^{-2} day^{-1}`
.. |DD| replace:: :math:`Decimal Degree`