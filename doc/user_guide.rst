.. include:: abbreviations.txt

**************
Design of PCSE
**************


The Python Crop Simulation Environment builds on the heritage
provided by the earlier approaches developed in Wageningen,
notably the Fortran Simulation Environment. The `FSE manual`
(van Kraalingen, 1995) provides
a very good overview on the principles of Euler integration
and its application to crop simulation models. Therefore,
we will not discuss this aspect here.

Nevertheless, PCSE also tries to improve on these approaches
by separating the simulation logic into a number of
distinct components that play a
role in the implementation of (crop) simulation models:

 1. The dynamic part of the simulation is taken care of by a
    dedicated simulation engine which handles the initialization,
    the ordering of rate/state updates as well as keeping
    track of time, retrieving weather data and calling the
    agromanagement module.
 2. Solving the differential equations for crop growth and updating
    the crop state is deferred to SimulationObjects that
    implement dedicated growth or development processes such as
    phenology or |CO2| assimilation.
 3. An AgroManagement module is included which takes care of
    signalling agricultural management actions such as sowing, harvesting,
    irrigation, etc.
 4. Several tools are available for providing weather data and
    reading parameter values from files or databases.

Next, an overview of the different components in PCSE will be provided.

SimulationObjects
=================

PCSE  uses SimulationObjects to group parts of the crop simulation model
that form a logical entity into separate program code sections. In this
way the crop simulation model is grouped into sections that implement certain
biophysical processes such as phenology, assimilation, respiration, etc.
Simulation objects can be grouped to form components that perform the simulation
of an entire crop or a soil profile.

This approach has several advantages:

* Model code with a certain purpose is grouped together, making it easier
  to read, understand and maintain.
* A SimulationObject contains only parameters, rate and state variables
  that are needed. In contrast, with monolythic code it is often unclear (at
  first glance at least) what biophysical process they belong to.
* Isolation of process implementations creates less dependencies, but more
  importantly, dependencies are evident from the code which makes it easier
  to modify individual SimulationObjects.
* SimulationObjects can be tested individually by comparing output vs the
  expected output (e.g. unit testing).
* SimulationObjects can be exchanged for other objects with the same purpose
  but a different biophysical approach. For example, the WOFOST assimilation
  approach could be easily replaced by a more simple Light Use Efficiency or
  Water Use Efficiency approach, by replacing the SimulationObject that
  handles the |CO2| assimilation.

Characteristics of SimulationObjects
------------------------------------

Each SimulationObject is defined in the same way and has a couple of standard
sections and methods which facilitates understanding and readability.
Each SimulationObject has parameters to define the mathematical relationships,
it has state variables to define the state of the system and it has rate
variables that describe the rate of change from one time step to the next.
Moreover, a SimulationObject may contain other SimulationObjects that
together form a logical structure. Finally, the SimulationObject must implemented
separate code sections for initialization, rate calculation and integration
of the rates of change. A finalization step which is called at the end of the simulation
can be added optionally. The skeleton of a SimulationObject looks like this:

.. code-block:: python

    class CropProcess(SimulationObject):

        class Parameters(ParamTemplate):
            PAR1 = Float()
            # more parameters defined here

        class StateVariables(StatesTemplate):
            STATE1 = Float()
            # more state variables defined here

        class RateVariables(RatesTemplate):
            RATE1 = Float()
            # more rate variables defined here

        def initialize(day, kiosk, parametervalues):
            """Initializes the SimulationObject with given parametervalues."""

        @prepare_rates
        def calc_rates(day, drv):
            """Calculate the rates of change given the current states and driving
            variables (drv)."""

        @prepare_states
        def integrate(day, delt):
            """Integrate the rates of change on the current state variables
            multiplied by the time-step
            """

        @prepare_states
        def finalize(day):
            """do some final calculations when the simulation is finishing."""


The strict separation of program logic was copied from the Fortran Simulation
Environment (FSE, Rappoldt and Van Kraalingen 1996. http://edepot.wur.nl/4411
and Van Kraalingen 1995, http://edepot.wur.nl/35555) and
is critical to ensure that the simulation results are correct.
The different calculations types (integration, driving variables and
rate calculations) should be strictly separated. In other words, first all
states should be updated, subsequently all driving variables should be calculated,
after which all rates of change should be calculated. If this rule is not
applied rigorously, some rates may pertain to states at
the current time whereas others will pertain to states from the previous time
step (Van Kraalingen, 1995).

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
routine calls with a large risk of mapping names onto the wrong variable
in the calling routine. Moreover, the variables in the call are often a
mixture of model parameters, driving variables, rate/state variables and variables that are
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


In PCSE the communication between
SimulationObjects is taken care of by the so-called `VariableKiosk`. The
metaphore kiosk is used because the SimulationObjects publish
their rate and/or state variables (or a subset) into the kiosk, other
SimulationObjects can subsequently request the variable value from the kiosk
without any knowledge about the SimulationObject that published it.
Therefore, the VariableKiosk is shared by all SimulationObjects and must
be provided when SimulationObjects initialize.

The VariableKiosk
-----------------

The VariableKiosk is an integral part of PCSE and it has several
important functions. First of all it takes
care of communicating variables between SimulationObjects.
When a variable has been published,
PCSE updates the value of that variables into the `VariableKiosk` at each
simulation cycle. Under the hood, this functionality is provided by the
`traitlets` module from which all SimulationObjects are derived.

To avoid that state or rate variables of previous time steps keep
lagging in the VariableKiosk, the kiosk contents kiosk are flushed during each
time step. After all rates have been calculated, the values of all state
variables are flushed from the VariableKiosk. Similarly, after the update of the
state variables, the values of rate variables are flushed from the kiosk.

Second, to enforce that variable names are unique across the entire model
the VariableKiosk registers all variables across all SimulationObjects.
Variable names that are defined in a `RateVariables` or
`StateVariables` class definition are automatically
registered in the VariableKiosk and the name is checked for uniqueness.
Moreover, it is tracked which SimulationObject registers which variable and
only the SimulationObject that registers and publishes a variable can change
its value. All other objects can retrieve that value, but not change it.
Uniqueness of variable names is important for retrieving their value; if two
SimulationObject would define the same variable than PCSE would only find the
first one. It is not even guaranteed that this would be the same
variable between sessions.

Finally, the VariableKiosk serves as a unique identifier for a PCSE model instances.
Since the VariableKiosk is shared across all SimulationObjects within the
hierarchy, it can be used as a unique identifier of a PCSE model instance
by retrieving its identifier with python's `id()` function. This aspect is
important when running ensembles of PCSE models in combination with
sending and receiving signals (see :ref:`EventsAndSignals`).

Logging
-------
Besides the features described above, a SimulationObject has a standard
interface for sending log messages through `self.logger`. This works
using the standard python logging facility, so `self.logger.info(<msg>)`
sends a message with level INFO to the log file.

The Engine
==========

The PCSE Engine provides the environment where SimulationObjects are 'living'.
The engine takes care of reading the model configuration, initializing model
components (e.g. groups of SimulationObjects), driving the simulation
forward by calling the SimulationObjects, calling the agromanagement
unit, keeping track of time and providing the weather data needed.

Configuration files
-------------------

The configuration of a model in PCSE is read from a configuration file
by the Engine. This configuration defines the following aspects of
a simulation:

* the component to be used for simulation the soil dynamics;
* the component to be used for simulating the crop dynamics;
* the component to be used for agromanagement actions;
* the names of the variables to be stored for output;
* the frequency when output is generated;
* the names of the variables to be stored for summary output. Summary
  output is only generated when the crop simulation is finished.

.. _ContinuousSimulation:

Continuous simulation in PCSE
-----------------------------

To implement continuous simulation, the engine in PCSE uses the same approach as
FSE: Euler integration with a fixed time step of one day.  The following
figures shows the principle of continuous simulation
and the execution order of various steps.

.. figure:: continuous_simulation.png
    :align: center
    :scale: 50%

    Order of calculations for continuous simulation using Euler integration
    (after Van Kraalingen, 1995).

The steps in the proces cycle that are evident from the figure above are
implemented in the simulation `Engine` which is completely separated
from the model logistic itself. Therefore the Engine is generic and can be
used for any model that is defined in PCSE.

From the previous figure it is clear that before the simulation can start the
system has to be initialized:

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

Tracking time
-------------

PCSE contains a dedicated timer module which keeps track of time during the
simulation. Since PCSE runs on fixed timesteps of one day, this is a fairly
simple module which uses python `datetime.date` and `datetime.timedelta`
objects to save and update the simulation time. Additionally, the timer
module generates output signals that indicate when simulation results must be
stored. Depending on the configuration, output signals can be generated
every *n* days or at the last day of each month/dekad.

Retrieving weather data
-----------------------

Weather data must be provided to the `Engine` using a special component called
a `WeatherDataProvider`. The WeatherDataProvider encapsulates the weather
data and provides an interface for the Engine to retrieve the weather data
for a given day. Within PCSE several WeatherDataProviders are available
that can retrieve data from different sources.

Agromanagement
==============

Agromanagement is a specific part of PCSE because it does not deal
with continuous processed but rather with events that signal dedicated actions
such as sowing, harvesting, irrigation, etc. Such events are a kind of
discontinuities in the simulation and they are difficult to implement within
the normal processing loop for continuous simulation
(see :ref:`ContinuousSimulation`). Therefore, PCSE uses a different for
dealing with agromanagement actions; it uses signals that are broadcasted
through the entire simulation environment.

The purpose of the agromanagement module is to check at each time-step
if management actions are scheduled. If so, the signals that are defined
for these actions are broadcasted by the agromanagement module together
with the parameters that are associated with the particular action.

Currently, PCSE contains one agromanagement component `AgroManagementSingleCrop`
which is used for simulating a single cropping season. New `AgroManagement`
modules will be developed which will allow continuous simulation with
rotation of crops and the support for irrigation and nutrient applications.

.. _EventsAndSignals:

Sending and handling signals
============================

As already note before, during the simulation, events may occur that
have to be notified to one or more parts of the model. Clear examples of
such events are described above in the section on agromanagement. However,
agromanagent is not the only component where signals and events play an
important role. For example OUTPUT signals are generated by the `timer`
module to indicate that the current state of variables should be saved
for later use. Moreover, ending the crop simulation can be signalled by
several modules:
1) the harvest date is reached (the `AgroManagement` module), 2) the crop
reaches physiological maturity (the `Phenology` module) and 3) all leaves die
(the `Leaf_dynamics` module). Moreover, several
simulation objects  may need to be notified of such an event and in
some cases variables need to be passed with a given event.
It is clear that implementing this kind of communication through classical
variable passing leads to large dependencies between modules which is
undesirable.

Therefore PCSE uses the pyDispatcher module available from
http://pydispatcher.sourceforge.net/. This module allows a SimulationObject
to send signals and to register one or more handlers which 'react' to a
given signal. Using this approach information can be shared efficiently
between SimulationObjects, the Engine and the AgroManagement module that
would otherwise be difficult to implement.

The signals that are used by PCSE are all defined in the `signals` module.
Currently, the number of distinct signals is still quite limited and consist of:

* CROP_START: which is used when a crop simulation is started. The most
  important parameter is the crop simulation component itself which send
  to the Engine in order become part of the simulation loop.
* CROP_FINISH: when the crop simulation is finished.
* TERMINATE":when the simulation needs to be halted.
* OUTPUT: when the current state of selected variables must be stored.

.. _DrivingVar:

Driving variables
=================

To run the crop simulation meteorological inputs are needed as
driving variables (e.g. the meteorology drives the system from one state to the
next) which are provided by a `WeatherDataProvider`.
PCSE uses the following daily meteorological variables:

========= ========================================================= ===============
Name        Description                                               Unit
========= ========================================================= ===============
TMAX      Daily maximum temperature                                  |C|
TMIN      Daily minimum temperature                                  |C|
VAP       Mean daily vapour pressure                                 |hPa|
WIND      Mean daily wind speed at 2 m above ground level            |msec-1|
RAIN      Precipitation (rainfall or water equivalent in case of
          snow or hail).                                             |cmday-1|
IRRAD     Daily global radiation                                     |Jm-2day-1|
SNOWDEPTH Depth of snow cover (optional)                             |cm|
========= ========================================================= ===============

The snow depth is an optional meteorological variable and is only used for
estimating the impact of frost damage on the crop (if enabled). Snow depth can
also be simulated by the `SnowMAUS` module if observations are not available
on a daily basis. Furthermore there are some meteorological variables which
are derived from the previous ones:

====== ========================================================= ===============
Name   Description                                                 Unit
====== ========================================================= ===============
E0     Penman potential evaporation for a free water surface      |cmday-1|
ES0    Penman potential evaporation for a bare soil surface       |cmday-1|
ET0    Penman or Penman-Monteith potential evaporation
       for a reference crop canopy                                |cmday-1|
TEMP   Mean daily temperature (TMIN + TMAX)/2                     |C|
DTEMP  Mean daytime temperature (TEMP + TMAX)/2                   |C|
TMINRA The 7-day running average of TMIN                          |C|
====== ========================================================= ===============


