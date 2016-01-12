.. include:: abbreviations.txt

******************
Understanding PCSE
******************

An overview of PCSE
===================

The Python Crop Simulation Environment builds on the heritage
provided by the earlier approaches developed in Wageningen,
notably the Fortran Simulation Environment. The `FSE manual <http://edepot.wur.nl/35555>`_
(van Kraalingen, 1995) provides a very good overview on the principles of Euler integration
and its application to crop simulation models. Therefore, we will not discuss this in detail
here.

Nevertheless, PCSE also tries to improve on these approaches
by separating the simulation logic into a number of
distinct components that play a role in the implementation of (crop)
simulation models:

 1. The dynamic part of the simulation is taken care of by a
    dedicated simulation `Engine` which handles the initialization,
    the ordering of rate/state updates for the soil and plant
    modules as well as keeping track of time, retrieving weather data and
    calling the agromanager module.
 2. Solving the differential equations for soil/plant system and updating
    the model state is deferred to SimulationObjects that
    implement (bio)physical processes such as phenology or |CO2| assimilation.
 3. An AgroManager module is included which takes care of
    signalling agricultural management actions such as sowing, harvesting,
    irrigation, etc.
 4. Communication between PCSE components is implemented by either exporting
    variables into a shared state object or by implementing signals that can be
    broadcasted and received by any PCSE object.
 5. Several tools are available for providing weather data and
    reading parameter values from files or databases.

Next, an overview of the different components in PCSE will be provided.

The Engine
==========

The PCSE Engine provides the environment where the simulation takes place.
The engine takes care of reading the model configuration, initializing model
components, driving the simulation
forward by calling the SimulationObjects, calling the agromanagement
unit, keeping track of time, providing the weather data needed and
storing the model variables during the simulation for later output.
The Engine itself is generic and can be used for any model that is defined
in PCSE.

.. _ContinuousSimulation:

.. Continuous simulation in PCSE
.. -----------------------------

To implement continuous simulation, the engine uses the same approach as
FSE: Euler integration with a fixed time step of one day.  The following
figures shows the principle of continuous simulation
and the execution order of various steps.

.. figure:: continuous_simulation.png
   :align: center
   :width: 500 px

   Order of calculations for continuous simulation using Euler integration
   (after Van Kraalingen, 1995).

The steps in the process cycle that are shown in the figure above are
implemented in the simulation `Engine` which is completely separated
from the model logic itself. Moreover, it demonstrates that before
the simulation can start the engine has to be initialized which involves
several steps:

1. The model configuration must be loaded;
2. The AgroManager module must be initialized and called to determine
   the first and last of the simulation sequence;
3. The timer must be initialized with the first and last day of the
   simulation sequence;
4. The soil component specified in the model configuration must be
   initialized.
5. The weather variables must be retrieved for the starting day;
6. The AgroManager must be called to trigger any management events that
   are scheduled for the starting day.
7. The initial rates of change based on the initial states and driving
   variables must be calculated;
8. Finally, output can be collected to save the initial states and rates of
   the simulation.

The next cycle in the simulation will now start with an update of the timer to
the next time step (e.g. day). Next, the rates of change of
the previous day will be integrated onto the state variables and the driving
variables for the current day will be retrieved. Finally, the rates of change
will be recalculated based on the new driving variables and updated model
states and so forth.

The simulation loop will terminate when some finish condition has been reached.
Usually, the `AgroManager` module will encounter the end of the agricultural
campaign and will issue a terminate signal that terminates the entire simulation.

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
together form a logical structure. Finally, the SimulationObject must implement
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
            self.params = self.Parameters(parametervalues)
            self.rates = self.RateVariables(kiosk)
            self.states = self.StateVariables(kiosk, STATE1=0.)

        @prepare_rates
        def calc_rates(day, drv):
            """Calculate the rates of change given the current states and driving
            variables (drv)."""

            # simple example of rate calculation using rainfall (drv.RAIN)
            self.rates.RATE1 = self.params.PAR1 * drv.RAIN

        @prepare_states
        def integrate(day, delt):
            """Integrate the rates of change on the current state variables
            multiplied by the time-step
            """
            self.states.STATE1 += self.rates.RATE1 * delt

        @prepare_states
        def finalize(day):
            """do some final calculations when the simulation is finishing."""


The strict separation of program logic was copied from the Fortran Simulation
Environment (FSE, `Rappoldt and Van Kraalingen 1996 <http://edepot.wur.nl/4411>`_
and `Van Kraalingen 1995 <http://edepot.wur.nl/35555>`_) and
is critical to ensure that the simulation results are correct.
The different calculations types (integration, driving variables and
rate calculations) should be strictly separated. In other words, first all
states should be updated, subsequently all driving variables should be calculated,
after which all rates of change should be calculated. If this rule is not
applied rigorously, some rates may pertain to states at
the current time whereas others will pertain to states from the previous time
step .

Compared to the FSE system and the
`FORTRAN implementation of WOFOST <https://github.com/ajwdewit/wofost>`_,
the `initialize()`, `calc_rates()`, `integrate()` and `finalize()` sections
match with the *ITASK* numbers 1, 2, 3, 4.

Communication between SimulationObjects
---------------------------------------

A complicating factor that arises when using modular code is how to arrange
the communication between SimulationObjects. For example, the `evapotranspiration`
SimulationObject will need information about the leaf area index from the
`leaf_dynamics` SimulationObject to calculate the crop transpiration
values.

In PCSE the communication between
SimulationObjects is taken care of by the so-called `VariableKiosk`. The
metaphore kiosk is used because the SimulationObjects publish
their rate and/or state variables (or a subset) into the kiosk, other
SimulationObjects can subsequently request the variable value from the kiosk
without any knowledge about the SimulationObject that published it.
Therefore, the VariableKiosk is shared by all SimulationObjects and must
be provided when SimulationObjects initialize. See the section on communication
between PCSE components for a detailed description of the variable kiosk.

The AgroManager
===============

Agromanagement is an intricate part of PCSE which is needed for
simulating the processes that are happening
on agriculture fields. In order for crops to growth, farmers must first plow the
fields and sow the crop. Next, they have to do proper management including
irrigation, weeding, nutrient application, pest control and finally harvesting.
All these actions have to be scheduled at specific dates, connected to certain
crop stages or in dependence of soil and weather conditions. Moreover specific
parameters such as the amount of irrigation or nutrients must be provided as well.

In previous versions of WOFOST, the options for agromanagement were limited to
sowing and harvesting. On the one had this was because agromangement was often assumed
to be optimal and thus there was little need for detailed agromanagement.
On the other hand, implementing agromanagement is relatively complex because
agromanagement consists of events that are happening rather then
continuous processes. As such, it does not fit well in the traditional simulation
cycle, see :ref:`ContinuousSimulation`.

Also from a technical point of view implementing such events through the traditional
function calls for rate calculation and state updates is not attractive. For
example, for indicating a nutrient application event several additional parameters
would have to be passed: e.g. the type of nutrient, the amount and its efficiency.
This has several drawbacks, first of all, only a limited number of SimulationObjects
will actually do something with this information while for all other objects, the
information is of no use. Second, nutrient application will usually happen only once
or twice in the growing cycle. So for a 200-day growing cycle there will be
198 days where the parameters do not carry any information. Nevertheless, they
would still be present in the function call, thereby decreasing the computational
efficiency and the readability of the code. Therefore, PCSE uses a very different
approach for agromanagement events which is based on signals (see XXX).

Defining agromanagement in PCSE
-------------------------------

Defining the agromanagement in PCSE is not very complicated and first starts with
defining a sequence of campaigns. Campaigns
start on a prescribed calendar date and finalize when the next campaign starts.
Each campaign is characterized by zero or one crop calendar, zero or more timed events and
zero or more state events. The crop calendar specifies the timing of the crop (sowing,
harvesting) while the timed and state events can be used to specify management actions
that are either dependent on time (a specific date) or a certain model state variable
such as crop development stage. Crop calendars and event definitions are only valid
for the campaign in which they are defined.

The data format used for defining the agromanagement in PCSE is called YAML. YAML is a
versatile format optimized for readability by humans while still having the power of XML.
However, the agromanagement definition in PCSE is by no means tied to YAML and can be
read from a database as well.

The structure of the data needed as input for the AgroManager is most easily understood
with an example (below). The example definition consists of three campaigns, the first
starting on 1999-08-01, the second starting on 2000-09-01 and the last campaign starting
on 2001-03-01. The first campaign consists of a crop calendar for winter-wheat starting
with sowing at the given crop_start_date. During the campaign there are timed events for
irrigation at 2000-05-25 and 2000-06-30. Moreover, there are state events for  fertilizer
application (event_signal: apply_npk) given by development stage (DVS) at DVS 0.3, 0.6 and 1.12.

The second campaign has no crop calendar, timed events or state events. This means that
this is a period of bare soil with only the water balance running. The third campaign is
for fodder maize sown at 2001-04-15 with two series of timed events (one for irrigation and
one for N/P/K application) and no state events. The end date of the simulation in this case
will be 2001-11-01 (2001-04-15 + 200 days).

An example of an agromanagement definition file::

    AgroManagement:
    - 1999-08-01:
        CropCalendar:
            crop_id: winter-wheat
            crop_start_date: 1999-09-15
            crop_start_type: sowing
            crop_end_date:
            crop_end_type: maturity
            max_duration: 300
        TimedEvents:
        -   event_signal: irrigate
            name:  Timed irrigation events
            comment: All irrigation amounts in cm
            events_table:
            - 2000-05-25: {irrigation_amount: 3.0}
            - 2000-06-30: {irrigation_amount: 2.5}
        StateEvents:
        -   event_signal: apply_npk
            event_state: DVS
            zero_condition: rising
            name: DVS-based N/P/K application table
            comment: all fertilizer amounts in kg/ha
            events_table:
            - 0.3: {N_amount : 1, P_amount: 3, K_amount: 4}
            - 0.6: {N_amount: 11, P_amount: 13, K_amount: 14}
            - 1.12: {N_amount: 21, P_amount: 23, K_amount: 24}
    - 2000-09-01:
        CropCalendar:
        TimedEvents:
        StateEvents
    - 2001-03-01:
        CropCalendar:
            crop_id: fodder-maize
            crop_start_date: 2001-04-15
            crop_start_type: sowing
            crop_end_date:
            crop_end_type: maturity
            max_duration: 200
        TimedEvents:
        -   event_signal: irrigate
            name:  Timed irrigation events
            comment: All irrigation amounts in cm
            events_table:
            - 2001-06-01: {irrigation_amount: 2.0}
            - 2001-07-21: {irrigation_amount: 5.0}
            - 2001-08-18: {irrigation_amount: 3.0}
            - 2001-09-19: {irrigation_amount: 2.5}
        -   event_signal: apply_npk
            name:  Timed N/P/K application table
            comment: All fertilizer amounts in kg/ha
            events_table:
            - 2001-05-25: {N_amount : 50, P_amount: 25, K_amount: 22}
            - 2001-07-05: {N_amount : 70, P_amount: 35, K_amount: 32}
        StateEvents:

Crop calendars
--------------

The crop calendar definition will be passed on to a `CropCalendar` object which is is
responsible for storing, checking, starting and ending the crop cycle during a PCSE simulation.
At each time step the instance of `CropCalendar` is called
and at the dates defined by its parameters it initiates the appropriate actions:

- sowing/emergence: A `crop_start` signal is dispatched including the parameters needed to
  start the new crop simulation object (crop_start_type and crop_id)
- maturity/harvest: the crop cycle is ended by dispatching a `crop_finish` signal with the
  appropriate parameters.

For a detailed description of a crop calendar see the code documentation on the CropCalendar in the
section on :ref:`Agromanagement`.

Timed events
------------

Timed events are management actions that are occurring on specific dates. As simulations in PCSE run
on daily time steps it is easy to schedule actions on dates. Time events are characterized by
an event signal, a name and comment that can be used to describe the event and finally an
events table that lists the dates for the events and the parameters that should be passed onward.

For a detailed description of a timed events see the code documentation on the TimedEventsDispatcher
in the section on :ref:`Agromanagement`.

State events
------------

State events are management actions that are tied to certain model states. Examples are actions such
as nutrient application that should be executed at certain crop stages, or irrigation application
that should occur only when the soil is dry. PCSE has a flexible definition of state events and an event
can be connected to any variable that is defined within PCSE.

Each state event is defined by an `event_signal`, an `event_state` (e.g. the model
state that triggers the event) and a `zero condition`. Moreover, an optional name and an
optional comment can be provided. Finally the events_table specifies at which model state values
the event occurs. The events_table is a list which provides for each state the parameters that
should be dispatched with the given event_signal.

Managing state events is more complicated than timed events because PCSE cannot determine beforehand at
which time step these events will trigger.
For finding the time step at which a state event occurs PCSE uses the concept of `zero-crossing`.
This means that a state event is triggered when (`model_state` - `event_state`) equals or
crosses zero. The `zero_condition` defines how this crossing should take place. The value of
`zero_condition` can be:

* `rising`: the event is triggered when (`model_state` - `event_state`) goes from a negative value towards
   zero or a positive value.
* `falling`: the event is triggered when (`model_state` - `event_state`) goes from a positive value towards
   zero or a negative value.
* `either`: the event is triggered when (`model_state` - `event_state`) crosses or reaches zero from any
   direction.


