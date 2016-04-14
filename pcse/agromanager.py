# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), Juli 2015
# from __future__ import print_function
"""Implementation of AgroManager and related classes for agromanagement actions in PCSE.

Available classes:

  * CropCalendar: A class for handling cropping calendars
  * TimedEventDispatcher: A class for handling timed events (e.g. events connected to a date)
  * StateEventDispatcher: A class for handling state events (e.g. events that happen when a state variable reaches
    a certain values.
  * AgroManager: A class for handling all agromanagement events which encapsulates
    the CropCalendar and Timed/State events.
"""

from datetime import date, timedelta
import logging
from collections import Counter

from .base_classes import DispatcherObject, VariableKiosk, SimulationObject, ParameterProvider, AncillaryObject
from .traitlets import HasTraits, Float, Int, Instance, Enum, Bool, List, Dict, Unicode
from . import exceptions as exc
from .util import ConfigurationLoader
from . import signals
from . import exceptions as exc

def cmp2(x, y):
    """
    Compare two values and return sign

    Surrogate for cmp() function in Python2
    """
    return (x > y) - (x < y)

def check_date_range(day, start, end):
    """returns True if start <= day < end

    Optionally, end may be None. in that case return True if start <= day

    :param day: the date that will be checked
    :param start: the start date of the range
    :param end: the end date of the range or None
    :return: True/False
    """

    if end is None:
        return start <= day
    else:
        return start <= day < end


class CropCalendar(HasTraits, DispatcherObject):
    """A crop calendar for managing the crop cycle.

    A `CropCalendar` object is responsible for storing, checking, starting and ending
    the crop cycle. The crop calendar is initialized by providing the parameters needed
    for defining the crop cycle. At each time step the instance of `CropCalendar` is called
    and at dates defined by its parameters it initiates the appropriate actions:

    - sowing/emergence: A `crop_start` signal is dispatched including the parameters needed to
      start the new crop simulation object
    - maturity/harvest: the crop cycle is ended by dispatching a `crop_finish` signal with the
      appropriate parameters.

    :param kiosk: The PCSE VariableKiosk instance
    :param crop_id: String identifying the crop
    :param crop_start_date: Start date of the crop simulation
    :param crop_start_type: Start type of the crop simulation ('sowing', 'emergence')
    :param crop_end_date: End date of the crop simulation
    :param crop_end_type: End type of the crop simulation ('harvest', 'maturity', 'earliest')
    :param max_duration: Integer describing the maximum duration of the crop cycle

    :return: A CropCalendar Instance
    """

    # Characteristics of the crop cycle
    crop_id = Unicode()
    crop_start_date = Instance(date)
    crop_start_type = Enum(["sowing", "emergence"])
    crop_end_date = Instance(date)
    crop_end_type = Enum(["maturity", "harvest", "earliest"])
    max_duration = Int()

    # system parameters
    kiosk = Instance(VariableKiosk)
    parameterprovider = Instance(ParameterProvider)
    mconf = Instance(ConfigurationLoader)
    logger = Instance(logging.Logger)

    # Counter for duration of the crop cycle
    duration = Int(0)
    in_crop_cycle = Bool(False)

    def __init__(self, kiosk, crop_id=None, crop_start_date=None,
                 crop_start_type=None, crop_end_date=None, crop_end_type=None, max_duration=None):

        # set up logging
        loggername = "%s.%s" % (self.__class__.__module__,
                                self.__class__.__name__)

        self.logger = logging.getLogger(loggername)
        self.kiosk = kiosk
        self.crop_id = crop_id
        self.crop_start_date = crop_start_date
        self.crop_start_type = crop_start_type
        self.crop_end_date = crop_end_date
        self.crop_end_type = crop_end_type
        self.max_duration = max_duration

        self._connect_signal(self._on_CROP_FINISH, signal=signals.crop_finish)

    def validate(self, campaign_start_date, next_campaign_start_date):
        """Validate the crop calendar internally and against the interval for
        the agricultural campaign.

        :param campaign_start_date: start date of this campaign
        :param next_campaign_start_date: start date of the next campaign
        """

        # Check that crop_start_date is before crop_end_date
        crop_end_date = self.crop_end_date
        if self.crop_end_type == "maturity":
            crop_end_date = self.crop_start_date + timedelta(days=self.max_duration)
        if self.crop_start_date >= crop_end_date:
            msg = "crop_end_date before or equal to crop_start_date for crop '%s'!"
            raise exc.PCSEError(msg % (self.crop_start_date, self.crop_end_date))

        # check that crop_start_date is within the campaign interval
        r = check_date_range(self.crop_start_date, campaign_start_date, next_campaign_start_date)
        if r is not True:
            msg = "Start date (%s) for crop '%s' not within campaign window (%s - %s)." % \
                  (self.crop_start_date, self.crop_id, campaign_start_date, next_campaign_start_date)
            raise exc.PCSEError(msg)

    def __call__(self, day):
        """Runs the crop calendar to determine if any actions are needed.

        :param day:  a date object for the current simulation day
        :param drv: the driving variables at this day
        :return: None
        """

        if self.in_crop_cycle:
            self.duration += 1

        # Start of the crop cycle
        if day == self.crop_start_date:  # Start a new crop
            self.duration = 0
            self.in_crop_cycle = True
            msg = "Starting crop (%s) on day %s" % (self.crop_id, day)
            self.logger.info(msg)
            self._send_signal(signal=signals.crop_start, day=day,
                              crop_id=self.crop_id, crop_start_type=self.crop_start_type,
                              crop_end_type=self.crop_end_type)

        # end of the crop cycle
        finish_type = None
        # Check if crop_end_date is reached for CROP_END_TYPE harvest/earliest
        if self.crop_end_type in ["harvest", "earliest"]:
            if day == self.crop_end_date:
                finish_type = "harvest"

        # Check for forced stop because maximum duration is reached
        if self.in_crop_cycle and self.duration == self.max_duration:
            finish_type = "max_duration"

        # If finish condition is reached send a signal to finish the crop
        if finish_type is not None:
            self.in_crop_cycle = False
            self._send_signal(signal=signals.crop_finish, day=day,
                              finish=finish_type, crop_delete=True)

    def _on_CROP_FINISH(self):
        """Register that crop has reached the end of its cycle.
        """
        self.in_crop_cycle = False

    def get_end_date(self):
        """Return the end date of the crop cycle.

        This is either given as the harvest date or calculated as
        crop_start_date + max_duration

        :return: a date object
        """
        if self.crop_end_type in ["harvest", 'earliest']:
            return self.crop_end_date
        else:
            return self.crop_start_date + timedelta(days=self.max_duration)

    def get_start_date(self):
        """Returns the start date of the cycle. This is always self.crop_start_date

        :return: the start date
        """
        return self.crop_start_date


class TimedEventsDispatcher(HasTraits, DispatcherObject):
    """Takes care handling events that are connected to a date.

    Events are handled by dispatching a signal (taken from the `signals` module)
    and providing the relevant parameters with the signal. TimedEvents can be
    most easily understood when looking at the definition in the agromanagement
    file. The following section (in YAML) provides the definition of two instances
    of TimedEventsDispatchers::

        TimedEvents:
        -   event_signal: irrigate
            name:  Timed irrigation events
            comment: All irrigation amounts in mm
            events_table:
            - 2000-01-01: {irrigation_amount: 20}
            - 2000-01-21: {irrigation_amount: 50}
            - 2000-03-18: {irrigation_amount: 30}
            - 2000-03-19: {irrigation_amount: 25}
        -   event_signal: apply_npk
            name:  Timed N/P/K application table
            comment: All fertilizer amounts in kg/ha
            events_table:
            - 2000-01-10: {N_amount : 10, P_amount: 5, K_amount: 2}
            - 2000-01-31: {N_amount : 30, P_amount: 15, K_amount: 12}
            - 2000-03-25: {N_amount : 50, P_amount: 25, K_amount: 22}
            - 2000-04-05: {N_amount : 70, P_amount: 35, K_amount: 32}

    Each TimedEventDispatcher is defined by an `event_signal`, an optional name,
    an optional comment and the events_table. The events_table is list which provides
    for each date the parameters that should be dispatched with the given
    event_signal.
    """
    event_signal = None
    events_table = List()
    days_with_events = Instance(Counter)
    kiosk = Instance(VariableKiosk)
    logger = Instance(logging.Logger)
    name = Unicode()
    comment = Unicode()

    def __init__(self, kiosk, event_signal, name, comment, events_table):
        """Initialising a TimedEventDispatcher

        :param kiosk: an instance of the VariableKiosk
        :param event_signal: the signal to be dispatched when the event occurs (from pcse.signals)
        :param name: the name of the event dispatcher
        :param comment: A comment that will be used in log message
        :param events_table: The events table, the structure here is a list of dicts, with each dict having only
            one key/value with the key being the date of the event and the value a dict of parameter values
            that should be dispatched with the signal.
        """

        # set up logging
        loggername = "%s.%s" % (self.__class__.__module__,
                                self.__class__.__name__)
        self.logger = logging.getLogger(loggername)

        self.kiosk = kiosk
        self.events_table = events_table
        self.name = name
        self.comment = comment

        # get signal from signals module
        if not hasattr(signals, event_signal):
            msg = "Signal '%s'  not defined in pcse.signals module."
            raise exc.PCSEError(msg % event_signal)
        # self.event_signal = getattr(signals, event_signal)
        self.event_signal = getattr(signals, event_signal)

        # Build a counter for the days with events.
        self.days_with_events = Counter()
        for ev in self.events_table:
            self.days_with_events.update(ev.keys())

        # Check if there are days with two or more events under the
        # same signal which is not allowed.
        multi_days = []
        for day, count in self.days_with_events.items():
            if count > 1:
                multi_days.append(day)
        if multi_days:
            msg = "Found days with more than 1 event for events table '%s' on days: %s"
            raise exc.PCSEError(msg % (self.name, multi_days))

    def validate(self, campaign_start_date, next_campaign_start_date):
        """Validates the timed events given the campaign window

        :param campaign_start_date: Start date of the campaign
        :param next_campaign_start_date: Start date of the next campaign, can be None
        """
        for event in self.events_table:
            day = event.keys()[0]
            r = check_date_range(day, campaign_start_date, next_campaign_start_date)
            if r is not True:
                msg = "Timed event at day %s not in campaign interval (%s - %s)" %\
                      (day, campaign_start_date, next_campaign_start_date)
                raise exc.PCSEError(msg)

    def __call__(self, day):
        """Runs the TimedEventDispatcher to determine if any actions are needed.

        :param day: a date object for the current simulation day
        :return: None
        """
        if day not in self.days_with_events:
            return

        for event in self.events_table:
            if day in event:
                msg = "Time event dispatched from '%s' at day %s" % (self.name, day)
                self.logger.info(msg)
                kwargs = event[day]
                self._send_signal(signal=self.event_signal, **kwargs)

    def get_end_date(self):
        """Returns the last date for which a timed event is given
        """
        return max(self.days_with_events)


class StateEventsDispatcher(HasTraits, DispatcherObject):
    """Takes care handling events that are connected to a model state variable.

    Events are handled by dispatching a signal (taken from the `signals` module)
    and providing the relevant parameters with the signal. StateEvents can be
    most easily understood when looking at the definition in the agromanagement
    file. The following section (in YAML) provides the definition of two instances
    of StateEventsDispatchers::

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
        -   event_signal: irrigate
            event_state: SM
            zero_condition: falling
            name: Soil moisture driven irrigation scheduling
            comment: all irrigation amounts in cm of water
            events_table:
            - 0.15: {irrigation_amount: 20}


    Each StateEventDispatcher is defined by an `event_signal`, an `event_state` (e.g. the model
    state that triggers the event) and a `zero condition`. Moreover, an optional name and an
    optional comment can be provided. Finally the events_table specifies at which model state values
    the event occurs. The events_table is a list which provides for each state the parameters that
    should be dispatched with the given event_signal.

    For finding the time step at which a state event occurs PCSE uses the concept of `zero-crossing`.
    This means that a state event is triggered when (`model_state` - `event_state`) equals or
    crosses zero. The `zero_condition` defines how this crossing should take place. The value of `zero_condition`
    can be:

    * `rising`: the event is triggered when (`model_state` - `event_state`) goes from a negative value towards
       zero or a positive value.
    * `falling`: the event is triggered when (`model_state` - `event_state`) goes from a positive value towards
       zero or a negative value.
    * `either`: the event is triggered when (`model_state` - `event_state`) crosses or reaches zero from any
       direction.

    The impact of the zero_condition can be illustrated using the example definitions above.
    The development stage of the crop (DVS) only increases from 0 at emergence to 2 at maturity. A StateEvent
    set on the DVS (first example) will therefore logically have a zero_condition 'rising' although 'either'
    could be used as well. A DVS-based event will not occur with zero_condition set to 'falling' as the value
    of DVS will not decrease.

    The soil moisture (SM) however can both increase and decrease. A StateEvent for applying irrigation (second
    example) will therefore be specified with a zero_condition 'falling' because the event must be triggered
    when the soil moisture level reaches or crosses the minimum level specified by the events_table. Note that
    if we set the zero_condition to 'either' the event would probably occur again the next time-step because
    the irrigation amount increase the soil moisture and (`model_state` - `event_state`) crosses zero again
    but from the other direction.
    """
    event_signal = None
    event_state = Unicode()
    zero_condition = Enum(['rising', 'falling', 'either'])
    events_table = List()
    kiosk = Instance(VariableKiosk)
    logger = Instance(logging.Logger)
    name = Unicode()
    comment = Unicode()
    previous_signs = List()

    def __init__(self, kiosk, event_signal, event_state, zero_condition, name,
                 comment, events_table):
        """Initialising a StateEventDispatcher

        :param kiosk: an instance of the VariableKiosk
        :param event_signal: the signal to be dispatched when the event occurs (from pcse.signals)
        :param event_state: the name of the state variable that should trigger the event
        :param zero_condition: the zero_condition, one of 'rising'|'falling'|'either'
        :param name: the name of the event dispatcher
        :param comment: A comment that will be used in log message
        :param events_table: The events table, the structure here is a list of dicts, with each dict having only
               one key/value with the key being the value of the state that should trigger the event and the
               value a dict of parameter values that should be dispatched with the signal.
        """

        # set up logging
        loggername = "%s.%s" % (self.__class__.__module__,
                                self.__class__.__name__)
        self.logger = logging.getLogger(loggername)

        self.kiosk = kiosk
        self.events_table = events_table
        self.zero_condition = zero_condition
        self.event_state = event_state
        self.name = name
        self.comment = comment

        # assign evaluation function for states
        if self.zero_condition == 'falling':
            self._evaluate_state = self._zero_condition_falling
        elif self.zero_condition == 'rising':
            self._evaluate_state = self._zero_condition_rising
        elif self.zero_condition == 'either':
            self._evaluate_state = self._zero_condition_either

        # assign Nones to self.zero_condition_signs to signal
        # that the sign have not yet been evaluated
        self.previous_signs = [None]*len(self.events_table)

        # get signal from signals module
        if not hasattr(signals, event_signal):
            msg = "Signal '%s' not defined in pcse.signals module."
            raise exc.PCSEError(msg % event_signal)
        self.event_signal = getattr(signals, event_signal)

        # Build a counter for the state events.
        self.states_with_events = Counter()
        for ev in self.events_table:
            self.states_with_events.update(ev.keys())

        # Check if there are days with two or more events under the
        # same signal which is not allowed.
        multi_states = []
        for state, count in self.states_with_events.items():
            if count > 1:
                multi_states.append(state)
        if multi_states:
            msg = "Found states with more than 1 event for events table '%s' for state: %s"
            raise exc.PCSEError(msg % (self.name, multi_states))

    def __call__(self, day):
        """Runs the TimedEventDispatcher to determine if any actions are needed.

        :param day: a date object for the current simulation day
        :return: None
        """
        if not self.event_state in self.kiosk:
            msg = "State variable '%s' not (yet) available in kiosk!" % self.event_state
            self.logger.warning(msg)
            return

        # Determine if any event should be trigger based on the current state and
        # the event_condition.
        current_state = self.kiosk[self.event_state]
        zero_condition_signs = []
        for event, zero_condition_sign in zip(self.events_table, self.previous_signs):
            state, keywords = event.items()[0]
            zcs = self._evaluate_state(current_state, state, keywords, zero_condition_sign)
            zero_condition_signs.append(zcs)
        self.previous_signs = zero_condition_signs


    def _zero_condition_falling(self, current_state, state, keywords, zero_condition_sign):
        sign = cmp2(current_state - state, 0)

        # is None: e.g. called the first time and zero_condition_sign is not yet calculated
        if zero_condition_sign is None:
            return sign

        if zero_condition_sign == 1 and sign in [-1, 0]:
            msg = "State event dispatched from '%s' at event_state %s" % (self.name, state)
            self.logger.info(msg)
            self._send_signal(signal=self.event_signal, **keywords)

        return sign

    def _zero_condition_rising(self, current_state, state, kwargs, zero_condition_sign):
        sign = cmp2(current_state - state, 0)

        # is None: e.g. called the first time and zero_condition_sign is not yet calculated
        if zero_condition_sign is None:
            return sign

        if zero_condition_sign == -1 and sign in [0, 1]:
            msg = "State event dispatched from '%s' at model state %s" % (self.name, current_state)
            self.logger.info(msg)
            self._send_signal(signal=self.event_signal, **kwargs)

        return sign

    def _zero_condition_either(self, current_state, state, keywords, zero_condition_sign):
        sign = cmp2(current_state - state, 0)

        # is None: e.g. called the first time and zero_condition_sign is not yet calculated
        if zero_condition_sign is None:
            return sign

        if (zero_condition_sign == 1 and sign in [-1, 0]) or \
           (zero_condition_sign == -1 and sign in [0, 1]):
            msg = "State event dispatched from %s at event_state %s" % (self.name, state)
            self.logger.info(msg)
            self._send_signal(signal=self.event_signal, **keywords)

        return sign


class AgroManager(AncillaryObject):
    """Class for continuous AgroManagement actions including crop rotations and events.

    See also the documentation for the classes `CropCalendar`, `TimedEventDispatcher` and `StateEventDispatcher`.

    The AgroManager takes care of executing agromanagent actions that typically occur on agricultural
    fields including planting and harvesting of the crop, as well as management actions such as fertilizer
    application, irrigation, mowing and spraying.

    The agromanagement during the simulation is implemented as a sequence of campaigns. Campaigns start on a
    prescribed calendar date and finalize when the next campaign starts. The simulation ends either explicitly by
    provided a trailing empty campaign or by deriving the end date from the crop calendar and timed events in the
    last campaign. See also the section below on `end_date` property.

    Each campaign is characterized by zero or one crop calendar, zero or more timed events and zero or more
    state events.
    The structure of the data needed as input for AgroManager is most easily understood with the example
    (in YAML) below. The definition consists of three campaigns, the first starting on 1999-08-01, the second
    starting on 2000-09-01 and the last campaign starting on 2001-03-01. The first campaign consists of a crop
    calendar for winter-wheat starting with sowing at the given crop_start_date. During the campaign there are
    timed events for irrigation at 2000-05-25 and 2000-06-30. Moreover, there are state events for  fertilizer
    application (event_signal: apply_npk) given by development stage (DVS) at DVS 0.3, 0.6 and 1.12.

    The second campaign has no crop calendar, timed events or state events. This means that this is a period of
    bare soil with only the water balance running. The third campaign is for fodder maize sown at 2001-04-15
    with two series of timed events (one for irrigation and one for N/P/K application) and no state events.
    The end date of the simulation in this case will be 2001-11-01 (2001-04-15 + 200 days).

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

    """

    # campaign start dates
    campaign_start_dates = List()

    # Overall engine start date and end date
    _start_date = Instance(date)
    _end_date = Instance(date)

    # campaign definitions
    crop_calendars = List()
    timed_event_dispatchers = List()
    state_event_dispatchers = List()

    _tmp_date = None  # Helper variable
    _icampaign = 0  # count the campaigns

    def initialize(self, kiosk, agromanagement):
        """Initialize the AgroManager.

        :param kiosk: A PCSE variable Kiosk
        :param agromanagement: the agromanagement definition, see the example above in YAML.
        """


        self.kiosk = kiosk
        self.crop_calendars = []
        self.timed_event_dispatchers = []
        self.state_event_dispatchers = []
        self.campaign_start_dates = []

        # Connect CROP_FINISH signal with handler
        self._connect_signal(self._on_CROP_FINISH, signals.crop_finish)

        # First get and validate the dates of the different campaigns
        for campaign in agromanagement:
            # Check if campaign start dates is in chronological order
            campaign_start_date = campaign.keys()[0]
            self._check_campaign_date(campaign_start_date)
            self.campaign_start_dates.append(campaign_start_date)

        # Add None to the list of campaign dates to signal the end of the
        # number of campaigns.
        self.campaign_start_dates.append(None)

        # Walk through the different campaigns and build crop calendars and
        # timed/state event dispatchers
        for campaign, campaign_start, next_campaign in \
                zip(agromanagement, self.campaign_start_dates[:-1], self.campaign_start_dates[1:]):

            # Get the campaign definition for the start date
            campaign_def = campaign[campaign_start]

            if self._is_empty_campaign(campaign_def):  # no campaign definition for this campaign, e.g. fallow
                self.crop_calendars.append(None)
                self.timed_event_dispatchers.append(None)
                self.state_event_dispatchers.append(None)
                continue

            # get crop calendar definition for this campaign
            cc_def = campaign_def['CropCalendar']
            if cc_def is not None:
                cc = CropCalendar(kiosk, **cc_def)
                cc.validate(campaign_start, next_campaign)
                self.crop_calendars.append(cc)
            else:
                self.crop_calendars.append(None)

            # Get definition of timed events and build TimedEventsDispatchers
            te_def = campaign_def['TimedEvents']
            if te_def is not None:
                te_dsp = self._build_TimedEventDispatchers(kiosk, te_def)
                for te in te_dsp:
                    te.validate(campaign_start, next_campaign)
                self.timed_event_dispatchers.append(te_dsp)
            else:
                self.timed_event_dispatchers.append(None)

            # Get definition of state events and build StateEventsDispatchers
            se_def = campaign_def['StateEvents']
            if se_def is not None:
                se_dsp = self._build_StateEventDispatchers(kiosk, se_def)
                self.state_event_dispatchers.append(se_dsp)
            else:
                self.state_event_dispatchers.append(None)

    def _is_empty_campaign(self, campaign_def):
        """"Check if the campaign definition is empty"""

        if campaign_def is None:
            return True

        attrs = ["CropCalendar", "TimedEvents", "StateEvents"]
        r = []
        for attr in attrs:
            if attr in campaign_def:
                if campaign_def[attr] is None:
                    r.append(True)
                else:
                    r.append(False)
        if r == [True]*3:
            return True

        return False

    @property
    def start_date(self):
        """Retrieves the start date of the agromanagement sequence, e.g. the first simulation date

        :return: a date object
        """
        if self._start_date is None:
            self._start_date = self.campaign_start_dates[0]

        return self._start_date

    @property
    def end_date(self):
        """Retrieves the end date of the agromanagement sequence, e.g. the last simulation date.

        :return: a date object

        Getting the last simulation date is more complicated because there are two options.

        **1. Adding an explicit trailing empty campaign**

        The first option is to explicitly define the end date of the simulation by adding a
        'trailing empty campaign' to the agromanagement definition.
        An example of an agromanagement definition with a 'trailing empty campaigns' (YAML format) is
        given below. This example will run the simulation until 2001-01-01::

            Version: 1.0
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
                StateEvents:
            - 2001-01-01:

        Note that in configurations where the last campaign contains a definition for state events, a trailing
        empty campaign *must* be provided because the end date cannot be determined. The following campaign
        definition will therefore lead to an error::

            Version: 1.0
            AgroManagement:
            - 2001-01-01:
                CropCalendar:
                    crop_id: fodder-maize
                    crop_start_date: 2001-04-15
                    crop_start_type: sowing
                    crop_end_date:
                    crop_end_type: maturity
                    max_duration: 200
                TimedEvents:
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


        **2. Without an explicit trailing campaign**

        The second option is that there is no trailing empty campaign and in that case the end date of the simulation
        is retrieved from the crop calendar and/or the timed events that are scheduled. In the example below, the
        end date will be 2000-08-05 as this is the harvest date and there are no timed events scheduled after this
        date::

            Version: 1.0
            AgroManagement:
            - 1999-09-01:
                CropCalendar:
                    crop_id: winter-wheat
                    crop_start_date: 1999-10-01
                    crop_start_type: sowing
                    crop_end_date: 2000-08-05
                    crop_end_type: harvest
                    max_duration: 330
                TimedEvents:
                -   event_signal: irrigate
                    name:  Timed irrigation events
                    comment: All irrigation amounts in cm
                    events_table:
                    - 2000-05-01: {irrigation_amount: 2, efficiency: 0.7}
                    - 2000-06-21: {irrigation_amount: 5, efficiency: 0.7}
                    - 2000-07-18: {irrigation_amount: 3, efficiency: 0.7}
                StateEvents:

        In the case that there is no harvest date provided and the crop runs till maturity, the end date from
        the crop calendar will be estimated as the crop_start_date plus the max_duration.

        """
        if self._end_date is None:

            # First check if the last campaign definition is an empty trailing campaign and use that date.
            if self.crop_calendars[-1] is None and \
               self.timed_event_dispatchers[-1] is None and \
               self.state_event_dispatchers[-1] is None:
                self._end_date = self.campaign_start_dates[-2]  # use -2 here because None is
                                                                # appended to campaign_start_dates
                return self._end_date

            # Check if there are state events defined in the last campaign without specifying the end date
            # explicitly with an trailing empty campaign
            if self.state_event_dispatchers[-1] is not None:
                msg = "In the AgroManagement definition, the last campaign with start date '%s' contains StateEvents. " \
                      "When specifying StateEvents, the end date of the campaign must be explicitly" \
                      "given by a trailing empty campaign."
                raise exc.PCSEError(msg)

            # Walk over the crop calendars and timed events to get the last date.
            cc_dates = []
            te_dates = []
            for cc, teds in zip(self.crop_calendars, self.timed_event_dispatchers):
                if cc is not None:
                    cc_dates.append(cc.get_end_date())
                if teds is not None:
                    te_dates.extend([t.get_end_date() for t in teds])

            # If no end dates can be found raise an error because the agromanagement sequence
            # consists only of empty campaigns
            if not cc_dates and not te_dates:
                msg = "Empty agromanagement definition: no campaigns with crop calendars or timed events provided!"
                raise exc.PCSEError(msg)

            end_date = date(1, 1, 1)
            if cc_dates:
                end_date = max(max(cc_dates), end_date)
            if te_dates:
                end_date = max(max(te_dates), end_date)
            self._end_date = end_date

        return self._end_date

    def _check_campaign_date(self, campaign_start_date):
        """
        :param campaign_start_date: Start date of the agricultural campaign
        :return: None
        """
        if not isinstance(campaign_start_date, date):
            msg = "Campaign start must be given as a date."
            raise exc.PCSEError(msg)

        if self._tmp_date is None:
            self._tmp_date = campaign_start_date
        else:
            if campaign_start_date <= self._tmp_date:
                msg = "The agricultural campaigns are not sequential " \
                      "in the agromanagement definition."
                raise exc.PCSEError(msg)

    def _build_TimedEventDispatchers(self, kiosk, event_definitions):
        r = []
        for ev_def in event_definitions:
            ev_dispatcher = TimedEventsDispatcher(kiosk, **ev_def)
            r.append(ev_dispatcher)
        return r

    def _build_StateEventDispatchers(self, kiosk, event_definitions):
        r = []
        for ev_def in event_definitions:
            ev_dispatcher = StateEventsDispatcher(kiosk, **ev_def)
            r.append(ev_dispatcher)
        return r

    def __call__(self, day, drv):
        """Calls the AgroManager to execute and crop calendar actions, timed or state events.

        :param day: The current simulation date
        :param drv: The driving variables for the current day
        :return: None
        """

        # Check if the agromanager should switch to a new campaign
        if day == self.campaign_start_dates[self._icampaign+1]:
            self._icampaign += 1
            # if new campaign, throw out the previous campaign definition
            self.crop_calendars.pop(0)
            self.timed_event_dispatchers.pop(0)
            self.state_event_dispatchers.pop(0)

        # call handlers for the crop calendar, timed and state events
        if self.crop_calendars[0] is not None:
            self.crop_calendars[0](day)

        if self.timed_event_dispatchers[0] is not None:
            for ev_dsp in self.timed_event_dispatchers[0]:
                ev_dsp(day)

        if self.state_event_dispatchers[0] is not None:
            for ev_dsp in self.state_event_dispatchers[0]:
                ev_dsp(day)

    def _on_CROP_FINISH(self, day):
        """Send signal to terminate after the crop cycle finishes.

        The simulation will be terminated when the following conditions are met:
        1. There are no campaigns defined after the current campaign
        2. There are no StateEvents active
        3. There are no TimedEvents scheduled after the current date.
        """


        if self.campaign_start_dates[self._icampaign+1] is not None:
            return  #  e.g. There is a next campaign defined

        if self.state_event_dispatchers[0] is not None:
            return  # there are state events active that may trigger in the future

        if self.timed_event_dispatchers[0] is not None:
            a = 1
            end_dates = [t.get_end_date() for t in self.timed_event_dispatchers[0]]
            if end_dates:
                if max(end_dates) > day:  # There is at least one scheduled event after the current day
                    return
        self._send_signal(signal=signals.terminate)
