from __future__ import print_function
from datetime import date, timedelta

import yaml
import logging
from collections import Counter
import numpy as np

from pcse.base_classes import DispatcherObject, VariableKiosk, SimulationObject
from pcse.traitlets import HasTraits, Float, Int, Instance, Enum, Bool, List, Dict, Unicode
from pcse import exceptions as exc

class TestSignals(object):
    APPLY_NPK = "APPLY_NPK"
    IRRIGATE = "IRRIGATE"
signals = TestSignals()

class MyModel(SimulationObject):
    def initialize(self, day, kiosk, *args):

        self._connect_signal(self._on_SIGNAL, signals.APPLY_NPK)
        self._connect_signal(self._on_SIGNAL, signals.IRRIGATE)

    def _on_SIGNAL(self, signal, sender, **kwargs):
        msg = "signal %s received with args: %s" % (signal, kwargs)
        print(msg)


class TimedEventsDispatcher(HasTraits, DispatcherObject):
    event_signal = None
    events_table = List()
    days_with_events = Instance(Counter)
    kiosk = Instance(VariableKiosk)
    logger = Instance(logging.Logger)
    name = Unicode()
    comment = Unicode()

    def __init__(self, kiosk, event_signal, name, comment, events_table):

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

    def __call__(self, day):
        if not day in self.days_with_events:
            return

        for event in self.events_table:
            if day in event:
                msg = "Time event dispatched from '%s' at day %s" % (self.name, day)
                print(msg)
                kwargs = event[day]
                self._send_signal(signal=self.event_signal, **kwargs)


class StateEventsDispatcher(HasTraits, DispatcherObject):
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

    def __call__(self):
        if not self.event_state in self.kiosk:
            msg = "State variable '%s' not available in kiosk!" % self.event_state
            raise exc.PCSEError(msg)

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
        sign = cmp(current_state - state, 0)

        # is None: e.g. called the first time and zero_condition_sign is not yet calculated
        if zero_condition_sign is None:
            return sign

        if zero_condition_sign == 1 and sign in [-1, 0]:
            msg = "State event dispatched from '%s' at event_state %s" % (self.name, state)
            print(msg)
            self._send_signal(signal=self.event_signal, **keywords)

        return sign

    def _zero_condition_rising(self, current_state, state, kwargs, zero_condition_sign):
        sign = cmp(current_state - state, 0)

        # is None: e.g. called the first time and zero_condition_sign is not yet calculated
        if zero_condition_sign is None:
            return sign

        if zero_condition_sign == -1 and sign in [0, 1]:
            msg = "State event dispatched from '%s' at model state %s" % (self.name, current_state)
            print(msg)
            self._send_signal(signal=self.event_signal, **kwargs)

        return sign

    def _zero_condition_either(self, current_state, state, keywords, zero_condition_sign):
        sign = cmp(current_state - state, 0)

        # is None: e.g. called the first time and zero_condition_sign is not yet calculated
        if zero_condition_sign is None:
            return sign

        if (zero_condition_sign == 1 and sign in [-1, 0]) or \
           (zero_condition_sign == -1 and sign in [0, 1]):
            msg = "State event dispatched from %s at event_state %s" % (self.name, state)
            print(msg)
            self._send_signal(signal=self.event_signal, **keywords)

        return sign


def main():
    with open("events.yaml") as fp:
        event_definitions = yaml.load(fp)

    startdate = date(2000, 1, 1)
    kiosk = VariableKiosk()
    kiosk.register_variable(1, "DVS", type="S", publish=True)
    kiosk.register_variable(1, "SM", type="S", publish=True)

    # register a model that listens to signals sent by events.
    model = MyModel(startdate, kiosk)

    # Build the dispatchers for time and state events
    timed_event_dispatchers = []
    for event_def in event_definitions['TimedEvents']:
        timed_event_dispatchers.append(TimedEventsDispatcher(kiosk, **event_def))
    state_event_dispatchers = []
    for event_def in event_definitions['StateEvents']:
        state_event_dispatchers.append(StateEventsDispatcher(kiosk, **event_def))

    # The 'simulation loop'
    doys = range(100)
    development = np.arange(0,2,0.02)
    np.random.seed(1000)
    soil_moisture = np.random.rand(200)
    for doy, dvs, sm in zip(doys, development, soil_moisture):

        day = startdate + timedelta(days=doy)
        kiosk.set_variable(1, "DVS", dvs)
        kiosk.set_variable(1, "SM", sm)
        print("day: %s, DVS: %5.2f, SM: %5.3f" % (day, dvs, sm))

        # Fire timed events
        for ed in timed_event_dispatchers:
            ed(day)
        # Fire state events
        for ed in state_event_dispatchers:
            ed()

if __name__ == "__main__":
    main()