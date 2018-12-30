# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
from __future__ import print_function
import datetime

from .pydispatch import dispatcher
from .base import AncillaryObject, VariableKiosk
from .traitlets import HasTraits, Instance, Bool, Int, Enum
from . import signals
from .util import is_a_dekad, is_a_month, is_a_week


class Timer(AncillaryObject):
    """This class implements a basic timer for use with the WOFOST crop model.
    
    This object implements a simple timer that increments the current time with
    a fixed time-step of one day at each call and returns its value. Moreover,
    it generates OUTPUT signals in daily, dekadal or monthly time-steps that
    can be caught in order to store the state of the simulation for later use.
        
    Initializing the timer::

        timer = Timer(start_date, kiosk, final_date, mconf)
        CurrentDate = timer()
        
    **Signals sent or handled:**
 
        * "OUTPUT": sent when the condition for generating output is True
          which depends on the output type and interval.
 

  """

    start_date = Instance(datetime.date)
    end_date = Instance(datetime.date)
    current_date = Instance(datetime.date)
    time_step = Instance(datetime.timedelta)
    interval_type = Enum(["daily", "weekly", "dekadal", "monthly"])
    output_weekday = Int()
    interval_days = Int()
    generate_output = Bool(False)
    day_counter = Int(0)
    first_call = Bool(True)
    _in_crop_cycle = Bool()

    def initialize(self, kiosk, start_date, end_date, mconf):
        """
        :param day: Start date of the simulation
        :param kiosk: Variable kiosk of the PCSE instance
        :param end_date: Final date of the simulation. For example, this date
            represents (START_DATE + MAX_DURATION) for a single cropping season.
            This date is *not* the harvest date because signalling harvest is taken
            care of by the `AgroManagement` module.
        :param mconf: A ConfigurationLoader object, the timer needs access to the
            configuration attributes mconf.OUTPUT_INTERVAL, mconf.OUTPUT_VARS and
            mconf.OUTPUT_INTERVAL_DAYS

        """
        
        self.kiosk = kiosk
        self.start_date = start_date
        self.end_date = end_date
        self.current_date = start_date
        # self.day_counter = 0
        # Settings for generating output. Note that if no OUTPUT_VARS are listed
        # in that case no OUTPUT signals will be generated.
        self.generate_output = bool(mconf.OUTPUT_VARS)
        self.interval_type = mconf.OUTPUT_INTERVAL.lower()
        self.output_weekday = mconf.OUTPUT_WEEKDAY
        self.interval_days = mconf.OUTPUT_INTERVAL_DAYS
        self.time_step = datetime.timedelta(days=1)
        # self.first_call = True

    def __call__(self):
        
        # On first call only return the current date, do not increase time
        if self.first_call is True:
            self.first_call = False
            self.logger.debug("Model time at first call: %s" % self.current_date)
        else:
            self.current_date += self.time_step
            self.day_counter += 1
            self.logger.debug("Model time updated to: %s" % self.current_date)

        # Check if output should be generated
        output = False
        if self.generate_output:
            if self.interval_type == "daily":
                if (self.day_counter % self.interval_days) == 0:
                    output = True
            elif self.interval_type == "weekly":
                if is_a_week(self.current_date, self.output_weekday):
                    output = True 
            elif self.interval_type == "dekadal":
                if is_a_dekad(self.current_date):
                    output = True
            elif self.interval_type == "monthly":
                if is_a_month(self.current_date):
                    output = True

        # Send output signal if True
        if output:
            self._send_signal(signal=signals.output)
            
        # If end date is reached send the terminate signal
        if self.current_date >= self.end_date:
            msg = "Reached end of simulation period as specified by END_DATE."
            self.logger.info(msg)
            self._send_signal(signal=signals.terminate)
            
        return self.current_date, float(self.time_step.days)


def simple_test():
    "Only used for testing timer routine"

    class Container(object):
        pass

    def on_OUTPUT():
        print("Output generated.")
    
    Start = datetime.date(2000, 1, 1)
    End = datetime.date(2000, 2, 1)
    kiosk = VariableKiosk()
    dispatcher.connect(on_OUTPUT, signal=signals.output,
                       sender=dispatcher.Any)

    mconf = Container()
    mconf.OUTPUT_INTERVAL = "dekadal"
    mconf.OUTPUT_INTERVAL_DAYS = 4
    mconf.OUTPUT_VARS = ["dummy"]

    print("-----------------------------------------")
    print("Dekadal output")
    print("-----------------------------------------")
    timer = Timer(Start, kiosk, End, mconf)
    for i in range(100):
        today = timer()

    print("-----------------------------------------")
    print("Monthly output")
    print("-----------------------------------------")
    mconf.OUTPUT_INTERVAL = "monthly"
    timer = Timer(Start, kiosk, End, mconf)
    for i in range(150):
        today = timer()

    print("-----------------------------------------")
    print("daily output with 4 day intervals")
    print("-----------------------------------------")
    mconf.OUTPUT_INTERVAL = "daily"
    timer = Timer(Start, kiosk, End, mconf)
    for i in range(150):
        today = timer()

if __name__ == '__main__':
    simple_test()
