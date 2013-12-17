import datetime

from .pydispatch import dispatcher
from .base_classes import AncillaryObject
from .traitlets import HasTraits, Instance, Bool, Int, Enum
from . import signals
#from base_classes import VariableKiosk, SimulationObject, StatesTemplate
#from traitlets import Instance, Bool, Int, Enum
#import signals


class Timer(AncillaryObject):
    """This class implements a basic timer for use with the WOFOST crop model.
    
    This object implements a simple timer that increments the current time with
    a fixed time-step of one day at each call and returns its value. Moreover,
    it generates OUTPUT signals in daily, dekadal or monthly time-steps that
    can be caught in order to store the state of the simulation for later use.
        
    Initializing the timer::

        timer = Timer(start_date, kiosk, final_date, interval_type,
                      interval_days)
        CurrentDate = timer()
        
    **Signals sent or handled:**
 
        * "OUTPUT": sent when the condition for generating output is True
          which depends on the output type and interval.
 

  """

    start_date   = Instance(datetime.date)
    final_date   = Instance(datetime.date)
    current_date = Instance(datetime.date)
    time_step    = Instance(datetime.timedelta)
    interval_type = Enum(["daily", "dekadal", "monthly"])
    interval_days = Int
    day_counter   = Int
    first_call = Bool

    def initialize(self, start_date, kiosk, final_date, interval_type="daily",
                    interval_days=1):
        """
        :param day: Start date of the simulation
        :param kiosk: Variable kiosk of the PyWOFOST instance
        :param final_date: Final date of the simulation. For example, this date
            represents (START_DATE + MAX_DURATION) for a single cropping season.
            This date is *not* the harvest date because signalling harvest is taken
            care of by the `AgroManagement` module.
        :param interval_type: Interval type for storing simulation results through
            OUTPUT signals. Can be one of "daily"|"dekadal"|"monthly" defaults
            to "daily".
        :param interval_days: Number of days between daily output, defaults to 1.
            Is ignored in case of dekadal or monthly output.
        """
        
        self.kiosk = kiosk
        self.start_date = start_date
        self.final_date = final_date
        self.current_date = start_date
        self.day_counter = 0
        self.interval_type = interval_type.lower()
        self.interval_days = interval_days
        self.time_step = datetime.timedelta(days=1)
        self.first_call = True

    def _is_a_month(self, day):
        """Returns True if the date is on the last day of a month."""

        if day.month==12:
            if (day == datetime.date(day.year, day.month, 31)):
                return True
        else:
            if (day == datetime.date(day.year, day.month+1, 1) - \
                       datetime.timedelta(days=1)):
                return True
        
        return False

    def _is_a_dekad(self, day):
        """Returns True if the date is on a dekad boundary, i.e. the 10th,
        the 20th or the last day of each month"""
        if day.month==12:
            if (day == datetime.date(day.year, day.month, 10)):
                return True
            elif (day == datetime.date(day.year, day.month, 20)):
                return True
            elif (day == datetime.date(day.year, day.month, 31)):
                return True
        else:
            if (day == datetime.date(day.year, day.month, 10)):
                return True
            elif (day == datetime.date(day.year, day.month, 20)):
                return True
            elif (day == datetime.date(day.year, day.month+1, 1) - \
                       datetime.timedelta(days=1)):
                return True
        
        return False
        
    def __call__(self):
        
        # On first call only return the current date, do not increase time
        if self.first_call is True:
            self.first_call = False
            self.logger.info("Model time at first call: %s" % self.current_date)
        else:
            self.current_date += self.time_step
            self.day_counter += 1
            self.logger.info("Model time updated to: %s" % self.current_date)

        # Check if output should be generated
        output = False
        if self.interval_type == "daily":
            if (self.day_counter % self.interval_days) == 0:
                output = True
        elif self.interval_type == "dekadal":
            if self._is_a_dekad(self.current_date):
                output = True
        elif self.interval_type == "monthly":
            if self._is_a_month(self.current_date):
                output = True

        # Send output signal if True
        if output:
            self._send_signal(signal=signals.output)
            
        # If final date is reached send the terminate signal
        if (self.current_date >= self.final_date):
            self._send_signal(signal=signals.terminate)
            
        return self.current_date

def simple_test():
    "Only used for testing timer routine"

    def on_OUTPUT():
        print "Output generated."
    
    Start = datetime.date(2000,1,1)
    End = datetime.date(2000,2,1)
    dispatcher.connect(on_OUTPUT, signal=signals.output,
                       sender=dispatcher.Any)
    timer = Timer(Start, End, "Dekadal")
    print "-----------------------------------------"
    print "Dekadal output"
    print "-----------------------------------------"
    for i in range(100):
        today = timer()
    print "-----------------------------------------"
    print "Monthly output"
    print "-----------------------------------------"
    timer = Timer(Start, End, "Monthly")
    for i in range(150):
        today = timer()
    print "-----------------------------------------"
    print "daily output with 4 day intervals"
    print "-----------------------------------------"
    timer = Timer(Start, End, "daily", 4)
    for i in range(150):
        today = timer()


if __name__ == '__main__':
    simple_test()
