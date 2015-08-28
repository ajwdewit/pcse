# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""This module defines and describes the signals used by PCSE

Signals are used by PCSE to notify components of events such as sowing,
harvest and termination. Events can be send by any SimulationObject through
its `SimulationObject._send_signal()` method. Similarly, any SimulationObject
can receive signals by registering a handler through the
`SimulationObject._connect_signal()` method.
Variables can be passed to the handler of the signal through
positional or keyword arguments. However, it is highly discouraged to use
positional arguments when sending signals in order to avoid conflicts between
positional and keyword arguments.

An example can help to clarify how signals are used in PCSE but check also the
documentation of the PyDispatcher_ package for more information::

    import sys, os
    import math
    sys.path.append('/home/wit015/Sources/python/pcse/')
    import datetime as dt
    
    import pcse
    from pcse.base_classes import SimulationObject, VariableKiosk
    
    mysignal = "My first signal"
    
    class MySimObj(SimulationObject):
        
        def initialize(self, day, kiosk):
            self._connect_signal(self.handle_mysignal, mysignal)
    
        def handle_mysignal(self, arg1, arg2):
            print "Value of arg1,2: %s, %s" % (arg1, arg2)
    
        def send_signal_with_exact_arguments(self):
            self._send_signal(signal=mysignal, arg2=math.pi, arg1=None)
    
        def send_signal_with_more_arguments(self):
            self._send_signal(signal=mysignal, arg2=math.pi, arg1=None, 
                              extra_arg="extra")
    
        def send_signal_with_missing_arguments(self):
            self._send_signal(signal=mysignal, arg2=math.pi, extra_arg="extra")
    
            
    # Create an instance of MySimObj
    day = dt.date(2000,1,1)
    k = VariableKiosk()
    mysimobj = MySimObj(day, k)
    
    # This sends exactly the right amount of keyword arguments
    mysimobj.send_signal_with_exact_arguments()
    
    # this sends an additional keyword argument 'extra_arg' which is ignored.
    mysimobj.send_signal_with_more_arguments()
    
    # this sends the signal with a missing 'arg1' keyword argument which the handler
    # expects and thus causes an error, raising a TypeError
    try:
        mysimobj.send_signal_with_missing_arguments()
    except TypeError, exc:
        print "TypeError occurred: %s" % exc

Saving this code as a file `test_signals.py` and importing it gives the
following output::

    >>> import test_signals
    Value of arg1,2: None, 3.14159265359
    Value of arg1,2: None, 3.14159265359
    TypeError occurred: handle_mysignal() takes exactly 3 non-keyword arguments (1 given)

Currently the following signals are used within PCSE with the following
keywords.*

**CROP_START**

 Indicates that a new crop cycle will start.
 
 self._send_signal(signal=signals.crop_start, day=<date>,
                   crop_id=<string>, crop_start_type=<string>,
                   crop_end_type=<string>)

 keyword arguments with signals.crop_start:
    
    * day: Current date
    * crop_id: a string identifying the crop

**CROP_FINISH**

 Indicates that the current crop cycle is finished
 
 self._send_signal(signal=signals.crop_finish, day=<date>,
                   finish_type=<string>, crop_delete=<True|False>)

keyword arguments with signals.crop_finish:
    
    * day: Current date
    * finish_type: string describing the reason for finishing the simulation, e.g.
      maturity, harvest, all leaves died, maximum duration reached, etc.
    * crop_delete: Set to True when the CropSimulation object must be deleted
      from the system, for example for the implementation of crop rotations.
      Defaults to False.

**TERMINATE**
 
 Indicates that the entire system should terminate (crop & soil water balance) and
 that terminal output should be collected.

 self._send_signal(signal=signals.terminate)

 No keyword arguments are defined for this signal

**OUTPUT**

 Indicates that the model state should be saved for later use.

 self._send_signal(signal=signals.output)
 
 No keyword arguments are defined for this signal

**SUMMARY_OUTPUT**

 Indicates that the model state should be saved for later use,
 SUMMARY_OUTPUT is only generated when a CROP_FINISH signal is
 received indicating that the crop simulation must finish.

 self._send_signal(signal=signals.output)

 No keyword arguments are defined for this signal

**APPLY_NPK**

Is used for application of Nitrate/Phosphate/Potassium (N/P/K) fertilizer

self._send_signal(signal=signals.apply_npk, N_amount=<float>, P_amount=<float>, K_amount=<float>,
                  N_recovery<float>, P_recovery=<float>, K_recovery=<float>)

Keyword arguments with signals.apply_npk:

    * N/P/K_amount: Amount of fertilizer in kg/ha applied on this day.
    * N/P/K_recovery: Recovery fraction for the given type of fertilizer


**IRRIGATE**

Is used for sending irrigation events.

self._send_signal(signal=signals.irrigate, amount=<float>, efficiency=<float>)

Keyword arguments with signals.irrigate:

    * amount: Amount of irrigation in cm water applied on this day.
    * efficiency: efficiency of irrigation, meaning that the total amount of water that
      is added to the soil reservoir equals amount * efficiency


.. _PyDispatcher: http://pydispatcher.sourceforge.net/
"""

crop_start = "CROP_START"
crop_emerged = "CROP_EMERGED"
crop_finish = "CROP_FINISH"
terminate = "TERMINATE"
output = "OUTPUT"
summary_output = "SUMMARY_OUTPUT"
apply_npk = "APPLY_NPK"
apply_n = "APPLY_N"
irrigate = "IRRIGATE" 
