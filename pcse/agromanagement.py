# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Classes that simulate AgroManagement actions in PCSE.

Agromanagement actions are those actions that are originating from
human actors rather than biophysical principles. Examples are sowing
and harvesting, etc.
"""

import datetime

from .base_classes import ParameterProvider, AncillaryObject, ParamTemplate
from .traitlets import HasTraits, Float, Int, Instance, Enum, Bool
from . import signals
from . import exceptions as exc
from .util import ConfigurationLoader


class AgroManagementSingleCrop(AncillaryObject):
    """Agromanagement for running a single crop season.
    
    Simple agromanagement class which takes care of:
    
    * signalling the start of the crop cycle and initializing the crop
      simulation object
    * signalling the end of the crop cycle (if the crop is to stop
      at a certain harvest date)
    * signalling the maximum crop cycle duration.

    **Simulation parameters:**
    
    ================ =========================================== =======  ======
     Name            Description                                  Type     Unit
    ================ =========================================== =======  ======
    CROP_START_DATE  Date object describing start date of           STi     -
                     the crop simulation.
    CROP_START_TYPE  string representing start type: 'sowing'|      STi     -
                     'emergence'
    CROP_END_DATE    date representing end  date of the crop        STi     -
                     simulation. This date will not be used
                     in case CROP_END_TYPE=='maturity'.
    CROP_END_TYPE    String representing end type: 'harvest'|       STi     -
                     'maturity'|'earliest'
    MAX_DURATION     Integer describing maximum duration of         STi    days
                     the crop simulation
    ================ =========================================== =======  ======
    
    **Signals sent or handled:**

    * "CROP_START": sent when `day == CROP_START_DATE`
    * "CROP_FINISH": sent when `day == CROP_END_DATE`
    * "TERMINATE": sent when a "CROP_FINISH" signal is received.
    """
    # system configuration
    mconf = Instance(ConfigurationLoader)

    # Placeholders for the parameters that are needed to start the crop
    parameterprovider = Instance(ParameterProvider)
    timerdata = Instance(dict)
    soildata = Instance(dict)
    cropdata = Instance(dict)
    sitedata = Instance(dict)
    duration = Int(0)
    in_crop_cycle = Bool(False)
    
    class Parameters(ParamTemplate):
        MAX_DURATION = Int(-99)
        CROP_START_DATE = Instance(datetime.date)
        CROP_START_TYPE = Enum(["sowing", "emergence"])
        CROP_END_DATE = Instance(datetime.date)
        CROP_END_TYPE = Enum(["maturity", "harvest", "earliest"])
                 
    def initialize(self, day, kiosk, mconf, parameterprovider):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE model instance
        :param mconf: ConfigurationLoader instance
        :param parameterprovider: `ParameterProvider` object providing parameters as
                key/value pairs
        """
        self.params = self.Parameters(parameterprovider)
        self.kiosk = kiosk
        self.mconf = mconf
        self.duration = 0
        self.parameterprovider = parameterprovider

        # TERMINATE signal should be issued directly after signal CROP_FINISH
        # This handler takes care of that
        self._connect_signal(self._on_CROP_FINISH, signal=signals.crop_finish)
        
        # Check if sequence of dates is OK, otherwise the crop
        # will never start or finish
        if day > self.params.CROP_START_DATE:
            msg = ("CROP_START_DATE before simulation start day: "
                   "crop simulation will never start.")
            raise exc.PCSEError(msg)
        if self.params.CROP_END_TYPE in ("harvest","earliest"):
            if self.params.CROP_END_DATE <= self.params.CROP_START_DATE:
                msg = ("CROP_END_DATE <= CROP_START_DATE: " +
                       "crop simulation will never finish!")
                raise exc.PCSEError(msg)
                

    def __call__(self, day, drv):
        
        self.duration += 1
        
        # Check if crop sowing/emergence date is reached.
        if day == self.params.CROP_START_DATE:
            if self.in_crop_cycle:
                msg = ("Crop sowing/emergence date reached while existing " +
                       "crop still active!")
                raise exc.PCSEError(msg)

            # Initialize the crop simulation object and send it to the
            # combined soil/crop system using the CROP_START signal.
            cropsimulation = self.mconf.CROP(day, self.kiosk, self.parameterprovider)
            self._start_new_crop(day, cropsimulation)
        
        finish_cropsimulation = False
        # Check if CROP_END_DATE is reached for CROP_END_TYPE harvest/earliest
        if self.params.CROP_END_TYPE in ["harvest", "earliest"]:
            if day >= self.params.CROP_END_DATE:
                finish_cropsimulation = True
                finish_type = "harvest"

        # Check for forced stop because maximum duration is reached
        if self.in_crop_cycle and self.duration >= self.params.MAX_DURATION:
            finish_cropsimulation = True
            finish_type = "max_duration"
        
        # If finish condition is reached send a signal to finish the crop
        if finish_cropsimulation is True:
            self.in_crop_cycle = False
            self._send_signal(signal=signals.crop_finish, day=day,
                              finish_type=finish_type)
    
    def _start_new_crop(self, day, cropsimulation):
        """Starts a new simulation by sending the apropriate signal and
        variables.
        """
        self.duration = 0
        self.in_crop_cycle = True
        self._send_signal(signal=signals.crop_start, day=day,
                          cropsimulation=cropsimulation)
                
    def _on_CROP_FINISH(self):
        "Send signal to terminate system after real crop is finished."

        self.in_crop_cycle = False
        self._send_signal(signal=signals.terminate)
