# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import datetime

from .base_classes import AncillaryObject, ParamTemplate
from .traitlets import HasTraits, Float, Int, Instance, Enum, Bool
from . import signals
from . import exceptions as exc

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
    START_DATE       Date object describing start date of           STi     -
                     the entire simulation period.
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
    * "CROP_FINISH": sent when `day == CROP_END_DATE` and CROP_END_TYPE in ['harvest','earliest']
    * "TERMINATE": sent when a "CROP_FINISH" signal is received.
    """

    # Placeholders for the parameters that are needed to start the crop
    duration = Int(0)
    in_crop_cycle = Bool(False)

    class Parameters(ParamTemplate):
        MAX_DURATION = Int(-99)
        START_DATE = Instance(datetime.date)
        CROP_START_DATE = Instance(datetime.date)
        CROP_START_TYPE = Enum(["sowing", "emergence"])
        CROP_END_DATE = Instance(datetime.date)
        CROP_END_TYPE = Enum(["maturity", "harvest", "earliest"])

    def initialize(self, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: a ParameterProvider object providing parameters as key/value pairs
        """
        self.params = self.Parameters(parvalues)
        self.kiosk = kiosk
        self.duration = 0

        # TERMINATE signal should be issued directly after signal CROP_FINISH
        # This handler takes care of that
        self._connect_signal(self._on_CROP_FINISH, signal=signals.crop_finish)

        # Check if sequence of dates is OK, otherwise the crop
        # will never start or finish
        if self.params.START_DATE > self.params.CROP_START_DATE:
            msg = ("CROP_START_DATE before simulation start day: "
                   "crop simulation will never start.")
            raise exc.PCSEError(msg)
        if self.params.CROP_END_TYPE in ("harvest", "earliest"):
            if self.params.CROP_END_DATE <= self.params.CROP_START_DATE:
                msg = ("CROP_END_DATE <= CROP_START_DATE: "
                       "crop simulation will never finish!")
                raise exc.PCSEError(msg)

    def __call__(self, day, drv):

        self.duration += 1

        # Check if crop sowing/emergence date is reached.
        if day == self.params.CROP_START_DATE:
            if self.in_crop_cycle:
                msg = ("Crop sowing/emergence date reached while existing "
                       "crop still active!")
                raise exc.PCSEError(msg)
            self.duration = 0
            self.in_crop_cycle = True
            self._send_signal(signal=signals.crop_start, day=day,
                              crop_start_type=self.params.CROP_START_TYPE,
                              crop_end_type=self.params.CROP_END_TYPE)

        # Check if CROP_END_DATE is reached for CROP_END_TYPE harvest/earliest
        finish_type = None
        if self.params.CROP_END_TYPE in ["harvest", "earliest"]:
            if day >= self.params.CROP_END_DATE:
                finish_type = "harvest"

        # Check for forced stop because maximum duration is reached
        if self.in_crop_cycle and self.duration >= self.params.MAX_DURATION:
            finish_type = "max_duration"

        # If finish condition is reached send a signal to finish the crop
        if finish_type is not None:
            self.in_crop_cycle = False
            self._send_signal(signal=signals.crop_finish, day=day,
                              finish=finish_type)

    def _on_CROP_FINISH(self):
        """Send signal to terminate system after real crop is finished.
        """
        self.in_crop_cycle = False
        self._send_signal(signal=signals.terminate)
