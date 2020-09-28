# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import datetime

from ..traitlets import Float, Instance
from ..decorators import prepare_rates, prepare_states
from ..base import StatesTemplate, SimulationObject
from .. import signals

from .evapotranspiration import Simple_Evapotranspiration as Evapotranspiration
from .root_dynamics import Simple_Root_Dynamics as Root_Dynamics
from .leaf_dynamics import CSDM_Leaf_Dynamics as Leaf_Dynamics


class GreenLayerCrop(SimulationObject):
    """Top level object organizing the different components of the crop
    simulation for a crop model which only simulates water use but does not
    simulate growth of biomass, etc. The approach used here is very similar to
    the FAO Water Requirement Satisfaction Index (WRSI).

    The processes that are implemented as embedded simulation objects consist of:
    
        1. Evapotranspiration taken from the WOFOST model
        2. Leaf dynamics as defined by the CSDM model (a logistic/exponential LAI curve)
        3. Root dynamics taken from the WOFOST model

    **Simulation parameters:**
    
    None in this class, but see classes for evapotranspiration, leaf dynamics and
    root dynamics.
    
    **State variables:**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    CTRAT    Total crop transpiration                           N    cm
    DOF      Date representing the day of finish of the crop    N    -
             simulation.
    SumAET   Sum of actual crop + soil evapotranspiration       N    cm
    SumPET   Sum of potential crop + soil evapotranspiration    N    cm
    FINISH   String representing the reason for finishing the   N    -
             simulation: maturity, harvest, leave death, etc.
    WRSI     Water Requirement Satisfaction Index computed
             as SumAET/SumPET * 100                             N     %
    =======  ================================================= ==== ============

 
     **Rate variables:**

    None
    """
    
    # sub-model components for crop simulation
    evtra = Instance(SimulationObject)
    lv_dynamics = Instance(SimulationObject)
    ro_dynamics = Instance(SimulationObject)
    
    class StateVariables(StatesTemplate):
        CTRAT = Float(-99.) # Crop total transpiration
        DOF = Instance(datetime.date)
        FINISH = Instance(str)
        WRSI = Float()
        SumPET = Float() # Sum of potential crop evapotranspiration
        SumAET = Float() # Sum of actual crop evapotranspiration

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """

        self.kiosk = kiosk
        
        # Initialize components of the crop
        self.evtra = Evapotranspiration(day, kiosk, parvalues)
        self.ro_dynamics = Root_Dynamics(day, kiosk, parvalues)
        self.lv_dynamics = Leaf_Dynamics(day, kiosk, parvalues)

        self.states = self.StateVariables(kiosk, CTRAT=0.0, DOF=None, FINISH=None, WRSI=100,
                                          SumPET=0., SumAET=0.)
            
        # assign handler for CROP_FINISH signal
        self._connect_signal(self._on_CROP_FINISH, signal=signals.crop_finish)

    @prepare_rates
    def calc_rates(self, day, drv):
        states = self.states

        # (evapo)transpiration rates
        self.evtra(day, drv)

        # Root growth
        self.ro_dynamics.calc_rates(day, drv)
        # leaf growth
        self.lv_dynamics.calc_rates(day, drv)

    @prepare_states
    def integrate(self, day, delt=1.0):
        states = self.states
        
        # Integrate states on leaves, storage organs, stems and roots
        self.ro_dynamics.integrate(day, delt)
        self.lv_dynamics.integrate(day, delt)

        # total crop transpiration (CTRAT)
        states.CTRAT += self.kiosk["TRA"] * delt

        # Sum of total and actual evapotranspiration
        states.SumPET += (self.kiosk["TRAMX"] + self.kiosk["EVS"]) * delt
        states.SumAET += (self.kiosk["TRA"] + self.kiosk["EVS"]) * delt
        # compute Water Requirements Satisfaction Index according to FAO
        states.WRSI = states.SumAET/states.SumPET * 100
        
    def _on_CROP_FINISH(self, day, finish_type, *args, **kwargs):
        """Handler for setting day of finish (DOF) and reason for
        crop finishing (FINISH).
        """
        self._for_finalize["DOF"] = day
        self._for_finalize["FINISH"]= finish_type
