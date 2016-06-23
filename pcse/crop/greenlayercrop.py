# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import datetime

from ..traitlets import Float, Int, Instance, Enum, Unicode
from ..decorators import prepare_rates, prepare_states
from ..util import limit
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
     SimulationObject
from .. import signals
from .. import exceptions as exc


from .evapotranspiration import Simple_Evapotranspiration as Evapotranspiration
from .root_dynamics import Simple_Root_Dynamics as Root_Dynamics
from .leaf_dynamics import CSDM_Leaf_Dynamics as Leaf_Dynamics

class GreenLayerCrop(SimulationObject):
    """Top level object organizing the different components of the crop
    simulation for a fake crop which only simulates water use but does not
    simulate growth of biomass, etc. All crop parameters are hard-coded into the
    respective modules, only soil parameters can be passed.
    
    The FakeCropSimulation is used to put before the real crop starts in order
    to reduce the effect of water balance initialization.
            
    The processes that are implemented as embedded simulation objects consist of:
    
        1. Evapotranspiration (self.evtra)
        2. Leaf dynamics (self.lv_dynamics)
        3. Root dynamics (self.ro_dynamics)

    **Simulation parameters:**
    
    None
    
    **State variables:**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    CTRAT    Total crop transpiration                           N    cm
    DOF      Date representing the day of finish of the crop    N    -
             simulation. 
    FINISH   String representing the reason for finishing the   N    -
             simulation: maturity, harvest, leave death, etc.
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

        self.states = self.StateVariables(kiosk, CTRAT=0.0, DOF=None, FINISH=None)
            
        # assign handler for CROP_FINISH signal
        self._connect_signal(self._on_CROP_FINISH, signal=signals.crop_finish)

    #---------------------------------------------------------------------------
    @prepare_rates
    def calc_rates(self, day, drv):
        states = self.states

        # (evapo)transpiration rates
        self.evtra(day, drv)

        # Root growth
        self.ro_dynamics.calc_rates(day, drv)
        # leaf growth
        self.lv_dynamics.calc_rates(day, drv)

    #---------------------------------------------------------------------------
    @prepare_states
    def integrate(self, day, delt=1.0):
        states = self.states
        
        # Integrate states on leaves, storage organs, stems and roots
        self.ro_dynamics.integrate(day, delt)
        self.lv_dynamics.integrate(day, delt)

        # total crop transpiration (CTRAT)
        states.CTRAT += self.kiosk["TRA"]
        
    #---------------------------------------------------------------------------
    def _on_CROP_FINISH(self, day, finish):
        """Handler for setting day of finish (DOF) and reason for
        crop finishing (FINISH).
        """
        self._for_finalize["DOF"] = day
        self._for_finalize["FINISH"]= finish
