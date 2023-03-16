# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit and Iwan Supit (allard.dewit@wur.nl), July 2015
# Approach based on LINTUL N/P/K made by Joost Wolf
import numpy as np
from pcse.traitlets import Float
from pcse.decorators import prepare_rates, prepare_states
from pcse.base import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject
from pcse import signals
from ..traitlets import Float, Int, Instance, Bool

class N_soil_dynamics_layered(SimulationObject):
    """Provides unlimited soil N/P/K for potential production simulations.

    NAVAIL just remains 100 kg/ha whatever the crop takes.
    """

    class StateVariables(StatesTemplate):
        ORGMAT = Instance(np.ndarray)
        CORG   = Instance(np.ndarray)
        NORG   = Instance(np.ndarray)
        NAVAIL = Float(-99.)  # total mineral N from soil and fertiliser  kg N ha-1

    def initialize(self, day, kiosk, parvalues):
        self.states = self.StateVariables(kiosk, publish=["NAVAIL"], NAVAIL=100.)

    def calc_rates(self, day, drv):
        pass

    @prepare_states
    def integrate(self, day, delt=1.0):
        self.touch()
