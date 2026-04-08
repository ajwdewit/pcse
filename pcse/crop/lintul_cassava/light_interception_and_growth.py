# -*- coding: utf-8 -*-
# Herman Berghuijs (herman.berghuijs@wur.nl), Allard de Wit (allard.dewit@wur.nl), Tom Schut (tom.schut@wur.nl)
# February 2026

import numpy as np
from pcse.base import ParamTemplate, RatesTemplate, SimulationObject, StatesTemplate
from pcse.traitlets import Float
from pcse.util import AfgenTrait

J_to_MJ = 1e-6

class light_interception_and_growth(SimulationObject):
    """
    Class to simulate the light interception and subsequent production of new dry matter

    This class calculates the amount of light that is intercepted using Lambert-Beer's law. The daily potential dry
    matter production is calculated as the product of the amount of intercepted light and the light use efficiency.
    The actual dry matter production is also calculated. This can be smaller than the potential dry matter production
    in case of water stress, nutrient stress or suboptimal temperatures.

    **Simulation parameters**

    =================  ================================================  ======  ===========================
    Name               Description                                       Type    Unit
    =================  ================================================  ======  ===========================
    K_EXT              Extinction coefficient of light in the canopy      SCr     m2 ground m-2 leaf
    FPAR               Fraction of PAR in global radiation                SCr     MJ PAR MJ-1 radiation
    LUE_OPT            Light use efficiency in absence of water or
                       nutrient stress and under optimal temperatures     SCr     g DM MJ-1 PAR
    TTB                Reduction rate of light use efficiency due to
                       suboptimal temperatures as a table function
                       of temperature                                     TCr     -
    =================  ================================================  ======  ===========================

    **State variables**

    =================  ==============================================  ======  ===========================
    Name                Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    PAR                 Total amount of intercepted PAR                 Y       MJ PAR m-2 ground
    =================  ==============================================  ======  ===========================

    **Rate variables**

    =================  ==============================================  ======  ===========================
    Name                Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    RPAR                Daily rate of PAR interception                  N       MJ PAR m-2 ground d-1
    =================  ==============================================  ======  ===========================

    **Auxillary variables**

    =================  ==============================================  ======  ===========================
    Name                Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    GTOTAL              Rate of total dry matter production             Y       g DM m-2 ground d-1
    =================  ==============================================  ======  ===========================

    """
    class Parameters(ParamTemplate):
        FPAR = Float()
        K_EXT = Float()
        LUE_OPT = Float()
        TTB = AfgenTrait()

    class RateVariables(RatesTemplate):
        GTOTAL = Float()
        RPAR = Float()

    class StateVariables(StatesTemplate):
        PAR = Float()

    def initialize(self, day, kiosk, parvalues, delt = 1.):
        PAR = 0.
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk,
                                        publish = ["GTOTAL"])
        self.states = self.StateVariables(kiosk,
                                          publish = ["PAR"],
                                          PAR = PAR)

    def calc_rates(self,  day, drv, delt = 1.):
        p = self.params
        r = self.rates
        k = self.kiosk

        DTR = drv.IRRAD * J_to_MJ
        RPAR = p.FPAR * DTR

        # Light interception and total crop growth rate.
        PARINT = RPAR * (1 - np.exp(-p.K_EXT * k.LAI))  # MJ m-2 d-1

        TEMPRF = p.TTB(drv.TEMP)
        LUE = p.LUE_OPT * TEMPRF # g DM m-2 d-1

        # When water stress is more severe or nutrient is stress is more severe
        if k.RFTRA <= k.NPKI:
            GTOTAL = LUE * PARINT * k.RFTRA * (1 - k.DORMANCY)  # g DM m-2 d-1
        else:
            GTOTAL = LUE * PARINT * k.NPKI * (1 - k.DORMANCY)  # g DM m-2 d-1

        r.RPAR = RPAR
        r.GTOTAL = GTOTAL

    def integrate(self, day, drv, delt = 1):
        r = self.rates
        s = self.states
        s.PAR += r.RPAR