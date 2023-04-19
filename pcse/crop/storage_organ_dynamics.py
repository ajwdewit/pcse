# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014

from ..traitlets import Float
from ..decorators import prepare_rates, prepare_states
from ..base import ParamTemplate, StatesTemplate, RatesTemplate, SimulationObject


class WOFOST_Storage_Organ_Dynamics(SimulationObject):
    """Implementation of storage organ dynamics.
    
    Storage organs are the simplest component of the plant in WOFOST and
    consist of a static pool of biomass. Growth of the storage organs is the
    result of assimilate partitioning. Death of storage organs is not
    implemented and the corresponding rate variable (DRSO) is always set to
    zero.
    
    Pods are green elements of the plant canopy and can as such contribute
    to the total photosynthetic active area. This is expressed as the Pod
    Area Index which is obtained by multiplying pod biomass with a fixed
    Specific Pod Area (SPA).

    **Simulation parameters**
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    TDWI     Initial total crop dry weight                  SCr      |kg ha-1|
    SPA      Specific Pod Area                              SCr      |ha kg-1|
    =======  ============================================= =======  ============    

    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    PAI      Pod Area Index                                     Y     -
    WSO      Weight of living storage organs                    Y     |kg ha-1|
    DWSO     Weight of dead storage organs                      N     |kg ha-1|
    TWSO     Total weight of storage organs                     Y     |kg ha-1|
    =======  ================================================= ==== ============

    **Rate variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    GRSO     Growth rate storage organs                         N   |kg ha-1 d-1|
    DRSO     Death rate storage organs                          N   |kg ha-1 d-1|
    GWSO     Net change in storage organ biomass                N   |kg ha-1 d-1|
    =======  ================================================= ==== ============
    
    **Signals send or handled**
    
    None
    
    **External dependencies**
    
    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    ADMI     Above-ground dry matter             CropSimulation     |kg ha-1 d-1|
             increase
    FO       Fraction biomass to storage organs  DVS_Partitioning    - 
    FR       Fraction biomass to roots           DVS_Partitioning    - 
    =======  =================================== =================  ============
    """

    class Parameters(ParamTemplate):      
        SPA  = Float(-99.)
        TDWI = Float(-99.)

    class StateVariables(StatesTemplate):
        WSO  = Float(-99.) # Weight living storage organs
        DWSO = Float(-99.) # Weight dead storage organs
        TWSO = Float(-99.) # Total weight storage organs
        PAI  = Float(-99.) # Pod Area Index

    class RateVariables(RatesTemplate):
        GRSO = Float(-99.)
        DRSO = Float(-99.)
        GWSO = Float(-99.)
        
    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.rates  = self.RateVariables(kiosk, publish = ["GRSO"])
        self.kiosk = kiosk
        
        p = self.params
        k = self.kiosk

        WSO  = (p.TDWI * (1 - k.FR)) * k.FO
        s = dict(
            WSO = WSO,
            DWSO = 0.,
            TWSO = WSO,
            # Initial Pod Area Index
            PAI = WSO * p.SPA
        )
        self.states = self.StateVariables(kiosk, publish=["TWSO","WSO","PAI"], **s)

    @prepare_rates
    def calc_rates(self, day, drv):
        r  = self.rates
        k = self.kiosk
        
        # Growth/death rate organs
        r.GRSO = k.ADMI * k.FO
        r.DRSO = 0.0
        r.GWSO = r.GRSO - r.DRSO + k.REALLOC_SO

    @prepare_states
    def integrate(self, day, delt=1.0):
        p = self.params
        r = self.rates
        s = self.states

        s.WSO += r.GWSO
        s.DWSO += r.DRSO
        s.TWSO = s.WSO + s.DWSO
        s.PAI = s.WSO * p.SPA

    @prepare_states
    def _set_variable_WSO(self, nWSO):
        s = self.states
        p = self.params
        oWSO, oTWSO, oPAI = s.WSO, s.TWSO, s.PAI
        s.WSO = nWSO
        s.TWSO = s.DWSO + nWSO
        s.PAI = s.WSO * p.SPA

        increments = {"WSO": s.WSO - oWSO,
                      "PAI": s.PAI - oPAI,
                      "TWSO": s.TWSO - oTWSO}
        return increments