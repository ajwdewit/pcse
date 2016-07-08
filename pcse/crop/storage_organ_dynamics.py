# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014

from ..traitlets import Float, Int, Instance
from ..decorators import prepare_rates, prepare_states
from ..util import limit
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject, VariableKiosk

class WOFOST_Storage_Organ_Dynamics(SimulationObject):
    """Implementation of storage organ dynamics.
    
    Storage organs are the most simple component of the plant in WOFOST and
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
        self.rates  = self.RateVariables(kiosk)
        self.kiosk = kiosk
        
        # INITIAL STATES
        params = self.params
        # Initial storage organ biomass
        FO = self.kiosk["FO"]
        FR = self.kiosk["FR"]
        WSO  = (params.TDWI * (1-FR)) * FO
        DWSO = 0.
        TWSO = WSO + DWSO
        # Initial Pod Area Index
        PAI = WSO * params.SPA

        self.states = self.StateVariables(kiosk, publish=["TWSO","WSO","PAI"],
                                          WSO=WSO, DWSO=DWSO, TWSO=TWSO,
                                          PAI=PAI)

    @prepare_rates
    def calc_rates(self, day, drv):
        rates  = self.rates
        states = self.states
        params = self.params
        
        FO = self.kiosk["FO"]
        ADMI = self.kiosk["ADMI"]

        # Growth/death rate organs
        rates.GRSO = ADMI * FO
        rates.DRSO = 0.0
        rates.GWSO = rates.GRSO - rates.DRSO

    @prepare_states
    def integrate(self, day, delt=1.0):
        params = self.params
        rates = self.rates
        states = self.states

        # Stem biomass (living, dead, total)
        states.WSO += rates.GWSO
        states.DWSO += rates.DRSO
        states.TWSO = states.WSO + states.DWSO

        # Calculate Pod Area Index (SAI)
        states.PAI = states.WSO * params.SPA
