# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
from copy import deepcopy

from ..traitlets import Float, Int, Instance, AfgenTrait
from ..decorators import prepare_rates, prepare_states
from ..util import limit, merge_dict
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject, VariableKiosk

#-------------------------------------------------------------------------------
class WOFOST_Root_Dynamics(SimulationObject):
    """Root biomass dynamics and rooting depth.
    
    Root growth and root biomass dynamics in WOFOST are separate processes,
    with the only exception that root growth stops when no more biomass is sent
    to the root system.
    
    Root biomass increase results from the assimilates partitioned to
    the root system. Root death is defined as the current root biomass
    multiplied by a relative death rate (`RDRRTB`). The latter as a function
    of the development stage (`DVS`).
    
    Increase in root depth is a simple linear expansion over time unti the
    maximum rooting depth (`RDM`) is reached.
    
    **Simulation parameters**
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    RDI      Initial rooting depth                          SCr      cm
    RRI      Daily increase in rooting depth                SCr      |cm day-1|
    RDMCR    Maximum rooting depth of the crop              SCR      cm
    RDMSOL   Maximum rooting depth of the soil              SSo      cm
    TDWI     Initial total crop dry weight                  SCr      |kg ha-1|
    IAIRDU   Presence of air ducts in the root (1) or       SCr      -
             not (0)
    RDRRTB   Relative death rate of roots as a function     TCr      -
             of development stage
    =======  ============================================= =======  ============
    

    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    RD       Current rooting depth                              Y     cm
    RDM      Maximum attainable rooting depth at the minimum    N     cm
             of the soil and crop maximum rooting depth
    WRT      Weight of living roots                             Y     |kg ha-1|
    DWRT     Weight of dead roots                               N     |kg ha-1|
    TWRT     Total weight of roots                              Y     |kg ha-1|
    =======  ================================================= ==== ============

    **Rate variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    RR       Growth rate root depth                             N    cm
    GRRT     Growth rate root biomass                           N   |kg ha-1 d-1|
    DRRT     Death rate root biomass                            N   |kg ha-1 d-1|
    GWRT     Net change in root biomass                         N   |kg ha-1 d-1|
    =======  ================================================= ==== ============
    
    **Signals send or handled**
    
    None
    
    **External dependencies:**
    
    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    DVS      Crop development stage              DVS_Phenology       -
    DMI      Total dry matter                    CropSimulation     |kg ha-1 d-1|
             increase
    FR       Fraction biomass to roots           DVS_Partitioning    - 
    =======  =================================== =================  ============
    """

    class Parameters(ParamTemplate):
        """Traits-based class for storing rooting depth parameters
        """
        RDI    = Float(-99.)    
        RRI    = Float(-99.)
        RDMCR  = Float(-99.)
        RDMSOL = Float(-99.)
        TDWI   = Float(-99.)
        IAIRDU = Float(-99)
        RDRRTB = AfgenTrait()
                    
    class RateVariables(RatesTemplate):
        RR   = Float(-99.)
        GRRT = Float(-99.)
        DRRT = Float(-99.)
        GWRT = Float(-99.)

    class StateVariables(StatesTemplate):
        RD   = Float(-99.)
        RDM  = Float(-99.)
        WRT  = Float(-99.)
        DWRT = Float(-99.)
        TWRT = Float(-99.)
        
    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish="DRRT")
        self.kiosk = kiosk
        
        # INITIAL STATES
        params = self.params
        # Initial root depth states
        rdmax = max(params.RDI, min(params.RDMCR, params.RDMSOL))
        RDM = rdmax
        RD = params.RDI
        # initial root biomass states
        FR = self.kiosk["FR"]
        WRT  = params.TDWI * FR
        DWRT = 0.
        TWRT = WRT + DWRT

        self.states = self.StateVariables(kiosk, publish=["RD","WRT","TWRT"],
                                          RD=RD, RDM=RDM, WRT=WRT, DWRT=DWRT,
                                          TWRT=TWRT)
    @prepare_rates
    def calc_rates(self, day, drv):
        params = self.params
        rates = self.rates
        states = self.states

        # Increase in root biomass
        DMI = self.kiosk["DMI"]
        DVS = self.kiosk["DVS"]
        FR = self.kiosk["FR"]
        rates.GRRT = FR * DMI
        rates.DRRT = states.WRT * params.RDRRTB(DVS)
        rates.GWRT = rates.GRRT - rates.DRRT
        
        # Increase in root depth
        rates.RR = min((states.RDM - states.RD), params.RRI)
        # Do not let the roots growth if partioning to the roots
        # (variable FR) is zero.
        if FR == 0.:
            rates.RR = 0.
    
    @prepare_states
    def integrate(self, day, delt=1.0):
        rates = self.rates
        states = self.states

        # Dry weight of living roots
        states.WRT += rates.GWRT
        # Dry weight of dead roots
        states.DWRT += rates.DRRT
        # Total weight dry + living roots
        states.TWRT = states.WRT + states.DWRT

        # New root depth
        states.RD += rates.RR

class Simple_Root_Dynamics(SimulationObject):
    """Simple class for linear root growth.
    
    Increase in root depth is a simple linear expansion over time until the
    maximum rooting depth (`RDM`) is reached.
    
    **Simulation parameters**
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    RDI      Initial rooting depth                          SCr      cm
    RRI      Daily increase in rooting depth                SCr      |cm day-1|
    RDMCR    Maximum rooting depth of the crop              SCR      cm
    RDMSOL   Maximum rooting depth of the soil              SSo      cm
    =======  ============================================= =======  ============
    

    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    RD       Current rooting depth                              Y     cm
    RDM      Maximum attainable rooting depth at the minimum    N     cm
             of the soil and crop maximum rooting depth
    =======  ================================================= ==== ============

    **Rate variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    RR       Growth rate root depth                             N    cm
    =======  ================================================= ==== ============
    
    **Signals send or handled**
    
    None
    
    **External dependencies:**
    
    None
    """

    class Parameters(ParamTemplate):
        """Traits-based class for storing rooting depth parameters
        """
        RDI    = Float(-99.)    
        RRI    = Float(-99.)
        RDMCR  = Float(-99.)
        RDMSOL = Float(-99.)
                    
    class RateVariables(RatesTemplate):
        RR   = Float(-99.)

    class StateVariables(StatesTemplate):
        RD   = Float(-99.)
        RDM  = Float(-99.)
        
    def initialize(self, day, kiosk, cropdata, soildata):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        :param soildata: dictionary with WOFOST soildata key/value pairs
        """
        
        # Merge cropdata and soildata dictionaries in order to pass them to
        # the Parameters class.
        parvalues = merge_dict(cropdata, soildata)
        
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk)
        self.kiosk = kiosk
        
        # INITIAL STATES
        params = self.params
        # Initial root depth states
        rdmax = max(params.RDI, min(params.RDMCR, params.RDMSOL))
        RDM = rdmax
        RD = params.RDI

        self.states = self.StateVariables(kiosk, publish=["RD"],
                                          RD=RD, RDM=RDM)
    @prepare_rates
    def calc_rates(self, day, drv):
        params = self.params
        rates = self.rates
        states = self.states
        
        # Increase in root depth
        rates.RR = min((states.RDM - states.RD), params.RRI)
    
    @prepare_states
    def integrate(self, day, delt=1.0):
        rates = self.rates
        states = self.states

        # New root depth
        states.RD += rates.RR
