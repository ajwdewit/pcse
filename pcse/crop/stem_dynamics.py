# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014

from ..traitlets import Float
from ..decorators import prepare_rates, prepare_states
from ..util import AfgenTrait
from ..base import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject

class WOFOST_Stem_Dynamics(SimulationObject):
    """Implementation of stem biomass dynamics.
    
    Stem biomass increase results from the assimilates partitioned to
    the stem system. Stem death is defined as the current stem biomass
    multiplied by a relative death rate (`RDRSTB`). The latter as a function
    of the development stage (`DVS`).
    
    Stems are green elements of the plant canopy and can as such contribute
    to the total photosynthetic active area. This is expressed as the Stem
    Area Index which is obtained by multiplying stem biomass with the
    Specific Stem Area (SSATB), which is a function of DVS.

    **Simulation parameters**:
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    TDWI     Initial total crop dry weight                  SCr       |kg ha-1|
    RDRSTB   Relative death rate of stems as a function     TCr       -
             of development stage
    SSATB    Specific Stem Area as a function of            TCr       |ha kg-1|
             development stage
    =======  ============================================= =======  ============
    

    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    SAI      Stem Area Index                                    Y     -
    WST      Weight of living stems                             Y     |kg ha-1|
    DWST     Weight of dead stems                               N     |kg ha-1|
    TWST     Total weight of stems                              Y     |kg ha-1|
    =======  ================================================= ==== ============

    **Rate variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    GRST     Growth rate stem biomass                           N   |kg ha-1 d-1|
    DRST     Death rate stem biomass                            N   |kg ha-1 d-1|
    GWST     Net change in stem biomass                         N   |kg ha-1 d-1|
    =======  ================================================= ==== ============
    
    **Signals send or handled**
    
    None
    
    **External dependencies:**
    
    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    DVS      Crop development stage              DVS_Phenology       -
    ADMI     Above-ground dry matter             CropSimulation     |kg ha-1 d-1|
             increase
    FR       Fraction biomass to roots           DVS_Partitioning    - 
    FS       Fraction biomass to stems           DVS_Partitioning    - 
    =======  =================================== =================  ============
    """

    class Parameters(ParamTemplate):      
        RDRSTB = AfgenTrait()
        SSATB  = AfgenTrait()
        TDWI   = Float(-99.)

    class StateVariables(StatesTemplate):
        WST  = Float(-99.)
        DWST = Float(-99.)
        TWST = Float(-99.)
        SAI  = Float(-99.) # Stem Area Index

    class RateVariables(RatesTemplate):
        GRST = Float(-99.)
        DRST = Float(-99.)
        GWST = Float(-99.)
        
    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """
        
        self.params = self.Parameters(parvalues)
        self.rates  = self.RateVariables(kiosk, publish=["DRST", "GRST"])
        self.kiosk  = kiosk

        p = self.params
        k = self.kiosk

        WST = (p.TDWI * (1 - k.FR)) * k.FS
        s = dict(
            WST  = WST,
            DWST = 0.,
            TWST = WST,
            SAI = WST * p.SSATB(k.DVS)
        )

        self.states = self.StateVariables(kiosk, publish=["TWST","WST","SAI"], **s)

    @prepare_rates
    def calc_rates(self, day, drv):
        r  = self.rates
        s = self.states
        p = self.params
        k = self.kiosk
        
        # Growth/death rate stems
        r.GRST = k.ADMI * k.FS
        r.DRST = p.RDRSTB(k.DVS) * s.WST
        r.GWST = r.GRST - r.DRST - k.REALLOC_ST

    @prepare_states
    def integrate(self, day, delt=1.0):
        p = self.params
        r = self.rates
        s = self.states
        k = self.kiosk

        # Stem biomass (living, dead, total)
        s.WST += r.GWST
        s.DWST += r.DRST
        s.TWST = s.WST + s.DWST

        # Calculate Stem Area Index (SAI)
        s.SAI = s.WST * p.SSATB(k.DVS)

    @prepare_states
    def _set_variable_WST(self, nWST):
        s = self.states
        p = self.params
        k = self.kiosk

        oWST = s.WST
        oTWST = s.TWST
        oSAI = s.SAI
        s.WST = nWST
        s.TWST = s.DWST + nWST
        s.SAI = s.WST * p.SSATB(k.DVS)

        increments = {"WST": s.WST - oWST,
                      "SAI": s.SAI - oSAI,
                      "TWST": s.TWST - oTWST}
        return increments
