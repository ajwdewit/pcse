# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
# Adapted for translocation for stem dry matter by
# Iwan Supit (iwan.supit@wur.nl), June 2016

from ..traitlets import Float, Int, Instance, AfgenTrait
from ..decorators import prepare_rates, prepare_states
from ..util import limit
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject, VariableKiosk

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
    FRTRL    fraction of stem dry weight that can be        SCr        -
             remobilized for grain growth
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
    TRANSL   Amount of stem biomass that can be tarnslocated    Y   |kg ha-1 d-1|
    DVR      Development rate                                   N     -
    =======  ================================================= ==== ============
    
    **Signals send or handled**
    
    None
    
    **External dependencies:**
    
    =======  =================================== =================       ============
     Name     Description                         Provided by             Unit
    =======  =================================== =================       ============
    DVS      Crop development stage              DVS_Phenology             -
    DVR      Crop development rate               DVS_Phenology             -
    ADMI     Above-ground dry matter             CropSimulation          |kg ha-1 d-1|
             increase
    FR       Fraction biomass to roots           DVS_Partitioning          -
    FS       Fraction biomass to stems           DVS_Partitioning          -
    REDUCTL  Translocation reduction             storage_organ_dynamics  |kg ha-1 d-1|
    =======  =================================== ======================  ============
    """

    class Parameters(ParamTemplate):      
        RDRSTB = AfgenTrait()
        SSATB  = AfgenTrait()
        ISINK  = Float(-99.)
        FRTRL  = Float(-99.)
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
        TRANSL = Float(-99.)
        
    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.rates  = self.RateVariables(kiosk, publish=["DRST","TRANSL"])
        self.kiosk  = kiosk

        # INITIAL STATES
        params = self.params
        # Set initial stem biomass
        FS = self.kiosk["FS"]
        FR = self.kiosk["FR"]
        WST  = (params.TDWI * (1-FR)) * FS
        DWST = 0.
        TWST = WST + DWST
        # Initial Stem Area Index
        DVS = self.kiosk["DVS"]
        SAI = WST * params.SSATB(DVS)

        self.states = self.StateVariables(kiosk, publish=["TWST","WST","SAI"],
                                          WST=WST, DWST=DWST, TWST=TWST, SAI=SAI)

    @prepare_rates
    def calc_rates(self, day, drv):
        rates  = self.rates
        states = self.states
        params = self.params
        
        DVS  = self.kiosk["DVS"]
        DVR  = self.kiosk["DVR"]
        FS   = self.kiosk["FS"]
        ADMI = self.kiosk["ADMI"]

        # Calcualte the amount of stem dry matter that can be translocated
        if DVS > 1.0:
            rates.TRANSL = (states.WST + states.DWST) * DVR * params.FRTRL
        else:
            rates.TRANSL = 0.

        # Growth/death rate stems
        rates.GRST = ADMI * FS
        rates.DRST = params.RDRSTB(DVS) * states.WST
        rates.GWST = rates.GRST - rates.DRST - rates.TRANSL


    @prepare_states
    def integrate(self, day):
        params = self.params
        rates = self.rates
        states = self.states

        # in case of sink limitation not all assimilates
        # can be translocated to the storage organs
        # consequently the stem growth rate has to be increased
        # with the translocation reduction REDUCTL
        if params.ISINK ==1:
            if self.kiosk["REDUCTL"] is None:
                REDUCTL = 0
            else:
                REDUCTL = self.kiosk["REDUCTL"]
        else:
            REDUCTL = 0.

        # Stem biomass (living, dead, total)
        states.WST += (rates.GWST + REDUCTL)
        states.DWST += rates.DRST
        states.TWST = states.WST + states.DWST

        # Calculate Stem Area Index (SAI)
        DVS = self.kiosk["DVS"]
        states.SAI = states.WST * params.SSATB(DVS)
