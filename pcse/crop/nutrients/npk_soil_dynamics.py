# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit and Iwan Supit (allard.dewit@wur.nl), July 2015
# Approach based on LINTUL N/P/K made by Joost Wolf

from ...traitlets import Float
from ...decorators import prepare_rates, prepare_states
from ...base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject
from ... import signals

class NPK_Soil_Dynamics(SimulationObject):
    """A simple module for soil N/P/K dynamics.

    This modules represents the soil as a bucket for available N/P/K consisting
    of two components: 1) a native soil supply which consists of an initial
    amount of N/P/K which will become available with a fixed fraction every day
    and 2) an external supply which is computed as an amount of N/P/K supplied
    and multiplied by a recovery fraction in order to have an effective amount of
    N/P/K that is available for crop growth.

    This module does not simulate any soil physiological processes and is only
    a book-keeping approach for N/P/K availability. On the other hand, it
    requires no detailed soil parameters. Only an initial soil amount, the
    fertilizer inputs, a recovery fraction and a background supply.

    **Simulation parameters**

    ============  ============================================= =======  ==============
     Name          Description                                   Type     Unit
    ============  ============================================= =======  ==============
    NSOILBASE     Base soil supply of N available through        SSi      |kg ha-1|
                  mineralisation
    NSOILBASE_FR  Fraction of base soil N that comes available   SSi        -
                  every day
    PSOILBASE     Base soil supply of N available through        SSi      |kg ha-1|
                  mineralisation
    PSOILBASE_FR  Fraction of base soil N that comes available             -
                  every day
    KSOILBASE     Base soil supply of N available through        SSi      |kg ha-1|
                  mineralisation
    KSOILBASE_FR  Fraction of base soil N that comes available   SSi        -
                  every day
    DVSNPK_STOP   Development stage after which no nutrients
                  are taken up from the soil by the crop.
    BG_N_SUPPLY   Background supply of N through atmospheric     SSi      |kg ha-1 d-1|
                  deposition.
    BG_P_SUPPLY   Background supply of P through atmospheric     SSi      |kg ha-1 d-1|
                  deposition.
    BG_K_SUPPLY   Background supply of K through atmospheric     SSi      |kg ha-1 d-1|
                  deposition.
    ============  ============================================= =======  ==============


    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
     NSOIL    total mineral soil N available at start of         N    [kg ha-1]
              growth period
     PSOIL    total mineral soil P available at start of         N    [kg ha-1]
              growth period
     KSOIL    total mineral soil K available at start of         N    [kg ha-1]
              growth period
     NAVAIL   Total mineral N from soil and fertiliser           Y    |kg ha-1|
     PAVAIL   Total mineral N from soil and fertiliser           Y    |kg ha-1|
     KAVAIL   Total mineral N from soil and fertiliser           Y    |kg ha-1|
    =======  ================================================= ==== ============

    **Rate variables**

    ==============  ================================================= ==== =============
     Name            Description                                       Pbl      Unit
    ==============  ================================================= ==== =============
    RNSOIL           Rate of change on total soil mineral N            N   |kg ha-1 d-1|
    RPSOIL           Rate of change on total soil mineral P            N   |kg ha-1 d-1|
    RKSOIL           Rate of change on total soil mineral K            N   |kg ha-1 d-1|

    RNAVAIL          Total change in N availability                    N   |kg ha-1 d-1|
    RPAVAIL          Total change in P availability                    N   |kg ha-1 d-1|
    RKAVAIL          Total change in K availability                    N   |kg ha-1 d-1|

    # Rate of fertilizer supply for N/P/K [kg/ha/day]
    FERT_N_SUPPLY    Supply of fertilizer N. This will be supplied     N   |kg ha-1 d-1|
                     by the AgroManager module through the event
                     mechanism. See the section on signals below.
    FERT_P_SUPPLY    As previous for P                                 N   |kg ha-1 d-1|
    FERT_K_SUPPLY    As previous for K                                 N   |kg ha-1 d-1|
    ==============  ================================================= ==== =============

    **Signals send or handled**

    `NPK_Soil_Dynamics` receives the following signals:
        * APPLY_NPK: Is received when an external input from N/P/K fertilizer
          is provided. See `_on_APPLY_NPK()` for details.

    **External dependencies:**

    =========  =================================== ===================  ==============
     Name       Description                         Provided by          Unit
    =========  =================================== ===================  ==============
    DVS        Crop development stage              DVS_Phenology           -
    TRA        Actual crop transpiration           Evapotranspiration     |cm|
               increase
    TRAMX      Potential crop transpiration        Evapotranspiration     |cm|
               increase
    RNUPTAKE   Rate of N uptake by the crop        NPK_Demand_Uptake     |kg ha-1 d-1|
    RPUPTAKE   Rate of P uptake by the crop        NPK_Demand_Uptake     |kg ha-1 d-1|
    RKUPTAKE   Rate of K uptake by the crop        NPK_Demand_Uptake     |kg ha-1 d-1|
    =========  =================================== ===================  ==============
    """

    NSOILI = Float(-99.) # initial soil N amount
    PSOILI = Float(-99.) # initial soil P amount
    KSOILI = Float(-99.) # initial soil K amount
    
    class Parameters(ParamTemplate):      
        NSOILBASE = Float(-99.)  # total mineral soil N available at start of growth period [kg N/ha]
        NSOILBASE_FR = Float(-99.)  # fraction of soil mineral N coming available per day [day-1]

        PSOILBASE = Float(-99.)  # total mineral soil P available at start of growth period [kg N/ha]
        PSOILBASE_FR = Float(-99.)  # fraction of soil mineral P coming available per day [day-1]
        
        KSOILBASE = Float(-99.)  # total mineral soil K available at start of growth period [kg N/ha]
        KSOILBASE_FR = Float(-99.)  # fraction of soil mineral K coming available per day [day-1]
        
        DVSNPK_STOP = Float(-99.)  # Development stage after which no nutrients are absorbed

        # Background rates of N/P/K supply [kg/ha/day]
        BG_N_SUPPLY = Float()
        BG_P_SUPPLY = Float()
        BG_K_SUPPLY = Float()

    class StateVariables(StatesTemplate):
        NSOIL = Float(-99.)  # mineral N available from soil for crop    kg N ha-1
        PSOIL = Float(-99.)  # mineral N available from soil for crop    kg N ha-1
        KSOIL = Float(-99.)  # mineral N available from soil for crop    kg N ha-1

        NAVAIL = Float(-99.)  # total mineral N from soil and fertiliser  kg N ha-1
        PAVAIL = Float(-99.)  # total mineral P from soil and fertiliser  kg N ha-1
        KAVAIL = Float(-99.)  # total mineral K from soil and fertiliser  kg N ha-1
      
    class RateVariables(RatesTemplate):
        RNSOIL = Float(-99.)
        RPSOIL = Float(-99.)
        RKSOIL = Float(-99.)
        
        RNAVAIL = Float(-99.)
        RPAVAIL = Float(-99.)
        RKAVAIL = Float(-99.)

        # Rate of fertilizer supply for N/P/K [kg/ha/day]
        FERT_N_SUPPLY = Float()
        FERT_P_SUPPLY = Float()
        FERT_K_SUPPLY = Float()

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk)
        self.kiosk = kiosk
        
        # INITIAL STATES
        p = self.params
        self.NSOILI = p.NSOILBASE
        self.PSOILI = p.PSOILBASE
        self.KSOILI = p.KSOILBASE
        
        self.states = self.StateVariables(kiosk,
            publish=["NAVAIL", "PAVAIL", "KAVAIL"],
            NSOIL=p.NSOILBASE, PSOIL=p.PSOILBASE, KSOIL=p.KSOILBASE,
            NAVAIL=0., PAVAIL=0., KAVAIL=0.)

        self._connect_signal(self._on_APPLY_NPK, signals.apply_npk)
        
    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states
        p = self.params
        
        TRA   = self.kiosk["TRA"]
        TRAMX = self.kiosk["TRAMX"]
        DVS   = self.kiosk["DVS"]
        
        NUPTR = self.kiosk["RNUPTAKE"]
        PUPTR = self.kiosk["RPUPTAKE"]
        KUPTR = self.kiosk["RKUPTAKE"]
        
        TRANRF = TRA/TRAMX
        
        if DVS < p.DVSNPK_STOP and TRANRF > 0.01 :
            NutrientLIMIT = 1.0
        else:
            NutrientLIMIT = 0.
                    
        r.RNSOIL = -max(0., min(p.NSOILBASE_FR * self.NSOILI * NutrientLIMIT, s.NSOIL))
        r.RPSOIL = -max(0., min(p.PSOILBASE_FR * self.PSOILI * NutrientLIMIT, s.PSOIL))
        r.RKSOIL = -max(0., min(p.KSOILBASE_FR * self.KSOILI * NutrientLIMIT, s.KSOIL))
               
        r.RNAVAIL = r.FERT_N_SUPPLY + p.BG_N_SUPPLY - NUPTR - r.RNSOIL
        r.RPAVAIL = r.FERT_P_SUPPLY + p.BG_P_SUPPLY - PUPTR - r.RPSOIL
        r.RKAVAIL = r.FERT_K_SUPPLY + p.BG_K_SUPPLY - KUPTR - r.RKSOIL
        
    @prepare_states
    def integrate(self, day, delt=1.0):
        rates = self.rates
        states = self.states

        # mineral NPK amount in the soil
        states.NSOIL += rates.RNSOIL
        states.PSOIL += rates.RPSOIL
        states.KSOIL += rates.RKSOIL
        
        # total (soil + fertilizer) NPK amount in soil
        states.NAVAIL += rates.RNAVAIL
        states.PAVAIL += rates.RPAVAIL
        states.KAVAIL += rates.RKAVAIL

    def _on_APPLY_NPK(self, N_amount=None, P_amount=None, K_amount=None, N_recovery=None,
                      P_recovery=None, K_recovery=None):

        r = self.rates
        r.unlock()
        r.FERT_K_SUPPLY = N_amount * N_recovery
        r.FERT_P_SUPPLY = P_amount * P_recovery
        r.FERT_K_SUPPLY = K_amount * K_recovery
        r.lock()