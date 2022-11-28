# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit and Iwan Supit (allard.dewit@wur.nl), July 2015
# Approach based on LINTUL N/P/K made by Joost Wolf

from pcse.traitlets import Float
from pcse.decorators import prepare_rates, prepare_states
from pcse.base import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject
from pcse import signals


class NPK_PotentialProduction(SimulationObject):
    """Provides unlimited soil N/P/K for potential production simulations.

    NAVAIL, KAVAIL and PAVAIL just remain 100 kg/ha whatever the crop takes.
    """

    class StateVariables(StatesTemplate):
        NAVAIL = Float(-99.)  # total mineral N from soil and fertiliser  kg N ha-1
        PAVAIL = Float(-99.)  # total mineral P from soil and fertiliser  kg N ha-1
        KAVAIL = Float(-99.)  # total mineral K from soil and fertiliser  kg N ha-1

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        """
        self.states = self.StateVariables(kiosk, publish=["NAVAIL", "PAVAIL", "KAVAIL"],
                                          NAVAIL=100., PAVAIL=100., KAVAIL=100.)

    def calc_rates(self, day, drv):
        pass

    @prepare_states
    def integrate(self, day, delt=1.0):
        self.touch()


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
    NAVAILI       Initial N available in the N pool              SSi      |kg ha-1|
    PAVAILI       Initial P available in the P pool              SSi      |kg ha-1|
    KAVAILI       Initial K available in the K pool              SSi      |kg ha-1|
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
    RNuptake   Rate of N uptake by the crop        NPK_Demand_Uptake     |kg ha-1 d-1|
    RPuptake   Rate of P uptake by the crop        NPK_Demand_Uptake     |kg ha-1 d-1|
    RKuptake   Rate of K uptake by the crop        NPK_Demand_Uptake     |kg ha-1 d-1|
    =========  =================================== ===================  ==============
    """

    NSOILI = Float(-99.) # initial soil N amount
    PSOILI = Float(-99.) # initial soil P amount
    KSOILI = Float(-99.) # initial soil K amount

    # placeholders for FERT_N/P/K_SUPPLY
    _FERT_N_SUPPLY = Float(0.)
    _FERT_P_SUPPLY = Float(0.)
    _FERT_K_SUPPLY = Float(0.)

    class Parameters(ParamTemplate):
        NSOILBASE = Float(-99.)  # total mineral soil N available at start of growth period [kg N/ha]
        NSOILBASE_FR = Float(-99.)  # fraction of soil mineral N coming available per day [day-1]

        PSOILBASE = Float(-99.)  # total mineral soil P available at start of growth period [kg N/ha]
        PSOILBASE_FR = Float(-99.)  # fraction of soil mineral P coming available per day [day-1]
        
        KSOILBASE = Float(-99.)  # total mineral soil K available at start of growth period [kg N/ha]
        KSOILBASE_FR = Float(-99.)  # fraction of soil mineral K coming available per day [day-1]

        # Initial values of available nutrients which is different from the previous ones
        #
        NAVAILI = Float()
        PAVAILI = Float()
        KAVAILI = Float()

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
            NAVAIL=p.NAVAILI, PAVAIL=p.PAVAILI, KAVAIL=p.KAVAILI)

        self._connect_signal(self._on_APPLY_NPK, signals.apply_npk)
        
    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states
        p = self.params
        k = self.kiosk

        # Rate of supplied N/P/K
        r.FERT_N_SUPPLY = self._FERT_N_SUPPLY
        r.FERT_P_SUPPLY = self._FERT_P_SUPPLY
        r.FERT_K_SUPPLY = self._FERT_K_SUPPLY
        self._FERT_N_SUPPLY = 0.
        self._FERT_P_SUPPLY = 0.
        self._FERT_K_SUPPLY = 0.

        r.RNSOIL = -max(0., min(p.NSOILBASE_FR * self.NSOILI, s.NSOIL))
        r.RPSOIL = -max(0., min(p.PSOILBASE_FR * self.PSOILI, s.PSOIL))
        r.RKSOIL = -max(0., min(p.KSOILBASE_FR * self.KSOILI, s.KSOIL))

        # Check uptake rates from crop, if a crop is actually growing
        RNuptake = k.RNuptake if "RNuptake" in self.kiosk else 0.
        RPuptake = k.RNuptake if "RPuptake" in self.kiosk else 0.
        RKuptake = k.RNuptake if "RKuptake" in self.kiosk else 0.

        r.RNAVAIL = r.FERT_N_SUPPLY + p.BG_N_SUPPLY - RNuptake - r.RNSOIL
        r.RPAVAIL = r.FERT_P_SUPPLY + p.BG_P_SUPPLY - RPuptake - r.RPSOIL
        r.RKAVAIL = r.FERT_K_SUPPLY + p.BG_K_SUPPLY - RKuptake - r.RKSOIL
        
    @prepare_states
    def integrate(self, day, delt=1.0):
        rates = self.rates
        states = self.states

        # mineral NPK amount in the soil
        states.NSOIL += rates.RNSOIL * delt
        states.PSOIL += rates.RPSOIL * delt
        states.KSOIL += rates.RKSOIL * delt
        
        # total (soil + fertilizer) NPK amount in soil
        states.NAVAIL += rates.RNAVAIL * delt
        states.PAVAIL += rates.RPAVAIL * delt
        states.KAVAIL += rates.RKAVAIL * delt

    def _on_APPLY_NPK(self, N_amount=None, P_amount=None, K_amount=None, N_recovery=None,
                      P_recovery=None, K_recovery=None):

        self._FERT_N_SUPPLY = N_amount * N_recovery
        self._FERT_P_SUPPLY = P_amount * P_recovery
        self._FERT_K_SUPPLY = K_amount * K_recovery