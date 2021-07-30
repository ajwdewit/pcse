# -*- coding: utf-8 -*-
# Copyright (c) 2021 Wageningen Environmental Research
# Allard de Wit and Iwan Supit (allard.dewit@wur.nl), July 2015
# Approach based on LINTUL N made by Joost Wolf

from pcse.traitlets import Float
from pcse.decorators import prepare_rates, prepare_states
from pcse.base import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject
from pcse import signals


class N_PotentialProduction(SimulationObject):
    """Provides unlimited soil N/ for potential production simulations.

    NAVAIL just remain 100 kg/ha whatever the crop takes.
    """

    class StateVariables(StatesTemplate):
        NAVAIL = Float(-99.)  # total mineral N from soil and fertiliser  kg N ha-1

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        """
        self.states = self.StateVariables(kiosk, publish=["NAVAIL"], NAVAIL=100.)

    def calc_rates(self, day, drv):
        pass

    @prepare_states
    def integrate(self, day, delt=1.0):
        self.touch()


class N_Soil_Dynamics(SimulationObject):
    """A simple module for soil N dynamics.

    This modules represents the soil as a bucket for available N consisting
    of two components: 1) a native soil supply which consists of an initial
    amount of N which will become available with a fixed fraction every day
    and 2) an external supply which is computed as an amount of N supplied
    and multiplied by a recovery fraction in order to have an effective amount of
    N that is available for crop growth.

    This module does not simulate any soil physiological processes and is only
    a book-keeping approach for N availability. On the other hand, it
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
    BG_N_SUPPLY   Background supply of N through atmospheric     SSi      |kg ha-1 d-1|
                  deposition.
    ============  ============================================= =======  ==============


    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
     NSOIL    total mineral soil N available at start of         N    [kg ha-1]
              growth period
     NAVAIL   Total mineral N from soil and fertiliser           Y    |kg ha-1|
    =======  ================================================= ==== ============

    **Rate variables**

    ==============  ================================================= ==== =============
     Name            Description                                       Pbl      Unit
    ==============  ================================================= ==== =============
    RNSOIL           Rate of change on total soil mineral N            N   |kg ha-1 d-1|
    RNAVAIL          Total change in N availability                    N   |kg ha-1 d-1|
    # Rate of fertilizer supply for N/P/K [kg/ha/day]
    FERT_N_SUPPLY    Supply of fertilizer N. This will be supplied     N   |kg ha-1 d-1|
                     by the AgroManager module through the event
                     mechanism. See the section on signals below.
    ==============  ================================================= ==== =============

    **Signals send or handled**

    `N_Soil_Dynamics` receives the following signals:
        * APPLY_N: Is received when an external input from N fertilizer
          is provided. See `_on_APPLY_N()` for details.

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
    =========  =================================== ===================  ==============
    """

    NSOILI = Float(-99.) # initial soil N amount

    class Parameters(ParamTemplate):      
        NSOILBASE = Float(-99.)  # total mineral soil N available at start of growth period [kg N/ha]
        NSOILBASE_FR = Float(-99.)  # fraction of soil mineral N coming available per day [day-1]
        BG_N_SUPPLY = Float() # Background rates of N/P/K supply [kg/ha/day]

    class StateVariables(StatesTemplate):
        NSOIL = Float(-99.)  # mineral N available from soil for crop    kg N ha-1
        NAVAIL = Float(-99.)  # total mineral N from soil and fertiliser  kg N ha-1

    class RateVariables(RatesTemplate):
        RNSOIL = Float(-99.)
        RNAVAIL = Float(-99.)
        # Rate of fertilizer supply for N [kg/ha/day]
        FERT_N_SUPPLY = Float()

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: parameter values
        """

        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk)
        self.kiosk = kiosk
        
        # INITIAL STATES
        p = self.params
        self.NSOILI = p.NSOILBASE
        self.states = self.StateVariables(kiosk, publish=["NAVAIL"], NSOIL=p.NSOILBASE, NAVAIL=0.)
        self._connect_signal(self._on_APPLY_N, signals.apply_n)
        
    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states
        p = self.params
        k = self.kiosk

        r.RNSOIL = max(0., min(p.NSOILBASE_FR * self.NSOILI, s.NSOIL))
        r.RNAVAIL = r.FERT_N_SUPPLY + p.BG_N_SUPPLY - k.RNuptake + r.RNSOIL

    @prepare_states
    def integrate(self, day, delt=1.0):
        rates = self.rates
        states = self.states

        # mineral NPK amount in the soil
        states.NSOIL += rates.RNSOIL
        # total (soil + fertilizer) NPK amount in soil
        states.NAVAIL += rates.RNAVAIL

    @prepare_rates
    def _on_APPLY_N(self, N_amount=None, N_recovery=None, *args, **kwargs):
        r = self.rates
        r.FERT_N_SUPPLY = N_amount * N_recovery
