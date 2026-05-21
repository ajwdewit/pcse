# -*- coding: utf-8 -*-
# Herman Berghuijs (herman.berghuijs@wur.nl), Allard de Wit (allard.dewit@wur.nl), Tom Schut (tom.schut@wur.nl)
# March 2026

from pcse.base import ParamTemplate, RatesTemplate, SimulationObject, StatesTemplate
from pcse.traitlets import Float

import pcse.signals

kg_to_g = 1000
ha_to_m2 = 10000

class soil_nutrient_dynamics_PP(SimulationObject):
    """
    Class to simulate N, P, and K dynamics in soil under non-nutrient limited conditions

    This class simulates a fake water balance of the dynamics of the amounts of N, P, and K that are available for
    root uptake. Arbitrary values are assigned to these amount at each simulated day. This class is used to run LINTUL
    Cassava under non-nutrient limited conditions.

    ** Auxillary variables **

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    NAVAIL             Amount of available nitrogen in soil            Y       g N m-2 ground d-1
    KAVAIL             Amount of available potassium in soil           Y       g K m-2 ground d-1
    PAVAIL             Amount of available phosphorus in soil          Y       g P m-2 ground d-1
    =================  ==============================================  ======  ===========================
    """
    class Parameters(ParamTemplate):
        pass

    class StateVariables(StatesTemplate):
        NAVAIL = Float()
        PAVAIL = Float()
        KAVAIL = Float()

    class RateVariables(RatesTemplate):
        pass

    def initialize(self, day, kiosk, parameters):
        self.kiosk = kiosk
        self.rates = self.RateVariables(kiosk, publish = [])
        self.params = self.Parameters(parameters)
        p = self.params
        NAVAIL = 1000.
        PAVAIL = 1000.
        KAVAIL = 1000.
        self.states = self.StateVariables(kiosk,
                                          NAVAIL = NAVAIL,
                                          PAVAIL = PAVAIL,
                                          KAVAIL = KAVAIL,
                                          publish = ["NAVAIL",
                                                     "PAVAIL",
                                                     "KAVAIL"]
                                          )
    def calc_rates(self, day, drv):
        pass

    def integrate(self, day, drv):
        self.touch()

class soil_nutrient_dynamics(SimulationObject):
    """
    Class to simulate N, P, and K dynamics in soil under nutrient limited conditions

    This class simulates the dynamics of the amount of N, P, and K. It distiguishes between the amounts of N, P, and K
    that are available for plant uptake and those that are not available for plant uptake but become available over
    time. The latter amounts of N, P, and K are further subdivided into N, P, and K in the soil supply and in the
    applied fertilizers.

    ** Simulation parameters **

    =================  ==============================================  ======  ===========================
    Name               Description                                     Type     Unit
    =================  ==============================================  ======  ===========================
    NMINI              Initial total amount of nitrogen in soil        SCr      g N m-2 ground
    KMINI              Initial total amount of potassium in soil       SCr      g K m-2 ground
    K_WATER            Transpiration reduction factor at which the
                       nutrient uptake rate is reduced by a half due
                       to water stress                                 SCr     -
    PMINI              Initial total amount of phosphorus in soil      SCr      g P m-2 ground
    RTNMINF            Relative rate at which nitrogen becomes
                       available from applied fertilizers              SCr      d-1
    RTNMINS            Relative rate at which nitrogen becomes
                       available from the soil supply                  SCr      d-1
    RTKMINF            Relative rate at which potassium becomes
                       available from applied fertilizers              SCr      d-1
    RTKMINS            Relative rate at which potassium becomes
                       available from the soil supply                  SCr      d-1
    RTPMINF            Relative rate at which phosphorus becomes
                       available from applied fertilizers              SCr      d-1
    RTPMINS            Relative rate at which phosphorus becomes
                       available from the soil supply                  SCr      d-1
    =================  ==============================================  ======  ===========================

    ** State variables **

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    NMINF              Amount of nitrogen in applied fertilizers
                       and not yet available for root uptake           N       g N m-2 ground
    NMINS              Amount of nitrogen in the soil supply and not
                       yet available for root uptake                   N       g N m-2 ground
    NMINT              Amount of nitrogen in the soil that is
                       available for root uptake.                      N       g N m-2 ground
    KMINF              Amount of potassium in applied fertilizers
                       and not yet available for root uptake           N       g K m-2 ground
    KMINS              Amount of potassium in the soil supply and not
                       yet available for root uptake                   N       g K m-2 ground
    KMINT              Amount of potassium in the soil that is
                       available for root uptake.                      N       g K m-2 ground
    PMINF              Amount of phosphorus in applied fertilizers
                       and not yet available for root uptake           N       g P m-2 ground
    PMINS              Amount of phosphorus in the soil supply and not
                       yet available for root uptake                   N       g P m-2 ground
    PMINT              Amount of phosphorus in the soil that is
                       available for root uptake.                      N       g P m-2 ground
    =================  ==============================================  ======  ===========================

    ** Rate variables **

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    RNMINF             Rate of change of nitrogen amount in applied
                       fertilizers                                     N       g N m-2 ground d-1
    RNMINS             Rate of change of nitrogen in the soil supply   N       g N m-2 ground d-1
    RNMINT             Rate of change of available nitrogen amount     N       g N m-2 ground d-1
    RKMINF             Rate of change of potassium amount in applied
                       fertilizers                                     N       g K m-2 ground d-1
    RKMINS             Rate of change of potassium in the soil supply  N       g K m-2 ground d-1
    RKMINT             Rate of change of available potassium amount    N       g K m-2 ground d-1
    RPMINF             Rate of change of phosphorus amount in applied
                       fertilizers                                     N       g P m-2 ground d-1
    RPMINS             Rate of change of phosphorus in the soil supply N       g P m-2 ground d-1
    RPMINT             Rate of change of available phosphorus amount   N       g P m-2 ground d-1
    =================  ==============================================  ======  ===========================

    This class is a Python implementation of the calculations related to the amounts of N, P, and K in the
    soil in the R function nutrientdyn in the R version of the model LINTUL Cassava NPK (Adiele et al., 2022; Ezui et
    al., 2018).

    Authors nutrientdyn:             Rob van den Beuken, Tom Schut
    Authors Python implementation:   Herman Berghuijs, Allard de Wit, Tom Schut

    References:
    Adiele J.G., Schut A.G.T., Ezui K.S., Giller K.E. (2022) LINTUL-Cassava-NPK: A simulation
    model for nutrient-limited cassava growth. Field Crops Research 281: ARTN 108488.
    https://doi.org/10.1007/s13593-020-00649-w

    Ezui K.S., Leffelaar P.A., Franke A.C., Mando A., Giller K.E. (2018) Simulating drought impact
    and mitigation in cassava using the LINTUL model. Field Crops Research 219: 256-272.
    https://doi.org/10.1016/j.fcr.2018.01.033
    """

    # Placeholders
    _RFERTN = None
    _RFERTP = None
    _RFERTK = None
    _NRECOV = None
    _PRECOV = None
    _KRECOV = None

    class Parameters(ParamTemplate):
        RTNMINS = Float()
        RTPMINS = Float()
        RTKMINS = Float()
        RTNMINF = Float()
        RTPMINF = Float()
        RTKMINF = Float()
        NMINI = Float()
        PMINI = Float()
        KMINI = Float()
        K_WATER = Float()

    class StateVariables(StatesTemplate):
        NMINF = Float()
        PMINF = Float()
        KMINF = Float()
        NMINS = Float()
        PMINS = Float()
        KMINS = Float()
        NMINT = Float()
        PMINT = Float()
        KMINT = Float()

        NAVAIL = Float()
        PAVAIL = Float()
        KAVAIL = Float()

    class RateVariables(RatesTemplate):
        RNMINF = Float()
        RPMINF = Float()
        RKMINF = Float()
        RNMINS = Float()
        RPMINS = Float()
        RKMINS = Float()
        RNMINT = Float()
        RPMINT = Float()
        RKMINT = Float()

    def initialize(self, day, kiosk, parameters):
        self.kiosk = kiosk
        self.rates = self.RateVariables(kiosk, publish = [])
        self.params = self.Parameters(parameters)
        p = self.params
        KMINT = 0.25 * p.KMINI
        NMINT = 0.25 * p.NMINI
        PMINT = 0.25 * p.PMINI
        KMINS = 0.75 * p.KMINI
        NMINS = 0.75 * p.NMINI
        PMINS = 0.75 * p.PMINI
        NMINF = 0.
        PMINF = 0.
        KMINF = 0.

        self.states = self.StateVariables(kiosk,
                                          NMINT = NMINT,
                                          PMINT = PMINT,
                                          KMINT = KMINT,
                                          NMINS = NMINS,
                                          PMINS = PMINS,
                                          KMINS = KMINS,
                                          NMINF = NMINF,
                                          PMINF = PMINF,
                                          KMINF = KMINF,
                                          NAVAIL = NMINT,
                                          PAVAIL = PMINT,
                                          KAVAIL = KMINT,
                                          publish = [
                                              "NAVAIL",
                                              "PAVAIL",
                                              "KAVAIL"
                                          ]
                                          )
        self._connect_signal(self._on_APPLYNPK, pcse.signals.apply_npk)

    def calc_rates(self, day, drv, delt=1):
        k = self.kiosk
        p = self.params
        r = self.rates
        s = self.states
        Time = drv.DAY.timetuple().tm_yday

        RFERTN = (self._RFERTN * self._NRECOV * kg_to_g) / ha_to_m2  # g N m-2 d-1
        RFERTP = (self._RFERTP * self._PRECOV * kg_to_g) / ha_to_m2  # g N m-2 d-1
        RFERTK = (self._RFERTK * self._KRECOV * kg_to_g) / ha_to_m2  # g N m-2 d-1

        self._RFERTN = 0
        self._RFERTP = 0
        self._RFERTK = 0
        self._NRECOV = 0
        self._PRECOV = 0
        self._KRECOV = 0

        if s.NMINS / delt < p.RTNMINS:
            RNMINS = -s.NMINS / delt  # g N m-2 d-1
        else:
            RNMINS = -p.RTNMINS  # g N m-2 d-1

        if s.PMINS / delt  < p.RTPMINS:
            RPMINS = -s.PMINS / delt  # g P m-2 d-1
        else:
            RPMINS  = -p.RTPMINS  # g P m-2 d-1

        if s.KMINS / delt < p.RTKMINS:
            RKMINS = -s.KMINS / delt  # g K m-2 d-1
        else:
            RKMINS = -p.RTKMINS  # g K m-2 d-1

        # ------------ Fertilizer supply
        # Fertilizer nutrient supply
        # Pool in the soil which is not yet available for plant uptake
        # supply rate rate that becomes available for uptake
        WLIMIT = k.RFTRA / (p.K_WATER + k.RFTRA)

        RNMINF = RFERTN - p.RTNMINF * s.NMINF * WLIMIT  # g N m-2 d-1
        RPMINF = RFERTP - p.RTPMINF * s.PMINF * WLIMIT  # g P m-2 d-1
        RKMINF = RFERTK - p.RTKMINF * s.KMINF * WLIMIT  # g K m-2 d-1

        # Change in total inorganic N/P/K in soil as function of fertilizer input,
        # soil N/P/K mineralization and crop uptake.
        RNMINT = p.RTNMINF * s.NMINF * WLIMIT + (-RNMINS) - k.RNUPTR  # g N m-2 d-1
        RPMINT = p.RTPMINF * s.PMINF * WLIMIT + (-RPMINS) - k.RPUPTR  # g P m-2 d-1
        RKMINT = p.RTKMINF * s.KMINF * WLIMIT + (-RKMINS) - k.RKUPTR  # g K m-2 d-1        #

        r.RNMINT = RNMINT
        r.RPMINT = RPMINT
        r.RKMINT = RKMINT
        r.RNMINS = RNMINS
        r.RPMINS = RPMINS
        r.RKMINS = RKMINS
        r.RNMINF = RNMINF
        r.RPMINF = RPMINF
        r.RKMINF = RKMINF

    def integrate(self, day, drv, delt = 1):
        r = self.rates
        s = self.states

        NMINT = s.NMINT + delt * r.RNMINT
        PMINT = s.PMINT + delt * r.RPMINT
        KMINT = s.KMINT + delt * r.RKMINT

        s.NAVAIL = NMINT
        s.PAVAIL = PMINT
        s.KAVAIL = KMINT
        s.NMINT += delt * r.RNMINT
        s.PMINT += delt * r.RPMINT
        s.KMINT += delt * r.RKMINT
        s.NMINS += delt * r.RNMINS
        s.PMINS += delt * r.RPMINS
        s.KMINS += delt * r.RKMINS
        s.NMINF += delt * r.RNMINF
        s.PMINF += delt * r.RPMINF
        s.KMINF += delt * r.RKMINF

    def _on_APPLYNPK(self,
                    N_amount = None, N_recovery = None,
                    P_amount = None, P_recovery = None,
                    K_amount = None, K_recovery = None):
        r = self.rates
        self._RFERTN = N_amount
        self._RFERTP = P_amount
        self._RFERTK = K_amount
        self._NRECOV = N_recovery
        self._PRECOV = P_recovery
        self._KRECOV = K_recovery

