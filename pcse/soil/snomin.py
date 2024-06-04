# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Herman Berghuijs (herman.berghuijs@wur.nl) and Allard de Wit (allard.dewit@wur.nl), January 2024

import numpy as np
from .. import exceptions as exc
from pcse.decorators import prepare_rates, prepare_states
from pcse.base import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject
from pcse import signals
from ..traitlets import Float, Int, Instance, Bool
# from .soiln_profile import SoilNProfile


class SNOMIN(SimulationObject):
    """
    SNOMIN (Soil Nitrogen module for Mineral and Inorganic Nitrogen) is a layered soil nitrogen balance. A
    full mathematical description of the model is given by Berghuijs et al (2024).

    Berghuijs HNC, Silva JV, Reidsma P, De Wit AJW (2024) Expanding the WOFOST crop model
    to explore options for sustainable nitrogen management: A study for winter wheat in
    the Netherlands. European Journal of Agronomy 154 ARTN 127099. https://doi.org/10.1016/j.eja.2024.127099

    **Simulation parameters:**

    ========== ====================================================  ====================
     Name      Description                                           Unit
    ========== ====================================================  ====================
    A0SOM      Initial age of soil organic material                  y
    CNRatioBio C:N ratio of microbial biomass                        kg C kg-1 N
    FASDIS     Fraction of assimilation to dissimilation             -
    KDENIT_REF Reference first order denitrification rate constant   d-1
    KNIT_REF   Reference first order nitrification rate constant     d-1
    KSORP      Sorption coefficient ammonium (m3 water kg-1 soil)    m3 soil kg-1 soil
    MRCDIS     Michaelis Menten constant for response factor
               denitrification to soil respiration                   kg C m-2 d-1
    NO3ConcR   NO3-N concentration in rain water                     mg NO3--N L- water
    NH4ConcR   NH4-N concentration in rain water                     mg NH4+-N L-1 water
    NO3I       Initial amount of NO3-N :sup:`1`                      kg NO3--N ha-1
    NH4I       Initial amount of NH4-N :sup:`1`                      kg NH4+-N ha-1)
    WFPS_CRIT  Critical water filled pore space fraction             m3 water m-3 pore
    ========== ====================================================  ====================

    :sup:`1` This state variable is defined for each soil layer

    **State variables**

    ========== ====================================================  ==============
     Name      Description                                            Unit
    ========== ====================================================  ==============
    AGE        Appearant age of amendment (d) :sup:`1`                 d
    ORGMAT     Amount of organic matter (kg ORG ha-1) :sup:`1`         kg OM m-2
    CORG       Amount of C in organic matter (kg C ha-1) :sup:`1`      kg C m-2
    NORG       Amount of N in organic matter (kg N ha-1) :sup:`1`      kg N m-2
    NH4        Amount of NH4-N (kg N ha-1) :sup:`2`                    kg NH4-N m-2
    NO3        Amount of NO3-N (kg N ha-1) :sup:`2`                    kg NO3-N m-2
    ========== ====================================================  ==============

    | :sup:`1` This state variable is defined for each combination of soil layer and amendment
    | :sup:`2` This state variable is defined for each soil layer

    **Rate variables**

    ========== ==================================================  ====================
     Name      Description                                          Unit
    ========== ==================================================  ====================
    RAGE       Rate of change of apparent age :sup:`2`              d d-1
    RAGEAM     Initial apparent age :sup:`2`                        d d-1
    RAGEAG     Rate of ageing of amendment :sup:`2`                 d d-1
    RCORG      Rate of change of organic C :sup:`2`                 kg C m-2 d-1
    RCORGAM    Rate pf application organic C :sup:`2`               kg C m-2 d-1
    RCORGDIS   Dissimilation rate of organic C :sup:`2`             kg C m-2 d-1
    RNH4       Rate of change amount of NH4+-N :sup:`1`             kg NH4+-N m-2 d-1
    RNH4AM     Rate of NH4+-N application :sup:`1`                  kg NH4+-N m-2 d-1
    RNH4DEPOS  Rate of NH4-N deposition :sup:`1`                    kg NH4+-N m-2 d-1
    RNH4IN     Rate of NH4+-N inflow from adjacent layer :sup:`1`   kg NH4+-N m-2 d-1
    RNH4MIN    Net rate of mineralization :sup:`1`                  kg NH4+-N m-2 d-1
    RNH4NITR   Rate of nitrification :sup:`1`                       kg NH4+-N m-2 d-1
    RNH4OUT    Rate of NH4+-N outflow to adjacent layer :sup:`1`    kg NH4+-N m-2 d-1
    RNH4UP     Rate of NH4+-N root uptake :sup:`1`                  kg NH4+-N m-2 d-1
    RNO3       Rate of change amount of NO3--N :sup:`1`             kg NO3--N m-2 d-1
    RNO3AM     Rate of NO3--N application :sup:`1`                  kg NO3--N m-2 d-1
    RNO3DENITR Rate of denitrification :sup:`1`                     kg NO3--N m-2 d-1
    RNO3DEPOS  Rate of NO3--N deposition :sup:`1`                   kg NO3--N m-2 d-1
    RNO3IN     Rate of NH4+-N inflow from adjacent layer :sup:`1`   kg NO3+-N m-2 d-1
    RNO3NITR   Rate of nitrification :sup:`1`                       kg NO3--N m-2 d-1
    RNO3OUT    Rate of NO3--N outflow to adjacent layer :sup:`1`    kg NO3--N m-2 d-1
    RNO3UP     Rate of NO3--N root uptake :sup:`1`                  kg NO3--N m-2 d-1
    RNORG      Rate of change of organic N :sup:`2`                 kg N m-2 d-1
    RNORGAM    Rate pf application organic N :sup:`2`               kg N m-2 d-1
    RNORGDIS   Dissimilation rate of organic matter :sup:`2`        kg N m-2 d-1
    RORGMAT    Rate of change of organic material :sup:`2`          kg OM m-2 d-1
    RORGMATAM  Rate of application organic matter :sup:`2`          kg OM m-2 d-1
    RORGMATDIS Dissimilation rate of organic matter :sup:`2`        kg OM m-2 d-1
    ========== ==================================================  ====================

    | :sup:`1` This state variable is defined for each soil layer
    | :sup:`2` This state variable is defined for each combination of soil layer and amendment


    **Signals send or handled**

    `SNOMIN` receives the following signals:
        * APPLY_N_SNOMIN: Is received when an external input from N fertilizer
          is provided. See `_on_APPLY_N_SNOMIN()` and `signals.apply_n_snomin` for details.
    """

    # Placeholders initial values
    _ORGMATI = None
    _CORGI = None
    _NORGI = None
    _NH4I = None
    _NO3I = None

    # Placeholders
    _RNO3AM = None
    _RNH4AM = None
    _RAGEAM = None
    _RORGMATAM = None
    _RCORGAM = None
    _RNORGAM = None

    # Unit conversions
    g_to_kg = 1e-3
    cm_to_m = 1e-2
    cm2_to_ha = 1e-8
    cm3_to_m3 = 1e-6
    ha_to_m2 = 1e-4
    m2_to_ha = 1e-4
    y_to_d = 365.25

    # placeholder for soil object
    soiln_profile = None

    class StateVariables(StatesTemplate):
        AGE0   = Instance(np.ndarray) # Initial age of material (d)
        AGE    = Instance(np.ndarray) # Appearant age of material (d)
        ORGMAT = Instance(np.ndarray) # Amount of organic matter (kg ORG ha-1)
        CORG   = Instance(np.ndarray) # Amount of C in organic matter (kg C ha-1)
        NORG   = Instance(np.ndarray) #
        NH4    = Instance(np.ndarray) # Amount of NH4-N (kg N ha-1)
        NO3    = Instance(np.ndarray) # Amount of NO3-N (kg N ha-1)
        NAVAIL = Float()  # total mineral N from soil and fertiliser  kg N ha-1
        NDENITCUM = Float()
        NO3LEACHCUM = Float()
        NH4LEACHCUM = Float()
        NLOSSCUM = Float()

        RORGMATDISTT = Float()
        RORGMATAMTT = Float()
        RCORGDISTT = Float()
        RCORGAMTT = Float()
        RNORGDISTT = Float()
        RNORGAMTT = Float()

        RNO3NITRTT = Float()
        RNO3DENITRTT = Float()
        RNO3UPTT = Float()
        RNO3INTT = Float()
        RNO3OUTTT = Float()
        RNO3AMTT = Float()
        RNO3DEPOSTT = Float()

        RNH4MINTT = Float()
        RNH4NITRTT = Float()
        RNH4UPTT = Float()
        RNH4INTT = Float()
        RNH4OUTTT = Float()
        RNH4AMTT = Float()
        RNH4DEPOSTT = Float()

        ORGMATT = Float()
        CORGT = Float()
        NORGT = Float()
        RMINT = Float()
        NH4T = Float()
        NO3T = Float()

    class RateVariables(RatesTemplate):
        RAGE = Instance(np.ndarray)
        RORGMAT = Instance(np.ndarray)
        RCORG = Instance(np.ndarray)
        RNORG = Instance(np.ndarray)

        RAGEAG = Instance(np.ndarray)

        RORGMATDIS = Instance(np.ndarray)
        RCORGDIS = Instance(np.ndarray)
        RNORGDIS = Instance(np.ndarray)

        RAGEAM = Instance(np.ndarray)
        RORGMATAM = Instance(np.ndarray)
        RCORGAM = Instance(np.ndarray)
        RNORGAM = Instance(np.ndarray)

        RNH4 = Instance(np.ndarray)
        RNH4MIN = Instance(np.ndarray)
        RNH4NITR = Instance(np.ndarray)
        RNH4UP = Instance(np.ndarray)
        RNH4IN = Instance(np.ndarray)
        RNH4OUT = Instance(np.ndarray)
        RNH4AM = Instance(np.ndarray)
        RNH4DEPOS = Instance(np.ndarray)

        RNO3 = Instance(np.ndarray)
        RNO3NITR = Instance(np.ndarray)
        RNO3DENITR = Instance(np.ndarray)
        RNO3UP = Instance(np.ndarray)
        RNO3IN = Instance(np.ndarray)
        RNO3OUT = Instance(np.ndarray)
        RNO3AM = Instance(np.ndarray)
        RNO3DEPOS = Instance(np.ndarray)

        RNH4LEACHCUM = Float()
        RNO3LEACHCUM = Float()
        RNDENITCUM = Float()
        RNLOSS = Float()

    class Parameters(ParamTemplate):
        A0SOM = Float()             # Initial age of humus (y)
        CNRatioBio = Float()        # C:N ratio of microbial biomass (kg C kg-1 N)
        FASDIS = Float()            # Fraction of assimilation to dissimilation (kg ORG kg-1 ORG)
        KDENIT_REF = Float()        # Reference first order denitrification rate constant (d-1)
        KNIT_REF = Float()          # Reference first order nitrification rate constant (d-1)
        KSORP = Float()             # Sorption coefficient ammonium (m3 water kg-1 soil)
        MRCDIS = Float()            # Michaelis Menten constant for response factor denitrification to soil respiration
        NO3ConcR = Float()          # NO3-N concentration in rain water (mg N L-1)
        NH4ConcR = Float()          # NH4-N concentration in rain water (mg N L-1)
        NO3I = Instance(list)       # Initial amount of NO3-N (kg N ha-1)
        NH4I = Instance(list)       # Initial amount of NH4-N (kg N ha-1)
        WFPS_CRIT = Float()         # Critical water filled pore space fraction (m3 water m-3 pore) for denitrification

    def initialize(self, day, kiosk, parvalues):
        self.kiosk  = kiosk
        self.params = self.Parameters(parvalues)
        if "soil_profile" not in parvalues:
            msg = "Cannot find 'soil_profile' object in `parvalues`. The 'soil_profile' object should be " \
                  "instantiated by the multi-layer waterbalance before SNOMIN can run. It looks like SNOMIN " \
                  "was started before the waterbalance."
            raise exc.PCSEError(msg)
        self.soiln_profile = parvalues["soil_profile"]

        # Initialize module
        sinm = self.SoilInorganicNModel()

        # Initialize state variables
        # parvalues._soildata["soil_profile"] = self.soiln_profile
        NH4 = np.zeros(len(self.soiln_profile))
        NO3 = np.zeros_like(NH4)
        AGE =  np.zeros((1, len(self.soiln_profile)))
        AGE0 = np.zeros_like(AGE)
        ORGMAT = np.zeros_like(AGE)
        CORG =  np.zeros_like(AGE)
        NORG =  np.zeros_like(AGE)
        minip_C = self.SoilOrganicNModel.MINIP_C()
        for il, layer in enumerate(self.soiln_profile):
            NH4[il] = self.params.NH4I[il] * self.m2_to_ha
            NO3[il] = self.params.NO3I[il] * self.m2_to_ha
            AGE0[0,il] = self.params.A0SOM * self.y_to_d
            AGE[0,il] = self.params.A0SOM * self.y_to_d
            ORGMAT[0,il] = layer.RHOD_kg_per_m3 * layer.FSOMI * layer.Thickness_m
            CORG[0,il] = minip_C.calculate_organic_C(ORGMAT[0,il])
            NORG[0,il] = CORG[0, il] / layer.CNRatioSOMI

        states = dict(
            NH4 = NH4,
            NO3 = NO3,
            AGE = AGE,
            AGE0 = AGE0,
            ORGMAT = ORGMAT,
            CORG =  CORG,
            NORG =  NORG,

            # Initialize state variables to check the mass balance of organic matter
            RORGMATDISTT = 0.,
            RORGMATAMTT = 0.,

            # Initialize state variables to check the mass balance of organic C
            RCORGDISTT = 0.,
            RCORGAMTT = 0.,

            # Initialize state variables to check the mass balance of organic N
            RNORGDISTT = 0.,
            RNORGAMTT = 0.,

            # Initialize state variables to check the mass balance of NH4-N
            RNH4MINTT = 0.,
            RNH4NITRTT = 0.,
            RNH4UPTT = 0.,
            RNH4INTT = 0.,
            RNH4OUTTT = 0.,
            RNH4AMTT = 0.,
            RNH4DEPOSTT = 0.,

            # Initialize state variables to check the mass balance of NO3-N
            RNO3NITRTT = 0.,
            RNO3DENITRTT = 0.,
            RNO3UPTT = 0.,
            RNO3INTT = 0.,
            RNO3OUTTT = 0.,
            RNO3AMTT = 0.,
            RNO3DEPOSTT = 0.,

            # Initialize output variables
            CORGT = np.sum(CORG) / self.m2_to_ha,
            NORGT = np.sum(NORG) / self.m2_to_ha,
            ORGMATT = np.sum(ORGMAT) / self.m2_to_ha,
            RMINT = 0.,
            NH4T = np.sum(NH4)/ self.m2_to_ha,
            NO3T = np.sum(NO3)/ self.m2_to_ha,
            NAVAIL = 0.,
            NH4LEACHCUM = 0.,
            NO3LEACHCUM = 0.,
            NDENITCUM = 0.,
            NLOSSCUM = 0.,
        )

        self.states = self.StateVariables(kiosk, publish=["NAVAIL", "ORGMATT", "CORGT", "NORGT"], **states)
        self.rates = self.RateVariables(kiosk)

        self._RAGEAM = np.zeros_like(AGE)
        self._RORGMATAM = np.zeros_like(ORGMAT)
        self._RCORGAM = np.zeros_like(CORG)
        self._RNORGAM = np.zeros_like(NORG)
        self._RNH4AM = np.zeros_like(NH4)
        self._RNO3AM = np.zeros_like(NO3)

        self._ORGMATI = ORGMAT
        self._CORGI = CORG
        self._NORGI = NORG
        self._NH4I = NH4
        self._NO3I = NO3

        # Connect module to signal AgroManager
        self._connect_signal(self._on_APPLY_N_SNOMIN, signals.apply_n_snomin)

    @prepare_rates
    def calc_rates(self, day, drv):
        k = self.kiosk
        r = self.rates
        s = self.states
        p = self.params

        delt = 1.0

        # Initialize model components
        sonm = self.SoilOrganicNModel()
        sinm = self.SoilInorganicNModel()

        # Get external variables in the right unit
        infiltration_rate_m_per_d = self.get_infiltration_rate(k)
        flow_m_per_d = self.get_water_flow_rates(k)
        N_demand_soil = self.get_N_demand(k)
        RD_m = self.get_root_length(k)
        SM = self.get_soil_moisture_content(k)
        pF = self.get_pF(self.soiln_profile, SM)
        T = drv.TEMP

        # Get soil pH
        pH = self.get_pH(self.soiln_profile)

        # Collect ammendment rates
        r.RAGEAM = self._RAGEAM
        r.RORGMATAM = self._RORGMATAM
        r.RCORGAM = self._RCORGAM
        r.RNORGAM = self._RNORGAM
        r.RNH4AM = self._RNH4AM
        r.RNO3AM = self._RNO3AM

        # Reset placeholders for ammendment rates
        self._RAGEAM = np.zeros_like(s.AGE)
        self._RORGMATAM = np.zeros_like(r.RORGMATAM)
        self._RCORGAM = np.zeros_like(r.RCORGAM)
        self._RNORGAM = np.zeros_like(r.RNORGAM)
        self._RNH4AM = np.zeros_like(r.RNH4)
        self._RNO3AM = np.zeros_like(r.RNH4)

        # Calculate increase aparent age of each ammendment
        r.RAGEAG = sonm.calculate_apparent_age_increase_rate(s.AGE, delt, pF, pH, T)

        # Calculate dissimilation rates of each ammendment
        r.RORGMATDIS, r.RCORGDIS, r.RNORGDIS = \
            sonm.calculate_dissimilation_rates(s.AGE, p.CNRatioBio, p.FASDIS, s.NORG, s.ORGMAT, pF, pH, T)

        # Calculate rates of change apparent age, organic matter, organic C, and organic N
        r.RAGE = r.RAGEAG + r.RAGEAM
        r.RORGMAT = r.RORGMATAM - r.RORGMATDIS
        r.RCORG = r.RCORGAM - r.RCORGDIS
        r.RNORG = r.RNORGAM - r.RNORGDIS

        # Calculate N uptake rates
        r.RNH4UP, r.RNO3UP = sinm.calculate_N_uptake_rates(self.soiln_profile, delt, p.KSORP, N_demand_soil,
                                                           s.NH4, s.NO3, RD_m, SM)

        # Calculate remaining amounts of NH4-N and NO3-N after uptake and calculate chemical conversion
        NH4PRE = s.NH4 - r.RNH4UP * delt
        NO3PRE = s.NO3 - r.RNO3UP * delt
        r.RNH4MIN, r.RNH4NITR, r.RNO3NITR, r.RNO3DENITR = \
            sinm.calculate_reaction_rates(self.soiln_profile, p.KDENIT_REF, p.KNIT_REF, p.KSORP, p.MRCDIS,
                                          NH4PRE, NO3PRE, r.RCORGDIS, r.RNORGDIS, SM, T, p.WFPS_CRIT)

        # For each layer: if the net mineralization rate is negative (net immobilization) and there is not enough
        # NH4-N in the layer to be taken up by the organic matter, only the amount of N that remains after
        # denitrification is available for immobilization. RNORGDIS is updated as well to close the N balance again.
        for il, layer in enumerate(self.soiln_profile):
            if NH4PRE[il] + (r.RNH4MIN[il] - r.RNH4NITR[il]) * delt < 0:
                r.RNH4MIN[il] = NH4PRE[il] - r.RNH4NITR[il]
                RNORDIST = r.RNORGDIS[:,il].sum()
                for iam in range(0,len(r.RNORGDIS[:,il])):
                    r.RNORGDIS[iam,il] = (r.RNH4MIN[il] / RNORDIST) * r.RNORGDIS[iam,il]
            r.RNORG = r.RNORGAM - r.RNORGDIS

        # Calculate deposition rates
        r.RNH4DEPOS, r.RNO3DEPOS = sinm.calculate_deposition_rates(self.soiln_profile, infiltration_rate_m_per_d,
                                                                   s.NH4, p.NH4ConcR, s.NO3, p.NO3ConcR)

        # Calculate remaining amounts of NH4-N and NO3-N after uptake and reaction and calculate inorganic
        # N flow rates between layers
        NH4PRE2 = NH4PRE + (r.RNH4AM + r.RNH4MIN + r.RNH4DEPOS - r.RNH4NITR) * delt
        NO3PRE2 = NO3PRE + (r.RNO3AM + r.RNO3NITR + r.RNO3DEPOS  - r.RNO3DENITR) * delt
        r.RNH4IN, r.RNH4OUT, r.RNO3IN, r.RNO3OUT = \
            sinm.calculate_flow_rates(self.soiln_profile, flow_m_per_d, p.KSORP, NH4PRE2, NO3PRE2, SM)

        # Calculate rates of change NH4-N and NO3-N
        r.RNH4 = r.RNH4AM + r.RNH4MIN + r.RNH4DEPOS - r.RNH4NITR - r.RNH4UP + r.RNH4IN - r.RNH4OUT
        r.RNO3 = r.RNO3AM + r.RNO3NITR + r.RNO3DEPOS - r.RNO3DENITR - r.RNO3UP + r.RNO3IN - r.RNO3OUT

        # Calculate output rate variables for N loss
        r.RNH4LEACHCUM =  (1/self.m2_to_ha) * r.RNH4OUT[-1]
        r.RNO3LEACHCUM =  (1/self.m2_to_ha) * r.RNO3OUT[-1]
        r.RNDENITCUM = (1/self.m2_to_ha) * r.RNO3DENITR.sum()

    @prepare_states
    def integrate(self, day, delt=1.0):
        k = self.kiosk
        p = self.params
        r = self.rates
        s = self.states

        # Initialize soil module
        sinm = self.SoilInorganicNModel()

        # Calculate apparent ages and amounts of organic matter, organic C, and organic N in next time step
        AGE = s.AGE + r.RAGE * delt
        ORGMAT = s.ORGMAT + r.RORGMAT * delt
        CORG = s.CORG + r.RCORG * delt
        NORG= s.NORG + r.RNORG * delt

        # Calculate amounts of NH4-N and NO3-N in next time step
        NH4 = s.NH4 + r.RNH4 * delt
        NO3 = s.NO3 + r.RNO3 * delt

        # Update state variables
        s.AGE = AGE
        s.ORGMAT = ORGMAT
        s.CORG = CORG
        s.NORG = NORG
        s.NH4 = NH4
        s.NO3 = NO3

        # Get external state variable
        RD_m = self.get_root_length(k)
        SM = self.get_soil_moisture_content(k)

        # Calculate the amount of N available for root uptake in the next time step
        s.NAVAIL = sinm.calculate_NAVAIL(self.soiln_profile, p.KSORP, s.NH4, s.NO3, RD_m, SM) / self.m2_to_ha
        self.check_mass_balances(day, delt)

        # Set output variables
        s.ORGMATT = np.sum(s.ORGMAT)  * (1/self.m2_to_ha)
        s.CORGT = np.sum(s.CORG)  * (1/self.m2_to_ha)
        s.NORGT = np.sum(s.NORG)  * (1/self.m2_to_ha)
        s.RMINT += np.sum(r.RNORGDIS) * (1/self.m2_to_ha)
        s.NH4T = np.sum(s.NH4) * (1/self.m2_to_ha)
        s.NO3T = np.sum(s.NO3) * (1/self.m2_to_ha)
        NH4LEACHCUM = s.NH4LEACHCUM + r.RNH4LEACHCUM * delt
        NO3LEACHCUM = s.NO3LEACHCUM + r.RNO3LEACHCUM * delt
        NDENITCUM = s.NDENITCUM + r.RNDENITCUM * delt
        NLOSSCUM =  NH4LEACHCUM + NO3LEACHCUM + NDENITCUM
        s.NH4LEACHCUM = NH4LEACHCUM
        s.NO3LEACHCUM = NO3LEACHCUM
        s.NDENITCUM = NDENITCUM
        s.NLOSSCUM = NLOSSCUM

    def _on_APPLY_N_SNOMIN(self, amount=None, application_depth = None, cnratio=None, f_orgmat=None,
                           f_NH4N = None, f_NO3N = None, initial_age =None):
        """This function calculates the application rates of organic matter, organic C, organic N, NH4-N, NO3-N
        and the initial apparent age of the applied material at the date of application.

        For each amendment, the following variables need to be provided in the AgroManagement
        file of the simulation:

        **Amendment properties**
        ================== ======================================================    =========================
         Name              Description                                               Unit
        ================== ======================================================    =========================
        amount             Amount of material in amendment                           kg material ha-1
        application_depth  Depth over which the amendment is applied in the soil     cm
        cnratio            C:N ratio of organic matter in material                   kg C kg-1 N
        initial_age        Initial apparent age of organic matter in material        y
        f_NH4N             Fraction of NH4+-N in material                            kg NH4+-N kg-1 material
        f_NO3N             Fraction of NO3--N in material                            kg NO3--N kg-1 material
        f_orgmat           Fraction of organic matter in amendment                   kg OM kg-1 material
        ================== ======================================================    =========================
        """

        r = self.rates
        s = self.states
        delt = 1.

        # Create model components
        sinm = self.SoilInorganicNModel()
        sonm = self.SoilOrganicNModel()

        # Initialize amendment rates
        RAGE_am = np.zeros((1, len(self.soiln_profile)))
        AGE0_am = np.zeros_like(RAGE_am)
        RORGMAT_am = np.zeros_like(RAGE_am)
        RCORG_am = np.zeros_like(RAGE_am)
        RNORG_am = np.zeros_like(RAGE_am)
        RNH4_am = np.zeros_like(s.NH4)
        RNO3_am = np.zeros_like(s.NO3)

        # Prevents that a part of the N is not applied if the application depth is less than thickness of the upper layer
        if application_depth < self.soiln_profile[0].Thickness:
            application_depth = self.soiln_profile[0].Thickness

        AGE0_am[0, :] = initial_age * self.y_to_d
        RAGE_am[0, :] = initial_age * self.y_to_d
        RNH4_am, RNO3_am = np.array(sinm.calculate_N_application_amounts(self.soiln_profile, amount, application_depth,
                                                                          f_NH4N, f_NO3N)) * self.m2_to_ha
        RORGMAT_am, RCORG_am, RNORG_am = np.array(sonm.calculate_application_rates(self.soiln_profile, amount,
                                                                                   application_depth, cnratio, f_orgmat)
                                                  ) * self.m2_to_ha

        # Add a new column to the state variables for organic amendments to add the new amendment.
        s.AGE0 = np.concatenate((s.AGE0, AGE0_am), axis = 0)
        s.ORGMAT = np.concatenate((s.ORGMAT, np.zeros((1, len(self.soiln_profile)))), axis = 0)
        s.CORG = np.concatenate((s.CORG, np.zeros((1, len(self.soiln_profile)))), axis = 0)
        s.NORG = np.concatenate((s.NORG, np.zeros((1, len(self.soiln_profile)))), axis = 0)
        s.AGE = np.concatenate((s.AGE, np.zeros((1, len(self.soiln_profile)))), axis = 0)

        # Store the amendment rates
        self._RAGEAM = np.concatenate((self._RAGEAM, RAGE_am), axis = 0)
        self._RORGMATAM = np.concatenate((self._RORGMATAM, RORGMAT_am), axis = 0)
        self._RCORGAM = np.concatenate((self._RCORGAM, RCORG_am), axis = 0)
        self._RNORGAM = np.concatenate(( self._RNORGAM, RNORG_am), axis = 0)
        self._RNH4AM = RNH4_am
        self._RNO3AM = RNO3_am

    def get_infiltration_rate(self, k):
        infiltration_rate_m_per_d = k.RIN * self.cm_to_m
        return infiltration_rate_m_per_d

    def get_pF(self, soiln_profile, SM):
        pF = np.zeros_like(SM)
        for il, layer in enumerate(soiln_profile):
            pF[il] = layer.PFfromSM(SM[il])
        return pF

    def get_pH(self, soiln_profile):
        pH = np.zeros(len(soiln_profile))
        for il, layer in enumerate(soiln_profile):
            pH[il] = layer.Soil_pH
        return pH

    def get_soil_moisture_content(self, k):
        SM = k.SM
        return SM

    def get_water_flow_rates(self, k):
        flow_m_per_d = k.Flow * self.cm_to_m
        return flow_m_per_d

    def get_N_demand(self, k):
        if "RNuptake" in k:
            N_demand_soil = k.RNuptake * self.m2_to_ha
        else:
            N_demand_soil = 0.
        return N_demand_soil

    def get_root_length(self, k):
        if "RD" in k:
            RD_m = k.RD * self.cm_to_m
        else:
            RD_m = 0.
        return RD_m

    def check_mass_balances(self, day, delt):
        s = self.states
        r = self.rates

        s.RORGMATAMTT += delt * r.RORGMATAM.sum()
        s.RORGMATDISTT += delt * r.RORGMATDIS.sum()
        s.RCORGAMTT += delt * r.RCORGAM.sum()
        s.RCORGDISTT += delt * r.RCORGDIS.sum()
        s.RNORGAMTT += delt * r.RNORGAM.sum()
        s.RNORGDISTT += delt * r.RNORGDIS.sum()

        s.RNH4MINTT += delt * r.RNH4MIN.sum()
        s.RNH4NITRTT += delt * r.RNH4NITR.sum()
        s.RNH4UPTT += delt * r.RNH4UP.sum()
        s.RNH4INTT += delt * r.RNH4IN.sum()
        s.RNH4OUTTT += delt * r.RNH4OUT.sum()
        s.RNH4AMTT += delt * r.RNH4AM.sum()
        s.RNH4DEPOSTT += delt * r.RNH4DEPOS.sum()

        s.RNO3NITRTT += delt * r.RNO3NITR.sum()
        s.RNO3DENITRTT += delt * r.RNO3DENITR.sum()
        s.RNO3UPTT += delt * r.RNO3UP.sum()
        s.RNO3INTT += delt * r.RNO3IN.sum()
        s.RNO3OUTTT += delt * r.RNO3OUT.sum()
        s.RNO3AMTT += delt * r.RNO3AM.sum()
        s.RNO3DEPOSTT += delt * r.RNO3DEPOS.sum()

        ORGMATBAL = self._ORGMATI.sum() - s.ORGMAT.sum() + s.RORGMATAMTT - s.RORGMATDISTT
        if abs(ORGMATBAL) > 0.0001:
            msg = "Organic matter mass balance is not closing on %s with checksum: %f" % (day, ORGMATBAL)
            raise exc.SoilOrganicMatterBalanceError(msg)

        CORGBAL = self._CORGI.sum() - s.CORG.sum() + s.RCORGAMTT - s.RCORGDISTT
        if abs(CORGBAL) > 0.0001:
            msg = "Organic carbon mass balance is not closing on %s with checksum: %f" % (day, CORGBAL)
            raise exc.SoilOrganicCarbonBalanceError(msg)

        NORGBAL = self._NORGI.sum() - s.NORG.sum() + s.RNORGAMTT - s.RNORGDISTT
        if abs(NORGBAL) > 0.0001:
            msg = "Organic carbon mass balance is not closing on %s with checksum: %f" % (day, NORGBAL)
            raise exc.SoilOrganicNitrogenBalanceError(msg)

        NH4BAL = self._NH4I.sum() - s.NH4.sum() + s.RNH4AMTT + s.RNH4INTT + s.RNH4MINTT + s.RNH4DEPOSTT - s.RNH4NITRTT - s.RNH4OUTTT - s.RNH4UPTT
        if abs(NH4BAL) > 0.0001:
            msg = "NH4-N mass balance is not closing on %s with checksum: %f" % (day, NH4BAL)
            raise exc.SoilAmmoniumBalanceError(msg)

        NO3BAL = self._NO3I.sum() - s.NO3.sum() + s.RNO3AMTT + s.RNO3NITRTT + s.RNO3INTT + s.RNO3DEPOSTT - s.RNO3DENITRTT - s.RNO3OUTTT - s.RNO3UPTT
        if abs(NO3BAL) > 0.0001:
            msg = "NO3-N mass balance is not closing on %s with checksum: %f" % (day, NO3BAL)
            raise exc.SoilNitrateBalanceError(msg)

    class SoilInorganicNModel:
        def calculate_N_application_amounts(self, soiln_profile, amount, application_depth, f_NH4N, f_NO3N):
            samm = self.SoilAmmoniumNModel()
            sni = self.SoilNNitrateModel()
            RNH4_am = np.zeros(len(soiln_profile))
            RNO3_am = np.zeros_like(RNH4_am)
            zmin = 0
            for il, layer in enumerate(soiln_profile):
                zmax = zmin + soiln_profile[il].Thickness
                RNH4_am[il] = samm.calculate_NH4_application_amount(amount, application_depth, f_NH4N, layer.Thickness, zmax, zmin)
                RNO3_am[il] = sni.calculate_NO3_application_amount(amount, application_depth, f_NO3N, layer.Thickness, zmax, zmin)
                zmin = zmax
            return RNH4_am, RNO3_am

        def calculate_flow_rates(self, soiln_profile, flow_m_per_d, KSORP, NH4, NO3, SM):
            samm = self.SoilAmmoniumNModel()
            sni = self.SoilNNitrateModel()
            RNH4IN, RNH4OUT = samm.calculate_NH4_flow_rates(soiln_profile, flow_m_per_d, KSORP, NH4, SM)
            RNO3IN, RNO3OUT = sni.calculate_NO3_flow_rates(soiln_profile, flow_m_per_d, NO3, SM)
            return RNH4IN, RNH4OUT, RNO3IN, RNO3OUT

        def calculate_deposition_rates(self,soiln_profile,infiltration_rate_m_per_d, NH4, NH4ConcR, NO3, NO3ConcR):
            samm = self.SoilAmmoniumNModel()
            sni = self.SoilNNitrateModel()
            RNH4DEPOS = samm.calculate_NH4_deposition_rates(soiln_profile, infiltration_rate_m_per_d, NH4, NH4ConcR)
            RNO3DEPOS = sni.calculate_NO3_deposition_rates(soiln_profile, infiltration_rate_m_per_d, NO3, NO3ConcR)
            return RNH4DEPOS, RNO3DEPOS

        def calculate_reaction_rates(self, soiln_profile, KDENIT_REF, KNIT_REF, KSORP, MRCDIS, NH4, NO3, RCORGDIS, RNORGDIS, SM, T, WFPS_CRIT):
            samm = self.SoilAmmoniumNModel()
            sni = self.SoilNNitrateModel()
            RNH4MIN, RNH4NITR = samm.calculate_NH4_reaction_rates(soiln_profile, KNIT_REF, KSORP, NH4, RNORGDIS, SM, T)
            RNO3NITR, RNO3DENITR = sni.calculate_NO3_reaction_rates(soiln_profile, KDENIT_REF, MRCDIS, NO3, RCORGDIS, RNH4NITR, SM, T, WFPS_CRIT)
            return RNH4MIN, RNH4NITR, RNO3NITR, RNO3DENITR

        def calculate_NAVAIL(self, soiln_profile, KSORP, NH4, NO3, RD_m, SM):
            samm = self.SoilAmmoniumNModel()
            sni = self.SoilNNitrateModel()
            zmin = 0.
            NAVAIL = 0.
            for il, layer in enumerate(soiln_profile):
                zmax = zmin + layer.Thickness_m
                NH4_avail_layer = samm.calculate_available_NH4(KSORP, NH4[il], RD_m, layer.RHOD_kg_per_m3, SM[il], zmax, zmin)
                NO3_avail_layer = sni.calculate_available_NO3(NO3[il], RD_m, SM[il], zmax, zmin)
                NAVAIL += (NH4_avail_layer + NO3_avail_layer)
                zmin = zmax
            return NAVAIL

        def calculate_N_uptake_rates(self, soiln_profile, delt, KSORP, N_demand_soil, NH4, NO3, RD_m, SM):
            RNH4UP = np.zeros_like(NH4)
            RNO3UP = np.zeros_like(NO3)
            samm = self.SoilAmmoniumNModel()
            sni = self.SoilNNitrateModel()
            zmin = 0.
            for il, layer in enumerate(soiln_profile):
                zmax = zmin + layer.Thickness_m
                RNH4UP[il] = samm.calculate_NH4_plant_uptake_rate(KSORP, N_demand_soil, NH4[il], RD_m, layer.RHOD_kg_per_m3, SM[il], zmax, zmin)
                N_demand_soil -= RNH4UP[il] * delt
                RNO3UP[il] = sni.calculate_NO3_plant_uptake_rate(N_demand_soil, NO3[il], RD_m, SM[il], zmax, zmin)
                N_demand_soil -= RNO3UP[il] * delt
                zmin = zmax
            return RNH4UP, RNO3UP

        class SoilAmmoniumNModel:
            def calculate_NH4_deposition_rates(self, soiln_profile, infiltration_rate_m_per_d, NH4, NH4ConcR):
                RNH4DEPOS = np.zeros_like(NH4)
                mg_to_kg = 1e-6
                L_to_m3 = 1e-3
                for il, layer in enumerate(soiln_profile):
                    if il == 0:
                        RNH4DEPOS[il] = (mg_to_kg / L_to_m3) * NH4ConcR * infiltration_rate_m_per_d
                    else:
                        RNH4DEPOS[il] = 0.
                return RNH4DEPOS

            def calculate_NH4_reaction_rates(self, soiln_profile, KNITREF, KSORP, NH4, RNORGDIS, SM, T):
                RNH4MIN = np.zeros_like(NH4)
                RNH4NITR = np.zeros_like(NH4)
                for il, layer in enumerate(soiln_profile):
                    RNMIN_kg_per_m2 = RNORGDIS[:,il].sum()
                    RNH4MIN[il] = self.calculate_mineralization_rate(RNMIN_kg_per_m2)
                    RNH4NITR[il] = self.calculate_nitrification_rate(KNITREF, KSORP, layer.Thickness_m, NH4[il], layer.RHOD_kg_per_m3, SM[il], layer.SM0, T)
                return RNH4MIN, RNH4NITR

            def calculate_NH4_flow_rates(self, soiln_profile, flow_m_per_d, KSORP, NH4, SM):
                cNH4Kwel = 0.

                RNH4IN = np.zeros_like(NH4)
                RNH4OUT = np.zeros_like(NH4)
                RHOD = np.zeros_like(NH4)
                cNH4 = np.zeros_like(NH4)
                dz = np.zeros_like(NH4)

                # Downward flow
                for il, layer in enumerate(soiln_profile):
                    RHOD[il] = layer.RHOD_kg_per_m3
                    dz[il] = layer.Thickness_m
                    cNH4[il] = self.calculate_NH4_concentration(KSORP, dz[il], NH4[il], RHOD[il], SM[il])

                for il in range(0,len(soiln_profile)):
                    if flow_m_per_d[il] >= 0.:
                        if il == 0:
                            RNH4IN[il] += 0.
                        else:
                            RNH4IN[il] += flow_m_per_d[il] * cNH4[il - 1]
                            RNH4OUT[il-1] += flow_m_per_d[il] * cNH4[il - 1]
                if flow_m_per_d[len(NH4) - 1] >= 0.:
                    RNH4OUT[len(NH4) - 1] += flow_m_per_d[len(NH4)] * cNH4[len(NH4) - 1]

                ## Upward flow
                for il in reversed(range(0,len(soiln_profile))):
                    if flow_m_per_d[il + 1] < 0.:
                        if il == len(NH4) - 1:
                            RNH4IN[il] += - flow_m_per_d[il + 1] * cNH4Kwel
                        else:
                            RNH4IN[il] += - flow_m_per_d[il + 1] * cNH4[il + 1]
                            RNH4OUT[il + 1] += - flow_m_per_d[il + 1] * cNH4[il + 1]
                if flow_m_per_d[0] < 0.:
                    RNH4OUT[0] += - flow_m_per_d[0] * cNH4[0]
                else:
                    RNH4OUT[0] += 0.
                return RNH4IN, RNH4OUT

            def calculate_NH4_application_amount(self, amount, application_depth, f_NH4N, layer_thickness, zmax, zmin):
                if application_depth > zmax:
                    NH4_am = (layer_thickness / application_depth) * f_NH4N *  amount
                elif zmin <= application_depth <= zmax:
                    NH4_am = ((application_depth - zmin) / application_depth) * f_NH4N  * amount
                else:
                    NH4_am = 0.
                return NH4_am

            def calculate_NH4_concentration(self, KSORP, layer_thickness, NH4, RHOD_kg_per_m3, SM):
                cNH4 = (1 / ( KSORP * RHOD_kg_per_m3 + SM)) * NH4 / layer_thickness
                return cNH4

            def calculate_available_NH4(self, KSORP, NH4, RD, RHOD_kg_per_m3, SM, zmax, zmin):
                layer_thickness = zmax - zmin
                cNH4 = self.calculate_NH4_concentration(KSORP, layer_thickness, NH4, RHOD_kg_per_m3, SM)

                if RD <= zmin:
                    NH4_avail = 0.
                elif RD > zmax:
                    NH4_avail = (SM / ( KSORP * RHOD_kg_per_m3 + SM)) * NH4
                else:
                    NH4_avail = ((RD - zmin)/ layer_thickness) * (SM / ( KSORP * RHOD_kg_per_m3 + SM)) * NH4
                return NH4_avail

            def calculate_mineralization_rate(self, rNMINs_layer):
                RNH4MIN = rNMINs_layer.sum()
                return RNH4MIN

            def calculate_nitrification_rate(self, KNIT_REF, KSORP, layer_thickness, NH4, RHOD_kg_per_m3, SM, SM0, T):
                cNH4 = self.calculate_NH4_concentration(KSORP, layer_thickness, NH4, RHOD_kg_per_m3, SM)
                fWNIT = self.calculate_soil_moisture_response_nitrification_rate_constant(SM, SM0)
                fT = self.calculate_temperature_response_nitrification_rate_constant(T)
                RNH4NIT = fWNIT * fT *  KNIT_REF * SM * cNH4 * layer_thickness
                return RNH4NIT

            def calculate_NH4_plant_uptake_rate(self, KSORP, N_demand_soil, NH4, RD_m, RHOD_kg_per_m3, SM, zmax, zmin):
                NH4_av = self.calculate_available_NH4(KSORP, NH4, RD_m, RHOD_kg_per_m3, SM, zmax, zmin)
                RNH4UP = min(N_demand_soil, NH4_av)
                return RNH4UP

            def calculate_soil_moisture_response_nitrification_rate_constant(self, SM, SM0):
                WFPS = SM / SM0
                fWNIT = 0.9 / (1. + np.exp(-15 *(WFPS - 0.45))) + 0.1 - 1/(1+np.exp(-50. * (WFPS - 0.95)))
                return fWNIT

            def calculate_temperature_response_nitrification_rate_constant(self, T):
                fT = 1/(1+np.exp(-0.26*(T-17.)))-1/(1+np.exp(-0.77*(T-41.9)))
                return fT

        class SoilNNitrateModel:
            def calculate_NO3_deposition_rates(self, soiln_profile, infiltration_rate_m_per_d, NO3, NO3ConcR):
                mg_to_kg = 1e-6
                L_to_m3 = 1e-3
                RNO3DEPOS = np.zeros_like(NO3)
                for il, layer in enumerate(soiln_profile):
                    if il == 0:
                        RNO3DEPOS[il] = (mg_to_kg / L_to_m3) * NO3ConcR * infiltration_rate_m_per_d
                    else:
                        RNO3DEPOS[il] = 0.
                return RNO3DEPOS

            def calculate_NO3_flow_rates(self, soiln_profile, flow_m_per_d, NO3, SM):
                cNO3Kwel = 0.

                RNO3IN = np.zeros_like(NO3)
                RNO3OUT = np.zeros_like(NO3)
                cNO3 = np.zeros_like(NO3)
                dz = np.zeros_like(NO3)

                # Downward flow
                for il, layer in enumerate(soiln_profile):
                    dz[il] = layer.Thickness_m
                    cNO3[il] =  self.calculate_NO3_concentration(dz[il], NO3[il], SM[il])

                for il in range(0,len(soiln_profile)):
                    if flow_m_per_d[il] >= 0.:
                        if il == 0:
                            RNO3IN[il] += 0.
                        else:
                            RNO3IN[il] += flow_m_per_d[il] * cNO3[il - 1]
                            RNO3OUT[il-1] += flow_m_per_d[il] * cNO3[il - 1]
                if flow_m_per_d[len(NO3) - 1] >= 0.:
                    RNO3OUT[len(NO3) - 1] += flow_m_per_d[len(NO3)] * cNO3[len(NO3) - 1]
                else:
                    RNO3OUT[len(NO3) - 1] += 0

                # Upward flow
                for il in reversed(range(0,len(soiln_profile))):
                    if flow_m_per_d[il + 1] < 0.:
                        if il == len(NO3) - 1:
                            RNO3IN[il] += - flow_m_per_d[il + 1] * cNO3Kwel
                        else:
                            RNO3IN[il] += - flow_m_per_d[il + 1] * cNO3[il + 1]
                            RNO3OUT[il + 1] += - flow_m_per_d[il + 1] * cNO3[il + 1]
                if flow_m_per_d[0] < 0.:
                    RNO3OUT[0] += - flow_m_per_d[0] * cNO3[0]
                else:
                    RNO3OUT[0] += 0.
                return RNO3IN, RNO3OUT

            def calculate_NO3_reaction_rates(self, soiln_profile, KDENIT_REF, MRCDIS, NO3, RCORGDIS,
                                             RNH4NITR, SM, T, WFPS_CRIT):
                RNO3NITR = np.zeros_like(NO3)
                RNO3DENITR = np.zeros_like(NO3)
                for il, layer in enumerate(soiln_profile):
                    RNO3NITR[il] = RNH4NITR[il]
                    RNO3DENITR[il] = \
                        self.calculate_denitrification_rate(layer.Thickness_m, NO3[il], KDENIT_REF, MRCDIS,
                                                            RCORGDIS[:,il].sum(), SM[il], layer.SM0, T, WFPS_CRIT)
                return RNO3NITR, RNO3DENITR

            def calculate_NO3_application_amount(self, amount, application_depth, f_NO3N, layer_thickness, zmax, zmin):
                if application_depth > zmax:
                    NO3_am = (layer_thickness / application_depth) * f_NO3N *  amount
                elif zmin <= application_depth <= zmax:
                    NO3_am = ((application_depth - zmin) / application_depth) * f_NO3N  * amount
                else:
                    NO3_am = 0.
                return NO3_am

            def calculate_NO3_concentration(self, layer_thickness, NO3, SM):
                cNO3 = NO3 / (layer_thickness * SM)
                return cNO3

            def calculate_available_NO3(self, NO3, RD, SM, zmax, zmin):
                layer_thickness = zmax - zmin
                if RD <= zmin:
                    NO3_avail_layer = 0.
                elif RD > zmax:
                    NO3_avail_layer = NO3
                else:
                    NO3_avail_layer = ((RD - zmin)/ layer_thickness) * NO3
                return NO3_avail_layer

            def calculate_denitrification_rate(self, layer_thickness, NO3, KDENIT_REF, MRCDIS, RCORGT_kg_per_m2,
                                               SM, SM0, T, WFPS_CRIT):
                cNO3 = self.calculate_NO3_concentration(layer_thickness, NO3, SM)
                fR = self.calculate_soil_respiration_response_denitrifiation_rate_constant(RCORGT_kg_per_m2, MRCDIS)
                fW = self.calculate_soil_moisture_response_denitrification_rate_constant(SM, SM0, WFPS_CRIT)
                fT = self.calculate_temperature_response_denitrification_rate_constant(T)
                RNO3DENIT = fW * fT * fR * KDENIT_REF * SM * cNO3 * layer_thickness
                return RNO3DENIT

            def calculate_NO3_plant_uptake_rate(self, N_demand_soil, NO3, RD_m, SM, zmax, zmin):
                NO3_av = self.calculate_available_NO3(NO3, RD_m, SM, zmax, zmin)
                RNO3UP = min(N_demand_soil, NO3_av)
                return RNO3UP

            def calculate_soil_moisture_response_denitrification_rate_constant(self, SM, SM0, WFPS_CRIT):
                WFPS = SM / SM0
                if WFPS < WFPS_CRIT:
                    fW = 0.
                else:
                    fW = np.power((WFPS - WFPS_CRIT)/(1 - WFPS_CRIT),2)
                return fW

            def calculate_soil_respiration_response_denitrifiation_rate_constant(self, RCORGT, MRCDIS):
                fR = RCORGT / (MRCDIS + RCORGT)
                return fR

            def calculate_temperature_response_denitrification_rate_constant(self, T):
                fT = 1/(1+np.exp(-0.26*(T-17)))-1/(1+np.exp(-0.77*(T-41.9)))
                return fT

    class SoilOrganicNModel:
        def calculate_apparent_age_increase_rate(self, AGE, delt, pF, pH, T):
            RAGEAG = np.zeros_like(AGE)
            janssen = self.Janssen()
            for am in range(0, AGE.shape[0]):
                for il in range(0, AGE.shape[1]):
                    RAGEAG[am,il] = janssen.calculate_increase_apparent_age_rate(delt, pF[il], pH[il], T)
            return RAGEAG

        def calculate_application_rates(self, soiln_profile, amount, application_depth, cnratio, f_orgmat):
            RORGMAT_am = np.zeros((1, len(soiln_profile)))
            RCORG_am = np.zeros_like(RORGMAT_am)
            RNORG_am = np.zeros_like(RORGMAT_am)

            zmin = 0.
            for il, layer in enumerate(soiln_profile):
                zmax = zmin + soiln_profile[il].Thickness
                RORGMAT_am[0, il] = \
                    self.calculate_organic_material_application_amount(amount, application_depth, f_orgmat,
                                                                       layer.Thickness, zmax, zmin)
                RCORG_am[0, il] = \
                    self.calculate_organic_carbon_application_amount(amount, application_depth, f_orgmat,
                                                                     layer.Thickness, zmax, zmin)
                RNORG_am[0, il] = \
                    self.calculate_organic_nitrogen_application_amount(amount, application_depth, cnratio, f_orgmat,
                                                                       layer.Thickness, zmax, zmin)
                zmin = zmax

            return RORGMAT_am, RCORG_am, RNORG_am

        def calculate_dissimilation_rates(self, AGE, CNRatioBio, FASDIS, NORG, ORGMAT, pF, pH, T):
            RORGMATDIS = np.zeros_like(AGE)
            RCORGDIS = np.zeros_like(AGE)
            RNORGDIS = np.zeros_like(AGE)
            janssen = self.Janssen()
            minip_c = self.MINIP_C()
            minip_n = self.MINIP_N()

            for am in range(0, AGE.shape[0]):
                for il in range(0, AGE.shape[1]):
                    if ORGMAT[am, il] > 0:
                        RORGMATDIS[am,il] = \
                            janssen.calculate_dissimilation_rate_OM_T(ORGMAT[am,il], AGE[am,il], pF[il], pH[il], T)
                        RCORGDIS[am,il] = \
                            minip_c.calculate_dissimilation_rate_C(janssen, ORGMAT[am,il], AGE[am,il], pF[il],
                                                                   pH[il], T)
                        RNORGDIS[am,il] = \
                            minip_n.calculate_dissimilation_rate_N(janssen, minip_c, ORGMAT[am,il], NORG[am,il],
                                                                   FASDIS, CNRatioBio, AGE[am,il], pF[il], pH[il], T)
                    else:
                        RORGMATDIS[am,il] = 0.
                        RCORGDIS[am,il] = 0.
                        RNORGDIS[am,il] = 0.
            return RORGMATDIS, RCORGDIS, RNORGDIS

        def calculate_organic_carbon_application_amount(self, amount, application_depth, f_orgmat, layer_thickness, zmax, zmin):
            minip_C = self.MINIP_C()
            ORGMAT_am = \
                self.calculate_organic_material_application_amount(amount, application_depth, f_orgmat,
                                                                   layer_thickness, zmax, zmin)
            CORG_am = minip_C.calculate_organic_C(ORGMAT_am)
            return CORG_am

        def calculate_organic_nitrogen_application_amount(self, amount, application_depth, cnratio, f_orgmat, layer_thickness, zmax, zmin):
            minip_C = self.MINIP_C()
            ORGMAT_am = self.calculate_organic_material_application_amount(amount, application_depth, f_orgmat, layer_thickness, zmax, zmin)
            CORG_am = minip_C.calculate_organic_C(ORGMAT_am)
            if cnratio == 0:
                NORG_am = 0.
            else:
                NORG_am = CORG_am / cnratio
            return NORG_am

        def calculate_organic_material_application_amount(self, amount, application_depth, f_orgmat, layer_thickness, zmax, zmin):
            if application_depth > zmax:
                ORGMAT_am = (layer_thickness / application_depth) * f_orgmat * amount
            elif zmin <= application_depth <= zmax:
                ORGMAT_am = ((application_depth - zmin) / application_depth) * f_orgmat * amount
            else:
                ORGMAT_am = 0
            return ORGMAT_am

        class Janssen:
            m = 1.6
            b = 2.82
            y_to_d = 365.25

            def calculate_increase_apparent_age_rate(self, dt, pF, pH, T):
                f_pH = self.calculate_pH_response_dissimilation_rate(pH)
                f_T = self.calculate_temperature_response_dissimilation_rate_Yang(T)
                f_SM = self.calculate_soil_moisture_response_dissimilation_rate(pF)
                dA = f_pH * f_T * f_SM * dt
                return dA

            def calculate_relative_dissimilation_rate_OM_T(self, t, pF, pH, T):
                m = self.m
                b = self.b
                f_pH = self.calculate_pH_response_dissimilation_rate(pH)
                f_T = self.calculate_temperature_response_dissimilation_rate_Yang(T)
                f_SM = self.calculate_soil_moisture_response_dissimilation_rate(pF)
                k = f_pH * f_T * f_SM * b *  pow(t/self.y_to_d, -m) / self.y_to_d
                return k

            def calculate_dissimilation_rate_OM_T(self, OM, t, pF, pH, T):
                k = self.calculate_relative_dissimilation_rate_OM_T(t, pF, pH, T)
                rate = k * OM
                return rate

            def calculate_soil_moisture_response_dissimilation_rate(self, pF):
                if pF < 2.7:
                    f_SM = 1.0
                elif pF < 4.2:
                    f_SM = 1.0 * (4.2 - pF) / (4.2 - 2.7)
                else:
                    f_SM = 0.
                return f_SM

            def calculate_pH_response_dissimilation_rate(self, pH):
                f_pH = 1 / (1 + np.exp(-1.5 * (pH - 4)))
                return f_pH

            def calculate_temperature_response_dissimilation_rate(self, T):
                f_T = pow(2, (T-9)/9)
                return f_T

            def calculate_temperature_response_dissimilation_rate_Yang(self, T):
                if T < -1:
                    f_T = 0.
                elif T < 9.:
                    f_T = 0.09 * (T + 1)
                elif T < 27.:
                    f_T = 0.88 * pow(2, (T-9)/9)
                else:
                    f_T = 3.5
                return f_T

        class MINIP_C:
            OM_to_C = 0.58
            y_to_d = 365.

            def calculate_assimilation_rate(self, janssen, OM, f_ass_dis, t, pF, pH, T):
                r_disc = self.calculate_dissimilation_rate_C(janssen, OM, t, pF, pH, T)
                r_ass = r_disc * f_ass_dis
                return r_ass

            def calculate_dissimilation_rate_C(self, janssen, OM, t, pF, pH, T):
                k = janssen.calculate_relative_dissimilation_rate_OM_T(t, pF, pH, T)
                Corg = self.calculate_organic_C(OM)
                rate = k * Corg
                return rate

            def calculate_total_conversion_rate_C(self, janssen, OM, f_ass_dis, t, pF, pH, T):
                r_dis_C = self.calculate_dissimilation_rate_C(janssen, OM, t, pF, pH, T)
                r_ass_C = self.calculate_assimilation_rate(janssen, OM, f_ass_dis, t, pF, pH, T)
                r_conv_C = r_dis_C + r_ass_C
                return r_conv_C

            def calculate_organic_C(self, OM):
                Corg = OM * self.OM_to_C
                return Corg

        class MINIP_N:

            def calculate_total_conversion_rate_N(self, janssen, minip_c, OM, Norg, f_ass_dis, t, pF, pH, T):
                r_conv_C = minip_c.calculate_total_conversion_rate_C(janssen, OM, f_ass_dis, t, pF, pH, T)
                C = minip_c.calculate_organic_C(OM)
                r_conv_N = r_conv_C * (Norg/C)
                return r_conv_N

            def calculate_assimilation_rate_N(self, janssen, minip_c, OM, f_ass_dis, f_C_N_microbial, t, pF, pH, T):
                r_ass_C = minip_c.calculate_assimilation_rate(janssen, OM, f_ass_dis, t, pF, pH, T)
                r_ass_N = r_ass_C/f_C_N_microbial
                return r_ass_N

            def calculate_dissimilation_rate_N(self, janssen, minip_c, OM, Norg, f_ass_dis, f_C_N_microbial, t, pF, pH, T):
                r_ass_N = self.calculate_assimilation_rate_N(janssen, minip_c, OM, f_ass_dis, f_C_N_microbial, t, pF, pH, T)
                r_conv_N = self.calculate_total_conversion_rate_N(janssen, minip_c, OM, Norg, f_ass_dis, t, pF, pH, T)
                r_diss_N = r_conv_N - r_ass_N
                return r_diss_N
