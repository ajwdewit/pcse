# -*- coding: utf-8 -*-
# Copyright (c) 2021 Wageningen Environmental Research
# Allard de Wit (allard.dewit@wur.nl), July 2021
# Approach based on LINGRA N made by Joost Wolf
from collections import namedtuple

import pcse
from pcse import exceptions as exc
from pcse.traitlets import Float, Int, Instance, Bool
from pcse.decorators import prepare_rates, prepare_states
from pcse.base import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject
from pcse.util import AfgenTrait, limit

MaxNutrientConcentrations = namedtuple("MaxNutrientConcentrations", ["NMAXLV", "NMAXRT",])


class N_Demand_Uptake(SimulationObject):
    """Calculates the crop N demand and its uptake from the soil.

    Crop N demand is calculated as the difference between the
    actual N concentration (kg N per kg biomass) in the
    vegetative plant organs (leaves, stems and roots) and the maximum
    N concentration for each organ. N uptake is then estimated
    as the minimum of supply from the soil and demand from the crop.

    **Simulation parameters**

    ============  ============================================= =======  ======================
     Name          Description                                   Type     Unit
    ============  ============================================= =======  ======================
    NMAXLV_TB      Maximum N concentration in leaves as          TCr     kg N kg-1 dry biomass
                   function of DVS
    NMAXRT_FR      Maximum N concentration in roots as fraction  SCr     -
                   of maximum N concentration in leaves
    ============  ============================================= ======= =======================


    **Rate variables**

    ===========  ================================================= ==== ================
     Name         Description                                      Pbl      Unit
    ===========  ================================================= ==== ================
    RNuptakeLV     Rate of N uptake in leaves                        Y   |kg N ha-1 d-1|
    RNuptakeRT     Rate of N uptake in roots                         Y   |kg N ha-1 d-1|

    RNuptake       Total rate of N uptake                            Y   |kg N ha-1 d-1|
    NdemandLV      Ndemand of leaves based on current growth rate    N   |kg N ha-1|
                   and deficienties from previous time steps
    NdemandRT      N demand of roots, idem as leaves                 N   |kg N ha-1|

    Ndemand        Total N demand (leaves + roots)                   N   |kg N ha-1|
    ===========  ================================================= ==== ================

    **Signals send or handled**

    None

    **External dependencies**

    ================  =================================== ====================  ===========
     Name              Description                         Provided by            Unit
    ================  =================================== ====================  ===========
    DVS               Crop development stage              DVS_Phenology              -
    NAVAIL            Total available N from soil         NPK_Soil_Dynamics      |kg ha-1|
    ================  =================================== ====================  ===========

    """

    class Parameters(ParamTemplate):
        NMAXLV_TB = AfgenTrait()  # maximum N concentration in leaves as function of dvs
        NMAXRT_FR = Float(-99.)  # maximum N concentration in roots as fraction of maximum N concentration in leaves
        NUPTAKE_MAX = Float(-99)

    class RateVariables(RatesTemplate):
        RNuptakeLV = Float(-99.)  # N uptake rate [kg ha-1 d -1]
        RNuptakeRT = Float(-99.)
        RNuptake = Float(-99.)  # Total N uptake rate [kg ha-1 d -1]

        NdemandLV = Float(-99.)
        NdemandRT = Float(-99.)
        Ndemand = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: a ParameterProvider with parameter key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.kiosk = kiosk

        self.rates = self.RateVariables(kiosk, publish=["RNuptakeLV", "RNuptakeRT", "RNuptake", ])

    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states
        p = self.params
        k = self.kiosk

        delt = 1.0
        mc = self._compute_N_max_concentrations()

        # Total NPK demand of leaves, stems, roots and storage organs
        # Demand consists of a demand carried over from previous timesteps plus a demand from new growth
        # Note that we are pre-integrating here, so a multiplication with time-step delt is required

        # N demand [kg ha-1]
        r.NdemandLV = max(mc.NMAXLV * k.WeightLVgreen - k.NamountLV, 0.) + max(k.LVgrowth * mc.NMAXLV, 0) * delt
        r.NdemandRT = max(mc.NMAXRT * k.WeightRT - k.NamountRT, 0.) + max(k.dWeightRT * mc.NMAXRT, 0) * delt
        r.Ndemand = r.NdemandLV + r.NdemandRT

        # No nutrients are absorbed
        # when severe water shortage occurs i.e. RFTRA <= 0.01
        NutrientLIMIT = 1.0 if k.RFTRA > 0.01 else 0.

        # N uptake rate
        # if no demand then uptake rate = 0.
        if r.Ndemand == 0.:
            r.RNuptake = r.RNuptakeLV = r.RNuptakeRT = 0.
        else:
            # N uptake rate from soil, with a maximum equal to NUPTAKE_MAX
            RNuptake = (max(0., min(r.Ndemand, k.NAVAIL)) * NutrientLIMIT)
            r.RNuptake = min(RNuptake, p.NUPTAKE_MAX)
            # Distribution over roots/leaves
            r.RNuptakeLV = (r.NdemandLV / r.Ndemand) * r.RNuptake
            r.RNuptakeRT = (r.NdemandRT / r.Ndemand) * r.RNuptake

    @prepare_states
    def integrate(self, day, delt=1.0):
        pass

    def _compute_N_max_concentrations(self):
        """Computes the maximum N concentrations in leaves, stems, roots and storage organs.

        Note that max concentrations are first derived from the dilution curve for leaves.
        Maximum concentrations for stems and roots are computed as a fraction of the
        concentration for leaves. Maximum concentration for storage organs is directly taken from
        the parameters NMAXSO.
        """

        p = self.params
        k = self.kiosk
        NMAXLV = p.NMAXLV_TB(k.DVS)
        max_NPK_conc = MaxNutrientConcentrations(
            # Maximum NPK concentrations in leaves [kg N kg-1 DM]
            NMAXLV=NMAXLV,
            NMAXRT=p.NMAXRT_FR * NMAXLV,
        )

        return max_NPK_conc


class N_Stress(SimulationObject):
    """Implementation of N stress calculation through nitrogen nutrition index.

    Stress factors are calculated based on the mass concentrations of N in
    the vegetative biomass of the plant. For each pool of nutrients, four
    concentrations are calculated based on the biomass for leaves and stems:
    - the actual concentration based on the actual amount of nutrients
      divided by the vegetative biomass.
    - The maximum concentration, being the maximum that the plant can absorb
      into its leaves and stems.
    - The critical concentration, being the concentration that is needed to
      maintain growth rates that are not limited by N (regulated by NCRIT_FR).
      For N, the critical concentration can be lower than the maximum
      concentration. This concentration is sometimes called 'optimal
      concentration'.
    - The residual concentration which is the amount that is locked
      into the plant structural biomass and cannot be mobilized anymore.

    The stress index (SI) is determined as a simple ratio between those
    concentrations according to:

    :math:`SI = (C_{a) - C_{r})/(C_{c} - C_{r})`

    with subscript `a`, `r` and `c` being the actual, residual and critical
    concentration for the nutrient. This results in the nitrogen nutrition index
    (NNI). Finally, the reduction factor for assimilation (RFNUTR) is calculated using the
    reduction factor for light use efficiency (NLUE).

    **Simulation parameters**

    ============  ============================================= =======  ======================
     Name          Description                                   Type     Unit
    ============  ============================================= =======  ======================
    NMAXLV_TB      Maximum N concentration in leaves as          TCr     kg N kg-1 dry biomass
                   function of DVS
    NMAXRT_FR      Maximum N concentration in roots as fraction  SCr     -
                   of maximum N concentration in leaves

    NCRIT_FR       Critical N concentration as fraction of       SCr     -
                   maximum N concentration for vegetative
                   plant organs as a whole (leaves + stems)
    NRESIDLV       Residual N fraction in leaves                 SCr     kg N kg-1 dry biomass
    NLUE           Impact of N stress on Light use efficiency    SCr     -
    ============  ============================================= ======= =======================

    **Rate variables**

    The rate variables here are not real rate variables in the sense that they are derived
    state variables and do not represent a rate. However, as they are directly used
    in the rate variable calculation it is logical to put them here.

    =======  ================================================= ==== ==============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ==============
    NNI       Nitrogen nutrition index                          Y     -
    RFNUTR    Reduction factor for light use efficiency         Y     -
    =======  ================================================= ==== ==============


    **External dependencies:**

    ==============  =================================== =====================  ==============
     Name            Description                         Provided by            Unit
    ==============  =================================== =====================  ==============
    DVS              Crop development stage              DVS_Phenology            -
    WST              Dry weight of living stems          WOFOST_Stem_Dynamics   |kg ha-1|
    WeightLVgreen    Dry weight of living leaves         WOFOST_Leaf_Dynamics   |kg ha-1|
    NamountLV        Amount of N in leaves               N_Crop_Dynamics        |kg ha-1|
    ==============  =================================== =====================  ==============
    """

    class Parameters(ParamTemplate):
        NMAXLV_TB = AfgenTrait()  # maximum N concentration in leaves as function of dvs
        NCRIT_FR = Float(-99.)  # optimal N concentration as fraction of maximum N concentration
        NRESIDLV = Float(-99.)  # residual N fraction in leaves [kg N kg-1 dry biomass]
        NLUE = Float()

    class RateVariables(RatesTemplate):
        NNI = Float()
        RFNUTR = Float()

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: current date
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: ParameterProvider with parameter key/value pairs
        """

        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish=["NNI", "RFNUTR"])

    @prepare_rates
    def __call__(self, day, drv):

        p = self.params
        r = self.rates
        k = self.kiosk

        # Maximum N concentrations in leaves (kg N kg-1 DM)
        NMAXLV = p.NMAXLV_TB(k.DVS)

        # Total vegetative living above-ground biomass (kg DM ha-1)
        VBM = k.WeightLVgreen

        # Critical (Optimal) N amount in vegetative above-ground living biomass
        # and its N concentration
        NcriticalLV = p.NCRIT_FR * NMAXLV * VBM

        # if above-ground living biomass = 0 then optimum = 0
        if VBM > 0.:
            NcriticalVBM = NcriticalLV / VBM
            NconcentrationVBM = (k.NamountLV) / VBM
            NresidualVBM = (k.WeightLVgreen * p.NRESIDLV) / VBM
        else:
            NcriticalVBM = 0.
            NconcentrationVBM = 0.
            NresidualVBM = 0.

        if (NcriticalVBM - NresidualVBM) > 0.:
            r.NNI = limit(0.001, 1.0, (NconcentrationVBM - NresidualVBM) / (NcriticalVBM - NresidualVBM))
        else:
            r.NNI = 0.001

        r.RFNUTR = limit(0., 1.0, 1. - p.NLUE * (1.0 - r.NNI) ** 2)

        return r.NNI


class N_Crop_Dynamics(SimulationObject):
    """Implementation of overall N crop dynamics.

    NPK_Crop_Dynamics implements the overall logic of N book-keeping within the
    crop.

    **Simulation parameters**

    =============  ================================================ =======  ======================
     Name           Description                                      Type     Unit
    =============  ================================================ =======  ======================
    NMAXLV_TB      Maximum N concentration in leaves as            TCr     kg N kg-1 dry biomass
                   function of dvs
    NMAXRT_FR      Maximum N concentration in roots as fraction    SCr     -
    NRESIDLV       Residual N fraction in leaves                   SCr     kg N kg-1 dry biomass
    NRESIDRT       Residual N fraction in roots                    SCr     kg N kg-1 dry biomass
    =============  ================================================ =======  ======================

    **State variables**

    ==========  ================================================= ==== ============
     Name        Description                                      Pbl      Unit
    ==========  ================================================= ==== ============
    NamountLV    Actual N amount in living leaves                   Y   |kg N ha-1|
    NamountRT    Actual N amount in living roots                    Y   |kg N ha-1|
    Nuptake_T    total absorbed N amount                            N   |kg N ha-1|
    Nlosses_T    Total N amount lost due to senescence              N   |kg N ha-1|
    ==========  ================================================= ==== ============

    **Rate variables**

    ===========  ================================================= ==== ============
     Name         Description                                      Pbl      Unit
    ===========  ================================================= ==== ============
    RNamountLV     Weight increase (N) in leaves                    N   |kg ha-1 d-1|
    RNamountRT     Weight increase (N) in roots                     N   |kg ha-1 d-1|
    RNdeathLV      Rate of N loss in leaves                         N   |kg ha-1 d-1|
    RNdeathRT      Rate of N loss in roots                          N   |kg ha-1 d-1|
    RNloss         N loss due to senescence                         N   |kg ha-1 d-1|
    ===========  ================================================= ==== ============

    **Signals send or handled**

    None

    **External dependencies**

    =======  =================================== ====================  ============
     Name     Description                         Provided by            Unit
    =======  =================================== ====================  ============
    LVdeath     Death rate of leaves                WOFOST_Leaf_Dynamics  |kg ha-1 d-|
    =======  =================================== ====================  ============
    """

    WeightLV_remaining = Float()
    _flag_MOWING = Bool(False)

    demand_uptake = Instance(SimulationObject)
    NamountLVI = Float(-99.)  # initial soil N amount in leaves
    NamountRTI = Float(-99.)  # initial soil N amount in roots

    class Parameters(ParamTemplate):
        NMAXLV_TB = AfgenTrait()
        NMAXRT_FR = Float(-99.)
        NRESIDLV = Float(-99.)  # residual N fraction in leaves [kg N kg-1 dry biomass]
        NRESIDRT = Float(-99.)  # residual N fraction in roots [kg N kg-1 dry biomass]

    class StateVariables(StatesTemplate):
        NamountLV = Float(-99.)  # N amount in leaves [kg N ha-1]
        NamountRT = Float(-99.)  # N amount in roots [kg N ]
        Nuptake_T = Float(-99.)  # total absorbed N amount [kg N ]
        Nlosses_T = Float(-99.)

    class RateVariables(RatesTemplate):
        RNamountLV = Float(-99.)  # Net rates of NPK in different plant organs
        RNamountRT = Float(-99.)
        RNdeathLV = Float(-99.)  # N loss rate leaves [kg ha-1 d-1]
        RNharvestLV = Float()  # N loss due to harvesting [kg ha-1 d-1]
        RNdeathRT = Float(-99.)  # N loss rate roots  [kg ha-1 d-1]
        RNloss = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: dictionary with parameters as key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk)
        self.kiosk = kiosk

        # Initialize components of the npk_crop_dynamics
        self.demand_uptake = N_Demand_Uptake(day, kiosk, parvalues)

        # INITIAL STATES
        params = self.params
        k = kiosk

        # Initial amounts
        self.NamountLVI = NamountLV = k.WeightLVgreen * params.NMAXLV_TB(k.DVS)
        self.NamountRTI = NamountRT = k.WeightRT * params.NMAXLV_TB(k.DVS) * params.NMAXRT_FR

        self.states = self.StateVariables(kiosk, publish=["NamountLV", "NamountRT"],
                                          NamountLV=NamountLV,NamountRT=NamountRT,Nuptake_T=0., Nlosses_T=0.)
        self._connect_signal(self._on_MOWING, signal=pcse.signals.mowing)

    @prepare_rates
    def calc_rates(self, day, drv):
        rates = self.rates
        params = self.params
        k = self.kiosk
        s = self.states

        self.demand_uptake.calc_rates(day, drv)

        if self._flag_MOWING is True:
            # Compute loss of N due to harvesting
            rates.RNharvestLV = k.dWeightHARV/k.WeightLVgreen * s.NamountLV
            rates.RNdeathLV = 0.0
        else:
            # Compute loss of N due to death of plant material
            rates.RNdeathLV = params.NRESIDLV * k.LVdeath
            rates.RNharvestLV = 0.0
        rates.RNdeathRT = 0.0
        rates.RNloss = rates.RNdeathLV + rates.RNdeathRT + rates.RNharvestLV

        # N rates in leaves and root computed as uptake - death - harvesting.
        rates.RNamountLV = k.RNuptakeLV - rates.RNdeathLV - rates.RNharvestLV
        rates.RNamountRT = k.RNuptakeRT - rates.RNdeathRT

        self._check_N_balance(day)
        self._flag_MOWING = False

    @prepare_states
    def integrate(self, day, delt=1.0):
        rates = self.rates
        states = self.states
        k = self.kiosk

        # N amount in leaves, stems, root and storage organs
        states.NamountLV += rates.RNamountLV
        states.NamountRT += rates.RNamountRT

        self.demand_uptake.integrate(day, delt)

        # total N uptake from soil
        states.Nuptake_T += k.RNuptake
        # total N losses from dying material
        states.Nlosses_T += rates.RNloss

    def _check_N_balance(self, day):
        s = self.states
        checksum = abs(s.Nuptake_T + (self.NamountLVI + self.NamountRTI) -
                       (s.NamountLV + s.NamountRT + s.Nlosses_T))

        if abs(checksum) >= 1.0:
            msg = "N flows not balanced on day %s\n" % day
            msg += "Checksum: %f, Nuptake_T: %f\n" % (checksum, s.Nuptake_T)
            msg += "NamountLVI: %f, NamountRTI: %f\n" % \
                   (self.NamountLVI, self.NamountRTI)
            msg += "NamountLV: %f, NamountRT: %f, \n" % \
                   (s.NamountLV, s.NamountRT)
            msg += "NLOSST: %f\n" % (s.Nlosses_T)
            raise exc.NutrientBalanceError(msg)

    def _on_MOWING(self, biomass_remaining):
        """Handler for grass mowing events
        """
        self.WeightLV_remaining = biomass_remaining
        self._flag_MOWING = True